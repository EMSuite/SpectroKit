//
//  LeastSquaresFit.swift
//  Role: Unweighted least squares with non-negative line areas: the NAMED DEFAULT of the
//        spectrum fit. Every result carries `LeastSquaresFit.label`.
//
//  Lawson and Hanson, "Solving Least Squares Problems" (1974), Algorithm 23.10 (NNLS), with
//  the passive-set subproblem solved by Accelerate's QR-based `dgelsy_` (see
//  FitLinearAlgebra). Columns of a design that may take either sign (a polynomial or
//  continuum background) are eliminated exactly first (variable projection): with
//  F the free block, P_F the projector onto its span, NNLS runs on ((I - P_F) N, (I - P_F) y)
//  and the free coefficients follow by ordinary least squares from the residual. This is the same
//  minimiser as a single NNLS with F split into +F and -F, without the singular pair.
//
//  DEVIATION (eXSpy/hyperspy): `Model1D.fit` with its default "ls" loss is an UNBOUNDED
//  Levenberg-Marquardt on the line areas (`ext_force_positive` acts only with `bounded=True`),
//  and it stops short of the least-squares minimum: on EDS_TEM_FePt_nanoparticles it returns
//  2800.0213 / 14921.183 where the exact minimiser is 2800.2296 / 14920.805 (7e-5 and 3e-5
//  relative). Here the minimum is exact and areas cannot go negative.
//

import Foundation

package nonisolated enum LeastSquaresFit {
    /// Printed beside every least-squares result.
    package static let label = "least squares, unweighted"

    package struct Solution: Sendable, Equatable {
        /// Coefficients of the non-negative columns (>= 0).
        package var nonNegative: [Double]
        /// Coefficients of the free-sign columns.
        package var free: [Double]
        package var iterations: Int
        package var converged: Bool
        /// False when the data or the answer is not finite or a free-block solve failed (then `converged` is false too).
        package var answered: Bool = true
    }

    /// `min || sqrt(w) (y - N x - F b) ||` over `x >= 0`, `b` free. `rowWeights` nil is the unweighted fit.
    package static func solve(
        nonNegative n: ColumnMatrix, free f: ColumnMatrix = ColumnMatrix(rows: 0, cols: 0),
        y: [Double], rowWeights: [Double]? = nil
    ) -> Solution {
        let m = y.count
        precondition(n.rows == m && (f.cols == 0 || f.rows == m))
        let sw = rowWeights?.map { $0.squareRoot() }
        let nW = sw.map { n.scalingRows($0) } ?? n
        let fW = f.cols > 0 ? (sw.map { f.scalingRows($0) } ?? f) : f
        var yW = y
        if let sw { for i in 0..<m { yW[i] *= sw[i] } }

        // Column scaling of N to unit norm; a column with no support stays at 0 and never enters.
        var norms = [Double](repeating: 0, count: n.cols)
        var active: [Int] = []
        var scaled = nW
        for j in 0..<n.cols {
            var s = 0.0
            for i in 0..<m { let v = nW.data[j * m + i]; s += v * v }
            let nj = s.squareRoot()
            norms[j] = nj
            if nj > 1e-300 {
                active.append(j)
                for i in 0..<m { scaled.data[j * m + i] /= nj }
            }
        }

        // Eliminate the free block.
        // Every column active (the usual case) needs no second copy of the design.
        var work = active.count == n.cols ? scaled : scaled.selecting(active)
        var target = yW
        // A free block the solver cannot answer, or a NaN in the data, is not a converged fit (it used to report converged zeros).
        var answered = yW.allSatisfy(\.isFinite)
        guard answered else {
            return Solution(nonNegative: [Double](repeating: 0, count: n.cols), free: [Double](repeating: 0, count: f.cols),
                            iterations: 0, converged: false, answered: false)
        }
        if fW.cols > 0 {
            var rhs = ColumnMatrix(rows: m, cols: work.cols + 1)
            rhs.data = work.data + yW
            if let (x, _) = FitLinearAlgebra.leastSquares(fW, rhs) {
                let proj = fW.times1(x)
                for k in 0..<rhs.data.count { rhs.data[k] -= proj[k] }
            } else { answered = false }
            work = ColumnMatrix(rows: m, cols: work.cols, data: Array(rhs.data[0..<(work.cols * m)]))
            target = Array(rhs.data[(work.cols * m)...])
        }

        let (xa, iterations, converged) = lawsonHanson(work, target)
        var x = [Double](repeating: 0, count: n.cols)
        for (k, j) in active.enumerated() { x[j] = xa[k] / norms[j] }

        // Free coefficients from the residual of the non-negative part.
        var b = [Double](repeating: 0, count: f.cols)
        if fW.cols > 0 {
            let nx = nW.times(x)
            var r = yW
            for i in 0..<m { r[i] -= nx[i] }
            if let (sol, _) = FitLinearAlgebra.leastSquares(fW, ColumnMatrix(rows: m, cols: 1, data: r)) { b = sol.data } else { answered = false }
        }
        if !(x.allSatisfy(\.isFinite) && b.allSatisfy(\.isFinite)) { answered = false }
        return Solution(nonNegative: x, free: b, iterations: iterations, converged: converged && answered, answered: answered)
    }

    /// Algorithm 23.10 on unit-norm columns. Returns the solution on those columns.
    ///
    /// Every buffer of the loop (gradient, model, passive index list, passive solution and the LAPACK workspace) is made once here
    /// and reused by each iteration (macEM#641: allocating them per iteration made 5.35 GB of transient arrays per EDS Auto ID).
    /// The arithmetic is unchanged: the same passive columns go to the same `dgelsy_` call in the same order.
    private static func lawsonHanson(_ a: ColumnMatrix, _ y: [Double]) -> ([Double], Int, Bool) {
        let m = a.rows, n = a.cols
        var x = [Double](repeating: 0, count: n)
        if n == 0 || m == 0 { return (x, 0, true) }
        var passive = [Bool](repeating: false, count: n)
        let yNorm = y.reduce(0) { $0 + $1 * $1 }.squareRoot()
        let tol = 1e-12 * max(yNorm, 1e-300)
        var residual = y
        var w = [Double](repeating: 0, count: n)
        var ax = [Double](repeating: 0, count: m)
        var idx: [Int] = []
        idx.reserveCapacity(n)
        var s = [Double](repeating: 0, count: n)
        let workspace = LeastSquaresWorkspace(rows: m, maxCols: n)
        var outer = 0
        let maxOuter = 3 * n + 100
        var converged = false
        while outer < maxOuter {
            outer += 1
            // Gradient on the zero set.
            a.transposeTimes(residual, into: &w)
            var t = -1; var best = tol
            for j in 0..<n where !passive[j] && w[j] > best { best = w[j]; t = j }
            if t < 0 { converged = true; break }
            passive[t] = true
            // Inner loop: keep the passive solution positive.
            var guardInner = 0
            while guardInner < 3 * n + 10 {
                guardInner += 1
                idx.removeAll(keepingCapacity: true)
                for j in 0..<n where passive[j] { idx.append(j) }
                let k = idx.count
                guard FitLinearAlgebra.leastSquares(a, columns: idx, y, solution: &s, workspace: workspace) else {
                    passive[t] = false; break
                }
                var allPositive = true
                for c in 0..<k where !(s[c] > 0) { allPositive = false; break }
                if allPositive {
                    for (c, j) in idx.enumerated() { x[j] = s[c] }
                    break
                }
                // Step towards s as far as feasibility allows.
                var alpha = Double.infinity
                for (c, j) in idx.enumerated() where s[c] <= 0 {
                    let denom = x[j] - s[c]
                    if denom > 0 { alpha = min(alpha, x[j] / denom) }
                }
                if !alpha.isFinite { alpha = 0 }
                for (c, j) in idx.enumerated() { x[j] += alpha * (s[c] - x[j]) }
                for j in idx where x[j] <= 1e-14 { x[j] = 0; passive[j] = false }
            }
            a.times(x, into: &ax)
            for i in 0..<m { residual[i] = y[i] - ax[i] }
        }
        return (x, outer, converged)
    }
}

extension ColumnMatrix {
    /// `A X` for a coefficient matrix `X` (cols x k), column-major.
    fileprivate nonisolated func times1(_ x: ColumnMatrix) -> [Double] {
        precondition(x.rows == cols)
        var out = [Double](repeating: 0, count: rows * x.cols)
        for c in 0..<x.cols {
            let col = Array(x.data[(c * cols)..<((c + 1) * cols)])
            let y = times(col)
            for i in 0..<rows { out[c * rows + i] = y[i] }
        }
        return out
    }
}
