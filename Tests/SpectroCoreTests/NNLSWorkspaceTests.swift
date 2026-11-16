//
//  NNLSWorkspaceTests.swift
//  Role: macEM#641. The NNLS active-set loop of `LeastSquaresFit` now reuses one set of buffers per solve instead of allocating
//        the passive design, the LAPACK buffers and the gradient on every iteration. The numbers must not move: these tests
//        hold the new solver to a frozen copy of the solver before the change (`Before` below), bit for bit, on designs shaped
//        like an EDS fit (overlapping Gaussian lines, a free polynomial background, noisy counts) that drive the inner loop
//        through dropped columns; and they check that a workspace reused across different column sets answers as a fresh solve.
//

import Foundation
import Testing
import SpectroCore

struct NNLSWorkspaceTests {

    private struct SplitMix64: RandomNumberGenerator {
        var state: UInt64
        mutating func next() -> UInt64 {
            state &+= 0x9E3779B97F4A7C15
            var z = state
            z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9
            z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
            return z ^ (z >> 31)
        }
        mutating func uniform() -> Double { Double(next() >> 11) / Double(1 << 53) }
    }

    /// A design of `lines` Gaussians (some nearly on top of each other) on `m` channels, a degree-`poly` polynomial as the free
    /// block, and counts from a few of the lines plus noise that can push a fitted area below zero.
    private static func problem(seed: UInt64, m: Int, lines: Int, poly: Int)
        -> (n: ColumnMatrix, f: ColumnMatrix, y: [Double], weights: [Double]) {
        var rng = SplitMix64(state: seed)
        var cols: [[Double]] = []
        for _ in 0..<lines {
            let c = rng.uniform() * Double(m)
            let s = 2 + rng.uniform() * 12
            cols.append((0..<m).map { i in let d = (Double(i) - c) / s; return exp(-0.5 * d * d) })
        }
        // A column with no support never enters (the zero-norm path).
        if lines > 3 { cols[lines / 2] = [Double](repeating: 0, count: m) }
        let free = (0...poly).map { p in (0..<m).map { i in pow(Double(i) / Double(m), Double(p)) } }
        var y = [Double](repeating: 0, count: m)
        for (j, c) in cols.enumerated() where j % 3 != 1 {
            let amp = 50 + rng.uniform() * 2000
            for i in 0..<m { y[i] += amp * c[i] }
        }
        for i in 0..<m { y[i] += 30 + 20 * Double(i) / Double(m) + (rng.uniform() - 0.5) * 60 }
        let weights = (0..<m).map { _ in 0.5 + rng.uniform() }
        return (ColumnMatrix(columns: cols), ColumnMatrix(columns: free), y, weights)
    }

    @Test func testSolverMatchesTheSolverBeforeTheChangeBitForBit() {
        var droppedSomewhere = false
        for seed in UInt64(1)...60 {
            let m = 64 + Int(seed % 7) * 97
            let lines = 3 + Int(seed % 11) * 2
            let (n, f, y, w) = Self.problem(seed: seed, m: m, lines: lines, poly: Int(seed % 4))
            for free in [ColumnMatrix(rows: 0, cols: 0), f] {
                for weights in [nil, w] {
                    let now = LeastSquaresFit.solve(nonNegative: n, free: free, y: y, rowWeights: weights)
                    let before = Before.solve(nonNegative: n, free: free, y: y, rowWeights: weights)
                    #expect(Before.Answer(now) == before, "seed \(seed), m \(m), \(lines) lines, free \(free.cols), weighted \(weights != nil)")
                    if now.iterations > now.nonNegative.filter({ $0 > 0 }).count + 1 { droppedSomewhere = true }
                }
            }
        }
        // The cases reach the inner loop's step-and-drop branch, not only the plain add-one-column path.
        #expect(droppedSomewhere)
    }

    @Test func testReusedWorkspaceAnswersAsAFreshSolve() throws {
        let (n, _, y, _) = Self.problem(seed: 7, m: 300, lines: 12, poly: 0)
        let workspace = LeastSquaresWorkspace(rows: n.rows, maxCols: n.cols)
        var x = [Double](repeating: .nan, count: n.cols)
        // Wider, then narrower, then wider sets: stale pivots or right-hand-side rows from the last call must not leak in.
        for js in [[0, 1, 2, 3, 4, 5, 6, 7], [2, 9], [], [11, 3, 0, 8, 5], Array(0..<12), [4]] {
            #expect(FitLinearAlgebra.leastSquares(n, columns: js, y, solution: &x, workspace: workspace))
            if js.isEmpty { continue }
            let fresh = try #require(FitLinearAlgebra.leastSquares(n.selecting(js), ColumnMatrix(rows: n.rows, cols: 1, data: y)))
            #expect(Array(x[0..<js.count]) == fresh.x.data, "columns \(js)")
        }
    }

    @Test func testInPlaceProductsMatchTheAllocatingOnes() {
        let (n, _, y, _) = Self.problem(seed: 3, m: 129, lines: 9, poly: 0)
        var x = (0..<n.cols).map { Double($0 % 3) * 1.5 }
        x[1] = 0
        var ax = [Double](repeating: .nan, count: n.rows)
        n.times(x, into: &ax)
        #expect(ax == n.times(x))
        var g = [Double](repeating: .nan, count: n.cols)
        n.transposeTimes(y, into: &g)
        #expect(g == n.transposeTimes(y))
    }
}

/// `LeastSquaresFit.solve` as it was before macEM#641 (SpectroKit 0.14.3), kept verbatim as the reference the reuse must not move.
private enum Before {
    struct Answer: Equatable {
        var nonNegative: [Double], free: [Double], iterations: Int, converged: Bool, answered: Bool
        init(nonNegative: [Double], free: [Double], iterations: Int, converged: Bool, answered: Bool) {
            self.nonNegative = nonNegative; self.free = free; self.iterations = iterations; self.converged = converged; self.answered = answered
        }
        init(_ s: LeastSquaresFit.Solution) {
            self.init(nonNegative: s.nonNegative, free: s.free, iterations: s.iterations, converged: s.converged, answered: s.answered)
        }
    }

    static func solve(nonNegative n: ColumnMatrix, free f: ColumnMatrix, y: [Double], rowWeights: [Double]?)
        -> Answer {
        let m = y.count
        let sw = rowWeights?.map { $0.squareRoot() }
        let nW = sw.map { n.scalingRows($0) } ?? n
        let fW = f.cols > 0 ? (sw.map { f.scalingRows($0) } ?? f) : f
        var yW = y
        if let sw { for i in 0..<m { yW[i] *= sw[i] } }
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
        var work = scaled.selecting(active)
        var target = yW
        var answered = yW.allSatisfy(\.isFinite)
        guard answered else {
            return Answer(nonNegative: [Double](repeating: 0, count: n.cols), free: [Double](repeating: 0, count: f.cols),
                          iterations: 0, converged: false, answered: false)
        }
        if fW.cols > 0 {
            var rhs = ColumnMatrix(rows: m, cols: work.cols + 1)
            rhs.data = work.data + yW
            if let (x, _) = FitLinearAlgebra.leastSquares(fW, rhs) {
                var proj = [Double](repeating: 0, count: m * x.cols)
                for c in 0..<x.cols {
                    let col = Array(x.data[(c * fW.cols)..<((c + 1) * fW.cols)])
                    let yc = fW.times(col)
                    for i in 0..<m { proj[c * m + i] = yc[i] }
                }
                for k in 0..<rhs.data.count { rhs.data[k] -= proj[k] }
            } else { answered = false }
            work = ColumnMatrix(rows: m, cols: work.cols, data: Array(rhs.data[0..<(work.cols * m)]))
            target = Array(rhs.data[(work.cols * m)...])
        }
        let (xa, iterations, converged) = lawsonHanson(work, target)
        var x = [Double](repeating: 0, count: n.cols)
        for (k, j) in active.enumerated() { x[j] = xa[k] / norms[j] }
        var b = [Double](repeating: 0, count: f.cols)
        if fW.cols > 0 {
            let nx = nW.times(x)
            var r = yW
            for i in 0..<m { r[i] -= nx[i] }
            if let (sol, _) = FitLinearAlgebra.leastSquares(fW, ColumnMatrix(rows: m, cols: 1, data: r)) { b = sol.data } else { answered = false }
        }
        if !(x.allSatisfy(\.isFinite) && b.allSatisfy(\.isFinite)) { answered = false }
        return Answer(nonNegative: x, free: b, iterations: iterations, converged: converged && answered, answered: answered)
    }

    static func lawsonHanson(_ a: ColumnMatrix, _ y: [Double]) -> ([Double], Int, Bool) {
        let m = a.rows, n = a.cols
        var x = [Double](repeating: 0, count: n)
        if n == 0 || m == 0 { return (x, 0, true) }
        var passive = [Bool](repeating: false, count: n)
        let yNorm = y.reduce(0) { $0 + $1 * $1 }.squareRoot()
        let tol = 1e-12 * max(yNorm, 1e-300)
        var residual = y
        var outer = 0
        let maxOuter = 3 * n + 100
        var converged = false
        while outer < maxOuter {
            outer += 1
            let w = a.transposeTimes(residual)
            var t = -1; var best = tol
            for j in 0..<n where !passive[j] && w[j] > best { best = w[j]; t = j }
            if t < 0 { converged = true; break }
            passive[t] = true
            var guardInner = 0
            while guardInner < 3 * n + 10 {
                guardInner += 1
                let idx = (0..<n).filter { passive[$0] }
                guard let (s, _) = FitLinearAlgebra.leastSquares(a.selecting(idx), ColumnMatrix(rows: m, cols: 1, data: y)) else {
                    passive[t] = false; break
                }
                if s.data.allSatisfy({ $0 > 0 }) {
                    for (k, j) in idx.enumerated() { x[j] = s.data[k] }
                    break
                }
                var alpha = Double.infinity
                for (k, j) in idx.enumerated() where s.data[k] <= 0 {
                    let denom = x[j] - s.data[k]
                    if denom > 0 { alpha = min(alpha, x[j] / denom) }
                }
                if !alpha.isFinite { alpha = 0 }
                for (k, j) in idx.enumerated() { x[j] += alpha * (s.data[k] - x[j]) }
                for j in idx where x[j] <= 1e-14 { x[j] = 0; passive[j] = false }
            }
            let ax = a.times(x)
            for i in 0..<m { residual[i] = y[i] - ax[i] }
        }
        return (x, outer, converged)
    }
}
