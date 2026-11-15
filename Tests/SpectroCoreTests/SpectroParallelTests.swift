//
//  SpectroParallelTests.swift
//  The shared pixel loop stops inside a wave, not only between waves (#106), and the long SpectroCore calls that had no
//  `cancelled` closure take one: true throws `SpectroCancelled`, nothing is returned. `IncompleteGamma.logGamma` (lgamma_r) is
//  `lgamma`'s value.
//

import Foundation
import Testing
import os
import SpectroCore

struct SpectroParallelTests {

    /// The first ask between waves says go, every later ask says stop: no chunk runs. Red under: the loop asking only between
    /// waves (the whole first wave runs, 2 x cores chunks).
    @Test func aCancelInsideAWaveStopsItsChunks() throws {
        let asks = OSAllocatedUnfairLock(initialState: 0)
        let ran = OSAllocatedUnfairLock(initialState: 0)
        let cancelled: @Sendable () -> Bool = { asks.withLock { $0 += 1; return $0 > 1 } }
        #expect(throws: SpectroCancelled.self) {
            try SpectroParallel.forEachPixel(count: 100_000, chunk: 1, cancelled: cancelled) { _ in ran.withLock { $0 += 1 } }
        }
        #expect(ran.withLock { $0 } == 0)
    }

    /// Never cancelled, every pixel runs once, in any order. Red under: a chunk skipped or run twice by the per-chunk check.
    @Test func everyPixelRunsOnceWhenNotCancelled() throws {
        let slots = SpectroSlots<Int>(count: 10_007, initial: 0)
        try SpectroParallel.forEachPixel(count: slots.count, chunk: 64, cancelled: { false }) { p in slots.set(p, slots[p] + 1) }
        #expect(slots.values.allSatisfy { $0 == 1 })
    }

    @Test func theLongCoreCallsTakeACancel() throws {
        let stop: @Sendable () -> Bool = { true }
        #expect(throws: SpectroCancelled.self) {
            try SpectrumSmoothing.gaussian(sigma: 1).smoothBatch([Double](repeating: 1, count: 64), spectrumLength: 8, cancelled: stop)
        }
        let r = try GridRegistration(matrix: [1, 0, 0, 0, 1, 0], scanWidth: 2, scanHeight: 2, spectrumWidth: 2, spectrumHeight: 2)
        let mask = [true, false, true, false]
        #expect(throws: SpectroCancelled.self) { try MaskTransport.transport(mask: mask, registration: r, cancelled: stop) }
        #expect(throws: SpectroCancelled.self) {
            try MaskTransport.regions([SparseRegion(pixelIndices: [0, 2])], registration: r, cancelled: stop)
        }
        let image = DenseSpectrumImage(ny: 2, nx: 2, channels: 2, counts: [1, 2, 3, 4, 5, 6, 7, 8])
        #expect(throws: SpectroCancelled.self) {
            try PoolBuilder.pool(image: image, scanMask: mask, registration: r, name: "left", kind: "drawn", cancelled: stop)
        }
        let t = try MaskTransport.transport(mask: mask, registration: r)
        #expect(throws: SpectroCancelled.self) {
            try PoolBuilder.pools(image: image, transported: t, names: ["left"], kind: "phase", validation: nil, cancelled: stop)
        }
        let path = [PixelPosition(row: 7, column: 2), PixelPosition(row: 7, column: 20)]
        for weighting in [LineProfileWeighting.hard, .bilinear] {
            #expect(throws: SpectroCancelled.self) {
                try SpectrumLineProfile.compute(image: SpectrumLineProfileTests.cube(), path: path, width: 3, weighting: weighting,
                                                cancelled: stop)
            }
        }
    }

    /// The same calls without a cancel answer as before (the default closure reads task cancellation, false here).
    @Test func theDefaultCancelLetsTheCallsRun() throws {
        let r = try GridRegistration(matrix: [1, 0, 0, 0, 1, 0], scanWidth: 2, scanHeight: 2, spectrumWidth: 2, spectrumHeight: 2)
        let image = DenseSpectrumImage(ny: 2, nx: 2, channels: 2, counts: [1, 2, 3, 4, 5, 6, 7, 8])
        let pool = try PoolBuilder.pool(image: image, scanMask: [true, false, true, false], registration: r, name: "left", kind: "drawn")
        #expect(pool.counts == [6, 8])
        let regions = try MaskTransport.regions([SparseRegion(pixelIndices: [0, 2])], registration: r)
        #expect(regions.entries.count == 1)
    }

    /// The reentrant form is the same number. Red under: a sign or branch change in the swap.
    @Test func logGammaIsLgamma() {
        for a in [0.5, 1, 2.5, 10, 37.25, 170] { #expect(IncompleteGamma.logGamma(a) == lgamma(a)) }
    }
}
