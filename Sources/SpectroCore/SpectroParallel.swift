//
//  SpectroParallel.swift
//  Role: the one pixel loop of the per-pixel spectroscopy calls (EELS fits and maps, EDX quantification maps): pixels are independent,
//        so the loop runs in chunks on all cores, a wave at a time, and each chunk writes only its own slots (`SpectroSlots`), so an
//        answer never depends on the schedule. Cancellation is asked between waves on the calling thread and again before every
//        chunk: a chunk that sees `cancelled()` true stops the wave's remaining chunks. The default closure
//        (`SpectroCancellation.task`) reads Swift task cancellation, which only the calling thread carries (it runs a share of
//        every wave's chunks), so a cancel lands within about one chunk of that thread's work rather than a whole wave (#106).
//        CPU float64 only (no Metal: the work per pixel is smaller than a copy).
//

import Foundation
import os

package nonisolated enum SpectroParallel {
    /// Runs `body(p)` for every pixel `0..<count`, in chunks of `chunk` pixels, several chunks at a time. Throws `SpectroCancelled`
    /// when `cancelled()` is true between waves or before a chunk; nothing partial is kept by the callers. `body` must write only the
    /// slot of its own pixel.
    package static func forEachPixel(
        count: Int, chunk: Int = 256,
        cancelled: @Sendable () -> Bool, progress: (@Sendable (Double) -> Void)? = nil,
        body: @Sendable (Int) -> Void
    ) throws {
        try forEachChunk(count: count, chunk: chunk, cancelled: cancelled, progress: progress) { range in
            for p in range { body(p) }
        }
    }

    /// `forEachPixel` a chunk at a time, for loops that keep scratch per chunk: `body` gets each chunk's index range once.
    package static func forEachChunk(
        count: Int, chunk: Int,
        cancelled: @Sendable () -> Bool, progress: (@Sendable (Double) -> Void)? = nil,
        body: @Sendable (Range<Int>) -> Void
    ) throws {
        guard count > 0 else { progress?(1); return }
        let chunk = max(chunk, 1)
        let chunks = (count + chunk - 1) / chunk
        let wave = max(ProcessInfo.processInfo.activeProcessorCount, 1) * 2
        let stopped = OSAllocatedUnfairLock(initialState: false)
        var next = 0
        while next < chunks {
            try SpectroCancellation.check(cancelled)
            progress?(Double(next * chunk) / Double(count))
            let end = min(next + wave, chunks)
            let first = next
            DispatchQueue.concurrentPerform(iterations: end - first) { i in
                if stopped.withLock({ $0 }) { return }
                if cancelled() { stopped.withLock { $0 = true }; return }
                let ci = first + i
                body((ci * chunk)..<min(count, (ci + 1) * chunk))
            }
            if stopped.withLock({ $0 }) { throw SpectroCancelled() }
            next = end
        }
        progress?(1)
    }
}

/// A buffer several chunks write to at disjoint indices (the pixel loop's output): raw storage, so no array uniqueness check races
/// between threads. `@unchecked Sendable` because every writer owns a disjoint set of indices (one pixel's slot) and `values` is read
/// only after the loop has joined; that invariant is the callers' (`SpectroParallel.forEachPixel`'s `body` writes its own pixel).
package final class SpectroSlots<T>: @unchecked Sendable {
    private let storage: UnsafeMutablePointer<T>
    package let count: Int

    package init(count: Int, initial: T) {
        self.count = count
        storage = UnsafeMutablePointer<T>.allocate(capacity: max(count, 1))
        storage.initialize(repeating: initial, count: count)
    }

    deinit {
        storage.deinitialize(count: count)
        storage.deallocate()
    }

    /// Writes slot `i` (callers write only their own pixel's slots).
    package func set(_ i: Int, _ value: T) { storage[i] = value }
    package subscript(i: Int) -> T { storage[i] }
    package var values: [T] { Array(UnsafeBufferPointer(start: storage, count: count)) }
}
