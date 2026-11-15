//
//  EELSParallel.swift
//  Role: the EELS pixel loop, which is SpectroCore's `SpectroParallel` (shared with the EDX quantification maps since #106): chunks on
//        all cores, a wave at a time, cancellation asked between waves and before every chunk, results in pixel order.
//

import Foundation
import SpectroCore

package nonisolated enum EELSParallel {
    /// `SpectroParallel.forEachPixel`: runs `body(p)` for every pixel, throws `SpectroCancelled` when `cancelled()` is true.
    package static func forEachPixel(
        count: Int, chunk: Int = 256,
        cancelled: @Sendable () -> Bool, progress: (@Sendable (Double) -> Void)?,
        body: @Sendable (Int) -> Void
    ) throws {
        try SpectroParallel.forEachPixel(count: count, chunk: chunk, cancelled: cancelled, progress: progress, body: body)
    }
}

/// The pixel loop's output buffer, SpectroCore's `SpectroSlots`.
package typealias EELSSlots<T> = SpectroSlots<T>
