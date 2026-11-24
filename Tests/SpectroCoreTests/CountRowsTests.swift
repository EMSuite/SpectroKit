//
//  CountRowsTests.swift
//  SpectroKit#134: a spectrum image hands out the counts of a band of scan rows in C order (row, column, channel), so a writer can
//  stream the cube a band at a time (Electron Microscopy Studio's HyperSpy export) without holding it twice. Dense at 2 and 4 bytes
//  a count, sparse (absent channels 0) and loaded all give the same exact counts; rows outside the image, and a conformer that
//  does not say how to read its rows, are refused by name.
//  Red before: `countRows` did not exist.
//

import Foundation
import Testing
import SpectroCore

struct CountRowsTests {
    private static let ny = 4, nx = 3, channels = 5

    /// A planted cube: every count names its row, column and channel; every third entry 0, so the sparse store has gaps.
    private static let counts: [UInt32] = (0 ..< ny * nx * channels).map { (i: Int) -> UInt32 in
        if i % 3 == 0 { return 0 }
        let row: Int = i / (nx * channels)
        let column: Int = (i / channels) % nx
        let channel: Int = i % channels
        let named = UInt32(1000 * row + 10 * column + channel)
        let wide: UInt32 = i % 2 == 0 ? 0 : 70_000
        return named + wide
    }

    /// The same cube as a CSR store holding only the nonzero entries.
    private static func sparse() throws -> SparseSpectrumImage {
        var offsets = [0], index: [UInt16] = [], values: [UInt32] = []
        for p in 0 ..< ny * nx {
            for c in 0 ..< channels where counts[p * channels + c] != 0 {
                index.append(UInt16(c))
                values.append(counts[p * channels + c])
            }
            offsets.append(values.count)
        }
        return try SparseSpectrumImage(ny: ny, nx: nx, channels: channels, rowOffsets: offsets, channelIndex: index, counts: values)
    }

    private static func band(_ rows: Range<Int>) -> [UInt32] {
        Array(counts[(rows.lowerBound * nx * channels) ..< (rows.upperBound * nx * channels)])
    }

    /// Mutation: the band's start not at its first row, the sparse store's entries put in the wrong channel or a gap not 0 - red.
    @Test func everyStoreHandsOutTheSameBand() throws {
        let narrow = DenseSpectrumImage(ny: Self.ny, nx: Self.nx, channels: Self.channels, counts16: Self.counts.map { UInt16($0 % 65_536) })
        let wide = DenseSpectrumImage(ny: Self.ny, nx: Self.nx, channels: Self.channels, counts: Self.counts)
        #expect(wide.bytesPerCount == 4)
        let meta = SpectrumImageMetadata(fileName: "planted", filePath: "planted", scanWidth: Self.nx, scanHeight: Self.ny,
                                         channelCount: Self.channels, energyOffsetEV: 0, energyDispersionEV: 10)
        let loaded = LoadedSpectrumImage(image: try Self.sparse(), metadata: meta,
                                         energyAxis: EnergyAxis(offset: 0, scale: 0.01, size: Self.channels))
        for rows in [0 ..< 4, 1 ..< 3, 3 ..< 4, 2 ..< 2] {
            #expect(try wide.countRows(rows) == Self.band(rows))
            #expect(try Self.sparse().countRows(rows) == Self.band(rows))
            #expect(try loaded.countRows(rows) == Self.band(rows))
            #expect(try narrow.countRows(rows) == Self.band(rows).map { $0 % 65_536 })
        }
    }

    /// Mutation: the rows not checked, or a bare conformer handing out zeros - red.
    @Test func rowsOutsideTheImageAndABareConformerAreRefused() throws {
        let wide = DenseSpectrumImage(ny: Self.ny, nx: Self.nx, channels: Self.channels, counts: Self.counts)
        #expect(throws: SpectrumRowsError.rowsOutsideImage(lower: 3, upper: 5, ny: 4)) { try wide.countRows(3 ..< 5) }
        #expect(throws: SpectrumRowsError.rowsOutsideImage(lower: -1, upper: 1, ny: 4)) { try Self.sparse().countRows(-1 ..< 1) }
        #expect(throws: SpectrumRowsError.notReadableByRows) { try Bare().countRows(0 ..< 1) }
    }
}

/// A conformer that keeps its pixels to itself.
private struct Bare: SpectrumImage {
    var ny: Int { 1 }
    var nx: Int { 1 }
    var channels: Int { 1 }
    func sum(mask: PixelMask?) -> [UInt64] { [0] }
    func sum(mask: PixelMask?, cancelled: @Sendable () -> Bool, progress: (@Sendable (Double) -> Void)?) throws -> [UInt64] { [0] }
    func windowSums(_ ranges: [Range<Int>]) -> [[UInt64]] { ranges.map { _ in [0] } }
    func windowSums(_ ranges: [Range<Int>], cancelled: @Sendable () -> Bool, progress: (@Sendable (Double) -> Void)?) throws -> [[UInt64]] {
        ranges.map { _ in [0] }
    }
}
