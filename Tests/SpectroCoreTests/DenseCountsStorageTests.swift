//
//  DenseCountsStorageTests.swift
//  SpectroKit#95: a dense EDS spectrum image whose counts all fit 16 bits is held at two bytes a count, as HyperSpy keeps the file's
//  dtype, not four. Sums, window sums, the line profile's reads and the smoothing input stay the same exact numbers.
//  Red before: `bytesPerCount` did not exist and every dense cube held UInt32.
//

import Foundation
import Testing
import SpectroCore

struct DenseCountsStorageTests {

    /// A seeded 7 x 5 x 64 cube of counts below `limit`.
    private func cube(limit: UInt32, seed: UInt64) -> [UInt32] {
        var state = seed
        return (0 ..< 7 * 5 * 64).map { _ in
            state = state &* 6364136223846793005 &+ 1442695040888963407
            return UInt32(truncatingIfNeeded: state >> 33) % limit
        }
    }

    /// The exact sums, written out plainly.
    private func reference(_ counts: [UInt32], mask: PixelMask?, channels: Int) -> [UInt64] {
        var out = [UInt64](repeating: 0, count: channels)
        for p in 0 ..< counts.count / channels where mask?[p] ?? true {
            for c in 0 ..< channels { out[c] += UInt64(counts[p * channels + c]) }
        }
        return out
    }

    @Test func sixteenBitCountsAreHeldAtTwoBytes() {
        let counts = cube(limit: 65_536, seed: 95)
        let image = DenseSpectrumImage(ny: 5, nx: 7, channels: 64, counts: counts)
        #expect(image.bytesPerCount == 2)
        #expect(image.counts == counts)
        let direct = DenseSpectrumImage(ny: 5, nx: 7, channels: 64, counts16: counts.map { UInt16($0) })
        #expect(direct.bytesPerCount == 2)
        #expect(direct.counts == counts)
    }

    @Test func aCountAboveSixteenBitsKeepsFourBytes() {
        var counts = cube(limit: 1000, seed: 96)
        counts[17] = 65_536
        let image = DenseSpectrumImage(ny: 5, nx: 7, channels: 64, counts: counts)
        #expect(image.bytesPerCount == 4)
        #expect(image.count(at: 17) == 65_536)
    }

    @Test func theSumsAreTheSameExactNumbers() throws {
        let mask = (0 ..< 35).map { $0 % 3 != 0 }
        let ranges = [0 ..< 10, 5 ..< 40, 60 ..< 70, 30 ..< 30]
        for limit: UInt32 in [65_536, 4_000_000] {
            let counts = cube(limit: limit, seed: 97)
            let image = DenseSpectrumImage(ny: 5, nx: 7, channels: 64, counts: counts)
            #expect(image.sum(mask: nil) == reference(counts, mask: nil, channels: 64))
            #expect(image.sum(mask: mask) == reference(counts, mask: mask, channels: 64))
            let windows = image.windowSums(ranges)
            for (r, range) in ranges.enumerated() {
                for p in 0 ..< 35 {
                    let lo = max(range.lowerBound, 0), hi = min(range.upperBound, 64)
                    let expected = lo < hi ? (lo ..< hi).reduce(UInt64(0)) { $0 + UInt64(counts[p * 64 + $1]) } : 0
                    #expect(windows[r][p] == expected)
                }
            }
            #expect(image.doubles() == counts.map { Double($0) })
            var acc = [Double](repeating: 0, count: 64)
            image.accumulate(pixel: 3, weight: 0.5, into: &acc)
            #expect(acc == (0 ..< 64).map { 0.5 * Double(counts[3 * 64 + $0]) })
        }
    }
}
