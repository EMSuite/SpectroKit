//
//  FloatCubeMethodsTests.swift
//  Role: SpectroKit#97. Quantification, the edge model (joint, windowed, background only, two-area), MLLS and the thickness with a zero-loss
//        spectrum image read a float32 cube a pixel at a time instead of widening the whole cube to Double for the call. Every read widens
//        exactly, so each method gives the same numbers, bit for bit, from the Floats as from the same values held as Doubles.
//

import Foundation
import Testing
import SpectroTestSupport
import SpectroCore
import SpectroEELS

struct FloatCubeMethodsTests {

    /// The case's core loss rounded to Float, held once as Floats and once as those Floats widened to Double.
    private func pair(_ core: CoreLossSpectrumImage) throws -> (single: CoreLossSpectrumImage, double: CoreLossSpectrumImage) {
        let floats = core.storedCounts.doubles.map { Float($0) }
        let single = try CoreLossSpectrumImage(nx: core.nx, ny: core.ny, channels: core.channels, counts: floats,
                                               energyOffsetEV: core.energyOffsetEV, energyScaleEV: core.energyScaleEV)
        let double = try CoreLossSpectrumImage(nx: core.nx, ny: core.ny, channels: core.channels, counts: floats.map { Double($0) },
                                               energyOffsetEV: core.energyOffsetEV, energyScaleEV: core.energyScaleEV)
        #expect(single.storedCounts.precision == .single)
        return (single, double)
    }

    private func pair(_ low: LowLossSpectrumImage) throws -> (single: LowLossSpectrumImage, double: LowLossSpectrumImage) {
        let floats = low.storedCounts.doubles.map { Float($0) }
        let single = try LowLossSpectrumImage(nx: low.nx, ny: low.ny, channels: low.channels, counts: floats,
                                              energyOffsetEV: low.energyOffsetEV, energyScaleEV: low.energyScaleEV)
        let double = try LowLossSpectrumImage(nx: low.nx, ny: low.ny, channels: low.channels, counts: floats.map { Double($0) },
                                              energyOffsetEV: low.energyOffsetEV, energyScaleEV: low.energyScaleEV)
        return (single, double)
    }

    private func model(_ name: String) throws -> ModelCase {
        let all = ModelFixtures.pins["cases"] as! [String: [String: Any]]
        return try ModelFixtures.make(name, all[name]!)
    }

    @Test func testQuantificationGivesTheSameMap() throws {
        let core = try pair(try EdgeFixtures.core("edge-core-R-ss-noiseless.f64"))
        let low = try pair(try EdgeFixtures.low("edge-low-R-ss-noiseless.f64"))
        let m = EdgeFixtures.plantedMicroscope
        let a = try EELSEdgeQuantification.ratioMethod(coreLoss: core.single, microscope: m, lowLoss: .integratedToWidth(low.single))
        let b = try EELSEdgeQuantification.ratioMethod(coreLoss: core.double, microscope: m, lowLoss: .integratedToWidth(low.double))
        #expect(a.pixels == b.pixels)
    }

    @Test func testTheEdgeModelGivesTheSameFits() throws {
        for name in ["cno", "fel-conv"] {
            let c = try model(name)
            let core = try pair(c.core)
            let low = try c.low.map { try pair($0) }
            let a = try EELSModelFit.fit(coreLoss: core.single, edges: c.edges, microscope: c.microscope, lowLoss: low?.single)
            let b = try EELSModelFit.fit(coreLoss: core.double, edges: c.edges, microscope: c.microscope, lowLoss: low?.double)
            #expect(a.pixels == b.pixels, Comment(rawValue: name))
            #expect(a.backgrounds == b.backgrounds, Comment(rawValue: name))
            #expect(a.residualSumOfSquares == b.residualSumOfSquares, Comment(rawValue: name))
            #expect(try EELSModelFit.fitBackground(coreLoss: core.single, edges: c.edges)
                    == EELSModelFit.fitBackground(coreLoss: core.double, edges: c.edges), Comment(rawValue: name))
            #expect(try EELSModelFit.twoAreaBackground(coreLoss: core.single, edges: c.edges)
                    == EELSModelFit.twoAreaBackground(coreLoss: core.double, edges: c.edges), Comment(rawValue: name))
            let s = try EELSMLLS.fit(coreLoss: core.single, edges: c.edges, microscope: c.microscope, lowLoss: low?.single)
            let d = try EELSMLLS.fit(coreLoss: core.double, edges: c.edges, microscope: c.microscope, lowLoss: low?.double)
            #expect(s.pixels == d.pixels, Comment(rawValue: name))
        }
    }

    /// The per-edge windowed fit (`edgeWindows`), on a planted C, N, O image of two pixels.
    @Test func testTheWindowedFitGivesTheSameFits() throws {
        let axis = EnergyAxis(offset: 150, scale: 0.75, size: 2020)
        let energies = (0 ..< axis.size).map { axis.energy(ofChannel: $0) }
        let edges: [CoreLossEdge] = [.carbonK, .nitrogenK, .oxygenK]
        var counts = [Double]()
        for (pixel, ratio) in [[3.0, 1.0, 2.0], [1.0, 2.0, 0.5]].enumerated() {
            var spectrum = energies.map { 3.0e10 * pow($0, -3) * Double(pixel + 1) }
            for (g, edge) in edges.enumerated() {
                let xs = EdgeCrossSection(table: GOSTable.bundled(edge)!, beamEnergyKeV: 200, convergenceMrad: 10, collectionMrad: 20)!
                let col = xs.differential(energies).map { $0 * axis.scale }
                for k in 0 ..< axis.size { spectrum[k] += ratio[g] * 1e6 * col[k] }
            }
            counts += spectrum
        }
        let core = try pair(try CoreLossSpectrumImage(nx: 2, ny: 1, channels: axis.size, counts: counts, energyOffsetEV: 150, energyScaleEV: 0.75))
        let c = EdgeWindows.defaults(onsetEV: 284, family: "K", precedingSignalEndEV: nil, nextOnsetEV: 401, axis: axis)
        let n = EdgeWindows.defaults(onsetEV: 401, family: "K", precedingSignalEndEV: c.signalEV.upperBound, nextOnsetEV: 532, axis: axis)
        let o = EdgeWindows.defaults(onsetEV: 532, family: "K", precedingSignalEndEV: n.signalEV.upperBound, nextOnsetEV: nil, axis: axis)
        let mic = EELSMicroscope(beamEnergyKeV: 200, convergenceMrad: 10, collectionMrad: 20)
        let a = try EELSModelFit.fit(coreLoss: core.single, edges: edges, microscope: mic, edgeWindows: [c, n, o])
        let b = try EELSModelFit.fit(coreLoss: core.double, edges: edges, microscope: mic, edgeWindows: [c, n, o])
        #expect(a.pixels == b.pixels)
        #expect(a.backgrounds == b.backgrounds)
    }

    @Test func testThicknessWithAZeroLossSpectrumImageGivesTheSameMap() throws {
        let image = try pair(try EELSFixtures.image("exspy-linescan"))
        let zeroLoss = try pair(try EELSFixtures.lineScanZeroLoss())
        let a = try EELSThickness.logRatio(image.single, zeroLoss: zeroLoss.single)
        let b = try EELSThickness.logRatio(image.double, zeroLoss: zeroLoss.double)
        #expect(a.pixels.map { $0.tOverLambda?.bitPattern } == b.pixels.map { $0.tOverLambda?.bitPattern })
        #expect(a.pixels.map { $0.i0?.bitPattern } == b.pixels.map { $0.i0?.bitPattern })
        #expect(a.pixels.contains { $0.i0 != nil })
    }
}
