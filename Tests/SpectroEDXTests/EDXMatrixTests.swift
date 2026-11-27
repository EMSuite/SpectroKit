import Foundation
import Testing
import SpectroTestSupport
import SpectroCore
import SpectroEDX

/// The EDX reference matrix (`tools/edx-matrix/pins.py`): exSpy 0.3.2's window route from spectrum to composition on
/// 288 planted spectra (12 materials, thickness 10, 50 and 150 nm, two detector set-ups, 80 and 300 kV, two doses).
/// The kit runs the same chain end to end on the same counts: the take-off angle, exSpy's lines' FWHM, the integration and
/// background windows (with the overlap merge), the net counts, Cliff-Lorimer thin and with the absorption iteration, weight
/// and atomic percent, the mass thickness.
///
/// Tolerances are the ones exSpy's own tests use for each function (never stricter, D17; `docs/parity.md` rows):
/// windows 1e-3 keV absolute (`test_background_substraction`), net counts 0.1 absolute (`test_intensity`), percent 1e-3
/// absolute (`test_quant_lorimer*`), take-off and mass thickness 1e-7 relative (allclose defaults). The iteration count of
/// the absorption loop is an index and must be equal. The loop runs exSpy's signed stopping rule (`.exspySigned`), the one
/// these pins record.
///
/// The planted composition is in the fixture for the record. The window route's own error against it reaches 11 wt% where a
/// line of another family falls in a window (Ti L beside O K in SrTiO3, Ti K-beta on V K-alpha); the kit is held to exSpy here.
struct EDXMatrixTests {

    struct Case: Decodable {
        let name: String
        let material: String
        let resolutionMnKaEV: Double
        let tiltDegrees: Double
        let elevationDegrees: Double
        let takeOffDegrees: Double
        let thicknessNm: Double
        let lines: [String]
        let elements: [String]
        let lineEnergiesKeV: [Double]
        let countsOffset: Int
        let integrationWindows: [[Double]]
        let backgroundWindows: [[Double]]
        let intensities: [Double]
        let kFactors: [Double]
        let weightPercent: [Double]
        let atomicPercent: [Double]
        let massThickness: Double
        let iterations: Int
        let thinWeightPercent: [Double]
        let plantedWeightPercent: [Double]
    }
    struct Axis: Decodable { let size: Int; let offset: Double; let scale: Double }
    struct Pins: Decodable { let exspy: String; let axis: Axis; let countsFile: String; let cases: [Case] }

    private static let pins: Pins = {
        let data = try! Data(contentsOf: EDXFixtures.url("edx-matrix.json"))
        return try! JSONDecoder().decode(Pins.self, from: data)
    }()
    private static let counts: [UInt32] = {
        let data = try! Data(contentsOf: EDXFixtures.url(pins.countsFile))
        return data.withUnsafeBytes { raw in
            (0..<(data.count / 4)).map { UInt32(littleEndian: raw.loadUnaligned(fromByteOffset: $0 * 4, as: UInt32.self)) }
        }
    }()
    private static let table: MassAbsorptionTable = {
        let url = EDXFixtures.repoRoot.appendingPathComponent("Sources/SpectroEDX/Resources/Spectroscopy/FFastMAC.csv")
        return try! MassAbsorptionTable(contentsOf: url)
    }()

    private static let rtol = 1e-7, windowKeV = 1e-3, countsAtol = 0.1, percent = 1e-3

    private func close(_ a: Double, _ b: Double) -> Bool {
        abs(a - b) <= Self.rtol * abs(b)
    }
    private func expectClose(_ a: [Double], _ b: [Double], atol: Double, _ what: String, _ c: Case,
                             sourceLocation: SourceLocation = #_sourceLocation) {
        #expect(a.count == b.count, "\(c.name) \(what): \(a.count) vs \(b.count) values", sourceLocation: sourceLocation)
        for (i, (x, y)) in zip(a, b).enumerated() {
            #expect(abs(x - y) <= atol, "\(c.name) \(what)[\(i)] \(x) vs exSpy \(y)", sourceLocation: sourceLocation)
        }
    }

    /// The whole chain on one case; every difference is recorded as a failure.
    private func check(_ c: Case) throws {
        let axis = EnergyAxis(offset: Self.pins.axis.offset, scale: Self.pins.axis.scale, size: Self.pins.axis.size)
        let spectrum = Self.counts[c.countsOffset..<(c.countsOffset + axis.size)].map(UInt64.init)

        let toa = EDXQuantification.takeOffAngle(tiltAlpha: c.tiltDegrees, azimuth: 0, elevation: c.elevationDegrees)
        #expect(close(toa, c.takeOffDegrees), "\(c.name) take-off \(toa) vs \(c.takeOffDegrees)")

        let lines = try c.lines.map { try SpectralLine(id: $0, resolutionMnKaEV: c.resolutionMnKaEV) }
        #expect(lines.map(\.energy) == c.lineEnergiesKeV, "\(c.name) line energies")
        let iw = WindowIntensity.integrationWindows(lines: lines)
        let bw = BackgroundWindows.estimate(lines: lines)
        for i in lines.indices {
            expectClose(iw[i], c.integrationWindows[i], atol: Self.windowKeV, "integration window \(c.lines[i])", c)
            expectClose(bw[i], c.backgroundWindows[i], atol: Self.windowKeV, "background window \(c.lines[i])", c)
        }
        let resolved = try lines.indices.map { try WindowIntensity.resolve(integration: iw[$0], background: bw[$0], axis: axis) }
        let net = WindowIntensity.netCounts(spectrum: spectrum, windows: resolved)
        expectClose(net, c.intensities, atol: Self.countsAtol, "net counts", c)

        let thin = try CliffLorimer.weightFractions(intensities: net, kFactors: c.kFactors).map { $0 * 100 }
        expectClose(thin, c.thinWeightPercent, atol: Self.percent, "thin wt%", c)

        let setup = AbsorptionSetup(elements: c.elements, lineEnergiesKeV: c.lineEnergiesKeV,
                                    geometry: SingleDetectorGeometry(takeOffDegrees: toa), table: Self.table)
        let res = try AbsorptionCorrection.cliffLorimer(intensities: net, kFactors: c.kFactors, setup: setup,
                                                        thicknessNm: c.thicknessNm, options: .init(criterion: .exspySigned))
        #expect(res.iterations == c.iterations, "\(c.name) iterations \(res.iterations) vs exSpy \(c.iterations)")
        expectClose(res.composition, c.weightPercent, atol: Self.percent, "wt%", c)
        let atomic = try CompositionConversion.weightToAtomic(res.composition, elements: c.elements)
        expectClose(atomic, c.atomicPercent, atol: Self.percent, "at%", c)
        #expect(close(res.massThickness, c.massThickness), "\(c.name) mass thickness \(res.massThickness) vs \(c.massThickness)")
    }

    @Test func fixtureCoversTheMatrix() {
        let cases = Self.pins.cases
        #expect(Self.pins.exspy == "0.3.2")
        #expect(cases.count == 288)
        #expect(Set(cases.map(\.material)).count == 12)
        #expect(Self.counts.count == cases.count * Self.pins.axis.size)
        #expect(cases.contains { $0.lines.contains { $0.hasSuffix("_La") } }, "an L line")
        #expect(cases.contains { c in zip(c.backgroundWindows, c.backgroundWindows.dropFirst()).contains { $0 == $1 } }, "merged windows")
    }

    /// One case per material, the high dose at set-up A, 300 kV, 50 nm.
    @Test func matrixOnePerMaterial() throws {
        let picked = Self.pins.cases.filter { $0.name.hasSuffix(" A 300kV 50nm high") }
        #expect(picked.count == 12)
        for c in picked { try check(c) }
    }

    @Test func matrixAll() throws {
        for c in Self.pins.cases { try check(c) }
    }
}
