//
//  SEMEDSParityTests.swift
//  exspy's window path (`EDSWindowPath`, eXSpy 0.3.2's `get_lines_intensity`, `estimate_integration_windows`,
//  `estimate_background_windows`, `quantification_cliff_lorimer` and `weight_to_atomic` on the kit's own line table, windows and channel
//  rule) against the bundles of the app's tools/parity/sem on the BAM TM002 spectrum: `m7/eds-lines-tm002` (figure F16: the lines exspy
//  picks, the integration and background windows, the raw and background-subtracted intensities), `m7/eds-quant-tm002` (figure F17:
//  Cliff-Lorimer weight and atomic percent with unit and planted k-factors, the take-off angle, each line's FWHM, the default background
//  windows, a 2.1-FWHM window, the lines at 5, 20 and 30 kV) with exspy's element table (its `db_` arrays, exact). Tolerances are the
//  bundles' (the line table and the channel rule exact, float64 reductions rtol 1e-7, the app's docs/parity.md rows "SEM-EDS (exspy's
//  window path)").
//  Moved from Electron Microscopy Studio's EMStudioTests/SEMEDSParityTests.swift (at e7212f9) by SpectroKit#131 with the two bundles and
//  their input, byte for byte (Fixtures/sem); no number changed. The app read the spectrum through its HyperSpy door; here the hspy file
//  is read with FormatoHDF5's plain HDF5 access (counts, axis, `Sample.elements`, the SEM conditions), and the app's own table
//  (`XrayLineDatabase`) is checked through the bundle's `db_` arrays it was generated from.
//  Red before the merge: no `EDSWindowPath`. Red after a change to the channel rounding, the window merge or the Cliff-Lorimer chain:
//  the intensities move by whole channels' counts and the composition by percent.
//

import CryptoKit
import Foundation
import FormatoCore
import FormatoHDF5
import SpectroCore
import SpectroEDX
import XCTest

nonisolated final class SEMEDSParityTests: XCTestCase {

    /// The TM002 spectrum as the hspy file stores it.
    private struct Spectrum {
        let counts: [Double]
        let axis: EnergyAxis
        let units: String
        let elements: [String]
        let conditions: EDSWindowPath.Conditions
    }

    private struct Unreadable: Error, CustomStringConvertible {
        let description: String
    }

    private func tm002(_ ref: SEMParityReference) throws -> Spectrum {
        let url = try XCTUnwrap(ref.inputURL(), "the bundle names the TM002 fixture")
        let input = try XCTUnwrap(ref.manifest.inputs.first)
        let digest = try SHA256.hash(data: Data(contentsOf: url)).map { String(format: "%02x", $0) }.joined()
        XCTAssertEqual(digest, input.sha256, "the fixture is the bundle's input byte for byte")

        let reader = try HDF5Reader(url: url)
        let name = try XCTUnwrap(reader.children(of: "/Experiments").first, "the file holds one experiment")
        let experiment = "/Experiments/" + name
        func attribute(_ key: String, _ group: String) throws -> HDF5AttributeValue? {
            let path = experiment + group
            guard try reader.attributeNames(of: path).contains(key) else { return nil }
            return try reader.attribute(key, of: path)
        }
        func number(_ key: String, _ group: String) throws -> Double? {
            switch try attribute(key, group) {
            case .real(let v)?: return v
            case .integer(let v)?: return Double(v)
            default: return nil
            }
        }
        guard let offset = try number("offset", "/axis-0"), let scale = try number("scale", "/axis-0") else {
            throw Unreadable(description: "the axis has no offset or scale")
        }
        guard case .text(let units)? = try attribute("units", "/axis-0") else { throw Unreadable(description: "the axis has no units") }
        let counts = try reader.readDataset(experiment + "/data").doubles()
        let elements = try reader.readStrings(experiment + "/metadata/Sample/_list_elements")
        let sem = "/metadata/Acquisition_instrument/SEM", eds = sem + "/Detector/EDS", stage = sem + "/Stage"
        let conditions = try EDSWindowPath.Conditions(
            beamEnergyKeV: number("beam_energy", sem), energyResolutionMnKaEV: number("energy_resolution_MnKa", eds),
            elevationAngleDeg: number("elevation_angle", eds), azimuthAngleDeg: number("azimuth_angle", eds),
            tiltAlphaDeg: number("tilt_alpha", stage), tiltBetaDeg: number("tilt_beta", stage),
            liveTimeS: number("live_time", eds), realTimeS: number("real_time", eds))
        return Spectrum(counts: counts, axis: EnergyAxis(offset: offset, scale: scale, size: counts.count), units: units,
                        elements: elements, conditions: conditions)
    }

    private static func element(_ line: String) -> String { XRayLines.line(line)?.element ?? "" }

    // MARK: - One line table

    /// exspy's element table (`exspy._misc.elements.elements`, the app's generated `XrayLineDatabase`, held exact by its `db_` arrays)
    /// and the kit's `XRayLines` (from `xray_lines.json`) with `QuantElementData`: every line the app's table holds is in the kit's
    /// at the same energy and weight, with no line more; every element has the same Z and atomic weight and the same lines.
    func testTheAppsLineTableIsTheKits() throws {
        let ref = try SEMParityReference.load("m7/eds-quant-tm002")
        let elements = try XCTUnwrap(ref.values["db_elements"] as? [String])
        let lines = try XCTUnwrap(ref.values["db_xray_lines"] as? [String])
        let z = try ref.array("db_z").values, atomicWeight = try ref.array("db_atomic_weight").values
        let energy = try ref.array("db_line_energy_keV").values, weight = try ref.array("db_line_weight").values
        XCTAssertEqual(elements.count, 95)
        XCTAssertEqual(lines.count, 947)
        XCTAssertEqual(energy.count, lines.count)
        XCTAssertEqual(weight.count, lines.count)
        XCTAssertEqual(XRayLines.all.count, lines.count, "the kit's table has no line the app's lacks")
        XCTAssertEqual(Set(XRayLines.all.map(\.id)), Set(lines))

        var zOf: [String: Double] = [:]
        for (i, symbol) in elements.enumerated() {
            zOf[symbol] = z[i]
            let entry = try XCTUnwrap(QuantElementData.entry(symbol), symbol)
            XCTAssertEqual(Double(entry.z), z[i], symbol)
            XCTAssertEqual(entry.atomicWeight, atomicWeight[i], symbol)
            let own = Set(lines.filter { $0.hasPrefix(symbol + "_") })
            XCTAssertEqual(Set(XRayLines.all.filter { $0.element == symbol }.map(\.id)), own, symbol)
        }
        for (i, id) in lines.enumerated() {
            let line = try XCTUnwrap(XRayLines.line(id), id)
            XCTAssertEqual(line.energy, energy[i], id)
            XCTAssertEqual(line.weight, weight[i], id)
            XCTAssertEqual(Double(line.atomicNumber), zOf[line.element], id)
        }
        XCTAssertEqual(XRayLines.line("Mn_La")?.energy, 0.63316)
        XCTAssertNil(XRayLines.line("Mn_Xa"))
        XCTAssertEqual(XRayLines.lines(of: "Li"), [], "exspy lists Li without X-ray lines")
    }

    // MARK: - The conditions the file states

    func testTM002ConditionsAreTheFilesAndTheTakeOffAngleIsExspys() throws {
        let ref = try SEMParityReference.load("m7/eds-quant-tm002")
        let spectrum = try tm002(ref)
        XCTAssertEqual(spectrum.units, "keV")
        XCTAssertEqual(spectrum.elements, ["Al", "C", "Cu", "Mn", "Zr"])
        let c = spectrum.conditions
        XCTAssertEqual(c.beamEnergyKeV, 10.0)
        XCTAssertEqual(c.energyResolutionMnKaEV, 130.0)
        XCTAssertEqual(c.elevationAngleDeg, 37.0)
        XCTAssertEqual(c.azimuthAngleDeg, 0.0)
        XCTAssertEqual(c.tiltAlphaDeg, 0.0)
        XCTAssertEqual(c.tiltBetaDeg, 0.0)
        XCTAssertEqual(try XCTUnwrap(c.liveTimeS), 19.997292, accuracy: 1e-12)
        SEMAssertParity(c.takeOffAngleDeg, ref, scalar: "take_off_angle_deg")
        // exspy's docstring: 20° of stage tilt adds 20° to the take-off angle of this geometry.
        XCTAssertEqual(EDXQuantification.takeOffAngle(tiltAlpha: 20, azimuth: 0, elevation: 37), 57.0, accuracy: 1e-9)
    }

    /// What the file leaves out takes exspy's defaults (`EDSConfig`).
    func testConditionsTakeExspysDefaults() {
        let bare = EDSWindowPath.Conditions(beamEnergyKeV: nil)
        XCTAssertNil(bare.beamEnergyKeV)
        XCTAssertEqual(bare.energyResolutionMnKaEV, 130)
        XCTAssertEqual(bare.elevationAngleDeg, 35)
        XCTAssertEqual(bare.azimuthAngleDeg, 0)
        XCTAssertEqual(bare.tiltAlphaDeg, 0)
        XCTAssertEqual(bare.tiltBetaDeg, 0)
        XCTAssertNil(bare.liveTimeS)
        XCTAssertNil(bare.realTimeS)
        XCTAssertEqual(bare.takeOffAngleDeg, EDXQuantification.takeOffAngle(tiltAlpha: 0, azimuth: 0, elevation: 35))
    }

    // MARK: - Figure F16

    func testFigureF16LinesWindowsAndIntensitiesMatchExspy() throws {
        let ref = try SEMParityReference.load("m7/eds-lines-tm002")
        let spectrum = try tm002(ref)
        let resolution = spectrum.conditions.energyResolutionMnKaEV
        let lines = EDSWindowPath.lines(forElements: spectrum.elements, axis: spectrum.axis, beamEnergyKeV: spectrum.conditions.beamEnergyKeV)
        XCTAssertEqual(lines, ref.values["xray_lines"] as? [String])
        let energies = try XCTUnwrap(ref.values["line_energies_keV"] as? [String: Double])
        for line in lines { XCTAssertEqual(XRayLines.line(line)?.energy, energies[line], line) }

        let integration = try EDSWindowPath.integrationWindows(lines, width: 2.0, resolutionMnKaEV: resolution)
        try SEMAssertParity(integration.flatMap { $0 }, ref.array("integration_windows_keV"), shape: [lines.count, 2])
        let background = try EDSWindowPath.backgroundWindows(lines, lineWidth: (5.0, 2.0), resolutionMnKaEV: resolution)
        try SEMAssertParity(background.flatMap { $0 }, ref.array("background_windows_keV"), shape: [lines.count, 4])

        let raw = try EDSWindowPath.intensities(counts: spectrum.counts, axis: spectrum.axis, integrationWindows: integration)
        try SEMAssertParity(raw, ref.array("intensity_raw"))
        let net = try EDSWindowPath.intensities(counts: spectrum.counts, axis: spectrum.axis, integrationWindows: integration,
                                                backgroundWindows: background)
        try SEMAssertParity(net, ref.array("intensity_background_subtracted"))

        // One implementation: the kit's integer window sums (what `ElementWindows` reads) give the same numbers on these counts.
        let integerCounts = spectrum.counts.map { UInt64($0) }
        let resolved = try zip(integration, background).map { window, bw in
            try WindowIntensity.resolve(integration: window, background: bw, axis: spectrum.axis)
        }
        XCTAssertEqual(WindowIntensity.netCounts(spectrum: integerCounts, windows: resolved), net)
    }

    // MARK: - Figure F17 and the rest of the window path

    /// F17 pins exspy's raw arithmetic (`.exspyHundredPercent`); on TM002 every line is far above the floor, so the kit's refusing
    /// route (`.refuse`, and `EDXQuantification.quantify`, its default) gives the same numbers.
    func testFigureF17CliffLorimerMatchesExspy() throws {
        let ref = try SEMParityReference.load("m7/eds-quant-tm002")
        let spectrum = try tm002(ref)
        let resolution = spectrum.conditions.energyResolutionMnKaEV
        let lines = EDSWindowPath.lines(forElements: spectrum.elements, axis: spectrum.axis, beamEnergyKeV: spectrum.conditions.beamEnergyKeV)
        XCTAssertEqual(lines, ref.values["xray_lines"] as? [String])
        var fwhm: [Double] = []
        for line in lines {
            let energy = try XCTUnwrap(XRayLines.line(line)?.energy, line)
            try fwhm.append(XCTUnwrap(XRayLines.fwhm(resolutionMnKaEV: resolution, atEnergy: energy), line))
        }
        try SEMAssertParity(fwhm, ref.array("line_fwhm_keV"))

        // The intensities the EDS room quantifies: background-subtracted with exspy's tutorial windows (5, 2).
        let net = try EDSWindowPath.intensities(
            counts: spectrum.counts, axis: spectrum.axis,
            integrationWindows: EDSWindowPath.integrationWindows(lines, resolutionMnKaEV: resolution),
            backgroundWindows: EDSWindowPath.backgroundWindows(lines, lineWidth: (5.0, 2.0), resolutionMnKaEV: resolution))
        let elements = lines.map { Self.element($0) }
        let planted = try XCTUnwrap(ref.values["kfactors_planted"] as? [Double])
        for (name, k) in [("unit", [Double](repeating: 1, count: lines.count)), ("planted", planted)] {
            let expectedWeight = try ref.array("weight_percent_\(name)"), expectedAtomic = try ref.array("atomic_percent_\(name)")
            for rule in [SingleLineRule.exspyHundredPercent, .refuse] {
                let weight = try EDSWindowPath.cliffLorimer(intensities: net, kFactors: k, singleLine: rule).map { $0 * 100.0 }
                SEMAssertParity(weight, expectedWeight)
                try SEMAssertParity(EDXQuantification.weightToAtomic(weight, elements: elements), expectedAtomic)
                XCTAssertEqual(weight.reduce(0, +), 100, accuracy: 1e-9, "\(name) \(rule)")
            }
            let setup = EDXQuantificationSetup(elements: elements, method: .cliffLorimer(kFactors: k, thicknessNm: nil))
            let pixel = try EDXQuantification.quantify(intensities: net, setup: setup)
            SEMAssertParity(pixel.weightPercent, expectedWeight)
            SEMAssertParity(pixel.atomicPercent, expectedAtomic)
        }
    }

    func testDefaultBackgroundWindowsAndAWiderWindowMatchExspy() throws {
        let ref = try SEMParityReference.load("m7/eds-quant-tm002")
        let spectrum = try tm002(ref)
        let resolution = spectrum.conditions.energyResolutionMnKaEV
        let lines = try XCTUnwrap(ref.values["xray_lines"] as? [String])
        let integration = try EDSWindowPath.integrationWindows(lines, resolutionMnKaEV: resolution)
        let background = try EDSWindowPath.backgroundWindows(lines, resolutionMnKaEV: resolution)
        try SEMAssertParity(background.flatMap { $0 }, ref.array("background_windows_default_keV"), shape: [lines.count, 4])
        let net = try EDSWindowPath.intensities(counts: spectrum.counts, axis: spectrum.axis, integrationWindows: integration,
                                                backgroundWindows: background)
        try SEMAssertParity(net, ref.array("intensity_background_default"))
        let wider = try EDSWindowPath.integrationWindows(lines, width: 2.1, resolutionMnKaEV: resolution)
        let raw = try EDSWindowPath.intensities(counts: spectrum.counts, axis: spectrum.axis, integrationWindows: wider)
        try SEMAssertParity(raw, ref.array("intensity_raw_width_2_1"))
    }

    func testLinesFollowTheBeamEnergyAsExspyPicksThem() throws {
        let ref = try SEMParityReference.load("m7/eds-quant-tm002")
        let spectrum = try tm002(ref)
        let byBeam = try XCTUnwrap(ref.values["xray_lines_by_beam_energy"] as? [String: [String]])
        XCTAssertEqual(byBeam.count, 3)
        for (kv, expected) in byBeam {
            let beam = try XCTUnwrap(Double(kv))
            XCTAssertEqual(EDSWindowPath.lines(forElements: spectrum.elements, axis: spectrum.axis, beamEnergyKeV: beam), expected, "\(kv) kV")
        }
    }

    // MARK: - Rules the bundles do not reach

    /// HyperSpy's channel rule, the one the kit's windows use (`EnergyAxis`): half a channel rounds towards zero for a non-negative
    /// energy, away from zero for a negative one; off the axis there is no channel.
    func testChannelRoundingIsHyperSpys() {
        let axis = EnergyAxis(offset: -0.1, scale: 0.01, size: 1024)
        XCTAssertEqual(try axis.index(of: -0.1), 0)
        XCTAssertEqual(try axis.index(of: 0.0), 10)
        XCTAssertEqual(try axis.index(of: 0.105), 20, "20.5 rounds towards zero")
        XCTAssertEqual(try axis.index(of: 0.106), 21)
        XCTAssertEqual(try axis.index(of: -0.045), 6, "5.5 at a negative energy rounds away from zero")
        XCTAssertThrowsError(try axis.index(of: -0.2))
        XCTAssertThrowsError(try axis.index(of: 10.5))
        XCTAssertEqual(try axis.channelRange(from: -0.5, to: 0.0), 0 ..< 10, "a start below the axis starts at 0")
        XCTAssertEqual(try axis.channelRange(from: 10.0, to: 20.0), 1010 ..< 1024, "a stop above the axis runs to the end")
        XCTAssertThrowsError(try axis.channelRange(from: 11.0, to: 12.0))
    }

    /// exspy's floor: a line at or below 0.1 counts never serves as the reference. One line above it is the whole composition on
    /// exspy's route and a refusal on the kit's.
    func testCliffLorimerFloorAndSingleLine() throws {
        XCTAssertEqual(try EDSWindowPath.cliffLorimer(intensities: [0.05, 200, 0], kFactors: [1, 1, 1], singleLine: .exspyHundredPercent), [0, 1, 0])
        XCTAssertEqual(try EDSWindowPath.cliffLorimer(intensities: [0, 0], kFactors: [1, 1], singleLine: .exspyHundredPercent), [0, 0])
        XCTAssertThrowsError(try EDSWindowPath.cliffLorimer(intensities: [0.05, 200, 0], kFactors: [1, 1, 1], singleLine: .refuse)) { error in
            XCTAssertEqual(error as? QuantError, .needTwoLines(linesAbove: 1))
        }
        for rule in [SingleLineRule.exspyHundredPercent, .refuse] {
            let two = try EDSWindowPath.cliffLorimer(intensities: [100, 300], kFactors: [1, 1], singleLine: rule)
            XCTAssertEqual(two[0], 0.25, accuracy: 1e-15)
            XCTAssertEqual(two[1], 0.75, accuracy: 1e-15)
            XCTAssertThrowsError(try EDSWindowPath.cliffLorimer(intensities: [1, 2], kFactors: [1], singleLine: rule))
        }
        XCTAssertEqual(try EDXQuantification.weightToAtomic([0, 0], elements: ["Cu", "Sn"]), [0, 0])
    }

    /// exspy's own docstring examples on this spectrum (`get_lines_intensity`): Mn La with the default window, and Mn Ka with
    /// background windows when Mn is the only element.
    func testExspyDocstringIntensities() throws {
        let ref = try SEMParityReference.load("m7/eds-lines-tm002")
        let spectrum = try tm002(ref)
        let resolution = spectrum.conditions.energyResolutionMnKaEV
        let la = try EDSWindowPath.intensities(counts: spectrum.counts, axis: spectrum.axis,
                                               integrationWindows: EDSWindowPath.integrationWindows(["Mn_La"], resolutionMnKaEV: resolution))
        XCTAssertEqual(la, [96700])
        let ka = ["Mn_Ka"]
        let ka21 = try EDSWindowPath.intensities(counts: spectrum.counts, axis: spectrum.axis,
                                                 integrationWindows: EDSWindowPath.integrationWindows(ka, width: 2.1, resolutionMnKaEV: resolution))
        XCTAssertEqual(ka21, [53597])
        let net = try EDSWindowPath.intensities(counts: spectrum.counts, axis: spectrum.axis,
                                                integrationWindows: EDSWindowPath.integrationWindows(ka, resolutionMnKaEV: resolution),
                                                backgroundWindows: EDSWindowPath.backgroundWindows(ka, resolutionMnKaEV: resolution))
        XCTAssertEqual(net, [46716])
    }

    /// The lines outside the spectral range are the rest, in order; a spectrum and its axis must agree, and a window off the axis is
    /// refused where HyperSpy raises.
    func testRangeAndRefusals() {
        let axis = EnergyAxis(offset: 0, scale: 0.01, size: 1024)
        let split = EDSWindowPath.linesInSpectralRange(["Mn_La", "Mn_Ka", "Xx_Ka", "Al_Ka"], axis: axis, beamEnergyKeV: 5)
        XCTAssertEqual(split.inRange, ["Mn_La", "Al_Ka"])
        XCTAssertEqual(split.outOfRange, ["Mn_Ka", "Xx_Ka"])
        XCTAssertThrowsError(try EDSWindowPath.intensities(counts: [1, 2], axis: EnergyAxis(offset: 0, scale: 1, size: 3),
                                                           integrationWindows: [[0, 1]])) { error in
            XCTAssertEqual(error as? EDSWindowPathError, .channelCount(counts: 2, axis: 3))
        }
        let small = EnergyAxis(offset: 0, scale: 0.01, size: 100)
        let counts = [Double](repeating: 1, count: 100)
        XCTAssertThrowsError(try EDSWindowPath.intensities(counts: counts, axis: small, integrationWindows: [[1.5, 2]])) { error in
            XCTAssertEqual(error as? WindowIntensityError, .axis(.sliceBeyondAxis(start: 1.5, stop: 2)))
        }
        XCTAssertThrowsError(try EDSWindowPath.intensities(counts: counts, axis: small, integrationWindows: [[0.1, 0.2]],
                                                           backgroundWindows: [])) { error in
            XCTAssertEqual(error as? EDSWindowPathError, .windowCount(background: 0, integration: 1))
        }
        XCTAssertThrowsError(try EDSWindowPath.integrationWindows(["Mn_Xa"], resolutionMnKaEV: 130)) { error in
            XCTAssertEqual(error as? SpectralLineError, .unknownLine("Mn_Xa"))
        }
    }
}
