//
//  EDSWindowPath.swift
//  Role: eXSpy's window path for one spectrum, public: the lines `add_lines()` picks for a list of elements at a beam energy, their
//        integration and background windows, the line intensities with and without background, the acquisition conditions with
//        eXSpy's defaults, and Cliff-Lorimer on the intensities with the single-line rule named by the caller. A front over the
//        kit's own functions (`XRayLines.defaultLines`, `WindowIntensity`, `BackgroundWindows`, `CliffLorimer`, `TakeOff`): no
//        second implementation. Merged from Electron Microscopy Studio's `SEMEDSLines` (at e7212f9, SpectroKit#131).
//
//  Ported from eXSpy 0.3.2 (the bundles) and 7185a4d1 (the kit's pins), the same functions in both: `signals/eds.py`
//  (`_get_lines_from_elements`, `_get_xray_lines_in_spectral_range`, `estimate_integration_windows`,
//  `estimate_background_windows`, `get_lines_intensity`, `get_take_off_angle`), `_misc/eds/utils.py` (`get_FWHM_at_Energy`,
//  `take_off_angle`, `quantification_cliff_lorimer`), `_misc/material.py` (`_weight_to_atomic`), and HyperSpy's
//  `UniformDataAxis.value2index` for the channel rule (`EnergyAxis`). Line width: Fiori and Newbury 1978; k-factor ratio: Cliff and
//  Lorimer 1975.
//
//  Numbers: float64, eXSpy's order of operations. A spectrum is a few thousand channels: plain loops on the CPU (D16). Parity: the
//  bundles m7/eds-lines-tm002 and m7/eds-quant-tm002 (BAM TM002), the line table and the channel rule exact, the rest rtol 1e-7
//  (SEMEDSParityTests).
//
//  DEVIATIONS from eXSpy (each the kit's, unchanged here): H and He get no line (`XRayLines.notRealLines`); a line whose width the
//  Fiori-Newbury law cannot give throws `SpectralLineError.widthUndefined` where eXSpy's `math.sqrt` raises; a NaN window edge throws
//  `EnergyAxisError.notANumber`.
//

import Foundation
import SpectroCore

/// Why the window path refuses a spectrum or its windows.
public nonisolated enum EDSWindowPathError: Error, Equatable, Sendable {
    /// The spectrum holds `counts` channels and its axis `axis`.
    case channelCount(counts: Int, axis: Int)
    /// One background window pair per integration window: `background` pairs for `integration` lines.
    case windowCount(background: Int, integration: Int)
    /// An integration window has two edges, a background window pair four; this one has `edges`.
    case windowEdges(expected: Int, edges: Int)
}

/// eXSpy's window path on one spectrum: lines, windows, intensities, conditions and Cliff-Lorimer.
///
/// Energies are keV, the axis is the spectrum's `EnergyAxis`. The functions are the kit's own (the line table of `XRayLines`, the
/// windows and channel rule of `WindowIntensity`), so a window or an intensity here is the one `ElementWindows` builds from the same
/// line and rules.
public nonisolated enum EDSWindowPath {
    /// What the acquisition says about the detector and the stage, as eXSpy reads `Acquisition_instrument.SEM`; a value the file does
    /// not state takes eXSpy's preference default (`_set_default_param`, `EDSConfig`). Unlike `EDXQuantification.takeOffAngle`,
    /// which takes every angle from its caller, these are eXSpy's defaults, stated here by name.
    public struct Conditions: Equatable, Sendable {
        public var beamEnergyKeV: Double?
        /// `Detector.EDS.energy_resolution_MnKa`, eV (default 130).
        public var energyResolutionMnKaEV: Double
        /// `Detector.EDS.elevation_angle`, degrees (default 35).
        public var elevationAngleDeg: Double
        /// `Detector.EDS.azimuth_angle`, degrees (default 0).
        public var azimuthAngleDeg: Double
        /// `Stage.tilt_alpha`, degrees (default 0).
        public var tiltAlphaDeg: Double
        /// `Stage.tilt_beta`, degrees (default 0, as `get_take_off_angle` reads it).
        public var tiltBetaDeg: Double
        public var liveTimeS: Double?
        public var realTimeS: Double?

        /// eXSpy's defaults (`EDSConfig`): 130 eV at Mn K alpha, 35 degrees elevation, 0 degrees azimuth and stage tilt.
        public static let defaultResolutionMnKaEV = 130.0
        public static let defaultElevationDeg = 35.0
        public static let defaultAzimuthDeg = 0.0
        public static let defaultTiltDeg = 0.0

        /// A nil value takes eXSpy's default.
        public init(beamEnergyKeV: Double?, energyResolutionMnKaEV: Double? = nil, elevationAngleDeg: Double? = nil,
                    azimuthAngleDeg: Double? = nil, tiltAlphaDeg: Double? = nil, tiltBetaDeg: Double? = nil,
                    liveTimeS: Double? = nil, realTimeS: Double? = nil) {
            self.beamEnergyKeV = beamEnergyKeV
            self.energyResolutionMnKaEV = energyResolutionMnKaEV ?? Self.defaultResolutionMnKaEV
            self.elevationAngleDeg = elevationAngleDeg ?? Self.defaultElevationDeg
            self.azimuthAngleDeg = azimuthAngleDeg ?? Self.defaultAzimuthDeg
            self.tiltAlphaDeg = tiltAlphaDeg ?? Self.defaultTiltDeg
            self.tiltBetaDeg = tiltBetaDeg ?? Self.defaultTiltDeg
            self.liveTimeS = liveTimeS
            self.realTimeS = realTimeS
        }

        /// `get_take_off_angle()`: the take-off angle from the stage tilts and the detector's azimuth and elevation, degrees
        /// (`EDXQuantification.takeOffAngle`).
        public var takeOffAngleDeg: Double {
            TakeOff.angle(tiltAlpha: tiltAlphaDeg, azimuth: azimuthAngleDeg, elevation: elevationAngleDeg, tiltBeta: tiltBetaDeg)
        }
    }

    // MARK: - Lines

    /// `_get_lines_from_elements(elements, only_one, only_lines)`, which `add_lines()` runs for `Sample.elements`: per element its
    /// lines of `onlyLines` (default `["a"]`: Ka, La, Ma) inside the spectral range; with `onlyOne` the first of them, in Python's
    /// string order, below half the beam energy, else the last; all lines sorted at the end. Without a beam energy the axis' last
    /// channel stands in, as eXSpy's fallback. An element without lines, or a symbol the table does not know, adds none.
    public static func lines(forElements elements: [String], axis: EnergyAxis, beamEnergyKeV: Double?, onlyOne: Bool = true,
                             onlyLines: [String]? = ["a"]) -> [String] {
        XRayLines.defaultLines(elements: elements, axis: axis, beamEnergy: beamEnergyKeV, onlyOne: onlyOne, onlyLines: onlyLines)
    }

    /// `_get_xray_lines_in_spectral_range`: the lines strictly between the axis' first channel and the smaller of its last channel and
    /// the beam energy; the rest (eXSpy warns about them), each list in the order given. A line the table lacks is out of range.
    public static func linesInSpectralRange(_ lines: [String], axis: EnergyAxis,
                                            beamEnergyKeV: Double?) -> (inRange: [String], outOfRange: [String]) {
        var inside: [String] = []
        var outside: [String] = []
        for line in lines {
            if XRayLines.linesInRange([line], axis: axis, beamEnergy: beamEnergyKeV).isEmpty { outside.append(line) } else { inside.append(line) }
        }
        return (inside, outside)
    }

    // MARK: - Windows

    /// `estimate_integration_windows(windows_width)`: per line `[E - width FWHM / 2, E + width FWHM / 2]`, keV, the FWHM from the
    /// detector's resolution at Mn K alpha (`XRayLines.fwhm`).
    public static func integrationWindows(_ lines: [String], width: Double = 2.0, resolutionMnKaEV: Double) throws -> [[Double]] {
        WindowIntensity.integrationWindows(lines: try spectralLines(lines, resolutionMnKaEV: resolutionMnKaEV), windowsWidth: width)
    }

    /// `estimate_background_windows(line_width, windows_width)`: per line a window `lineWidth.left` FWHM left of the line and one
    /// `lineWidth.right` FWHM right of it, each `windowsWidth` FWHM wide, `[l0, l1, r0, r1]` keV; then eXSpy's merge of overlapping
    /// windows (`BackgroundWindows.merge`). eXSpy's defaults are (2, 2) and 1; its EDS tutorial uses (5, 2).
    public static func backgroundWindows(_ lines: [String], lineWidth: (left: Double, right: Double) = (2, 2), windowsWidth: Double = 1,
                                         resolutionMnKaEV: Double) throws -> [[Double]] {
        BackgroundWindows.estimate(lines: try spectralLines(lines, resolutionMnKaEV: resolutionMnKaEV),
                                   lineWidth: [lineWidth.left, lineWidth.right], windowsWidth: windowsWidth)
    }

    private static func spectralLines(_ lines: [String], resolutionMnKaEV: Double) throws -> [SpectralLine] {
        try lines.map { try SpectralLine(id: $0, resolutionMnKaEV: resolutionMnKaEV) }
    }

    // MARK: - Intensities

    /// `get_lines_intensity(integration_windows=…, background_windows=…)` on one spectrum: the counts summed over each integration
    /// window (the EDS axis is binned, so `integrate1D` is a sum) and, with background windows, less the sum of the two background
    /// windows times the integration window's channel count over theirs. The channels are `WindowIntensity.resolve`'s (HyperSpy's
    /// `value2index`: half a channel goes towards zero at a non-negative energy); an edge off the axis throws where eXSpy raises
    /// (`WindowIntensityError`). A window narrower than one channel sums no channel, as eXSpy's slice does.
    public static func intensities(counts: [Double], axis: EnergyAxis, integrationWindows: [[Double]],
                                   backgroundWindows: [[Double]]? = nil) throws -> [Double] {
        guard counts.count == axis.size else { throw EDSWindowPathError.channelCount(counts: counts.count, axis: axis.size) }
        if let bw = backgroundWindows, bw.count != integrationWindows.count {
            throw EDSWindowPathError.windowCount(background: bw.count, integration: integrationWindows.count)
        }
        func sum(_ range: Range<Int>) -> Double {
            var s = 0.0
            for c in range { s += counts[c] }
            return s
        }
        var out: [Double] = []
        out.reserveCapacity(integrationWindows.count)
        for (i, window) in integrationWindows.enumerated() {
            guard window.count == 2 else { throw EDSWindowPathError.windowEdges(expected: 2, edges: window.count) }
            let background = backgroundWindows?[i]
            if let bw = background, bw.count != 4 { throw EDSWindowPathError.windowEdges(expected: 4, edges: bw.count) }
            let resolved = try WindowIntensity.resolve(integration: window, background: background, axis: axis)
            let left = resolved.background.map { sum($0.left) } ?? 0
            let right = resolved.background.map { sum($0.right) } ?? 0
            out.append(resolved.net(sum(resolved.signal), left, right))
        }
        return out
    }

    // MARK: - Cliff-Lorimer

    /// `quantification_cliff_lorimer(intensities, kfactors)` on one spectrum, as weight fractions (sum 1), with no absorption
    /// correction. The two routes differ only where fewer than two lines lie above eXSpy's 0.1-count floor, and the caller names
    /// the route:
    /// - `.refuse`, the kit's route and its default everywhere else (`EDXQuantification.quantify`): such a spectrum throws
    ///   `QuantError.needTwoLines`.
    /// - `.exspyHundredPercent`, eXSpy's raw arithmetic: the one line above the floor is 1 and the others 0 (all 0 with none).
    /// With two or more lines above the floor both routes run the same arithmetic. Times 100 is eXSpy's `quantification(method="CL")`;
    /// `EDXQuantification.weightToAtomic` gives the atomic percent.
    public static func cliffLorimer(intensities: [Double], kFactors: [Double], singleLine: SingleLineRule) throws -> [Double] {
        try CliffLorimer.weightFractions(intensities: intensities, kFactors: kFactors, singleLine: singleLine)
    }
}
