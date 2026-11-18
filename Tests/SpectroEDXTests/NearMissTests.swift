//
//  NearMissTests.swift
//  Role: SpectroKit#124, Auto ID's near misses. Paul's choice (2026-11-18): Auto ID suggests them and never adds them, so each case
//        lands in `ProposalResult.possible` (a consumer offers a one-click Add) and never in `proposed`. The three cases were measured on
//        SpectroKit 0.15.0 through the room's own Auto ID path (EMSuite/EMSuiteUI#143); the files stay out of the repo, so carbon and
//        gallium are planted here (noise-free: the fit's own line model plus a continuum) and the Mg case is the result tier on its
//        measured numbers.
//

import Foundation
import Testing
import SpectroCore
import SpectroEDX

struct NearMissTests {

    private static let axis = EnergyAxis(offset: 0, scale: 0.01, size: 1100)   // 0 to 10.99 keV
    private static let energies = (0..<axis.size).map { axis.energy(ofChannel: $0) }

    /// Kramers through a window that cuts the continuum below about 0.4 keV, as `ProposerPolicyTests`.
    private static let kramers: @Sendable (Double) -> Double = { e in
        let x = max(e, 0.05)
        return 60 * max(200 - e, 0) / x * exp(-pow(0.30 / x, 3)) * 0.02 + 0.5
    }

    /// Expected counts (not drawn): the continuum x dose plus the named line groups' areas x dose.
    private static func expected(_ areas: [String: Double], dose k: Double = 10) -> [Double] {
        var out = energies.map { $0 >= 0.2 ? k * kramers($0) : 0 }
        for (id, area) in areas {
            let el = String(id.prefix { $0 != "_" })
            let g = EDSLineModel.build(elements: [el], axis: axis, beamEnergy: 200, resolutionMnKaEV: 130).groups.first { $0.id == id }!
            let col = EDSLineModel.column(of: g.lines, energies: energies, scale: axis.scale)
            for i in out.indices { out[i] += k * area * col[i] }
        }
        return out
    }

    private static func settings(_ current: [String]) -> FitSettings {
        FitSettings.standard(elements: current, axis: axis, resolutionMnKaEV: 130, beamEnergy: 200)
    }

    private static func candidate(_ r: ProposalResult, _ el: String) -> ElementCandidate? { r.candidates.first { $0.element == el } }

    // MARK: Carbon below the 0.45 keV floor

    /// Measured: on the YAlO3 DM4 (Zenodo 10609594) C was refused before any test (its line at 0.28 keV lies below 0.45 keV); with the floor at
    /// 0.2 keV it read net 512.4 ± 58.3 against L_D 110.5 (4.6 L_D). Now C is measured once with the settled model and shown as possible, never
    /// proposed and never refused. Red before: C was refused and absent from the candidates.
    @Test func carbonBelowTheFloorIsPossibleNeverProposed() throws {
        let r = try ElementProposer().propose(counts: Self.expected(["Al_Ka": 9000, "C_Ka": 600]), axis: Self.axis, settings: Self.settings(["Al"]))
        let c = try #require(Self.candidate(r, "C"), "C is measured below the floor")
        #expect(c.group == "C_Ka" && c.aboveDetectionLimit, "C net \(c.net) L_D \(c.detectionLimit)")
        #expect(r.possible.contains(c) && !c.isProposed && !r.proposed.contains(c) && !r.sumPeakQuestions.contains(c))
        #expect(c.possibleReason?.hasPrefix("Its K\u{03B1} (0.28 keV) lies below 0.45 keV, the lowest line Auto ID adds") == true, "\(c.possibleReason ?? "nil")")
        #expect(!r.refused.contains { $0.element == "C" })
        // The unlisted-line check leaves it out, as before it was measured.
        #expect(c.measuredBelowFloor)
        #expect(r.notes.first?.contains("C K\u{03B1} is measured") == true, "\(r.notes.first ?? "")")
        // Never added: the reported Al numbers are those of the run that cannot see carbon.
        let without = try ElementProposer(pool: ElementProposer.defaultPool.filter { $0 != "C" })
            .propose(counts: Self.expected(["Al_Ka": 9000, "C_Ka": 600]), axis: Self.axis, settings: Self.settings(["Al"]))
        #expect(without.candidates.filter { $0.element != "C" } == r.candidates.filter { $0.element != "C" })
    }

    /// The floor still holds for everything else (N K-alpha stays refused) and for C when the fitted range cannot hold its line.
    @Test func carbonOutsideTheFittedRangeIsStillRefused() throws {
        var s = Self.settings(["Al"])
        s.fitFrom = 0.35
        let r = try ElementProposer().propose(counts: Self.expected(["Al_Ka": 9000, "C_Ka": 600]), axis: Self.axis, settings: s)
        #expect(Self.candidate(r, "C") == nil)
        #expect(r.refused.contains { $0.element == "C" } && r.refused.contains { $0.element == "N" })
        #expect(r.notes.first?.contains("C K\u{03B1}") == false)
    }

    // MARK: Gallium on a pile-up energy, with its L-alpha

    /// Measured: on Velox SI 1456 Ga Kα (9.252 keV, 545.9 L_D) was held as a sum-peak question (Cu+Ge at 9.236 keV, Ni+Si), though Ga Lα
    /// (1.098 keV, 107.5 L_D) was there too and no pile-up makes it. Planted: Ni and Si listed and strong, so Ni+Si (9.218 keV) sits 34 eV from
    /// Ga Kα. With its Lα, Ga is possible (the sum-peak conflict kept), never proposed; without it, Ga stays a sum-peak question.
    /// Red before: Ga was a sum-peak question in both.
    @Test func galliumWithItsLAlphaIsPossibleNotAPileUpQuestion() throws {
        let withL = try ElementProposer().propose(counts: Self.expected(["Ni_Ka": 20000, "Si_Ka": 20000, "Ga_Ka": 4000, "Ga_La": 1000]),
                                                  axis: Self.axis, settings: Self.settings(["Ni", "Si"]))
        let ga = try #require(Self.candidate(withL, "Ga"))
        #expect(ga.group == "Ga_Ka" && ga.hasSumPeakQuestion, "\(ga.group) \(ga.conflicts.map(\.question))")
        #expect(withL.possible.contains(ga) && !withL.proposed.contains(ga) && !withL.sumPeakQuestions.contains(ga))
        let reason = try #require(ga.possibleReason)
        #expect(reason.contains("L\u{03B1} at 1.10 keV is not a pile-up energy") && reason.contains("Ni+Si sum peak"), "\(reason)")

        let kOnly = try ElementProposer().propose(counts: Self.expected(["Ni_Ka": 20000, "Si_Ka": 20000, "Ga_Ka": 4000]),
                                                  axis: Self.axis, settings: Self.settings(["Ni", "Si"]))
        let held = try #require(Self.candidate(kOnly, "Ga"))
        #expect(kOnly.sumPeakQuestions.contains(held) && held.possibleReason == nil, "\(held.possibleReason ?? "nil")")
    }

    // MARK: Mg beside Al: no kit change

    /// Measured: on AlMgSi_4D_EDX.dm4 (whole map) Mg Kα read net 544.1 ± 60.0 against L_D 172.0 (3.16 L_D) and was held back only by the gain
    /// bar (81 against 100). The bar stays (Paul, 2026-11-18); the candidate is already in `possible`, which the room shows with Add.
    @Test func magnesiumHeldByTheGainBarIsPossible() {
        let sigmaZero = (172.0 - 2.71) / 3.29
        let mg = ElementCandidate(element: "Mg", group: "Mg_Ka", energyKeV: 1.2536, net: 544.1, sigma: 60.0, sigmaZero: sigmaZero,
                                  criticalLevel: 1.645 * sigmaZero, detectionLimit: 172.0, conflicts: [], suggestedRole: .fitOnly,
                                  holeRegionNote: nil, misfit: 1,
                                  possibleReason: "Adding it improves the fit by 81, below the 100 a new line needs.")
        let r = ProposalResult(candidates: [mg], sumPeaks: [], refused: [], currie: .standard, notes: [], passes: 1, settled: true)
        #expect(mg.aboveDetectionLimit && !mg.isProposed)
        #expect(r.possible == [mg] && r.proposed.isEmpty && r.sumPeakQuestions.isEmpty)
    }
}
