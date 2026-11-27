# Status: SpectroKit 0.17.3 <!-- x-release-please-version -->

It builds on FormatoKit `"0.12.0"..<"1.0.0"` and ImagoKit `"0.6.1"..<"1.0.0"`, against HyperSpy 2.4.0 and
exSpy 0.3.2 / 7185a4d1. macOS 15 and later (iOS was dropped from `Package.swift` in 0.13.0; it never built there).

Live:
- **EDX**: the line table, element identification (the published identification loop, `ElementProposer`, which refines the energy
  axis on its own K picks, from 0.6 keV up so C and O K-alpha do not tilt it; the opt-in downward confirmation rule and the uncorroborated tier; carbon below 0.45 keV and a sum-peak
  question that another line answers are shown as possible, never added), net counts, the EDS fit with
  escape peaks and, by default where the data support them, sum peaks; the continuum step at absorption edges broadened by the
  detector by default, an optional Si dead-layer term, the opt-in NIST EPQ line weights; Cliff-Lorimer, zeta and cross-section
  quantification with the absorption iteration, per pixel over a spectrum image, k-factor import, sigma on wt% and at% (with
  take-off angle, thickness and dose terms) and per-pixel sigma maps, the per-region result; exSpy's window path for one SEM spectrum
  (`EDSWindowPath`: the lines `add_lines()` picks at a beam energy, integration and background windows, line intensities, the
  acquisition conditions with exSpy's defaults and their take-off angle, Cliff-Lorimer with the single-line rule named), on the app's
  TM002 bundles.
- **EELS**: t/lambda thickness, zero-loss alignment, power-law extrapolation, Fourier-log, Fourier-ratio and Richardson-Lucy
  deconvolution, the edge model fit (joint, or per-edge windows with a 1σ per edge intensity) on 383 bundled K, L and M edges, the MLLS
  fit, composition in at% and wt%, region sums.
- **SpectroCore**: line-profile pooling, display smoothing, sparse regions; a spectrum image's counts a band of scan rows at a time
  (`countRows`, exact, dense or sparse), for a writer that streams the cube.

Since 0.16.0: 0.17.0 adds `countRows`, which hands out a spectrum image's counts a band of scan rows at a time; 0.17.1 refuses a
power-law extrapolation past half the machine's memory before allocating it and takes FormatoKit and ImagoKit from the EMSuite org. What each
version changed, with the migration notes: `CHANGELOG.md`.

Parity row by row: `docs/parity.md` (92 rows: 59 present, 10 partial, 23 missing);
timings against Python: `docs/benchmarks.md`; what is open: the GitHub issues; known limits: `docs/parity.md`; fixtures and their generators: `Tests/Fixtures/README.md`.

Tests: the fast targets run per change (`swift test -c release --parallel --skip SpectroKitParityTests`); the pinned Python comparisons
(golden baseline, exSpy/HyperSpy/pyEELSMODEL pins, the planted EDX runs) are `SpectroKitParityTests`, run in CI on every push and before a tag.
