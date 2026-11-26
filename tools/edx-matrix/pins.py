"""The EDX reference matrix: exSpy's window route from spectrum to composition on many materials.

For each case a spectrum is planted with exSpy's own line model (`EDSTEMSpectrum.create_model`, Gaussian lines at exSpy's
table weights and its FWHM law) plus a Kramers continuum, then exSpy 0.3.2 runs what a user runs:

    s.set_elements(elements); s.add_lines()                       # exSpy picks the lines
    iw = s.estimate_integration_windows()                         # 2 FWHM
    bw = s.estimate_background_windows()                          # line_width [2, 2], windows_width 1, merged
    I  = s.get_lines_intensity(background_windows=bw, integration_windows=iw)
    s.quantification(I, "CL", k, absorption_correction=True, thickness=t)   # weight and atomic percent

Every expected number in the fixture is exSpy's output. The planted composition is stored beside it as truth for the
record; the Swift test holds the kit to exSpy, not to the truth.

Matrix: 12 materials x thickness {10, 50, 150} nm x two detector set-ups x {80, 300} kV x two doses = 288 cases.
  set-up A: tilt 0, elevation 18 deg, 130 eV at Mn Ka;  set-up B: tilt 15, elevation 30 deg, 136 eV.
  dose "high": the expected spectrum scaled to 1e7 counts and rounded;  dose "low": Poisson at 1e5 counts (seeded).
Axis: 2560 channels, offset -0.2 keV, 10 eV.

k-factors are planted values, not physical ones (a user brings their own): K lines 1 + 0.04 (Z - 13), at least 0.7;
L lines 1.5 + 0.03 (Z - 40). The planted line areas are C / k times exSpy's own thin-film transmission (1 / acf of
`get_abs_corr_zeta` at the planted composition and `CL_get_mass_thickness`), so CL with absorption returns the
planted composition up to window spill and noise.

    References/python/bin/python -I -B tools/edx-matrix/pins.py      # writes Tests/SpectroEDXTests/Fixtures/edx-matrix.{json,u32}
"""
import hashlib, json, logging, pathlib, sys, warnings
warnings.filterwarnings("ignore"); logging.getLogger("hyperspy").setLevel(logging.ERROR)
import numpy as np
import hyperspy, exspy
from exspy.signals import EDSTEMSpectrum
from exspy._misc.eds import utils as eds_utils
from exspy._misc.elements import elements as edb

import exspy.signals.eds_tem as eds_tem_module


class _Grab(logging.Handler):            # exSpy logs "Convergence reached after N interations" at INFO: N is pinned
    def __init__(self): super().__init__(); self.n = None
    def emit(self, rec):
        msg = rec.getMessage()
        if msg.startswith("Convergence reached after"): self.n = int(msg.split()[3])


_grab = _Grab()
eds_tem_module._logger.addHandler(_grab); eds_tem_module._logger.setLevel(logging.INFO); eds_tem_module._logger.propagate = False

assert exspy.__version__ == "0.3.2" and hyperspy.__version__ == "2.4.0" and np.__version__ == "2.5.4", (exspy.__version__, hyperspy.__version__, np.__version__)

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "Tests/SpectroEDXTests/Fixtures"
SIZE, OFFSET, SCALE = 2560, -0.2, 0.01
SEED = 20261126

# name -> {element: atoms per formula unit} (converted to weight percent with exSpy's atomic weights)
MATERIALS = {
    "Al-4Cu":    {"wt": {"Al": 96.0, "Cu": 4.0}},
    "SrTiO3":    {"at": {"Sr": 1, "Ti": 1, "O": 3}},
    "FePt":      {"at": {"Fe": 1, "Pt": 1}},
    "NiAl":      {"at": {"Ni": 1, "Al": 1}},
    "Ti-6Al-4V": {"wt": {"Ti": 90.0, "Al": 6.0, "V": 4.0}},
    "316 steel": {"wt": {"Fe": 68.0, "Cr": 17.0, "Ni": 12.0, "Mo": 3.0}},
    "GaAs":      {"at": {"Ga": 1, "As": 1}},
    "SiGe":      {"at": {"Si": 1, "Ge": 1}},
    "ZrO2":      {"at": {"Zr": 1, "O": 2}},
    "MgAl2O4":   {"at": {"Mg": 1, "Al": 2, "O": 4}},
    "CuZn":      {"wt": {"Cu": 70.0, "Zn": 30.0}},
    "AuAg":      {"wt": {"Au": 50.0, "Ag": 50.0}},
}
THICKNESS = (10.0, 50.0, 150.0)
SETUPS = {"A": dict(tilt=0.0, elevation=18.0, resolution=130.0), "B": dict(tilt=15.0, elevation=30.0, resolution=136.0)}
VOLTAGES = (80.0, 300.0)
DOSES = {"high": 1e7, "low": 1e5}


def weight_percent(m):
    if "wt" in m: return dict(m["wt"])
    mass = {e: n * edb[e]["General_properties"]["atomic_weight"] for e, n in m["at"].items()}
    tot = sum(mass.values())
    return {e: 100 * v / tot for e, v in mass.items()}


def k_factor(line):
    el, fam = line.split("_")
    z = edb[el]["General_properties"]["Z"]
    if fam.startswith("K"): return max(0.7, 1 + 0.04 * (z - 13))
    if fam.startswith("L"): return max(0.9, 1.5 + 0.03 * (z - 40))
    return 3.0                                                        # M lines


def spectrum(elements, kv, setup):
    s = EDSTEMSpectrum(np.zeros(SIZE))
    ax = s.axes_manager.signal_axes[0]; ax.scale = SCALE; ax.offset = OFFSET; ax.units = "keV"; ax.name = "Energy"
    s.set_microscope_parameters(beam_energy=kv, tilt_stage=setup["tilt"], azimuth_angle=0.0,
                                elevation_angle=setup["elevation"], energy_resolution_MnKa=setup["resolution"],
                                live_time=1.0, beam_current=1.0)
    s.set_elements(elements); s.add_lines()
    return s


def continuum(energy, kv, zbar):
    e = np.clip(energy, 0.05, None)
    k = np.where(energy > 0.1, zbar * np.clip(kv - e, 0, None) / e, 0.0)
    eff = np.exp(-0.25 / e ** 3) * (1 - np.exp(-60.0 / e ** 2.5))   # window cut-off and a thin Si crystal
    return k * eff


def plant(s, wt, thickness, kv, rng, dose):
    lines = list(s.metadata.Sample.xray_lines)
    elements = [l.split("_")[0] for l in lines]
    w = [wt[e] for e in elements]
    toa = s.get_take_off_angle()
    wsig = [hyperspy.signals.BaseSignal(np.array([v])) for v in w]
    for sig, e, l in zip(wsig, elements, lines):
        sig.metadata.set_item("Sample.elements", [e]); sig.metadata.set_item("Sample.xray_lines", [l])
    mt = s.CL_get_mass_thickness(wsig, float(thickness))
    mts = hyperspy.signals.BaseSignal(np.asarray(mt, dtype=float).reshape(1))
    acf = np.asarray(eds_utils.get_abs_corr_zeta(wsig, mts, toa)).ravel()
    assert np.all(np.isfinite(acf)), (acf, mt, toa)
    m = s.create_model(auto_background=False)
    areas = {}
    for comp in m.xray_lines:
        if comp.name in lines:                                        # the measured line: absorbed as exSpy models it
            i = lines.index(comp.name)
            areas[comp.name] = w[i] / 100 / k_factor(comp.name) / acf[i]
        else:                                                         # the element's other families: present, not absorbed
            areas[comp.name] = wt[comp.name.split("_")[0]] / 100 / k_factor(comp.name)
        comp.A.value = areas[comp.name]
    m.assign_current_values_to_all()
    peaks = np.asarray(m.as_signal().data, dtype=float)
    energy = OFFSET + SCALE * np.arange(SIZE)
    zbar = sum(wt[e] / 100 * edb[e]["General_properties"]["Z"] for e in wt)
    bg = continuum(energy, kv, zbar)
    bg *= 2.0 * peaks.sum() / bg.sum()                                # continuum = 2 x the line counts
    mean = peaks + bg
    assert np.all(np.isfinite(mean)) and mean.min() >= 0, (peaks.min(), bg.min(), np.isfinite(peaks).all(), np.isfinite(bg).all())
    mean *= dose[1] / mean.sum()
    counts = np.rint(mean) if dose[0] == "high" else rng.poisson(mean).astype(float)
    return counts.astype(np.uint32), [float(a) for a in acf], float(np.asarray(mt).ravel()[0])


def run_exspy(s, counts, lines, thickness):
    s.data = counts.astype(float)
    iw = s.estimate_integration_windows()
    bw = s.estimate_background_windows()
    I = s.get_lines_intensity(background_windows=bw, integration_windows=iw)
    names = [i.metadata.Sample.xray_lines[0] for i in I]
    assert names == lines, (names, lines)
    k = [k_factor(l) for l in lines]
    out = dict(integrationWindows=np.asarray(iw).tolist(), backgroundWindows=np.asarray(bw).tolist(),
               intensities=[float(np.asarray(i.data).ravel()[0]) for i in I], kFactors=k)
    for units, key in (("weight", "weightPercent"), ("atomic", "atomicPercent")):
        _grab.n = None
        comp, mt = s.quantification(I, "CL", k, composition_units=units, absorption_correction=True,
                                    thickness=float(thickness), show_progressbar=False)
        out[key] = [float(np.asarray(c.data).ravel()[0]) for c in comp]
        out["massThickness"] = float(np.asarray(mt.data).ravel()[0])
        out["iterations"] = _grab.n
    comp = s.quantification(I, "CL", k, composition_units="weight", show_progressbar=False)
    out["thinWeightPercent"] = [float(np.asarray(c.data).ravel()[0]) for c in comp]
    return out


def main():
    rng = np.random.default_rng(SEED)
    cases, blobs = [], []
    for name, m in MATERIALS.items():
        wt = weight_percent(m)
        for setup_name, setup in SETUPS.items():
            for kv in VOLTAGES:
                for t in THICKNESS:
                    for dose in DOSES.items():
                        s = spectrum(sorted(wt), kv, setup)
                        lines = list(s.metadata.Sample.xray_lines)
                        counts, acf, mt = plant(s, wt, t, kv, rng, dose)
                        r = run_exspy(s, counts, lines, t)
                        planted = [wt[l.split("_")[0]] for l in lines]
                        cases.append(dict(
                            name=f"{name} {setup_name} {kv:.0f}kV {t:.0f}nm {dose[0]}", material=name, setup=setup_name,
                            beamEnergyKeV=kv, thicknessNm=t, dose=dose[0], resolutionMnKaEV=setup["resolution"],
                            tiltDegrees=setup["tilt"], elevationDegrees=setup["elevation"],
                            takeOffDegrees=float(s.get_take_off_angle()), lines=lines,
                            elements=[l.split("_")[0] for l in lines],
                            lineEnergiesKeV=[float(eds_utils._get_energy_xray_line(l)) for l in lines],
                            countsOffset=len(blobs) * SIZE, plantedWeightPercent=planted,
                            plantedMassThickness=mt, plantedAbsorptionFactors=acf, **r))
                        blobs.append(counts)
    data = np.concatenate(blobs).astype("<u4").tobytes()
    (OUT / "edx-matrix.u32").write_bytes(data)
    doc = dict(generator="tools/edx-matrix/pins.py", exspy=exspy.__version__, hyperspy=hyperspy.__version__, numpy=np.__version__,
               axis=dict(size=SIZE, offset=OFFSET, scale=SCALE), seed=SEED,
               countsFile="edx-matrix.u32", countsSha256=hashlib.sha256(data).hexdigest(), cases=cases)
    (OUT / "edx-matrix.json").write_text(json.dumps(doc, indent=1) + "\n")
    err = [max(abs(a - b) for a, b in zip(c["weightPercent"], c["plantedWeightPercent"])) for c in cases]
    for dose in DOSES:
        e = [x for x, c in zip(err, cases) if c["dose"] == dose]
        print(f"{dose}: {len(e)} cases, max |exSpy - planted| wt% median {np.median(e):.3f}, max {max(e):.3f}")
    worst = sorted(zip(err, [c["name"] for c in cases]))[-6:]
    print("worst:", worst)
    print(len(cases), "cases,", len(data) / 1e6, "MB counts")


if __name__ == "__main__":
    main()
