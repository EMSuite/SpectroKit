# tools/edx-matrix: exSpy's EDX window route end to end on twelve materials

`pins.py` plants 288 spectra with exSpy 0.3.2's own line model and a Kramers continuum. It runs exSpy's window route on each:
lines, integration and background windows, net counts, then Cliff-Lorimer thin and with the absorption iteration. It writes
`Tests/SpectroEDXTests/Fixtures/edx-matrix.json` (exSpy's answers, the planted truth, the layout) and `edx-matrix.u32` (the counts,
2560 channels per case, uint32 little-endian). The script's docstring defines the matrix: materials, thicknesses, detector set-ups,
voltages, doses and k-factors.

    tools/reference/setup.sh
    References/python/bin/python -I -B tools/edx-matrix/pins.py      # about 5 minutes

`EDXMatrixTests` runs the kit's chain on the same counts and holds it to exSpy at exSpy's own test tolerances. The
planted composition is a record, not a test. exSpy's window route misses it by a median of 0.3 wt% at 1e7 counts and
0.6 wt% at 1e5, but by up to 11 wt% where a line of another family falls in a window. Examples: Ti L in the O K windows
of SrTiO3, and Ti K-beta on V K-alpha in Ti-6Al-4V. That is the window method's limit, and the reason the app fits the
spectrum.
