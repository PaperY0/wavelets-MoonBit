"""End-to-end differential test: wavelets CLI vs a live PyWavelets 1.8.0.

The baked golden suite (golden_data_test.mbt) already compares the *library*
against PyWavelets exhaustively. This script additionally drives the *CLI*
(`moon run cmd/main`) so the JSON plumbing, argument parsing and serialization
are verified end-to-end. It writes DIFF-REPORT.md.

Setup (local, does not touch the global env):
    python -m pip install --target tools/ref/pylibs -r tools/ref/requirements.txt
Run:
    python tools/diff_test.py
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "ref", "pylibs"))
import numpy as np
import pywt

TOL = 1e-9
MODES = ["zero", "constant", "symmetric", "periodic", "smooth",
         "periodization", "reflect", "antisymmetric", "antireflect"]
# Representative wavelets spanning every family (the golden suite covers all 106).
DW = ["haar", "db2", "db4", "db8", "sym4", "sym8", "coif1", "coif3",
      "bior3.5", "rbio3.5", "dmey"]
CW = ["mexh", "morl", "gaus1", "gaus2", "gaus4"]

stats = {"total": 0, "match": 0, "mismatch": 0}
failures = []


def cli(*args):
    r = subprocess.run(["moon", "run", "cmd/main", "--", *args],
                       cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"CLI failed ({args[:3]}): {r.stderr[-300:]}")
    return json.loads(r.stdout.strip())


def nums(a):
    """CLI JSON numbers: null -> NaN."""
    return np.array([np.nan if v is None else float(v) for v in a], dtype=float)


def close(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return False
    an, bn = np.isnan(a), np.isnan(b)
    if not np.array_equal(an, bn):
        return False
    return bool(np.allclose(a[~an], b[~bn], rtol=TOL, atol=TOL))


def check(label, got, ref):
    stats["total"] += 1
    if close(got, ref):
        stats["match"] += 1
    else:
        stats["mismatch"] += 1
        if len(failures) < 40:
            failures.append(label)


def sig(n, seed=0):
    return np.arange(1.0, n + 1.0) + 0.25 * np.sin(np.arange(n) + seed)


def main():
    # Warm up the build so per-case `moon run` only executes the cached artifact.
    subprocess.run(["moon", "build", "cmd/main"], cwd=ROOT,
                   capture_output=True, text=True)

    # ---- dwt ----
    for w in DW:
        for mode in MODES:
            for n in [1, 3, 8, 9, 16, 31, 64]:
                x = sig(n)
                try:
                    a, d = pywt.dwt(x, w, mode=mode)
                except Exception:
                    # upstream rejects e.g. n=1 for [anti]reflect; skip likewise
                    continue
                got = cli("dwt", w, mode, *[repr(float(v)) for v in x])
                check(f"dwt {w}/{mode}/n={n} cA", nums(got["cA"]), a)
                check(f"dwt {w}/{mode}/n={n} cD", nums(got["cD"]), d)

    # ---- wavedec (valid levels only) ----
    for w in ["haar", "db2", "sym4", "coif1"]:
        for mode in ["symmetric", "periodization", "zero"]:
            x = sig(64)
            level = min(3, pywt.dwt_max_level(64, pywt.Wavelet(w).dec_len))
            ref = pywt.wavedec(x, w, mode=mode, level=level)
            got = cli("wavedec", w, mode, str(level), *[repr(float(v)) for v in x])
            for k, (g, r) in enumerate(zip(got["coeffs"], ref)):
                check(f"wavedec {w}/{mode}/L={level} coeff[{k}]", nums(g), r)

    # ---- cwt (CLI defaults: sampling_period=1, precision=12) ----
    t = np.linspace(0.0, 1.0, 32, endpoint=False)
    chirp = np.cos(2 * np.pi * (3 * t + 8 * t * t))
    for w in CW:
        for scales in [[1.0, 2.0, 4.0, 8.0], [1.0, 3.0, 7.0]]:
            for x in [chirp, sig(16)]:
                coefs, freqs = pywt.cwt(x, scales, w, sampling_period=1.0, precision=12)
                got = cli("cwt", w, ",".join(repr(s) for s in scales),
                          *[repr(float(v)) for v in x])
                check(f"cwt {w}/{scales} coefs", np.array(got["coefs"], dtype=float), coefs)
                check(f"cwt {w}/{scales} freqs", nums(got["freqs"]), freqs)

    # ---- wavefun ----
    for w in ["haar", "db2", "sym4", "bior3.5", "mexh", "morl", "gaus2"]:
        for level in [1, 3, 5]:
            if w in CW:
                psi, x = pywt.ContinuousWavelet(w).wavefun(level)
                got = cli("wavefun", w, str(level))
                check(f"wavefun {w}/L={level} psi", nums(got["psi"]), psi)
                check(f"wavefun {w}/L={level} x", nums(got["x"]), x)
            else:
                wv = pywt.Wavelet(w)
                got = cli("wavefun", w, str(level))
                out = wv.wavefun(level)
                if len(out) == 3:
                    phi, psi, x = out
                    check(f"wavefun {w}/L={level} phi", nums(got["phi"]), phi)
                    check(f"wavefun {w}/L={level} psi", nums(got["psi"]), psi)
                    check(f"wavefun {w}/L={level} x", nums(got["x"]), x)
                else:
                    phi_d, psi_d, phi_r, psi_r, x = out
                    check(f"wavefun {w}/L={level} phi", nums(got["phi"]), phi_d)
                    check(f"wavefun {w}/L={level} psi", nums(got["psi"]), psi_d)
                    check(f"wavefun {w}/L={level} phi_r", nums(got["phi_r"]), phi_r)
                    check(f"wavefun {w}/L={level} psi_r", nums(got["psi_r"]), psi_r)
                    check(f"wavefun {w}/L={level} x", nums(got["x"]), x)

    # ---- threshold ----
    xs = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    for mode in ["soft", "hard", "greater", "less", "garrote"]:
        for tval in [0.0, 2.0]:
            with np.errstate(invalid="ignore", divide="ignore"):
                ref = pywt.threshold(xs, tval, mode=mode)
            got = cli("threshold", mode, repr(tval), *[repr(v) for v in xs])
            check(f"threshold {mode}/t={tval}", nums(got["r"]), ref)

    # ---- pad ----
    x = sig(5)
    for mode in MODES:
        for left, right in [(0, 0), (3, 4), (7, 2)]:
            ref = pywt.pad(x, (left, right), mode)
            got = cli("pad", mode, str(left), str(right), *[repr(float(v)) for v in x])
            check(f"pad {mode}/({left},{right})", nums(got["r"]), ref)

    # ---- report ----
    lines = [
        "# Differential report — CLI vs PyWavelets 1.8.0",
        "",
        f"Generated by `tools/diff_test.py` on {os.environ.get('USERNAME', 'local')}.",
        f"Tolerance: rtol=atol={TOL}; NaN matched against JSON null.",
        "",
        f"- Total checks: **{stats['total']}**",
        f"- Matched: **{stats['match']}**",
        f"- Mismatched: **{stats['mismatch']}**",
        "",
        "The CLI uses default `sampling_period=1`, `precision=12` for CWT; the",
        "reference calls match those defaults. Over-level decomposition is not",
        "exercised here (the library raises by design; see docs/known-deviations.md).",
        "",
    ]
    if failures:
        lines += ["## Mismatches (first 40)", ""] + [f"- {f}" for f in failures]
    else:
        lines += ["**0 unexpected mismatches.** Every CLI result agrees with PyWavelets.", ""]
    report = "\n".join(lines)
    with open(os.path.join(ROOT, "DIFF-REPORT.md"), "w", encoding="utf-8") as f:
        f.write(report)
    print(report)
    return 1 if stats["mismatch"] else 0


if __name__ == "__main__":
    sys.exit(main())