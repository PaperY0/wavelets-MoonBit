"""Generate golden test vectors from PyWavelets and embed them in golden_data.mbt.

Run from the project root:  python tools/gen_golden.py

Every vector is computed by the reference implementation (tools/ref/pylibs) so
the MoonBit port can be checked against exact upstream behaviour, including
edge cases (odd lengths, periodization, boundary thresholds).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "ref", "pylibs"))
import numpy as np
import pywt

rng = np.random.default_rng(42)

MODES = ["zero", "constant", "symmetric", "periodic", "smooth",
         "periodization", "reflect", "antisymmetric", "antireflect"]


def _clean(v):
    v = float(v)
    # NaN (e.g. garrote at t=0, x=0) is encoded as null; tests assert is_nan.
    return None if v != v else v


def L(a):
    return [_clean(v) for v in np.asarray(a).ravel()]


def M(a):
    return [[_clean(v) for v in row] for row in np.asarray(a)]


def main():
    g = {}

    signals = {
        "ramp1": [1.0],
        "ramp2": [1.0, 2.0],
        "ramp3": [1.0, 2.0, 3.0],
        "ramp7": list(np.arange(1.0, 8.0)),
        "ramp8": list(np.arange(1.0, 9.0)),
        "ramp9": list(np.arange(1.0, 10.0)),
        "ramp16": list(np.arange(1.0, 17.0)),
        "rand32": L(rng.standard_normal(32)),
    }

    # ---- dwt ----
    dwt_cases = []
    for w in ["haar", "db1", "db2", "db4", "db10", "sym4", "sym8", "coif1",
              "coif3", "bior1.3", "bior3.5", "rbio2.2", "dmey"]:
        for mode in MODES:
            for key in ["ramp1", "ramp2", "ramp3", "ramp7", "ramp8", "ramp9",
                        "ramp16", "rand32"]:
                x = signals[key]
                try:
                    a, d = pywt.dwt(x, w, mode=mode)
                except Exception:
                    continue
                dwt_cases.append({"w": w, "mode": mode, "x": x,
                                  "a": L(a), "d": L(d)})
    g["dwt"] = dwt_cases

    # ---- idwt (roundtrip + single-branch) ----
    idwt_cases = []
    for w in ["haar", "db2", "bior3.5", "dmey", "sym4"]:
        for mode in MODES:
            for key in ["ramp8", "ramp9", "rand32"]:
                x = signals[key]
                a, d = pywt.dwt(x, w, mode=mode)
                r = pywt.idwt(a, d, w, mode=mode)
                idwt_cases.append({"w": w, "mode": mode, "a": L(a), "d": L(d),
                                   "r": L(r)})
                ra = pywt.idwt(a, None, w, mode=mode)
                rd = pywt.idwt(None, d, w, mode=mode)
                idwt_cases.append({"w": w, "mode": mode, "a": L(a), "d": None,
                                   "r": L(ra)})
                idwt_cases.append({"w": w, "mode": mode, "a": None, "d": L(d),
                                   "r": L(rd)})
    g["idwt"] = idwt_cases

    # ---- wavedec / waverec ----
    ml_cases = []
    for w in ["haar", "db2", "sym4", "coif2", "bior3.5"]:
        for mode in ["symmetric", "periodization", "smooth", "reflect"]:
            for n in [16, 31, 32]:
                x = L(rng.standard_normal(n))
                maxlvl = pywt.dwt_max_level(n, pywt.Wavelet(w).dec_len)
                for level in range(1, min(3, maxlvl) + 1):
                    c = pywt.wavedec(x, w, mode=mode, level=level)
                    r = pywt.waverec(c, w, mode=mode)
                    ml_cases.append({"w": w, "mode": mode, "level": level,
                                     "x": x, "coeffs": [L(a) for a in c],
                                     "rec": L(r)})
    g["wavedec"] = ml_cases

    # ---- upcoef / downcoef ----
    ud_cases = []
    for w in ["haar", "db2", "sym4"]:
        x = signals["ramp8"]
        a, d = pywt.dwt(x, w, mode="symmetric")
        for part, c in [("a", a), ("d", d)]:
            for level in [1, 2]:
                for take in [0, 8]:
                    r = pywt.upcoef(part, c, w, level=level, take=take)
                    ud_cases.append({"op": "up", "w": w, "part": part,
                                     "level": level, "take": take,
                                     "c": L(c), "r": L(r)})
            for level in [1, 2]:
                r = pywt.downcoef(part, x, w, mode="symmetric", level=level)
                ud_cases.append({"op": "down", "w": w, "part": part,
                                 "level": level, "take": 0,
                                 "c": x, "r": L(r)})
    g["updown"] = ud_cases

    # ---- dwt2 / idwt2 ----
    mat_cases = []
    mats = {
        "m44": M(np.arange(1.0, 17.0).reshape(4, 4)),
        "m57": M(rng.standard_normal((5, 7))),
        "m86": M(rng.standard_normal((8, 6))),
    }
    for w in ["haar", "db2", "bior3.5"]:
        for mode in ["symmetric", "periodization", "zero"]:
            for key, m in mats.items():
                a, (h, v, d) = pywt.dwt2(m, w, mode=mode)
                r = pywt.idwt2((a, (h, v, d)), w, mode=mode)
                mat_cases.append({"w": w, "mode": mode, "x": m,
                                  "a": M(a), "h": M(h), "v": M(v), "d": M(d),
                                  "rec": M(r)})
    g["dwt2"] = mat_cases

    # ---- wavedec2 / waverec2 ----
    # 仅选用层数有效的矩阵/小波组合（本库对超层报错，见 known-deviations）。
    ml2_cases = []
    ml2_mats = {
        "m86": mats["m86"],
        "m1616": M(rng.standard_normal((16, 16))),
        "m1513": M(rng.standard_normal((15, 13))),
    }
    for w in ["haar", "db2"]:
        for mode in ["symmetric", "periodization"]:
            for key, level in [("m86", 2 if w == "haar" else 1),
                               ("m1616", 2), ("m1513", 2)]:
                m = ml2_mats[key]
                c = pywt.wavedec2(m, w, mode=mode, level=level)
                r = pywt.waverec2(c, w, mode=mode)
                ml2_cases.append({"w": w, "mode": mode, "level": level,
                                  "x": m, "approx": M(c[0]),
                                  "details": [[M(t) for t in tup] for tup in c[1:]],
                                  "rec": M(r)})
    g["wavedec2"] = ml2_cases

    # ---- cwt ----
    cwt_cases = []
    t = np.linspace(0.0, 1.0, 32, endpoint=False)
    chirp = L(np.cos(2 * np.pi * (3 * t + 8 * t * t)))
    for w in ["mexh", "morl", "gaus1", "gaus2", "gaus4"]:
        for scales in [[1.0, 2.0, 4.0, 8.0], [1.0, 3.0, 7.0]]:
            for x in [chirp, signals["ramp16"]]:
                coefs, freqs = pywt.cwt(x, scales, w, sampling_period=0.5,
                                        precision=12)
                cwt_cases.append({"w": w, "scales": scales, "x": x,
                                  "dt": 0.5, "precision": 12,
                                  "coefs": M(coefs), "freqs": L(freqs)})
    g["cwt"] = cwt_cases

    # ---- wavefun ----
    wf_cases = []
    for w in ["haar", "db2", "sym4", "coif1", "bior3.5", "rbio2.2", "dmey"]:
        levels = [1, 2, 3] if w != "dmey" else [1, 2]
        if w in ("haar", "db2"):
            levels.append(5)
        for level in levels:
            out = pywt.Wavelet(w).wavefun(level)
            if len(out) == 3:
                phi, psi, x = out
                wf_cases.append({"w": w, "level": level, "kind": "orthogonal",
                                 "phi": L(phi), "psi": L(psi), "x": L(x)})
            else:
                phi_d, psi_d, phi_r, psi_r, x = out
                wf_cases.append({"w": w, "level": level, "kind": "biorthogonal",
                                 "phi": L(phi_d), "psi": L(psi_d),
                                 "phi_r": L(phi_r), "psi_r": L(psi_r), "x": L(x)})
    for w in ["mexh", "morl", "gaus1", "gaus3", "gaus8"]:
        for level in [3, 8]:
            psi, x = pywt.ContinuousWavelet(w).wavefun(level)
            wf_cases.append({"w": w, "level": level, "kind": "continuous",
                             "psi": L(psi), "x": L(x)})
    g["wavefun"] = wf_cases

    # ---- central_frequency / integrate_wavelet ----
    cf_cases = []
    for w in ["mexh", "morl", "gaus1", "gaus2", "gaus3", "gaus4", "gaus5",
              "gaus6", "gaus7", "gaus8"]:
        for precision in [8, 12]:
            cf_cases.append({"w": w, "precision": precision,
                             "f": float(pywt.central_frequency(w, precision=precision))})
    g["centrfrq"] = cf_cases

    iw_cases = []
    for w in ["mexh", "morl", "gaus1", "haar", "db2"]:
        out = pywt.integrate_wavelet(w, precision=8)
        if len(out) == 2:
            int_psi, x = out
            iw_cases.append({"w": w, "precision": 8, "int_psi": L(int_psi), "x": L(x)})
    g["intwave"] = iw_cases

    # ---- pad ----
    pad_cases = []
    for mode in MODES:
        for x in [[1.0, 2.0, 3.0], [5.0], list(np.arange(1.0, 9.0))]:
            for left, right in [(4, 4), (0, 3), (2, 0)]:
                r = pywt.pad(x, (left, right), mode)
                pad_cases.append({"mode": mode, "x": x, "left": left,
                                  "right": right, "r": L(r)})
    g["pad"] = pad_cases

    # ---- threshold ----
    th_cases = []
    xs = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    for mode in ["soft", "hard", "greater", "less", "garrote"]:
        for tval in [0.0, 2.0]:
            r = pywt.threshold(xs, tval, mode=mode)
            th_cases.append({"mode": mode, "t": tval, "x": xs, "r": L(r)})
    g["threshold"] = th_cases

    # ---- max_level / coeff_len ----
    mlvl = []
    for n in [0, 1, 2, 3, 4, 7, 8, 16, 100, 1000]:
        for f in [2, 4, 10, 20]:
            mlvl.append({"n": n, "f": f, "r": int(pywt.dwt_max_level(n, f))})
    g["maxlevel"] = mlvl
    cl = []
    for n in [1, 2, 3, 7, 8, 9, 16, 31]:
        for f in [2, 4, 12, 62]:
            for mode in ["symmetric", "periodization"]:
                cl.append({"n": n, "f": f, "mode": mode,
                           "r": int(pywt.dwt_coeff_len(n, f, mode))})
    g["coefflen"] = cl

    # ---- qmf / orthogonal_filter_bank ----
    qmf_cases = []
    for x in [signals["ramp8"], pywt.Wavelet("db2").rec_lo]:
        qmf_cases.append({"x": L(x), "r": L(pywt.qmf(x))})
    g["qmf"] = qmf_cases
    ofb = []
    for w in ["db2", "sym4", "coif1"]:
        h = pywt.Wavelet(w).rec_lo
        dl, dh, rl, rh = pywt.orthogonal_filter_bank(h)
        ofb.append({"x": L(h), "dec_lo": L(dl), "dec_hi": L(dh),
                    "rec_lo": L(rl), "rec_hi": L(rh)})
    g["orthfilt"] = ofb

    payload = json.dumps(g, separators=(",", ":"), ensure_ascii=True)
    # embed as a MoonBit string literal (pure ASCII, JSON-escaped == MoonBit-escaped)
    mbt = ("// GENERATED by tools/gen_golden.py from PyWavelets "
           f"{pywt.__version__} — do not edit by hand.\n"
           "// Golden vectors for black-box differential tests.\n\n"
           "///|\n"
           "let golden_json : String = " + json.dumps(payload) + "\n")
    out = os.path.join(ROOT, "golden_data_test.mbt")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write(mbt)
    with open(os.path.join(HERE, "golden.json"), "w", encoding="utf-8") as f:
        f.write(payload)
    counts = {k: len(v) for k, v in g.items()}
    print("wrote", out, f"({len(payload)} bytes JSON)")
    print("cases:", counts)


if __name__ == "__main__":
    main()