"""Reverse-engineer exact pywt semantics for the MoonBit port.

Run: python tools/probe.py   (uses tools/ref/pylibs)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "ref", "pylibs"))
import numpy as np
import pywt

print("== 1. discrete wavefun cascade ==")
for name in ["haar", "db2"]:
    w = pywt.Wavelet(name)
    for lvl in [1, 2, 3]:
        phi, psi, x = w.wavefun(lvl)
        print(name, lvl, "len", len(phi), "x", x[0], x[-1], "phi[:6]", np.round(phi[:6], 6))

print("== 2. gaus normalization constants ==")
for n in range(1, 9):
    w = pywt.ContinuousWavelet(f"gaus{n}")
    # analytic nth derivative of exp(-t^2) via physicists' Hermite: (-1)^n H_n(t) e^{-t^2}
    def ana(x, n):
        Hprev, H = np.ones_like(x), 2 * x
        if n == 0:
            return np.exp(-x * x)
        for k in range(1, n):
            Hprev, H = H, 2 * x * H - 2 * k * Hprev
        return ((-1) ** n) * H * np.exp(-x * x)
    consts = []
    for lvl in [8, 12]:
        psi, x = w.wavefun(lvl)
        a = ana(x, n)
        c = float(np.sum(psi * a) / np.sum(a * a))
        consts.append(c)
        err = float(np.max(np.abs(psi - c * a)))
        print(f"gaus{n} lvl{lvl}: c={c:.17g} maxerr={err:.3e}")
    print(f"  level-independent: {abs(consts[0]-consts[1]) < 1e-14}")

print("== 3. periodization idwt ==")
x8 = np.arange(1.0, 9.0)
a, d = pywt.dwt(x8, "db2", mode="periodization")
r = pywt.idwt(a, d, "db2", mode="periodization")
print("db2 per idwt len:", len(r), "roundtrip err:", float(np.max(np.abs(r - x8))))
x7 = np.arange(7.0)
a7, d7 = pywt.dwt(x7, "haar", mode="periodization")
r7 = pywt.idwt(a7, d7, "haar", mode="periodization")
print("haar per N=7 rec:", r7)

print("== 4. dwt2 orientation ==")
m = np.array([[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0], [9.0, 10.0, 11.0, 12.0], [13.0, 14.0, 15.0, 16.0]])
cA, (cH, cV, cD) = pywt.dwt2(m, "haar")
print("cA:", cA.tolist())
print("cH:", cH.tolist())
print("cV:", cV.tolist())
print("cD:", cD.tolist())
# manual: rows first (axis=1) then cols (axis=0)?
lo_r, hi_r = pywt.dwt(m, "haar", axis=1)
ll, lh = pywt.dwt(lo_r, "haar", axis=0)
hl, hh = pywt.dwt(hi_r, "haar", axis=0)
print("rows-first LL==cA:", np.allclose(ll, cA), "LH==cH:", np.allclose(lh, cH), "HL==cV:", np.allclose(hl, cV), "HH==cD:", np.allclose(hh, cD))

print("== 5. dwt_max_level edges ==")
for n, f in [(1, 2), (2, 4), (3, 4), (4, 4), (1000, 10), (0, 2)]:
    print(n, f, pywt.dwt_max_level(n, f))

print("== 6. threshold boundaries ==")
xs = np.array([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0])
for mode in ["soft", "hard", "greater", "less", "garrote"]:
    print(mode, pywt.threshold(xs, 2.0, mode=mode).tolist())

print("== 7. pad antisymmetric / smooth / antireflect ==")
x = np.array([1.0, 2.0, 3.0])
for mode in ["antisymmetric", "smooth", "antireflect", "constant", "periodic"]:
    print(mode, pywt.pad(x, (4, 4), mode).tolist())
print("smooth len1:", pywt.pad(np.array([5.0]), (2, 2), "smooth").tolist())
print("symmetric len1:", pywt.pad(np.array([5.0]), (2, 2), "symmetric").tolist())

print("== 8. wavedec default level ==")
c = pywt.wavedec(np.arange(16.0), "db2")
print("N=16 db2 default level:", len(c) - 1, "lens:", [len(a) for a in c])

print("== 9. idwt cA longer by 1 ==")
try:
    r = pywt.idwt(np.zeros(6), np.zeros(5), "db2")
    print("idwt(6,5) len:", len(r))
except Exception as e:
    print("idwt(6,5) error:", e)

print("== 10. continuous wavefun x grid ==")
psi, x = pywt.ContinuousWavelet("mexh").wavefun(3)
print("mexh lvl3 len:", len(psi), "x:", x.tolist())