# Known deviations from PyWavelets

This port targets **numerical equivalence** with PyWavelets 1.8.0 for the
ported surface. Every deliberate difference is listed here; anything not listed
is intended to match upstream bit-for-bit (verified by the golden differential
tests in `golden_data_test.mbt`).

## Scope reductions (not implemented in v0.1)

1. **Complex continuous wavelets.** The `cgau`, `shan`, `fbsp`, `cmor` families
   are complex-valued. v0.1 ports only the 10 real continuous wavelets (`mexh`,
   `morl`, `gaus1`–`gaus8`). `Wavelet::new("cgau1")` raises `UnknownWavelet`.
   Rationale: the CWT/`wavefun`/`central_frequency` machinery is real-valued;
   complex support needs a separate coefficient type and is roadmap.

2. **Stationary wavelet transform (SWT/`swt`), wavelet packets (`WaveletPacket`),
   the FFT-based CWT method (`method="fft"`), and custom filter banks
   (`Wavelet` from user filters).** Not in v0.1. The CWT here is the convolution
   method only.

## Behavioral differences (ported surface, intentional)

3. **Over-level decomposition raises instead of warning.** Upstream
   `wavedec`/`wavedec2` emit a `UserWarning` when `level > dwt_max_level(...)`
   and proceed (all coefficients suffer boundary effects). MoonBit has no
   warning channel; this port raises `InvalidLevel`. Callers should query
   `dwt_max_level`/`wmaxlev` first. The golden cases use only valid levels.

4. **`central_frequency` is continuous-only.** Upstream also accepts discrete
   wavelets (via `wavefun`). Here `central_frequency`/`scale2frequency`/
   `frequency2scale` raise `NotContinuous` for discrete wavelets. `wavefun`
   itself works for both kinds.

5. **`integrate_wavelet` for discrete biorthogonal wavelets raises
   `NotSupported`.** Upstream integrates the decomposition-side wavelet of a
   biorthogonal pair via a cascade that this port does not reproduce in v0.1.
   Orthogonal discrete and all continuous wavelets integrate normally.

6. **`wavefun` level cap = 20.** Upstream allows arbitrary precision; this port
   rejects `level > 20` with `InvalidArgument` to bound the `2^level` grid.

## Numerical quirks preserved (NOT deviations — matched on purpose)

- **`dmey` is an FIR approximation of the Meyer wavelet, not exactly
  orthogonal.** Its `dec_lo` energy is `1.0022448292411643` (not 1), and
  `idwt(dwt(x))` carries a reconstruction error up to ~0.034. This is true of
  upstream too; the port reproduces it exactly. The property tests exclude
  `dmey` from the exact perfect-reconstruction sweep and assert the approximate
  bound separately.
- **`waverec` can return more samples than the input** (e.g. N=9 with db2 → 10)
  because upstream trims only the `cA == cD + 1` case at each step. Preserved.
- **Odd-length periodization** copies the last sample to make the length even
  before the periodic wrap; the periodized extension pads `(F-2)/2` on each
  side. Matched for all 106 discrete wavelets.
- **`soft`/`garrote` thresholding at `x == 0, value == 0` yields NaN** (`0/0`),
  and negative inputs clipped to zero keep `-0.0`. Reproduced via the
  `x * clip(1 - t/|x|)` form rather than `sign(x) * max(|x| - t, 0)`.
- **`greater`/`less` thresholding compare the signed value**, not the magnitude.
- **`dwt` requires N≥2 for `reflect`/`antireflect`** (period `2N-2`), raising
  `InvalidArgument` otherwise — exactly upstream's restriction. `pad` and `idwt`
  have no such restriction (n=1 reflect pads to a constant), also matching upstream.
- **`upcoef` detail branch** uses `rec_hi` on the first refinement and `rec_lo`
  thereafter; **`downcoef`** keeps the approximation branch on intermediate
  levels and selects `part` only at the final level. Both match upstream.

## JSON / CLI serialization

- The CLI emits JSON; non-finite doubles (NaN/±Inf) are written as `null`
  (matching `JSON.stringify`). The MoonBit API returns the raw `Double`.
- The pinned wasm-gc host routes `eprintln` to stdout; machine consumers must
  check the process exit code before parsing output as JSON.