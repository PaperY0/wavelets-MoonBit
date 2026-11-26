# Project Agents.md Guide

This is a [MoonBit](https://docs.moonbitlang.com) project: a port of PyWavelets
(`PaperY0/wavelets`). See `README.md` for the API and `docs/known-deviations.md`
for the deliberate differences from upstream.

## Project Structure

- MoonBit packages are organized per directory; each directory contains a
  `moon.pkg` file listing its dependencies. Each package has its files and
  blackbox test files (ending in `_test.mbt`) and whitebox test files (ending in
  `_wbtest.mbt`).
- In the toplevel directory, there is a `moon.mod` file listing module metadata.
- The library is the root package; the CLI lives in `cmd/main`.
- `wavelet_data.mbt` and `golden_data_test.mbt` are **generated** — do not edit
  by hand. Regenerate with `python tools/gen_wavelets.py` and
  `python tools/gen_golden.py` (both need `tools/ref/pylibs`, see
  `tools/ref/requirements.txt`).

## Coding convention

- MoonBit code is organized in block style, each block is separated by `///|`,
  the order of each block is irrelevant.
- Public enums and the error type are `pub(all)` so downstream packages can
  construct/raise them; internal structs are `priv`.
- Float literals in scientific notation must carry a decimal point (`1.0e-9`,
  not `1e-9`). `@math` has no `sqrt`/`abs`; use `Double::sqrt(x)` / `x.sqrt()`.
- A package-level `using @json {type Json}` may appear only once across the
  package's files (it lives in `helpers_test.mbt` for the test package).

## Tooling

- `moon fmt` formats code; `moon info` regenerates `pkg.generated.mbti`.
- In the last step run `moon info && moon fmt`, then check the `.mbti` diff: if
  nothing changed, the refactor is invisible to external users.
- `moon test` runs the suite. The golden differential tests compare against
  PyWavelets 1.8.0 outputs baked into `golden_data_test.mbt`; the property tests
  verify perfect reconstruction across all wavelets/modes/lengths.
- Prefer assertion tests for these well-defined scientific results. Use
  `moon coverage analyze > uncovered.log` to inspect uncovered lines.
- After regenerating data modules, run `moon fmt` again (generators are not
  guaranteed fmt-clean).
- `moon publish` packs the working directory (not just git-tracked files), so the
  ~50MB `tools/ref/pylibs` reference install must be moved out of the module
  temporarily before publishing and restored afterwards.