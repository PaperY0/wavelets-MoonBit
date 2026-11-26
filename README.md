# PaperY0/wavelets

> [PyWavelets](https://github.com/PyWavelets/pywt) 的 MoonBit 移植。
> 输入信号，输出 **1D/2D 离散小波变换（DWT/IDWT、多层分解）+ 连续小波变换（CWT）+
> 小波函数/尺度函数 + 阈值化 + 9 种信号延拓模式**，内置 **106 个离散小波 + 10 个实值连续小波**。
> 纯计算、无 IO，`wasm` / `wasm-gc` / `js` 三后端可用。

[![CI](https://github.com/PaperY0/wavelets-MoonBit/actions/workflows/ci.yml/badge.svg)](https://github.com/PaperY0/wavelets-MoonBit/actions/workflows/ci.yml)
[English below](#papery0wavelets-en)

## 为什么是它

MoonBit 生态里**没有小波变换库**。对 mooncakes 注册表全量 3983 个模块（剔除测试/撤回后
3838 个）逐一扫描名称与描述，`wavelet` / `小波` / `DWT` / `CWT` / `PyWavelets` /
`daubechies` / `morlet` 等关键词**零命中**：现存 DSP 模块只做 FFT/FIR/IIR，唯一的
JPEG2000 解码器把逆变换封在 codec 内部、并不对外暴露小波 API。小波分析（时频局部化、
多分辨率分解、去噪、压缩、特征提取）在生态里是空白。

本库把 Python 事实标准 **PyWavelets 1.8.0** 完整移植过来，延续作者其它移植项目
（`zxcvbn`、`crf`、`conllu`、`pronouncing`）的范式：**成熟算法库移植 + 纯计算 +
对参考实现逐用例差分验证**。滤波器系数与上游**逐位一致**（17 位有效数字，脚本生成）。

## 能力范围

- **106 个离散小波**：`haar`、`db1`–`db38`、`sym2`–`sym20`、`coif1`–`coif17`、
  `bior`/`rbio` 双正交族、`dmey`；**10 个实值连续小波**：`mexh`、`morl`、`gaus1`–`gaus8`。
- **1D DWT/IDWT**：全部 **9 种延拓模式**（zero/constant/symmetric/periodic/smooth/
  periodization/reflect/antisymmetric/antireflect），忠实保留上游的奇数长度
  periodization 与“waverec 可能比输入多一个样本”等怪癖。
- **多层分解**：`wavedec`/`waverec`、`dwt_max_level`/`wmaxlev`、`dwt_coeff_len`；
  **部分变换** `upcoef`/`downcoef`（多层分支语义与上游一致）、`qmf`、`orthogonal_filter_bank`。
- **2D**：`dwt2`/`idwt2`/`wavedec2`/`waverec2`，轻量行主序 `Matrix` 类型。
- **CWT（卷积法）**：`cwt`、`central_frequency`、`scale2frequency`/`frequency2scale`、
  `integrate_wavelet`、`wavefun`（正交 + 双正交两侧）。
- **阈值化**：`threshold`（soft/hard/greater/less/garrote），复刻上游精确边界语义
  （含 `0/0 → NaN` 与负零保留）。

## API 设计

```moonbit nocheck
// 1D 单层 / 多层
pub fn dwt(Array[Double], Wavelet, mode? : Mode) -> (Array[Double], Array[Double]) raise WaveletError
pub fn idwt(Array[Double]?, Array[Double]?, Wavelet, mode? : Mode) -> Array[Double] raise WaveletError
pub fn wavedec(Array[Double], Wavelet, mode? : Mode, level? : Int?) -> Array[Array[Double]] raise WaveletError
pub fn waverec(Array[Array[Double]], Wavelet, mode? : Mode) -> Array[Double] raise WaveletError

// 2D
pub fn dwt2(Matrix, Wavelet, mode? : Mode) -> Coeffs2 raise WaveletError
pub fn idwt2(Coeffs2, Wavelet, mode? : Mode) -> Matrix raise WaveletError

// 连续小波
pub fn cwt(Array[Double], Array[Double], Wavelet, sampling_period? : Double, precision? : Int) -> (Matrix, Array[Double]) raise WaveletError
pub fn Wavelet::wavefun(Self, level? : Int) -> Wavefun raise WaveletError

// 构造与查询
pub fn Wavelet::new(String) -> Wavelet raise WaveletError   // "db2" / "sym4" / "mexh" ...
pub fn wavelist(kind? : WaveletKind?) -> Array[String]
pub fn threshold(Array[Double], Double, mode? : ThresholdMode, substitute? : Double) -> Array[Double] raise WaveletError
```

强类型枚举 `Mode` / `WaveletKind` / `Symmetry` / `CoeffPart` / `ThresholdMode` 均为
`pub(all)`，下游可直接构造；错误统一为 `WaveletError` 子错误分类
（`UnknownWavelet` / `NotDiscrete` / `InvalidMode` / `InvalidLevel` / `ShapeMismatch` /
`NotSupported` 等）。完整接口见 [`pkg.generated.mbti`](pkg.generated.mbti)。

## 三个使用场景

1. **信号去噪 / 压缩**：`wavedec` 分解 → `threshold` 软阈值细节系数 → `waverec` 重构。
2. **时频分析**：`cwt(signal, scales, Wavelet::new("mexh"))` 得到尺度–时间系数矩阵与对应频率。
3. **CLI / 差分对拍**：`moon run cmd/main -- dwt db2 symmetric 1 2 3 4 5 6 7 8` 输出 JSON；
   `tools/diff_test.py` 驱动 CLI 与本地 PyWavelets 逐用例对拍。

```bash
moon run cmd/main -- list discrete
moon run cmd/main -- info db2
moon run cmd/main -- dwt db2 symmetric 1 2 3 4 5 6 7 8
moon run cmd/main -- wavedec db2 symmetric 2 1 2 3 4 5 6 7 8
moon run cmd/main -- cwt mexh 1,2,4,8 1 2 3 4 5 6 7 8
moon run cmd/main -- wavefun db2 5
moon run cmd/main -- threshold soft 2 -3 -2 -1 0 1 2 3
moon run cmd/main -- pad antisymmetric 4 4 1 2 3
```

### 安装与可运行示例

```bash
moon add PaperY0/wavelets@0.1.0     # 从 mooncakes 安装
moon run examples/denoise           # 1D 小波软阈值去噪（打印去噪前后 MSE）
moon run examples/cwt_timefreq      # CWT 时频分析（打印各尺度频率与幅值）
```

[`examples/`](examples/) 提供两个最小可运行样例（纯计算、只依赖公共 API）；
`cmd/main` CLI 覆盖全部变换并输出 JSON，供 `tools/diff_test.py` 差分对拍。

## 与上游的刻意差异

逐条见 [`docs/known-deviations.md`](docs/known-deviations.md)。要点：

- **复值连续小波族**（cgau/shan/fbsp/cmor）不在 v0.1；只移植 10 个实值连续小波。
- **超层分解报错而非警告**：`wavedec`/`wavedec2` 在 `level > dwt_max_level` 时
  raise `InvalidLevel`（MoonBit 无警告通道）；上游是 warning 后继续。
- `central_frequency` 仅连续小波；双正交离散小波的 `integrate_wavelet` raise `NotSupported`；
  `wavefun` 精度上限 20。
- **SWT、小波包、FFT 法 CWT、自定义滤波器组**为 roadmap，不在 v0.1。

数值怪癖则**刻意保留**（非偏差）：`dmey` 是 Meyer 的 FIR 近似、非严格正交（重构误差
~0.03，上游同样如此）；`waverec` 可能多返回样本；soft/garrote 在 `x=0,t=0` 得 NaN；
`dwt` 对 reflect/antireflect 要求 N≥2（与上游一致，pad/idwt 无此限制）。

## 测试策略

**43 个 MoonBit 测试全部通过**，三层来源：

1. **黄金差分**（`golden_data_test.mbt`，由 `tools/gen_golden.py` 从 PyWavelets 1.8.0
   烘焙）：**1771 个用例**覆盖 dwt/idwt/wavedec/updown/dwt2/wavedec2/cwt/wavefun/
   centrfrq/intwave/pad/threshold/maxlevel/coefflen/qmf/orthfilt，容差 `1e-9`
   （NaN 编码为 JSON `null`）。
2. **性质扫描**：105 个正交/双正交小波 × 9 模式 × 9 长度的**完美重构**；多层与 2D 往返；
   滤波器组正交性；CWT 线性；分支可加性；dmey 近似界。
3. **白盒**：upsample/conv/imod/extend_sample/FFT（对朴素 DFT）/Hermite/linspace。

另有 **CLI 端到端差分**：[`DIFF-REPORT.md`](DIFF-REPORT.md)——`tools/diff_test.py`
驱动 CLI 对本地 PyWavelets 跑 **1527 项检查，0 处不一致**。

```bash
moon check   # 0 error
moon test    # 43/43（默认 wasm 后端）
moon test --target wasm-gc
moon test --target js
```

## 非目标（明确不做）

复值连续小波、SWT/小波包、FFT 法 CWT、自定义滤波器组、GPU/SIMD 加速（roadmap）。
本库是**纯计算移植**，不做 plotting、不做文件 IO、不读运行时钟。

## 开发

工具链版本固定在 [`.moonbit-version`](.moonbit-version)：

```bash
moon check && moon info && moon fmt && moon test
python -m pip install --target tools/ref/pylibs -r tools/ref/requirements.txt
python tools/gen_wavelets.py   # 重新生成 wavelet_data.mbt（滤波器系数）
python tools/gen_golden.py     # 重新生成 golden_data_test.mbt + tools/golden.json
python tools/diff_test.py      # CLI 对 PyWavelets 差分，重写 DIFF-REPORT.md
# 发布到 mooncakes：先临时移走 tools/ref/pylibs（约 50MB 参考依赖，避免打包上传超限），
# 发布完成后移回。pylibs 已在 .gitignore 中，不进入版本库。
moon publish
# 生成器产物需再跑一次 moon fmt
```

> 生成/对拍脚本需要 Python 3 与本地安装的 `PyWavelets==1.8.0` / `numpy==2.5.5`
> （`pip --target tools/ref/pylibs`，不污染全局环境）。Windows 上命令名通常是 `python`。

## 数据与许可

- 滤波器系数、延拓模式语义、变换算法与黄金向量均来自
  [PyWavelets](https://github.com/PyWavelets/pywt)（**MIT**，(c) Filip Wasilewski 及
  PyWavelets Developers），对 1.8.0 验证。
- 本项目以 MIT 发布，[`LICENSE`](LICENSE) 保留上游版权声明。

## 当前状态（2026-11-26）

✅ 106 离散 + 10 实值连续小波，系数逐位一致 · ✅ 1D/2D DWT/IDWT、多层、CWT、wavefun、
阈值化、9 延拓模式 · ✅ 43/43 测试（含 1771 黄金差分用例）· ✅ CLI 差分 1527/1527 一致 ·
✅ `moon check` 0 error · ✅ 纯计算无 IO，wasm/wasm-gc/js 三后端 ·
✅ CI：check/build/test 三后端 + 接口门禁 + 上游差分 · ✅ 两个可运行示例（`examples/`）·
✅ 已发布 mooncakes：`moon add PaperY0/wavelets@0.1.0`。

---

<a id="papery0wavelets-en"></a>

## PaperY0/wavelets (EN)

A MoonBit port of [PyWavelets](https://github.com/PyWavelets/pywt) 1.8.0: 1D/2D
discrete wavelet transforms (DWT/IDWT, multilevel), continuous wavelet transform
(CWT, convolution method), scaling/wavelet functions, thresholding, and 9 signal
extension modes. **106 discrete + 10 real continuous wavelets**, filter
coefficients bit-identical to upstream. Pure computation, no IO; wasm / wasm-gc /
js targets.

```moonbit nocheck
let w = @wavelets.Wavelet::new("db2")
let (cA, cD) = @wavelets.dwt([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], w)
let rec = @wavelets.idwt(Some(cA), Some(cD), w)
let coeffs = @wavelets.wavedec(signal, w, level=Some(3))
let (coefs, freqs) = @wavelets.cwt(signal, [1.0, 2.0, 4.0], @wavelets.Wavelet::new("mexh"))
```

A full mooncakes registry scan (3983 modules) found zero wavelet/DWT/CWT modules;
the niche is empty. Verification: 43 MoonBit tests including 1771 golden
differential cases baked from PyWavelets 1.8.0, plus property sweeps (perfect
reconstruction across 105 wavelets × 9 modes × 9 lengths); `tools/diff_test.py`
drives the CLI against a live PyWavelets with 1527 checks and 0 mismatches.

Deliberate deviations (see `docs/known-deviations.md`): complex continuous
wavelet families, SWT, wavelet packets, FFT-based CWT and custom filter banks are
out of scope for v0.1; over-level decomposition raises instead of warning;
`central_frequency` is continuous-only. Numerical quirks (dmey's non-exact
orthogonality, `waverec` length, soft/garrote NaN at 0/0, the N≥2 reflect/
antireflect restriction) are preserved to match upstream.

Filter coefficients, mode semantics, algorithms and golden vectors are from
PyWavelets (MIT, (c) Filip Wasilewski and the PyWavelets Developers). Released
under MIT; LICENSE retains upstream attribution. Install: `moon add PaperY0/wavelets@0.1.0`.