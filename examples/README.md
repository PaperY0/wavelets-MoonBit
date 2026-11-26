# 可运行示例

两个最小使用样例，均只依赖本库公共 API，纯计算、可直接运行：

```bash
moon run examples/denoise        # 1D 小波软阈值去噪，打印去噪前后 MSE
moon run examples/cwt_timefreq   # CWT 时频分析，打印各尺度频率与最大幅值
```

另有 JSON CLI（`cmd/main`）覆盖全部变换，见根 README「使用场景」。

- `denoise/`：含噪信号 → `wavedec` → 细节支 `threshold(Soft)` → `waverec`，
  验证 MSE 下降。对应申报书场景 1（信号去噪）。
- `cwt_timefreq/`：chirp 信号 → `cwt(mexh, scales)` → 读取系数矩阵与频率向量。
  对应申报书场景 3（时频分析）。