# Radius mechanism ablation

固定矩阵：三个 `data/*.ply` 输入 × 五个模式；每组参数与任务要求一致，不启用 Chen。

在仓库根目录运行：

```powershell
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=E:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DRAR_BUILD_FIELD_VIEWER=OFF -DRAR_BUILD_ABLATION_VALIDATOR=ON
cmake --build build --config Release --target rar_cgal rar_radius_ablation_test rar_smoke_test rar_chen_continuation_test rar_validate_volume -j 2
ctest --test-dir build -C Release --output-on-failure
python experiments/radius_ablation/run.py
python experiments/radius_ablation/analyze.py experiments/radius_ablation/results/<时间戳>
python experiments/radius_ablation/validate_volumes.py experiments/radius_ablation/results/<时间戳> build/Release/rar_validate_volume.exe
```

`run.py` 默认使用新的时间戳目录，也可指定 `--exe` 和 `--results`。结果目录已存在时拒绝覆盖。
运行前要求跟踪的源文件已提交；已有未跟踪输入允许保留，其状态及输入 SHA-256 会写入清单。
其他机器需要按实际位置设置 CMake toolchain。Python 依赖及本次版本见结果目录的 `analysis_environment.json`。
`mesh_metrics.py` 复用了现有 `tools/compare_wine.py` 的质量和 VTK locator 函数，避免导入原脚本的顶层批处理。

每组保留 `command.json`、`command.txt`、`stdout.log`、`stderr.log`、`status.json`、`output.ply`、`initial.csv`、`metrics.json` 及采样距离数组。
初始 CSV 在 remeshing 之前由输入网格副本导出；分析脚本核验顶点坐标与原始输入一致，以及五组初始诊断完全一致。
`--no-export-field` 仅关闭原有常规场导出；本实验的 `--radius-diagnostics` 独立启用。
`runtime_seconds` 是程序内部的场构造及 remeshing 计时，不含独立诊断、输入和输出 I/O；`wall_seconds` 包含完整进程。

误差参考为相同的输入 PLY，不是额外的解析曲面或 ground truth。
每个方向 80,000 个面积均匀随机样本；参考种子 731，输出种子 343。
距离由 VTK 最近三角形查询获得；CD² 为两个方向平均平方距离的均值，另存未平方平均距离。
P99 为双向合并样本的百分位，同时保存各方向 P99。最大距离仅是最大采样距离，不能称为精确 Hausdorff。

`validate_volumes.py` 额外检查输入和输出是否闭合、有无自交、是否满足 CGAL `does_bound_a_volume`，并检查有向体积为正。
仅在两端均通过时填入 `(V_output - V_input) / V_input`。闭合但自交的网格仍将相对体积变化标为缺失；`vol` 仅保留代数有向体积，不解释为有效实体体积。

局部区域用初始顶点的 Cap/Recovery/Unchanged 标签定义，在三角形内做重心线性插值。
区域面积严格等于对应顶点 `vertex_area` 之和。输出面重心映射到输入的最近三角形，在最近点插值三个区域权重；每个输出面按权重分摊为“等效面数”。
密度为等效面数除以固定输入区域面积，并与 Pmap 对比。该最近面映射不约束拓扑对应，薄结构附近可能跨层；记录反向法线比例和映射距离 P99 用于诊断，但这些指标不能排除同向近邻层的错配。

新增 CLI 模式只选择初始尺寸：Pmap=`h_raw`，Cap=`min(h_raw,h_radius)`，Recovery=`max(...)`。
所有新模式直接调用 CGAL `at(v,mesh)` 获取 raw；完整 Radius 和原生 CGAL 的旧实现未改动。
Cap/Recovery 的数值场使用严格 min/max；诊断标签以 `64*机器epsilon*max(h_raw,h_radius)` 排除算术舍入噪声。
`clipped_unchanged` 要求 Guard 分支实际触发、两种尺寸严格相等且等于配置的上下界之一。
非有限 raw 会明确报错；不会擅自把 CGAL 的 NaN 改成其他尺寸。Radius 保留旧代码对非有限曲率返回 max-edge 的语义。

可选旧版本回归：在独立 worktree 构建原提交，再运行：

```powershell
python experiments/radius_ablation/check_legacy.py experiments/radius_ablation/results/<时间戳> build/radius-baseline/Release/rar_cgal.exe
```

原始大文件保存在本地 `results/`（已忽略）；紧凑汇总、报告和脚本保留在 Git 中。


# 任务：Radius Guard 的机制分离与消融实验

## 研究目的

我们正在研究 CGAL Adaptive Remeshing 的 Radius Guard，想验证它是否分别解决了两类初始尺寸异常：

- **H1：异常细化恢复（Recovery）**。当 `t = epsilon * curvature` 较大时，原始曲率—尺寸映射可能给出异常小的目标尺寸，产生过度细化。Radius Guard 可能对此有所缓解。
- **H2：过大尺寸限制（Cap）**。某些高曲率区域原始目标尺寸大于局部曲率半径，Radius Guard 可能限制其尺寸，从而改善局部几何逼近。

此前 Wine Glass、Camel2、Shark13 三个模型的结果不同。我们希望通过机制分离，判断上述假设是否成立。

**请直接在当前仓库实现、编译、测试并运行实验。不要只给出实现方案。**

## 1. 首先检查现有实现

仓库：`RedHe-rring/RAR_CGAL`

基于当前 `exp/chen-element-count-objective` 分支。

主要检查：

- `include/rar/CGALAdaptiveSizing.h`
- `include/rar/CGALAdaptiveRemesher.h`
- `include/rar/PropertyMapSizingField.h`
- `include/rar/RemeshConfig.h`
- `src/main.cpp`
- `tests/smoke.cpp`

先检查 `git status`，不要覆盖未提交修改。建议在安全的前提下创建独立实验分支 `exp/radius-mechanism-ablation`。

## 2. 实现机制分离

对输入网格的每个顶点，计算：

- `h_raw`：由 CGAL `PMP::Adaptive_sizing_field` 提供的原始目标尺寸，优先直接调用其 `at(v, mesh)`，而不是重复实现其内部公式。
- `h_radius`：由现有 `cgal_adaptive_radius_target_length()` 得到的目标尺寸。
- `kappa`：CGAL 的最大绝对主曲率。
- `t = epsilon * kappa`。

在这两个目标尺寸基础上定义：

**Cap-only：**

`h_cap(v) = min(h_raw(v), h_radius(v))`

它只保留 Radius Guard 使目标尺寸减小的部分。

**Recovery-only：**

`h_recovery(v) = max(h_raw(v), h_radius(v))`

它只保留 Radius Guard 使目标尺寸增大的部分。

**Full Radius：**

`h_full(v) = h_radius(v)`

注意这是一种以实际尺寸变化方向定义的机制分离，而不是仅仅根据理论阈值划分。

所有目标尺寸应使用同样的 min/max-edge 设置。对于零曲率、非有限曲率、数值裁剪和阈值附近的情况，需要按照现有代码的语义正确处理。

增加以下 CLI 模式：

- `cgal-adaptive-pmap`：CGAL 原始顶点尺寸 + PropertyMapSizingField
- `cgal-adaptive-radius-cap`：Cap-only
- `cgal-adaptive-radius-recovery`：Recovery-only

保留现有的：

- `cgal-adaptive`
- `cgal-adaptive-radius`

新模式应与现有 Radius 模式使用同一个 `PropertyMapSizingField` 和相同的 CGAL Remeshing 参数。

特别注意，原始 `cgal-adaptive` 直接使用 CGAL 自带的 Sizing Field，而 Radius 使用 `PropertyMapSizingField`。必须保留 `cgal-adaptive-pmap` 作为控制组，避免将 Adapter 差异误判为 Radius 公式的效果。

不要修改已有 Chen 优化代码，也不要重写 Split / Collapse / Flip / Relaxation / Projection。

## 3. 增加初始场诊断导出

在 Remeshing 开始之前，使用原始输入网格导出 CSV。

至少包含以下字段：

`vertex_id, x, y, z, vertex_area, kappa, epsilon_kappa, h_raw, h_radius, h_cap, h_recovery, delta, ratio, mechanism`

其中：

- `vertex_area` 为关联三角形面积的三分之一之和。
- `delta = h_radius - h_raw`
- `ratio = h_radius / h_raw`
- `mechanism` 根据实际尺寸变化标记为 Cap、Recovery 或 Unchanged。

确保 CSV 的数据来自 Remeshing 前的原始网格，不能使用已经变化的输出网格计算。

统计以下指标：

1. Radius 触发面积比例：`t > 0.183`
2. 实际发生 Cap 的面积比例
3. 实际发生 Recovery 的面积比例
4. `t > 1.8165` 和 `t > 2` 的面积比例
5. `h_raw` 与 `h_radius` 的分布、修改幅度和极值
6. 修改前后的元素数量代理量，例如 `sum(A_i / h_i^2)`

需要额外记录由于 min/max-edge 裁剪而触发分支但目标尺寸没有发生变化的情况。

所有面积比例应使用输入网格面积权重，而不是单纯的顶点数量比例。

## 4. 先做单元测试

在 `tests/smoke.cpp` 中或者新建测试文件，验证：

- `min(h_raw, h_radius)` 不大于 `h_raw`
- `max(h_raw, h_radius)` 不小于 `h_raw`
- `h_cap` 和 `h_recovery` 分别只保留一个方向的变化
- `h_raw` 和 `h_radius` 相等时，两种机制均不改变尺寸
- 处理零曲率、非有限曲率以及 min/max-edge 裁剪
- 当 `t > 2` 时，正确复现 CGAL 原始目标尺寸与 Radius 目标尺寸的实际差异

同时检查：现有 `cgal-adaptive` 与新的 `cgal-adaptive-pmap` 即便顶点目标尺寸相同，最终 Remeshing 行为是否仍有差异。如有差异，必须在报告中说明。

编译并运行 CTest，记录结果。

## 5. 执行批量实验

使用以下三个输入模型：

| Model | epsilon |
|---|---|
| wine_glass24.ply | 0.0214378565 |
| camel2.ply | 0.001 |
| shark13.ply | 0.003106642746 |

先定位仓库中的真实模型路径。如某个输入缺失，请明确指出，不要以其他模型替代。

每个模型运行五种模式：

1. cgal-adaptive
2. cgal-adaptive-pmap
3. cgal-adaptive-radius-cap
4. cgal-adaptive-radius-recovery
5. cgal-adaptive-radius

固定参数：

- `--min-edge 0.001`
- `--max-edge 0.5`
- `--iterations 5`
- `--relax-steps 3`
- `--no-preserve-features`
- Projection 开启
- 使用上述明确指定的 epsilon

本阶段不启用 Chen。

在 `experiments/radius_ablation/` 下创建实验脚本和结果目录。

保存完整命令、标准输出、错误日志、运行状态、输出 PLY 和初始场诊断。

如果运行失败，记录明确的失败原因，不要静默忽略，也不要自动更换参数。

## 6. 计算实验指标

对于每个输出网格，统计：

- Vertex / Face Count
- MinAngle P05
- Tri.Q P05
- Aspect Ratio P95
- Watertight / Component Count
- Runtime

如果原始参考网格可用，继续计算：

- 双向采样 Chamfer Distance，明确说明是否平方
- 双向采样 P99 距离
- 最大采样距离（注明不是精确 Hausdorff）
- Relative Volume Change（仅在体积定义有效时）

几何误差比较需要使用相同的采样数量、固定随机种子、统一距离定义。

建议每个方向至少 80,000 个采样点，并使用三角形表面最近距离，而不是只计算顶点到顶点的距离。

检查分析脚本已有的依赖和指标实现，尽量复用；如依赖缺失，不要伪造结果。

## 7. 重点分析两个假设

**Wine Glass：**

检验 Recovery-only 是否能够解释 Full Radius 的大量面数减少。

**Shark13：**

检验 Cap-only 是否能够解释 Full Radius 的误差或质量改善。

**Camel2：**

检验 Cap-only 与 Recovery-only 是否都只有很小影响。

此外，不能只根据全局指标判断。请尽可能生成：

- `epsilon * kappa` 的分布直方图
- `h_radius / h_raw` 的分布直方图
- 初始网格的 Cap/Recovery 区域可视化
- 各区域的局部三角形密度变化

对于局部对应关系，如果需要把输出面映射到原始网格，请说明映射算法及其局限，特别注意靠近薄结构的情况。

## 8. 最终输出

生成：

`experiments/radius_ablation/REPORT.md`

报告至少回答：

1. 原始 CGAL 和 PropertyMapSizingField 的差异是否明显？
2. Wine Glass 的主要变化来自 Cap 还是 Recovery？
3. Shark13 的主要变化来自 Cap 还是 Recovery？
4. 为什么 Camel2 的变化较小？
5. 哪些现象得到证实，哪些假设被否定？
6. Radius Guard 是否仍然有需要改进的地方？

同时输出：

- 每个模型的一张五方法指标对比表
- 所有模型的汇总 CSV
- 初始场诊断 CSV
- 必要的图片
- 可复现批处理脚本

**验收标准：**

- 现有模式的行为保持不变。
- 新模式通过编译和单元测试。
- 每个新模式只改变指定的初始尺寸场。
- 三个模型使用相同的对照参数。
- 实验结果可由脚本重现。
- 最终报告明确区分事实、机制假设与尚未验证的结论。

不要为了让假设成立而修改数据或调整实验参数。

完成后给出修改的文件清单、测试结果、实验结果、主要发现和未解决问题。不要自动推送 GitHub，也不要删除已有实验结果。