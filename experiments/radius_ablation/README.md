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
