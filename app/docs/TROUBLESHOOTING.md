# 排错手册

先看日志，不要先改参数：

- 建系：`tleap.out`、`pdb4amber.log`、`antechamber.log`、`sqm.out`、`parmchk2.log`
- MD：`<stage>.out`、`<stage>.mdinfo`、PBS `.o<jobid>`
- cpptraj：`*_prep.log`、`*_nastruct.log`
- MMPBSA：`ante.log`、`mmpbsa.log`、`mmpbsa.dat`

---

## 1. 建系 / tleap

### `Could not find bond parameter for ...`

- `lig.frcmod` 没加载，或加载了错的；
- GAFF/GAFF2 混用：`antechamber -at gaff2` 但 `parmchk2 -s gaff` 或 `source leaprc.gaff`；
- mol2 原子类型不在 GAFF 库；
- 非标准残基（修饰核苷酸、辅因子）需要自己的参数。

排查：

```bash
grep -n "loadAmberParams\|loadMol2\|source leaprc" systems/<label>/leap.in
grep -n "Could not find" systems/<label>/tleap.out
head -30 ligands/<label>/template/lig.frcmod
```

### `ATTN, need revision`

- 看 tleap 输出里的原子/残基；
- 常见：错误质子化/电荷、缺失 frcmod 条目、异常价态、金属/辅因子被当作标准残基；
- 不要盲目 `--allow-attn`；确认化学合理并记录原因。

### `addIons` 后净电荷不为 0

- 正体系需要 `Cl-`（`--counterion anion`），负体系需要 `K+`/`Na+`（`--counterion cation`）；
- 配体电荷错；
- 体系已含离子，`addIons ... 0` 可能加错数量；
- 看 `tleap.out` 的 `Total unperturbed charge:` 和 `build_qc.tsv` 的 `total_charge`。

### `solvateOct`/`solvateBox` 失败或盒子太小

- 增大 `PAD`（10–12 Å）；
- 细长溶质改 `BOX=rect`；
- 盒子小于两倍 cutoff 会导致 GPU PME 错误。

### `tleap Errors = 0` 但拓扑不对

- 看 Warnings/ATTN，不只看 Errors；
- 检查 `RADIUS_SET`、原子/残基计数、水模型、离子数、配体 resid；
- 建系成功 ≠ 科学正确：质子化、姿态、电荷都可能错。

### `G5 ... OP1/OP2/P does not have a type`

- 5′ 端带磷酸，但 `G5` 模板没有 `P/OP1/OP2`；
- 去掉第一个残基的 `P/OP1/OP2`（5′-OH），或自定义端基模板。

### 蛋白二硫键 / 非标准残基

- `pdb4amber` 会写 `*_sslink`；tleap 需要 `bond` 命令；
- 用 `--extra-leap` 文件：

```leap
bond com.26.SG com.84.SG
```

（holo 单位是 `com`，apo 是 `mol`。）

### pdb4amber 去掉了想要的结构氢

- 默认 `--nohyd`；如果配体/辅因子的结构氢重要，先手工处理，或不要用 `--nohyd`（需要改 `build_system.sh` 调用）；
- 配体 PDB 必须有显式氢，`param_ligand.sh` 会拒绝无氢输入。

---

## 2. 配体参数化

### `antechamber` 在 `sqm` 阶段失败

- 错误 `-nc`（净电荷）；
- 缺氢/质子化不对；
- 几何太差；
- 看 `sqm.out` 最后的 SCF/几何错误。

### mol2 电荷和 ≠ 请求电荷

- 质子化状态与形式电荷不一致；
- 不要用 `-nc` 硬凑；重新确定质子化。

### `parmchk2` 输出为空/很短

- mol2 没有可用原子类型；
- 检查 `-at gaff2` 和 mol2 的 `@<TRIPOS>ATOM` 类型列。

### 移植缺原子

- 姿态 PDB 原子名和模板 mol2 不一致；
- 对齐命名，或对该姿态单独参数化。

---

## 3. MD / GPU

### 启动 NaN / vlimit

- 初始几何重叠：多最小化、修姿态；
- 电荷/参数错：看 `build_qc.tsv`、`parmchk2`；
- 约束 mask 错（约束到水/离子）；
- `dt` 太大（没有 HMR 时不要超过 0.002 ps）；
- 盒子太小。

### GPU PME box 错误

`Periodic box dimensions have changed too much`

- 用更短的 NPT release stage（`RELEASE_NS=0.1`、`RELEASE_STAGES=100`）；
- 或先用 CPU `sander.MPI` 平衡密度；
- 增大 `PAD`、减小 `cut`。

### `gpu_neighbor_list_setup :: Small box detected`

`pmemd.cuda` 在盒子某个维度只有 ≤2 个 neighbor-list cell 时拒绝运行（cut=10 Å 时约对应
box < 35–40 Å）。只含几个残基的测试体系容易触发；正常蛋白体系 `PAD=10 Å` 一般不会。

- 增大 `PAD`（如 15–20 Å）重建体系；
- 或确认后用 `-AllowSmallBox`（Amber 明确不推荐，仅限测试）。

### `mpirun` / `sander.MPI` 挂起（min 阶段无输出）

继承的 `LD_LIBRARY_PATH` 里若混有另一套 OpenMPI（例如之前 source 过别的 `amber.sh`、
conda 的 openmpi，或系统 OpenMPI），`mpirun` 会加载版本不一致的 `libmpi`，作业卡在
`01_min_restraint` 且不写 `.out`。

PBS 模板现在会在 source `amber.sh` 前删除 `LD_LIBRARY_PATH` 中的 `openmpi` 路径。
手工运行时也应这样做，或在一个干净环境里 source 目标 `amber.sh`：

```bash
LD_LIBRARY_PATH=$(printf '%s' "${LD_LIBRARY_PATH:-}" | tr ':' '\n' \
  | grep -v -i openmpi | paste -sd: - || true)
export LD_LIBRARY_PATH
source /path/to/amber.sh
```

用 `ldd "$AMBERHOME/bin/sander.MPI" | grep mpi` 检查所有 MPI 库是否来自同一个 OpenMPI。

### `ntx=5` restart 错误

- `irest=1, ntx=5` 需要 restart 里有速度；
- 最小化 restart 没有速度：加热用 `irest=0, ntx=1`；
- restart 截断/损坏：看文件大小和上一 stage `.out`。

### 生产 chunk 顺序

- chunk N 读 chunk N−1 的 `.rst7`；
- 不要并发/乱序；
- 用 `bin/amber-kit status` 看 NSTEP/时间；
- 本工具包的 `run_prod.pbs` 会跳过已完成 chunk，未完成 chunk 及其下游重算。

### 轨迹分子跳盒

- `iwrap=0` 未 imaging 是正常的；
- cpptraj `autoimage` 后再分析/导出；
- 双链 RNA 用 `anchor/fixed`。

---

## 4. cpptraj

### mask 选错

- 含空格/`&` 的 mask 必须加引号：

```cpptraj
strip ':1-27 & !@H=' parmout dry.prmtop
```

- 用 `masks.env` 的 mask。

### `parm`/`trajin` 原子数不匹配

- 拓扑和轨迹必须原子一一对应；
- 溶剂化拓扑配 dry 轨迹会失败；
- 用 `make_cpptraj_inputs.py` 生成的 dry 拓扑/轨迹对。

### `nastruct` 报错

- `resrange` 只覆盖 RNA，不含配体/离子；
- 修饰核苷酸需要特殊处理；
- 用 stride 抽帧。

---

## 5. MMPBSA

### `ante-MMPBSA.py` FileExists

- 脚本默认拒绝覆盖；用 `--force`（只清它管理的输出）或删旧目录。

### 原子数不匹配 / `cpptraj cannot pair`

- 轨迹是 dry，拓扑是 solvated（或反过来）；
- 用溶剂化拓扑 + `ante-MMPBSA.py` 生成与 dry 轨迹匹配的 complex/receptor/ligand 拓扑；
- ligand mask 在 strip 后 resid 可能变；优先用 `:LIG`。

### `RADIUS_SET` 不匹配 / 正 ΔG

- `igb=8` 需要 `mbondi3`；`igb=5` 需要 `mbondi2`；
- 不要只改输入 `igb`；要重建拓扑；
- 大正 ΔG 常见于错误质子化/电荷、配体解离、radii 不一致；
- MM/PBSA 是相对 proxy，不是 Kd。

### 只选了 1 帧

- `startframe/endframe/interval` 是 1-based；
- 检查帧数和 interval；`mmpbsa_one.sh` 会在 <2 帧时警告。

---

## 6. Python / 环境

- `amber.sh` 可能把 Amber 自带 `parmed` 放到 `PYTHONPATH`；`python3 -I` 会忽略 `PYTHONPATH`。
- 本机建议：
  - QC/建系：`/opt/software/amber_general_workflow/envs/qc/bin/python`（parmed 4.2）；
  - 分析：`/opt/software/amber_general_workflow/envs/analysis/bin/python`（MDAnalysis 2.7 / scipy 1.13）。
  - 仅 RNA 堆积分析需要 `shapely`；本机尚未安装（离线环境），蛋白体系不受影响。
- 如果 `import parmed` 失败：

```bash
/opt/software/amber_general_workflow/envs/qc/bin/python -I -c "import parmed, numpy"
```

- 如果 `import MDAnalysis` 失败：检查 `PY_ANALYSIS`。
- `analyze_replica.py` 在 `python -I` 下会显式把脚本目录加入 `sys.path`，以便 import `analyze_ligand_hbonds`。

---

## 7. PBS / 集群

### 作业不在队列

- `qsub` 是否报错；
- `#PBS` 是否在可执行语句前；
- 队列名/资源是否合法：`qstat -q`、`qstat -Q`。

### 目录不对

- `PBS_O_WORKDIR` 是提交目录，不是脚本目录；
- 从副本目录 `qsub`，或用 `qsub -d /path`。

### 日志找不到

- `#PBS -o` 目录必须存在且可写；
- `qstat -f <jobid>` 看 `Output_Path`/`Error_Path`；
- `#PBS -j oe` 合并到 `.o`。

### GPU 不可见

- 登录节点没有 GPU；必须在计算节点；
- 看 `PBS_GPUFILE`、`CUDA_VISIBLE_DEVICES`；
- `nvidia-smi -L` 在作业日志里；
- 其他集群的 `PBS_GPUFILE` 格式可能不同；模板的正则提取末尾数字，必要时改。

### 双副本互相覆盖

- 两个 `REP_DIR` 相同；`run_two_gpu`/`run_multi_gpu` 会拒绝；
- 检查 `multi-gpu` 的 `REP_DIRS`。

### `qstat` 一直显示 `C`

- 已完成作业保留；`qdel -p <id>` 清理。

---

## 8. 排序 / 统计

- 每个副本是一个样本，不是每一帧；
- 报告 per-replica mean/SD/SEM；
- 少于 5 个分子不要宣称泛化；
- 3 个分子完美排序 p=1/6；
- 排序不符时先查 setup/采样，不要调参数凑排序；
- 如果报告 `igb=8` 敏感性，必须用 `mbondi3` 拓扑重建；复用 `mbondi2` 的 `igb=8` 结果是无效的。
