# Amber 通用 MD 流程教程（蛋白/RNA–配体）

> **日常蛋白–小分子/肽/抗体请优先使用一键入口 [`bin/amber-auto`](docs/AUTO.md)。** 本教程是高级路径，覆盖 RNA、多模型 NMR、逐阶段手动控制和 `amber-kit` 的完整配置。

本教程按“准备 → PBS → 建系 → 平衡 → 生产 → 分析 → 排错”组织。所有命令都在项目根目录执行；`bin/amber-kit` 是统一入口，底层脚本在 `scripts/`。

> 先读 [`README.md`](README.md) 的“重要原则”。PBS 头的逐行解释在 [`pbs/README.md`](pbs/README.md)，本教程只做摘要。

---

## 0. 这套流程的边界

适合：

- 蛋白–小分子、RNA–小分子、apo 蛋白/RNA 的显式溶剂 MD；
- 多模型 NMR 体系的模型选择、配体姿态移植；
- 3–10 个分子的相对 MM/GBSA 排序（仅作 pilot，不作筛选器）。

不适合（需要额外准备/其他方法）：

- 对接、缺失 loop/侧链建模、修饰核苷酸、共价配体、金属配位；
- FEP/ABFE、伞形采样、GaMD、REMD；
- 多轨迹 MM/PBSA、熵计算；
- 蛋白–RNA 复合物（需自定义 `leap.in`，见 5.6）。

---

## 1. 目录、配置与检查

### 1.1 复制一份项目

```bash
cp -r ~/Documents/amber_general_workflow ~/Documents/my_project
cd ~/Documents/my_project
bin/amber-kit init
```

`init` 会创建：

```
config/config.sh          # 从 config.example.sh 复制
config/systems.tsv        # 从 systems.example.tsv 复制
inputs/ ligands/ systems/ md/ analysis/ logs/
```

### 1.2 编辑 `config/config.sh`

最少改这些：

```bash
PROJECT_NAME="my_project"

AMBER_SH="/opt/software/amber22/amber.sh"        # 你的 Amber 初始化脚本
PY_QC="/opt/software/amber_general_workflow/envs/qc/bin/python"   # 有 parmed
PY_ANALYSIS="/opt/software/amber_general_workflow/envs/analysis/bin/python"      # 有 MDAnalysis

SYSTEM_TYPE="rna-ligand"      # rna-ligand|protein-ligand|rna-apo|protein-apo
RECEPTOR_PDB="inputs/receptor.pdb"
RECEPTOR_MODEL=""             # 多模型时填模型号；也可在 systems.tsv 的 models 列给

RECEPTOR_FF="RNA.YIL"         # RNA.YIL|RNA.OL3|ff14SB|ff19SB
SOLVENT="explicit"
WATER="tip3p"                 # RNA.YIL 配 tip3p；ff19SB 配 opc
BOX="oct"; PAD=10.0
CATION="K+"; ANION="Cl-"; COUNTERION="cation"; SALT_PAIRS=0
RADII="mbondi2"; IGB=5        # igb=5↔mbondi2；igb=8↔mbondi3

PBS_QUEUE="batch"
PBS_PPN=8
PBS_GPUS=1                    # 单个 pmemd.cuda 只用 1 张
PBS_GPUS_PER_NODE=4           # 本机每节点 4 张；其他集群改这里
PBS_TWO_GPU_PPN=48            # 整节点 CPU 核数
PBS_WALLTIME="24:00:00"
```

### 1.3 编辑 `config/systems.tsv`

每个配体（或每个质子化状态）一行，制表符分隔。表头必须是第一行（可以保留上面的 `#` 注释）：

```tsv
label	ligand_pdb	charge	resname	models	replicas	seed_base	control	experimental	experimental_type	notes
lig1	inputs/lig1.pdb	2	LIG	7,11,10	1,2,3	20261101	0	7.9	ic50	A1AZM
lig2	inputs/lig2.pdb	3	LIG	5,16,1	1,2,3	20261201	0	18.9	ic50	A1AZL
lig3	inputs/lig3.pdb	1	LIG	15,14,9	1,2,3	20261301	0	100.4	ic50	MQC +1
lig3_p2	inputs/lig3_p2.pdb	2	LIG	15,14,9	1,2,3	20261401	1	NA	NA	MQC +2 control
```

规则：

- `models` 为空 → 单模型；非空 → `systems/<label>/model_<NN>/`，与 `replicas` 按顺序对应。
- `replicas` 为空 → 默认 `1`。
- `control=1` 只进排名表，不参与与实验值的相关/置换统计。
- `experimental_type`：`ic50`/`kd`（越小越强）或 `dg`（越负越强）；`NA` 表示无。
- **空列必须保留制表符**；`amber-kit` 用 unit separator 解析，不会像 `read` 那样吞掉空字段。

### 1.4 检查

```bash
bin/amber-kit check
bin/amber-kit check --strict   # 输入 PDB 不存在也报错
```

检查内容：Amber/Python 环境、`SYSTEM_TYPE`、力场/水配对、`igb/radii`、PBS 资源格式、`systems.tsv` 列。

---

## 2. 输入准备

### 2.1 受体 PDB

- 可以是**只有受体**，也可以是**受体+配体复合物**；`amber-kit build` 会在受体 PDB 含配体 resname 时用 `scripts/split_complex.py` 拆分。
- 多模型 NMR PDB：用 `RECEPTOR_MODEL` 或 `systems.tsv` 的 `models` 列选模型（1-based）。
- `pdb4amber` 会：
  - `--dry` 去水；
  - `--most-populous` 保留主要 altloc；
  - `--nohyd` 去氢，tleap 按模板补氢。
- 非标准 HETATM（辅因子、修饰残基、金属）默认报错；确认属于体系后用 `--allow-hetatm` 重跑 `build_system.sh`，或先手工处理。

### 2.2 配体 PDB

要求：

- **显式氢**；
- **单残基**；
- 单模型（脚本保留第一个模型）；
- 原子名唯一；如果要做姿态移植，原子名在各姿态间一致。

`param_ligand.sh` 会拒绝无氢、多残基的 PDB。质子化状态必须你自己确定（pH、assay 条件、deposited state 都可能不同）；脚本不猜。

### 2.3 从复合物拆分配体

```bash
python -I scripts/split_complex.py \
  --pdb inputs/complex.pdb --model 7 \
  --ligand-resname MQC --out-resname LIG \
  --receptor-out inputs/receptor.pdb \
  --ligand-out inputs/lig.pdb \
  --complex-ref-out inputs/complex_ref.pdb
```

配体被重命名为 `LIG`、chain L、resseq 1；坐标/原子名保留。`complex_ref.pdb` 只用于肉眼核对，不参与 tleap。

### 2.4 多模型选择（可选）

```bash
python -I scripts/select_models.py \
  --pdb inputs/9CPG_allmodels.pdb --rna \
  --core-residues 1-6,8-13,14-19,21-26 \
  --loop-residues 7,20 --terminal-residues 1,13,14,26 \
  --ligand-resname LIG --n-select 3 \
  --out inputs/selected_models.tsv
```

输出：每个模型的 `core_rmsd_to_medoid_A`、`min_lig_loop_A`、`min_lig_terminal_A`、`end_stacked`、`selected`。规则：

- `medoid` = 到所有模型平均 core RMSD 最小；
- `diverse` = 距已选集合最远的模型；
- `end_stacked=1` = 配体更靠近末端而非 loop（>15 Å 或 terminal < loop 且 terminal < 4.5 Å），选模型时尽量避开。

把选中的模型号填进 `systems.tsv` 的 `models` 列。

### 2.5 一个常见的 RNA 建系坑：5′-磷酸

如果 PDB 的第一个核苷酸带 `P/OP1/OP2`，而 tleap 把它识别成 5′ 端 `G5`/`A5` 模板，可能报：

```
FATAL: Atom .R<G5 1>.A<OP1 34> does not have a type.
```

处理：

- 去掉第一个残基的 `P/OP1/OP2`（5′-OH），让 tleap 用标准 5′ 端模板；或
- 用自定义 `leap.in` 加载带 5′-磷酸的端基模板。

本工具包的测试夹具 `tests/data/tiny_rna.pdb` 就是去掉 5′-磷酸的版本。

---

## 3. 配体参数化

### 3.1 单配体

```bash
source /opt/software/amber22/amber.sh
export AMBER_MD_PYTHON=/opt/software/amber_general_workflow/envs/qc/bin/python

scripts/param_ligand.sh \
  --pdb inputs/lig.pdb --charge 2 --resname LIG \
  --gaff gaff2 --outdir ligands/lig1/template
```

输出：

- `lig.mol2`、`lig.frcmod`；
- `ligand_params.json`（电荷、原子数、ATTN、来源）；
- `antechamber.log`、`sqm.*` 等审计文件。

检查：

- mol2 电荷和 = `--charge`（±0.01）；
- `parmchk2` 非空；
- ATTN 默认失败；确认化学合理后才 `--allow-attn`。

### 3.2 多姿态移植

同一配体的不同对接/NMR 姿态，建议用模板电荷/类型移植，避免每个姿态重新 AM1-BCC 导致电荷噪声：

```bash
scripts/param_ligand.sh \
  --pdb inputs/lig_pose07.pdb --charge 2 --resname LIG --gaff gaff2 \
  --outdir ligands/lig1/model_07 \
  --transplant ligands/lig1/template/lig.mol2 --onto inputs/lig_pose07.pdb
```

`transplant_lig_coords.py` 按原子名匹配；模板缺少目标原子会硬报错。

### 3.3 `amber-kit build` 会自动做这些

- 单模型：`param_ligand.sh` 写入 `ligands/<label>/template/`；
- 多模型：第一个模型做模板，其余模型 `--transplant` 到 `ligands/<label>/model_<NN>/`；
- 然后 `build_system.sh` 用对应的 `lig.mol2/frcmod` 建系。

---

## 4. 建系与 QC

### 4.1 力场/水配对

| 溶质 | 推荐基线 | 常见替代 |
|---|---|---|
| RNA | `RNA.YIL` + TIP3P | `RNA.OL3` |
| 蛋白 | `ff19SB` + OPC | `ff14SB` + TIP3P |
| 配体 | GAFF2 + AM1-BCC | GAFF + AM1-BCC；RESP |
| RNA+配体 | RNA.YIL + GAFF2 + TIP3P | — |
| 蛋白+配体 | ff19SB + GAFF2 + OPC | ff14SB + GAFF2 + TIP3P |

不要混用 GAFF/GAFF2；不要只改 `igb` 而不改 `radii`。

### 4.2 `build_system.sh` 做了什么

1. 检查受体 PDB（水、HETATM、多模型、配体 resname 冲突）；
2. `pdb4amber --dry --most-populous --nohyd`；
3. 渲染 `leap.in`（`templates/build/*.in`）；
4. 加载配体 `frcmod/mol2`，`combine`；
5. `solvateOct`/`solvateBox`，`addIons` 中和，可选 `addIonsRand` 盐；
6. `set default PBRadii` 后 `saveAmberParm`；
7. 解析 `tleap.out` 的 `Errors/Warnings/ATTN`，检查净电荷/radii；
8. 运行 `qc_system.py`，写 `build_qc.tsv`/`build_manifest.json`。

### 4.3 直接调用（示例）

RNA–配体：

```bash
scripts/build_system.sh \
  --type rna-ligand \
  --receptor inputs/rna.pdb \
  --ligand ligands/lig1/template/lig.mol2 \
  --frcmod ligands/lig1/template/lig.frcmod \
  --lig-resname LIG --ligand-pdb inputs/lig.pdb --gaff gaff2 \
  --ff RNA.YIL --solvent explicit --water tip3p \
  --box oct --pad 10.0 --cation K+ --counterion cation \
  --radii mbondi2 --outdir systems/lig1
```

蛋白–配体：

```bash
scripts/build_system.sh \
  --type protein-ligand \
  --receptor inputs/protein.pdb \
  --ligand ligands/lig1/template/lig.mol2 \
  --frcmod ligands/lig1/template/lig.frcmod \
  --lig-resname LIG --ligand-pdb inputs/lig.pdb --gaff gaff2 \
  --ff ff19SB --solvent explicit --water opc \
  --box oct --pad 10.0 --cation Na+ --counterion none \
  --radii mbondi2 --outdir systems/lig1
```

只渲染不运行：

```bash
scripts/build_system.sh ... --render-only --force
# 检查 systems/lig1/leap.in
```

### 4.4 QC 门槛

`build_qc.tsv` 至少检查：

- `atoms/residues` 合理；
- `total_charge` ≈ 0（中和后）；
- `lig_charge` = 你填的电荷；
- `lig_pose_rmsd_A` 小（同一姿态）；
- `lig_min_dist_A` > 2.5 Å（无严重重叠）；
- `radii` = `mbondi2`/`mbondi3`；
- `restraint_atoms` 不为 0，且不包含水/离子。

`tleap` 必须 `Errors = 0`；`ATTN` 默认失败。

### 4.5 自动推导 masks

```bash
python -I scripts/residue_masks.py \
  --prmtop systems/lig1/complex.prmtop \
  --ligand-resname LIG \
  --json systems/lig1/masks.json \
  --env systems/lig1/masks.env --force
```

示例 `masks.env`：

```bash
RESTRAINT_MASK=':1-27 & !@H='
RECEPTOR_MASK=':1-26'
LIGAND_MASK=':27'
LIGAND_RESNAME_MASK=':LIG'
LIGAND_RESNAME='LIG'
LIGAND_RESID='27'
SOLUTE_RESIDUES='1-27'
LAST_SOLUTE_RESID='27'
RNA_RESRANGE='1-26'
```

`amber-kit inputs/analyze/mmpbsa` 会自动 source 它；不要再手写 `:1-27`。

### 4.6 蛋白–RNA 复合物（超出默认流程）

`build_system.sh` 只支持单一受体类型。蛋白–RNA 复合物需要：

1. 自定义 `leap.in`：同时 `source leaprc.protein.ff19SB` 和 `source leaprc.RNA.YIL`，分别 `loadPdb` 后 `combine`；
2. 手动检查原子类型冲突、末端、离子；
3. 用 `qc_system.py` + `residue_masks.py` 做 QC 和 masks；
4. 其余 MD/PBS/分析步骤可复用。

---

## 5. PBS 前缀（摘要）

完整逐行解释见 [`pbs/README.md`](pbs/README.md)。

标准单 GPU 头：

```bash
#!/bin/bash
#PBS -N myproject_lig1_r1_equil
#PBS -q batch
#PBS -l nodes=1:ppn=8:gpus=1
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o /path/to/project/logs/
#PBS -m n
```

必须记住：

- `#PBS` 必须在可执行语句之前；
- `cd "$PBS_O_WORKDIR"`，并从副本目录 `qsub`；
- GPU 绑定用 `PBS_GPUFILE`，保留 Torque 已设的 `CUDA_VISIBLE_DEVICES`；
- 单 `pmemd.cuda` 只申请 1 张 GPU；多副本用 `multi-gpu` 模板；
- `-o` 目录必须存在；
- 100 ns 约 7–8 h；4×25 ns 可放 24 h；超时用 chunk 或缩短 chunk；
- 本工具包的 `run_equil/run_prod` 会按 `.rst7 + Total wall time + NSTEP` 续跑，不会从头重算。

其他集群（GPU 数量不同）：

```bash
PBS_GPUS=1                 # 每个作业仍 1 张
PBS_GPUS_PER_NODE=4        # 节点实际 GPU 数
PBS_TWO_GPU_PPN=48         # 整节点 CPU 核数
```

```bash
amber-kit submit multi-gpu --pair lig1:1,lig1:2,lig2:1,lig2:2 --dry-run
```

---

## 6. 生成 MD 输入

```bash
bin/amber-kit inputs --label lig1 --rep 1
# 或全部
bin/amber-kit inputs
```

每个 `md/<label>/r<rep>/` 会得到：

```
01_min_restraint.in   5000 SD + 5000 CG，溶质重原子 10 kcal/mol/Å²
02_min_free.in        10000 SD，无约束
03_heat.in            NVT 0→300 K，2 ns，溶质重原子约束
04a..04j.in           NPT 300 K/1 atm，约束 10→0，10×1 ns
05_prod_c01..c04.in   NPT 生产，4×25 ns
bench.in              0.2 ns benchmark
run_config.json       机器可读配置
stage_table.tsv       阶段表
README.md             手跑命令
```

关键：

- `restraintmask` 来自 `masks.env`；不含水/离子；
- `ntx/irest`：加热 `irest=0, ntx=1`；后续 `irest=1, ntx=5`；
- `nstlim = ns × 1000 / dt`；`ntwx = traj_interval_ps / dt`；
- 种子：`seed_base + (rep-1)`，每个 stage 再递增；记录在 `run_config.json`；
- 生产 chunk 有状态，严格 c01→c04。

先检查：

```bash
bin/amber-kit inputs --label lig1 --rep 1 --print-only
```

---

## 7. 平衡与生产

### 7.1 提交平衡

```bash
bin/amber-kit submit equil --label lig1 --rep 1 --dry-run
bin/amber-kit submit equil --label lig1 --rep 1
```

渲染出的 `run_equil.pbs` 会：

- min1/min2 用 `sander.MPI -np 8`；
- heat/04a..j 用 `pmemd.cuda`；
- 每个 stage 完成检查 `.rst7 + Total wall time + NSTEP>=nstlim`；
- 已完成 stage 跳过；未完成 stage 及其下游重算；
- 最后 `cp 04j.rst7 npt.rst7`。

### 7.2 提交生产

`npt.rst7` 出现后：

```bash
bin/amber-kit submit prod --label lig1 --rep 1
```

或单 chunk：

```bash
bin/amber-kit submit prod-chunk --label lig1 --rep 1 --chunk c02
```

生产链最后 `cp 05_prod_c04.rst7 prod_final.rst7`。

### 7.3 进度

```bash
bin/amber-kit status --queue
qstat -a
```

`status` 显示每个副本的最后 stage、NSTEP、ns/day、是否完成/缺失/错误。

### 7.4 平衡 QC（进生产前）

- 温度 300 ± 5 K；
- 压力 1 ± 0.1 atm（NPT）；
- 密度 1.0 ± 0.02 g/cm³（水；蛋白/RNA 体系略低正常）；
- 盒子体积稳定，无塌缩/爆炸；
- 势能有限，无 NaN；
- 溶质 RMSD 不单调漂移；
- 配体仍在结合位点（未解离）；
- RNA 碱基配对/堆积基本保留。

不通过时先修 setup（姿态、质子化、电荷、最小化、盒子），不要直接延长生产。

### 7.5 断点续跑

- 直接重新 `qsub` 同一个 `run_equil.pbs` / `run_prod.pbs`；
- 已完成的 stage 不会重算；
- 未完成的 stage 从上一 stage 的 `.rst7` 重跑；
- 不要手动删除已完成的 `.rst7`；
- 不要同时投两个同副本作业（`flock` 会拒绝第二个）。

---

## 8. 分析

### 8.1 Dry trajectory

```bash
bin/amber-kit analyze --label lig1 --rep 1 --force
```

内部：

```bash
python -I scripts/make_cpptraj_inputs.py \
  --top md/lig1/r1/complex.prmtop \
  --traj 'md/lig1/r1/05_prod_c*.nc' \
  --strip-mask ':WAT,K+,Cl-,Na+' \
  --anchor ':1-13' --fixed ':14-26,27' \
  --outdir analysis/cpptraj --prefix lig1_r1 \
  --nastruct --resrange 1-26 --run
```

要点：

- 先 `autoimage` 再 `strip`；`iwrap=0` 的轨迹必须 imaging；
- 双链 RNA 用 `anchor/fixed`，避免一条链相对另一条跳盒；
- dry topology 必须和 dry trajectory 原子一一对应；
- `nastruct` 只对 RNA，`resrange` 不含配体/离子。

### 8.2 结构指标

`analyze_replica.py` 输出 `analysis/python/<label>_r<rep>.npz` + `.json`：

- 受体 backbone RMSD（蛋白 N,CA,C,O；RNA 骨架）；
- 配体重原子 RMSD（受体叠合后）、配体内部 RMSD、Rg、COM 漂移；
- 受体 Rg、SASA（`cpptraj surf` LCPO，stride 可调）；
- 配体–受体最小距离、接触原子/残基、每残基接触占有率、buried fraction；
- 双向 H 键；
- loop–配体接触（`LOOP_RESIDUES`）。

```bash
python -I scripts/analyze_replica.py --help
```

### 8.3 配体 H 键

```bash
python -I scripts/analyze_ligand_hbonds.py \
  --prmtop analysis/cpptraj/lig1_r1_dry.prmtop \
  --traj analysis/cpptraj/lig1_r1_dry.nc \
  --ligand-resname LIG --stride 2 \
  --out analysis/python/lig1_r1_hb.npz --force
```

判据：H···A < 2.5 Å、D···A < 3.5 Å、D-H···A > 120°；供体由 H 与 N/O/S 的距离 <1.3 Å 判定。

### 8.4 RNA 专有

**nastruct 汇总**：

```bash
python -I scripts/parse_nastruct.py \
  --prefix analysis/cpptraj/lig1_r1 \
  --label lig1 --rep 1 \
  --out analysis/tables/nastruct_lig1_r1.csv \
  --pairs-out analysis/tables/nastruct_pairs_lig1_r1.csv
```

输出每个碱基对的 `bp_occupancy`、`hb_mean`、buckle/propeller/opening，以及整体 rise/twist/roll/tilt/groove。

**堆积**：

```bash
python -I scripts/analyze_stacking.py \
  --prmtop analysis/cpptraj/lig1_r1_dry.prmtop \
  --traj analysis/cpptraj/lig1_r1_dry.nc \
  --ligand-resname LIG \
  --ligand-ring-atoms 'ring1:C1,C2,C3,C4,C5,C6;ring2:C9,C10,C11,C12,C13,C14' \
  --target-residues 7,20,6,8,19,21 \
  --out analysis/python/lig1_r1_stack.npz --stride 10
```

输出 `dist/angle/overlap`；常用判据：质心距离 < 4.5 Å、夹角 < 30°、投影重叠 > 0。阈值要在看结果前定好。

### 8.5 MM/GBSA

```bash
bin/amber-kit mmpbsa --label lig1 --rep 1 --check-only   # 只拆拓扑/写输入
bin/amber-kit mmpbsa --label lig1 --rep 1 --force        # 运行
```

默认：

- 单轨迹；
- `igb=5`、`saltcon=0.12 M`、`mbondi2`；
- 最后 50% 轨迹，每 100 帧取一帧（10 ps 帧 → 1 ns）；
- `idecomp=1` 每残基分解；
- 从**溶剂化拓扑** `ante-MMPBSA.py -s ':WAT,K+,Cl-,Na+' -n ':LIG' --radii mbondi2` 拆分。

输出：

```
analysis/mmpbsa/<label>_r<rep>/mmpbsa.dat
analysis/mmpbsa/<label>_r<rep>/mmpbsa.csv
analysis/mmpbsa/<label>_r<rep>/FINAL_DECOMP_MMPBSA.dat
analysis/mmpbsa/<label>_r<rep>/summary.tsv
```

注意：

- 不要用 dry prmtop 直接跑 MMPBSA；radii 不一致；
- `igb=8` 必须 `mbondi3` 拓扑；只改输入里的 `igb` 是无效的；
- MM/GBSA 没有熵，不是 Kd；只用于同一协议内的相对比较。

### 8.6 汇总与排序

```bash
bin/amber-kit report
```

输出：

- `analysis/tables/per_replica_metrics.csv`；
- `analysis/tables/system_summary.csv`（跨副本 mean/SD/SEM）；
- `analysis/tables/ranking.csv`（`rank_all` 含对照；`rank_calibration` 只含实验体系；bootstrap `p_rank1`；Spearman/置换 p）。

解读：

- 每个副本是一个独立样本；不要拿帧当独立样本；
- 少于 5 个分子不要宣称泛化；
- 3 个分子完美排序的随机概率 = 1/6；
- 如果 MD 排序与实验不符，先查 setup/采样，不要调介电常数/盐/窗口去凑排序。

---

## 9. 两个完整示例

### 9.1 蛋白–配体（ff19SB/OPC/Na+）

`config/config.sh`：

```bash
SYSTEM_TYPE="protein-ligand"
RECEPTOR_FF="ff19SB"
WATER="opc"
CATION="Na+"
COUNTERION="none"       # 体系净电荷为 0 时；否则 cation/anion
RADII="mbondi2"; IGB=5
```

`config/systems.tsv`：

```tsv
label	ligand_pdb	charge	resname	models	replicas	seed_base	control	experimental	experimental_type	notes
lig1	inputs/lig1.pdb	0	LIG		1,2,3	20261101	0	NA	NA	example
```

命令：

```bash
cp inputs/protein.pdb inputs/receptor.pdb
bin/amber-kit check
bin/amber-kit build
bin/amber-kit inputs
bin/amber-kit submit equil --label lig1 --rep 1 --dry-run
bin/amber-kit submit equil --label lig1 --rep 1
# npt.rst7 后
bin/amber-kit submit prod --label lig1 --rep 1
bin/amber-kit analyze --label lig1 --rep 1
bin/amber-kit mmpbsa --label lig1 --rep 1
bin/amber-kit report
```

### 9.2 RNA–配体（RNA.YIL/TIP3P/K+）

`config/config.sh`：

```bash
SYSTEM_TYPE="rna-ligand"
RECEPTOR_FF="RNA.YIL"
WATER="tip3p"
CATION="K+"
COUNTERION="cation"
BOX="oct"; PAD=10.0
RADII="mbondi2"; IGB=5
ANCHOR_MASK=":1-13"
FIXED_MASK=":14-26,27"
RNA_RESRANGE="1-26"
LOOP_RESIDUES="7,20"
LIGAND_RING_ATOMS="ring1:C1,C2,C3,C4,C5,C6;ring2:C9,C10,C11,C12,C13,C14"
STACK_TARGETS="7,20,6,8,19,21"
```

`config/systems.tsv`：

```tsv
label	ligand_pdb	charge	resname	models	replicas	seed_base	control	experimental	experimental_type	notes
lig1	inputs/lig1.pdb	2	LIG	7,11,10	1,2,3	20261101	0	7.9	ic50	example
```

命令同 9.1；`build` 会自动推导 `:1-N & !@H=` 的 restraint mask（26 nt + 配体）。

---

## 10. 排错与检查清单

### 10.1 建系

| 症状 | 处理 |
|---|---|
| `Could not find bond parameter` | 检查 `frcmod` 是否加载、GAFF/GAFF2 是否一致、配体原子类型 |
| `ATTN, need revision` | 看具体原子；常见于错误质子化/电荷、缺失参数；确认后再 `--allow-attn` |
| `addIons` 后净电荷不为 0 | 检查 counterion 方向（正体系加 `Cl-`，负体系加 `K+`/`Na+`）、配体电荷 |
| `G5 ... OP1 does not have a type` | 5′-磷酸问题；去掉端基 `P/OP1/OP2` 或用自定义端基模板 |
| 盒子太小/PME 错误 | 增大 `PAD`（10–12 Å）、改用 `rect`、减小 `cut` |
| `tleap Errors=0` 但拓扑不对 | 看 Warnings/ATTN、`build_qc.tsv`、`RADIUS_SET`、水/离子数、配体 resid |

### 10.2 MD

| 症状 | 处理 |
|---|---|
| 启动 NaN/vlimit | 姿态重叠、电荷/参数错、约束 mask 错、dt 太大；多最小化/修姿态 |
| `ntx=5` 报错 | restart 无速度；加热从 `irest=0, ntx=1`，之后才 `ntx=5` |
| GPU PME box 错误 | 缩短 release stage/增加段数；或 CPU 先平衡密度 |
| 生产 chunk 顺序错 | 严格 c01→c04；用 `status` 检查 NSTEP/时间 |
| 轨迹分子跳盒 | `iwrap=0` 正常；分析前 `autoimage` |
| 作业超时 | 用 `prod-chunk`；减小 `PROD_CHUNK_NS`；续跑 |

### 10.3 分析

| 症状 | 处理 |
|---|---|
| cpptraj mask 选错 | 含空格/`&` 的 mask 要加引号；用 `masks.env` |
| `parm/traj atom mismatch` | dry/solvated 拓扑与轨迹不匹配；用 `make_cpptraj_inputs.py` 生成的 dry 对 |
| `nastruct` 报错 | `resrange` 只覆盖 RNA；不要含配体/离子；stride 抽帧 |
| MMPBSA 原子数不匹配 | 用溶剂化拓扑 + `ante-MMPBSA.py`；dry 轨迹对 dry complex 拓扑 |
| MMPBSA 正 ΔG/异常 | 检查质子化/电荷、配体是否解离、radii/igb 是否匹配 |
| 只有 1 帧 | `startframe/endframe/interval` 1-based；检查轨迹帧数 |

### 10.4 PBS

| 症状 | 处理 |
|---|---|
| 作业不在队列 | 看 `qsub` 输出；`#PBS` 是否在可执行语句前 |
| 目录不对 | 从副本目录 `qsub`；`PBS_O_WORKDIR` 是提交目录 |
| 日志找不到 | `-o` 目录必须存在；`qstat -f` 看 `Output_Path` |
| GPU 不可见 | 在计算节点跑；检查 `PBS_GPUFILE`/`CUDA_VISIBLE_DEVICES`；`nvidia-smi -L` |
| 双副本互相覆盖 | 两个 `REP_DIR` 相同；模板会拒绝 |
| `qstat` 一直 `C` | 已完成；`qdel -p <id>` 清理 |

### 10.5 发表前检查清单

- [ ] 力场/水/离子/盒子/配体电荷/质子化都记录；
- [ ] `tleap Errors=0`，ATTN 已解释；
- [ ] `igb/radii` 匹配；
- [ ] 副本数 ≥3，种子记录；
- [ ] 平衡 QC 通过；
- [ ] 分析窗口/阈值在看结果前定好；
- [ ] MM/GBSA 用溶剂化拓扑，帧窗口/interval 记录；
- [ ] 报告 per-replica mean/SD，不只报一个数；
- [ ] n<5 时不宣称泛化；
- [ ] 没有为了排序去调参数。
