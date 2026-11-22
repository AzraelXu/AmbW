# amber_general_workflow

一套**与具体项目无关**的 Amber 显式溶剂 MD 工具包：体系准备 → 平衡/生产 → 分析，重点覆盖 **PBS 作业前缀、GPU 绑定、分段与续跑**。适用于蛋白–小分子、蛋白–肽、蛋白–抗体、RNA–配体、apo 蛋白/RNA。

**推荐入口：`bin/amber-auto`。** 给一个或多个复合物 PDB，一条命令完成建系、QC、输入生成和 PBS 依赖链提交：

```bash
bin/amber-auto run complex1.pdb complex2.pdb
```

- 一键流程教程：[`docs/AUTO.md`](docs/AUTO.md)
- 错误与逃生口：[`docs/ERRORS.md`](docs/ERRORS.md)
- 高级 `amber-kit` 教程：[`TUTORIAL.md`](TUTORIAL.md)、[`docs/ADVANCED_AMBER_KIT.md`](docs/ADVANCED_AMBER_KIT.md)
- PBS 前缀逐行详解：[`pbs/README.md`](pbs/README.md)
- 命令速查：[`docs/CHEATSHEET.md`](docs/CHEATSHEET.md)
- 排错：[`docs/TROUBLESHOOTING.md`](docs/TROUBLESHOOTING.md)

> 本工具包不替代科学判断：配体净电荷/质子化、力场/水模型、分析窗口、MM/GBSA 参数都必须由你在看到结果前确定。`amber-auto` 会明确打印它推断的电荷、方法和 SMILES；模糊时直接报错，不静默猜测。

---

## 1. 能做什么

| 模块 | 内容 |
|---|---|
| 一键流程 | `amber-auto`：复合物分类 → 配体补氢/电荷推断 → Reduce 质子化 → tleap 建系/QC → masks/inputs → PBS 依赖链 → cpptraj/分析/MM-GBSA |
| 体系准备 | `pdb4amber` + `tleap`；RNA/蛋白/apo/holo；蛋白–肽/抗体；GAFF2 + AM1-BCC；多模型 NMR 拆分/移植；QC；自动推导 masks |
| 模拟 | min1/min2 → NVT 加热 → 10 段 NPT 释放 → 4×25 ns 生产；Torque PBS；GPU 绑定；断点续跑 |
| 分析 | dry trajectory；RMSD/Rg/SASA/接触/H 键；RNA nastruct/堆积；MM/GBSA + 分解；多副本汇总 |

**不做**：对接、缺失 loop 建模、N-糖链、共价配体、金属配位、FEP/ABFE、伞形采样、GaMD、REMD、多轨迹 MM/PBSA、熵计算。蛋白–RNA 复合物需要自定义 `leap.in`，请使用 `amber-kit`。

---

## 2. 环境

本机示例（其他集群请改 `config/config.sh`）：

```bash
source /opt/software/amber22/amber.sh
# QC/建系：parmed
/opt/software/amber_general_workflow/envs/qc/bin/python -c "import parmed, numpy"
# 分析：MDAnalysis/scipy（仅 RNA 堆积分析还缺 shapely，见 docs/TROUBLESHOOTING.md）
/opt/software/amber_general_workflow/envs/analysis/bin/python -c "import MDAnalysis, numpy, scipy"
```

PBS：Torque（`qsub/qstat/qdel`），队列 `batch`。其他集群见 [`pbs/README.md`](pbs/README.md) 第 6 节。

---

## 3. 使用流程

### 3.1 一键流程（推荐）

```bash
# 建系 + 生成输入 + 自动提交 equil -> prod -> post
bin/amber-auto run complex1.pdb complex2.pdb

# 只看分类、链角色、配体候选和电荷推断，不建系
bin/amber-auto inspect complex1.pdb

# 只建系/生成输入，不提交
bin/amber-auto build complex1.pdb

# 先检查 PBS 依赖链，不真正 qsub
bin/amber-auto run complex1.pdb --dry-run

# 进度 / 汇总报告 / 断点续跑
bin/amber-auto status --queue
bin/amber-auto report
bin/amber-auto resume
```

默认值：`ff19SB + OPC`、显式溶剂、100 ns（4×25 ns）、1 个副本、生产结束后自动跑 cpptraj 分析 + MM/GBSA。常用覆盖：

```bash
# 3 个副本、200 ns
bin/amber-auto run complex.pdb --replicas 3 --prod-ns 200

# 手动指定配体和电荷（自动推断失败或需要特定质子化态时）
bin/amber-auto run complex.pdb --ligand 0QE:A:501 --charge -2

# 显式指定肽/抗体链角色
bin/amber-auto run complex.pdb --type protein-peptide \
  --receptor-chains A --partner-chains B

# 隐式溶剂快速筛查
bin/amber-auto run complex.pdb --solvent implicit --prod-ns 20
```

输出在 `./amber-campaign/`：

```
campaign.json                        # 全局状态、job id、依赖、warnings
complexes/<label>/prepare/           # 清洗后的蛋白/配体、Reduce 日志、二硫键
complexes/<label>/systems/<label>/   # prmtop/inpcrd、masks、QC、linkage
complexes/<label>/md/<label>/r1/     # MD 输入、run_equil/prod/post.pbs
complexes/<label>/analysis/          # dry trajectory、RMSD/接触/H 键、MM-GBSA
report/                              # summary.tsv / report.md / report.json
```

糖链、二价金属、共价配体、非 MSE 修饰残基、多个配体候选、电荷模糊时默认**硬报错**并列出具体残基；对应的显式逃生口见 [`docs/ERRORS.md`](docs/ERRORS.md)。

### 3.2 高级：amber-kit 旧流程（兼容）

> 以下步骤保留给需要精细控制、RNA 体系或已有项目的用户；普通蛋白–配体/肽/抗体请优先用 `amber-auto`。

```bash
# 0. 复制一份（每个体系/项目一份，自包含）
cp -r ~/Documents/amber_general_workflow ~/Documents/my_project
cd ~/Documents/my_project

# 1. 初始化（生成 config/config.sh、config/systems.tsv、目录）
bin/amber-kit init

# 2. 放输入
#    inputs/receptor.pdb  受体（或受体+配体复合物）
#    inputs/lig.pdb       配体（必须含显式氢、单残基）
#    编辑 config/config.sh：SYSTEM_TYPE、RECEPTOR_FF、WATER、LIGAND_CHARGE ...
#    编辑 config/systems.tsv：每个配体一行

# 3. 检查环境/配置
bin/amber-kit check

# 4. 配体参数化 + 建系 + QC + 推导 masks
bin/amber-kit build

# 5. 生成每个副本的 MD 输入
bin/amber-kit inputs

# 6. 先干跑检查 PBS 头和命令，再提交
bin/amber-kit submit equil --label lig1 --rep 1 --dry-run
bin/amber-kit submit equil --label lig1 --rep 1

# 7. npt.rst7 出现后提交生产
bin/amber-kit submit prod --label lig1 --rep 1

# 8. 进度
bin/amber-kit status --queue

# 9. 分析 + MM/GBSA
bin/amber-kit analyze --label lig1 --rep 1 --force
bin/amber-kit mmpbsa  --label lig1 --rep 1 --force
bin/amber-kit report
```

---

## 4. 目录

```
config/config.example.sh     全局配置（力场/水/离子/PBS/分析）
config/systems.example.tsv   每个配体/模型/副本一行
bin/amber-kit                统一入口
scripts/                     建系/模拟输入/分析/PBS 辅助
templates/                   tleap、Amber &cntrl、cpptraj、MMPBSA 模板
pbs/                         PBS 模板 + 前缀详解
docs/                        速查、排错
examples/                    protein-ligand / rna-ligand 示例配置
tests/                       smoke test
```

项目运行后还会有：

```
inputs/       你的 PDB
ligands/      GAFF2/AM1-BCC 参数
systems/      complex.prmtop / complex.inpcrd / masks.env / build_qc.tsv
md/<label>/r<rep>/  MD 输入、轨迹、restart
analysis/     cpptraj、python npz、mmpbsa、tables
logs/         PBS 日志
```

---

## 5. 常用命令

```bash
bin/amber-kit check
bin/amber-kit build --label lig1 --render-only     # 只渲染 leap.in
bin/amber-kit inputs --label lig1 --rep 1 --print-only
bin/amber-kit submit equil --label lig1 --rep 1 --dry-run
bin/amber-kit submit prod --label lig1 --rep 1
bin/amber-kit submit prod-chunk --label lig1 --rep 1 --chunk c02
bin/amber-kit submit multi-gpu --pair lig1:1,lig1:2 --dry-run
bin/amber-kit status --queue
bin/amber-kit analyze --label lig1 --rep 1 --force
bin/amber-kit mmpbsa --label lig1 --rep 1 --check-only
bin/amber-kit report
```

单独运行底层脚本：

```bash
source /opt/software/amber22/amber.sh
scripts/build_system.sh --help
scripts/param_ligand.sh --help
scripts/setup_md_inputs.py --help
scripts/analyze_replica.py --help
scripts/mmpbsa_one.sh --help
```

---

## 6. 重要原则

1. **配体电荷/质子化必须手填**；脚本不猜。
2. `antechamber -at gaff2` ⇔ `parmchk2 -s gaff2` ⇔ `leaprc.gaff2`，三者必须一致。
3. `tleap` 必须 `Errors = 0`；`ATTN` 默认失败，确认化学合理后才 `--allow-attn`。
4. `igb=5 ↔ mbondi2`，`igb=8 ↔ mbondi3`；改 `igb` 必须重建拓扑。
5. 生产 chunk 有状态，严格顺序；本工具包的 PBS 模板会自动跳过已完成 stage。
6. MM/GBSA 用**溶剂化拓扑** + `ante-MMPBSA.py`；不要用 dry prmtop。
7. 少于 5 个分子不要宣称泛化；3 个分子完美排序的随机概率是 1/6；MM/GBSA 不是 Kd。
8. 不要为了匹配实验排序去调介电常数、盐浓度、窗口或 QC 阈值。

---

## 7. 复制到其他集群

只改 `config/config.sh`：

```bash
AMBER_SH=/path/to/amber.sh
PY_QC=/path/to/python-with-parmed
PY_ANALYSIS=/path/to/python-with-MDAnalysis
PBS_QUEUE=batch
PBS_PPN=8
PBS_GPUS=1                 # 单个 pmemd.cuda 只用 1 张
PBS_GPUS_PER_NODE=4        # 节点实际 GPU 数（本机 4）
PBS_TWO_GPU_PPN=48         # 整节点 CPU 核数
PBS_LOG_DIR=$PROJECT_ROOT/logs
```

GPU 数量不同：

- 每个作业仍用 1 张 GPU：`PBS_GPUS=1` 不变，节点有几张只影响并发。
- 一个作业跑多个副本：用 `submit multi-gpu`，把 `PBS_GPUS_PER_NODE` 改成节点 GPU 数，`PBS_TWO_GPU_PPN` 改成节点核数。

详见 [`pbs/README.md`](pbs/README.md) 第 6 节。
