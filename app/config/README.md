# 配置说明

## `auto.example.json`（amber-auto）

一键入口 `bin/amber-auto` 使用 JSON 默认配置，不需要手改 shell 文件。可复制为 `config/auto.json` 后用 `--config` 指定，或直接命令行覆盖。字段包括工具路径、力场/水/盒子、MD 长度、PBS 资源和 MM-GBSA 参数。

```bash
bin/amber-auto run complex.pdb --config config/auto.json
```

## `config.sh`（amber-kit 高级路径）

由 `config.example.sh` 复制生成。所有值都可以被命令行/`systems.tsv` 覆盖（配体字段）。

### 路径与解释器

| 变量 | 说明 |
|---|---|
| `PROJECT_ROOT` | 自动设为 `config/..`；一般不要改 |
| `PROJECT_NAME` | 作业名前缀；不要有空格 |
| `AMBER_SH` | Amber 初始化脚本；其他集群改成自己的路径 |
| `PY_QC` | 有 `parmed` + `numpy` 的 Python（建系/QC） |
| `PY_ANALYSIS` | 有 `MDAnalysis`/`scipy`/`shapely`/`pandas` 的 Python（分析） |
| `MPIRUN` | `mpirun` 路径/命令 |

### 体系与力场

| 变量 | 取值 | 说明 |
|---|---|---|
| `SYSTEM_TYPE` | `rna-ligand` / `protein-ligand` / `rna-apo` / `protein-apo` | 决定 build/inputs/分析默认 |
| `RECEPTOR_PDB` | 路径 | 受体或受体+配体复合物 |
| `RECEPTOR_MODEL` | 整数/空 | 多模型 NMR 选第 N 个模型 |
| `LIGAND_RESNAME` | 1–4 字符 | 配体残基名 |
| `LIGAND_CHARGE` | 整数/空 | **必须手填**；空则 holo build 报错 |
| `GAFF` | `gaff2`/`gaff` | 必须与 `antechamber -at`、`parmchk2 -s`、`leaprc.gaff*` 一致 |
| `RECEPTOR_FF` | `RNA.YIL`/`RNA.OL3`/`RNA.ROC`/`RNA.LJbb`/`ff14SB`/`ff19SB` | 力场 |
| `SOLVENT` | `explicit`/`implicit` | 本工具包以 explicit 为主 |
| `WATER` | `tip3p`/`opc`/`spce`/`tip4pd`/`fb3` | 与力场配对 |
| `BOX` | `oct`/`rect` | `solvateOct`/`solvateBox` |
| `PAD` | Å | 溶质到盒边；10–12 |
| `CATION`/`ANION` | `K+`/`Na+`/`Cl-` | 中和 |
| `COUNTERION` | `cation`/`anion`/`none` | 中和方向 |
| `SALT_PAIRS` | 整数 | `addIonsRand` 盐对；仅 explicit |
| `RADII` | `mbondi2`/`mbondi3` | `igb=5`↔`mbondi2`；`igb=8`↔`mbondi3` |
| `IGB` | 5/8/... | MM/GBSA 的 GB 模型；与 `RADII` 匹配 |

### masks

| 变量 | 说明 |
|---|---|
| `RESTRAINT_MASK` | `auto` 或 Amber mask；`auto` 由 `residue_masks.py` 从拓扑推导 |
| `RECEPTOR_MASK`/`LIGAND_MASK` | 同上；供分析/MMPBSA |
| `STRIP_MASK` | 去水/离子；必须包含体系实际有的离子 |
| `ANCHOR_MASK`/`FIXED_MASK` | 双链 RNA 的 `autoimage` 锚定；蛋白留空 |

### MD

| 变量 | 说明 |
|---|---|
| `HEAT_NS` | NVT 加热长度，默认 2 ns |
| `RELEASE_STAGES` | NPT 约束释放段数，默认 10 |
| `RELEASE_NS` | 每段长度，默认 1 ns |
| `PROD_CHUNKS` | 生产 chunk 数，默认 4 |
| `PROD_CHUNK_NS` | 每 chunk 长度，默认 25 ns |
| `DT` | 积分步长，默认 0.002 ps |
| `TRAJ_INTERVAL_PS` | 轨迹写盘间隔，默认 10 ps |
| `HMR` / `HMR_MASS` | `1` 开启 HMR（非水氢 3.024 Da）；需 `DT=0.004`，与 `--render-only` 互斥 |
| `SEED_BASE` | 未在 `systems.tsv` 给 `seed_base` 时的默认种子 |

### PBS

| 变量 | 说明 |
|---|---|
| `PBS_QUEUE` | 队列 |
| `PBS_NODES`/`PBS_PPN`/`PBS_GPUS` | 单 GPU 作业资源 |
| `PBS_WALLTIME` | `HH:MM:SS` |
| `PBS_JOIN`/`PBS_MAIL` | `oe`/`n` |
| `PBS_LOG_DIR` | 日志目录；必须存在 |
| `PBS_TWO_GPU_PPN` | 整节点模板的 ppn |
| `PBS_GPUS_PER_NODE` | 节点实际 GPU 数；`multi-gpu` 使用 |

### 分析

| 变量 | 说明 |
|---|---|
| `ANALYSIS_WINDOW` | `last50`/`last25`/`all` |
| `CONTACT_CUTOFF` | 配体–受体接触，Å |
| `HBOND_HA_CUTOFF`/`HBOND_DA_CUTOFF`/`HBOND_ANGLE` | H 键判据 |
| `CONTACT_STRIDE`/`SASA_STRIDE` | 接触/SASA 抽帧 |
| `MMPBSA_METHOD` | `gb`/`pb` |
| `MMPBSA_IGB`/`MMPBSA_SALTCON` | GB 参数 |
| `MMPBSA_INTERVAL` | 每 N 帧取一帧 |
| `MMPBSA_DECOMP` | `auto` 或 `print_res` 字符串 |
| `MMPBSA_RADII` | 拆分拓扑的 radii |
| `RNA_RESRANGE` | RNA nastruct 范围；空=auto |
| `LOOP_RESIDUES` | loop 残基，逗号分隔 |
| `LIGAND_RING_ATOMS` | 配体环原子组，`name:atom1,atom2;...` |
| `STACK_TARGETS` | 堆积目标残基号，逗号分隔 |

---

## `systems.tsv`

制表符分隔；第一行表头（`label ...`）。空列要保留制表符。

| 列 | 说明 |
|---|---|
| `label` | 体系标签；用于 `systems/<label>`、`md/<label>/r<rep>`、分析文件名 |
| `ligand_pdb` | 配体 PDB（含显式 H、单残基）；apo 填 `-` |
| `charge` | 配体净电荷整数；apo 填 `-` |
| `resname` | 配体残基名 1–4 字符；apo 填 `-` |
| `models` | 逗号分隔模型号；空=单模型；与 `replicas` 顺序对应 |
| `replicas` | 逗号分隔副本号；空=1 |
| `seed_base` | r1 的种子；rN = seed_base + (N-1) |
| `control` | `1`=对照（进排名表，不进统计） |
| `experimental` | IC50/Kd（µM）或 ΔG（kcal/mol）；无则 `NA` |
| `experimental_type` | `ic50`/`kd`/`dg`/`NA` |
| `notes` | 备注；不要含制表符 |

### 示例

```tsv
label	ligand_pdb	charge	resname	models	replicas	seed_base	control	experimental	experimental_type	notes
lig1	inputs/lig1.pdb	2	LIG	7,11,10	1,2,3	20261101	0	7.9	ic50	A1AZM
lig3_p2	inputs/lig3_p2.pdb	2	LIG	15,14,9	1,2,3	20261401	1	NA	NA	MQC +2 control
```

---

## 复制到其他集群

只改 `config.sh` 的这几项：

```bash
AMBER_SH=...
PY_QC=...
PY_ANALYSIS=...
PBS_QUEUE=...
PBS_PPN=...
PBS_GPUS=1
PBS_GPUS_PER_NODE=<节点 GPU 数>
PBS_TWO_GPU_PPN=<节点 CPU 核数>
PBS_LOG_DIR=...
```

GPU 数量不同：

- 单个 `pmemd.cuda` 永远只要 1 张；`PBS_GPUS=1`。
- 一个作业跑 N 个副本：`submit multi-gpu`，`PBS_GPUS_PER_NODE=N`。
