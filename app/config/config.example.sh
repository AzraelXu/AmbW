#!/bin/bash
# =============================================================================
# amber_general_workflow 全局配置
#
# 用法：
#   cp config/config.example.sh config/config.sh
#   然后按本文件注释修改 config/config.sh
#
# systems.tsv 中每一行可以覆盖本文件的配体字段（ligand_pdb/charge/resname）。
# 本文件只放"全局"选择：体系类型、力场、水、离子、盒子、MD 长度、PBS 资源、
# 分析默认值。
# =============================================================================

# -----------------------------------------------------------------------------
# 项目路径
# -----------------------------------------------------------------------------
# config.sh 位于 <project>/config/config.sh，因此 PROJECT_ROOT 自动是 <project>
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_NAME="my_project"

# -----------------------------------------------------------------------------
# 可执行环境
# -----------------------------------------------------------------------------
AMBER_SH="/home/ylqx/amber22/amber.sh"
# QC/建系需要 parmed；本机 ambertools_env 有 parmed 4.2 + numpy 1.24
PY_QC="/opt/software/amber_general_workflow/envs/qc/bin/python"
# 分析需要 MDAnalysis/scipy + numpy；本机 prolif 环境（MDA 2.7 / scipy 1.13）。
# 仅 RNA 堆积分析 analyze_stacking.py 额外需要 shapely，本机该环境未安装。
PY_ANALYSIS="/opt/software/amber_general_workflow/envs/analysis/bin/python"
MPIRUN="mpirun"

# -----------------------------------------------------------------------------
# 体系类型
# -----------------------------------------------------------------------------
# rna-ligand | protein-ligand | rna-apo | protein-apo
SYSTEM_TYPE="rna-ligand"

# 受体 PDB。可以是：
#   - 只有受体（推荐）
#   - 受体 + 配体的复合物（amber-kit build 会用 split_complex.py 拆开）
#   - 多模型 NMR PDB（用 RECEPTOR_MODEL 或 systems.tsv 的 models 列选模型）
RECEPTOR_PDB="inputs/receptor.pdb"
RECEPTOR_MODEL=""          # 例如 7；空 = 使用第 1 个模型 / 单模型文件

# 配体默认值；systems.tsv 每行可覆盖
LIGAND_RESNAME="LIG"
LIGAND_CHARGE=""           # 必须手填（+2、+1、0、-1 ...）；空则 build 报错
GAFF="gaff2"               # gaff2 | gaff

# -----------------------------------------------------------------------------
# 力场 / 水 / 盒子 / 离子
# -----------------------------------------------------------------------------
# RNA:     RNA.YIL（推荐）、RNA.OL3、RNA.ROC、RNA.LJbb
# Protein: ff14SB、ff19SB
RECEPTOR_FF="RNA.YIL"

# explicit | implicit
SOLVENT="explicit"
# explicit 时的水模型：tip3p | opc | spce | tip4pd | fb3
# 推荐配对：RNA.YIL + tip3p；ff14SB + tip3p；ff19SB + opc
WATER="tip3p"

# oct（截角八面体）| rect（长方体）
BOX="oct"
PAD=10.0                   # 溶质到盒边距离，Å

# 中和离子：RNA 常用 K+，蛋白常用 Na+；ANION 用于带正电体系
CATION="K+"
ANION="Cl-"
COUNTERION="cation"        # cation | anion | none
SALT_PAIRS=0               # >0 时 addIonsRand CATION N ANION N（仅 explicit）

# PBRadii：igb=5 -> mbondi2；igb=8 -> mbondi3
RADII="mbondi2"
IGB=5

# -----------------------------------------------------------------------------
# 掩码
# -----------------------------------------------------------------------------
# auto = 建系后由 residue_masks.py 从 prmtop 自动推导
RESTRAINT_MASK="auto"      # 会变成 ":1-N & !@H="
RECEPTOR_MASK="auto"       # 会变成 ":1-..."（不含配体）
LIGAND_MASK="auto"         # 会变成 ":LIG" 或 ":27"
STRIP_MASK=":WAT,K+,Cl-,Na+"

# 双链 RNA 的 autoimage 锚定；蛋白/单链留空
ANCHOR_MASK=""             # 例如 ":1-13"
FIXED_MASK=""              # 例如 ":14-26,27"

# -----------------------------------------------------------------------------
# MD 长度与采样
# -----------------------------------------------------------------------------
HEAT_NS=2                  # NVT 加热
RELEASE_STAGES=10          # NPT 约束释放段数（默认 10 段）
RELEASE_NS=1               # 每段长度
PROD_CHUNKS=4              # 生产 chunk 数
PROD_CHUNK_NS=25           # 每 chunk 长度；4 × 25 = 100 ns
DT=0.002                   # ps；SHAKE 下 0.002
TRAJ_INTERVAL_PS=10        # 轨迹写盘间隔；10 ps -> ntwx=5000
SEED_BASE=20261001         # 每行 systems.tsv 的 seed_base 优先

# -----------------------------------------------------------------------------
# PBS 资源（Torque / cluster2）
# -----------------------------------------------------------------------------
PBS_QUEUE="batch"
PBS_NODES=1
PBS_PPN=8
PBS_GPUS=1
PBS_WALLTIME="24:00:00"
PBS_JOIN="oe"
PBS_MAIL="n"
PBS_LOG_DIR="$PROJECT_ROOT/logs"   # 必须存在；amber-kit init 会创建
PBS_TWO_GPU_PPN=48                 # 整节点多 GPU fallback 的 ppn
PBS_GPUS_PER_NODE=4                # 本机 node02–node08 每节点 4×RTX3090


# -----------------------------------------------------------------------------
# 分析默认值
# -----------------------------------------------------------------------------
ANALYSIS_WINDOW="last50"   # last50 | last25 | all
CONTACT_CUTOFF=4.0         # 配体-受体接触，Å
HBOND_HA_CUTOFF=2.5        # H···A，Å
HBOND_DA_CUTOFF=3.5        # D···A，Å
HBOND_ANGLE=120            # D-H···A，度
CONTACT_STRIDE=10          # 接触计算抽帧
SASA_STRIDE=10             # cpptraj surf 抽帧

MMPBSA_METHOD="gb"         # gb | pb
MMPBSA_IGB=5               # gb 时：5（mbondi2）或 8（mbondi3）
MMPBSA_SALTCON=0.12        # M
MMPBSA_INTERVAL=100        # 每 100 帧取一帧（10 ps 帧 -> 1 ns）
MMPBSA_DECOMP="auto"       # auto -> "1-<last solute residue>"
MMPBSA_RADII="$RADII"

# RNA 专有分析
RNA_RESRANGE=""            # 空 = auto（从拓扑识别 A/C/G/U）
LOOP_RESIDUES=""           # 例如 "7,20"
LIGAND_RING_ATOMS=""       # 例如 "ring1:C1,C2,C3,C4,C5,C6;ring2:C9,C10,C11,C12,C13,C14"
STACK_TARGETS=""           # 例如 "7,20,6,8,19,21"
