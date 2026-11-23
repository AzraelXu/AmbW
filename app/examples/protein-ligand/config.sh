#!/bin/bash
# 示例：蛋白–配体（ff19SB + OPC + Na+）
# 使用方法：复制到 <project>/config/config.sh，再按需修改。
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_NAME="protein_lig"

AMBER_SH="/home/ylqx/amber22/amber.sh"
PY_QC="/opt/software/amber_general_workflow/envs/qc/bin/python"
PY_ANALYSIS="/opt/software/amber_general_workflow/envs/analysis/bin/python"
MPIRUN="mpirun"

SYSTEM_TYPE="protein-ligand"
RECEPTOR_PDB="inputs/receptor.pdb"
RECEPTOR_MODEL=""
LIGAND_RESNAME="LIG"
LIGAND_CHARGE=""          # 必须手填
GAFF="gaff2"

RECEPTOR_FF="ff19SB"
SOLVENT="explicit"
WATER="opc"
BOX="oct"; PAD=10.0
CATION="Na+"; ANION="Cl-"; COUNTERION="cation"; SALT_PAIRS=0
RADII="mbondi2"; IGB=5

RESTRAINT_MASK="auto"; RECEPTOR_MASK="auto"; LIGAND_MASK="auto"
STRIP_MASK=":WAT,Na+,Cl-,K+"
ANCHOR_MASK=""; FIXED_MASK=""

HEAT_NS=2; RELEASE_STAGES=10; RELEASE_NS=1
PROD_CHUNKS=4; PROD_CHUNK_NS=25; DT=0.002; TRAJ_INTERVAL_PS=10
HMR=0; HMR_MASS=3.024     # HMR=1 needs DT=0.004 and a minimized start
SEED_BASE=20261001

PBS_QUEUE="batch"; PBS_NODES=1; PBS_PPN=8; PBS_GPUS=1
PBS_WALLTIME="24:00:00"; PBS_JOIN="oe"; PBS_MAIL="n"
PBS_LOG_DIR="$PROJECT_ROOT/logs"
PBS_TWO_GPU_PPN=48
PBS_GPUS_PER_NODE=2

ANALYSIS_WINDOW="last50"
CONTACT_CUTOFF=4.0; HBOND_HA_CUTOFF=2.5; HBOND_DA_CUTOFF=3.5; HBOND_ANGLE=120
CONTACT_STRIDE=10; SASA_STRIDE=10
MMPBSA_METHOD="gb"; MMPBSA_IGB=5; MMPBSA_SALTCON=0.12
MMPBSA_INTERVAL=100; MMPBSA_DECOMP="auto"; MMPBSA_RADII="$RADII"
RNA_RESRANGE=""; LOOP_RESIDUES=""
LIGAND_RING_ATOMS=""; STACK_TARGETS=""
