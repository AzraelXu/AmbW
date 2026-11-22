# 命令速查

## 环境

```bash
source /opt/software/amber22/amber.sh
which pmemd.cuda cpptraj tleap antechamber parmchk2
qstat -a
```

## 一键流程（amber-auto，推荐）

```bash
# 建系 + 自动提交 equil -> prod -> post
bin/amber-auto run complex1.pdb complex2.pdb

bin/amber-auto inspect complex1.pdb          # 只分类
bin/amber-auto build complex1.pdb            # 只建系/生成输入
bin/amber-auto run complex1.pdb --dry-run    # 只打印 PBS 依赖链
bin/amber-auto status --queue
bin/amber-auto report
bin/amber-auto resume

# 常用覆盖
--replicas 3 --prod-ns 200
--ligand 0QE:A:501 --charge -2
--type protein-peptide --receptor-chains A --partner-chains B
--solvent implicit --prod-ns 20
--no-mmpbsa
--drop-glycans / --drop-metals / --drop-other-ligands
```

详见 [`AUTO.md`](AUTO.md)、错误处理 [`ERRORS.md`](ERRORS.md)。

## 项目（amber-kit 高级路径）

```bash
cp -r ~/Documents/amber_general_workflow ~/Documents/my_project
cd ~/Documents/my_project
bin/amber-kit init
bin/amber-kit check
bin/amber-kit check --strict
```

## 输入准备

```bash
# 从复合物拆受体/配体
python -I scripts/split_complex.py --pdb complex.pdb --model 7 \
  --ligand-resname MQC --out-resname LIG \
  --receptor-out receptor.pdb --ligand-out lig.pdb

# 取一个 NMR 模型
python -I scripts/extract_model.py --pdb allmodels.pdb --model 7 --out model_07.pdb

# 选 medoid + diverse
python -I scripts/select_models.py --pdb allmodels.pdb --rna \
  --core-residues 1-6,8-13,14-19,21-26 --loop-residues 7,20 \
  --terminal-residues 1,13,14,26 --ligand-resname LIG --n-select 3 \
  --out selected.tsv
```

## 配体参数化

```bash
export AMBER_MD_PYTHON=/opt/software/amber_general_workflow/envs/qc/bin/python
scripts/param_ligand.sh --pdb lig.pdb --charge 2 --resname LIG \
  --gaff gaff2 --outdir ligands/lig1/template

# 姿态移植
scripts/param_ligand.sh --pdb pose07.pdb --charge 2 --resname LIG \
  --gaff gaff2 --outdir ligands/lig1/model_07 \
  --transplant ligands/lig1/template/lig.mol2 --onto pose07.pdb
```

## 建系

```bash
# 整套
bin/amber-kit build
bin/amber-kit build --label lig1
bin/amber-kit build --label lig1 --render-only

# 单步
scripts/build_system.sh --type rna-ligand \
  --receptor rna.pdb --ligand lig.mol2 --frcmod lig.frcmod \
  --lig-resname LIG --ff RNA.YIL --solvent explicit --water tip3p \
  --box oct --pad 10 --cation K+ --counterion cation \
  --radii mbondi2 --outdir systems/lig1

python -I scripts/qc_system.py --prmtop systems/lig1/complex.prmtop \
  --inpcrd systems/lig1/complex.inpcrd --ligand-resname LIG \
  --expect-lig-charge 2 --expect-net-charge 0 --expect-water yes \
  --expect-radii mbondi2 --restraint-mask ':1-27 & !@H=' --strict

python -I scripts/residue_masks.py --prmtop systems/lig1/complex.prmtop \
  --ligand-resname LIG --json masks.json --env masks.env --force
```

## MD 输入

```bash
bin/amber-kit inputs
bin/amber-kit inputs --label lig1 --rep 1
bin/amber-kit inputs --label lig1 --rep 1 --print-only

python -I scripts/setup_md_inputs.py --out md/lig1/r1 --seed 20261101 \
  --receptor-type rna --solvent explicit --restraint-mask ':1-27 & !@H=' \
  --prmtop systems/lig1/complex.prmtop --ligand-resname LIG \
  --heat-ns 2 --release-stages 10 --release-ns 1 \
  --prod-chunks 4 --prod-chunk-ns 25 --dt 0.002 --traj-interval-ps 10
```

## PBS

```bash
bin/amber-kit submit equil --label lig1 --rep 1 --dry-run
bin/amber-kit submit equil --label lig1 --rep 1
bin/amber-kit submit prod --label lig1 --rep 1
bin/amber-kit submit prod-chunk --label lig1 --rep 1 --chunk c02
bin/amber-kit submit two-gpu --pair lig1:1,lig1:2 --dry-run
bin/amber-kit submit multi-gpu --pair lig1:1,lig1:2,lig2:1,lig2:2 --dry-run
bin/amber-kit submit param --label lig1 --dry-run
bin/amber-kit submit analysis --label lig1 --rep 1 --dry-run
bin/amber-kit submit mmpbsa --label lig1 --rep 1 --dry-run

qstat -a
qstat -f <jobid>
qdel <jobid>
qdel -p <jobid>
```

## 分析

```bash
bin/amber-kit status --queue

bin/amber-kit analyze --label lig1 --rep 1 --force
bin/amber-kit analyze --label lig1 --rep 1 --rna --force
bin/amber-kit analyze --label lig1 --rep 1 --no-sasa --force

python -I scripts/make_cpptraj_inputs.py --top md/lig1/r1/complex.prmtop \
  --traj 'md/lig1/r1/05_prod_c*.nc' --strip-mask ':WAT,K+,Cl-,Na+' \
  --anchor ':1-13' --fixed ':14-26,27' \
  --outdir analysis/cpptraj --prefix lig1_r1 --nastruct --resrange 1-26 --run

python -I scripts/analyze_replica.py --prmtop analysis/cpptraj/lig1_r1_dry.prmtop \
  --traj analysis/cpptraj/lig1_r1_dry.nc --label lig1 --rep 1 \
  --ligand-resname LIG --rna --out analysis/python/lig1_r1.npz --force

python -I scripts/analyze_ligand_hbonds.py --prmtop analysis/cpptraj/lig1_r1_dry.prmtop \
  --traj analysis/cpptraj/lig1_r1_dry.nc --ligand-resname LIG \
  --out analysis/python/lig1_r1_hb.npz --force

python -I scripts/analyze_stacking.py --prmtop analysis/cpptraj/lig1_r1_dry.prmtop \
  --traj analysis/cpptraj/lig1_r1_dry.nc --ligand-resname LIG \
  --ligand-ring-atoms 'ring1:C1,C2,C3,C4,C5,C6' --target-residues 7,20 \
  --out analysis/python/lig1_r1_stack.npz --force

python -I scripts/parse_nastruct.py --prefix analysis/cpptraj/lig1_r1 \
  --label lig1 --rep 1 --out analysis/tables/nastruct_lig1_r1.csv --force
```

## MM/GBSA

```bash
bin/amber-kit mmpbsa --label lig1 --rep 1 --check-only
bin/amber-kit mmpbsa --label lig1 --rep 1 --force
bin/amber-kit mmpbsa --label lig1 --rep 1 --method pb --force

scripts/mmpbsa_one.sh --top md/lig1/r1/complex.prmtop \
  --traj analysis/cpptraj/lig1_r1_dry.nc --ligand-mask :LIG \
  --strip-mask ':WAT,K+,Cl-,Na+' --radii mbondi2 --method gb --igb 5 \
  --saltcon 0.12 --interval 100 --decomp-residues 1-27 \
  --outdir analysis/mmpbsa/lig1_r1
```

## 汇总

```bash
bin/amber-kit report

python -I scripts/aggregate_results.py --root . --systems config/systems.tsv \
  --window last50 --out analysis/tables --force

python -I scripts/rank_ligands.py --root . --systems config/systems.tsv \
  --out analysis/tables/ranking.csv --force
```

## 结果文件

```
systems/<label>/complex.prmtop
systems/<label>/complex.inpcrd
systems/<label>/masks.env
systems/<label>/build_qc.tsv
md/<label>/r<rep>/05_prod_c*.nc
analysis/cpptraj/<label>_r<rep>_dry.nc
analysis/python/<label>_r<rep>.npz
analysis/mmpbsa/<label>_r<rep>/mmpbsa.dat
analysis/tables/ranking.csv
```
