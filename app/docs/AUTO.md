# amber-auto 一键流程

`bin/amber-auto` 面向日常使用：给一个或多个复合物 PDB，自动完成

```
分类 → 清洗/补氢/电荷推断 → Reduce 质子化 → 二硫键 → tleap 建系 + QC
     → masks → MD 输入 → qsub equil → prod → post(cpptraj + 分析 + MM-GBSA)
```

## 1. 环境要求

```bash
source /opt/software/amber22/amber.sh          # Amber22
/opt/software/amber_general_workflow/envs/analysis/bin/python -c "import rdkit"   # 电荷推断
/opt/software/amber_general_workflow/envs/qc/bin/python -c "import parmed, numpy"
/opt/software/amber_general_workflow/envs/analysis/bin/python -c "import MDAnalysis, scipy"
command -v obabel || ls /opt/schrodinger2021-2/utilities/obabel
qstat -q                                    # Torque batch 队列
```

`amber-auto` 不要求复制仓库、不要求手改 `config/config.sh` 或 `systems.tsv`。默认参数在 `config/auto.example.json`，可用 `--config` 指向自己的 JSON。
共享安装可把站点配置放在任意位置，用环境变量 `AMBER_AUTO_CONFIG` 指定（`--config` 优先于环境变量），这样只读的共享代码不需要每个用户改。

## 2. 基本用法

```bash
# 一个或多个复合并行建系，全部成功后自动提交
bin/amber-auto run complex1.pdb complex2.pdb

# 只分类，不建系
bin/amber-auto inspect complex1.pdb

# 只建系 + 生成输入，不提交
bin/amber-auto build complex1.pdb

# 渲染 PBS 并打印依赖链，不 qsub
bin/amber-auto run complex1.pdb --dry-run

# 进度、报告、断点续跑
bin/amber-auto status --queue
bin/amber-auto report
bin/amber-auto resume
```

输出目录默认 `./amber-campaign`，可用 `--campaign DIR` 指定。

## 3. 支持的体系

| 类型 | 自动识别依据 | 建系方式 |
|---|---|---|
| 蛋白–小分子 | 非水/非离子/非添加剂的 HETATM 候选，按接触和埋藏评分 | 受体 ff19SB + 配体 GAFF2/AM1-BCC |
| 蛋白–肽 | ≥2 条蛋白链，其中短链 ≤50 残基且与长链接触 | 整条复合物按蛋白力场建系 |
| 蛋白–抗体 | 90–500 残基、含链内二硫键的 Ig 样链，与非 Ig 链接触 | 同上 |
| 蛋白–蛋白 | 其他多链蛋白复合物 | 同上 |
| apo 蛋白 | 无配体、单链或链断裂片段 | 同上 |

可覆盖：

```bash
--type protein-ligand|protein-peptide|protein-antibody|protein-protein|apo-protein
--ligand RESNAME[:CHAIN[:RESID]]
--receptor-chains A --partner-chains B
--charge -2                 # 配体净电荷
--ligand-smiles 'O=C([O-])c1ccccc1'   # 校验键级/电荷
```

## 4. 配体识别与电荷

- 配体候选必须含 ≥3 个重原子、与蛋白有接触、非水/离子/结晶添加剂。
- 多个候选分数接近、同 resname 多拷贝、共价配体（<1.9 Å）默认硬报错。
- 配体已有氢：RDKit 枚举键级/电荷；无氢：OpenBabel 先按 pH 7.4 补氢并感知键级，RDKit 校验。
- 电荷模糊时硬报错，列出候选 q/SMILES/score；使用 `--charge` 或 `--ligand-smiles` 覆盖。
- 所有推断结果写入 `prepare/ligand_charge.json` 和 `campaign.json`。

## 5. 蛋白质子化、二硫键、链断裂

- 蛋白部分先运行 `pdb4amber --reduce`（不加 `--nohyd`），保留 HID/HIE/HIP 与 Asn/Gln 翻转，再去氢交给 tleap 补氢。
- SSBOND + SG–SG <2.5 Å 自动生成 `CYX` 和 `bond` 命令。
- C–N >2.0 Å 或编号跳变处插入 `TER`；两段各自带电荷末端，并在报告中警告。
- MSE 自动转 MET/SE→SD；其他修饰残基默认报错。

## 6. 默认 MD 设置

| 项 | 默认 |
|---|---|
| 力场/水 | ff19SB + OPC（蛋白）；配体 GAFF2/AM1-BCC |
| 溶剂 | 显式，截角八面体，10 Å padding |
| 离子 | Na+/Cl− 中和，不加额外盐 |
| 平衡 | 2 ns 加热 + 10×1 ns NPT 释放 |
| 生产 | 4×25 ns = 100 ns |
| 副本 | 1 |
| 分析与 MM-GBSA | cpptraj autoimage/strip → RMSD/Rg/接触/H 键/SASA → MM/GBSA (igb=5, 0.12 M) |

细长蛋白、长肽或多链抗体建议用 `--box rect`，否则 `solvateOct` 会按最长边扩成立方体，产生大量水。球状蛋白用默认 `oct` 即可。

常用覆盖：

```bash
--replicas 3 --prod-ns 200
--solvent implicit --prod-ns 20
--ff ff14SB --water tip3p
--box rect --pad 9          # 细长/多链体系
--prod-chunks 8 --prod-ns 200
--no-mmpbsa
--mmpbsa-method pb --mmpbsa-igb 5 --mmpbsa-saltcon 0.15
```

## 7. PBS 依赖链

每个副本提交三个 Torque 作业：

```bash
qsub run_equil.pbs
qsub -W depend=afterok:<equil_id> run_prod.pbs
qsub -W depend=afterok:<prod_id> run_post.pbs
```

- `run_equil` / `run_prod` 复用现有模板，自带 stage 跳过/重算和 `flock`。
- `run_post` 是 CPU 作业，运行 `scripts/postprocess_one.sh`。
- `campaign.json` 记录 job id、依赖、提交时间；`status` 合并 qstat 与文件证据。
- 生产被 walltime 杀掉后，已完成 chunk 会保留；`amber-auto resume` 从下一个缺失阶段重投。
- 需要固定或避开节点时用 `--pbs-nodes node02`（或配置 `"pbs_nodes": "node02"`）；默认 `1` 由 Torque 调度。

## 8. 目录与状态

```
amber-campaign/
  campaign.json
  complexes/<label>/
    components.json
    prepare/{protein_prepared.pdb,ligand_prepared.pdb,extra_leap.in,
             residue_map.tsv,prepare_manifest.json,ligand_charge.json}
    ligands/LIG/{lig.mol2,lig.frcmod}
    systems/<label>/{complex.prmtop,complex.inpcrd,masks.env,
                     build_manifest.json,build_qc.tsv,build_linkage.json}
    md/<label>/r1/{01_min_restraint.in,...,run_equil.pbs,run_prod.pbs,run_post.pbs}
    analysis/{cpptraj,python,mmpbsa,tables}
  report/{summary.tsv,report.json,report.md}
```

`amber-auto status` 输出每个副本的最后 stage、NSTEP、ns/day、作业状态和 ΔG；`amber-auto report` 汇总每个复合物的 QC、MD 进度、RMSD/接触和 MM-GBSA。

## 9. 示例

```bash
# 蛋白–小分子，3 副本，200 ns
bin/amber-auto run kinase.pdb --replicas 3 --prod-ns 200

# 蛋白–肽，显式指定链
bin/amber-auto run fab_peptide.pdb --type protein-peptide \
  --receptor-chains A,B --partner-chains C

# 先用隐式溶剂快速筛查，再跑显式
bin/amber-auto run candidate.pdb --solvent implicit --prod-ns 20
bin/amber-auto run candidate.pdb --campaign ./explicit --prod-ns 100

# 配体电荷有歧义时
bin/amber-auto inspect complex.pdb           # 查看候选
bin/amber-auto run complex.pdb --charge -1   # 明确指定
```

## 10. 限制

- 不支持 N-糖链、二价金属、共价配体、非 MSE 修饰残基、蛋白–RNA/DNA。
- 仅支持单残基小分子配体；多拷贝/多候选需显式选择。
- MM-GBSA 是相对指标，不含熵；蛋白–蛋白/抗体 ΔG 只能作界面相对比较，不跨靶点排名。
- 默认只提交到 cluster2 的 Torque `batch` 队列；其他调度器需改 PBS 模板。
