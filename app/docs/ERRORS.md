# amber-auto 错误与逃生口

`amber-auto` 默认 fail-closed：宁可停下并列出具体残基，也不静默产生错误拓扑。下面是常见错误、原因和对应操作。

## 1. 输入分类

### `glycan_unsupported`

```
error: [glycan_unsupported] N-glycans are not supported: NAG A:301, ...
```

- 默认不支持糖链。
- 如果糖链是结晶添加物、确定要删除：`amber-auto run complex.pdb --drop-glycans`（显式、会记录警告）。
- 如果糖链是生物学相关：需要自定义 GLYCAM 参数，当前流程不覆盖。

### `metal_unsupported`

```
error: [metal_unsupported] divalent metals are not supported: ZN A:401
```

- 二价金属可能参与配位，默认不自动删除。
- 确认只是结晶离子：`--drop-metals`。
- 催化/结构金属：需要手动建系并补金属参数。

### `modified_unsupported`

```
error: [modified_unsupported] modified residues are not supported: SEP A:12
```

- 仅 MSE→MET 自动转换。
- 确定要去掉修饰：`--drop-modified`；否则先手动转成标准残基。

### `unknown_hetatm`

```
error: [unknown_hetatm] unknown HETATM residues ...: UNL A:501
```

- 用 `--ligand UNL:A:501` 指定为配体；或 `--drop-unknown-hetatm` 删除。
- 常见结晶添加剂（SO4/GOL/EDO/...）会自动丢弃并记录。

### `ligand_ambiguous` / `multiple_ligand_copies`

```
error: [ligand_ambiguous] multiple plausible ligand candidates: ...
```

- 用 `--ligand RESNAME[:CHAIN[:RESID]]` 选择。
- 同 resname 多拷贝时还要加 `--drop-other-ligands`。
- 共价配体（<1.9 Å）无法用此流程处理。

### `no_interface` / `ambiguous_chains`

- 多链体系无法自动判断受体/伙伴，或两条链不接触。
- 显式指定：

```bash
amber-auto run complex.pdb --type protein-peptide \
  --receptor-chains A --partner-chains B
```

- 链断裂导致的同一链多片段会按 apo/单链处理，并插入 TER。

## 2. 配体补氢与电荷

### RDKit / OpenBabel 不可用

```
error: RDKit is required for automatic ligand charge inference
error: ligand has no hydrogens and OpenBabel could not ...
```

- 在 `config/auto.example.json` 里设置 `rdkit_python`、`obabel`，或用 `--rdkit-python/--obabel` 覆盖。
- 也可以直接给含显式氢的配体 PDB，并用 `--charge` 指定净电荷。

### 电荷/键级模糊

```
error: ligand charge/bond-order assignment is ambiguous ...
error: inferred charge -2 is not confidently separated ...
```

- 查看 `prepare/ligand_charge.json` 的 `candidates`。
- 用 `--charge N` 强制；或 `--ligand-smiles` 给出正确质子化态模板：

```bash
amber-auto run complex.pdb \
  --ligand-smiles 'O=C([O-])c1ccccc1' --charge -1
```

### OpenBabel 推断与 `--charge` 不一致

- OpenBabel 根据几何和 pH 7.4 补氢；如果与你的预期不同，提供含所需氢的配体 PDB。
- `--ph none` 改用普通 `-h` 补氢；`--ph 5.5` 可模拟酸性条件。

## 3. 建系/QC

### `tleap finished with Errors = N` / `ATTN`

- 查看 `systems/<label>/tleap.out` 和 `build_qc.tsv`。
- 常见原因：配体电荷/参数不对、缺失原子、链断裂处理不当、二硫键未配对。
- `--allow-attn` 只用于确认过的无害警告；`amber-auto` 当前不默认暴露该开关，必要时直接调用 `scripts/build_system.sh`。

### 中和失败 / `net charge after neutralization`

- `amber-auto` 会自动换反号 counterion 重试一次。
- 仍失败通常是配体电荷或末端电荷估计错误；查看 `prepare/prepare_manifest.json` 的 `estimated_net_charge`，用 `--charge` 修正配体电荷。

### `prmtop residue linkage check failed`

- 说明 tleap 的残基顺序与 `residue_map.tsv` 不一致，partner mask/MM-GBSA 会错位。
- 检查 `build_linkage.json`；通常是 prepared PDB 含非标准残基或 tleap 重排。

## 4. MD/PBS

### `qsub did not return a job id`

- 检查 qsub 是否可用、队列是否存在、`PBS_O_WORKDIR` 是否可写。
- 用 `amber-auto run ... --dry-run` 查看实际命令和目录。

### 作业一直 `Q`

- `qstat -q` / `pbsnodes -a`；14 张 RTX 3090 是否都在忙。
- `amber-auto status --queue` 会同时显示文件进度和队列。

### 生产被 walltime 杀掉

- `run_prod.pbs` 自带 chunk 跳过；已完成 chunk 不会重算。
- 直接 `amber-auto resume`，它会从 `prod_final.rst7`/`npt.rst7` 判断缺失阶段并重投。
- 减小 `--prod-ns` 或增加 `--prod-chunks` 可缩短单作业时长。

### `MM/GBSA failed`

- 查看 `analysis/mmpbsa/<prefix>/mmpbsa.log`。
- 蛋白 partner 必须使用显式 `ligand_mask/receptor_mask`；`amber-auto` 已自动传入。
- 配体/partner 原子数很大时减少帧数：`--mmpbsa-interval 500`，或 `--no-mmpbsa`。
- 隐式溶剂使用 `--dry-top`；`amber-auto` 自动处理。

## 5. 常用诊断命令

```bash
amber-auto inspect complex.pdb --json
amber-auto run complex.pdb --dry-run
amber-auto status --queue
qstat -f <jobid>
tail -50 amber-campaign/logs/<jobname>.o<jobid>
cat amber-campaign/complexes/<label>/systems/<label>/build_qc.tsv
cat amber-campaign/complexes/<label>/prepare/prepare_manifest.json
```
