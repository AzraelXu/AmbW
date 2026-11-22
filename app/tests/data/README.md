# 测试夹具

这些文件只用于 smoke test，不是科学示例：

- `tiny_peptide.pdb`：tleap 生成的 Ala–Val 二肽（ff19SB）。
- `tiny_rna.pdb`：tleap 生成的 G–C 二核苷酸（RNA.YIL），已去掉 5′ 端 `P/OP1/OP2`，避免 `G5` 模板缺原子类型。
- `ethanol.pdb`：Open Babel 生成的乙醇（显式氢）。
- `lig.mol2` / `lig.frcmod`：`param_ligand.sh` 对乙醇做的 GAFF2/AM1-BCC 参数（resname `LIG`）。

它们用于验证：`pdb4amber`/`tleap` 建系、QC、`residue_masks.py`、`build_system.sh --render-only`，以及 PBS 模板渲染。不用于生产。

`amber-auto` 额外夹具：

- `auto_pl.pdb`：Ala–Val 二肽 + 乙醇配体（`LIG`），蛋白–小分子。
- `auto_pep.pdb`：60 残基 ALA 受体 + 8 残基 ALA 肽，蛋白–肽。
- `auto_antibody.pdb`：两条 100 残基 Ig 样链（含链内二硫键）+ 20 残基抗原，仅用于分类测试。
- `auto_glycan.pdb` / `auto_zn.pdb`：糖链/二价金属错误路径。
- `auto_mse.pdb`：MSE→MET 转换。
- `auto_ss.pdb`：CYS–CYS 二硫键。
- `auto_gap.pdb`：同一链的 C–N 断裂。
- `auto_benzoate_noH.pdb` / `auto_neutral_noH.pdb`：无氢配体电荷推断（−1 / 0）。

这些夹具只验证逻辑和错误路径，几何/化学并不完整，不能用于生产模拟。
