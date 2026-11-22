# RNA–配体示例

1. 复制 `config.sh` 到项目 `config/config.sh`。
2. 把 RNA（或 RNA+配体复合物）PDB 放到 `inputs/receptor.pdb`，配体 PDB（显式氢、单残基）放到 `inputs/lig1.pdb`。
3. 填 `LIGAND_CHARGE`，编辑 `systems.tsv`。
4. 运行：

```bash
bin/amber-kit check
bin/amber-kit build
bin/amber-kit inputs
bin/amber-kit submit equil --label lig1 --rep 1 --dry-run
```

要点：

- `RNA.YIL` 通常配 `TIP3P`；`K+` 中和；默认不加额外盐、不加 Mg²⁺。
- 双链 RNA 用 `ANCHOR_MASK`/`FIXED_MASK`（例如 `:1-13` / `:14-26,27`）做 `autoimage`。
- `RNA_RESRANGE` 只覆盖 RNA，不含配体/离子；`LOOP_RESIDUES` 用于 loop–配体接触。
- `LIGAND_RING_ATOMS`/`STACK_TARGETS` 用于堆积分析；阈值在看结果前定好。
- 5′ 端带磷酸的 PDB 可能触发 `G5 ... OP1 does not have a type`；去掉端基 `P/OP1/OP2`。
