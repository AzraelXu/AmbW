# 蛋白–配体示例

1. 复制 `config.sh` 到项目 `config/config.sh`。
2. 把蛋白 PDB 放到 `inputs/receptor.pdb`，配体 PDB（显式氢、单残基）放到 `inputs/lig1.pdb`。
3. 填 `LIGAND_CHARGE`，编辑 `systems.tsv`。
4. 运行：

```bash
bin/amber-kit check
bin/amber-kit build
bin/amber-kit inputs
bin/amber-kit submit equil --label lig1 --rep 1 --dry-run
```

要点：

- `ff19SB` 推荐 `OPC`；`ff14SB` 通常配 `TIP3P`。
- 中性蛋白 + 中性配体可用 `COUNTERION=none`；带电体系按净电荷方向选 `cation`/`anion`。
- `IGB=5` 配 `mbondi2`；若改 `IGB=8`，必须改 `RADII=mbondi3` 并重建拓扑。
- 配体电荷/质子化必须手填，脚本不猜。
