# 高级：amber-kit 旧流程

`amber-auto` 覆盖日常的蛋白–小分子/肽/抗体体系。以下情况继续使用 `amber-kit`：

- RNA–配体、apo RNA、RNA 的 nastruct/堆积分析；
- 多模型 NMR 的模型选择、配体姿态移植；
- 需要逐步检查/修改 `leap.in`、`mask`、每阶段输入；
- 已有按 `config/config.sh` + `config/systems.tsv` 组织的老项目。

## 快速入口

```bash
cp -r ~/Documents/amber_general_workflow ~/Documents/my_project
cd ~/Documents/my_project
bin/amber-kit init
# 编辑 config/config.sh 和 config/systems.tsv
bin/amber-kit check --strict
bin/amber-kit build
bin/amber-kit inputs
bin/amber-kit submit equil --label lig1 --rep 1 --dry-run
bin/amber-kit submit equil --label lig1 --rep 1
# npt.rst7 出现后
bin/amber-kit submit prod --label lig1 --rep 1
bin/amber-kit status --queue
bin/amber-kit analyze --label lig1 --rep 1 --force
bin/amber-kit mmpbsa  --label lig1 --rep 1 --force
bin/amber-kit report
```

完整教程：[`../TUTORIAL.md`](../TUTORIAL.md)。PBS 细节：[`../pbs/README.md`](../pbs/README.md)。

## 与 amber-auto 的关系

- `amber-auto` 生成的每个复合物目录仍保持 `amber-kit` 兼容布局；对 protein-ligand / apo，可以进入该目录用 `amber-kit status/analyze/mmpbsa` 做高级分析。
- `amber-auto` 不修改 `amber-kit` 的任何行为；两套入口可以并存。
- 如果 `amber-auto` 的分类/配体选择不满足需求，用 `amber-kit` 手动建系是受支持的退路。
