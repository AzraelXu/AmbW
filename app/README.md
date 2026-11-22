# app 包索引

本目录是 `amber_general_workflow` 的代码、脚本、模板、配置与测试。**用户教程已迁移到仓库根目录 [TUTORIAL.md](../TUTORIAL.md)，管理员/部署文档在根目录 [README.md](../README.md)。**

## 文档地图

| 文档 | 读者 | 内容 |
|---|---|---|
| [../TUTORIAL.md](../TUTORIAL.md) | 使用者 | 输入准备、分步操作、单体系/批量示例、参数、结果校验、排错、FAQ |
| [../README.md](../README.md) | 管理员/维护者 | 架构、依赖与部署、配置项、批调度、输出布局、运维、扩展规范 |
| [../README-SHARED.md](../README-SHARED.md) | 共享安装管理员 | `/opt` 安装布局、权限、更新流程、集群注意事项 |
| [TUTORIAL.md](TUTORIAL.md) | 高级用户 | `amber-kit` 旧流程完整教程（RNA、NMR 多模型、精细控制） |
| [docs/AUTO.md](docs/AUTO.md) | 使用者/维护者 | `amber-auto` 设计细节与能力边界 |
| [docs/ERRORS.md](docs/ERRORS.md) | 使用者 | `amber-auto` 错误码与逃生口 |
| [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) | 使用者/维护者 | 症状→原因→修复 |
| [docs/CHEATSHEET.md](docs/CHEATSHEET.md) | 使用者 | 命令速查 |
| [docs/CLI_STYLE.md](docs/CLI_STYLE.md) | 维护者 | 统一终端输出规范（颜色/流/日志/Next steps） |
| [pbs/README.md](pbs/README.md) | 维护者 | PBS 前缀逐行说明与移植 |
| [config/README.md](config/README.md) | 维护者 | 配置字段与 `systems.tsv` 语义 |

## 本目录结构

```
bin/{amber-auto,amber-kit}     # 编排器
scripts/                       # 可独立调用的原子步骤
templates/                     # tleap / mdin / cpptraj / MMPBSA 模板
pbs/                           # Torque 作业模板
config/                        # 默认配置、站点配置、amber-kit 示例
docs/                          # 设计与排错文档
examples/                      # amber-kit 示例项目
tests/                         # smoke test 与检测单测
logs/                          # 旧流程日志占位
```

## 能力摘要

- `amber-auto`：复合物 PDB → 分类/准备/配体参数化/tleap+QC/masks/MD 输入 → `equil → prod → post` Torque 依赖链 → cpptraj/分析/MM-GBSA。
- `amber-kit`：项目式旧流程，支持 RNA–配体、NMR 多模型、已有项目与精细控制。

## 科学判断警告

本工具包不替代科学判断：配体净电荷/质子化、力场/水模型、分析窗口、MM/GBSA 参数都必须由使用者在看到结果前确定。八条核心原则见 [../TUTORIAL.md](../TUTORIAL.md) 第 7.6 节。
