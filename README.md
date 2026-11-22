# amber_general_workflow（管理员/维护者文档）

> 使用者请从 [TUTORIAL.md](TUTORIAL.md) 开始。本文只面向部署、运维与二次开发，不包含面向用户的上手步骤。

## 1. 项目定位与边界

`amber_general_workflow` 是一套与具体靶点无关的 Amber 显式/隐式溶剂 MD 工具包，覆盖体系准备、平衡/生产、分析三段。它有两种入口：

| 入口 | 定位 | 输入 | 适用 |
|---|---|---|---|
| `bin/amber-auto` | 推荐的一键流程 | 一个或多个复合物 PDB | 蛋白–小分子、蛋白–肽、蛋白–抗体、蛋白–蛋白、apo 蛋白 |
| `bin/amber-kit` | 旧版项目式流程 | 项目目录 + `config/systems.tsv` | RNA–配体、NMR 多模型、需要精细控制或已有项目 |

**不做**：对接、缺失 loop 建模、N-糖链、共价配体、金属配位、FEP/ABFE、伞形采样、GaMD、REMD、多轨迹 MM/PBSA、熵计算；蛋白–RNA/DNA 复合物需走 `amber-kit` 并自定义 `leap.in`。

本工具包不替代科学判断：配体净电荷/质子化、力场/水模型、分析窗口、MM/GBSA 参数必须由使用者在看到结果前确定。

## 2. 架构与目录结构

### 2.1 入口与分层

```
bin/amber-auto                 # 仓库/共享安装的统一 wrapper（解析 Python 与 site config）
bin/amber-kit                  # 旧流程 wrapper（解析 PYTHON_RENDER）
app/bin/amber-auto             # 主编排器（Python，argparse）
app/bin/amber-kit              # 旧流程编排器（Bash）
app/scripts/*.py|*.sh          # 可独立调用的原子步骤
app/scripts/sch_ui.py          # Python 统一 CLI 输出模块
app/scripts/lib/ui.sh          # Bash 统一 CLI 输出模块
app/docs/CLI_STYLE.md          # 统一输出规范
app/templates/                 # tleap / mdin / cpptraj / MMPBSA 模板
app/pbs/*.pbs.in               # Torque 作业模板
app/config/                    # 默认配置、站点配置、amber-kit 示例
app/tests/                     # smoke test 与检测单测
```

`amber-auto` 的主链路（`app/bin/amber-auto`）：

```
inspect_complex.py
  → prepare_complex.py（pdb4amber/Reduce、配体提取、ligand_charge.py）
  → param_ligand.sh（antechamber AM1-BCC + parmchk2）
  → build_system.sh（tleap + qc_system.py）
  → residue_masks.py / verify_linkage.py
  → setup_md_inputs.py（渲染 min/heat/release/prod）
  → render_template.py（渲染 run_equil/run_prod/run_post PBS）
  → qsub equil → prod → post（afterok 依赖链）
```

### 2.2 仓库目录

```
README.md                      # 本文（管理员）
TUTORIAL.md                    # 用户教程
README-SHARED.md               # /opt 共享安装的部署附录
VERSION                        # 共享安装版本元数据
config/.gitkeep                # 部署时可放 auto.site.json
app/README.md                  # app 包索引（指向根文档）
app/TUTORIAL.md                # 旧版 amber-kit 高级教程
app/docs/                      # AUTO/ERRORS/TROUBLESHOOTING/CHEATSHEET 等
app/examples/                  # amber-kit 示例项目
app/logs/                      # 旧流程日志占位
```

### 2.3 运行时 campaign 目录

`amber-auto` 的默认 campaign 为 `./amber-campaign`，可用 `--campaign DIR` 指定。布局见第 6 节。

## 3. 依赖与安装部署

### 3.1 依赖矩阵

| 组件 | 用途 | 站点默认 |
|---|---|---|
| Amber 22 / AmberTools 22 | `tleap`、`pdb4amber`、`antechamber`、`parmchk2`、`pmemd.cuda`、`sander.MPI`、`cpptraj`、`MMPBSA.py` | `/opt/software/amber22` |
| conda env `qc` | parmed、numpy；建系 QC、masks、linkage | `/opt/software/amber_general_workflow/envs/qc` |
| conda env `analysis` | RDKit、MDAnalysis、scipy、netCDF4、pandas；电荷推断、分析、汇总 | `/opt/software/amber_general_workflow/envs/analysis` |
| OpenBabel | 无氢配体补氢/键级感知 | `/opt/schrodinger2021-2/utilities/obabel` |
| Torque 6.1 | `qsub`/`qstat`，队列 `batch` | `/opt/software/torque-6.1/bin` |
| GPU | `pmemd.cuda` | node02–node08，各 4×RTX 3090 |

`app/templates` 与 `app/scripts` 不依赖第三方 Python 包；只有 QC/分析脚本依赖上述两个 env。

### 3.2 源码 checkout 运行

当前 checkout 不含 `envs/`（`.gitignore` 已忽略），但 wrapper 会按以下顺序解析：

1. `AMBER_AUTO_PYTHON`（命令名或绝对路径）；
2. `$HERE/envs/analysis/bin/python`；
3. `/opt/software/amber_general_workflow/envs/analysis/bin/python`；
4. `command -v python3`。

`bin/amber-kit` 同理，变量名为 `PYTHON_RENDER`。因此从源码目录直接运行 `./bin/amber-auto` 时，会复用共享安装的 analysis env。

### 3.3 `/opt` 共享安装

共享安装由管理员维护，目录 root 所有、只读：

```
/opt/software/amber_general_workflow/
├── bin/{amber-auto,amber-kit}
├── app/
├── config/                    # 站点配置（可能为空；wrapper 会回退到 app/config/auto.site.json）
├── envs/{qc,analysis}/
├── README-SHARED.md
└── VERSION
```

注意：共享安装的 `config/auto.site.json` 可能不存在；`app/bin/amber-auto` 的 `load_config()` 会回退到 `app/config/auto.example.json` 或内置 `DEFAULTS`。本仓库的 `auto.example.json` 与 `auto.site.json` 当前内容一致，且为真实 `/opt` 绝对路径。

### 3.4 从源码部署到共享安装

```bash
sudo rsync -a --delete \
  --exclude '__pycache__' --exclude 'amber-campaign' --exclude '.git' \
  /home/user1/src/amber_general_workflow/app/ /opt/software/amber_general_workflow/app/
sudo rsync -a /home/user1/src/amber_general_workflow/bin/ /opt/software/amber_general_workflow/bin/
sudo chown -R root:root /opt/software/amber_general_workflow
sudo chmod -R a+rX /opt/software/amber_general_workflow
sudo chmod -R go-w /opt/software/amber_general_workflow
```

部署后必须保留 `app/bin/amber-auto` 的 shebang 能解析到 analysis env，或依赖 wrapper 的显式解释器。不要用 `--delete` 删除 `envs/`、`config/`、`VERSION`、`README-SHARED.md`。

移植到其他集群时，优先改 site config 或环境变量：`amber_sh`、`rdkit_python`、`py_qc`、`py_analysis`、`obabel`、`pbs_queue`、`pbs_nodes`、`ppn`、`gpus`、`equil_walltime`、`prod_walltime`、`post_walltime`。PBS 头语法与 GPU 绑定逻辑见 [app/pbs/README.md](app/pbs/README.md) 第 6 节。

## 4. 环境与配置项

### 4.1 解析优先级

`amber-auto` 的配置来源从高到低：

1. 命令行参数（如 `--prod-ns`、`--pbs-nodes`）；
2. `--config PATH`；
3. `AMBER_AUTO_CONFIG` 环境变量；
4. wrapper 自动选择的真实 site config：`<repo>/config/auto.site.json` → `<repo>/app/config/auto.site.json`；
5. `app/config/auto.example.json`；
6. 内置 `DEFAULTS`（`app/bin/amber-auto`）。

wrapper **不会**把 `auto.example.json` 作为自动回退项；若没有真实 site config，就让 `load_config()` 走它自己的既有默认。

### 4.2 关键配置项

`app/config/auto.example.json` 与 `DEFAULTS` 中的键：

| 键 | 默认 | 说明 |
|---|---|---|
| `amber_sh` | `/opt/software/amber22/amber.sh` | Amber 环境脚本 |
| `rdkit_python` | analysis env | 配体电荷推断 |
| `py_qc` | qc env | 建系 QC/masks/linkage |
| `py_analysis` | analysis env | 轨迹分析/汇总 |
| `obabel` | Schrödinger OpenBabel | 无氢配体补氢 |
| `campaign` | `./amber-campaign` | 输出根 |
| `ff` / `water` / `solvent` | `ff19SB` / `opc` / `explicit` | 力场/水/溶剂 |
| `box` / `pad` | `oct` / `10.0` | 盒型/缓冲 |
| `cation` / `anion` | `Na+` / `Cl-` | 中和离子 |
| `radii` / `igb` | `mbondi2` / `5` | GB 参数配对 |
| `heat_ns` / `release_stages` / `release_ns` | `2.0` / `10` / `1.0` | 平衡协议 |
| `prod_ns` / `prod_chunks` | `100.0` / `4` | 生产长度/分块 |
| `dt` / `traj_interval_ps` | `0.002` / `10.0` | 步长/轨迹间隔 |
| `replicas` / `seed_base` | `1` / `20261001` | 副本/随机种子 |
| `strip_mask` | `:WAT,Na+,Cl-,K+` | 分析剥离 |
| `ph` | `7.4` | 配体补氢 pH |
| `no_reduce` | `false` | 跳过 Reduce |
| `mmpbsa` / `mmpbsa_method` / `mmpbsa_igb` / `mmpbsa_saltcon` / `mmpbsa_interval` | `true` / `gb` / `5` / `0.12` / `100` | MM/GBSA |
| `param_pbs` / `param_ppn` / `param_walltime` / `param_amber_sh` / `param_nodes` | `false`（site 为 `true`）/ `40` / `04:00:00` / `null`（site 为 conda ABT 的 `amber.sh`）/ `null`（site 为 `1`） | 配体 AM1-BCC 参数化是否提交 PBS、申请 CPU 数、walltime、AmberTools 环境、节点 |
| `min_engine` | `sander` | 最小化引擎；`pmemd.cuda` 为实验选项，未通过本示例验证（见下） |
| `pbs_queue` / `pbs_nodes` / `ppn` / `gpus` | `batch` / `1` / `8` / `1` | Torque 资源 |
| `equil_walltime` / `prod_walltime` / `post_walltime` | `24:00:00` / `24:00:00` / `12:00:00` | walltime |

`--ff ff14SB` 单独使用时自动切 `water=tip3p`；`--igb 8` 单独使用时自动切 `radii=mbondi3`。

### 4.3 环境变量

| 变量 | 作用 |
|---|---|
| `AMBER_AUTO_PYTHON` | wrapper 解释器（命令名或绝对路径） |
| `AMBER_AUTO_CONFIG` | site config 路径；低于 `--config`，高于自动回退 |
| `AMBER_AUTO_QSUB` / `AMBER_AUTO_QSTAT` | 覆盖 `qsub`/`qstat` 命令（测试或换调度器） |
| `AMBER_MD_AMBER_SH` | 子脚本使用的 amber.sh；由 `amber-auto` 自动注入 |
| `AMBER_MD_PYTHON` | 子脚本 QC Python；由 `amber-auto` 从 `py_qc` 注入，且仅在可执行时注入 |
| `PYTHON_RENDER` | `amber-kit` 的模板渲染 Python |
| `SCH_COLOR` | 颜色模式 `auto`/`always`/`never`；子进程与 PBS 日志固定 `never` |

所有入口支持 `--color=auto|always|never`，优先级高于 `SCH_COLOR`。`auto` 只在对应流是 TTY 时着色；重定向到日志或管道会自动降级为纯文本。`--json`、JSON/TSV 报告、PBS 日志、生成的配置与模板始终不含 ANSI。

### 4.4 配对约束

- `antechamber -at gaff2` ⇔ `parmchk2 -s gaff2` ⇔ `leaprc.gaff2`，三者必须一致。
- `igb=5 ↔ mbondi2`，`igb=8 ↔ mbondi3`；改 `igb` 必须重建拓扑。
- `tleap` 必须 `Errors = 0`；`ATTN` 默认失败，只有化学审查后才允许放行。
- 显式溶剂用 `solvateOct/solvateBox` + `addIons`；隐式溶剂用 `igb`/`saltcon`，不加水。

## 5. 批量调度与资源管理

### 5.1 批量入口语义

`amber-auto run pdb1 pdb2 ... pdbN` 是唯一的原生多体系入口：

1. **串行**建系：逐个 `inspect → prepare → param → tleap → QC → masks → inputs → render`；
2. 全部建系完成后，才进入提交阶段；
3. 每个体系每个副本渲染三个 PBS 脚本，提交为 `equil → prod → post` 的 `afterok` 依赖链；
4. 不同体系/副本的链彼此独立，由 Torque 决定并发。

`--no-submit`：只建系/渲染，不 qsub；`--dry-run`：建系后打印提交计划并写 `submission_plan.json`，不 qsub。`--force` 会重建指定体系（删除该体系的 `prepare/`、`systems/<label>/`、`ligands/<resname>/`、`md/<label>/` 与 analysis 子目录），不会动其他体系。

配体 AM1-BCC 参数化默认可在计算节点以 PBS 作业运行（`param_pbs`，站点默认开启）：`amber-auto` 渲染 `run_param_ligand.pbs`，用 `param_ppn`（默认 40）个 CPU 提交，登录节点轮询等待 `param_status.txt`，成功后继续 tleap 建系。若没有 PBS，使用 `--no-param-pbs` 回退到本地运行。

`param_amber_sh` 决定参数化使用的 AmberTools 环境，默认回退到 `amber_sh`。**关键差异**：`/opt/software/amber22` 的独立 `sqm` 链接静态 BLAS，是串行的（实测 1 核）；conda `ambertools_env` 的 `sqm` 链接 OpenBLAS pthreads（`libopenblasp-r0.3.27.so`），设置 `OPENBLAS_NUM_THREADS=40` 后 BLAS 部分可用满 40 核。本站点 `auto.site.json` 已把 `param_amber_sh` 指向 conda ABT，因此 `param_ppn=40` 是实际有效的。PBS 作业日志会打印 `sqm` 路径、`OPENBLAS_NUM_THREADS` 与 BLAS 链接，便于核对。

最小化默认用 `sander.MPI`。`min_engine=pmemd.cuda` 时会渲染 `run_equil_gpu_min.pbs.in`，min1/min2 也走 GPU（`pmemd.cuda` 支持 `imin=1`）。但**该模式在本示例 86,964 atoms 体系上实测发散**：`01_min_restraint` 的 `BOND` 从 222 升到 69,496 kcal/mol，`02_min_free` 继承坏结构。因此站点默认保持 `sander`；`pmemd.cuda` 仅作为实验选项，使用前必须检查 `BOND`、`FINAL RESULTS` 和下一阶段能量，确认未发散。CPU 最小化较慢，但稳定；24h equil walltime 应据此评估。

### 5.2 部分失败策略（fail-closed）

- 任一体系建系失败且未给 `--submit-partial`：**不提交任何作业**；
- 给了 `--submit-partial`：只提交建系成功的体系；
- 退出码：部分失败 `5`，全部失败 `3`，提交失败 `5`（以 `cmd_build_or_run` 的实现为准）。

`build_one` 捕获 `BuildError` 与未预期异常，记录 `unexpected <Type>: <msg>` 并把该体系标记为 `failed`，然后继续处理下一个体系；不会让一个坏 PDB 中断整批建系。

### 5.3 GPU/节点

- 默认 `nodes=1:ppn=8:gpus=1`，每个 MD 作业一张 GPU；
- `--pbs-nodes node02` 会把节点名写进渲染后的 PBS（建系时固化，`resume` 不会改）；
- 模板优先使用 `PBS_GPUFILE` 推导 `CUDA_VISIBLE_DEVICES`，若 Torque 分配多张则取第一张；
- **node03 已知会让 `sander.MPI` 卡在 `01_min_restraint`**，运维上建议固定 node02 或 node04–node08；
- `run_two_gpu` / `run_multi_gpu` 是 `amber-kit` 的整节点多副本模板，`amber-auto` 不走这两条路径。

### 5.4 并发与节流现状

`amber-auto` 没有并发上限、作业数组或节流参数；并发度由 Torque 队列和节点资源决定。批量建系本身是串行的。若需要控制并发，应在站点层配置队列限制，或在外部编排器里分批调用。

## 6. 数据与输出布局

```
<campaign>/
  campaign.json                     # 状态、job id、依赖、warnings、config 快照
  .campaign.lock                    # 状态写入 flock
  submission_plan.json              # run --dry-run
  resume_plan.json                  # resume --dry-run
  inspections/<label>.json
  complexes/<label>/
    inputs/original.pdb             # 输入副本；源文件不被修改
    components.json
    prepare/                        # protein/ligand 准备产物、Reduce 日志、二硫键
    ligands/<RESNAME>/              # lig.mol2、lig.frcmod、antechamber/parmchk2 日志
    systems/<label>/                # prmtop/inpcrd、masks.env、QC、linkage
    md/<label>/r<rep>/              # MD 输入、run_*.pbs、轨迹、restart、锁文件
    analysis/                       # cpptraj、python npz/json、mmpbsa、tables
    logs/                           # 每作业自写日志 <jobname>.o<jobid>
  report/{summary.tsv,report.json,report.md}
  logs/
```

关键完成证据：

- 每个 stage：`<stage>.rst7` 非空 + `<stage>.out` 含 `Total wall time` + NSTEP ≥ nstlim；
- `npt.rst7`、`prod_final.rst7` 为最后一段 restart 的副本；
- `analysis/post_status_<prefix>.json` 的 `steps`（0=成功，1=失败，2=跳过）；
- `.run_equil.lock` / `.run_prod.lock` / `.run_prod_cNN.lock` 防止同副本并发。

## 7. 运维与故障排查

### 7.1 自检

```bash
./bin/amber-auto --help
/opt/software/amber_general_workflow/envs/qc/bin/python -I -c 'import parmed,numpy'
/opt/software/amber_general_workflow/envs/analysis/bin/python -I -c 'import rdkit,MDAnalysis,netCDF4'
source /opt/software/amber22/amber.sh
command -v tleap pdb4amber antechamber parmchk2 cpptraj pmemd.cuda sander.MPI MMPBSA.py
qstat -q
pbsnodes -a node02 | grep -E 'state =|np =|gpus =|availmem'
SKIP_BUILD=1 bash app/tests/auto_smoke.sh
```

### 7.2 作业与续跑

```bash
./bin/amber-auto status --campaign DIR --queue
qstat -u "$USER"
./bin/amber-auto resume --campaign DIR --dry-run   # 只写 resume_plan.json，不改 campaign 状态
./bin/amber-auto resume --campaign DIR             # 按文件证据从 equil/prod/post 续跑
```

`resume` 会跳过 qstat 中仍处于 `Q/R/H/W/S` 的作业；没有 `npt.rst7` 从 `equil` 重投，有 `npt.rst7` 无 `prod_final.rst7` 从 `prod` 重投，两者都有但 post 未完成则从 `post` 重投。`--dry-run` 不修改 `campaign.json`。

### 7.3 常见故障

| 现象 | 根因 | 处理 |
|---|---|---|
| wrapper 报找不到 Python | `envs/` 不存在且 `AMBER_AUTO_PYTHON` 未设 | 设 `AMBER_AUTO_PYTHON`，或安装共享 env |
| QC 报 parmed 缺失 | `py_qc` 不可执行，`pyenv.sh` 回退到 PATH python | 检查 `py_qc`；`amber-auto` 会注入 `AMBER_MD_PYTHON` |
| tleap `Errors != 0` | 缺参数/残基名不合法/电荷异常 | 查 `systems/<label>/tleap.out` 与 `build_qc.tsv` |
| tleap `ATTN > 0` | 参数近似或化学不合理 | 逐条审查；不要默认放行 |
| `Small box detected` | 盒子太小（cut=10 Å 时约 <35–40 Å） | `--pad 15` 或更大，重建 |
| node03 卡在 `01_min_restraint` | 已知节点问题 | `--pbs-nodes node02`（建系时固化） |
| 作业被杀 | walltime 不足/节点故障 | `resume`；检查 `#PBS -l walltime` |
| 续跑后重复建系报 duplicate | 同一 campaign 再次 `run` 同一 PDB | 用 `resume`，或新建 campaign，或 `--force` |
| dry-run 后状态异常 | 旧版本会写 DRY job | 本版本 dry-run 不写 campaign 状态 |

## 8. 扩展开发规范

- **脚本边界**：每个 `app/scripts/*` 应可独立运行，提供 `--help`、`--force`、JSON/TSV 结构化输出；失败用非零退出码与明确 stderr。
- **错误处理**：不要 `except: pass`；批处理边界记录异常并继续，原子步骤快速失败。
- **模板**：新增 stage 时同时更新 `app/templates/md/...`、`setup_md_inputs.py` 的 stage 表和 `run_equil/run_prod.pbs.in` 的 `stage_done` 逻辑。
- **QC 门槛**：`tleap Errors=0`、`ATTN=0`（默认）、净电荷 0、`RADIUS_SET` 匹配、mask 推导通过。新增体系类型必须给出对应的 QC 断言。
- **向后兼容**：新增能力用可选参数暴露；不改已有参数默认值；不改 JSON schema 的既有键。
- **测试**：
  - `bash app/tests/smoke.sh`：语法、tleap 渲染/小型真实建系、mask、PBS 模板渲染、通用性 grep；
  - `bash app/tests/auto_smoke.sh`：`amber-auto` 检测、电荷推断、MSE/CYX、glycan 拒绝、真实 build（可 `SKIP_BUILD=1`）；
  - `python app/tests/test_detection.py`：`inspect_complex.py` 单测。
- **输出格式**：统一使用 `app/scripts/sch_ui.py` / `app/scripts/lib/ui.sh`，不要在业务代码里写 ANSI；FATAL/ERROR/WARNING 走 stderr，OK/INFO/HEADER 走 stdout；日志、JSON/TSV、模板和生成配置保持纯文本；新脚本在 argparse 前调用 `sch_ui.extract_color_arg` / `ui_parse_color_args`。完整规范见 [app/docs/CLI_STYLE.md](app/docs/CLI_STYLE.md)。
- **文档**：改动配置/CLI 时同步更新根 `README.md`（管理员）与 `TUTORIAL.md`（用户），不要在 `app/README.md` 重复教程。

## 9. 变更记录

### 2026-11-22

- 修复仓库入口：`bin/amber-auto` / `bin/amber-kit` 增加 Python 与 site config 回退；不再把 `auto.example.json` 作为自动回退项。
- 修复 QC 环境：`amber-auto` 向子脚本注入 `AMBER_MD_PYTHON=$py_qc`（仅当可执行），避免 `build_system.sh` 回退到无 parmed 的 PATH python。
- 批量失败隔离：`build_one` / `cmd_build_or_run` 捕获未预期异常，单体系失败不再中断整批；保持 fail-closed 提交策略。
- 修复 `resume --dry-run`：dry-run 不再写 `campaign.json` 的 job/status，只写 `resume_plan.json`。
- 配体参数化支持 PBS：新增 `param_pbs` / `param_ppn` / `param_walltime` / `param_amber_sh` / `param_nodes` 与对应 CLI；`run_param_ligand.pbs.in` 增加 `OMP/MKL/OPENBLAS_NUM_THREADS` 与 `param_status.txt`。站点默认使用 conda ABT 的 OpenBLAS sqm，40 核实测有效。
- 新增实验性 GPU 最小化：`min_engine` 与 `--min-engine`；新增 `app/pbs/run_equil_gpu_min.pbs.in`，min1/min2 走 `pmemd.cuda`。**默认仍为 `sander`**；本示例 87k 体系实测 `BOND` 发散，未启用。
- 统一终端输出：新增 `app/scripts/sch_ui.py` / `app/scripts/lib/ui.sh` 与 `app/docs/CLI_STYLE.md`；`--color=auto|always|never`，日志/JSON/PBS 保持纯文本。
- 文档拆分：新增根 `README.md`（管理员）与根 `TUTORIAL.md`（用户）；`app/README.md` 收缩为索引。

## 10. 文档地图

| 文档 | 读者 | 内容 |
|---|---|---|
| [TUTORIAL.md](TUTORIAL.md) | 使用者 | 输入准备、分步操作、单体系/批量示例、参数、结果校验、排错、FAQ |
| [README-SHARED.md](README-SHARED.md) | 共享安装管理员 | `/opt` 安装布局、权限、更新流程、集群注意事项 |
| [app/TUTORIAL.md](app/TUTORIAL.md) | 高级用户 | `amber-kit` 旧流程完整教程（RNA、NMR 多模型、精细控制） |
| [app/docs/AUTO.md](app/docs/AUTO.md) | 使用者/维护者 | `amber-auto` 设计细节与能力边界 |
| [app/docs/ERRORS.md](app/docs/ERRORS.md) | 使用者 | `amber-auto` 错误码与逃生口 |
| [app/docs/TROUBLESHOOTING.md](app/docs/TROUBLESHOOTING.md) | 使用者/维护者 | 症状→原因→修复 |
| [app/docs/CHEATSHEET.md](app/docs/CHEATSHEET.md) | 使用者 | 命令速查 |
| [app/docs/CLI_STYLE.md](app/docs/CLI_STYLE.md) | 维护者 | 统一终端输出规范（颜色/流/日志/Next steps） |
| [app/pbs/README.md](app/pbs/README.md) | 维护者 | PBS 前缀逐行说明与移植 |
| [app/config/README.md](app/config/README.md) | 维护者 | 配置字段与 `systems.tsv` 语义 |
