# amber_general_workflow 使用教程

> 部署、配置与二次开发请看 [README.md](README.md)。本文只讲怎么用。

本教程以 `amber-auto` 为主；`amber-kit` 旧流程放在最后一节并链接完整教程。

## 1. 前置条件

- 能登录集群并在 `batch` 队列提交 Torque 作业（`qsub`/`qstat` 可用）。
- 有可写的项目目录；轨迹和日志建议放在 `/home`（本机约 7 TB 可用）。
- 不需要手动 `source amber.sh`，也不需要 `conda activate`；入口脚本会自己处理。
- 推荐直接使用共享安装入口：

```bash
export AMBER_AUTO=/opt/software/amber_general_workflow/bin/amber-auto
$AMBER_AUTO --help
```

从源码 checkout 运行时把 `$AMBER_AUTO` 换成 `/home/user1/src/amber_general_workflow/bin/amber-auto`；wrapper 会自动使用共享安装的 analysis 环境。

所有命令支持 `--color=auto|always|never`，默认 `auto`；输出重定向到日志/管道时自动使用纯文本。JSON、TSV、PBS 日志与生成的配置文件永远不含 ANSI 颜色码。

## 2. 输入准备与目录规范

### 2.1 输入是什么

每个体系只需要一个**复合物 PDB**：

- 蛋白链（一条或多条）；
- 小分子配体作为 HETATM 残基嵌在复合物里（推荐单残基、单拷贝、带显式氢）；
- 不要预先加水、离子或去氢；`amber-auto` 会自己处理。

目录可以是平铺的：

```
my_project/
├── HLD0117_WT.pdb
└── HLD0117_Y537S.pdb
```

也可以是每个体系一个子目录；批量运行时用 shell glob 指定即可，脚本不强制目录层级。

### 2.2 命名与配体

- 文件名会成为 campaign 里的 label（`HLD0117_WT.pdb` → `HLD0117_WT`），建议只用字母、数字、`_`、`-`、`.`。
- 配体 resname 可以是常规三字母名，也可以是本示例的 `" 01"`；`inspect` 会归一化为 `01`。
- 配体带显式氢时，`amber-auto` 会保留其质子化状态并用 RDKit 枚举键级/电荷；无氢时用 OpenBabel 按 `--ph`（默认 7.4）补氢，再用 RDKit 校验。
- **配体净电荷必须由你确认**。自动枚举给出候选后，建议用 `--charge` 显式钉住；电荷有歧义时工具会硬报错，不会静默猜。

本教程示例配体为 `C51H55F3N7O3S+`（质子化哌嗪，净电荷 **+1**），因此所有命令都带 `--charge 1`。

### 2.3 先看分类

```bash
$AMBER_AUTO inspect HLD0117_WT.pdb HLD0117_Y537S.pdb
```

重点看：

- `type` 应为 `protein-ligand`；
- `ligand` 是否为预期残基（本示例为 `01 L:1`）；
- 是否有 `errors` 或 `warnings`。

需要 JSON 时加 `--json`；需要把检查结果落到指定目录时加 `--campaign DIR`。

## 3. 从零跑通

### 3.1 推荐节奏

```bash
# 1) 分类
$AMBER_AUTO inspect complex.pdb

# 2) 只建系 + 渲染 PBS，不提交；先检查 QC 和输入
$AMBER_AUTO run complex.pdb --no-submit --campaign ./my-campaign

# 3) 通过检查后提交 equil -> prod -> post
$AMBER_AUTO resume --campaign ./my-campaign

# 4) 进度 / 汇总 / 断点续跑
$AMBER_AUTO status --campaign ./my-campaign --queue
$AMBER_AUTO report --campaign ./my-campaign
$AMBER_AUTO resume --campaign ./my-campaign
```

也可以一条命令直接建系并提交：

```bash
$AMBER_AUTO run complex.pdb --campaign ./my-campaign
```

`run` 会先串行建系；只有全部建系成功才提交。任一体系建系失败时默认**一个作业都不提交**，修好后重跑；确实只想推进成功体系时才加 `--submit-partial`。

### 3.2 demo 模式

`--dry-run` 会真实建系、渲染 PBS、写出 `submission_plan.json`，但**不 qsub**；`--no-submit` 只建系/渲染。两者都不会占用 GPU。

## 4. 单体系完整示例

以 `/home/user1/Documents/ER_BRD4/HLD-0117/amber_test/HLD0117_WT.pdb` 为例：

```bash
export AMBER_AUTO=/opt/software/amber_general_workflow/bin/amber-auto
WT=/home/user1/Documents/ER_BRD4/HLD-0117/amber_test/HLD0117_WT.pdb
RUN=/home/user1/Documents/ER_BRD4/HLD-0117/amber_runs_demo
mkdir -p "$RUN/logs"

# 分类
$AMBER_AUTO inspect "$WT" --campaign "$RUN/inspect"

# 建系 + 渲染，不提交（先过 QC）
$AMBER_AUTO run "$WT" --campaign "$RUN/wt" \
  --ligand 01:L:1 --charge 1 \
  --prod-ns 10 --prod-chunks 4 --mmpbsa-interval 20 \
  --pbs-nodes node02 --no-submit \
  2>&1 | tee "$RUN/logs/wt_build.log"

# 检查通过后提交
$AMBER_AUTO resume --campaign "$RUN/wt" 2>&1 | tee "$RUN/logs/wt_submit.log"

# 轮询
$AMBER_AUTO status --campaign "$RUN/wt" --queue
```

默认平衡段是 2 ns 加热 + 10×1 ns NPT 释放；生产 10 ns 分成 4 个 2.5 ns chunk。单副本总 MD 约 22 ns。10 ns 生产适合流程验证；科学结论请按体系需要加长 `--prod-ns` 和/或 `--replicas`。

## 5. 批量完整示例

对 `amber_test` 下全部两个体系，一条命令建系并提交：

```bash
export AMBER_AUTO=/opt/software/amber_general_workflow/bin/amber-auto
DATA=/home/user1/Documents/ER_BRD4/HLD-0117/amber_test
RUN=/home/user1/Documents/ER_BRD4/HLD-0117/amber_runs_20261122
CAMP=$RUN/campaign_10ns
mkdir -p "$RUN/logs"

# 1) 建系 + 渲染（不提交）
$AMBER_AUTO run "$DATA"/HLD0117_WT.pdb "$DATA"/HLD0117_Y537S.pdb \
  --campaign "$CAMP" \
  --ligand 01:L:1 --charge 1 \
  --prod-ns 300 --prod-chunks 12 \
  --prod-walltime 72:00:00 --post-walltime 24:00:00 \
  --mmpbsa-interval 100 \
  --pbs-nodes node02 --no-submit \
  2>&1 | tee "$RUN/logs/full_build.log"

# 2) 检查 campaign.json / build_manifest.json / masks.env / run_config.json

# 3) 提交
$AMBER_AUTO resume --campaign "$CAMP" 2>&1 | tee "$RUN/logs/full_submit.log"

# 4) 监控
$AMBER_AUTO status --campaign "$CAMP" --queue
qstat -u "$USER"
tail -n 50 "$CAMP"/complexes/*/logs/*.o* 2>/dev/null
```

两个体系会各自提交一条 `equil → prod → post` 链；它们可以并行，具体取决于 Torque 资源和节点占用。本示例固定 `node02` 以避开已知有问题的 node03；node02 有 4 张 GPU，两个单 GPU 作业可以并存。

### 5.1 日志与产物

- 建系/提交的终端输出：`$RUN/logs/*.log`（上面的 `tee`）；
- 每个 PBS 作业：`$CAMP/complexes/<label>/logs/<jobname>.o<jobid>`；
- 配体参数化：`$CAMP/complexes/<label>/ligands/LIG/{antechamber.log,parmchk2.log}`；
- tleap：`$CAMP/complexes/<label>/systems/<label>/tleap.out`；
- 轨迹/restart：`$CAMP/complexes/<label>/md/<label>/r1/`；
- 分析/MM-GBSA：`$CAMP/complexes/<label>/analysis/`；
- 汇总：`$CAMP/report/summary.tsv`。

### 5.2 耗时预期

以本示例规模（蛋白 354 残基 + 配体 + 显式水，实际原子数见 `build_qc.tsv` 或 prmtop `NATOM`），RTX 3090 上大致：

- 平衡 12 ns：约 2–6 h；
- 生产 300 ns（12×25 ns chunk）：约 1.5–3 天；
- post（cpptraj + 分析 + 约 150 帧 MM/GBSA）：约 2–8 h；
- 两体系并行时总墙钟约 2–4 天，以 `status` 的 `ns/day` 为准。

`run_prod.pbs` 在一个 72 h 作业里顺序跑完所有 chunk；如果作业被 walltime 或节点故障中断，直接 `resume`，已完成的 chunk 会保留，从第一个未完成 chunk 继续。

## 6. 关键参数说明

| 参数 | 作用 | 建议 |
|---|---|---|
| `--type` | 覆盖体系分类 | 自动识别不确定时显式指定 |
| `--ligand RESNAME[:CHAIN[:RESID]]` | 指定配体 | 本示例用 `01:L:1`；多候选/多拷贝时必填 |
| `--charge INT` | 配体净电荷 | **建议总是显式给出**；本示例 `1` |
| `--ligand-smiles SMILES` | 用模板校验键级/电荷 | 需要锁定质子化/互变异构时用 |
| `--solvent explicit\|implicit` | 显式/隐式溶剂 | 生产用 explicit；快速筛查可用 implicit |
| `--ff ff19SB\|ff14SB` | 蛋白力场 | 默认 ff19SB；配体固定 GAFF2 |
| `--water opc\|tip3p` | 水模型 | ff19SB 配 OPC；ff14SB 自动配 TIP3P |
| `--box oct\|rect` / `--pad` | 盒型/缓冲 | 球状蛋白 oct+10；细长/多链考虑 rect |
| `--heat-ns` / `--release-stages` / `--release-ns` | 平衡协议 | 默认 2 ns + 10×1 ns |
| `--prod-ns` / `--prod-chunks` | 生产长度/分块 | 每 chunk 建议 ≤25 ns；walltime 默认 24 h |
| `--replicas` | 副本数 | 科学结论建议 ≥3；流程验证 1 即可 |
| `--pbs-nodes` | 节点名或 `1` | 本集群建议 `node02`（避开 node03） |
| `--ppn` / `--gpus` | MD 作业 CPU/GPU 请求 | 默认 8/1 |
| `--param-pbs` / `--no-param-pbs` | 配体参数化是否提交 PBS | 站点默认开启；无 PBS 时用 `--no-param-pbs` |
| `--param-ppn` | 参数化 PBS 作业申请的 CPU 数 | 默认 40 |
| `--param-walltime` | 参数化作业 walltime | 默认 `04:00:00` |
| `--param-amber-sh` | 参数化使用的 AmberTools `amber.sh` | 站点默认 conda ABT；`/opt` Amber22 的 sqm 是串行的 |
| `--param-nodes` | 参数化 PBS 作业节点 | 默认跟 `pbs_nodes`；站点设为 `1` 让调度器挑空闲节点 |
| `--min-engine` | 最小化引擎 `sander` / `pmemd.cuda` | 默认 `sander`；`pmemd.cuda` 为实验选项，本示例大体系实测 BOND 发散，不建议直接用 |
| `--mmpbsa-interval` | MM/GBSA 取帧间隔 | 10 ns 轨迹建议 20（后 50% 约 25 帧） |
| `--no-mmpbsa` | 关闭 MM/GBSA | 只做结构分析时用 |
| `--no-submit` / `--dry-run` | 不提交 | 先检查再提交 |
| `--force` | 重建指定体系 | 会清该体系产物，不动其他体系 |

注意：`--pbs-nodes`、`--mmpbsa-interval` 等参数在**建系时**写进渲染后的 PBS；事后改参数需要重建（或新建 campaign）。`resume` 不会改变已渲染脚本。

## 7. 结果查看与正确性校验

### 7.1 先看汇总

```bash
column -t -s $'\t' "$CAMP/report/summary.tsv"
```

`summary.tsv` 的关键列：

- `status`：`done` 表示三个作业都成功且 post 完成；
- `qc`：建系 QC 状态，应为 `ok`；
- `charge` / `charge_method`：配体电荷与来源；
- `stage` / `nstep` / `ns_day`：最后完成的 stage、步数、速度；
- `dg` / `dg_sd` / `frames`：MM/GBSA 结果与帧数；
- `warnings`：非致命告警，逐条看。

### 7.2 建系正确性

```bash
cat "$CAMP"/complexes/<label>/prepare/prepare_manifest.json
cat "$CAMP"/complexes/<label>/systems/<label>/build_manifest.json
cat "$CAMP"/complexes/<label>/systems/<label>/build_qc.tsv
cat "$CAMP"/complexes/<label>/systems/<label>/masks.env
```

应满足：

- `prepare_manifest.json`：`type=protein-ligand`，配体电荷与 `--charge` 一致，`reduce_used=true`，二硫键列表符合预期；
- `build_manifest.json`：`status=ok`，tleap `Errors=0`，`ATTN=0`，中和后净电荷 0；
- `build_qc.tsv`：电荷、水、离子、`RADIUS_SET`、mask 检查通过；
- `masks.env`：`LIGAND_MASK` 为配体残基号（如 `:355`），`RECEPTOR_MASK` 覆盖蛋白。

### 7.3 MD 完成度

```bash
ls "$CAMP"/complexes/<label>/md/<label>/r1/ | sort
grep -l 'Total wall time' "$CAMP"/complexes/<label>/md/<label>/r1/*.out
cat "$CAMP"/complexes/<label>/analysis/post_status_<label>_r1.json
```

每个 stage 的 `.out` 应有 `Total wall time`，且最后 `NSTEP` 达到 `nstlim`；`npt.rst7`、`prod_final.rst7` 应存在且非空；`post_status` 中 `cpptraj`、`analyze`、`mmpbsa` 应为 0。

### 7.4 轨迹与分析

```bash
ls -lh "$CAMP"/complexes/<label>/md/<label>/r1/*.nc
ls -lh "$CAMP"/complexes/<label>/analysis/cpptraj/
ls -lh "$CAMP"/complexes/<label>/analysis/python/
column -t -s $'\t' "$CAMP"/complexes/<label>/analysis/mmpbsa/<label>_r1/summary.tsv
```

MM/GBSA 的 `frames` 应与预期一致（10 ns、10 ps 间隔、后 50% 窗口、`interval=20` → 约 25 帧）。帧数过少时工具会警告，必要时减小 `--mmpbsa-interval` 并重建。

### 7.5 输入未被修改

```bash
sha256sum "$DATA"/*.pdb
```

应与运行前记录一致；所有产物都在 `$CAMP` 下，源目录只读使用。

### 7.6 八条核心原则

1. **配体电荷/质子化必须手填**；脚本不猜，模糊时报错。
2. `antechamber -at gaff2` ⇔ `parmchk2 -s gaff2` ⇔ `leaprc.gaff2`，三者必须一致。
3. `tleap` 必须 `Errors = 0`；`ATTN` 默认失败，确认化学合理后才允许放行。
4. `igb=5 ↔ mbondi2`，`igb=8 ↔ mbondi3`；改 `igb` 必须重建拓扑。
5. 生产 chunk 有状态，严格顺序；PBS 模板会自动跳过已完成 stage。
6. MM/GBSA 用**溶剂化拓扑** + `ante-MMPBSA.py`；不要用 dry prmtop。
7. 少于 5 个分子不要宣称泛化；3 个分子完美排序的随机概率是 1/6；MM/GBSA 不是 Kd。
8. 不要为了匹配实验排序去调介电常数、盐浓度、窗口或 QC 阈值。

## 8. 常见报错与解决办法

| 报错/现象 | 原因 | 解决 |
|---|---|---|
| `ligand charge/bond-order assignment is ambiguous` | 配体电荷/键级多解 | 用 `inspect` 看候选，显式 `--charge` 或 `--ligand-smiles` |
| `--charge N is not chemically valid` | 指定电荷与结构不符 | 检查质子化态；换 `--ligand-smiles` |
| `unknown HETATM residues` | 有未识别的 HETATM | 指定 `--ligand`；或按需 `--drop-unknown-hetatm` |
| `glycan_unsupported` / `metal_unsupported` | 含糖链/二价金属 | 当前版本不支持；先手工去除或改用其他流程 |
| `tleap finished with Errors = N` | 缺参数/残基名不合法/电荷异常 | 看 `tleap.out` 与 `build_qc.tsv` |
| `tleap emitted N ATTN warning(s)` | 参数近似或化学不合理 | 逐条审查；不要默认放行 |
| `Small box detected` | 盒子太小 | 增大 `--pad`（如 15）并重建 |
| 作业卡在 `01_min_restraint` | node03 已知问题 | 用 `--pbs-nodes node02` 重建/重投 |
| 作业 `exit 1` / walltime 被杀 | 资源不足/超时 | `resume`；检查 PBS 日志与 `#PBS -l walltime` |
| `not submitting any job because at least one complex failed` | fail-closed | 修复失败体系后重跑；或明确 `--submit-partial` |
| 同一 campaign 再 `run` 报 duplicate | 输入 SHA 已存在 | 用 `resume`；新建 campaign；或 `--force` |
| wrapper 报找不到 Python | 源码 checkout 无 `envs/` | 设 `AMBER_AUTO_PYTHON`，或使用共享安装入口 |
| MM/GBSA 帧数很少 | 默认 `interval=100` + 后 50% 窗口 | 用 `--mmpbsa-interval 20` 或更小，重建 post PBS |

更多错误码见 [app/docs/ERRORS.md](app/docs/ERRORS.md)，系统排错见 [app/docs/TROUBLESHOOTING.md](app/docs/TROUBLESHOOTING.md)。

## 9. FAQ

**Q：跑一半断了怎么办？**
A：`$AMBER_AUTO resume --campaign "$CAMP"`。已完成的 stage 会跳过；未完成的 stage 及其下游会重算。

**Q：能换节点吗？**
A：节点在**建系时**写进 PBS。已建系后换节点需要重建或手工改 PBS；建议一开始就加 `--pbs-nodes node02`。

**Q：为什么默认只跑 100 ns / 1 副本？**
A：这是流程默认值，不是科学建议。生产长度和副本数应按体系和问题设定。

**Q：两个体系会并行吗？**
A：每个体系/副本是一条独立依赖链，Torque 有资源就会并行；工具本身不设并发上限。

**Q：一个体系失败，成功的会提交吗？**
A：默认不会（fail-closed）；加 `--submit-partial` 才只提交成功的体系。

**Q：链 B 为什么被当作受体的一部分？**
A：自动分类把两条蛋白链都视为受体（蛋白–蛋白/多链模式）。如果生物学上链 B 是 partner，用 `--receptor-chains A --partner-chains B` 重建。

**Q：MM/GBSA 能当 Kd 用吗？**
A：不能。它不含熵，是相对指标；不要跨靶点排名，也不要用它替代实验亲和力。

**Q：可以只做建系不做 MD 吗？**
A：可以，`run ... --no-submit` 或 `build ...`；后者只建系不渲染提交计划。

**Q：能不能把最小化放到 GPU？**
A：可以试 `--min-engine pmemd.cuda`，但**默认不用**。本工具包在 86,964 atoms 示例上实测 `pmemd.cuda` 最小化会让 `01_min_restraint` 的 `BOND` 从 222 升到 69,496 kcal/mol，结构发散，后续 `02_min_free` 继承坏结构。当前默认仍是 `sander.MPI`，较慢但稳定；若要用 GPU，必须检查 `BOND`、`FINAL RESULTS` 和下一阶段初始能量，确认没有发散。

**Q：配体参数化为什么申请 40 CPU？`sqm` 真的会用满吗？**
A：取决于 AmberTools 环境。`/opt/software/amber22` 的独立 `sqm` 链接静态 BLAS，是串行的，只占 1 核；conda `ambertools_env` 的 `sqm` 链接 OpenBLAS pthreads，设置 `OPENBLAS_NUM_THREADS=40` 后 BLAS 部分可用满 40 核。本站点 `auto.site.json` 的 `param_amber_sh` 指向 conda ABT，所以默认能利用 40 核；PBS 日志会打印 `sqm` 路径与 BLAS 链接供核对。若改用 `/opt` Amber22，40 核只是分配配额，实际 sqm 仍 1 核。

## 10. 高级：amber-kit 旧流程

需要 RNA–配体、NMR 多模型、已有项目或精细控制时，使用 `bin/amber-kit`。完整教程见：

- [app/TUTORIAL.md](app/TUTORIAL.md)
- [app/docs/ADVANCED_AMBER_KIT.md](app/docs/ADVANCED_AMBER_KIT.md)
- [app/pbs/README.md](app/pbs/README.md)
- [app/config/README.md](app/config/README.md)

快速入口：

```bash
cp -r /opt/software/amber_general_workflow ~/my_project   # 或从源码复制
cd ~/my_project
bin/amber-kit init
$EDITOR config/config.sh config/systems.tsv
bin/amber-kit check
bin/amber-kit build
bin/amber-kit inputs
bin/amber-kit submit equil --label lig1 --rep 1 --dry-run
bin/amber-kit submit equil --label lig1 --rep 1
bin/amber-kit status --queue
```
