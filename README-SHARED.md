# amber_general_workflow 共享安装

本目录是集群共享安装，所有用户可读/可执行；代码和数据目录只读，**不要直接修改**。
每个用户在自己的项目目录运行，`amber-campaign/`、轨迹、PBS 日志都写在自己的目录里。

```
/opt/software/amber_general_workflow/
├── bin/amber-auto        # 统一入口（推荐）
├── bin/amber-kit         # 旧流程入口（高级用户）
├── app/                  # 只读代码/模板/文档/测试
├── config/auto.site.json # 站点配置（Amber/Python/OpenBabel 路径）
├── envs/qc/              # parmed 环境（建系/QC）
├── envs/analysis/        # MDAnalysis/scipy/RDKit/netCDF4 环境（分析/电荷推断）
├── README-SHARED.md
└── VERSION
```

不需要 `source amber.sh`，也不需要 `conda activate`；入口脚本会自己处理。
Amber 使用 `/home/ylqx/amber22`，OpenBabel 使用
`/opt/schrodinger2021-2/utilities/obabel`。

## 普通用户：三步

```bash
mkdir -p ~/my_project && cd ~/my_project
cp /path/to/your_complex.pdb .

/opt/software/amber_general_workflow/bin/amber-auto inspect your_complex.pdb
/opt/software/amber_general_workflow/bin/amber-auto run your_complex.pdb
```

`run` 会建系、QC、生成 MD 输入，并自动提交 `equil → prod → post` 的 Torque
依赖链。输出在 `./amber-campaign/`。

可选：把入口加入 `PATH`，之后可直接写 `amber-auto`：

```bash
export PATH=/opt/software/amber_general_workflow/bin:$PATH
```

```bash
/opt/software/amber_general_workflow/bin/amber-auto status --queue
/opt/software/amber_general_workflow/bin/amber-auto report
/opt/software/amber_general_workflow/bin/amber-auto resume          # 断点续跑
```

如果作业落在 node03 并卡在 `01_min_restraint`（见下方“本集群节点注意”），直接固定节点：

```bash
/opt/software/amber_general_workflow/bin/amber-auto run complex.pdb --pbs-nodes node02
```

## 配置覆盖（每个用户可不同）

站点配置是全集群共用的。需要改力场、水模型、PBS 资源、MD 长度、节点等时，复制一份到
自己的目录再改：

```bash
cp /opt/software/amber_general_workflow/config/auto.site.json ~/my-auto.json
$EDITOR ~/my-auto.json
export AMBER_AUTO_CONFIG=~/my-auto.json     # 之后所有调用都用这份
```

也可以临时用命令行覆盖（优先级最高）：

```bash
.../bin/amber-auto run complex.pdb --prod-ns 200 --replicas 3 --pad 12
.../bin/amber-auto run complex.pdb --pbs-nodes node02
```

## PBS 说明

- 队列 `batch`，单 GPU 作业头为 `nodes=1:ppn=8:gpus=1`；每个用户的作业以自己的
  PBS 账号运行，互不影响。
- 可用 `--pbs-nodes node02`（或配置 `"pbs_nodes": "node02"`）把作业固定到指定节点；
  默认 `1` 让 Torque 调度。
- 生产默认 4×25 ns；分析/MM-GBSA 在 `post` 作业里自动运行。
- 本集群 Torque 不会把 `#PBS -o` 的输出写回 NFS，模板会自己在
  `amber-campaign/complexes/<label>/logs/<jobname>.o<jobid>` 落一份日志。
- 作业中断后直接 `resume`；已完成的 stage 会跳过。
- 模板会在 `source amber.sh` 前清理继承的 OpenMPI 库路径，避免多版本 `libmpi`
  混用导致 `mpirun` 挂起。

## 本集群节点注意（2026-10-31）

- **node03 目前会让 `sander.MPI` 卡在 `01_min_restraint`**：8 个 rank 100% CPU、
  不写 `.out`；用 `/opt/software/amber22` 和 `/home/ylqx/amber22` 两套 Amber 都能复现。
  基础 MPI（`mpirun hostname`、`MPI_Allreduce`）在 node03 正常，属于节点问题。
- 在管理员修复 / drain node03 之前，建议加 `--pbs-nodes node02`（或 node04–node08）。
- 该问题与共享安装、脚本版本无关；同一作业在 node02 上正常完成。

## 已验证

- 共享入口在 node02 完整跑通 `equil → prod → post`（作业 733/734/735，全部
  `exit 0`），生产 365.76 ns/day，MM/GBSA ΔG = 0.4541 ± 0.2806 kcal/mol。
- 共享目录的 `tests/smoke.sh` 和 `tests/auto_smoke.sh` 全部通过（16 项）。
- 以 `nobody` 用户验证了读取/执行共享代码与两个环境，并能成功 `inspect`。
- 验证 campaign（可删除）：`/home/user1/agw_shared_validate.CuE9Jq`

## 已知限制

- 仅 RNA 堆积分析需要 `shapely`，共享环境未安装（集群无法访问 PyPI）；蛋白体系不受影响。
- `pmemd.cuda` 对超小盒子（neighbor-list cell ≤2，cut=10 Å 时约 box < 35–40 Å）
  会报 `Small box detected`。测试小体系请用 `--pad 15` 或更大。
- `run_two_gpu` / `run_multi_gpu` 模板未在本集群实测；单 GPU 依赖链已验证。
- 本工具包不替代科学判断：配体电荷/质子化、力场/水模型、MM/GBSA 参数需自行确认。

## 管理员：更新与权限

共享目录应由 root 拥有并只读：

```bash
sudo chown -R root:root /opt/software/amber_general_workflow
sudo chmod -R a+rX /opt/software/amber_general_workflow
sudo chmod -R go-w /opt/software/amber_general_workflow
```

更新代码时，把新版本 rsync 到 `app/`（排除 `__pycache__`、campaign），并把代码/配置里的
`/home/user1/Software/anaconda3/envs/...` 替换为
`/opt/software/amber_general_workflow/envs/{qc,analysis}`。
默认 Amber 是 `/home/ylqx/amber22/amber.sh`（本集群 `/opt/software/amber22` 的
sander 串行/MPI 有问题，故默认改用前者）；如果目标共享安装的 Amber 位于
`/opt/software/amber22`，再把 `amber_sh`/`param_amber_sh` 相应替换。再把 `app/bin/amber-auto` 的
shebang 改为 `/opt/software/amber_general_workflow/envs/analysis/bin/python`，
并用 patched `config/auto.example.json` 生成 `config/auto.site.json`。
本安装的完整部署步骤可用 `history`/本 README 复现。

快速自检：

```bash
/opt/software/amber_general_workflow/envs/qc/bin/python -I -c "import parmed, numpy"
/opt/software/amber_general_workflow/envs/analysis/bin/python -I -c "import MDAnalysis, scipy, rdkit, netCDF4"
/opt/software/amber_general_workflow/bin/amber-auto --help
sudo -u nobody /opt/software/amber_general_workflow/bin/amber-auto --help
```
