# PBS 前缀与作业脚本指南

本目录的 `*.pbs.in` 是**模板**：`{{KEY}}` 由 `amber-kit submit` 用 `scripts/render_template.py` 严格替换后生成项目内的 `.pbs` 文件，再从副本目录 `qsub`。

> 关键：`#PBS` 指令必须写在**可执行语句之前**、每行以 `#PBS` 开头；Torque 在 shell 运行前就解析它们。模板里的 `{{...}}` 不能留在最终 `.pbs` 里。

---

## 1. 本集群（cluster2）作业类型

| 作业 | 模板 | `nodes` | `ppn` | `gpus` | `walltime` | 用途 |
|---|---|---|---:|---:|---:|---|
| 平衡 | `run_equil.pbs.in` | 1 | 8 | 1 | 24 h | min1→min2→heat→04a..04j |
| 生产 | `run_prod.pbs.in` | 1 | 8 | 1 | 24 h | c01→c04（4×25 ns） |
| 单 chunk | `run_prod_chunk.pbs.in` | 1 | 8 | 1 | 24 h | 只跑一个 25 ns chunk |
| 双副本 | `run_two_gpu.pbs.in` | 1 | 48 | 无 | 24 h | 整节点 2×GPU |
| 多副本 | `run_multi_gpu.pbs.in` | 1 | 48 | 无 | 24 h | 整节点 N×GPU（N 可配） |
| 配体参数化 | `run_param_ligand.pbs.in` | 1 | 4 | 无 | 2 h | antechamber/parmchk2 |
| 分析 | `run_analysis.pbs.in` | 1 | 8 | 无 | 4 h | cpptraj + MDAnalysis |
| MM/GBSA | `run_mmpbsa.pbs.in` | 1 | 8 | 无 | 4 h | MMPBSA.py |
| 一键后处理 | `run_post.pbs.in` | 1 | 8 | 无 | 12 h | cpptraj + 分析 + MM/GBSA（`amber-auto` 使用） |

本集群：Torque 6.1.1.1，队列 `batch`，节点 `node02`–`node08`，每节点 48 核 + 4×RTX 3090。
`#PBS -l nodes=1:ppn=8:gpus=1` 已验证可用，Torque 会写 `$PBS_GPUFILE` 并设置 `CUDA_VISIBLE_DEVICES`。

`amber-auto` 用 Torque 依赖自动串联三个作业，无需手动等待：

```bash
eq=$(cd replica && qsub run_equil.pbs)
prod=$(cd replica && qsub -W depend=afterok:$eq run_prod.pbs)
post=$(cd replica && qsub -W depend=afterok:$prod run_post.pbs)
```

---

## 2. 一个标准 GPU 头，逐行解释

```bash
#!/bin/bash
#PBS -N myproject_lig1_r1_equil
#PBS -q batch
#PBS -l nodes=1:ppn=8:gpus=1
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o /home/user1/projects/myproject/logs/
#PBS -m n
```

| 指令 | 含义 | 建议 |
|---|---|---|
| `#PBS -N` | 作业名 | `<project>_<label>_r<rep>_<stage>`；`qstat` 第一列可见；不要用空格/中文 |
| `#PBS -q` | 队列 | 本集群只有 `batch` |
| `#PBS -l nodes=1:ppn=8:gpus=1` | 1 节点、每节点 8 个 CPU 进程、1 张 GPU | 单个 `pmemd.cuda` 只要 1 张 GPU；8 核给 min 的 `sander.MPI`/PME 用 |
| `#PBS -l walltime=24:00:00` | 墙钟上限 | 100 ns（18–27k 原子）约 7–8 h；4×25 ns 可放进 24 h |
| `#PBS -j oe` | stdout/stderr 合并 | 看一个 `.o<jobid>` 即可 |
| `#PBS -o /path/logs/` | 输出目录 | **必须已存在**；`amber-kit init` 创建 `logs/`；末尾斜杠 |
| `#PBS -m n` | 不邮件 | 需要邮件用 `-m abe` 并加 `#PBS -M you@example.com` |

### 不要这样写

- `#PBS -l nodes=1:ppn=8:gpus=2` 再跑一个 `pmemd.cuda`：第二张 GPU 闲置；双副本请用 `run_multi_gpu.pbs.in` / `run_two_gpu.pbs.in`，整节点 `ppn=<核数>`。
- 把 `#PBS` 写在 `cd`/`source` 之后：Torque 不解析，资源按默认值分配。
- `-o` 指向不存在的目录：作业直接失败，日志找不到。

---

## 3. GPU 绑定

Torque 在作业启动前写 `$PBS_GPUFILE`，每行一个分配到的 GPU，例如：

```
node02-gpu0
node02-gpu1
```

模板里的标准写法：

```bash
export CUDA_DEVICE_ORDER=PCI_BUS_ID
if [ -z "${CUDA_VISIBLE_DEVICES:-}" ] && [ -n "${PBS_GPUFILE:-}" ] && [ -f "$PBS_GPUFILE" ]; then
  export CUDA_VISIBLE_DEVICES=$(sed -E 's/.*[^0-9]([0-9]+)$/\1/' "$PBS_GPUFILE" | paste -sd, -)
fi
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-${GPU_ID:-0}}"
# 单个 pmemd.cuda 只用一张；若分配了多张，取第一张
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES%%,*}"
```

要点：

- 优先保留 Torque 已设的 `CUDA_VISIBLE_DEVICES`，不要无条件覆盖。
- `sed -E 's/.*[^0-9]([0-9]+)$/\1/'` 同时兼容 `node02-gpu0`、`/dev/nvidia0`、`GPU0` 等格式。
- `nvidia-smi` 只在计算节点有；登录节点没有。
- 双/多副本模板用**显式** `CUDA_VISIBLE_DEVICES=0/1/...`，并要求整节点 `ppn=<核数>`，避免和其他作业共享 GPU。

---

## 4. 工作目录与日志

```bash
cd "$PBS_O_WORKDIR" || exit 1
```

- `$PBS_O_WORKDIR` 是**提交目录**，不是脚本所在目录。
- 因此 `amber-kit submit` 会 `cd` 到 `md/<label>/r<rep>/` 再 `qsub`。
- 手动提交也要这样：

```bash
cd /path/to/md/lig1/r1
qsub run_equil.pbs
```

- 如果输出目录用 `#PBS -o`，Torque 写 `<jobname>.o<jobid>`；`#PBS -j oe` 时错误也在这个文件。
- **本集群实测**：Torque 作业结束后没有把 `#PBS -o` 的输出写回 NFS（`Output_Path` 指向的文件不存在）。
  因此模板会在 `LOG_DIR` 内自行 `exec >> {{JOB_NAME}}.o$PBS_JOBID` 落一份日志，优先看这个文件。
- 模板还会把 Amber 自己的 `.out/.log/.mdinfo` 写在副本目录；先看 `.out`，再看 PBS `.o`。

---

## 5. 墙钟、分段与续跑

`pmemd.cuda` 不会因为墙钟到点而优雅退出；Torque 先发 `SIGTERM`，随后 kill。因此：

- **生产分段**：默认 4×25 ns；每个 chunk 是一个有状态的 restart 链。
- **平衡分段**：min1→min2→heat→04a..04j；每段都写 `.rst7`。
- **续跑**：`run_equil.pbs` / `run_prod.pbs` 的 `stage_done()` 检查：
  - `.rst7` 存在且非空；
  - `.out` 含 `Total wall time`；
  - `.out` 最后的 `NSTEP >= nstlim`（min 没有 `nstlim`，只检查前两条）。
- 已完成的 stage **不重算、不覆盖**；某个 stage 未完成时，它和所有下游都重算，避免“旧 c03 基于旧 c02”。
- `flock` 防止同一副本并发；锁文件 `.run_equil.lock` / `.run_prod.lock`。
- 墙钟不够时：
  - 用 `run_prod_chunk.pbs.in` 一次提交一个 chunk；
  - 或减小 `PROD_CHUNK_NS`、增大 `PROD_CHUNKS`（例如 10×10 ns）。
- GPU PME box 错误（`Periodic box dimensions have changed too much`）：
  - 把 `RELEASE_NS` 改短、`RELEASE_STAGES` 增多（例如 100×0.1 ns），由逐 stage 续跑处理；
  - 或先用 CPU `sander.MPI` 跑密度平衡，再回 GPU。

---

## 6. 复制到其他集群（GPU 数量不同）

本工具包只把“集群相关”的值放在 `config/config.sh`。复制到其他集群时改这些：

| 变量 | 作用 | 其他集群怎么改 |
|---|---|---|
| `AMBER_SH` | Amber 初始化脚本 | 改成该集群的 `amber.sh` 路径 |
| `PY_QC` | 有 `parmed` 的 Python | 改成该集群的 conda/venv 解释器 |
| `PY_ANALYSIS` | 有 MDAnalysis/scipy/shapely 的 Python | 同上 |
| `PBS_QUEUE` | 队列名 | `batch` 或该集群队列 |
| `PBS_PPN` | 单 GPU 作业的 CPU 核数 | 按节点核数/排队策略改（常见 4/8/16） |
| `PBS_GPUS` | 单个 `pmemd.cuda` 作业申请的 GPU 数 | **保持 1**；`pmemd.cuda` 只用一张 |
| `PBS_GPUS_PER_NODE` | 整节点模板可用 GPU 数 | 改成节点实际 GPU 数（2/4/8） |
| `PBS_TWO_GPU_PPN` | 整节点模板的 ppn | 改成节点 CPU 核数（或留空让 Torque 分配） |
| `PBS_LOG_DIR` | 日志目录 | 改成该集群可写目录 |

### GPU 数量不同的两种情况

1. **每个作业仍然只用 1 张 GPU**（最常见）：
   - `PBS_GPUS=1` 不变；
   - 节点有 4/8 张 GPU 只是能同时跑更多作业，不影响单副本模板；
   - `#PBS -l nodes=1:ppn=...:gpus=1` 即可。

2. **一个作业里跑多个副本，每个副本一张 GPU**：
   - 用 `run_multi_gpu.pbs.in`；
   - 把 `PBS_GPUS_PER_NODE` 改成节点 GPU 数；
   - 把 `PBS_TWO_GPU_PPN` 改成节点 CPU 核数（整节点分配，避免共享）；
   - 提交：

```bash
amber-kit submit multi-gpu --pair lig1:1,lig1:2,lig2:1,lig2:2 --dry-run
# 确认 REP_DIRS/GPUS_PER_NODE 后去掉 --dry-run
```

   - 模板按 `REP_DIRS` 个数启动，最多 `GPUS_PER_NODE` 个；每个副本用 `CUDA_VISIBLE_DEVICES=0..N-1`。

### Torque → SLURM 对照（只改头）

| Torque | SLURM |
|---|---|
| `#PBS -N job` | `#SBATCH -J job` |
| `#PBS -q batch` | `#SBATCH -p batch` |
| `#PBS -l nodes=1:ppn=8:gpus=1` | `#SBATCH --nodes=1 --ntasks-per-node=8 --gres=gpu:1` |
| `#PBS -l walltime=24:00:00` | `#SBATCH --time=24:00:00` |
| `#PBS -j oe` | `#SBATCH --output=...`（默认合并） |
| `#PBS -o dir/` | `#SBATCH --output=dir/%x.o%j` |
| `$PBS_O_WORKDIR` | `$SLURM_SUBMIT_DIR` |
| `$PBS_GPUFILE` | `$CUDA_VISIBLE_DEVICES`（SLURM 通常直接设） |
| `qsub` / `qstat` / `qdel` | `sbatch` / `squeue` / `scancel` |

本工具包默认只提供 Torque 模板；SLURM 用户可复制 `pbs/*.pbs.in` 改头，`amber-kit submit` 的渲染逻辑不变（但 `qsub` 调用需改）。

---

## 7. 排错速查

| 症状 | 原因/处理 |
|---|---|
| `qstat` 里没有作业 | `qsub` 报错被忽略；看终端/PBS `.e`；检查 `#PBS` 是否在可执行语句前 |
| 作业跑在错误目录 | `PBS_O_WORKDIR` 是提交目录；从副本目录 `qsub` |
| 输出文件找不到 | `#PBS -o` 目录不存在或无写权限；`qstat -f <id>` 看 `Output_Path` |
| `CUDA_VISIBLE_DEVICES` 为空 | Torque 未设；检查 `PBS_GPUFILE`；模板回退 `GPU_ID` |
| `pmemd.cuda` 找不到 GPU | 登录节点没有 GPU；必须在计算节点；看 `nvidia-smi -L` |
| 作业被墙钟杀掉 | 用 `run_prod_chunk.pbs` 或减小 chunk；看 `.rst7` 时间戳判断从哪续跑 |
| 重投后从头发车 | 检查 `stage_done` 的 `.rst7/.out/Total wall time/NSTEP`；不要手动删 `.rst7` |
| 双副本互相覆盖 | 两个 `REP_DIR` 相同；模板会拒绝 |
| 整节点模板只跑一个 GPU | `GPUS_PER_NODE` 设成 2/4/8；`REP_DIRS` 个数不超过它 |
| `qstat` 显示 `C` 很久 | 已完成作业保留；`qdel -p <id>` 清理 |

---

## 8. 手工提交示例

```bash
source /opt/software/amber22/amber.sh

# 平衡
cd /path/to/project/md/lig1/r1
qsub run_equil.pbs

# 生产（npt.rst7 存在后）
qsub run_prod.pbs

# 单 chunk
qsub run_05_prod_c02.pbs

# 分析 / MM-GBSA
qsub analysis_lig1_r1.pbs
qsub mmpbsa_lig1_r1.pbs

qstat -a
qstat -f <jobid>
```

更推荐用统一入口（会先渲染、再 `qsub`）：

```bash
amber-kit submit equil --label lig1 --rep 1 --dry-run
amber-kit submit equil --label lig1 --rep 1
amber-kit submit prod  --label lig1 --rep 1
amber-kit submit multi-gpu --pair lig1:1,lig1:2 --dry-run
```
