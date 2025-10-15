#!/bin/bash
#SBATCH --job-name=quality-filtering
#SBATCH --partition=cu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=120:00:00
#SBATCH --mem=150G
#SBATCH --output=/data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow/log/humsub_%A.out
#SBATCH --error=/data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow/log/humsub_%A.err
#SBATCH --mail-user=xg362573@163.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
date

# =============================
# 1) 环境准备
# =============================
mkdir -p /data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow/log
export PATH="/home/szhang/anaconda3/bin:$PATH"
source /home/szhang/anaconda3/etc/profile.d/conda.sh

# 设置临时目录与 Conda 缓存
export TMPDIR="${TMPDIR:-/tmp}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-/data1/projects/liullab/.conda_pkgs}"
mkdir -p "$CONDA_PKGS_DIRS"

# 禁用 Gurobi 调度警告
# export PULP_CBC_CMD=/usr/bin/false
# 禁用 libmamba 求解器，使用 classic
# export CONDA_SOLVER=classic #### 这个任务用这个会非常慢16个小时都下载不下来

# =============================
# 2) 进入工作目录
# =============================
cd /data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow

# =============================
# 3) 解锁（防止上次崩溃锁文件）
# =============================
/home/szhang/anaconda3/envs/snakemake/bin/snakemake \
  -s ./Snakefile \
  --configfile ../config/default_config.yaml \
  --use-conda --unlock || true

# =============================
# 4) 继续执行未完成的任务（不重跑已完成）
# =============================
/home/szhang/anaconda3/envs/snakemake/bin/snakemake \
  -s ./Snakefile \
  --configfile ../config/default_config.yaml \
  --use-conda --conda-frontend conda \
  --cores "${SLURM_CPUS_PER_TASK:-4}" \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 60 \
  --scheduler greedy \
  --printshellcmds \
  --show-failed-logs

# /home/szhang/anaconda3/envs/snakemake/bin/snakemake \
#   -s ./Snakefile \
#   --configfile ../config/default_config.yaml \
#   --use-conda --conda-frontend conda \
#   --cores "${SLURM_CPUS_PER_TASK:-4}" \
#   --rerun-incomplete --keep-going --latency-wait 60 \
#   --scheduler greedy --printshellcmds --show-failed-logs \
#   --allowed-rules prodigal

date
