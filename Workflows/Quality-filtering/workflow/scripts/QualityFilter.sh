#!/bin/bash
#SBATCH --job-name=quality-filtering
#SBATCH --partition=cu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4                 # 与 --cores 对齐
#SBATCH --time=120:00:00
#SBATCH --mem=100G                        # 留意分区上限；不满足会被拒
#SBATCH --output=/data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow/log/humsub_%A.out
#SBATCH --error=/data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow/log/humsub_%A.err
#SBATCH --mail-user=xg362573@163.com
#SBATCH --mail-type=BEGIN,END,FAIL        # 常用即可，ALL 也行

set -euo pipefail
date

# 1) 准备目录与环境变量（加速/避免权限问题）
mkdir -p /data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow/log
export PATH="/home/hjlu/miniconda3/bin:$PATH"
source /home/hjlu/miniconda3/etc/profile.d/conda.sh
export TMPDIR="${TMPDIR:-/tmp}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-/data1/projects/liullab/.conda_pkgs}"
# 如遇 libmamba/ABI 问题，启用 classic solver：
export CONDA_SOLVER=classic

# 2) 进入工作目录
cd /data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow

# 3) 先解锁（安全）
/home/hjlu/miniconda3/envs/humsub/bin/snakemake \
  --snakefile /data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow/Snakefile \
  --configfile /data1/projects/liullab/humsub/Workflows/Quality-filtering/config/default_config.yaml \
  --use-conda --unlock || true

# 4) 正式跑：只依赖 rule all（不要手动列目标！）
/home/hjlu/miniconda3/envs/humsub/bin/snakemake \
  -s /data1/projects/liullab/humsub/Workflows/Quality-filtering/workflow/Snakefile \
  --configfile /data1/projects/liullab/humsub/Workflows/Quality-filtering/config/default_config.yaml \
  --use-conda --conda-frontend conda \
  --cores "${SLURM_CPUS_PER_TASK:-4}" \
  --rerun-incomplete \
  --latency-wait 60 \
  --keep-going \
  --printshellcmds \
  --show-failed-logs

date
