# scripts/prodigal.py
from glob import glob
import os, sys
import logging, traceback
from itertools import repeat
import subprocess
from multiprocessing import Pool

logging.captureWarnings(True)

# ----------------- exception -> log -----------------
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.error("".join(["Uncaught exception: ",
                           *traceback.format_exception(exc_type, exc_value, exc_traceback)]))
sys.excepthook = handle_exception

# ----------------- progress (optional) --------------
try:
    from tqdm import tqdm
except ImportError:
    tqdm = list
    logging.warning("optional package tqdm is not installed")

# ----------------- helpers --------------------------
def run_shell(cmd: str):
    """Run a shell command; raise on failure with captured stderr."""
    cp = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if cp.returncode != 0:
        raise RuntimeError(f"CMD failed ({cp.returncode}): {cmd}\nSTDERR:\n{cp.stderr}")

def strip_suffix(s: str, suffix: str) -> str:
    """Remove exact suffix if present (safe alternative to rstrip)."""
    return s[:-len(suffix)] if suffix and s.endswith(suffix) else s

def splitext_ignore_gz(path: str):
    base, ext = os.path.splitext(path)
    if ext == ".gz":
        base, ext = os.path.splitext(base)
    return base, ext  # ext: .fna/.fa/.fasta/.gff/...

def norm_from_path(p: str, fmt1: str, fmt2: str) -> str:
    """统一样本名：去掉 .gz、再去掉 config 中的 format1/format2 或常见扩展."""
    b = os.path.basename(p.strip())
    if b.endswith(".gz"):
        b = b[:-3]
    # 优先按配置的两种格式剥离
    b = strip_suffix(b, fmt1) if fmt1 else b
    b = strip_suffix(b, fmt2) if fmt2 else b
    # 再兜底去一层常见扩展
    b, _ = os.path.splitext(b)
    return b

# ----------------- build dictionary -----------------
def prepare_faa_paths(paths, fmt1: str, fmt2: str):
    """
    返回 {sample_name -> subdir_relative_path}
    subdir 逻辑：沿用原脚本的 “最后3个目录（不含文件名）拼接”
    若路径层级不足则尽量取能取到的部分。
    """
    d = {}
    for raw in paths:
        path = raw.strip()
        if not path:
            continue
        sample = norm_from_path(path, fmt1, fmt2)
        parts = path.split("/")
        # 取 -4:-1（不足时自动兼容）
        sub = "/".join(parts[-4:-1]) if len(parts) >= 4 else "/".join(parts[:-1])
        d[sample] = sub
    return d

# ----------------- any2fasta ------------------------
def any2fasta(path_list_file: str, out_fasta_dir: str, fmt1: str, fmt2: str):
    """把输入清单统一转成 .fasta 到 out_fasta_dir；返回样本名列表。"""
    with open(path_list_file, "r") as fh:
        paths = [ln.strip() for ln in fh if ln.strip()]

    os.makedirs(out_fasta_dir, exist_ok=True)
    sample_names = []

    for file in tqdm(paths):
        filename, ext = splitext_ignore_gz(file)
        sample = norm_from_path(file, fmt1, fmt2)
        sample_names.append(sample)

        dst = f"{out_fasta_dir}/{sample}.fasta"
        # 压缩/gff/或非常见扩展 -> any2fasta
        if file.endswith(".gz") or ext.lower() == ".gff":
            logging.info(f"{sample}: not plain FASTA -> any2fasta")
            run_shell(f"any2fasta {file} > {dst}")
        else:
            logging.info(f"{sample}: FASTA detected -> symlink")
            src = os.path.abspath(file)
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(src, dst)

    return sample_names

# ----------------- prodigal runners -----------------
def run_prodigal_one(input_fasta: str, parameters: str, out_base: str, dct: dict):
    """
    运行单个样本的 prodigal。
    若样本在字典里：写入 out_base/{faa,fna}/{subdir}/{sample}.*
    否则：回退到 out_base/{faa,fna}/{sample}.*
    """
    sample = os.path.splitext(os.path.basename(input_fasta))[0]
    sub = dct.get(sample, None)

    if sub:
        faa_dir = f"{out_base}/faa/{sub}"
        fna_dir = f"{out_base}/fna/{sub}"
    else:
        logging.warning(f"[dictionary-miss] {sample} not in dictionary, flatten to top-level.")
        faa_dir = f"{out_base}/faa"
        fna_dir = f"{out_base}/fna"

    os.makedirs(faa_dir, exist_ok=True)
    os.makedirs(fna_dir, exist_ok=True)

    faa = f"{faa_dir}/{sample}.faa"
    fna = f"{fna_dir}/{sample}.fna"

    logging.info(f"Running prodigal on {sample}")
    cmd = f"prodigal {parameters} -i {input_fasta} -a {faa} -d {fna}"
    run_shell(cmd)  # 失败会抛异常并写入 snakemake 的 log

def run_multiple_prodigal(in_dir: str, threads: int, parameters: str, out_base: str, dct: dict, ext: str = ".fasta"):
    fastas = glob(os.path.join(in_dir, "*" + ext))
    assert fastas, f"No input *{ext} found in {in_dir}"
    with Pool(int(threads)) as pool:
        pool.starmap(run_prodigal_one, zip(fastas, repeat(parameters), repeat(out_base), repeat(dct)))

# ----------------- main -----------------------------
if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    cfg = snakemake.config
    fmt1 = cfg.get("format1", "")
    fmt2 = cfg.get("format2", "")
    tmp_dir = cfg.get("temporary_dir", "/tmp")
    prodigal_output = cfg.get("directory_faa", "quality_filtering/faa_files")

    # 保证顶层目录存在
    for sub in ("faa", "fna"):
        os.makedirs(f"{prodigal_output}/{sub}", exist_ok=True)

    # 每批次的中间 fasta 目录
    intermediate_fasta_dir = f"{tmp_dir}/intermediate_results/fasta/{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}"
    os.makedirs(intermediate_fasta_dir, exist_ok=True)

    # ---- 构建字典（保留映射逻辑，但更健壮）----
    dictionary = {}
    if cfg.get("input_text_file", ""):
        with open(cfg["input_text_file"], "r") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        dictionary = prepare_faa_paths(lines, fmt1, fmt2)
    elif cfg.get("genomes", ""):
        paths = glob(f"{cfg['genomes']}/*{fmt1}") if fmt1 else []
        paths += glob(f"{cfg['genomes']}/*{fmt2}") if fmt2 else []
        paths = [p.strip() for p in paths]
        dictionary = prepare_faa_paths(paths, fmt1, fmt2)

    # ---- 1) 统一 any2fasta 到中间目录 ----
    sample_names = any2fasta(snakemake.input[0], intermediate_fasta_dir, fmt1, fmt2)

    # ---- 2) 多进程跑 prodigal ----
    run_multiple_prodigal(
        intermediate_fasta_dir,
        threads=snakemake.threads,
        parameters=snakemake.params.parameters,
        out_base=prodigal_output,
        dct=dictionary,
        ext=".fasta",
    )

    # ---- 3) 输出目录放 .faa 的软链（相对路径）----
    os.makedirs(snakemake.output[0], exist_ok=True)
    faa_src_base = os.path.relpath(f"{prodigal_output}/faa", snakemake.output[0])

    for s in sample_names:
        sub = dictionary.get(s, "")
        src = f"{faa_src_base}/{sub}/{s}.faa" if sub else f"{faa_src_base}/{s}.faa"
        dst = f"{snakemake.output[0]}/{s}.faa"
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(src, dst)

    logging.info("Prodigal finished successfully.")
