from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import pysam

DEFAULT_MINIMAP2_THREADS = 16
MINIMAP2_PRESET = "asm20"


def run_minimap2_alignment(
    assembly_fasta: str,
    patch_fasta: str,
    out_bam: str,
    threads: int = DEFAULT_MINIMAP2_THREADS,
    drop_unmapped: bool = True,
    tmp_dir: Optional[str] = None,
    verbose: bool = False,
) -> str:
    """
    Align assembly (query) against patch sequences (reference), write sorted and indexed BAM.
    """
    minimap2_path = shutil.which("minimap2")
    if minimap2_path is None:
        raise RuntimeError("minimap2 not found in PATH")
    if verbose:
        print(f"[align] minimap2: {minimap2_path}")
        print(f"[align] assembly(query): {assembly_fasta}")
        print(f"[align] patch(reference): {patch_fasta}")
        print(f"[align] out_bam: {out_bam}")
        print(f"[align] tmp_dir: {tmp_dir or '<system-temp>'}")

    out_path = Path(out_bam)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".sam",
        delete=False,
        dir=tmp_dir,
        encoding="utf-8",
    ) as tmp_handle:
        sam_path = tmp_handle.name

    cmd = [
        minimap2_path,
        "-x",
        MINIMAP2_PRESET,
        "-t",
        str(threads),
        "-a",
        "--eqx",
    ]
    cmd.extend([patch_fasta, assembly_fasta])
    if verbose:
        print(f"[align] running: {' '.join(cmd)}")

    filtered_bam_path: Optional[str] = None
    try:
        with open(sam_path, "w", encoding="utf-8") as sam_out:
            if verbose:
                subprocess.run(cmd, stdout=sam_out, check=True)
            else:
                subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, check=True, text=True)

        sort_input = sam_path
        if drop_unmapped:
            if verbose:
                print("[align] dropping unmapped reads (flag 0x4)")
            with tempfile.NamedTemporaryFile(
                mode="wb",
                suffix=".bam",
                delete=False,
                dir=tmp_dir,
            ) as filtered_bam:
                filtered_bam_path = filtered_bam.name
            pysam.view("-b", "-F", "4", "-o", filtered_bam_path, sam_path, catch_stdout=False)
            sort_input = filtered_bam_path

        sort_threads = max(1, int(threads))
        if verbose:
            print(f"[align] sorting BAM with threads={sort_threads}")
        pysam.sort("-@", str(sort_threads), "-o", str(out_path), sort_input)
        pysam.index(str(out_path))
        if verbose:
            print(f"[align] wrote sorted BAM: {out_path}")
            print(f"[align] wrote BAM index: {out_path}.bai")
    except FileNotFoundError as exc:
        raise RuntimeError(f"required executable not found: {exc.filename}") from exc
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        if stderr:
            raise RuntimeError(f"minimap2 failed (exit {exc.returncode}): {stderr}") from exc
        raise RuntimeError(f"minimap2 failed (exit {exc.returncode})") from exc
    finally:
        if os.path.exists(sam_path):
            os.remove(sam_path)
        if filtered_bam_path and os.path.exists(filtered_bam_path):
            os.remove(filtered_bam_path)

    return str(out_path)
