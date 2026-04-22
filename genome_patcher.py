"""
Genome patcher CLI.

Pipeline:
1. Align assembly contigs to patch sequences with minimap2.
2. Parse alignments and patch assembly sequences using aligned patch segments.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from gp_align import run_minimap2_alignment
from gp_io import build_output_paths, resolve_tmp_dir
from gp_patching import run_patch_stage
from gp_visualization import generate_patch_plots


def _add_io_label_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument("-l", "--label", required=True, help="Output label prefix")


def _add_verbose_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose progress logs",
    )


def _add_alignment_stage_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="minimap2 threads (default: 16)",
    )
    parser.add_argument(
        "--tmp-dir",
        default=None,
        help="Temporary directory for intermediate files (default: <out-dir>/tmp)",
    )
    parser.add_argument(
        "--keep-unmapped",
        dest="drop_unmapped",
        action="store_false",
        default=True,
        help="Keep unmapped query sequences in output BAM (default is to drop them)",
    )


def _add_patching_filter_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--min-pid", type=float, default=0.999, help="Minimum identity filter")
    parser.add_argument("--min-aln-len", type=int, default=1000, help="Minimum aligned length")
    parser.add_argument(
        "--border-tolerance",
        type=int,
        default=10,
        help="Bases from patch borders considered open for chain linking (default: 10)",
    )


def _add_patching_mode_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-p",
        "--preserve-mode",
        action="store_true",
        dest="preserve_patch_seq",
        help=(
            "Keep patch-sequence intervals in chain bodies (do not substitute with assembly "
            "segments), while still linking chains through open endpoints and removing "
            "assembly contigs implied by selected hits"
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Patch genome assemblies using patch sequences")
    subparsers = parser.add_subparsers(dest="command", required=True)

    align_parser = subparsers.add_parser("align", help="Run stage 1 (minimap2 alignment)")
    align_parser.add_argument("--assembly", required=True, help="Assembly FASTA(.gz), used as query")
    align_parser.add_argument("--patch", required=True, help="Patch FASTA(.gz), used as reference")
    _add_io_label_args(align_parser)
    _add_alignment_stage_args(align_parser)
    _add_verbose_arg(align_parser)

    patch_parser = subparsers.add_parser("patch", help="Run stage 2 (patch from BAM)")
    patch_parser.add_argument("--assembly", required=True, help="Assembly FASTA(.gz)")
    patch_parser.add_argument("--patch", required=True, help="Patch FASTA(.gz)")
    patch_parser.add_argument(
        "--bam",
        default=None,
        help="Input BAM; default is <out-dir>/<label>.sorted.bam",
    )
    _add_io_label_args(patch_parser)
    _add_patching_filter_args(patch_parser)
    _add_patching_mode_args(patch_parser)
    _add_verbose_arg(patch_parser)

    run_parser = subparsers.add_parser("run", help="Run stage 1 + stage 2 end-to-end")
    run_parser.add_argument("--assembly", required=True, help="Assembly FASTA(.gz), query in minimap2")
    run_parser.add_argument("--patch", required=True, help="Patch FASTA(.gz), reference in minimap2")
    _add_io_label_args(run_parser)
    _add_alignment_stage_args(run_parser)
    _add_patching_filter_args(run_parser)
    _add_patching_mode_args(run_parser)
    _add_verbose_arg(run_parser)

    return parser


def _resolve_patch_bam_path(user_bam: Optional[str], default_bam: str) -> str:
    bam_path = user_bam if user_bam else default_bam
    if not Path(bam_path).exists():
        raise FileNotFoundError(
            f"BAM not found: {bam_path}. Run 'align' first or pass --bam explicitly."
        )
    return bam_path


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "align":
        paths = build_output_paths(args.out_dir, args.label)
        tmp_dir = resolve_tmp_dir(args.tmp_dir, args.out_dir)
        if args.verbose:
            print("[cli] command=align")
        run_minimap2_alignment(
            assembly_fasta=args.assembly,
            patch_fasta=args.patch,
            out_bam=paths.bam_path,
            threads=args.threads,
            drop_unmapped=args.drop_unmapped,
            tmp_dir=tmp_dir,
            verbose=args.verbose,
        )
        return

    if args.command == "patch":
        paths = build_output_paths(args.out_dir, args.label)
        bam_path = _resolve_patch_bam_path(args.bam, paths.bam_path)
        if args.verbose:
            print("[cli] command=patch")
        patch_result = run_patch_stage(
            assembly_fasta=args.assembly,
            patch_fasta=args.patch,
            bam_path=bam_path,
            out_fasta=paths.patched_fasta_path,
            out_changed=paths.changed_path,
            out_alignments_tsv=paths.alignments_tsv_path,
            min_pid=args.min_pid,
            min_aln_len=args.min_aln_len,
            border_tolerance=args.border_tolerance,
            preserve_patch_seq=args.preserve_patch_seq,
            verbose=args.verbose,
        )
        if args.verbose:
            print("[cli] stage=visualize")
        generate_patch_plots(
            context=patch_result.context,
            out_dir=args.out_dir,
            label=args.label,
            verbose=args.verbose,
        )
        return

    if args.command == "run":
        paths = build_output_paths(args.out_dir, args.label)
        tmp_dir = resolve_tmp_dir(args.tmp_dir, args.out_dir)
        if args.verbose:
            print("[cli] command=run")
            print("[cli] stage=align")
        run_minimap2_alignment(
            assembly_fasta=args.assembly,
            patch_fasta=args.patch,
            out_bam=paths.bam_path,
            threads=args.threads,
            drop_unmapped=args.drop_unmapped,
            tmp_dir=tmp_dir,
            verbose=args.verbose,
        )
        if args.verbose:
            print("[cli] stage=patch")
        patch_result = run_patch_stage(
            assembly_fasta=args.assembly,
            patch_fasta=args.patch,
            bam_path=paths.bam_path,
            out_fasta=paths.patched_fasta_path,
            out_changed=paths.changed_path,
            out_alignments_tsv=paths.alignments_tsv_path,
            min_pid=args.min_pid,
            min_aln_len=args.min_aln_len,
            border_tolerance=args.border_tolerance,
            preserve_patch_seq=args.preserve_patch_seq,
            verbose=args.verbose,
        )
        if args.verbose:
            print("[cli] stage=visualize")
        generate_patch_plots(
            context=patch_result.context,
            out_dir=args.out_dir,
            label=args.label,
            verbose=args.verbose,
        )
        return

    raise RuntimeError(f"unsupported command: {args.command}")


if __name__ == "__main__":
    main()
