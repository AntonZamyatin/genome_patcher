from __future__ import annotations

import gzip
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass(frozen=True)
class OutputPaths:
    bam_path: str
    patched_fasta_path: str
    changed_path: str
    alignments_tsv_path: str


@contextmanager
def open_text_auto(path: str) -> Iterator[TextIO]:
    if path.endswith(".gz"):
        handle = gzip.open(path, "rt")
    else:
        handle = open(path, "rt", encoding="utf-8")
    try:
        yield handle
    finally:
        handle.close()


def read_fasta_records(path: str) -> Dict[str, SeqRecord]:
    with open_text_auto(path) as handle:
        return {record.id: record for record in SeqIO.parse(handle, "fasta")}


def write_fasta(records: Dict[str, SeqRecord], out_fasta: str) -> None:
    with open(out_fasta, "w", encoding="utf-8") as handle:
        SeqIO.write(records.values(), handle, "fasta")


def write_changed_list(changed: List[str], out_path: str) -> None:
    with open(out_path, "w", encoding="utf-8") as handle:
        for contig in changed:
            handle.write(f"{contig}\n")


def build_output_paths(out_dir: str, label: str) -> OutputPaths:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    return OutputPaths(
        bam_path=str(out_path / f"{label}.sorted.bam"),
        patched_fasta_path=str(out_path / f"{label}.patched.fasta"),
        changed_path=str(out_path / f"{label}.changed.txt"),
        alignments_tsv_path=str(out_path / f"{label}.alignments.tsv"),
    )


def resolve_tmp_dir(tmp_dir: str | None, out_dir: str) -> str:
    if tmp_dir is None:
        tmp_path = Path(out_dir) / "tmp"
    else:
        tmp_path = Path(tmp_dir)
    tmp_path.mkdir(parents=True, exist_ok=True)
    return str(tmp_path)
