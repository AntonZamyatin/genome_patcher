from __future__ import annotations

import csv
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gp_io import read_fasta_records, write_changed_list, write_fasta


@dataclass(frozen=True)
class PatchHit:
    assembly_contig: str
    assembly_start: int
    assembly_end: int
    patch_contig: str
    patch_start: int
    patch_end: int
    patch_rc: bool
    pid: float
    aln_len: int


@dataclass(frozen=True)
class Segment:
    source: str  # "patch" | "assembly"
    contig: str
    start: int
    end: int
    reverse: bool = False

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)

    def __add__(self, other: Segment) -> Segment:
        if not isinstance(other, Segment):
            return NotImplemented
        if self.source != other.source or self.contig != other.contig or self.reverse != other.reverse:
            raise TypeError("segments are not mergeable")

        # Neighboring coordinates in either order represent zero-gap adjacency.
        if self.end == other.start or other.end == self.start:
            return Segment(
                source=self.source,
                contig=self.contig,
                start=min(self.start, other.start),
                end=max(self.end, other.end),
                reverse=self.reverse,
            )
        raise TypeError("segments are not mergeable")


@dataclass
class Chain:
    chain_id: str
    patch_ids: List[str]
    segments: List[Segment] = field(default_factory=list)
    chosen_hits: List[PatchHit] = field(default_factory=list)
    # (assembly_contig, coord, direction)
    left_open: Optional[Tuple[str, int, int]] = None
    right_open: Optional[Tuple[str, int, int]] = None

    @staticmethod
    def _append_merged(segments: List[Segment], segment: Segment) -> List[Segment]:
        if segment.length <= 0:
            return list(segments)
        if not segments:
            return [segment]
        try:
            merged = segments[-1] + segment
        except TypeError:
            return [*segments, segment]
        return [*segments[:-1], merged]

    def __add__(self, other: Segment | Chain) -> Chain:
        if isinstance(other, Segment):
            return Chain(
                chain_id=self.chain_id,
                patch_ids=list(self.patch_ids),
                segments=self._append_merged(self.segments, other),
                chosen_hits=list(self.chosen_hits),
                left_open=self.left_open,
                right_open=self.right_open,
            )
        if isinstance(other, Chain):
            out = Chain(
                chain_id=f"{self.chain_id}+{other.chain_id}",
                patch_ids=[*self.patch_ids, *other.patch_ids],
                segments=list(self.segments),
                chosen_hits=[*self.chosen_hits, *other.chosen_hits],
                left_open=self.left_open,
                right_open=other.right_open,
            )
            for segment in other.segments:
                out = out + segment
            return out
        return NotImplemented

    def __radd__(self, other: Segment) -> Chain:
        if not isinstance(other, Segment):
            return NotImplemented
        out = Chain(
            chain_id=self.chain_id,
            patch_ids=list(self.patch_ids),
            segments=[],
            chosen_hits=list(self.chosen_hits),
            left_open=self.left_open,
            right_open=self.right_open,
        )
        out = out + other
        for segment in self.segments:
            out = out + segment
        return out


Endpoint = Tuple[str, int, int, str]  # (assembly_contig, coord, direction, chain_id)
LinkMap = Dict[str, Tuple[str, Segment]]  # chain_id -> (other_chain_id, bridge_segment)
ChainGraph = Tuple[Dict[str, Chain], LinkMap, LinkMap, List[List[str]]]


@dataclass(frozen=True)
class PatchStageContext:
    patcher: "GenomePatcher"
    all_hits_by_patch: Dict[str, List[PatchHit]]
    chains: Dict[str, Chain]
    next_link: LinkMap
    prev_link: LinkMap
    paths: List[List[str]]


@dataclass(frozen=True)
class PatchStageResult:
    patched: Dict[str, SeqRecord]
    changed: List[str]
    context: PatchStageContext


class GenomePatcher:
    def __init__(
        self,
        assembly_fasta: str,
        patch_fasta: str,
        bam_path: str,
        min_pid: float = 0.999,
        min_aln_len: int = 10000,
        border_tolerance: int = 10,
        preserve_patch_seq: bool = False,
        verbose: bool = False,
    ) -> None:
        self.assembly_fasta = assembly_fasta
        self.patch_fasta = patch_fasta
        self.bam_path = bam_path
        self.min_pid = min_pid
        self.min_aln_len = min_aln_len
        self.border_tolerance = max(0, int(border_tolerance))
        self.preserve_patch_seq = bool(preserve_patch_seq)
        self.verbose = verbose

        self.assembly_records: Dict[str, SeqRecord] = {}
        self.patch_records: Dict[str, SeqRecord] = {}
        self._hits_cache: Optional[Dict[str, List[PatchHit]]] = None

    def load_fastas(self) -> None:
        self.assembly_records = read_fasta_records(self.assembly_fasta)
        self.patch_records = read_fasta_records(self.patch_fasta)
        if self.verbose:
            print(
                f"[patch] loaded assembly contigs: {len(self.assembly_records)}, "
                f"patch contigs: {len(self.patch_records)}"
            )

    def _iter_alignments(self) -> Iterable[pysam.AlignedSegment]:
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        try:
            for aln in bam.fetch(until_eof=True):
                yield aln
        finally:
            bam.close()

    @staticmethod
    def _pid_from_cigar(cigartuples: Optional[List[Tuple[int, int]]]) -> Tuple[float, int]:
        if not cigartuples:
            return 0.0, 0
        has_eqx = any(op in (7, 8) for op, _ in cigartuples)
        insertions = sum(length for op, length in cigartuples if op == 1)
        deletions  = sum(length for op, length in cigartuples if op == 2)
        if has_eqx:
            matches    = sum(length for op, length in cigartuples if op == 7)
            mismatches = sum(length for op, length in cigartuples if op == 8)
            aligned    = matches + mismatches
        else:
            matches = sum(length for op, length in cigartuples if op == 0)
            aligned = matches
        pid_denom = aligned + insertions + deletions
        pid = (float(matches) / float(pid_denom)) if pid_denom > 0 else 0.0
        return pid, aligned

    @staticmethod
    def _query_aligned_bounds(
        aln: pysam.AlignedSegment, query_len: int
    ) -> Optional[Tuple[int, int]]:
        cigartuples = aln.cigartuples or []
        if not cigartuples or query_len <= 0:
            return None

        left_clip = 0
        i = 0
        while i < len(cigartuples) and cigartuples[i][0] in (4, 5):
            left_clip += cigartuples[i][1]
            i += 1

        right_clip = 0
        j = len(cigartuples) - 1
        while j >= 0 and cigartuples[j][0] in (4, 5):
            right_clip += cigartuples[j][1]
            j -= 1

        if left_clip + right_clip > query_len:
            return None

        if aln.is_reverse:
            start = right_clip
            end = query_len - left_clip
        else:
            start = left_clip
            end = query_len - right_clip

        if start < 0 or end < 0 or start >= end or end > query_len:
            return None
        return start, end

    def collect_hits(self, use_cache: bool = True) -> Dict[str, List[PatchHit]]:
        if use_cache and self._hits_cache is not None:
            return self._hits_cache
        if not self.assembly_records or not self.patch_records:
            self.load_fastas()

        hits_by_assembly: Dict[str, List[PatchHit]] = {}
        total = 0
        unmapped = 0
        non_primary = 0
        missing_metrics = 0
        pid_or_len_filtered = 0
        contig_missing = 0
        bad_coords = 0
        kept = 0
        for aln in self._iter_alignments():
            total += 1
            if aln.is_unmapped:
                unmapped += 1
                continue
            # Keep non-primary alignments as candidates by design.
            if aln.is_secondary or aln.is_supplementary:
                non_primary += 1

            pid, aligned_len = self._pid_from_cigar(aln.cigartuples)
            if aligned_len <= 0:
                missing_metrics += 1
                continue

            assembly_contig = aln.query_name
            patch_contig = aln.reference_name
            if (
                assembly_contig not in self.assembly_records
                or patch_contig not in self.patch_records
            ):
                contig_missing += 1
                continue

            assembly_len = len(self.assembly_records[assembly_contig].seq)
            assembly_bounds = self._query_aligned_bounds(aln, assembly_len)
            if assembly_bounds is None:
                bad_coords += 1
                continue
            assembly_start, assembly_end = assembly_bounds
            patch_start = aln.reference_start
            patch_end = aln.reference_end
            patch_rc = aln.is_reverse

            if assembly_start >= assembly_end or patch_start >= patch_end:
                bad_coords += 1
                continue

            if pid < self.min_pid or aligned_len < self.min_aln_len:
                pid_or_len_filtered += 1
                continue

            hit = PatchHit(
                assembly_contig=assembly_contig,
                assembly_start=assembly_start,
                assembly_end=assembly_end,
                patch_contig=patch_contig,
                patch_start=patch_start,
                patch_end=patch_end,
                patch_rc=patch_rc,
                pid=pid,
                aln_len=aligned_len,
            )
            hits_by_assembly.setdefault(assembly_contig, []).append(hit)
            kept += 1

        for contig in hits_by_assembly:
            hits_by_assembly[contig].sort(key=lambda h: (h.assembly_start, h.assembly_end))

        self._hits_cache = hits_by_assembly
        if self.verbose:
            print(
                "[patch] alignments summary: "
                f"total={total}, kept={kept}, unmapped={unmapped}, non_primary_included={non_primary}, "
                f"missing_metrics={missing_metrics}, low_pid_or_short={pid_or_len_filtered}, "
                f"unknown_contig={contig_missing}, bad_coords={bad_coords}"
            )
            print(f"[patch] assembly contigs with candidate hits: {len(hits_by_assembly)}")
        return hits_by_assembly

    def write_alignment_table(self, out_path: str) -> int:
        if not self.assembly_records or not self.patch_records:
            self.load_fastas()
        hits_by_assembly = self.collect_hits()

        rows: List[Tuple[str, str, int, int, int, int, str, float, float]] = []
        for hits in hits_by_assembly.values():
            for hit in hits:
                contig_len = len(self.assembly_records[hit.assembly_contig].seq)
                aligned_span = max(0, hit.assembly_end - hit.assembly_start)
                aligned_pct = (100.0 * aligned_span / contig_len) if contig_len > 0 else 0.0
                rows.append(
                    (
                        hit.patch_contig,
                        hit.assembly_contig,
                        hit.patch_start,
                        hit.patch_end,
                        hit.assembly_start,
                        hit.assembly_end,
                        "-" if hit.patch_rc else "+",
                        hit.pid,
                        aligned_pct,
                    )
                )
        rows.sort(key=lambda r: (r[0], r[2], r[3], r[1], r[4], r[5]))

        with open(out_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "patch_seq_name",
                    "ctg_name",
                    "p_start",
                    "p_end",
                    "q_start",
                    "q_end",
                    "strand",
                    "pid",
                    "plen",
                ]
            )
            for row in rows:
                writer.writerow(
                    [
                        row[0],
                        row[1],
                        row[2],
                        row[3],
                        row[4],
                        row[5],
                        row[6],
                        f"{row[7]:.6f}",
                        f"{row[8]:.6f}",
                    ]
                )

        return len(rows)

    @staticmethod
    def _best_covering_hit_for_interval(
        hits: List[PatchHit], seg_start: int, seg_end: int
    ) -> Optional[PatchHit]:
        covering = [
            hit
            for hit in hits
            if hit.patch_start <= seg_start and hit.patch_end >= seg_end
        ]
        if not covering:
            return None
        covering.sort(
            key=lambda h: (
                h.pid,
                h.aln_len,
                -h.patch_start,
                h.patch_end,
                h.assembly_contig,
                h.assembly_start,
                h.assembly_end,
            ),
            reverse=True,
        )
        return covering[0]

    @staticmethod
    def _project_interval_to_assembly_segment(
        hit: PatchHit, seg_start: int, seg_end: int
    ) -> Segment:
        if seg_start < hit.patch_start or seg_end > hit.patch_end or seg_start >= seg_end:
            raise ValueError(
                "cannot project patch interval outside hit bounds: "
                f"patch={hit.patch_contig}:{seg_start}-{seg_end}, "
                f"hit={hit.patch_start}-{hit.patch_end}"
            )

        # Preserve exact aligned assembly bounds when taking the full hit span.
        if seg_start == hit.patch_start and seg_end == hit.patch_end:
            return Segment(
                source="assembly",
                contig=hit.assembly_contig,
                start=hit.assembly_start,
                end=hit.assembly_end,
                reverse=hit.patch_rc,
            )

        if hit.patch_rc:
            asm_start = hit.assembly_end - (seg_end - hit.patch_start)
            asm_end = hit.assembly_end - (seg_start - hit.patch_start)
        else:
            asm_start = hit.assembly_start + (seg_start - hit.patch_start)
            asm_end = hit.assembly_start + (seg_end - hit.patch_start)

        if asm_start < 0 or asm_end < 0 or asm_start >= asm_end:
            raise ValueError(
                "invalid projected assembly interval: "
                f"{hit.assembly_contig}:{asm_start}-{asm_end} "
                f"from patch interval {seg_start}-{seg_end}"
            )

        return Segment(
            source="assembly",
            contig=hit.assembly_contig,
            start=asm_start,
            end=asm_end,
            reverse=hit.patch_rc,
        )

    def _hits_by_patch(self) -> Dict[str, List[PatchHit]]:
        all_hits_by_patch: Dict[str, List[PatchHit]] = {}
        for assembly_hits in self.collect_hits().values():
            for hit in assembly_hits:
                all_hits_by_patch.setdefault(hit.patch_contig, []).append(hit)
        for patch_id in all_hits_by_patch:
            all_hits_by_patch[patch_id].sort(key=lambda h: (h.patch_start, h.patch_end))
        return all_hits_by_patch

    def _build_patch_chains(
        self, all_hits_by_patch: Optional[Dict[str, List[PatchHit]]] = None
    ) -> Dict[str, Chain]:
        if not self.assembly_records or not self.patch_records:
            self.load_fastas()
        if all_hits_by_patch is None:
            all_hits_by_patch = self._hits_by_patch()

        chains: Dict[str, Chain] = {}
        for patch_id, patch_rec in self.patch_records.items():
            patch_hits = list(all_hits_by_patch.get(patch_id, []))
            patch_hits.sort(key=lambda h: (h.patch_start, h.patch_end))
            if not patch_hits:
                continue

            patch_len = len(patch_rec.seq)
            boundaries = {0, patch_len}
            for hit in patch_hits:
                if hit.patch_start < patch_len and hit.patch_end > 0:
                    boundaries.add(max(0, min(patch_len, hit.patch_start)))
                    boundaries.add(max(0, min(patch_len, hit.patch_end)))
            sorted_bounds = sorted(boundaries)

            chain_id = f"patch::{patch_id}"
            chain = Chain(chain_id=chain_id, patch_ids=[patch_id])
            left_terminal_hit: Optional[PatchHit] = None
            right_terminal_hit: Optional[PatchHit] = None
            left_terminal_seg_start: Optional[int] = None
            right_terminal_seg_end: Optional[int] = None
            for seg_start, seg_end in zip(sorted_bounds, sorted_bounds[1:]):
                if seg_end <= seg_start:
                    continue
                best_hit = self._best_covering_hit_for_interval(
                    patch_hits, seg_start, seg_end
                )
                if best_hit is None:
                    chain = chain + Segment(
                        source="patch",
                        contig=patch_id,
                        start=seg_start,
                        end=seg_end,
                    )
                    continue
                chain.chosen_hits.append(best_hit)
                if left_terminal_seg_start is None or seg_start < left_terminal_seg_start:
                    left_terminal_hit = best_hit
                    left_terminal_seg_start = seg_start
                if right_terminal_seg_end is None or seg_end > right_terminal_seg_end:
                    right_terminal_hit = best_hit
                    right_terminal_seg_end = seg_end
                if self.preserve_patch_seq:
                    chain = chain + Segment(
                        source="patch",
                        contig=patch_id,
                        start=seg_start,
                        end=seg_end,
                    )
                else:
                    chain = chain + (
                        self._project_interval_to_assembly_segment(
                            best_hit,
                            seg_start,
                            seg_end,
                        )
                    )

            left_open: Optional[Tuple[str, int, int]] = None
            if (
                left_terminal_hit is not None
                and left_terminal_seg_start is not None
                and left_terminal_seg_start <= self.border_tolerance
            ):
                left_open = self._endpoint_from_border_hit(chain_id, "left", left_terminal_hit)

            right_open: Optional[Tuple[str, int, int]] = None
            if (
                right_terminal_hit is not None
                and right_terminal_seg_end is not None
                and right_terminal_seg_end >= patch_len - self.border_tolerance
            ):
                right_open = self._endpoint_from_border_hit(chain_id, "right", right_terminal_hit)

            # In substitution mode, open borders are resolved by assembly continuation/linking.
            # In preserve mode, keep full patch sequence segments.
            if not self.preserve_patch_seq:
                if left_open and chain.segments:
                    first = chain.segments[0]
                    if first.source == "patch" and first.contig == patch_id and first.start == 0:
                        chain.segments = chain.segments[1:]
                if right_open and chain.segments:
                    last = chain.segments[-1]
                    if (
                        last.source == "patch"
                        and last.contig == patch_id
                        and last.end == patch_len
                    ):
                        chain.segments = chain.segments[:-1]

            chain.left_open = left_open
            chain.right_open = right_open
            chains[chain_id] = chain

        return chains

    @staticmethod
    def _endpoint_from_border_hit(
        chain_id: str, side: str, hit: PatchHit
    ) -> Tuple[str, int, int]:
        del chain_id
        if side not in ("left", "right"):
            raise ValueError(f"unknown endpoint side: {side}")

        if side == "left":
            if hit.patch_rc:
                coord = hit.assembly_end
                direction = +1
            else:
                coord = hit.assembly_start
                direction = -1
        else:
            if hit.patch_rc:
                coord = hit.assembly_start
                direction = -1
            else:
                coord = hit.assembly_end
                direction = +1

        return hit.assembly_contig, coord, direction

    @staticmethod
    def _sorted_endpoints(
        chains: Dict[str, Chain],
    ) -> Tuple[List[Endpoint], List[Endpoint]]:
        left_endpoints = sorted(
            [
                (left[0], left[1], left[2], chain_id)
                for chain_id, chain in chains.items()
                for left in [chain.left_open]
                if left is not None
            ],
            key=lambda e: (e[0], e[1], e[3]),
        )
        right_endpoints = sorted(
            [
                (right[0], right[1], right[2], chain_id)
                for chain_id, chain in chains.items()
                for right in [chain.right_open]
                if right is not None
            ],
            key=lambda e: (e[0], e[1], e[3]),
        )
        return left_endpoints, right_endpoints

    def _resolve_chain_links(
        self, chains: Dict[str, Chain]
    ) -> Tuple[LinkMap, LinkMap]:
        left_endpoints, right_endpoints = self._sorted_endpoints(chains)

        left_by_asm: Dict[str, List[Endpoint]] = {}
        for endpoint in left_endpoints:
            left_by_asm.setdefault(endpoint[0], []).append(endpoint)

        next_link: LinkMap = {}
        prev_link: LinkMap = {}
        for right in right_endpoints:
            right_contig, right_coord, right_direction, right_chain = right
            if right_chain in next_link:
                continue

            same_asm = [
                left
                for left in left_by_asm.get(right_contig, [])
                if left[3] != right_chain and left[3] not in prev_link
            ]
            if not same_asm:
                continue

            valid = [left for left in same_asm if right_coord <= left[1]]
            valid = [left for left in valid if right_direction == +1 and left[2] == -1]
            if not valid:
                considered_right = [
                    ep for ep in right_endpoints if ep[0] == right_contig
                ]
                considered_right_fmt = ", ".join(
                    f"{ep[3]}:{ep[1]}" for ep in considered_right
                )
                considered_left_fmt = ", ".join(
                    f"{ep[3]}:{ep[1]}" for ep in same_asm
                )
                raise ValueError(
                    "invalid open-link partners: "
                    f"right {right_chain}:{right_contig}:{right_coord} "
                    "has same-contig left candidates but none satisfy right<=left and direction. "
                    f"assembly_contig={right_contig}; "
                    f"right_open_coords=[{considered_right_fmt}]; "
                    f"left_open_coords=[{considered_left_fmt}]"
                )

            min_gap = min(left[1] - right_coord for left in valid)
            closest = [left for left in valid if (left[1] - right_coord) == min_gap]
            if len(closest) != 1:
                raise ValueError(
                    "ambiguous closest partner for endpoint: "
                    f"right {right_chain}:{right_contig}:{right_coord}"
                )
            left = closest[0]
            left_chain = left[3]

            if left_chain in prev_link or right_chain in next_link:
                raise ValueError(
                    "endpoint reuse detected while linking chains: "
                    f"right={right_chain}, left={left_chain}"
                )

            bridge = Segment(
                source="assembly",
                contig=right_contig,
                start=right_coord,
                end=left[1],
            )
            next_link[right_chain] = (left_chain, bridge)
            prev_link[left_chain] = (right_chain, bridge)

        return next_link, prev_link

    @staticmethod
    def _build_chain_paths(
        chains: Dict[str, Chain],
        next_link: LinkMap,
        prev_link: LinkMap,
    ) -> List[List[str]]:
        starts = sorted(chain_id for chain_id in chains if chain_id not in prev_link)
        if not starts and chains:
            raise ValueError("no chain start found; cycle detected in link graph")

        paths: List[List[str]] = []
        visited: set[str] = set()
        for start in starts:
            path: List[str] = []
            seen_local: set[str] = set()
            curr = start
            while True:
                if curr in seen_local:
                    raise ValueError(f"cycle detected while walking chain path from {start}")
                seen_local.add(curr)
                path.append(curr)
                visited.add(curr)
                link = next_link.get(curr)
                if link is None:
                    break
                curr = link[0]
            paths.append(path)

        missing = [chain_id for chain_id in chains if chain_id not in visited]
        if missing:
            raise ValueError(f"unreachable chains after linking: {', '.join(sorted(missing))}")
        return paths

    def _build_chain_graph(
        self, all_hits_by_patch: Optional[Dict[str, List[PatchHit]]] = None
    ) -> ChainGraph:
        chains = self._build_patch_chains(all_hits_by_patch=all_hits_by_patch)
        next_link, prev_link = self._resolve_chain_links(chains)
        paths = self._build_chain_paths(chains, next_link, prev_link)
        return chains, next_link, prev_link, paths

    def build_stage_context(self) -> PatchStageContext:
        all_hits_by_patch = self._hits_by_patch()
        chains, next_link, prev_link, paths = self._build_chain_graph(
            all_hits_by_patch=all_hits_by_patch
        )
        return PatchStageContext(
            patcher=self,
            all_hits_by_patch=all_hits_by_patch,
            chains=chains,
            next_link=next_link,
            prev_link=prev_link,
            paths=paths,
        )

    def _assembly_len(self, contig_id: str) -> int:
        return len(self.assembly_records[contig_id].seq)

    @staticmethod
    def _chain_label(path: List[str], chains: Dict[str, Chain], chain_idx: int) -> str:
        del chain_idx
        patch_ids: List[str] = []
        for chain_id in path:
            patch_ids.extend(chains[chain_id].patch_ids)
        joined = "_".join(patch_ids)
        return f"patched_{joined}"

    def _segment_sequence(self, segment: Segment) -> str:
        if segment.end <= segment.start:
            return ""
        if segment.source == "patch":
            seq = self.patch_records[segment.contig].seq[segment.start : segment.end]
        elif segment.source == "assembly":
            seq = self.assembly_records[segment.contig].seq[segment.start : segment.end]
            if segment.reverse:
                seq = seq.reverse_complement()
        else:
            raise ValueError(f"unknown segment source: {segment.source}")
        return str(seq)

    def _emit_chain_record(
        self,
        path: List[str],
        chains: Dict[str, Chain],
        next_link: LinkMap,
        prev_link: LinkMap,
        chain_idx: int,
    ) -> Tuple[SeqRecord, set[str]]:
        used_assembly: set[str] = set()
        seq_parts: List[str] = []

        first = chains[path[0]]
        if self.preserve_patch_seq:
            used_assembly.update(hit.assembly_contig for hit in first.chosen_hits)
        if first.left_open and first.chain_id not in prev_link:
            left_contig, left_coord, left_direction = first.left_open
            if left_coord < 0 or left_coord > self._assembly_len(left_contig):
                raise ValueError(
                    f"invalid left open coordinate: {left_contig}:{left_coord}"
                )
            asm_len = self._assembly_len(left_contig)
            if left_direction == -1:
                open_left_segment = Segment(
                    source="assembly",
                    contig=left_contig,
                    start=0,
                    end=left_coord,
                    reverse=False,
                )
            else:
                open_left_segment = Segment(
                    source="assembly",
                    contig=left_contig,
                    start=left_coord,
                    end=asm_len,
                    reverse=True,
                )
            seq_parts.append(self._segment_sequence(open_left_segment))
            used_assembly.add(left_contig)

        for segment in first.segments:
            seq_parts.append(self._segment_sequence(segment))
            if segment.source == "assembly":
                used_assembly.add(segment.contig)

        for chain_id in path[:-1]:
            link = next_link.get(chain_id)
            if link is None:
                break
            to_chain, bridge = link
            if bridge.end < bridge.start:
                raise ValueError(
                    "invalid link coordinates: "
                    f"{bridge.contig}:{bridge.start}>{bridge.end} "
                    f"(from={chain_id}, to={to_chain})"
                )
            if bridge.end > bridge.start:
                seq_parts.append(self._segment_sequence(bridge))
                used_assembly.add(bridge.contig)

            nxt = chains[to_chain]
            if self.preserve_patch_seq:
                used_assembly.update(hit.assembly_contig for hit in nxt.chosen_hits)
            for segment in nxt.segments:
                seq_parts.append(self._segment_sequence(segment))
                if segment.source == "assembly":
                    used_assembly.add(segment.contig)

        last = chains[path[-1]]
        if last.right_open and last.chain_id not in next_link:
            right_contig, right_coord, right_direction = last.right_open
            asm_len = self._assembly_len(right_contig)
            if right_coord < 0 or right_coord > asm_len:
                raise ValueError(
                    f"invalid right open coordinate: {right_contig}:{right_coord}"
                )
            if right_direction == +1:
                open_right_segment = Segment(
                    source="assembly",
                    contig=right_contig,
                    start=right_coord,
                    end=asm_len,
                    reverse=False,
                )
            else:
                open_right_segment = Segment(
                    source="assembly",
                    contig=right_contig,
                    start=0,
                    end=right_coord,
                    reverse=True,
                )
            seq_parts.append(self._segment_sequence(open_right_segment))
            used_assembly.add(right_contig)

        contig_id = self._chain_label(path, chains, chain_idx)
        seq = "".join(seq_parts)
        patch_labels = ",".join(pid for cid in path for pid in chains[cid].patch_ids)
        rec = SeqRecord(Seq(seq), id=contig_id, description=f"patched_chain patches={patch_labels}")
        return rec, used_assembly

    def patch_all(self, chain_graph: Optional[ChainGraph] = None) -> Tuple[Dict[str, SeqRecord], List[str]]:
        if not self.assembly_records or not self.patch_records:
            self.load_fastas()

        if chain_graph is None:
            chains, next_link, prev_link, paths = self._build_chain_graph()
        else:
            chains, next_link, prev_link, paths = chain_graph

        out_records: Dict[str, SeqRecord] = {}
        used_assembly: set[str] = set()
        for idx, path in enumerate(paths, start=1):
            rec, touched = self._emit_chain_record(path, chains, next_link, prev_link, idx)
            out_records[rec.id] = rec
            used_assembly.update(touched)

        for contig_id, rec in self.assembly_records.items():
            if contig_id not in used_assembly:
                out_records[contig_id] = rec

        changed = sorted(used_assembly)
        if self.verbose:
            print(
                f"[patch] emitted chains={len(paths)}, "
                f"removed_assembly_contigs={len(changed)} (used in emitted chains), "
                f"kept_assembly_contigs={len(self.assembly_records) - len(changed)}, "
                f"assembly_used_in_emitted_chains={len(used_assembly)}"
            )
        return out_records, changed


def run_patch_stage(
    assembly_fasta: str,
    patch_fasta: str,
    bam_path: str,
    out_fasta: str,
    out_changed: Optional[str] = None,
    out_alignments_tsv: Optional[str] = None,
    min_pid: float = 0.999,
    min_aln_len: int = 1000,
    border_tolerance: int = 10,
    preserve_patch_seq: bool = False,
    verbose: bool = False,
) -> PatchStageResult:
    patcher = GenomePatcher(
        assembly_fasta=assembly_fasta,
        patch_fasta=patch_fasta,
        bam_path=bam_path,
        min_pid=min_pid,
        min_aln_len=min_aln_len,
        border_tolerance=border_tolerance,
        preserve_patch_seq=preserve_patch_seq,
        verbose=verbose,
    )
    alignments_written: Optional[int] = None
    if out_alignments_tsv:
        alignments_written = patcher.write_alignment_table(out_alignments_tsv)
        if verbose:
            print(
                f"[patch] wrote alignment table: {out_alignments_tsv}; "
                f"rows={alignments_written}"
            )
    context = patcher.build_stage_context()
    chain_graph: ChainGraph = (context.chains, context.next_link, context.prev_link, context.paths)
    patched, changed = patcher.patch_all(chain_graph=chain_graph)
    write_fasta(patched, out_fasta)
    if out_changed:
        write_changed_list(changed, out_changed)
    if verbose:
        print(
            f"[patch] wrote patched FASTA: {out_fasta}; "
            f"removed assembly contigs: {len(changed)}; changed list: {out_changed or '<disabled>'}"
        )
    return PatchStageResult(
        patched=patched,
        changed=changed,
        context=context,
    )
