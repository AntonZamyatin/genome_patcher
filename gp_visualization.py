from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List

from gp_patching import Chain, GenomePatcher, LinkMap, PatchHit, PatchStageContext, Segment


def _safe_name(raw: str) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", raw).strip("._")
    return name or "unnamed"


def _plot_single_patch(
    patch_id: str,
    patch_len: int,
    hits: List[PatchHit],
    out_png: Path,
    title: str,
) -> None:
    import matplotlib.pyplot as plt

    hits = sorted(hits, key=lambda h: (h.patch_start, h.patch_end))
    unique_contigs = sorted({h.assembly_contig for h in hits})
    cmap = plt.get_cmap("tab20")
    contig_colors: Dict[str, tuple] = {
        contig: cmap(i % 20) for i, contig in enumerate(unique_contigs)
    }

    fig_h = max(2.8, 1.5 + 0.75 * max(1, len(hits)))
    fig, ax = plt.subplots(figsize=(16, fig_h))

    baseline_y = 0.0
    ax.hlines(
        y=baseline_y,
        xmin=0,
        xmax=max(1, patch_len),
        color="black",
        linewidth=3,
    )
    ax.text(
        0,
        baseline_y + 0.2,
        f"{patch_id} (len={patch_len})",
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    row_step = 1.05
    for i, hit in enumerate(hits):
        y = baseline_y - (i + 1) * row_step
        color = contig_colors[hit.assembly_contig]
        ax.hlines(y=y, xmin=hit.patch_start, xmax=hit.patch_end, color=color, linewidth=7.5, alpha=0.92)

        if hit.patch_rc:
            x0, x1 = hit.patch_end, hit.patch_start
        else:
            x0, x1 = hit.patch_start, hit.patch_end
        ax.annotate(
            "",
            xy=(x1, y),
            xytext=(x0, y),
            arrowprops=dict(
                arrowstyle="-|>",
                color=color,
                lw=2.8,
                mutation_scale=24,
                shrinkA=0,
                shrinkB=0,
            ),
        )

        label = f"{hit.assembly_contig}:{hit.assembly_start}-{hit.assembly_end}"
        label_x = (hit.patch_start + hit.patch_end) / 2.0
        ax.text(label_x, y + 0.13, label, fontsize=8, color=color, ha="center", va="bottom")

    x_max = max(1, patch_len)
    ax.set_xlim(0, x_max)
    y_min = baseline_y - (len(hits) + 1.5) * row_step
    ax.set_ylim(y_min, baseline_y + 0.9)
    ax.set_xlabel("Patch sequence coordinate")
    ax.set_yticks([])
    ax.set_title(title)
    ax.grid(True, axis="x", linestyle="--", alpha=0.35)

    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def _collect_chain_segments(
    patcher: GenomePatcher,
    path: List[str],
    chains: Dict[str, Chain],
    next_link: LinkMap,
    prev_link: LinkMap,
) -> List[Segment]:
    segments: List[Segment] = []

    def append_segment(segment: Segment) -> None:
        if segment.end <= segment.start:
            return
        if not segments:
            segments.append(segment)
            return
        try:
            segments[-1] = segments[-1] + segment
        except TypeError:
            segments.append(segment)

    first = chains[path[0]]
    if first.left_open and first.chain_id not in prev_link:
        left_contig, left_coord, left_direction = first.left_open
        asm_len = patcher._assembly_len(left_contig)
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
        append_segment(open_left_segment)

    for segment in first.segments:
        append_segment(segment)

    for chain_id in path[:-1]:
        link = next_link.get(chain_id)
        if link is None:
            break
        to_chain, bridge = link
        append_segment(bridge)

        nxt = chains[to_chain]
        for segment in nxt.segments:
            append_segment(segment)

    last = chains[path[-1]]
    if last.right_open and last.chain_id not in next_link:
        right_contig, right_coord, right_direction = last.right_open
        asm_len = patcher._assembly_len(right_contig)
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
        append_segment(open_right_segment)

    return segments


def _segment_caption(segment: Segment) -> str:
    direction = "-" if segment.reverse else "+"
    seg_len = max(0, segment.end - segment.start)
    return f"{segment.contig}:{segment.start}-{segment.end} {direction} ({seg_len})"


def _plot_patched_chain(
    contig_id: str,
    segments: List[Segment],
    out_png: Path,
    title: str,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrow

    plot_segments = [s for s in segments if s.end > s.start]
    unique_contigs = sorted({s.contig for s in plot_segments if s.source == "assembly"})
    cmap = plt.get_cmap("tab20")
    contig_colors: Dict[str, tuple] = {
        contig: cmap(i % 20) for i, contig in enumerate(unique_contigs)
    }

    n = max(1, len(plot_segments))
    fig_h = max(4.0, 1.0 + 0.8 * n)
    fig, ax = plt.subplots(figsize=(16, fig_h))

    arrow_x = 0.30
    arrow_len = 0.72
    coord_x = 1.42
    coord_text_x = 1.48
    segment_start_coords: List[int] = []
    cursor = 0
    for idx, segment in enumerate(plot_segments):
        y = n - idx
        color = "black" if segment.source == "patch" else contig_colors[segment.contig]
        segment_start_coords.append(cursor)
        cursor += segment.end - segment.start

        if segment.reverse:
            start_y = y - (arrow_len / 2.0)
            dy = +arrow_len
        else:
            start_y = y + (arrow_len / 2.0)
            dy = -arrow_len
        ax.add_patch(
            FancyArrow(
                arrow_x,
                start_y,
                0.0,
                dy,
                width=0.14,
                head_width=0.28,
                head_length=0.16,
                length_includes_head=True,
                color=color,
                alpha=0.95,
            )
        )

        ax.text(
            0.50,
            y,
            _segment_caption(segment),
            fontsize=8,
            color=color,
            ha="left",
            va="center",
        )

    # Coordinate rail: segment starts in resulting patched-contig coordinates.
    ax.vlines(coord_x, ymin=0.5, ymax=n + 0.5, colors="black", linewidth=1.2, alpha=0.75)
    for idx, start_coord in enumerate(segment_start_coords):
        y = n - idx
        ax.hlines(y=y, xmin=coord_x - 0.02, xmax=coord_x + 0.02, colors="black", linewidth=1.2)
        ax.text(
            coord_text_x,
            y,
            f"{start_coord}",
            fontsize=8,
            color="black",
            ha="left",
            va="center",
        )

    total_len = max(0, cursor)

    ax.set_xlim(0, 2.2)
    ax.set_ylim(0, n + 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    ax.text(0.02, n + 0.35, f"{contig_id} (len={total_len})", fontsize=10, fontweight="bold")
    ax.text(coord_x, n + 0.35, "segment_start", fontsize=8, color="black", ha="left")
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def generate_patch_plots(
    context: PatchStageContext,
    out_dir: str,
    label: str,
    verbose: bool = False,
) -> List[str]:
    patcher = context.patcher
    all_by_patch = context.all_hits_by_patch

    plots_dir = Path(out_dir) / "patch_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    if verbose:
        print(f"[viz] writing per-patch plots to: {plots_dir}")
        print(f"[viz] patch sequences total: {len(patcher.patch_records)}")

    output_files: List[str] = []
    for patch_id, patch_rec in patcher.patch_records.items():
        all_hits = [
            hit
            for hit in all_by_patch.get(patch_id, [])
            if hit.pid >= patcher.min_pid and hit.aln_len >= patcher.min_aln_len
        ]
        out_all_png = plots_dir / f"{_safe_name(label)}.{_safe_name(patch_id)}.all_hits.png"

        _plot_single_patch(
            patch_id=patch_id,
            patch_len=len(patch_rec.seq),
            hits=all_hits,
            out_png=out_all_png,
            title=(
                f"All patch-stage hits on patch: {patch_id} "
                f"(pid>={patcher.min_pid}, aln_len>={patcher.min_aln_len})"
            ),
        )
        output_files.append(str(out_all_png))
        if verbose:
            print(
                f"[viz] {patch_id}: all_hits={len(all_hits)} -> {out_all_png}"
            )

    chains = context.chains
    next_link = context.next_link
    prev_link = context.prev_link
    paths = context.paths

    for idx, path in enumerate(paths, start=1):
        chain_id = patcher._chain_label(path, chains, idx)
        chain_segments = _collect_chain_segments(
            patcher=patcher,
            path=path,
            chains=chains,
            next_link=next_link,
            prev_link=prev_link,
        )
        out_chain_png = plots_dir / f"{_safe_name(label)}.{_safe_name(chain_id)}.patched_chain.png"
        patch_labels = ",".join(pid for cid in path for pid in chains[cid].patch_ids)
        _plot_patched_chain(
            contig_id=chain_id,
            segments=chain_segments,
            out_png=out_chain_png,
            title=f"Patched contig chain: {chain_id} (patches={patch_labels})",
        )
        output_files.append(str(out_chain_png))
        if verbose:
            print(
                f"[viz] {chain_id}: segments={len(chain_segments)} -> {out_chain_png}"
            )

    if verbose:
        print(f"[viz] generated {len(output_files)} plot(s)")
    return output_files
