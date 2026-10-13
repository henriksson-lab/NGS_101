"""Sparse diagram-first page renderer for branched-DNA RNA detection protocols."""
from __future__ import annotations

from branched_dna import target_scene, tree_panel
from chemdraw import Row, oligo, panel
from page import head, info
import seqprimers as sp


def render(module) -> str:
    model = module.MODEL
    lo, hi = model.labels_per_target
    total = f"{lo:,}" if lo == hi else f"{lo:,}–{hi:,}"
    boundary = getattr(
        module, "SEQUENCE_BOUNDARY",
        "Target, preamplifier, amplifier and label-probe sequences are proprietary. "
        "Role-labelled placeholders show only the published hybridization architecture.")
    chunks = [head(module.TITLE), '<div class="wrap">', f'<h1>{module.TITLE}</h1>',
              f'<p class="research-notes"><a href="{module.NOTES}">Research notes</a></p>',
              info(module.SOURCE), info(module.SUMMARY),
              f'<div class="caveat"><b>Sequence boundary.</b> {boundary}</div>',
              '<h2>(1) Prepare the specimen</h2>',
              panel([Row(chunks=[(module.PREPARATION, None, False)])],
                    caption=module.PREPARATION_CAPTION),
              '<h2>(2) Hybridize adjacent target-probe pairs</h2>',
              panel(target_scene().rows(), cls="small", caption=module.TARGET_CAPTION)]
    if hasattr(module, "sequence_example"):
        example = module.sequence_example()
        chunks.extend([
            f'<h3>{example.title}</h3>', info(example.note), '<seq>',
            *(oligo(name, segments) for name, segments in example.oligos), '</seq>',
            panel(example.scene.rows(), cls="small", caption=example.caption),
        ])
    chunks.extend([
              '<h2>(3) Assemble the branched-DNA signal tree</h2>',
              tree_panel(model,
                         f"Five representative amplifier arms are drawn without crossings; "
                         f"labels give the published multiplicities. One tree carries "
                         f"{model.labels_per_tree} label-probe sites, and the documented "
                         f"target-probe range gives {total} theoretical label sites per "
                         "target RNA."),
              '<h2>(4) Detect the signal</h2>',
              panel([Row(chunks=[(module.DETECTION, "umi", True)])],
                    caption=module.DETECTION_CAPTION),
              '<h2>Measurement, not sequencing</h2>', info(module.READOUT),
              sp.read_spec_table([
                  sp.read_spec("Sequencing output", None,
                               "No sequencing library or sequencing reads", "not_applicable")
              ], "This is a branched-DNA imaging/flow-cytometry assay, not an NGS library."),
              '</div>'])
    return "\n".join(chunks)
