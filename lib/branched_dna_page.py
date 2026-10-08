"""Sparse diagram-first page renderer for branched-DNA RNA detection protocols."""
from __future__ import annotations

from branched_dna import target_scene, tree_rows
from chemdraw import Row, panel
from page import head, info


def render(module) -> str:
    model = module.MODEL
    lo, hi = model.labels_per_target
    total = f"{lo:,}" if lo == hi else f"{lo:,}–{hi:,}"
    chunks = [head(module.TITLE), '<div class="wrap">', f'<h1>{module.TITLE}</h1>',
              f'<p class="research-notes"><a href="{module.NOTES}">Research notes</a></p>',
              info(module.SOURCE), info(module.SUMMARY),
              '<div class="caveat"><b>Sequence boundary.</b> Target, preamplifier, amplifier '
              'and label-probe sequences are proprietary. Role-labelled placeholders show '
              'only the published hybridization architecture.</div>',
              '<h2>(1) Prepare the specimen</h2>',
              panel([Row(chunks=[(module.PREPARATION, None, False)])],
                    caption=module.PREPARATION_CAPTION),
              '<h2>(2) Hybridize adjacent target-probe pairs</h2>',
              panel(target_scene().rows(), cls="small", caption=module.TARGET_CAPTION),
              '<h2>(3) Assemble the branched-DNA signal tree</h2>',
              panel(tree_rows(model), cls="small",
                    caption=(f"One tree carries {model.labels_per_tree} label-probe sites; "
                             f"the documented target-probe range gives {total} theoretical "
                             "label sites per target RNA.")),
              '<h2>(4) Detect the signal</h2>',
              panel([Row(chunks=[(module.DETECTION, "umi", True)])],
                    caption=module.DETECTION_CAPTION),
              '<h2>Measurement, not sequencing</h2>', info(module.READOUT), '</div>']
    return "\n".join(chunks)
