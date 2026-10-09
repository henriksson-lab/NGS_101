"""Diagram-first renderer for NEB commercial protocol modules."""
from __future__ import annotations

import seqprimers as sp
from chemdraw import Construct, Scene, annotation_rows, duplex_rows, panel, strand_row
from page import head, info


def _rows(obj):
    if isinstance(obj, Scene):
        return obj.rows()
    if isinstance(obj, Construct):
        return [strand_row(obj), *annotation_rows(obj)]
    return obj


def render(module) -> str:
    chunks = [head(module.TITLE), '<div class="wrap">', f'<h1>{module.TITLE}</h1>',
              f'<p class="research-notes"><a href="{module.NOTES}">Research notes</a></p>',
              info(module.SOURCE), info(module.SUMMARY)]
    caveat = getattr(module, "CAVEAT", "")
    if caveat:
        chunks.append(f'<div class="caveat"><b>INFERRED boundary.</b> {caveat}</div>')
    for i, (heading, obj, caption) in enumerate(module.STEPS, 1):
        chunks.extend([f'<h2>({i}) {heading}</h2>', panel(_rows(obj), cls="small", caption=caption)])
    lib = module.final_library()
    chunks.extend(['<h2>Final sequencing molecule</h2>',
                   panel(duplex_rows(lib, label="library"),
                         cls="long", caption=module.FINAL_CAPTION),
                   sp.section(lib, module.SEQ_PRIMERS, intro=module.READOUT,
                              required_roles=module.REQUIRED_ROLES,
                              read_lengths=getattr(module, "READ_LENGTHS", None)), '</div>'])
    return "\n".join(chunks)
