"""Diagram-first renderer shared by the eight protocol-class gap pages."""
from __future__ import annotations

import seqprimers as sp
from chemdraw import Scene, annotation_rows, panel
from page import head, info


def render(m) -> str:
    parts = [head(m.TITLE), '<div class="wrap">', f'<h1>{m.TITLE}</h1>',
             f'<p class="research-notes"><a href="{m.NOTES}">Research notes</a></p>',
             info(m.SOURCE), info(m.SUMMARY)]
    if getattr(m, "CAVEAT", ""):
        parts.append(f'<div class="caveat">{m.CAVEAT}</div>')
    for heading, rows, caption in m.sections():
        parts.extend([f'<h2>{heading}</h2>', panel(rows, cls="long", caption=caption)])
    lib = getattr(m, "FINAL_LIBRARY", None)
    primers = getattr(m, "SEQ_PRIMERS", ())
    if lib is not None:
        # Bracketed role tokens deliberately stand for unavailable molecular sequence.
        # Scene still owns placement; excluding those prose columns from base-pair
        # validation avoids pretending that a literal "[adapter]" has a complement.
        unknown = tuple(s.name + "'" for s in lib if s.is_role_token())
        parts.extend(['<h2>Final sequencing library</h2>',
                      panel([*Scene.duplex(list(lib), label="library",
                                           unpaired=unknown).rows(),
                             *annotation_rows(lib)], cls="long",
                            caption=m.FINAL_CAPTION)])
    if primers:
        parts.append(sp.section(lib, primers, intro=m.SEQUENCING_INTRO,
                                required_roles=tuple(p.role for p in primers)))
    else:
        ending = getattr(m, "SEQUENCING_ENDING", "")
        if not ending:
            raise ValueError(f"{m.TITLE}: page has neither sequencing primers nor a declared endpoint")
        parts.extend(['<h2>Sequencing entry</h2>', info(ending)])
        unavailable = getattr(m, "SEQUENCING_UNAVAILABLE", "")
        if unavailable:
            parts.append(sp.unavailable_diagram(lib, unavailable))
    parts.extend(['</div>'])
    return "\n".join(parts)
