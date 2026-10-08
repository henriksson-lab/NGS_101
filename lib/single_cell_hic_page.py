"""Compact renderer shared by paper-specific single-cell 3C modules."""
from __future__ import annotations

import seqprimers as sp
from chemdraw import Scene, annotation_rows, junction_row, panel
from page import head, info


def render(module) -> str:
    w = module.WORKFLOW
    contact = module.contact_product()
    junctions = getattr(module, "JUNCTIONS",
                        (("locus A", "contact junction", "ligation"),))
    contact_rows = [*Scene.duplex(list(contact), label="one captured contact").rows(),
                    *(junction_row(contact, left, right, label)
                      for left, right, label in junctions),
                    *annotation_rows(contact)]
    final = module.final_library()
    final_rows = [*Scene.duplex(list(final), label="sequencing library").rows(),
                  *annotation_rows(final)]
    steps = "".join(f'<h3>({i}) {title}</h3><p>{text}</p>'
                    for i, (title, text) in enumerate(module.STEPS, 1))
    primers = getattr(module, "SEQ_PRIMERS", ())
    required = getattr(module, "REQUIRED_ROLES", ())
    if primers:
        sequencing = sp.section(final, primers, heading="Sequencing",
            intro=module.READOUT, required_roles=required)
    else:
        sequencing = (f'<h2>Sequencing</h2><p>{module.READOUT}</p>' +
            sp.unavailable_diagram(final, module.SEQUENCING_UNAVAILABLE))
    return "\n".join([
        head(module.TITLE), '<div class="wrap">', f'<h1>{module.TITLE}</h1>',
        f'<p class="research-notes"><a href="{module.NOTES}">Research notes</a></p>',
        info(module.CITATION), info(module.SUMMARY),
        info(module.LIBRARY_CAVEAT) if getattr(module, "LIBRARY_CAVEAT", "") else "",
        '<h2>Protocol path</h2>', steps,
        '<h2>The proximity-ligation product</h2>',
        panel(contact_rows, cls="small", caption=module.JUNCTION_CAPTION),
        '<h2>Final sequencing molecule</h2>',
        panel(final_rows, cls="small", caption=module.LIBRARY_CAPTION),
        sequencing, '</div>'])
