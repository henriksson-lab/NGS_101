"""Compact renderer shared by paper-specific single-cell 3C modules."""
from __future__ import annotations

import seqprimers as sp
from chemdraw import duplex_rows, junction_row, panel
from page import head, info


def render(module) -> str:
    w = module.WORKFLOW
    contact = module.contact_product()
    junctions = getattr(module, "JUNCTIONS",
                        (("locus A", "contact junction", "ligation"),))
    contact_label = "one captured contact"
    contact_duplex = duplex_rows(contact, label=contact_label)
    sequence_column = len(contact_label) + 1 + len("5'- ")
    contact_rows = [*contact_duplex[:2],
                    *(junction_row(contact, left, right, label,
                                   prefix_width=sequence_column)
                      for left, right, label in junctions),
                    *contact_duplex[2:]]
    final = module.final_library()
    final_rows = duplex_rows(final, label="sequencing library")
    steps = "".join(f'<h3>({i}) {title}</h3><p>{text}</p>'
                    for i, (title, text) in enumerate(module.STEPS, 1))
    primers = getattr(module, "SEQ_PRIMERS", ())
    required = getattr(module, "REQUIRED_ROLES", ())
    if primers:
        sequencing = sp.section(final, primers, heading="Sequencing",
            intro=module.READOUT, required_roles=required,
            read_lengths=getattr(module, "READ_LENGTHS", None))
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
