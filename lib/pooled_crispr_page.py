"""Small shared renderer for the four deliberately separate pooled-screen pages."""
from __future__ import annotations

import html
import pooled_crispr as M
import seqprimers as sp
from chemdraw import duplex_rows, oligo, panel
from page import head, info


def construct_panel(con, caption: str) -> str:
    return panel(duplex_rows(con, label=con.name),
                 cls="small", caption=caption)


def map_panel(version: str) -> str:
    parts = M.vector_map(version)
    text = "  ->  ".join(name for name, _tag in parts)
    return f'<pre class="panel small"><seq>{html.escape(text)}</seq></pre>'


def cloning() -> str:
    top, bottom = M.guide_oligos()
    return f'''<h2>Guide cloning</h2>
<seq>{oligo("top guide oligo", [M.seg("BsmBI end", top[:4], "me"), M.seg("+1 G + guide", top[4:], "cbc", placeholder=True)])}
{oligo("bottom guide oligo", [M.seg("BsmBI end", bottom[:4], "me"), M.seg("guide complement", bottom[4:], "cbc", placeholder=True)])}</seq>
{panel(M.cloning_scene().rows(), cls="small", caption="BsmBI-directed ligation. Asterisks mark the two junctions; the 1,885-bp stuffer is gone.")}'''


def one_pcr_readout(heading: str = "One-PCR screen readout") -> str:
    lib = M.one_pcr_library()
    return f'''<h2>{heading}</h2>
{construct_panel(lib, "One PCR adds the complete P5/P7 library around the integrated guide cassette.")}
{sp.section(lib, M.one_pcr_sequencing_primers(), heading="Sequencing primers", intro="The published run uses 80 Read 1 cycles and an 8-cycle i7 read; both primer sites are located below.", required_roles=("Read 1", "Index 1 (i7)"))}'''


def page(title: str, note: str, source: str, body: str) -> str:
    return "\n".join([head(title), '<div class="wrap">', f"<h1>{title}</h1>",
        f'<p class="research-notes"><a href="{note}">Research notes</a></p>', info(source),
        body, "</div>"])
