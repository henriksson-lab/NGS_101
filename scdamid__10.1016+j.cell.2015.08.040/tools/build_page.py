#!/usr/bin/env python3
"""Build the molecule-focused scDamID schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import scdamid as D
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel, protocol_graph_panels
from adapter_ligation import AdapterLigationGeometry
from endprep import dA_tailed_scene, repair_and_dA_tail
from page import head, info

OUT = HERE.parent / "scdamid.html"


def oligos() -> str:
    preamp = [D.predecessor("N4", "N" * D.RANDOM_NT, placeholder=True),
              D.predecessor("DamID stub", D.PREAMP_STUB),
              D.predecessor("restored junction", D.JUNCTION)]
    pairs = len(D.ADRB)
    return ("<h2>Key oligos</h2>"
            "<p><info>INFERRED — the paper's own oligo figure (S1E) is unavailable. These "
            "are the conventional-DamID oligos it says it reuses (AdRt, AdRb and AdR_PCR, "
            "published with Vogel et al. 2007), plus the NNNN 5' tail the figure legend "
            "does state.</info></p><seq>"
            + oligo("DamID AdRt", [D.predecessor("", D.ADRT)])
            + oligo("DamID AdRb", [D.predecessor("", D.ADRB)])
            + oligo("single pre-amplification primer", preamp)
            + "</seq>"
            + panel(D.damid_adaptor_scene().rows(), cls="small",
                    caption="INFERRED — AdRb pairs only the "
                            f"{pairs}-nt ligating end of AdRt, leaving a long 5' tail."))


def construction() -> str:
    final = D.final_library()
    dpni = Scene.duplex(D.dpni_fragment().segments)
    final_sc = Scene.duplex(final.segments)
    final_sc.strands["top"].label = final_sc.strands["bottom"].label = ""
    return f"""<h2>Library construction</h2>
<h3>(1) DpnI cuts Dam-methylated GATC sites</h3>
{panel(dpni.rows(), cls="small",
       caption="The selected blunt fragment lies between two Dam-methylated GATC sites.")}
<h3>(2) Ligate the DamID adaptor to both blunt ends</h3>
{panel(D.damid_adaptor_scene().rows(), cls="small",
       caption="INFERRED — the defining paper establishes adaptor ligation after DpnI but "
               "prints no adaptor; the bases shown are the reused conventional-DamID pair.")}
<h3>(3) Amplify with one N4-tailed primer</h3>
{panel(D.preamp_duplex().rows(), cls="small",
       caption="INFERRED — the same restored adaptor/GATC junction at both ends lets one primer "
               "amplify the selected fragment.")}
<h3>(4) End repair, dA-tail, ligate the indexed Illumina adaptor and PCR</h3>
{panel(dA_tailed_scene(repair_and_dA_tail(D.preamp_product(), inferred=True)).rows(), cls="small",
       caption="INFERRED — reconstructed end preparation produces the dA-tailed substrate before adapter ligation.")}
{panel(AdapterLigationGeometry.unavailable("The accessible defining paper establishes indexed Illumina adapter ligation but does not print that adapter's terminal bases, and names no kit.").rows(), cls="small",
       caption="The authoritative source does not disclose the indexed adapter overhang.")}
{panel([*final_sc.rows(), *annotation_rows(final)], cls="long",
       caption="INFERRED — every dotted region is unpublished for this paper: the outer "
               "Illumina arms are canonical constants and the index width a placeholder, "
               "because Table S2 is unavailable and no kit is named.")}"""


def sequencing() -> str:
    return ("<h2>Sequencing</h2>"
            + info("Single-end 51-nt Read 1 (stated by the Figure S1E legend) plus an i7 "
                   "index whose width is unpublished and drawn as a placeholder. Read 1 "
                   "begins with N4, the DamID stub and the restored GATC before entering "
                   "genomic DNA.")
            + sp.section(D.final_library(), D.SEQ_PRIMERS,
                         intro="Only the two reads actually used by this single-end, "
                               "single-index run are declared.",
                         required_roles=D.RUN_ROLES))


def render() -> str:
    graph = D.protocol_graph()
    return "\n".join([
        head("scDamID library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>scDamID &mdash; single-cell nuclear-lamina contact mapping</h1>',
        info('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.08.040">'
             'Kind et al., <i>Cell</i> (2015)</a>.'),
        '<div class="caveat"><b>Sequence boundary.</b> The accessible paper establishes DpnI '
        'selection, adaptor ligation, single-primer amplification and indexed Illumina '
        'sequencing, but its own oligo figure (S1E) and indexing-primer table (S2) are behind '
        'the PMC download gate. Every dotted base and every step labelled <b>INFERRED</b> is '
        'therefore unestablished for this paper: the DamID adaptor and pre-amplification primer '
        'are the conventional-DamID oligos it says it reuses (published with Vogel et al. '
        '2007), and the outer Illumina arms and index width are canonical parts and '
        'placeholders.</div>',
        oligos(), '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        sequencing(), '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
