#!/usr/bin/env python3
"""Build the diagram-first VASA-seq chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import vasa as V
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "vasa-seq.html"


def preamble() -> str:
    return f'''<div class="wrap">
<h1>VASA-seq &mdash; total RNA by fragmenting and poly(A)-tailing</h1>
{info('Defining source: <a href="https://doi.org/10.1038/s41587-022-01361-8">Salmen et al., Nature Biotechnology (2022)</a>.')}
<div class="caveat"><b>Inferred regions.</b> The paper cites earlier capture-oligo and VASA-plate PCR designs without printing them, and its VASA-drop adapter annotation does not establish a complete internal sequence. Bracketed dotted regions preserve those boundaries without importing secondary reconstructions.</div>
'''


def oligos() -> str:
    i5 = V.drop_i5_primer("N" * V.INDEX_LEN)
    i7 = V.drop_i7_primer("N" * V.INDEX_LEN)
    rows = [
        oligo("RA3", [V.seg("RA3", V.RA3, "r2")], mods="5rApp / 3SpC3"),
        oligo("RTP", V.rtp_segments()),
        oligo("VD_ILMN8_i5", list(i5)),
        oligo("VD_ILMN8_i7", list(i7)),
    ]
    return "<h2>Oligos printed by the paper</h2><seq>" + "\n".join(rows) + "</seq>"


def capture_scene(fmt: str) -> Scene:
    target = V.repaired_fragment()
    primer = V.capture_primer(fmt)
    sc = Scene()
    sc.strand("RNA", list(target), label="RNA")
    sc.anneal(f"{fmt} capture oligo", list(primer), to="RNA",
              pair=("oligo-dT", "added poly(A)"),
              unpaired=(f"{fmt} barcode / T7 region", "UFI", "cell barcode"))
    sc.arrow(f"{fmt} capture oligo", "reverse transcriptase")
    return sc


def rtp_scene() -> Scene:
    arna = V.ligated_arna()
    sc = Scene()
    sc.strand("aRNA + RA3", list(arna), label="aRNA")
    sc.anneal("RTP", V.rtp_segments(), to="aRNA + RA3", pair=("RTP core", "RA3"),
              unpaired=("RTP 5' G",))
    sc.junction("aRNA + RA3", "rRNA-depleted aRNA", "RA3")
    sc.arrow("RTP", "reverse transcriptase")
    return sc


def final_panel(con, caption: str) -> str:
    sc = V.schematic_duplex(con)
    return panel([*sc.rows(), *annotation_rows(con)], cls="small", caption=caption)


def steps() -> str:
    repaired = V.repaired_fragment()
    plate_ivt, drop_ivt = V.ivt_template("plate"), V.ivt_template("drop")
    return f'''<h2>Library generation</h2>
<h3>(1) Heat-fragment total RNA, repair the ends, and add poly(A)</h3>
{panel([strand_row(repaired), *annotation_rows(repaired)], cls="small", caption="Every repaired RNA fragment receives a poly(A) tail; its exact length is not specified.")}

<h3>(2) Barcode and reverse-transcribe each poly(A)-tailed fragment</h3>
{panel(capture_scene("plate").rows(), cls="small", caption="INFERRED — VASA-plate uses a cited CEL-seq2/SORT-seq capture oligo whose sequence is not printed in the VASA-seq paper.")}
{panel(capture_scene("drop").rows(), cls="small", caption="INFERRED — VASA-drop uses a cited inDrop v3 bead oligo whose sequence is not printed in the VASA-seq paper.")}

<h3>(3) Make double-stranded templates and amplify by T7 IVT</h3>
{final_panel(plate_ivt, "INFERRED — VASA-plate IVT template; the bracketed capture/T7 region is retained only as a role boundary.")}
{final_panel(drop_ivt, "INFERRED — VASA-drop IVT template; the bracketed capture/T7 region is retained only as a role boundary.")}

<h3>(4) Deplete rRNA from amplified RNA</h3>
{info('DNA probes hybridize to rRNA-derived aRNA; thermostable RNase H removes those molecules, then DNase removes the probes.')}

<h3>(5) Ligate a 3' adapter to aRNA, then reverse-transcribe</h3>
{panel(rtp_scene().rows(), cls="small", caption="The supplement prints this RA3/RTP pair, which anneals exactly. Its assignment to VASA-drop does not interlock with the printed drop PCR arm, so the drop adapter junction remains unresolved below.")}

<h3>(6a) VASA-plate: amplify with the cited plate PCR system</h3>
{final_panel(V.plate_library(), "INFERRED — final VASA-plate library. The paper establishes the UFI/barcode read boundary and RA3, but does not print the terminal PCR arms.")}

<h3>(6b) VASA-drop: amplify with one VD_ILMN8_i5 / VD_ILMN8_i7 pair</h3>
{final_panel(V.drop_library(), "INFERRED — final VASA-drop library. PCR-derived ends and paper-reported read boundaries are shown; the unpublished bead-linker/ligation geometry remains bracketed.")}
'''


def sequencing() -> str:
    return f'''<h2>Read layout</h2>
{sp.unavailable_diagram(V.plate_library(), "the VASA-plate library-primer sequences are not printed")}
{sp.unavailable_diagram(V.drop_library(), "the complete VASA-drop primer landing sites are not published")}
{table(("Format", "Read", "Content"), V.read_layouts())}
</div>'''


def main() -> None:
    OUT.write_text("\n".join([head("VASA-seq library chemistry"), preamble(), oligos(),
                              steps(), sequencing()]), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
