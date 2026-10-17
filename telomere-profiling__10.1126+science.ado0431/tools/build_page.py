#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telomere_profiling as M
from telomere_page import render_page
OUT=HERE.parent/"telomere-profiling.html"
protocol_graph=M.protocol_graph

CAVEAT=("Provenance: reaction conditions and oligos follow the published Science Supplementary Materials "
        "(Materials and Methods pp. 3–4; Table S3, fourteen TeloTag top strands with /5Phos/ and /3BioTEG/; "
        "fig. S1A, NB50 with its six phase splints). Only NB50 splints are printed; the model derives the splint "
        "families of the other TeloTags from their printed strands so that each duplex pairs, and draws NB50. "
        "Published discrepancies: " + "; ".join(M.PUBLISHED_DISCREPANCIES) + ". "
        "The complete ONT SQK-LSK110 motor-bearing adapter assembly is not disclosed. "
        "Model correction: earlier versions took conditions and splint rows from the superseded bioRxiv v2 preprint "
        "(internal enzymes ClaI/AsiSI/PmeI/BamHI, release by PvuI/PacI/EcoRI, kit SQK-LSK114); the published "
        "Methods specify ClaI, PvuI and SQK-LSK110. Fig. S5 additionally compares AsiSI/PvuI and BamHI/EcoRI "
        "digests to show that fragment length does not affect telomere length; those branches are not drawn.")

def render():
    return render_page(title=M.TITLE,notes="01_telomere-profiling.html",citation_html=M.SOURCE,panels=M.panels(),
        entry="motor-loaded ONT adapter ---> one enriched telomere strand through nanopore; no sequencing primer",
        caveat=CAVEAT,graph=protocol_graph(),oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
