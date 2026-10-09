"""Original Liu et al. DBiT-seq architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import nextera_library, seg, spatial_rt_scene
from chemdraw import Row, feature

TITLE = "DBiT-seq — deterministic barcoding in tissue"
NOTES = "01_dbit-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1016/j.cell.2020.10.026">Liu et al., <i>Cell</i> (2020)</a>.'
SUMMARY = "Two perpendicular 50-channel microfluidic flows install an 8-nt A barcode during reverse transcription and ligate an 8-nt B barcode plus UMI at each channel intersection, creating up to 2,500 addressable tissue pixels."
CAVEAT = "The paper publishes the barcode and linker tables in supplements. This page models their verified roles and lengths; the individual whitelist bases are represented by placeholders."
BARCODE_A=feature("spatial_barcode_a","spatial_barcode","combinatorial",group="spatial_barcode",part="flow A")
BARCODE_B=feature("spatial_barcode_b","spatial_barcode","combinatorial",group="spatial_barcode",part="flow B")
UMI_FEATURE=feature("umi","umi","random")

FINAL_LIBRARY, SEQ_PRIMERS = nextera_library([
    seg("cDNA / ADT insert", "X" * 34, placeholder=True),
    seg("poly(dT) junction", "T" * 10, placeholder=True),
    seg("barcode A", "A" * 8, "cbc", placeholder=True, feature=BARCODE_A),
    seg("ligation linker", "L" * 15, "me", placeholder=True),
    seg("barcode B", "B" * 8, "cbc", placeholder=True, feature=BARCODE_B),
    seg("UMI", "U" * 10, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("PCR handle", "H" * 22, placeholder=True)], "DBiT-seq library")
FINAL_CAPTION = "After tissue recovery, template switching and cDNA amplification, Nextera XT tagmentation produces the paired-end library. The internal A–B junction retains the row/column address."
SEQUENCING_INTRO = "The original 2×100 run reads the A/B spatial barcode and UMI from one end and transcript or antibody-tag sequence from the other."

def sections():
    return [
        ("Flow A stripes: reverse-transcribe with the first coordinate",
         spatial_rt_scene(barcode_parts=(("barcode A (8 nt)",8),), umi=0,
                          surface="15-nt ligation linker").rows(),
         "Each of 50 parallel channels supplies one A barcode attached to oligo-dT; in situ reverse transcription writes that stripe identity into cDNA."),
        ("Flow B stripes: ligate the perpendicular coordinate",
         [Row(chunks=[("A1  A2  A3  …  A50", "cbc", False)]),
          Row(chunks=[(" ╳   ╳   ╳       ╳    ← B1", "me", False)]),
          Row(chunks=[(" ╳   ╳   ╳       ╳    ← B2", "me", False)]),
          Row(chunks=[("each intersection: [A 8 nt] ** [B 8 nt][UMI][PCR handle]", "cbc", False)])],
         "A second chip is rotated 90°. T4 ligase joins the two 15-nt linker system at each intersection; ** marks the in-tissue ligation boundary."),
        ("Recover, template-switch, amplify and tagment",
         [Row(chunks=[("spatial cDNA → streptavidin recovery → template switch → PCR", None, False)]),
          Row(chunks=[("→ Nextera XT tagmentation → paired-end library", "me", False)])],
         "Biotin on barcode B enables recovery. The spatial address is installed before bulk extraction and survives the downstream library conversion."),
    ]
