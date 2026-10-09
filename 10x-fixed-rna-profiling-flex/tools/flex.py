"""10x Chromium Fixed RNA Profiling / Flex probe-pair chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import adjacent_probe_scene

TITLE="10x Fixed RNA Profiling (Flex)"
NOTES="01_flex.html"
SOURCE='Commercial protocol: 10x Genomics <a href="https://www.10xgenomics.com/support/cn/flex-gene-expression/documentation/steps/library-prep/chromium-fixed-rna-profiling-reagent-kits-for-singleplexed-samples-with-feature-barcode-technology-for-protein-using-barcode-oligo-capture">Fixed RNA Profiling user guide CG000674 Rev C</a>.'
SUMMARY="Pairs of gene-specific probes hybridize adjacently on formaldehyde-fixed RNA and are ligated. The ligated probe product—not the RNA itself—is then given a GEM cell barcode and UMI and converted into a sequencing library."
CAVEAT="Probe-pair and gel-bead sequences are commercial. Their published adjacency, ligation and barcode roles are drawn with inferred placeholders."

LIB, PRIMERS = truseq_library([
    seg("cell barcode", "B"*16, "cbc", placeholder=True, inferred=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="10x Flex whitelist")),
    seg("UMI", "U"*12, "umi", placeholder=True, inferred=True,
        feature=feature("umi", "umi", "random")),
    seg("ligated probe-pair identifier", "X"*30, placeholder=True, inferred=True)],
    "Flex gene-expression library", inferred_adapters=True)
FINAL_LIBRARIES=(("Final probe-derived gene-expression library", LIB, PRIMERS,
                  "INFERRED — Read 1 contains cell barcode and UMI; Read 2 reports the ligated probe-pair identifier, with undisclosed bases visibly inferred.",
                  "The declared Illumina primer sites are shown on the completed probe-derived library."),)

def sections():
    return [
        ("Hybridize adjacent probe halves to fixed RNA", adjacent_probe_scene().rows(),
         "A valid target is represented by two gene-specific probe halves bound next to one another."),
        ("Ligate the adjacent probe pair",
         [Row(chunks=[("[left probe barcode] ** [right probe barcode]", "cbc", False)]),
          Row(chunks=[("                 ** ligation", "me", False)])],
         "Only adjacent, correctly hybridized probe halves form the amplifiable reporter molecule."),
        ("Partition, barcode and amplify",
         [Row(chunks=[("ligated reporter + GEM bead → [cell barcode][UMI][probe-pair identifier]", "cbc", False)]),
          Row(chunks=[("→ pre-amplification → sample-index PCR", None, False)])],
         "INFERRED — The cell barcode and UMI are added after probe ligation; exact commercial transfer-oligo bases are placeholders."),
    ]
