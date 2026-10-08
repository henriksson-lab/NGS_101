"""10x Chromium Next GEM Single Cell Multiome ATAC + Gene Expression."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Construct, Row
import illumina as il
import nextera as nx
import seqprimers as sp
from multimodal_spatial import modality_split

TITLE = "10x Chromium Single Cell Multiome ATAC + Gene Expression"
NOTES = "01_multiome.html"
SOURCE = 'Commercial protocol: 10x Genomics <a href="https://www.10xgenomics.com/support/instruments/chromium-x-series/chromium-next-gem-single-cell-multiome-atac-plus-gene-expression-reagent-kits-user-guide">CG000338 Rev G</a>.'
SUMMARY = "Nuclei are tagmented before GEM formation. Inside each GEM, one gel-bead barcode is transferred to both accessible-DNA fragments and poly(A)-derived cDNA, yielding linked ATAC and gene-expression libraries from the same nucleus."
CAVEAT = "10x does not disclose the complete gel-bead oligos. Published barcode, UMI and adapter roles are shown as length-preserving placeholders; no proprietary bases are inferred."

GEX, GEX_PRIMERS = truseq_library([
    seg("cell barcode", "B" * 16, "cbc", placeholder=True, inferred=True),
    seg("UMI", "U" * 12, "umi", placeholder=True, inferred=True),
    seg("poly(dT) junction", "T" * 12, placeholder=True, inferred=True),
    seg("cDNA", "X" * 36, placeholder=True)], "Multiome gene-expression library",
    inferred_adapters=True)
ATAC = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True, inferred=True),
    seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
    seg("accessible genomic DNA", "X" * 34, placeholder=True),
    seg("right mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("i7 sample index", "I" * 8, "cbc", placeholder=True),
    seg("P7 reverse complement", il.P7_RC, "p7")], name="Multiome ATAC library")
ATAC_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
_atac_problems = sp.verify(ATAC, ATAC_PRIMERS,
                           required_roles=tuple(p.role for p in ATAC_PRIMERS))
if _atac_problems:
    raise ValueError("Multiome ATAC library: " + "; ".join(_atac_problems))

FINAL_LIBRARIES = (
    ("Final gene-expression library", GEX, GEX_PRIMERS,
     "INFERRED — Read 1 reports the GEM cell barcode and UMI; undisclosed adapter and bead-oligo bases remain visibly inferred.",
     "The declared Illumina primers are located on the final gene-expression duplex."),
    ("Final ATAC library", ATAC, ATAC_PRIMERS,
     "INFERRED — Paired genomic reads flank the accessible fragment; Index 2 reads the undisclosed 16-base GEM barcode and Index 1 reads the sample index.",
     "The declared Nextera sequencing primers are located on the final ATAC duplex."),
)

def sections():
    return [
        ("Tagment accessible chromatin before partitioning",
         [Row(chunks=[("nucleus chromatin → Tn5 cleavage + adapter transfer → tagged open-DNA fragments", "me", False)])],
         "Transposition occurs in bulk nuclei; intact nuclei retain both tagged chromatin and RNA."),
        ("Transfer one GEM barcode to RNA and ATAC molecules",
         [Row(chunks=[("poly(A) RNA + barcoded oligo-dT → cell barcode + UMI + cDNA", "cbc", False)]),
          Row(chunks=[("tagged genomic DNA + barcoded bridge → cell barcode + ATAC fragment", "cbc", False)])],
         "INFERRED — Both reactions occur inside the same GEM and therefore use the same published barcode role; undisclosed bases are placeholders."),
        ("Separate and amplify the two modalities", modality_split("gene-expression cDNA library", "ATAC DNA library"),
         "After GEM cleanup, modality-specific amplification preserves the common cell identity while producing independently sequenceable libraries."),
    ]
