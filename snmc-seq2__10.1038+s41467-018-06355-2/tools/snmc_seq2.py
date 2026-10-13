"""Luo et al. 2018 snmC-seq2 chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from base_conversion import bisulfite_path
from batch_ngs import seg, truseq_library
from chemdraw import Row, Scene, feature
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

TITLE = "snmC-seq2 — indexed single-nucleus methylomes"
NOTES = "01_snmc-seq2.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41467-018-06355-2">Luo et al., <i>Nature Communications</i> (2018)</a>.'
SUMMARY = "Bisulfite-converted nuclear DNA is copied with one of eight indexed random primers, completed by Adaptase, and amplified with unique dual indexes."
CAVEAT = ("The defining supplement prints the random-primer and dual-index families, but not the bases added by the proprietary Adaptase reaction. "
          "That short region is explicitly inferred below. The printed terminal PCR "
          "families have not yet been transcribed into this compact module, so their "
          "landing geometry remains unavailable.")
INLINE = "CGATGT"  # representative P5L_AD002 random-primer barcode


def library():
    inserts = [seg("6-nt inline cell barcode", INLINE, "cbc",
                   feature=feature("cell_inline", "cell_barcode", "whitelist", whitelist="published random-primer barcode set")),
               seg("random N9", "N" * 9, placeholder=True),
               seg("bisulfite-converted genomic insert", "X" * 36, placeholder=True),
               seg("Adaptase-added bases", "[Adaptase bases]", placeholder=True, inferred=True,
                   note="proprietary Adaptase sequence is not printed")]
    return truseq_library(inserts, "snmC-seq2 library", dual_index=True)


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = "INFERRED — dual-index endpoint with the source-printed 6-nt inline barcode and N9; the Adaptase-added bases are unavailable."
SEQUENCING_INTRO = "Read 1 begins with the inline cell barcode and random-primer sequence; i5 and i7 form the unique dual index."


def sections():
    return [
        ("Bisulfite-convert each nucleus", [
            Row(chunks=[(bisulfite_path(protected=False).text(), None, False)]),
            Row(chunks=[(bisulfite_path(protected=True).text(), "w1", False)]),
        ], "Unmethylated C is copied as T; 5mC remains C."),
        ("Indexed random priming", [
            Row(chunks=[("5′ spacer — Read 1 tail — 6-nt cell barcode — N9 → extension", "cbc", False)])
        ], "Eight source-printed random primers barcode plate quadrants before pooling."),
        ("INFERRED — Adaptase completion and indexed PCR", [
            Row(chunks=[("random-primer product → Adaptase tail → P5/i5 + P7/i7 PCR", None, False)])
        ], "The PCR primer families are printed; only the proprietary Adaptase-added bases are unknown."),
    ]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    nucleus=g.input("nucleus","one isolated nucleus containing genomic DNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
    bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
    converted=g.product("converted","bisulfite-converted nuclear DNA",MaterialKind.MOLECULAR,state=sections()[0][1])
    g.react("bisulfite_convert","Bisulfite-convert the nuclear DNA",StepKind.CHEMICAL_CONVERSION,((nucleus,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
    random_primer=g.input("random_primer","source-printed indexed N9 random primer",MaterialKind.OLIGO)
    polymerase=g.input("strand_polymerase","random-primer extension polymerase",MaterialKind.CATALYST)
    primed=g.product("primed","inline-cell-barcoded random-primer products",MaterialKind.MOLECULAR,state=sections()[1][1])
    g.react("random_prime","Copy converted DNA with an indexed N9 random primer",StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(primed,))
    adaptase=g.input("adaptase","proprietary Adaptase adapter-completion reagent",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (added bases not disclosed)")
    adaptase_enzyme=g.input("adaptase_enzyme","Adaptase completion enzyme",MaterialKind.CATALYST)
    adapted=g.product("adapted","Adaptase-completed dual-index PCR template",MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (Adaptase-added bases not disclosed)")
    g.react("adaptase_complete","Add the proprietary Adaptase completion end",StepKind.ADAPTER_ATTACHMENT,((primed,InputRole.TEMPLATE),(adaptase,InputRole.STRUCTURAL),(adaptase_enzyme,InputRole.CATALYST)),(adapted,))
    reviewed=["nucleus","bisulfite","bisulfite_convert","random_primer","strand_polymerase","random_prime","adaptase","adaptase_enzyme","adaptase_complete"]
    finish_pcr_branch(g,"final",adapted,FINAL_LIBRARY,SEQ_PRIMERS,reviewed,label="snmC-seq2 methylome",inferred_primers=True,limitation="The source prints the dual-index PCR families, but their exact sequences have not yet been transcribed into this compact module.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
