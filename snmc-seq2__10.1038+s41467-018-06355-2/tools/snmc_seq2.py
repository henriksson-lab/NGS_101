"""Luo et al. 2018 snmC-seq2 chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from base_conversion import bisulfite_path
from batch_ngs import seg
from chemdraw import Row, Scene, Segment, feature, oligo
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
import snmc2_core

TITLE = "snmC-seq2 — indexed single-nucleus methylomes"
NOTES = "01_snmc-seq2.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41467-018-06355-2">Luo et al., <i>Nature Communications</i> (2018)</a>.'
SUMMARY = "Bisulfite-converted nuclear DNA is copied with one of eight indexed random primers, completed by Adaptase, and amplified with unique dual indexes."
CAVEAT = ("The defining supplement prints the random-primer and 768 dual-index primer "
          "pairs. That sequence is exact below. The bases added by the proprietary "
          "Adaptase reaction are not printed, so only that short region remains inferred.")
INLINE = snmc2_core.INLINE
P5_INDEX_1 = snmc2_core.P5_INDEX_1
P7_INDEX_1 = snmc2_core.P7_INDEX_1
P5_INDEX_1_CODE = snmc2_core.P5_INDEX_1_CODE
P7_INDEX_1_CODE = snmc2_core.P7_INDEX_1_CODE


def library():
    return snmc2_core.library("bisulfite-converted genomic insert", insert_nt=36,
                              library_name="snmC-seq2 library")


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = "INFERRED — dual-index endpoint with the source-printed 6-nt inline barcode and H9; the Adaptase-added bases are unavailable."
SEQUENCING_INTRO = "Read 1 begins with the inline cell barcode and random-primer sequence; i5 and i7 form the unique dual index."


def sections():
    return [
        ("Bisulfite-convert each nucleus", [
            Row(chunks=[(bisulfite_path(protected=False).text(), None, False)]),
            Row(chunks=[(bisulfite_path(protected=True).text(), "w1", False)]),
        ], "Unmethylated C is copied as T; 5mC remains C."),
        ("Indexed random priming", [
            Row(chunks=[("5′ spacer — Read 1 tail — 6-nt cell barcode — H9 → extension", "cbc", False)])
        ], "Eight source-printed random primers barcode plate quadrants before pooling."),
        ("INFERRED — Adaptase completion; exact indexed PCR", [
            Row(chunks=[("random-primer product → Adaptase tail → P5/i5 + P7/i7 PCR", None, False)])
        ], "The representative PCR pair is source-exact; only the proprietary Adaptase-added bases are unknown."),
    ]


def oligos():
    yield oligo("P5L_AD002_H random primer", [
        Segment("Read 1 partial handle", "TTCCCTACACGACGCTCTTCCGATCT"),
        Segment("inline cell barcode", INLINE, "cbc",
                feature=feature("cell_inline", "cell_barcode", "whitelist",
                                whitelist="published random-primer barcode set")),
        Segment("random H9 (A/C/T only)", "H" * 9, placeholder=True),
    ], mods="/5SpC3/ before the printed sequence")
    yield oligo("dual-index P5 primer 1", [Segment("P5 primer", P5_INDEX_1)])
    yield oligo("dual-index P7 primer 1", [Segment("P7 primer", P7_INDEX_1)])


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    nucleus=g.input("nucleus","one isolated nucleus containing genomic DNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
    bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
    converted=g.product("converted","bisulfite-converted nuclear DNA",MaterialKind.MOLECULAR,state=sections()[0][1])
    g.react("bisulfite_convert","Bisulfite-convert the nuclear DNA",StepKind.CHEMICAL_CONVERSION,((nucleus,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
    random_primer=g.input("random_primer","source-printed indexed H9 random primer",MaterialKind.OLIGO)
    polymerase=g.input("strand_polymerase","random-primer extension polymerase",MaterialKind.CATALYST)
    primed=g.product("primed","inline-cell-barcoded random-primer products",MaterialKind.MOLECULAR,state=sections()[1][1])
    g.react("random_prime","Copy converted DNA with an indexed H9 random primer",StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(primed,))
    adaptase=g.input("adaptase","proprietary Adaptase adapter-completion reagent",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (added bases not disclosed)")
    adaptase_enzyme=g.input("adaptase_enzyme","Adaptase completion enzyme",MaterialKind.CATALYST)
    adapted=g.product("adapted","Adaptase-completed dual-index PCR template",MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (Adaptase-added bases not disclosed)")
    g.react("adaptase_complete","Add the proprietary Adaptase completion end",StepKind.ADAPTER_ATTACHMENT,((primed,InputRole.TEMPLATE),(adaptase,InputRole.STRUCTURAL),(adaptase_enzyme,InputRole.CATALYST)),(adapted,))
    reviewed=["nucleus","bisulfite","bisulfite_convert","random_primer","strand_polymerase","random_prime","adaptase","adaptase_enzyme","adaptase_complete"]
    pcr_geometry = snmc2_core.representative_pcr_geometry(FINAL_LIBRARY)
    finish_pcr_branch(g,"final",adapted,FINAL_LIBRARY,SEQ_PRIMERS,reviewed,
                      label="snmC-seq2 methylome", geometry=pcr_geometry)
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
