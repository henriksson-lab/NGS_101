"""Peterson et al. ddRAD-seq library chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Row, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from restriction import ECORI, MSPI

TITLE = "ddRAD-seq — double-digest RAD sequencing"
NOTES = "01_ddrad-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1371/journal.pone.0037135">'
          'Peterson et al., <i>PLOS ONE</i> (2012)</a>.')
SUMMARY = ("Digest genomic DNA with two enzymes, ligate end-selective adapters carrying "
           "an inline sample barcode, pool and size-select, then add the Illumina outer "
           "adapter and a second sample index by PCR.")
CAVEAT = ("The displayed GCATG inline barcode and ATCACG i7 are examples from published "
          "Sequences S1. This graph follows the route without optional biotin-P2 capture; "
          "same-end products are not physically removed at digestion. Protocol S1 says "
          "P1/P1 products can amplify but do not cluster. Their downstream recovery is unresolved here. "
          "The primer sites show compatibility with single- or paired-end sequencing, "
          "with or without an index read; no single acquired cycle recipe is assigned.")

INLINE = "GCATG"
I7 = "ATCACG"
P1_TOP = il.TRUSEQ_READ1 + INLINE
P1_BOTTOM = ECORI.overhang + "CATGC" + il.INDEX2_PRIMER_RC
P2_TOP = il.TRUSEQ_READ2
P2_BOTTOM = MSPI.overhang + "AGATCGGAAGAGCGAGAACAA"

# A source transcription is accepted only through the enzyme-derived end model.
P1_END = ECORI.adapter_end(P1_BOTTOM[:len(ECORI.overhang)])
P2_END = MSPI.adapter_end(P2_BOTTOM[:len(MSPI.overhang)])


def oligos():
    return [f"P1.1: 5′-{P1_TOP}-3′",
            f"P1.2: 5′-/5Phos/{P1_BOTTOM}-3′",
            f"P2.1: 5′-{P2_TOP}-3′",
            f"P2.2: 5′-/5Phos/{P2_BOTTOM}-3′"]


def adapter_geometry() -> AdapterLigationGeometry:
    p1=Scene(); p1.strand("P1",[
        Segment("P1 EcoRI overhang",P1_END.written_overhang),
        Segment("P1 adapter body",P1_BOTTOM[len(P1_END.written_overhang):])],
        label="P1.2 adapter strand", mod5="p")
    p1.ligation_end("P1","P1 EcoRI overhang",kind="cohesive",polarity="5-prime",
                    owner="adapter",text="EcoRI-compatible 5-prime overhang")
    p2=Scene(); p2.strand("P2",[
        Segment("P2 MspI overhang",P2_END.written_overhang),
        Segment("P2 adapter body",P2_BOTTOM[len(P2_END.written_overhang):])],
        label="P2.2 adapter strand", mod5="p")
    p2.ligation_end("P2","P2 MspI overhang",kind="cohesive",polarity="5-prime",
                    owner="adapter",text="MspI-compatible 5-prime overhang")
    return AdapterLigationGeometry.known(
        "cohesive",p1,p2,
        statement="The published P1 and P2 adapter overhangs are validated against the EcoRI and MspI cut-derived cohesive ends.")


def final_library() -> tuple[Construct, tuple]:
    lib = Construct([
        Segment("P5 before shared ACAC", il.P5[:-4], "p5"),
        Segment("Read 1 site (shares ACAC with P5)", il.TRUSEQ_READ1, "r1"),
        Segment("inline sample barcode", INLINE, "cbc",
                feature=feature("inline_sample_barcode", "sample_index", "whitelist",
                                whitelist="published P1 adapter barcode set")),
        Segment("EcoRI remnant", ECORI.overhang),
        Segment("size-selected genomic insert", "X" * 34, placeholder=True),
        Segment("MspI remnant", MSPI.overhang),
        Segment("Read 2 site start", "A", "r2"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("6-nt i7", I7, "cbc",
                feature=feature("sample_i7", "sample_index", "whitelist",
                                whitelist="published indexed PCR primer set")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="ddRAD-seq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid ddRAD-seq endpoint: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = final_library()
FINAL_CAPTION = ("A P1–insert–P2 molecule after PCR. Read 1 encounters the inline barcode "
                 "before the EcoRI-associated genomic end; Index 1 reports i7, and Read 2 "
                 "enters from the MspI-associated end.")
SEQUENCING_INTRO = ("Standard TruSeq Read 1, Index 1 and Read 2 primers bind. Demultiplexing "
                    "uses the inline P1 barcode together with the optional six-base i7.")


def sections():
    digest = [
        Row(chunks=[(f"EcoRI  G^{ECORI.overhang}C  →  5′-{ECORI.overhang} cohesive end", None, False)]),
        Row(chunks=[(f"MspI   C^{MSPI.overhang}G    →  5′-{MSPI.overhang} cohesive end", None, False)]),
    ]
    lib = FINAL_LIBRARY
    sc = Scene.duplex(list(lib), label="selected ddRAD molecule")
    sc.junction("top", "inline sample barcode", "EcoRI remnant", "P1 ligation")
    sc.junction("top", "MspI remnant", "Read 2 site start", "P2 ligation")
    sc.labels("top")
    return [
        ("Double digest", digest,
         "EcoRI-HF and MspI make distinguishable 5′ cohesive ends; their overhangs are derived from the cut sites."),
        ("Ligate end-selective adapters", sc.rows(),
         "P1 carries the inline barcode and fits the EcoRI end; P2 fits the MspI end. ** marks each ligation."),
        ("Pool, size-select and PCR", [Row(chunks=[
            ("barcoded ligation products → pooled size window → 8–12 PCR cycles → indexed library", None, False)])],
         "Pooling occurs after inline barcoding; a narrow size window defines the reduced genomic representation."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("genomic_dna","genomic DNA",MaterialKind.MOLECULAR,primary=True)
    ecori=graph.input("ecori","EcoRI-HF",MaterialKind.CATALYST)
    mspi=graph.input("mspi","MspI",MaterialKind.CATALYST)
    compatible=graph.product("compatible","double-digest fragments, including unlike- and same-enzyme ends",
                             MaterialKind.MOLECULAR)
    graph.react("double_digest","Digest genomic DNA with EcoRI-HF and MspI",
                StepKind.DIGESTION,
                ((dna,InputRole.TEMPLATE),(ecori,InputRole.CATALYST),
                 (mspi,InputRole.CATALYST)),(compatible,))
    adapters=graph.input("adapters","barcoded P1 and P2 cohesive-end adapters",
                         MaterialKind.ADAPTER,state=tuple(adapter_geometry().rows()))
    ligase=graph.input("ligase","T4 DNA ligase",MaterialKind.CATALYST)
    ligated=graph.product("ligated","P1/P2-ligated reduced-representation fragments",
                          MaterialKind.MOLECULAR)
    other=graph.product("other_ligated","same-end adapter-ligated fragment classes",
                        MaterialKind.MOLECULAR)
    graph.react("ligate","Ligate end-selective P1 and P2 adapters",StepKind.LIGATION,
                ((compatible,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,other),geometry=adapter_geometry())
    graph.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,
        "Protocol S1, optional streptavidin-removal section: P1/P1 products can amplify "
        "but do not cluster. Optional biotin-P2 capture is not selected here; "
        "same-end products may persist through downstream size selection/PCR. "
        "This compact graph does not reconstruct their recovery or amplification."))
    size_rule=graph.input("size_rule","narrow gel size-selection interval",
                          MaterialKind.SELECTION)
    selected=graph.product("selected","pooled size-selected P1/P2 fragments",
                           MaterialKind.MOLECULAR)
    offsize=graph.product("offsize","off-size ligation products",MaterialKind.WASTE)
    graph.react("size_select","Pool and select the intended fragment-size window",
                StepKind.SELECTION,
                ((ligated,InputRole.TEMPLATE),(size_rule,InputRole.SELECTION)),
                (selected,offsize))
    graph.discard(offsize,"off-size fragments are excised away")
    forward=graph.input("pcr_forward","P5-side library PCR primer",MaterialKind.OLIGO)
    reverse=graph.input("pcr_reverse","indexed P7-side library PCR primer",
                        MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","Phusion DNA polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("library_pcr","Complete the outer arms and add the i7 index",
                StepKind.PCR,
                ((selected,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The final primer sites are exact; the compact graph does not reconstruct every strand of the pre-PCR P1/P2 product."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "genomic_dna","ecori","mspi","double_digest","adapters","ligase",
        "ligate","size_rule","size_select","pcr_forward","pcr_reverse",
        "pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
