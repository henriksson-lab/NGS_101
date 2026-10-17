"""5-prime direct-capture Perturb-seq (standard sgRNA-CR1 branch)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind,
                            Terminal)

TITLE = "5′ direct-capture Perturb-seq — standard sgRNA-CR1 branch"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41587-020-0470-y">Replogle et al., <i>Nature Biotechnology</i> (2020)</a>.'

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-chromium-5prime-v1")
UMI_FEATURE = feature("umi", "umi", "random")
GUIDE_FEATURE = feature("guide", "guide_barcode", "unknown")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE,
            "protospacer": GUIDE_FEATURE, "copied protospacer": GUIDE_FEATURE,
            "i7": I7_FEATURE, "i7 reverse complement": I7_FEATURE}


def seg(name, top, tag=None, **kw):
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


TENX_RT_ADAPTER = "AAGCAGTGGTATCAACGCAGAGTAC"
GUIDE_ANNEAL = "CAAGTTGATAACGGACTAGCC"
OJR160 = TENX_RT_ADAPTER + GUIDE_ANNEAL

OJR163 = il.TRUSEQ_P5_FULL
I7 = "AGGAGTCC"                         # example index printed in the protocol
NESTED_ANNEAL = TENX_RT_ADAPTER[-5:] + GUIDE_ANNEAL
OJR165 = il.P7 + I7 + nx.ADAPTOR_S7 + NESTED_ANNEAL

# CR1, the parental constant region "as in pU6-sgRNA EF1Alpha-puro-T2A-BFP"
# (Supplementary Note 2).  Published Supplementary Table 1 prints the full
# guide scaffold of pBA900 (CR1 with cs2 appended at the 3' end, before the
# poly-T tract); CR1 is that printed scaffold upstream of cs2.  It carries A
# (not the T of the Chen 2013 F+E scaffold) at position 13 ("GAGCTAaGCTGG").
CS2 = "GCTCACCTATTAGCGGCTAAGG"
PBA900_SCAFFOLD = ("GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTT"
                   "GAAAAAGTGGCACCGAGTCGGTGCGCTCACCTATTAGCGGCTAAGGTTTTTTT")
SCAFFOLD = PBA900_SCAFFOLD[:PBA900_SCAFFOLD.index(CS2)]
SCAFFOLD_SITE = revcomp(GUIDE_ANNEAL)
SITE_OFFSET = SCAFFOLD.index(SCAFFOLD_SITE)

# The inspected sources state only that the bead carries the TruSeq Read 1
# site ("P5-TruSeq Read 1 / Index 2 (on bead sequence)", Supplementary Note 3)
# and that RT template-switches onto it.  The partial-Read-1 extent, switch
# spacer and rGrGrG come from the 10x 5' kit design (CG000086, not inspected)
# and are drawn as inferred.
PARTIAL_R1 = il.TRUSEQ_READ1[11:]
SWITCH_SPACER = "TTTCTTATAT"
BEAD_CAVEAT = ("INFERRED — the barcoded 10x 5′ bead oligo (partial Read 1 extent, "
               "16-nt cell barcode, 10-nt UMI, switch spacer and rGrGrG) is not printed in "
               "the paper or Supplementary Note 3; its bases follow the 10x 5′ kit design "
               "(User Guide CG000086, not inspected here). Only the TruSeq Read 1 handle on "
               "the bead and the 26-cycle Read 1 are stated by the source.")


def ojr160_segments():
    return [seg("10x RT adapter", TENX_RT_ADAPTER, "r1"),
            seg("guide annealing", GUIDE_ANNEAL, "tso")]


def ojr163_segments():
    return [seg("P5", il.P5, "p5"), seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1")]


def ojr165_segments():
    return [seg("P7", il.P7, "p7"), seg("i7", I7, "cbc"),
            seg("Read 2 arm", nx.ADAPTOR_S7, "r2"),
            seg("nested guide site", NESTED_ANNEAL, "tso")]


def bead_tso():
    return [
        seg("Partial Read 1", PARTIAL_R1, "r1", inferred=True),
        seg("cell barcode", "B" * 16, "cbc", placeholder=True, inferred=True),
        seg("UMI", "U" * 10, "umi", placeholder=True, inferred=True),
        seg("switch spacer", SWITCH_SPACER, "tso", inferred=True),
        seg("rGrGrG", "GGG", "tso", inferred=True),
    ]


def guide_rna():
    return [
        seg("protospacer", "N" * 20, "insert", placeholder=True),
        seg("CR1 before capture site", SCAFFOLD[:SITE_OFFSET], "insert"),
        seg("oJR160 site", SCAFFOLD_SITE, "tso"),
        seg("CR1 remainder", SCAFFOLD[SITE_OFFSET + len(SCAFFOLD_SITE):], "insert"),
    ]


def capture_scene():
    sc = Scene()
    sc.strand("guide", guide_rna(), label="sgRNA-CR1")
    sc.anneal(
        "oJR160",
        ojr160_segments(),
        to="guide", pair=("guide annealing", "oJR160 site"),
        label="oJR160", unpaired=("10x RT adapter",),
    )
    sc.mark("guide", "protospacer", "guide identity")
    sc.mark("guide", "oJR160 site", "guide-specific RT site")
    sc.arrow("oJR160", "reverse transcriptase")
    return sc


def first_strand():
    return Construct([
        *ojr160_segments(),
        seg("copied CR1", revcomp(SCAFFOLD[:SITE_OFFSET]), "insert"),
        seg("copied protospacer", "n" * 20, "insert", placeholder=True),
        seg("CCC", "CCC", "tso"),
    ], name="guide first strand")


def switch_scene():
    first = first_strand()
    sc = Scene()
    sc.strand("first", list(first), label="guide first strand")
    sc.anneal("bead TSO", bead_tso(), to="first", pair=("rGrGrG", "CCC"),
              label="barcoded bead TSO", above=True, mod5="bead")
    sc.mark("bead TSO", "cell barcode", "cell barcode + UMI", through="UMI")
    sc.arrow("first", "RT switches template")
    return sc


def amplified_guide_cdna():
    """Guide-bearing 10x cDNA before the guide-specific library PCR."""
    return Construct([
        *bead_tso(),
        seg("protospacer", "N" * 20, "insert", placeholder=True),
        seg("CR1 through capture site", SCAFFOLD[:SITE_OFFSET + len(SCAFFOLD_SITE)],
            "insert"),
        seg("RT adapter complement", revcomp(TENX_RT_ADAPTER), "r1"),
    ], name="barcoded guide cDNA")


def final_library():
    """The approximately 250-bp guide library produced by oJR163/oJR165 PCR."""
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        seg("cell barcode", "B" * 16, "cbc", placeholder=True, inferred=True),
        seg("UMI", "U" * 10, "umi", placeholder=True, inferred=True),
        seg("switch spacer", SWITCH_SPACER, "tso", inferred=True),
        seg("GGG", "GGG", "tso", inferred=True),
        seg("protospacer", "N" * 20, "insert", placeholder=True),
        seg("CR1 before capture site", SCAFFOLD[:SITE_OFFSET], "insert"),
        seg("oJR160 site", SCAFFOLD_SITE, "tso"),
        seg("nested adapter site", revcomp(TENX_RT_ADAPTER[-5:]), "r1"),
        seg("Index 1 / Read 2 arm", nx.INDEX1_PRIMER, "r2"),
        seg("i7 reverse complement", revcomp(I7), "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="5-prime direct-capture guide library")
    problems = sp.verify(lib, sequencing_primers(),
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid direct-capture library: " + "; ".join(problems))
    return lib


def final_scene():
    lib = final_library()
    sc = Scene.duplex(list(lib), label="guide library")
    sc.mark("top", "cell barcode", "Read 1: cell barcode", through="UMI")
    sc.mark("top", "protospacer", "guide identity reached by Read 2")
    sc.labels("top")
    return sc


def sequencing_primers():
    return (sp.TRUSEQ["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"])


def read_layout():
    lengths = {"Read 1": 26, "Index 1 (i7)": 8, "Read 2": 98}
    spans = sp.feature_spans(final_library(), sequencing_primers(), lengths)
    return [(x.read, x.cycles, x.feature.label) for x in spans]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    sgrna=graph.input("sgrna","sgRNA-CR1 in a captured single cell",MaterialKind.BIOLOGICAL,
                      state=tuple(capture_scene().rows()),primary=True)
    capture_primer=graph.input("ojr160","oJR160 guide-specific RT primer, 5 pmol spiked into 68.3 µL 10x RT master mix",
                               MaterialKind.OLIGO,state=tuple(ojr160_segments()))
    rt=graph.input("reverse_transcriptase","template-switching reverse transcriptase",
                   MaterialKind.CATALYST)
    first=graph.product("guide_first_strand","guide first-strand cDNA",MaterialKind.MOLECULAR,
                        state=first_strand())
    graph.react("guide_rt","Reverse-transcribe from the CR1 capture site",
                StepKind.REVERSE_TRANSCRIPTION,
                ((sgrna,InputRole.TEMPLATE),(capture_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    bead_tso_input=graph.input("bead_tso","INFERRED — barcoded 10x 5′ bead template-switch oligo (kit bases not in the inspected sources)",
                               MaterialKind.OLIGO,state=tuple(bead_tso()))
    switched=graph.product("switched_guide_cdna","cell-barcoded guide cDNA",
                           MaterialKind.MOLECULAR,state=tuple(switch_scene().rows()))
    graph.react("template_switch","INFERRED — bead bases: transfer cell barcode and UMI from the bead oligo",
                StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(bead_tso_input,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(switched,))
    pre_f=graph.input("preamp_forward","10x cDNA amplification forward primer (sequence not printed)",MaterialKind.OLIGO)
    pre_r=graph.input("preamp_reverse","10x Non-Poly(dT) primer PN-220106, binds the oJR160 adapter (sequence not printed)",MaterialKind.OLIGO)
    pre_pol=graph.input("preamp_polymerase","cDNA amplification polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_guide_cdna","amplified barcoded guide cDNA",
                            MaterialKind.MOLECULAR,state=amplified_guide_cdna())
    graph.react("preamp","Amplify the barcoded 10x cDNA (11 cycles)",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(pre_f,InputRole.STRUCTURAL),
                 (pre_r,InputRole.STRUCTURAL),(pre_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.unavailable(
                    amplified_guide_cdna(),reason="The compact page does not reconstruct the complete 10x preamplification duplex."))
    spri=graph.input("spri_selection","0.6× left-side SPRI; supernatant then 0.6×–1.2× double-sided SPRI",MaterialKind.SELECTION)
    guide_fraction=graph.product("guide_fraction","short guide-cDNA supernatant",
                                 MaterialKind.MOLECULAR,state=amplified_guide_cdna())
    gene_expression=graph.product("gene_expression_fraction","bead-retained gene-expression cDNA",
                                  MaterialKind.MOLECULAR)
    graph.react("size_select","Separate guide cDNA (~168 bp) from gene-expression cDNA",
                StepKind.SELECTION,
                ((amplified,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),
                (guide_fraction,gene_expression))
    graph.mark_terminal(Terminal(gene_expression.id,Disposition.STORED,
        "the bead-retained fraction proceeds through the separate gene-expression library workflow"))
    ojr163=graph.input("ojr163","oJR163 guide-library forward primer, 0.6 µM",MaterialKind.OLIGO,
                       state=tuple(ojr163_segments()))
    ojr165=graph.input("ojr165","indexed oJR165 nested reverse primer, 0.6 µM",MaterialKind.OLIGO,
                       state=tuple(ojr165_segments()))
    lib_pol=graph.input("library_polymerase","KAPA HiFi ReadyMix",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("guide_library_pcr","Build the nested guide library (4 × 5 ng, 12 cycles 98/70 °C)",StepKind.PCR,
                ((guide_fraction,InputRole.TEMPLATE),(ojr163,InputRole.STRUCTURAL),
                 (ojr165,InputRole.STRUCTURAL),(lib_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    final,reason="The final construct is exact, but the heterogeneous pre-PCR cDNA is represented as one strand."))
    graph.annotate(SourceDiscrepancy(
        "guide_library_pcr", "Guide-library primer concentration unit",
        (SourceClaim("Supplementary Note 3, step 5 (Supplementary Information)",
                     "0.6 uM of primers oJR163 and oJR165; 3 µL of 10 µM each primer per 50-µL reaction."),
         SourceClaim("Methods, ‘Sequencing library preparation’ (5′ sgRNA-CR1cs1 experiment)",
                     "guide molecules were amplified using 0.6 mM oJR163 and oJR166.")),
        "The graph uses 0.6 µM oJR163/oJR165 from Supplementary Note 3.",
        "Supplementary Note 3 is the step-by-step protocol for this standard-CR1 branch and its "
        "reaction table (3 µL of 10 µM in 50 µL) gives 0.6 µM; the Methods sentence concerns the "
        "cs1 branch's oJR166 reaction.",
        "Whether the Methods ‘mM’ is a unit misprint is not stated by the source."))
    spri2=graph.input("library_cleanup","0.8× SPRI cleanup, elute 40 µL",MaterialKind.SELECTION)
    purified=graph.product("purified_library",final.name,MaterialKind.LIBRARY,state=final)
    primer_waste=graph.product("library_cleanup_supernatant","SPRI supernatant (primers)",MaterialKind.WASTE)
    graph.react("library_cleanup_step","Purify the guide library",StepKind.CLEANUP,
                ((library,InputRole.TEMPLATE),(spri2,InputRole.SELECTION)),(purified,primer_waste))
    graph.discard(primer_waste,"Supplementary Note 3 step 5: 0.8× SPRI cleanup removes primers; supernatant discarded")
    graph.sequence(purified,SequencingPlan(library=final,primers=tuple(sequencing_primers()),
        read_lengths={"Read 1":26,"Index 1 (i7)":8,"Read 2":98}))
    graph.declare_reviewed_scope(("sgrna","ojr160","reverse_transcriptase","guide_rt",
      "bead_tso","template_switch","preamp_forward","preamp_reverse","preamp_polymerase",
      "preamp","spri_selection","size_select","ojr163","ojr165","library_polymerase",
      "guide_library_pcr","library","library_cleanup","library_cleanup_step",
      "purified_library"),SOURCE)
    return graph.require_publishable()
