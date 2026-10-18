"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role

UDITAS_TOP_N501 = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCNNNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"

UDITAS_BOTTOM = "CTGTCTCTTATACAC"

UDITAS_ANCHOR_CEP290 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGGACCATGGATGCACTCTGTAAATTCTCAT"

UDITAS_P7_N701 = "CAAGCAGAAGACGGCATACGAGATAGCGGAATGTGACTGGAGTTCAGACGTGT"


# Additional file 2 section 13 run sheet: MiSeq (a forward-strand instrument). Its
# Appendix A "i5 Barcode READ (Tx10 for UMI)" column prints the Index 2 read as the
# i5 barcode as written in the top oligo followed by the ten UMI cycles
# (e.g. TAGATCGC + T10 for N501), i.e. primed from the grafted P5 oligo.
INDEX2_FROM_P5 = sp.custom(
    "Index 2 (i5)", "flow-cell P5 oligo (MiSeq forward-strand workflow)", il.P5,
    "Giannoukos et al. 2018, Additional file 2 section 13 and Appendix A",
    "The 18-cycle i5 read reports the 8-nt pooling barcode, then the 10-nt UMI.")


def uditas_library() -> tuple[Construct, tuple]:
    """Published asymmetric UDiTaS library; primer landing sites validate at construction."""
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5 pooling barcode", "TAGATCGC", "cbc",
                feature=feature("pool_index_i5", "sample_index", "whitelist",
                                whitelist="UDiTaS N501-N508")),
        Segment("10-nt UMI", "U" * 10, "umi", placeholder=True,
                feature=feature("molecule_umi", "umi", "random")),
        Segment("Nextera S5", nx.S5, "s5"), Segment("mosaic end", nx.ME, "me"),
        role("genomic DNA from random Tn5 insertion to target anchor"),
        Segment("anchor-arm leading A", "A"),
        Segment("Index 1 + Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 sample index reverse complement",
                revcomp(UDITAS_P7_N701[len(il.P7):len(il.P7)+8]), "cbc",
                feature=feature("sample_index_i7", "sample_index", "whitelist",
                                whitelist="UDiTaS N701-N715 subset")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="UDiTaS library")
    primers = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"], INDEX2_FROM_P5, sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("UDiTaS final library: " + "; ".join(problems))
    return lib, primers


def adapter_scene():
    top=[*list(FINAL_LIBRARY)[:4],Segment("unpaired mosaic-end prefix",nx.ME[:-len(UDITAS_BOTTOM)]),
         Segment("paired mosaic-end suffix",nx.ME[-len(UDITAS_BOTTOM):])]
    sc=Scene();sc.strand("Tn5-A top",top)
    sc.anneal("Tn5-A bottom",[Segment("bottom /5Phos/ and terminal ddC",UDITAS_BOTTOM,
        note="5′ phosphate; terminal C is dideoxycytidine")],to="Tn5-A top",
        pair=("bottom /5Phos/ and terminal ddC","paired mosaic-end suffix"),mod5="p",mod3="ddC")
    return sc


def selected_template():
    return Construct([*list(FINAL_LIBRARY)[:5],role("genomic sequence between Tn5 insertion and the CEP290 anchor"),
        Segment("copied CEP290 anchor target",revcomp(UDITAS_ANCHOR_CEP290[34:])),
        Segment("copied anchor Read2 tail",revcomp(il.TRUSEQ_READ2),"r2")],
        name="round1 CEP290 target-anchored amplicon")


def workflow():
    return chain(molecule("edited genomic DNA", (role("wild-type, indel or rearranged target locus"),)), (
        ("Tagment with one custom UMI/index-bearing Tn5 adapter",
         molecule("asymmetrically tagmented genomic DNA", (role("P5 + i5 + 10-nt UMI + S5/ME"), role("random genomic fragment")),
                  junctions=(("P5 + i5 + 10-nt UMI + S5/ME", "random genomic fragment", "Tn5 strand transfer"),)),
         "The adapter marks a random fragment end before any target-selective amplification."),
        ("Round 1: amplify from P5 to one locus-specific anchor",
         molecule("unidirectional target amplicon", (role("P5/i5/UMI/S5/ME end"), role("genomic DNA spanning an editing outcome"), role("Read 2-tailed target anchor"))),
         "Only molecules containing the chosen anchor are recovered; their random Tn5 endpoints reduce length-dependent PCR bias."),
        ("Round 2: add the indexed P7 end",
         molecule("dual-index UDiTaS library", (role("P5"), role("i5 pooling barcode"), role("10-nt UMI"), role("S5/ME"), role("target-spanning insert"), role("Read 2 arm"), role("i7 sample index"), role("P7"))),
         "A common P5 primer and indexed P7 primer complete the sequencing library."),
    ))


TITLE = 'UDiTaS — Uni-Directional Targeted Sequencing'
NOTES = '01_uditas.html'
SOURCE = '<a href="https://doi.org/10.1186/s12864-018-4561-9">Giannoukos et al., BMC Genomics (2018)</a>.'
SUMMARY = 'A single custom Tn5 adapter gives every genomic fragment a random endpoint, a 10-nt UMI and an i5 pooling index; locus-specific PCR then captures indels, deletions, inversions and translocation junctions from one anchor.'
CAVEAT = ''
FINAL_CAPTION = 'N501/N701 example. The left custom Tn5 arm carries the UMI and pooling index; the right end is created from the locus-specific anchor by two PCR rounds.'
READ_LENGTHS = {'Read 1': 300, 'Index 1 (i7)': 8, 'Index 2 (i5)': 18, 'Read 2': 300}
READ_USES = {'Read 1': 'include', 'Index 1 (i7)': 'demultiplex',
             'Index 2 (i5)': 'extract_and_demultiplex', 'Read 2': 'include'}
SEQUENCING_INTRO = ('Additional file 2 (section 13) gives the MiSeq run sheet: Read 1 300 bases from the random Tn5 end, '
                    'Index 1 8 bases, Read 2 300 bases from the anchor-primer Read 2 tail, and Index 2 18 bases. '
                    'Its Appendix A prints each Index 2 read as the 8-nt i5 pooling barcode followed by the 10 UMI cycles, '
                    'the forward-strand order primed from the grafted P5 oligo; the article pipeline aligns the paired reads.')
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = uditas_library()
def oligos():
    return (f'Tn5-A top, N501 example — 5′-{UDITAS_TOP_N501}-3′', 'Tn5-A bottom — 5′-/5Phos/CTGTCTCTTATACA[ddC]-3′', f'CEP290 guide-323 anchor — 5′-{UDITAS_ANCHOR_CEP290}-3′', f'round-1/2 P5 primer — 5′-{il.P5}-3′', f'round-2 P7/i7 N701 example — 5′-{UDITAS_P7_N701}-3′')

def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    dna = graph.input("dna", "edited genomic DNA", MaterialKind.BIOLOGICAL,
                      state=molecule("edited genomic DNA", (role("edited target locus"),)),
                      primary=True)
    adapter = graph.input("tn5_adapter", "custom UMI/index-bearing Tn5 adapter",
                          MaterialKind.ADAPTER,state=adapter_scene())
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "asymmetrically tagmented genomic DNA",
                           MaterialKind.MOLECULAR, state=Construct([*list(FINAL_LIBRARY)[:5],
                             role("random genomic fragment; anchor coordinate not assigned")],name="Tn5-tagged genomic fragment"))
    waste = graph.product("tn5_waste", "Tn5 off-products and unusable fragments",
                          MaterialKind.MOLECULAR)
    graph.react("tagment", "Tagment with the custom UMI/index-bearing adapter",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.mark_terminal(Terminal(waste.id,Disposition.UNRESOLVED,
        "Lack of an anchored sequencing molecule is an amplification/analysis outcome; no physical end-class separation is established at tagmentation."))
    anchor = graph.input("anchor_primer", "CEP290 guide323 anchor primer (article Table2)", MaterialKind.OLIGO,
                         state=Construct([Segment("CEP290 anchor",UDITAS_ANCHOR_CEP290)]))
    p5 = graph.input("p5_primer", "P5-side amplification primer", MaterialKind.OLIGO,
                     state=Construct([Segment("P5 primer",il.P5)]))
    polymerase = graph.input("round1_polymerase", "2X Platinum SuperFi PCR Master Mix (round1 with TMAC)",
                            MaterialKind.CATALYST)
    selected = graph.product("selected", "target-anchored UDiTaS amplicon",
                             MaterialKind.MOLECULAR, state=selected_template())
    graph.react("anchor_pcr", "Amplify from the random Tn5 end to the target anchor",
                StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (anchor, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (selected,), geometry=PCRGeometry.unavailable(
                    tagged.state, reason="The printed locus-specific primer is declared, but its exact genomic site in the tagmented population is not supplied; final-library coordinates cannot substitute for the actual round1 template."))
    reviewed = ["dna", "tn5_adapter", "tn5", "tagment", "anchor_primer",
                "p5_primer", "round1_polymerase", "anchor_pcr"]
    forward=graph.input("index_forward","AppendixA P5/i5 amplification primer",MaterialKind.OLIGO,
        state=Construct([Segment("P5 primer",il.P5)]))
    reverse=graph.input("index_reverse","AppendixA i7_N701_SBS12 primer example",MaterialKind.OLIGO,
        state=Construct([Segment("N701 indexed primer",UDITAS_P7_N701)]))
    pol=graph.input("index_polymerase","2X Platinum SuperFi PCR Master Mix",MaterialKind.CATALYST)
    library=graph.product("index_library","UDiTaS indexing sequencing library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("index_pcr","Round2 indexing with P5 and i7_N701_SBS12",StepKind.PCR,
        ((selected,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
         (reverse,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),
        geometry=PCRGeometry.exact(selected.state,PCRPrimer("P5",il.P5),
            PCRPrimer("i7_N701_SBS12",UDITAS_P7_N701),
            statement="Additional file2 §7 and AppendixA primers bind the actual round1 amplicon; CEP290 anchor is the article Table2 example. Genomic extent remains molecule dependent."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS,read_uses=READ_USES))
    reviewed.extend(("index_forward","index_reverse","index_polymerase","index_pcr","index_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
