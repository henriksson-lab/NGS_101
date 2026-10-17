"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import nebnext as nb
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import truseq_library
from chemdraw import Construct, Segment, feature
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, closed_circle, da_tailed, molecule, role
from editing_outcomes import _finish_pcr

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
        Segment("i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i7", "sample_index", "whitelist",
                                whitelist="UDiTaS N701-N715 subset")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="UDiTaS library")
    primers = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"], INDEX2_FROM_P5, sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("UDiTaS final library: " + "; ".join(problems))
    return lib, primers


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
                          MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "asymmetrically tagmented genomic DNA",
                           MaterialKind.MOLECULAR, state=dna.state)
    waste = graph.product("tn5_waste", "Tn5 off-products and unusable fragments",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment with the custom UMI/index-bearing adapter",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not form a locus-anchored sequencing molecule")
    anchor = graph.input("anchor_primer", "locus-specific anchor primer", MaterialKind.OLIGO)
    p5 = graph.input("p5_primer", "P5-side amplification primer", MaterialKind.OLIGO)
    polymerase = graph.input("round1_polymerase", "round-1 PCR polymerase",
                            MaterialKind.CATALYST)
    selected = graph.product("selected", "target-anchored UDiTaS amplicon",
                             MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("anchor_pcr", "Amplify from the random Tn5 end to the target anchor",
                StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (anchor, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (selected,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The target-specific intermediate is not represented as an exact Construct."))
    reviewed = ["dna", "tn5_adapter", "tn5", "tagment", "anchor_primer",
                "p5_primer", "round1_polymerase", "anchor_pcr"]
    _finish_pcr(graph, "index", selected, FINAL_LIBRARY, SEQ_PRIMERS,
                SEQUENCING_INTRO, reviewed, label="UDiTaS indexing")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
