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

# Yin et al. 2019 Supplementary Table S11 has 14 random bases in four tracts
# separated by fixed T/A/T. Liu et al. 2022 STAR Protocols Table 2 prints a
# 10-random-base variant despite its own prose calling the RMB 14 bp. The
# defining 2019 table and both prose descriptions support this 14-base model.
PEM_BA_UP = "CCACGCGTGCTCTACANNNNTNNNANNNTNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGT"

PEM_BA_LOWER = "TGTAGAGCACGCGTGGNNNNNN"

PEM_P5_I5 = "AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGC"

PEM_P7 = il.P7

PEM_I7_B01 = "CAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGC"

PEM_BIO_RAG1AB = "GGACTGCTGGAGATTGCTCCAGAGAGGGTTTCCCCTCAAAGGAATCCTTC"

PEM_I5_BARCODES = ("AGCG", "GCCT", "AGGA", "TCAG")  # I5-R1..R4 / I5-Nested-1..4

PEM_NESTED_RAG1AB = "AGGATCTCACCCGGAACAGC"

PEM_I5_TAIL = il.TRUSEQ_READ1[2:]

PEM_I5_R1_RAG1AB = PEM_I5_TAIL + PEM_I5_BARCODES[0] + PEM_NESTED_RAG1AB


def pem_library() -> tuple[Construct, tuple]:
    """RAG1A/B, I5-R1 and I7-B01 example from the original PEM-seq table."""
    umi_kw = dict(tag="umi", placeholder=True)
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        Segment("4-nt inline I5-R1 sample index", PEM_I5_BARCODES[0], "cbc",
                feature=feature("sample_index_i5", "sample_index", "whitelist",
                                whitelist="PEM-seq I5-R1–R4")),
        role("bait-primer sequence and genomic DNA through the repair junction to the prey/shear end"),
        Segment("bridge fixed anchor", "CCACGCGTGCTCTACA"),
        Segment("RMB part 1", "N" * 4, feature=feature("rmb_1", "umi", "random", group="pem_rmb", part="1/4"), **umi_kw),
        Segment("RMB separator T1", "T"),
        Segment("RMB part 2", "N" * 3, feature=feature("rmb_2", "umi", "random", group="pem_rmb", part="2/4"), **umi_kw),
        Segment("RMB separator A", "A"),
        Segment("RMB part 3", "N" * 3, feature=feature("rmb_3", "umi", "random", group="pem_rmb", part="3/4"), **umi_kw),
        Segment("RMB separator T2", "T"),
        Segment("RMB part 4", "N" * 4, feature=feature("rmb_4", "umi", "random", group="pem_rmb", part="4/4"), **umi_kw),
        Segment("bridge adapter leading A", "A"),
        Segment("Index 1 + Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("I7-B01 sample index reverse complement", "ATCACG", "cbc",
                feature=feature("sample_index_i7", "sample_index", "whitelist",
                                whitelist="PEM-seq I7-B/R series")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="PEM-seq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("PEM-seq final library: " + "; ".join(problems))
    return lib, primers


def workflow():
    return chain(molecule("repaired genomic DNA from edited cells", (role("bait locus joined to germline, indel, deletion, inversion or translocation partner"),)), (
        ("Sonicate genomic DNA to 0.3–2 kb",
         molecule("randomly sheared repair products", (role("bait-primer site"), role("repair junction and prey DNA"), role("random shear end"))),
         "The bait-side primer site is fixed; the other end may contain local repair sequence or a genome-wide translocation partner."),
        ("Perform one extension from a 5-prime-biotinylated bait primer",
         molecule("single-copy biotinylated extension strand", (role("5-prime biotin bait primer"), role("one copied repair-product molecule"))),
         "One extension product—not an 80-cycle linear-amplification family—is made from each recoverable template before molecular barcoding."),
        ("Remove free primer and capture the extension strand on streptavidin",
         molecule("bead-bound bait extension", (role("streptavidin-bound 5-prime biotin"), role("single-stranded copied repair product"), role("free 3-prime end"))),
         "Affinity capture selects strands that contain the chosen bait-primer site."),
        ("Ligate the blocked bridge adapter to the strand's 3-prime end",
         molecule("RMB-tagged captured strand", (role("bait extension"), role("bridge fixed anchor"), role("14 random bases separated by T/A/T"), role("Read 2-side adapter")),
                  junctions=(("bait extension", "bridge fixed anchor", "single-strand ligation"),)),
         "The four random tracts jointly form the 14-base random molecular barcode; fixed T/A/T separators are not UMI bases."),
        ("Nested PCR adds the inline I5 barcode and I7 index arm",
         molecule("target-enriched PEM-seq amplicon", (role("Read 1 tail + 4-nt inline sample index + nested bait primer"), role("repair-product insert"), role("RMB bridge"), role("six-base i7 arm"))),
         "The nested primer moves closer to the bait DSB; the opposite primer copies the bridge and introduces an i7 sample index."),
        ("Tagging PCR completes P5 and P7",
         molecule("dual-ended PEM-seq library", (role("P5 and Read 1 site"), role("inline sample index and bait-to-prey insert"), role("14-base RMB"), role("Read 2 and Index 1 site"), role("i7 sample index and P7"))),
         "The published P5-I5 and P7 primers complete a 2×150-bp HiSeq library."),
    ))


def oligos():
    return (
        "BA-up — 5′-/5phos/CCACGCGTGCTCTACANNNNTNNNANNNTNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGT-NH2(C7)-3′",
        "BA-lower — 5′-TGTAGAGCACGCGTGGNNNNNN-NH2(C7)-3′",
        f"P5-I5 tagging primer — 5′-{PEM_P5_I5}-3′",
        f"P7 tagging primer — 5′-{PEM_P7}-3′",
        f"I7-B01 nested-PCR primer — 5′-{PEM_I7_B01}-3′",
        f"bioRAG1A/B extension primer — 5′-Biotin/{PEM_BIO_RAG1AB}-3′",
        f"I5-R1-RAG1A/B nested primer — 5′-{PEM_I5_R1_RAG1AB}-3′",
    )


TITLE = 'PEM-seq — primer-extension-mediated sequencing'
NOTES = '01_pem_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/s41421-019-0088-8">Yin et al., Cell Discovery (2019)</a>; exact universal oligos cross-checked against the authors’ detailed <a href="https://doi.org/10.1016/j.xpro.2021.101088">STAR Protocols procedure</a>.'
SUMMARY = 'A single biotinylated extension from one side of a bait DSB captures individual repair products; an on-bead bridge adapter assigns a 14-random-base molecular barcode before nested and indexing PCR.'
CAVEAT = ("<b>Published barcode conflict.</b> Yin et al. 2019 Supplementary Table S11 prints "
          "the bridge adapter random tract as NNNNTNNNANNNTNNNN (14 random bases, with "
          "fixed T/A/T separators). Liu et al. 2022 STAR Protocols Table 2 prints "
          "NNNTNNNNTNNN (10 random bases), although that protocol’s prose calls it a "
          "14-bp RMB. This model uses the defining 2019 Table S11 sequence because it "
          "agrees with both papers’ stated 14-base length. Whether the later printed "
          "10-base oligo was used in any experimental lot remains unresolved.")
FINAL_CAPTION = 'RAG1A/B with I5-R1 and I7-B01 is shown. Four random tracts separated by the fixed T/A/T bases form one 14-base RMB; Read 1 enters from the nested bait primer and Read 2 enters from the bridge side.'
READ_LENGTHS = {'Read 1': 150, 'Index 1 (i7)': 6, 'Read 2': 150}
SEQUENCING_INTRO = 'The defining protocol uses 2×150-bp HiSeq sequencing. Read 1 first reports the four-base inline I5 sample barcode and bait sequence; Read 2 approaches the multipart 14-base RMB from the bridge side, while Index 1 reports the six-base i7 barcode.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = pem_library()

def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    dna = graph.input("dna", "repaired genomic DNA from edited cells",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("repaired genomic DNA", (role("bait locus and repair junction"),)),
                      primary=True)
    shear = graph.input("shear_condition", "sonication condition", MaterialKind.CONDITION)
    fragments = graph.product("fragments", "0.3–2-kb repair products",
                              MaterialKind.MOLECULAR, state=dna.state)
    graph.react("sonicate", "Sonicate repaired genomic DNA to 0.3–2 kb", StepKind.FRAGMENTATION,
                ((dna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    bait = graph.input("bait_primer", "5-prime-biotinylated bait primer", MaterialKind.OLIGO)
    extension_enzyme = graph.input("extension_enzyme", "single-extension polymerase",
                                   MaterialKind.CATALYST)
    extended = graph.product("extended", "single-copy biotinylated extension strands",
                             MaterialKind.MOLECULAR, state=fragments.state)
    graph.react("extend", "Perform one extension from the biotinylated bait primer",
                StepKind.EXTENSION,
                ((fragments, InputRole.TEMPLATE), (bait, InputRole.STRUCTURAL),
                 (extension_enzyme, InputRole.CATALYST)), (extended,))
    beads = graph.input("streptavidin", "streptavidin beads", MaterialKind.SELECTION)
    captured = graph.product("captured", "bead-bound bait extensions",
                             MaterialKind.MOLECULAR, state=extended.state)
    rejected = graph.product("uncaptured", "uncaptured extension material", MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated extension strands", StepKind.AFFINITY_CAPTURE,
                ((extended, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "lacks a captured biotinylated bait extension")
    bridge = graph.input("bridge_adapter", "blocked random-molecular-barcode bridge adapter",
                         MaterialKind.ADAPTER)
    ligase = graph.input("bridge_ligase", "single-strand bridge ligase",
                         MaterialKind.CATALYST)
    tagged = graph.product("tagged", "RMB-tagged captured strands",
                           MaterialKind.MOLECULAR, state=captured.state)
    graph.react("ligate_bridge", "Ligate the bridge adapter to the free 3-prime end",
                StepKind.LIGATION,
                ((captured, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (tagged,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source gives the bridge oligos and blocking groups, but this graph does not yet reconstruct the on-bead single-strand ligation substrate as a Scene."))
    nested_f = graph.input("nested_forward", "nested bait-side PCR primer", MaterialKind.OLIGO)
    nested_r = graph.input("nested_reverse", "bridge-side indexed PCR primer", MaterialKind.OLIGO)
    nested_polymerase = graph.input("nested_polymerase", "nested-PCR polymerase",
                                    MaterialKind.CATALYST)
    nested = graph.product("nested", "nested target-enriched PEM-seq amplicons",
                           MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("nested_pcr", "Nested-PCR the bait and bridge ends", StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (nested_f, InputRole.STRUCTURAL),
                 (nested_r, InputRole.STRUCTURAL),
                 (nested_polymerase, InputRole.CATALYST)), (nested,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="Nested amplification precedes the exact final tagging-PCR construct."))
    reviewed = ["dna", "shear_condition", "sonicate", "bait_primer", "extension_enzyme",
                "extend", "streptavidin", "capture", "bridge_adapter", "bridge_ligase",
                "ligate_bridge", "nested_forward", "nested_reverse", "nested_polymerase",
                "nested_pcr"]
    _finish_pcr(graph, "tagging", nested, FINAL_LIBRARY, SEQ_PRIMERS,
                SEQUENCING_INTRO, reviewed, label="PEM-seq tagging")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
