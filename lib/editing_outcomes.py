"""Target-anchored genome-editing outcome libraries."""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
import nebnext as nb
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import truseq_library
from chemdraw import Construct, Segment, feature
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, closed_circle, da_tailed, molecule, role


UDITAS_TOP_N501 = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCNNNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
UDITAS_BOTTOM = "CTGTCTCTTATACAC"
UDITAS_ANCHOR_CEP290 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGGACCATGGATGCACTCTGTAAATTCTCAT"
UDITAS_P7_N701 = "CAAGCAGAAGACGGCATACGAGATAGCGGAATGTGACTGGAGTTCAGACGTGT"

# Yin et al. Supplementary Table S11; ordering modifications are retained by
# ``pem_oligos`` while these constants carry the molecular bases.
#
# PEM_BA_UP follows Table S11 (``...CTA CAN NNT NNN ANN NTN NNN AGA TCG...``), i.e.
# 14 random bases in four tracts separated by the fixed T/A/T bases, matching both
# papers' prose ("a 14-bp RMB"; Liu et al. 2022 STAR Protocols, and the defining
# paper's own RMB length comparison).  Do not "correct" it from the later protocol:
# STAR Protocols Table 2 prints the same oligo with only ``NNNTNNNNTNNN`` (10 random
# bases, no A separator), which contradicts Table S11 *and* that paper's own text.
# Table S11 is the defining source and is internally consistent, so it wins.
PEM_BA_UP = "CCACGCGTGCTCTACANNNNTNNNANNNTNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGT"
PEM_BA_LOWER = "TGTAGAGCACGCGTGGNNNNNN"
PEM_P5_I5 = "AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGC"
PEM_P7 = il.P7
PEM_I7_B01 = "CAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGC"
PEM_BIO_RAG1AB = "GGACTGCTGGAGATTGCTCCAGAGAGGGTTTCCCCTCAAAGGAATCCTTC"
# The inline I5 barcode is four bases, not five.  Table S11 lists I5-R1..R4 for three
# loci; only a 4-nt split leaves one shared nested primer per locus (RAG1A/B
# AGGATCTCACCCGGAACAGC, MYC1 AGAGCCATTCTCTGGCTCAG, RAG1C TATCAATATCCCACTGATGTATC).
# A 5-nt split would need the barcode's last base to be A at RAG1A/B and MYC1 but T
# at RAG1C.  Liu et al. 2022 STAR Protocols Table 2 independently colours the same
# four bases (AGCG/GCCT/AGGA/TCAG) as the index; that table is a vector object, so
# read it by rendering page 11 of the article PDF (see ref/MANIFEST.md).
PEM_I5_BARCODES = ("AGCG", "GCCT", "AGGA", "TCAG")  # I5-R1..R4 / I5-Nested-1..4
PEM_NESTED_RAG1AB = "AGGATCTCACCCGGAACAGC"
# The I5-Nested oligos begin two bases inside the canonical Read 1 arm.
PEM_I5_TAIL = il.TRUSEQ_READ1[2:]
PEM_I5_R1_RAG1AB = PEM_I5_TAIL + PEM_I5_BARCODES[0] + PEM_NESTED_RAG1AB

# Lazzarotto et al. Supplementary Protocol, reagent setup.  U is the internal
# deoxyuridine written as /ideoxyU/ in the ordering string returned below.
CHANGE_TRANSPOSON_TOP = "ACGUAGATGTGTATAAGAGACAG"
CHANGE_TRANSPOSON_BOTTOM = "CTGTCTCTTATACACATCTACGT"

# Illumina TruSeq ChIP Sample Preparation Guide 15023092 Rev. B, Table 5.
DISCOVER_AR001 = "ATCACG"


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
    primers = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"], sp.NEXTERA["I2"])
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


def pem_workflow():
    return chain(molecule("repaired genomic DNA from edited cells", (role("bait locus joined to germline, indel, deletion, inversion or translocation partner"),)), (
        ("Sonicate genomic DNA to 300–700 bp",
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


def pem_oligos():
    return (
        "BA-up — 5′-/5phos/CCACGCGTGCTCTACANNNNTNNNANNNTNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGT-NH2(C7)-3′",
        "BA-lower — 5′-TGTAGAGCACGCGTGGNNNNNN-NH2(C7)-3′",
        f"P5-I5 tagging primer — 5′-{PEM_P5_I5}-3′",
        f"P7 tagging primer — 5′-{PEM_P7}-3′",
        f"I7-B01 nested-PCR primer — 5′-{PEM_I7_B01}-3′",
        f"bioRAG1A/B extension primer — 5′-Biotin/{PEM_BIO_RAG1AB}-3′",
        f"I5-R1-RAG1A/B nested primer — 5′-{PEM_I5_R1_RAG1AB}-3′",
    )


def change_library() -> tuple[Construct, tuple]:
    """NEBNext dual-index product recovered outward from a Cas9-opened DNA circle."""
    insert = Construct([
        role("genomic DNA beginning exactly at an in-vitro Cas9 cleavage end"),
    ], name="Cas9-linearized CHANGE-seq insert")
    lib = nb.indexed_library(insert, "CHANGE-seq library")
    return lib, nb.STANDARD_PRIMERS


def change_workflow():
    return chain(molecule("purified high-molecular-weight genomic DNA", (role("genomic DNA"),)), (
        ("Tagment with the custom dU-bearing Tn5 transposome",
         molecule("approximately 400-bp tagmented fragments", (role("palindromic dU/Tn5 end"), role("genomic DNA"), role("palindromic dU/Tn5 end"))),
         "Tn5 installs the same 23-base duplex at both fragment ends; its 19-base mosaic end performs strand transfer and its four-base palindrome is the future ligation end."),
        ("Gap-repair the staggered transposition products",
         molecule("fully duplex tagmented fragments", (role("copied custom end"), role("duplex genomic insert"), role("copied custom end"))),
         "Uracil-tolerant polymerase and Taq ligase fill and seal the Tn5 gaps before the engineered ends are opened."),
        ("Open dU with USER and phosphorylate with T4 PNK",
         molecule("ligation-ready linear fragments", (role("four-base cohesive end"), role("duplex genomic insert"), role("four-base cohesive end"))),
         "USER cleaves at the internal deoxyuridine and PNK supplies ligatable 5-prime phosphates."),
        ("Circularize by dilute intramolecular ligation",
         closed_circle("covalently closed genomic circles", (role("genomic insert"), role("custom circle-closing junction")), closure_label="T4 DNA ligase closure"),
         "The complementary four-base ends close each fragment; low DNA concentration favours intramolecular ligation."),
        ("Digest every molecule that still has a free end",
         closed_circle("exonuclease-purified DNA circles", (role("closed genomic circle"),), closure_label="no free DNA end"),
         "Exonuclease I, lambda exonuclease and Plasmid-Safe DNase remove residual linear or incompletely closed DNA."),
        ("Cut circles with Cas9–sgRNA",
         molecule("target-containing circles linearized at Cas9 sites", (role("Cas9 cleavage end"), role("off-target or on-target genomic DNA"), role("Cas9 cleavage end"))),
         "Only a circle containing a sequence cleaved by the programmed RNP acquires the free ends needed for library construction."),
        ("Add one 3-prime dA to each Cas9-cleaved end",
         da_tailed("dA-tailed Cas9 cleavage product", "Cas9-to-Cas9 genomic fragment"),
         "The newly exposed Cas9 ends are shown before the hairpin adapter is introduced."),
        ("Ligate the NEBNext hairpin adapter",
         molecule("hairpin-adapter-ligated cleavage products", (role("NEBNext hairpin"), role("Cas9-to-Cas9 genomic fragment"), role("NEBNext hairpin"))),
         "A-tailing and the NEBNext adapter preserve each newly exposed Cas9 end as a sequencing boundary."),
        ("Open the hairpins with USER and PCR with dual indexes",
         molecule("dual-index CHANGE-seq library", (role("P5 + i5 + Read 1 arm"), role("genomic DNA starting at a Cas9 cut"), role("Read 2 arm + i7 + P7"))),
         "USER makes the commercial adapter PCR-addressable; twenty cycles with NEBNext i5 and i7 primers complete the Illumina library."),
    ))


def change_oligos():
    i5 = il.nebnext_i5_primer(il.NEBNEXT_I5_SET1["i501"])
    i7 = il.nebnext_i7_primer(il.NEBNEXT_I7_SET1["i701"])
    return (
        "oCRL225 — 5′-/5Phos/ACG/ideoxyU/AGATGTGTATAAGAGACAG-3′",
        f"oCRL226 — 5′-/5Phos/{CHANGE_TRANSPOSON_BOTTOM}-3′",
        "NEBNext Adapter for Illumina — 5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC[dU]ACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′",
        f"NEBNext i501 primer example — 5′-{i5}-3′",
        f"NEBNext i701 primer example — 5′-{i7}-3′",
    )


def discover_library() -> tuple[Construct, tuple]:
    """Historical single-index TruSeq ChIP product with a six-cycle index read."""
    return truseq_library([
        role("MRE11-bound genomic DNA from a Cas9 cut to a random sonication end"),
    ], "DISCOVER-Seq TruSeq ChIP library", dual_index=False, i7_nt=6)


def discover_workflow():
    return chain(molecule("living cells or tissue receiving a CRISPR nuclease", (role("chromosomal on-target and off-target sites"),)), (
        ("Let MRE11 bind newly made double-strand breaks",
         molecule("edited chromatin with repair complexes", (role("Cas9 cleavage end"), role("MRE11 at the DNA end"), role("surrounding chromatin"))),
         "The sampling time is optimized for the delivery system; the detailed protocol suggests approximately 8–16 hours after RNP delivery."),
        ("Formaldehyde-crosslink cells, lyse nuclei and sonicate chromatin",
         molecule("150–500-bp crosslinked chromatin fragments", (role("random sonication end"), role("MRE11-bound DNA reaching a nuclease cut"), role("Cas9 cleavage end"))),
         "Crosslinking preserves the transient MRE11–DNA complex while mechanical shearing makes ChIP-sized fragments."),
        ("Immunoprecipitate MRE11 on Protein A magnetic beads",
         molecule("MRE11-enriched chromatin", (role("antibody–MRE11 complex"), role("DNA adjacent to a cellular DSB"), role("Cas9 cleavage boundary"))),
         "An unedited or non-targeting sample measures background; ChIP enrichment is checked by qPCR before library construction."),
        ("Reverse crosslinks and purify the ChIP DNA",
         molecule("purified break-adjacent DNA", (role("random sonication end"), role("MRE11 footprint region"), role("native nuclease-cleavage end"))),
         "Proteinase K and RNase A remove the captured protein and RNA while retaining the DNA endpoints used by BLENDER."),
        ("End-repair and add one 3-prime dA",
         molecule("A-tailed ChIP fragments", (role("end-repaired random end"), role("break-adjacent insert"), role("end-repaired cleavage end"), role("3-prime dA"))),
         "The historical TruSeq ChIP kit blunts overhangs and then adenylates the repaired fragments."),
        ("Ligate a six-base indexed TruSeq ChIP adapter",
         molecule("adapter-ligated ChIP library", (role("TruSeq Read 1 side"), role("break-adjacent genomic insert"), role("TruSeq Read 2 side with six-base i7 index"))),
         "The discontinued kit supplies single-index adapters; Table 4/5 instructs recording six index bases in the sample sheet."),
        ("PCR-enrich, remove adapter dimers and pool",
         molecule("single-index DISCOVER-Seq library", (role("P5 and Read 1 arm"), role("MRE11-enriched insert"), role("Read 2 arm, i7 and P7"))),
         "Amplification completes the flow-cell ends; the detailed protocol adds a 0.8× bead cleanup when adapter dimers remain."),
    ))


def discover_oligos():
    return (
        f"TruSeq ChIP AR001 sample index — 5′-{DISCOVER_AR001}-3′",
        f"TruSeq universal PCR primer — 5′-{il.TRUSEQ_P5_FULL}-3′",
        f"TruSeq AR001 indexed P7 primer model — 5′-{il.P7}{DISCOVER_AR001}{il.TRUSEQ_READ2}-3′",
    )


def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label):
    forward = graph.input(prefix + "_forward", label + " forward primer", MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " reverse primer", MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                            MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", label + " PCR", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="Primer identities are represented in the final construct; the intermediate template is not an exact Construct."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers),
                                           limitation=""))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_pcr",
                     prefix + "_library"))
    return library


def _uditas_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
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
    graph.react("tagment", "Tagment with the custom UMi/index-bearing adapter",
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
                    spec.FINAL_LIBRARY, reason="The target-specific intermediate is not represented as an exact Construct."))
    reviewed = ["dna", "tn5_adapter", "tn5", "tagment", "anchor_primer",
                "p5_primer", "round1_polymerase", "anchor_pcr"]
    _finish_pcr(graph, "index", selected, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label="UDiTaS indexing")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _pem_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    dna = graph.input("dna", "repaired genomic DNA from edited cells",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("repaired genomic DNA", (role("bait locus and repair junction"),)),
                      primary=True)
    shear = graph.input("shear_condition", "sonication condition", MaterialKind.CONDITION)
    fragments = graph.product("fragments", "300–700-bp repair products",
                              MaterialKind.MOLECULAR, state=dna.state)
    graph.react("sonicate", "Sonicate repaired genomic DNA", StepKind.FRAGMENTATION,
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
                    spec.FINAL_LIBRARY, reason="Nested amplification precedes the exact final tagging-PCR construct."))
    reviewed = ["dna", "shear_condition", "sonicate", "bait_primer", "extension_enzyme",
                "extend", "streptavidin", "capture", "bridge_adapter", "bridge_ligase",
                "ligate_bridge", "nested_forward", "nested_reverse", "nested_polymerase",
                "nested_pcr"]
    _finish_pcr(graph, "tagging", nested, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label="PEM-seq tagging")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _change_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    dna = graph.input("dna", "purified high-molecular-weight genomic DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("purified HMW genomic DNA", (role("genomic DNA"),)),
                      primary=True)
    adapter = graph.input("tn5_adapter", "custom dU-bearing Tn5 transposon",
                          MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "custom-end tagmented genomic fragments",
                           MaterialKind.MOLECULAR, state=dna.state)
    tn5_waste = graph.product("tn5_waste", "Tn5 off-products and unusable fragments",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment with the custom dU-bearing transposome",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not yield a recoverable custom-end genomic circle")
    repair = graph.input("repair_mix", "gap-repair polymerase/ligase mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "fully duplex custom-end fragments",
                             MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("gap_repair", "Fill and seal the staggered transposition products",
                StepKind.END_REPAIR,
                ((tagged, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    user = graph.input("user_enzyme", "USER enzyme", MaterialKind.CATALYST)
    opened = graph.product("opened", "USER-opened cohesive-end fragments",
                           MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("open_du", "Open the internal dU sites", StepKind.DIGESTION,
                ((repaired, InputRole.TEMPLATE), (user, InputRole.CATALYST)), (opened,))
    pnk = graph.input("pnk", "T4 polynucleotide kinase", MaterialKind.CATALYST)
    phosphorylated = graph.product("phosphorylated", "phosphorylated cohesive-end fragments",
                                   MaterialKind.MOLECULAR, state=opened.state)
    graph.react("phosphorylate", "Phosphorylate the opened ligation ends",
                StepKind.PHOSPHORYLATION,
                ((opened, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (phosphorylated,))
    circle_ligase = graph.input("circle_ligase", "circle-closing DNA ligase",
                               MaterialKind.CATALYST)
    circles = graph.product("circles", "covalently closed genomic circles",
                            MaterialKind.MOLECULAR, state=phosphorylated.state)
    graph.react("circularize", "Circularize by dilute intramolecular ligation",
                StepKind.CIRCULARIZATION,
                ((phosphorylated, InputRole.TEMPLATE),
                 (circle_ligase, InputRole.CATALYST)), (circles,),
                geometry=CircularizationGeometry.known(
                    "double-strand intramolecular ligation",
                    statement="CHANGE-seq dilute-ligates the USER-opened ACGT cohesive ends of tagmented genomic fragments into covalently closed circles."))
    exonuclease = graph.input("exonuclease", "linear-DNA exonuclease mixture",
                              MaterialKind.CATALYST)
    purified = graph.product("purified_circles", "exonuclease-purified genomic circles",
                             MaterialKind.MOLECULAR, state=circles.state)
    graph.react("remove_linear", "Digest molecules that retain free ends",
                StepKind.DIGESTION,
                ((circles, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (purified,))
    sgrna = graph.input("sgrna", "programmed sgRNA", MaterialKind.OLIGO)
    cas9 = graph.input("cas9", "Cas9 nuclease", MaterialKind.CATALYST)
    linear = graph.product("linear", "Cas9-linearized target-containing circles",
                           MaterialKind.MOLECULAR, state=purified.state)
    graph.react("cas9_cut", "Reopen target-containing circles with Cas9–sgRNA",
                StepKind.GENOME_EDITING,
                ((purified, InputRole.TEMPLATE), (sgrna, InputRole.STRUCTURAL),
                 (cas9, InputRole.CATALYST)), (linear,))
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed Cas9-linearized products",
                           MaterialKind.MOLECULAR, state=linear.state)
    graph.react("tail", "Add one 3-prime dA to Cas9-cleaved ends", StepKind.TAILING,
                ((linear, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    hairpin = graph.input("nebnext_adapter", "NEBNext hairpin adapter", MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "NEBNext adapter ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated", "hairpin-adapter-ligated cleavage products",
                            MaterialKind.MOLECULAR, state=tailed.state)
    graph.react("ligate_adapter", "Ligate the NEBNext hairpin adapter", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (hairpin, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=reported_adapter_end(
                    "NEBNext hairpin adapter", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The NEBNext hairpin adapter presents the T end used with dA-tailed Cas9 products."))
    opener = graph.input("hairpin_user", "USER enzyme for the hairpin adapter",
                         MaterialKind.CATALYST)
    addressable = graph.product("addressable", "USER-opened PCR-addressable molecules",
                                MaterialKind.MOLECULAR, state=ligated.state)
    graph.react("open_hairpin", "Open adapter hairpins with USER", StepKind.DIGESTION,
                ((ligated, InputRole.TEMPLATE), (opener, InputRole.CATALYST)),
                (addressable,))
    reviewed = ["dna", "tn5_adapter", "tn5", "tagment", "repair_mix", "gap_repair",
                "user_enzyme", "open_du", "pnk", "phosphorylate", "circle_ligase",
                "circularize", "exonuclease", "remove_linear", "sgrna", "cas9",
                "cas9_cut", "tailing_enzyme", "tail", "nebnext_adapter",
                "adapter_ligase", "ligate_adapter", "hairpin_user", "open_hairpin"]
    _finish_pcr(graph, "index", addressable, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label="CHANGE-seq dual-index")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _discover_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    chromatin = graph.input("chromatin", "edited cells or tissue with nuclease-created breaks",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("edited chromatin", (role("chromosomal nuclease cuts"),)),
                            primary=True)
    mre11 = graph.input("mre11_antibody", "MRE11-binding antibody probe",
                        MaterialKind.SELECTION)
    bound = graph.product("mre11_bound", "MRE11-bound break chromatin",
                          MaterialKind.MOLECULAR, state=chromatin.state)
    unbound = graph.product("mre11_unbound", "chromatin without captured MRE11-bound ends",
                            MaterialKind.WASTE)
    graph.react("bind_mre11", "Recognize MRE11 at newly made DNA breaks", StepKind.BINDING,
                ((chromatin, InputRole.TEMPLATE), (mre11, InputRole.SELECTION)),
                (bound, unbound), changes_molecule=False)
    graph.discard(unbound, "does not carry the MRE11-bound break-end signal")
    formaldehyde = graph.input("formaldehyde", "formaldehyde cross-linking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed", "cross-linked MRE11–break chromatin",
                          MaterialKind.MOLECULAR, state=bound.state)
    graph.react("crosslink", "Cross-link MRE11–DNA complexes", StepKind.CROSS_LINKING,
                ((bound, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,))
    shear = graph.input("sonication", "chromatin sonication condition", MaterialKind.CONDITION)
    fragments = graph.product("fragments", "150–500-bp cross-linked chromatin fragments",
                              MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("sonicate", "Sonicate fixed chromatin", StepKind.FRAGMENTATION,
                ((fixed, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    beads = graph.input("protein_a_beads", "Protein A magnetic beads", MaterialKind.SELECTION)
    captured = graph.product("captured", "MRE11-enriched chromatin fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("chip_flowthrough", "MRE11-depleted ChIP flow-through",
                             MaterialKind.WASTE)
    graph.react("chip", "Immunoprecipitate MRE11-bound chromatin", StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "is not retained by MRE11 immunoprecipitation")
    reversal = graph.input("crosslink_reversal", "cross-link reversal and protein-removal condition",
                           MaterialKind.CONDITION)
    dna = graph.product("chip_dna", "purified break-adjacent ChIP DNA",
                        MaterialKind.MOLECULAR, state=captured.state)
    graph.react("purify_dna", "Reverse cross-links and remove protein",
                StepKind.DEPROTEINIZATION,
                ((captured, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (dna,))
    repair = graph.input("repair_enzyme", "TruSeq ChIP end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired ChIP fragments",
                             MaterialKind.MOLECULAR, state=dna.state)
    graph.react("repair", "End-repair ChIP DNA", StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed ChIP fragments",
                           MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("tail", "Add one 3-prime dA", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("truseq_adapter", "six-base indexed TruSeq ChIP adapter",
                          MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "TruSeq adapter ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated", "adapter-ligated ChIP fragments",
                            MaterialKind.MOLECULAR, state=tailed.state)
    graph.react("ligate", "Ligate the indexed TruSeq ChIP adapter", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=reported_adapter_end(
                    "TruSeq ChIP adapter", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The named TruSeq ChIP workflow ligates its indexed adapter to dA-tailed fragments."))
    reviewed = ["chromatin", "mre11_antibody", "bind_mre11", "formaldehyde", "crosslink",
                "sonication", "sonicate", "protein_a_beads", "chip", "crosslink_reversal",
                "purify_dna", "repair_enzyme", "repair", "tailing_enzyme", "tail",
                "truseq_adapter", "adapter_ligase", "ligate"]
    _finish_pcr(graph, "library", ligated, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label="DISCOVER-Seq library")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _attach_graph(key, spec):
    builders = {"uditas": _uditas_graph, "pem-seq": _pem_graph,
                "change-seq": _change_graph, "discover-seq": _discover_graph}
    spec.protocol_graph = lambda: builders[key](spec)
    return spec


def _spec(key, **fields):
    return _attach_graph(key, SimpleNamespace(**fields))


def get(key="uditas"):
    if key == "discover-seq":
        lib, primers = discover_library()
        return _spec(key,
            TITLE="DISCOVER-Seq — MRE11 ChIP-seq at CRISPR cleavage sites",
            NOTES="01_discover_seq.html",
            SOURCE='<a href="https://doi.org/10.1126/science.aav9023">Wienert et al., Science (2019)</a>; stepwise implementation from the authors’ <a href="https://doi.org/10.1038/s41596-020-0309-5">Nature Protocols procedure</a> and the named TruSeq ChIP kit guide.',
            SUMMARY="MRE11 recruited to newly cut chromosomal DNA is fixed in place and immunoprecipitated; ordinary paired-end ChIP-seq then preserves the sharp genomic endpoints that BLENDER recognizes as nuclease cuts.",
            CAVEAT="DISCOVER-Seq does not use a custom capture oligo or UMI. Its defining implementation delegates library construction to the now-discontinued single-index TruSeq ChIP kit; index identities are public, but the kit did not publish every supplied oligo as an ordering sequence.",
            oligos=discover_oligos, workflow=discover_workflow, FINAL_LIBRARY=lib,
            FINAL_CAPTION="Either genomic orientation can be captured. One library strand is shown from a random sonication boundary through MRE11-enriched DNA to the cellular nuclease-cleavage boundary; the six-base i7 sample index comes from the historical TruSeq ChIP adapter.",
            SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 100, "Index 1 (i7)": 6, "Read 2": 100},
            SEQUENCING_INTRO="The detailed protocol specifies 100-bp paired-end sequencing. Standard TruSeq Read 1, Index 1 and Read 2 primers bind the final library; paired fragment endpoints, not a barcode read, provide the cut-site signal.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
    if key == "change-seq":
        lib, primers = change_library()
        return _spec(key,
            TITLE="CHANGE-seq — circularization for high-throughput analysis of nuclease genome-wide effects",
            NOTES="01_change_seq.html",
            SOURCE='<a href="https://doi.org/10.1038/s41587-020-0555-7">Lazzarotto et al., Nature Biotechnology (2020)</a>; exact reaction order and custom transposon from the paper’s Supplementary Protocol.',
            SUMMARY="A dU-bearing custom Tn5 transposome converts genomic DNA into exonuclease-resistant circles; Cas9 selectively reopens circles at on- and off-target sites, which become the ends of a dual-index Illumina library.",
            CAVEAT="The two custom Tn5 oligos are paper-specific. The downstream NEBNext adapter and example Set 1 index primers are canonical vendor sequences for the commercial products named by the protocol.",
            oligos=change_oligos, workflow=change_workflow, FINAL_LIBRARY=lib,
            FINAL_CAPTION="One strand of the symmetric Cas9-linearized product is shown. Read 1 and Read 2 both enter genomic DNA from newly exposed nuclease-cleavage ends; the custom Tn5 circle-closing junction is internal to the selected insert rather than a sequencing adapter.",
            SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 150},
            SEQUENCING_INTRO="The defining protocol uses 150-bp paired-end MiSeq or NextSeq sequencing. Standard TruSeq primers bind the NEBNext-derived final construct; the paired genomic reads begin at the two Cas9-cleaved ends and point inward.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
    if key == "pem-seq":
        lib, primers = pem_library()
        return _spec(key,
            TITLE="PEM-seq — primer-extension-mediated sequencing", NOTES="01_pem_seq.html",
            SOURCE='<a href="https://doi.org/10.1038/s41421-019-0088-8">Yin et al., Cell Discovery (2019)</a>; exact universal oligos cross-checked against the authors’ detailed <a href="https://doi.org/10.1016/j.xpro.2021.101088">STAR Protocols procedure</a>.',
            SUMMARY="A single biotinylated extension from one side of a bait DSB captures individual repair products; an on-bead bridge adapter assigns a 14-random-base molecular barcode before nested and indexing PCR.", CAVEAT="",
            oligos=pem_oligos, workflow=pem_workflow, FINAL_LIBRARY=lib,
            FINAL_CAPTION="RAG1A/B with I5-R1 and I7-B01 is shown. Four random tracts separated by the fixed T/A/T bases form one 14-base RMB; Read 1 enters from the nested bait primer and Read 2 enters from the bridge side.",
            SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 6, "Read 2": 150},
            SEQUENCING_INTRO="The defining protocol uses 2×150-bp HiSeq sequencing. Read 1 first reports the four-base inline I5 sample barcode and bait sequence; Read 2 approaches the multipart 14-base RMB from the bridge side, while Index 1 reports the six-base i7 barcode.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
    if key != "uditas":
        raise KeyError(key)
    lib, primers = uditas_library()
    return _spec(key,
        TITLE="UDiTaS — Uni-Directional Targeted Sequencing", NOTES="01_uditas.html",
        SOURCE='<a href="https://doi.org/10.1186/s12864-018-4561-9">Giannoukos et al., BMC Genomics (2018)</a>.',
        SUMMARY="A single custom Tn5 adapter gives every genomic fragment a random endpoint, a 10-nt UMI and an i5 pooling index; locus-specific PCR then captures indels, deletions, inversions and translocation junctions from one anchor.", CAVEAT="",
        oligos=lambda: (
            f"Tn5-A top, N501 example — 5′-{UDITAS_TOP_N501}-3′",
            "Tn5-A bottom — 5′-/5Phos/CTGTCTCTTATACA[ddC]-3′",
            f"CEP290 guide-323 anchor — 5′-{UDITAS_ANCHOR_CEP290}-3′",
            f"round-1/2 P5 primer — 5′-{il.P5}-3′",
            f"round-2 P7/i7 N701 example — 5′-{UDITAS_P7_N701}-3′"),
        workflow=workflow, FINAL_LIBRARY=lib,
        FINAL_CAPTION="N501/N701 example. The left custom Tn5 arm carries the UMI and pooling index; the right end is created from the locus-specific anchor by two PCR rounds.",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Index 1 (i7)": 8, "Index 2 (i5)": 18},
        SEQUENCING_INTRO="The published MiSeq run uses a 300-base Read 1 from the random Tn5 end, an 8-base Index 1, and an 18-base reverse-complement Index 2 that traverses the 10-nt UMI followed by the 8-nt i5 pooling barcode.",
        SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
