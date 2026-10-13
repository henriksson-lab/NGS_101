"""Shared chemistry for structured and heavily modified RNA libraries.

The important invariant here is that end repair, demethylation and reverse-transcriptase
choice are separate operations.  A caller cannot describe AlkB treatment as generic
"library preparation", or silently turn a modification-preserving method into a
demethylating one.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from chemdraw import Construct, Scene, Segment, Workflow, feature
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role


def exact(name: str, seq: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, seq, tag, **kw)


@dataclass(frozen=True)
class StructuredRNAInput:
    name: str
    end_state: str
    demethylation: str = "none"
    preserve_modification_signatures: bool = True

    def __post_init__(self) -> None:
        if self.demethylation != "none" and self.preserve_modification_signatures:
            raise ValueError("demethylation cannot preserve the removed modification signatures")

    def scene(self) -> Scene:
        sc = molecule(self.name, (role("structured RNA body"), role(self.end_state)))
        sc.note("molecule", f"demethylation: {self.demethylation}")
        return sc


Y3 = exact("Y-3-prime-AD", "GTATCCAGTTGGAATTCTCGGGTGCCAAGG", "r2")
Y5_A = exact("Y-5-prime-AD-A", "GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGGA", "r1")
Y5_G = exact("Y-5-prime-AD-G", "GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGGG", "r1")
Y5_C = exact("Y-5-prime-AD-C", "GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGGC", "r1")
Y5_U = exact("Y-5-prime-AD-U", "GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGGT", "r1")

DM_TEMPLATE = exact("DM-TGIRT RNA template", "AGATCGGAAGAGCACACGTCTAGTTCTACAGTCCGACGATC", "r1")
DM_FORWARD = exact("DM-TGIRT multiplex PCR primer", "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC", "p5")
DM_REVERSE = exact("DM-TGIRT indexed PCR primer", "CAAGCAGAAGACGGCATACGAGAT[BARCODE]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", "p7", placeholder=True)

MIM_ADAPTER_I1 = exact("mim 3-prime adapter I1", "GATATCGTCAAGATCGGAAGAGCACACGTCTGAA", "r2")
MIM_RT = exact("mim primer-dependent RT oligo", "[5-prime-pRN]AGATCGGAAGAGCGTCGTGTAGGGAAAGAG[iSp18]GTGACTGGAGTTCAGACGTGTGCTC", "r2", placeholder=True)
MIM_FORWARD = exact("mim library forward primer", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC", "p5")
MIM_REVERSE = exact("mim indexed reverse primer", "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTG", "p7", placeholder=True,
                    feature=feature("sample_index", "sample_index", "random"))

PANDORA_RA3 = exact("NEBNext 3-prime SR adapter", "AGATCGGAAGAGCACACGTCT", "r2")


def _yamat() -> Workflow:
    source = StructuredRNAInput("mature aminoacylated tRNA", "3-prime NCCA end")
    return chain(source.scene(), (
        ("Deacylate mature tRNA at alkaline pH",
         molecule("deacylated mature tRNA", (role("paired acceptor stem"), role("free 3-prime NCCA"))),
         "The free 3′-OH and conserved NCCA overhang define the ligation substrate."),
        ("Anneal the four Y-5′ adapters with Y-3′-AD",
         molecule("Y-adapter–tRNA complex", (Y5_A, role("mature tRNA"), Y3)),
         "The four terminal adapter variants complement the four possible discriminator bases."),
        ("Seal both nicks with T4 RNA ligase 2",
         molecule("covalently Y-adapted tRNA", (Y5_A, role("mature tRNA"), Y3),
                  junctions=(("Y-5-prime-AD-A", "mature tRNA", "Rnl2 nick ligation"),
                             ("mature tRNA", "Y-3-prime-AD", "Rnl2 nick ligation"))),
         "Rnl2 seals the RNA/RNA and RNA/DNA nicks in the annealed acceptor-stem complex."),
        ("Reverse-transcribe, PCR and gel-select the tRNA product",
         molecule("YAMAT-seq library", (role("P5 and small-RNA read arm"),
                  role("mature tRNA cDNA"), role("indexed P7 arm"))),
         "The defining experiment used TruSeq Small RNA reagents and 100-nt single-read HiSeq sequencing."),
    ))


def _dm() -> Workflow:
    source = StructuredRNAInput("deacylated mature tRNA", "free 3-prime CCA",
                                "wild-type AlkB plus AlkB-D135S", False)
    return chain(source.scene(), (
        ("Demethylate Watson–Crick-face m1A, m3C and m1G",
         molecule("demethylated tRNA", (role("structured tRNA with relieved RT blocks"),
                                         role("free 3-prime CCA"))),
         "Wild-type AlkB targets m1A/m3C and D135S improves m1G removal."),
        ("Template-switch TGIRT from a blocked RNA/DNA starter duplex",
         molecule("full-length tRNA cDNA", (DM_TEMPLATE, role("tRNA-derived cDNA"))),
         "The DNA primer carries a one-base 3′ T overhang; TGIRT switches onto the tRNA CCA end."),
        ("Circularize the single-stranded cDNA",
         molecule("circular DM-TGIRT cDNA", (role("CircLigase-II-closed cDNA"),)),
         "Circularization places the two PCR sites around one continuous cDNA template."),
        ("Amplify with the multiplex and indexed primers",
         molecule("DM-tRNA-seq library", (DM_FORWARD, role("tRNA cDNA"), DM_REVERSE)),
         "The defining paper sequenced the indexed products on an Illumina HiSeq."),
    ))


def _mim() -> Workflow:
    source = StructuredRNAInput("60–100-nt native mature tRNA", "native 3-prime CCA")
    return chain(source.scene(), (
        ("Ligate one of four preadenylated barcoded 3′ adapters",
         molecule("adapter-ligated tRNA", (role("native modified tRNA"), MIM_ADAPTER_I1),
                  junctions=(("native modified tRNA", "mim 3-prime adapter I1", "Rnl2tr KQ ligation"),)),
         "The adapter barcode permits samples to be pooled before reverse transcription."),
        ("Reverse-transcribe for 16 h with TGIRT under low-salt conditions",
         molecule("full-length misincorporation-bearing cDNA", (MIM_RT,
                  role("tRNA cDNA retaining modification signatures"))),
         "Conditions favour readthrough rather than enzymatically removing the native modifications."),
        ("Hydrolyse RNA, size-select cDNA and circularize with CircLigase",
         molecule("circular mim-tRNAseq cDNA", (role("circular full-length tRNA cDNA"),)),
         "The RT primer's 5′ RN dinucleotide supports circularization."),
        ("PCR with common P5 and indexed P7 primers",
         molecule("mim-tRNAseq library", (MIM_FORWARD,
                  role("tRNA cDNA with misincorporations"), MIM_REVERSE)),
         "Five to six cycles precede 150-cycle Illumina NextSeq sequencing."),
    ))


def _pandora() -> Workflow:
    source = StructuredRNAInput("15–50-nt cellular small RNA", "heterogeneous blocked termini",
                                "AlkB demethylation", False)
    return chain(source.scene(), (
        ("Repair RNA termini with T4 polynucleotide kinase",
         molecule("end-repaired small RNA", (role("5-prime phosphate"),
                  role("small RNA"), role("3-prime hydroxyl"))),
         "PNK converts otherwise poorly ligatable terminal chemistries to 5′-P/3′-OH."),
        ("Demethylate internal RT-blocking bases with AlkB",
         molecule("PNK- and AlkB-treated small RNA", (role("ligatable demethylated small RNA"),)),
         "The two treatments address independent terminal and internal barriers."),
        ("Ligate the 3′ and 5′ small-RNA adapters sequentially",
         molecule("adapter-ligated small RNA", (role("5-prime SR adapter"),
                  role("repaired small RNA"), PANDORA_RA3),
                  junctions=(("5-prime SR adapter", "repaired small RNA", "RNA ligation"),
                             ("repaired small RNA", "NEBNext 3-prime SR adapter", "RNA ligation"))),
         "The article uses the NEBNext Small RNA Library Prep Set after repair."),
        ("Reverse-transcribe, PCR and PAGE-purify",
         molecule("PANDORA-seq library", (role("P5 small-RNA arm"),
                  role("small-RNA cDNA"), role("indexed P7 arm"))),
         "Libraries were sequenced single-end (50 or 75 cycles depending on sample)."),
    ))


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    summary: str
    notes: str
    endpoint: str
    caveat: str = ""


DEFS = {
    "yamat-seq": Definition("YAMAT-seq — Y-shaped adapter mature-tRNA sequencing",
        "10.1093/nar/gkx005", "Shigematsu et al., Nucleic Acids Research (2017)",
        "A four-way Y adapter recognizes the conserved mature-tRNA NCCA overhang and T4 RNA ligase 2 seals both acceptor-stem nicks in one reaction.",
        "01_yamat_seq.html", "A custom small-RNA sequencing primer reads through the adapter into mature-tRNA cDNA; the experiment used a 100-nt single read."),
    "dm-trna-seq": Definition("DM-tRNA-seq (DM-TGIRT-seq)",
        "10.1038/nmeth.3478", "Zheng et al., Nature Methods (2015)",
        "Wild-type and D135S AlkB remove major Watson–Crick-face methylations before TGIRT template switching, cDNA circularization and indexed PCR.",
        "01_dm_trna_seq.html", "The source gives both library PCR primers and states Illumina HiSeq sequencing, but not a separate sequencing-primer oligo."),
    "mim-trnaseq": Definition("mim-tRNAseq — modification-induced misincorporation tRNA sequencing",
        "10.1016/j.molcel.2021.01.028", "Behrens et al., Molecular Cell (2021)",
        "Barcoded 3′-adapter ligation and prolonged low-salt TGIRT reverse transcription preserve modification-dependent mismatches while producing full-length tRNA cDNA.",
        "01_mim_trnaseq.html", "The indexed library was sequenced for 150 cycles on an Illumina NextSeq; the source does not list a separate custom sequencing primer."),
    "pandora-seq": Definition("PANDORA-seq — pan-small-RNA sequencing",
        "10.1038/s41556-021-00652-7", "Shi et al., Nature Cell Biology (2021)",
        "T4 PNK repairs blocked RNA termini and AlkB removes internal methylation barriers before conventional small-RNA adapter ligation.",
        "01_pandora_seq.html", "NEBNext small-RNA libraries were read single-end for 50 or 75 cycles; kit sequencing-primer sequences are not reproduced in the paper."),
}


def _final(key: str) -> Construct:
    payload = {
        "yamat-seq": (role("Y-5-prime adapter-derived read entry"), role("mature tRNA cDNA"), role("Y-3-prime adapter-derived end")),
        "dm-trna-seq": (role("DM-TGIRT starter-derived entry"), role("tRNA cDNA")),
        "mim-trnaseq": (role("mim adapter-derived entry"), role("native-modification tRNA cDNA")),
        "pandora-seq": (role("small-RNA read entry"), role("repaired small-RNA cDNA")),
    }[key]
    return Construct((role("P5 flow-cell end"), *payload, role("sample index", tag="cbc", length_bp=6,
        molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7 flow-cell end")),
        name=DEFS[key].title)


def _finish_rt_pcr(graph: ProtocolGraph, prefix: str, template, final: Construct,
                   endpoint: str, reviewed: list[str], *, rt_primer_label: str):
    rt_primer = graph.input(prefix + "_rt_primer", rt_primer_label, MaterialKind.OLIGO)
    reverse_transcriptase = graph.input(prefix + "_reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", "structured-RNA-derived cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("structured-RNA-derived cDNA", (
                             role("structured RNA cDNA"),)))
    graph.react(prefix + "_rt", "Reverse-transcribe the adapted RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (cdna,))
    reviewed.extend((prefix + "_rt_primer", prefix + "_reverse_transcriptase",
                     prefix + "_rt"))
    return _finish_pcr(graph, prefix, cdna, final, endpoint, reviewed)


def _finish_pcr(graph: ProtocolGraph, prefix: str, template, final: Construct,
                endpoint: str, reviewed: list[str], *, sequence: bool = True):
    forward = graph.input(prefix + "_pcr_forward", "library PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_pcr_reverse", "indexed library PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", "indexed structured-RNA library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Index and amplify the library", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The source defines the PCR oligo roles but not a distinct sequencing-primer landing geometry on this role-level construct."))
    if sequence:
        graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((
        prefix + "_pcr_forward", prefix + "_pcr_reverse", prefix + "_pcr_polymerase",
        prefix + "_pcr", prefix + "_library",
    ))
    return library


def _yamat_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    trna = graph.input(
        "aminoacylated_trna", "mature aminoacylated tRNA",
        MaterialKind.BIOLOGICAL,
        state=StructuredRNAInput("mature aminoacylated tRNA", "3-prime NCCA end").scene(),
        primary=True)
    alkaline = graph.input("alkaline_condition", "alkaline deacylation condition",
                           MaterialKind.CONDITION)
    deacylated = graph.product(
        "deacylated_trna", "deacylated mature tRNA with free NCCA",
        MaterialKind.MOLECULAR,
        state=molecule("deacylated mature tRNA", (
            role("paired acceptor stem"), role("free 3-prime NCCA"))))
    graph.react("deacylate", "Deacylate mature tRNA at alkaline pH",
                StepKind.CHEMICAL_CONVERSION,
                ((trna, InputRole.TEMPLATE), (alkaline, InputRole.CONDITION)),
                (deacylated,))
    y_adapter = graph.input("y_adapter", "Y-3-prime-AD plus four Y-5-prime adapter variants",
                            MaterialKind.ADAPTER, state=molecule("Y adapter set", (
                                Y5_A, Y5_G, Y5_C, Y5_U, Y3)))
    complex_ = graph.product(
        "y_complex", "annealed Y-adapter–tRNA complex",
        MaterialKind.MOLECULAR,
        state=molecule("Y-adapter–tRNA complex", (Y5_A, role("mature tRNA"), Y3)))
    graph.react("anneal", "Anneal the discriminator-matched Y adapter",
                StepKind.HYBRIDIZATION,
                ((deacylated, InputRole.TEMPLATE), (y_adapter, InputRole.STRUCTURAL)),
                (complex_,))
    rnl2 = graph.input("rnl2", "T4 RNA ligase 2", MaterialKind.CATALYST)
    adapted = graph.product(
        "adapted_trna", "covalently Y-adapted tRNA", MaterialKind.MOLECULAR,
        state=molecule("covalently Y-adapted tRNA", (Y5_A, role("mature tRNA"), Y3),
            junctions=(("Y-5-prime-AD-A", "mature tRNA", "Rnl2 nick ligation"),
                       ("mature tRNA", "Y-3-prime-AD", "Rnl2 nick ligation"))))
    graph.react("seal_nicks", "Seal both Y-adapter nicks", StepKind.NICK_SEALING,
                ((complex_, InputRole.TEMPLATE), (y_adapter, InputRole.STRUCTURAL),
                 (rnl2, InputRole.CATALYST)), (adapted,))
    reviewed = [
        "aminoacylated_trna", "alkaline_condition", "deacylate", "y_adapter",
        "anneal", "rnl2", "seal_nicks",
    ]
    _finish_rt_pcr(graph, "library", adapted, final, d.endpoint, reviewed,
                   rt_primer_label="Y-adapter reverse-transcription primer")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _dm_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    trna = graph.input(
        "trna", "deacylated mature tRNA",
        MaterialKind.BIOLOGICAL,
        state=StructuredRNAInput(
            "deacylated mature tRNA", "free 3-prime CCA",
            "wild-type AlkB plus AlkB-D135S", False).scene(), primary=True)
    alkb = graph.input("alkb", "wild-type AlkB plus AlkB-D135S",
                       MaterialKind.CATALYST)
    demethylated = graph.product(
        "demethylated_trna", "demethylated structured tRNA",
        MaterialKind.MOLECULAR,
        state=molecule("demethylated tRNA", (
            role("structured tRNA with relieved RT blocks"), role("free 3-prime CCA"))))
    graph.react("demethylate", "Demethylate m1A, m3C and m1G",
                StepKind.CHEMICAL_CONVERSION,
                ((trna, InputRole.TEMPLATE), (alkb, InputRole.CATALYST)),
                (demethylated,))
    starter = graph.input("starter_duplex", "blocked DM-TGIRT RNA/DNA starter duplex",
                          MaterialKind.OLIGO, state=molecule("TGIRT starter duplex", (DM_TEMPLATE,)))
    tgirt = graph.input("tgirt", "thermostable group-II-intron reverse transcriptase",
                        MaterialKind.CATALYST)
    cdna = graph.product(
        "full_length_cdna", "full-length tRNA cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("full-length tRNA cDNA", (DM_TEMPLATE, role("tRNA-derived cDNA"))))
    graph.react("template_switch", "Template-switch TGIRT onto the tRNA CCA end",
                StepKind.TEMPLATE_SWITCHING,
                ((demethylated, InputRole.TEMPLATE), (starter, InputRole.STRUCTURAL),
                 (tgirt, InputRole.CATALYST)), (cdna,))
    circligase = graph.input("circligase", "CircLigase II", MaterialKind.CATALYST)
    circle = graph.product(
        "circular_cdna", "circular DM-TGIRT cDNA", MaterialKind.MOLECULAR,
        state=molecule("circular DM-TGIRT cDNA", (role("CircLigase-II-closed cDNA"),)))
    graph.react("circularize", "Circularize the single-stranded cDNA",
                StepKind.CIRCULARIZATION,
                ((cdna, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)),
                (circle,), geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The protocol circularizes single-stranded cDNA with CircLigase II."))
    reviewed = [
        "trna", "alkb", "demethylate", "starter_duplex", "tgirt",
        "template_switch", "circligase", "circularize",
    ]
    _finish_pcr(graph, "library", circle, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _mim_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    trna = graph.input(
        "native_trna", "60–100-nt native mature tRNA",
        MaterialKind.BIOLOGICAL,
        state=StructuredRNAInput("60–100-nt native mature tRNA",
                                 "native 3-prime CCA").scene(), primary=True)
    adapter = graph.input("barcoded_adapter", "preadenylated barcoded 3-prime adapter",
                          MaterialKind.ADAPTER, state=molecule("mim 3-prime adapter", (
                              MIM_ADAPTER_I1,)))
    ligase = graph.input("adapter_ligase", "truncated T4 RNA ligase 2 KQ",
                         MaterialKind.CATALYST)
    adapted = graph.product(
        "adapted_trna", "barcoded 3-prime-adapter-ligated tRNA",
        MaterialKind.MOLECULAR,
        state=molecule("adapter-ligated tRNA", (
            role("native modified tRNA"), MIM_ADAPTER_I1),
            junctions=(("native modified tRNA", "mim 3-prime adapter I1",
                        "Rnl2tr KQ ligation"),)))
    known_end = reported_adapter_end(
        "mim preadenylated 3-prime adapter", kind="single-strand",
        polarity="5-prime", terminal_sequence="G",
        statement="The source publishes the preadenylated, 3-prime-blocked adapter sequence and its RNA-ligation role.")
    graph.react("adapter_ligation", "Ligate a preadenylated barcoded 3-prime adapter",
                StepKind.LIGATION,
                ((trna, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=known_end)
    rt_oligo = graph.input("rt_oligo", "published mim primer-dependent RT oligo",
                           MaterialKind.OLIGO, state=molecule("mim RT oligo", (MIM_RT,)))
    tgirt = graph.input("tgirt", "TGIRT reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product(
        "full_length_cdna", "full-length modification-signature cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("full-length misincorporation-bearing cDNA", (
            MIM_RT, role("tRNA cDNA retaining modification signatures"))))
    graph.react("reverse_transcription", "Reverse-transcribe for 16 hours at low salt",
                StepKind.REVERSE_TRANSCRIPTION,
                ((adapted, InputRole.TEMPLATE), (rt_oligo, InputRole.STRUCTURAL),
                 (tgirt, InputRole.CATALYST)), (cdna,))
    hydrolysis = graph.input("rna_hydrolysis", "RNA hydrolysis condition",
                             MaterialKind.CONDITION)
    clean_cdna = graph.product("clean_cdna", "RNA-free full-length tRNA cDNA",
                               MaterialKind.MOLECULAR, state=cdna.state)
    degraded_rna = graph.product("degraded_rna", "hydrolysed RNA fragments",
                                 MaterialKind.WASTE)
    graph.react("hydrolyse_rna", "Hydrolyse the RNA template",
                StepKind.STRAND_SEPARATION,
                ((cdna, InputRole.TEMPLATE), (hydrolysis, InputRole.CONDITION)),
                (clean_cdna, degraded_rna))
    graph.discard(degraded_rna, "the full-length cDNA is the retained library substrate")
    size_rule = graph.input("cdna_size_selection", "full-length cDNA size-selection rule",
                            MaterialKind.SELECTION)
    sized = graph.product("sized_cdna", "size-selected full-length tRNA cDNA",
                          MaterialKind.MOLECULAR, state=clean_cdna.state)
    rejected = graph.product("short_cdna", "short or incomplete cDNA",
                             MaterialKind.WASTE)
    graph.react("select_cdna", "Size-select full-length cDNA", StepKind.SELECTION,
                ((clean_cdna, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (sized, rejected))
    graph.discard(rejected, "does not span the intended mature-tRNA product")
    circligase = graph.input("circligase", "CircLigase", MaterialKind.CATALYST)
    circle = graph.product(
        "circular_cdna", "circular full-length tRNA cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("circular mim-tRNAseq cDNA", (
            role("circular full-length tRNA cDNA"),)))
    graph.react("circularize", "Circularize full-length cDNA",
                StepKind.CIRCULARIZATION,
                ((sized, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)),
                (circle,), geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The protocol circularizes size-selected full-length cDNA with CircLigase."))
    reviewed = [
        "native_trna", "barcoded_adapter", "adapter_ligase", "adapter_ligation",
        "rt_oligo", "tgirt", "reverse_transcription", "rna_hydrolysis",
        "hydrolyse_rna", "cdna_size_selection", "select_cdna", "circligase",
        "circularize",
    ]
    _finish_pcr(graph, "library", circle, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _pandora_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "small_rna", "15–50-nt cellular small RNA with heterogeneous blocked termini",
        MaterialKind.BIOLOGICAL,
        state=StructuredRNAInput(
            "15–50-nt cellular small RNA", "heterogeneous blocked termini",
            "AlkB demethylation", False).scene(), primary=True)
    pnk = graph.input("pnk", "T4 polynucleotide kinase", MaterialKind.CATALYST)
    repaired = graph.product(
        "repaired_rna", "small RNA with 5-prime phosphate and 3-prime hydroxyl",
        MaterialKind.MOLECULAR,
        state=molecule("end-repaired small RNA", (
            role("5-prime phosphate"), role("small RNA"), role("3-prime hydroxyl"))))
    graph.react("end_repair", "Repair blocked RNA termini",
                StepKind.END_REPAIR,
                ((rna, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (repaired,))
    alkb = graph.input("alkb", "AlkB demethylase", MaterialKind.CATALYST)
    demethylated = graph.product(
        "demethylated_rna", "PNK- and AlkB-treated small RNA",
        MaterialKind.MOLECULAR,
        state=molecule("PNK- and AlkB-treated small RNA", (
            role("ligatable demethylated small RNA"),)))
    graph.react("demethylate", "Demethylate internal RT-blocking bases",
                StepKind.CHEMICAL_CONVERSION,
                ((repaired, InputRole.TEMPLATE), (alkb, InputRole.CATALYST)),
                (demethylated,))
    adapter3 = graph.input("three_adapter", "NEBNext 3-prime small-RNA adapter",
                           MaterialKind.ADAPTER, state=molecule("3-prime adapter", (
                               PANDORA_RA3,)))
    ligase3 = graph.input("three_ligase", "3-prime RNA adapter ligase",
                          MaterialKind.CATALYST)
    ligated3 = graph.product("three_ligated", "3-prime-adapted repaired small RNA",
                             MaterialKind.MOLECULAR, state=demethylated.state)
    graph.react("three_ligation", "Ligate the 3-prime small-RNA adapter",
                StepKind.LIGATION,
                ((demethylated, InputRole.TEMPLATE), (adapter3, InputRole.STRUCTURAL),
                 (ligase3, InputRole.CATALYST)), (ligated3,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source names the NEBNext kit and publishes the adapter sequence, but the page model does not reconstruct its complete modified ligation substrate."))
    adapter5 = graph.input("five_adapter", "NEBNext 5-prime small-RNA adapter",
                           MaterialKind.ADAPTER)
    ligase5 = graph.input("five_ligase", "5-prime RNA adapter ligase",
                          MaterialKind.CATALYST)
    ligated5 = graph.product(
        "adapted_rna", "two-adapter PANDORA small RNA",
        MaterialKind.MOLECULAR,
        state=molecule("adapter-ligated small RNA", (
            role("5-prime SR adapter"), role("repaired small RNA"), PANDORA_RA3),
            junctions=(("5-prime SR adapter", "repaired small RNA", "RNA ligation"),
                       ("repaired small RNA", "NEBNext 3-prime SR adapter", "RNA ligation"))))
    graph.react("five_ligation", "Ligate the 5-prime small-RNA adapter",
                StepKind.LIGATION,
                ((ligated3, InputRole.TEMPLATE), (adapter5, InputRole.STRUCTURAL),
                 (ligase5, InputRole.CATALYST)), (ligated5,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining article uses the NEBNext Small RNA kit; complete modified 5-prime adapter geometry is not reproduced."))
    reviewed = [
        "small_rna", "pnk", "end_repair", "alkb", "demethylate",
        "three_adapter", "three_ligase", "three_ligation", "five_adapter",
        "five_ligase", "five_ligation",
    ]
    rt_primer = graph.input("library_rt_primer", "NEBNext small-RNA reverse-transcription primer",
                            MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("library_reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "PANDORA small-RNA cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("PANDORA small-RNA cDNA", (
                             role("repaired small-RNA cDNA"),)))
    graph.react("library_rt", "Reverse-transcribe the adapted small RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated5, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (cdna,))
    reviewed.extend(("library_rt_primer", "library_reverse_transcriptase", "library_rt"))
    library = _finish_pcr(
        graph, "library", cdna, final, d.endpoint, reviewed, sequence=False)
    page = graph.input("page_selection", "PAGE library-size selection rule",
                       MaterialKind.SELECTION)
    selected_library = graph.product("selected_library", "PAGE-purified PANDORA-seq library",
                                     MaterialKind.LIBRARY, state=final)
    rejected_library = graph.product("rejected_library", "off-size amplified products",
                                     MaterialKind.WASTE)
    graph.react("page_select", "PAGE-purify the amplified library", StepKind.SELECTION,
                ((library, InputRole.TEMPLATE), (page, InputRole.SELECTION)),
                (selected_library, rejected_library))
    graph.discard(rejected_library, "falls outside the selected small-RNA library size")
    graph.sequence(selected_library, SequencingPlan(limitation=d.endpoint))
    reviewed.extend(("page_selection", "page_select", "selected_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def get(key: str):
    d = DEFS[key]
    workflows = {"yamat-seq": _yamat, "dm-trna-seq": _dm,
                 "mim-trnaseq": _mim, "pandora-seq": _pandora}
    oligo_map = {
        "yamat-seq": (
            "Y-3′-AD — 5′-/5Phos/GTATCCAGTTGGAATTCTCGGGTGCCAAGG/3ddC/-3′",
            "Y-5′-AD-A — 5′-GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGga-3′",
            "Y-5′-AD-G — 5′-GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGgg-3′",
            "Y-5′-AD-C — 5′-GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGgc-3′",
            "Y-5′-AD-U — 5′-GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTGgu-3′",
        ),
        "dm-trna-seq": (
            "TGIRT RNA template — 5′-AGAUCGGAAGAGCACACGUCUAGUUCUACAGUCCGACGAUC/3SpC3/-3′",
            "multiplex PCR primer — 5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC-3′",
            "indexed PCR primer — 5′-CAAGCAGAAGACGGCATACGAGAT[BARCODE]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-3′",
        ),
        "mim-trnaseq": (
            "3′ adapter I1 — 5′-/5Phos/GATATCGTCAAGATCGGAAGAGCACACGTCTGAA/ddC/-3′",
            "primer-dependent RT oligo — 5′-/5Phos/rNrNAGATCGGAAGAGCGTCGTGTAGGGAAAGAG/iSp18/GTGACTGGAGTTCAGACGTGTGCTC-3′",
            "library forward — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT*C-3′",
            "indexed reverse — 5′-CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGT*G-3′",
        ),
        "pandora-seq": (
            "NEBNext 3′ SR adapter — 5′-/rApp/AGATCGGAAGAGCACACGTCT/3AmMO/-3′",
        ),
    }
    final = _final(key)
    spec = SimpleNamespace(
        TITLE=d.title, SOURCE=f'{d.citation}; <a href="https://doi.org/{d.doi}">doi:{d.doi}</a>.',
        SUMMARY=d.summary, NOTES=d.notes, CAVEAT=d.caveat,
        oligos=lambda: oligo_map[key], workflow=workflows[key], FINAL_LIBRARY=final,
        FINAL_CAPTION="Role-level final library. Exact source-published oligos are listed above; undisclosed kit junctions remain role tokens.",
        SEQ_PRIMERS=(), READ_LENGTHS=None,
        SEQUENCING_INTRO="", SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE="The defining source does not publish a distinct sequencing-primer oligo or enough kit sequence to compute its exact landing site.")
    if key == "yamat-seq":
        spec.protocol_graph = lambda: _yamat_graph(d, final)
    elif key == "dm-trna-seq":
        spec.protocol_graph = lambda: _dm_graph(d, final)
    elif key == "mim-trnaseq":
        spec.protocol_graph = lambda: _mim_graph(d, final)
    else:
        spec.protocol_graph = lambda: _pandora_graph(d, final)
    return spec
