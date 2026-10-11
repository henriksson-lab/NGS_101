"""Shared chemistry for structured and heavily modified RNA libraries.

The important invariant here is that end repair, demethylation and reverse-transcriptase
choice are separate operations.  A caller cannot describe AlkB treatment as generic
"library preparation", or silently turn a modification-preserving method into a
demethylating one.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import Construct, Scene, Segment, Workflow, feature
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
    return SimpleNamespace(
        TITLE=d.title, SOURCE=f'{d.citation}; <a href="https://doi.org/{d.doi}">doi:{d.doi}</a>.',
        SUMMARY=d.summary, NOTES=d.notes, CAVEAT=d.caveat,
        oligos=lambda: oligo_map[key], workflow=workflows[key], FINAL_LIBRARY=_final(key),
        FINAL_CAPTION="Role-level final library. Exact source-published oligos are listed above; undisclosed kit junctions remain role tokens.",
        SEQ_PRIMERS=(), READ_LENGTHS=None,
        SEQUENCING_INTRO="", SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE="The defining source does not publish a distinct sequencing-primer oligo or enough kit sequence to compute its exact landing site.")
