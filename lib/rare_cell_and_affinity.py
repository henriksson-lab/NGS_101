"""Shared molecular models for rare-cell RNA capture, BAGs and CXXC CAP-seq.

The methods collected here share a useful design constraint: the named protocol is
primarily a *selection or partitioning operation*.  Library shells are therefore built
from canonical Illumina constructors where the source names a standard chemistry, and
otherwise remain explicit role tokens rather than guessed bases.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import illumina as il
import rt
import seqprimers as sp
from batch_ngs import nextera_library, seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature
from protocol_model import role, role_library


def _scene(label: str, *parts: Segment) -> Scene:
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    sc.labels("molecule")
    return sc


def _workflow(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        w.react(action, state.rows(), name=action, note=note)
    return w


TRANSCRIPT = role("transcript-derived insert")
GENOMIC = role("genomic insert")
RAISIN_CST = (("CHAPS", "0.49% w/v"), ("NaCl", "146 mM"),
              ("CaCl2", "1 mM"), ("MgCl2", "21 mM"),
              ("Tris pH 8.0", "10 mM"))
CAP_HUMAN_WASH_MM = 600
CAP_MOUSE_WASH_MM = 560


def raisin_particle() -> Scene:
    sc = Scene()
    sc.strand("nucleus", [role("intact nucleus")], label="RAISIN")
    sc.strand("ribosome", [role("rough ER / attached ribosomes"),
                            role("ribosome-bound poly(A) RNA")], label="retained exterior RNA")
    sc.note("ribosome", "CHAPS/salt extraction retains nucleus-associated ribosomes and mRNA")
    return sc


def smart_cdna_scene() -> Scene:
    sc = Scene()
    sc.strand("cDNA", [role("oligo-dT primer / SMART handle"), TRANSCRIPT,
                        role("copied locked TSO")], label="full-length cDNA")
    sc.note("cDNA", "poly(dT) reverse transcription + template switching; 21-cycle PCR")
    return sc


def raisin_workflow() -> Workflow:
    final, _ = raisin_library()
    return _workflow(_scene("frozen tissue", role("frozen tissue")), (
        ("Extract with CST (RAISIN) buffer", raisin_particle(),
         "CHAPS with salt, Ca²⁺ and Mg²⁺ preserves intact nuclei with attached ribosomes."),
        ("Sort one RAISIN per well, lyse and purify RNA", _scene(
            "single-RAISIN RNA", role("nuclear RNA"), role("ribosome-bound poly(A) RNA")),
         "The plate implementation sorts into TCL/β-mercaptoethanol and uses RNAClean XP."),
        ("Oligo-dT reverse-transcribe with a locked TSO, then amplify", smart_cdna_scene(),
         "Maxima reverse transcriptase is followed by 21 cycles of KAPA HiFi PCR."),
        ("Tagment amplified cDNA with Nextera XT and add sample indexes",
         Scene.duplex(list(final), label="RAISIN RNA-seq library"),
         "The custom barcode-adapter sequences were available only on request."),
    ))


def raisin_library() -> tuple[Construct, tuple]:
    return nextera_library([seg("RAISIN transcript-derived insert", "X" * 36,
                                placeholder=True)], "RAISIN RNA-seq library")


CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x GemCode v2")
UMI = feature("umi", "umi", "random")
SAMPLE_INDEX = feature("sample_index_i7", "sample_index", "unknown")


def gemcode_v2_library() -> tuple[Construct, tuple]:
    """The standard 10x v2 endpoint used by droplet RAISIN and MIRACL-seq."""
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        Segment("cell barcode", "B" * 16, "cbc", placeholder=True,
                feature=CELL_BARCODE),
        Segment("UMI", "U" * 10, "umi", placeholder=True, feature=UMI),
        Segment("poly(dT)30", "T" * 30),
        Segment("V", "V", placeholder=True), Segment("N", "N", placeholder=True),
        Segment("RAISIN transcript-derived insert", "X" * 36, placeholder=True),
        Segment("dA junction", "A"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
                feature=SAMPLE_INDEX),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="GemCode v2 RAISIN library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid GemCode v2 endpoint: " + "; ".join(problems))
    return lib, primers


@dataclass(frozen=True)
class MiraclLoading:
    coefficient: int

    def __post_init__(self) -> None:
        if self.coefficient < 2:
            raise ValueError("MIRACL-seq is defined by deliberate droplet overloading")

    def partition(self) -> Scene:
        return _scene(f"{self.coefficient}× overloaded GEM",
                      role("one 10x gel bead"), role("multiple RAISIN nuclei"),
                      role("oil-emulsion partition"))


def miracl_workflow(coefficient: int = 2) -> Workflow:
    loading = MiraclLoading(coefficient)
    final, _ = gemcode_v2_library()
    barcoded = _scene("barcoded first strand", role("10x Read 1 handle"),
                      role("16-nt cell barcode", molecular_feature=CELL_BARCODE),
                      role("10-nt UMI", molecular_feature=UMI), TRANSCRIPT)
    return _workflow(raisin_particle(), (
        (f"Overload the GemCode v2 channel ({coefficient}× shown)", loading.partition(),
         "The defining experiments tested 2×–8× loading; rare-cell homotypic multiplets remain unlikely."),
        ("Lyse in GEMs and perform barcoded reverse transcription", barcoded,
         "One bead barcode labels cDNA from every nucleus sharing that partition."),
        ("Amplify, shear, attach the 5′ adapter and sample index",
         Scene.duplex(list(final), label="MIRACL-seq library"),
         "This is the standard GemCode v2 molecular endpoint; heterotypic multiplets are removed computationally."),
    ))


# Verbatim capture oligos printed for the BAG platform.  The hybrid paper uses both.
BAG_DNA_CAPTURE_WRITTEN = "/5ACryd//iSp18/TGTGTTGGGTGTGTTTGGKKKKKKKGKKKKKKKKNN"
BAG_RNA_CAPTURE_WRITTEN = "/5ACryd//iSp18/AAGCAGTGGTATCAACGCAGAGTNNWNNNSTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

BAG_BARCODE = feature("bag_barcode", "cell_barcode", "combinatorial",
                      group="bag_barcode")
BAG_VARIETAL = feature("varietal_tag", "umi", "combinatorial",
                       group="varietal_tag")


def bag_final_library() -> Construct:
    return role_library("hybrid BAG-seq library", (
        role("Read 1 common sequence"),
        role("BAG barcode round 3", length_bp=6, molecular_feature=BAG_BARCODE),
        role("BAG barcode round 2", length_bp=6, molecular_feature=BAG_BARCODE),
        role("varietal tag round 2", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("BAG barcode round 1", length_bp=5, molecular_feature=BAG_BARCODE),
        role("varietal tag round 1", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("NlaIII-proximal genomic tag", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("CATG NlaIII site", length_bp=4),
        role("DNA- or RNA-derived insert"),
    ))


def bag_workflow() -> Workflow:
    captured = _scene("hybrid BAG", role("polyacrylamide gel matrix"),
                      role("5′-Acrydite random-TG DNA primer"),
                      role("5′-Acrydite poly(dT) RNA primer"),
                      role("single nucleus"))
    copies = _scene("gel-bound copies", role("Acrydite anchor"), GENOMIC,
                    role("Acrydite anchor"), TRANSCRIPT, role("TSO-derived NlaIII site"))
    cut = _scene("NlaIII-cleaved copies", role("Acrydite anchor"),
                 role("DNA- or RNA-derived copy"), role("CATG sticky end"))
    tagged = _scene("three-round tagged molecules", role("UP1"),
                    role("BAG barcode 1 + varietal tag 1"),
                    role("BAG barcode 2 + varietal tag 2"),
                    role("BAG barcode 3"), role("DNA- or RNA-derived copy"))
    final = bag_final_library()
    return _workflow(_scene("single nucleus", role("single sorted nucleus")), (
        ("Co-encapsulate at 85 °C and polymerize the hybrid BAG", captured,
         "Both 5′-Acrydite primers copolymerize into the gel; the lower denaturation temperature preserves RNA."),
        ("Reverse-transcribe RNA and Klenow-extend genomic DNA", copies,
         "Both copy types remain covalently gel-bound; the TSO supplies a double-stranded NlaIII site to cDNA."),
        ("Digest the bound copies with NlaIII", cut,
         "Cleavage exposes the sticky end used by the first split-pool extension."),
        ("Add two ligation/hybridization-extension tags, then a split PCR tag", tagged,
         "The three rounds form a 17-nt BAG barcode and a 12-nt varietal tag."),
        ("Pool and tagment the amplicons", _scene("hybrid BAG-seq library", *tuple(final)),
         "Paired-end 150-bp sequencing starts one read at the tag-bearing end and the other at a random tagmentation end."),
    ))


def cap_library() -> Construct:
    return role_library("CAP-seq library", (role("unmethylated CpG-rich insert"),))


def cap_workflow() -> Workflow:
    adapted = _scene("adapter-ligated genomic library", role("Solexa left adapter"),
                     role("genomic insert with CpG sites"), role("Solexa right adapter"))
    bound = _scene("CXXC-bound library", role("nickel-Sepharose-bound recombinant CXXC"),
                   role("adapter-ligated unmethylated CpG-rich fragment"))
    final = cap_library()
    return _workflow(adapted, (
        ("Bind adapter-ligated input to recombinant CXXC on nickel Sepharose", bound,
         "CXXC affinity selects clustered unmethylated CpGs; the adapters were installed before purification."),
        (f"Wash with {CAP_HUMAN_WASH_MM} mM NaCl (human) or "
         f"{CAP_MOUSE_WASH_MM} mM NaCl (mouse)", bound,
         "The lower-stringency mouse wash retains its relatively CpG-poor islands."),
        ("Elute at 1 M NaCl, pool and sequence",
         _scene("CAP-seq library", *tuple(final)),
         "The defining experiment used single-end Solexa sequencing."),
    ))


def get(key: str):
    if key == "raisin-rna-seq":
        final, primers = raisin_library()
        return SimpleNamespace(
            TITLE="RAISIN RNA-seq — intact nuclei with ribosome-bound RNA",
            NOTES="01_raisin_rna_seq.html",
            SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.08.003">Drokhlyansky et al., Cell (2020)</a>.',
            SUMMARY="CST extraction preserves an intact nucleus together with attached ribosomes and their poly(A) RNA; the plate implementation uses SMART-seq2 followed by Nextera XT.",
            CAVEAT="The paper does not print the locked TSO or custom barcode-adapter sequences; those workflow states are labelled by role rather than assigned bases.",
            FINAL_LIBRARY=final, SEQ_PRIMERS=primers, READ_LENGTHS=None,
            FINAL_CAPTION="A Nextera XT library made from plate-amplified RAISIN cDNA.",
            SEQUENCING_INTRO="Standard Nextera primers bind the completed library as computed below.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="",
            workflow=raisin_workflow, oligos=lambda: ())
    if key == "miracl-seq":
        final, primers = gemcode_v2_library()
        return SimpleNamespace(
            TITLE="MIRACL-seq — droplet overloading for rare-cell RNA-seq",
            NOTES="01_miracl_seq.html",
            SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.08.003">Drokhlyansky et al., Cell (2020)</a>.',
            SUMMARY="MIRACL-seq deliberately overloads GemCode v2 droplets with unsorted RAISIN nuclei, increasing rare-cell recovery while computationally removing heterotypic multiplets.",
            CAVEAT="MIRACL-seq changes partition loading, not the GemCode v2 oligo chemistry; the final construct therefore uses the standard v2 molecular endpoint.",
            FINAL_LIBRARY=final, SEQ_PRIMERS=primers, READ_LENGTHS={"Read 1": 26},
            FINAL_CAPTION="GemCode v2 3′ gene-expression library from overloaded RAISIN partitions.",
            SEQUENCING_INTRO="Read 1 records the 16-nt bead barcode and 10-nt UMI; Read 2 interrogates the transcript-derived insert.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="",
            workflow=miracl_workflow, oligos=lambda: ())
    if key == "hybrid-bag-seq":
        return SimpleNamespace(
            TITLE="Hybrid BAG-seq — joint single-nucleus DNA and RNA sequencing",
            NOTES="01_hybrid_bag_seq.html",
            SOURCE='Defining source: <a href="https://doi.org/10.1186/s13059-025-03790-5">Li et al., Genome Biology (2025)</a>; capture-primer sequences were checked against the original <a href="https://doi.org/10.1101/gr.253047.119">BAG-seq paper</a>.',
            SUMMARY="Acrydite random-TG and poly(dT) primers copolymerized into one gel capture genomic DNA and RNA-derived cDNA from the same nucleus before three rounds of split-pool tagging.",
            CAVEAT="The papers establish the barcode geometry and print the capture primers, but do not expose one complete base-level final adapter construct; unreleased regions are shown by role.",
            FINAL_LIBRARY=bag_final_library(), SEQ_PRIMERS=(), READ_LENGTHS=None,
            FINAL_CAPTION="Role-level paired-end library; Read 1 begins with the split-pool tags and Read 2 begins at a random tagmentation end.",
            SEQUENCING_INTRO="", SEQUENCING_ENDING="NovaSeq 6000 paired-end 150-bp reads were used. Read 1 begins with an 80-base tag/adapter block containing the 17-nt BAG barcode and 12-nt varietal tag.",
            SEQUENCING_UNAVAILABLE="The custom final-library sequencing-primer sites are not printed in the authoritative sources.",
            workflow=bag_workflow,
            oligos=lambda: (f"5′-{BAG_DNA_CAPTURE_WRITTEN}-3′ — Acrydite random-TG DNA capture primer",
                             f"5′-{BAG_RNA_CAPTURE_WRITTEN}-3′ — Acrydite poly(dT) RNA capture primer"))
    if key == "cap-seq":
        return SimpleNamespace(
            TITLE="CAP-seq — CXXC affinity purification sequencing",
            NOTES="01_cap_seq.html",
            SOURCE='Defining source: <a href="https://doi.org/10.1371/journal.pgen.1001134">Illingworth et al., PLOS Genetics (2010)</a>.',
            SUMMARY="Adapter-ligated genomic DNA is selected on immobilized recombinant CXXC for clusters of unmethylated CpGs, then eluted and sequenced.",
            CAVEAT="The defining paper refers to earlier Solexa library methods without printing the adapter or sequencing-primer bases; the library ends are therefore shown by role.",
            FINAL_LIBRARY=cap_library(), SEQ_PRIMERS=(), READ_LENGTHS=None,
            FINAL_CAPTION="Solexa library enriched for fragments carrying clustered unmethylated CpGs.",
            SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining study used single-end Solexa reads to map CXXC-enriched fragments.",
            SEQUENCING_UNAVAILABLE="Exact sequencing-primer sites cannot be drawn because the cited library protocol does not print the adapter sequences.",
            workflow=cap_workflow, oligos=lambda: ())
    raise KeyError(key)
