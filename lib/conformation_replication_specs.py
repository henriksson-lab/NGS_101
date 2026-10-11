"""Source-backed page specifications for targeted 3C and Repli-seq methods."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import seqprimers as sp
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature
from chromosome_capture import (Circular4CTemplate, ThreeCLibrary,
                                inverse_pcr_product, multiway_contact, region)
from replication_labeling import CopyNumberReplication, NascentPulse
from protocol_model import da_tailed


def role(name: str, *, tag=None, molecular_feature=None) -> Segment:
    return Segment(name, f"[{name}]", tag, placeholder=True,
                   feature=molecular_feature)


def molecule(label: str, names: tuple[str, ...], *,
             junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    parts = [region(name) for name in names]
    sc = Scene(); sc.strand("molecule", parts, label=label)
    for left, right, text in junctions:
        sc.junction("molecule", left, right, text)
    sc.labels("molecule")
    return sc


def chain(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        w.react(action, state.rows(), name=action, note=note)
    return w


def _three_c(enzyme: str, names: tuple[str, ...]) -> ThreeCLibrary:
    return ThreeCLibrary(enzyme, names)


def workflow_for(key: str) -> Workflow:
    if key == "4c-seq":
        contact = _three_c("DpnII", ("viewpoint fragment", "contacted reporter fragment"))
        circle = Circular4CTemplate(contact, "Csp6I", "viewpoint fragment")
        return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "proximal chromatin")), (
            ("Digest chromatin with the primary restriction enzyme and proximity-ligate",
             contact.scene(), "The first junction records proximity in fixed nuclei."),
            ("Reverse crosslinks, digest with a second enzyme and ligate dilute DNA",
             circle.scene(), "The second ligation recircularizes shorter 4C templates."),
            ("Inverse-PCR outward from the chosen viewpoint", inverse_pcr_product(circle),
             "A two-stage PCR adds sample and Illumina sequences to viewpoint-specific products."),
        ))
    if key == "ng-capture-c":
        contact = _three_c("DpnII or NlaIII", ("captured viewpoint fragment", "reporter fragment"))
        adapted = molecule("indexed 3C sequencing library",
                           ("left Illumina adapter", "viewpoint", "contact junction", "reporter", "right Illumina adapter"),
                           junctions=(("left Illumina adapter", "viewpoint", "adapter ligation"),
                                      ("viewpoint", "contact junction", "proximity ligation"),
                                      ("reporter", "right Illumina adapter", "adapter ligation")))
        return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "proximal chromatin")), (
            ("Digest and proximity-ligate to make a 3C library", contact.scene(),
             "No biotin fill-in is used to select the contact junction."),
            ("Sonicate to about 200 bp, end-repair and add indexed Illumina adapters", adapted,
             "Random sonication ends distinguish original molecules from PCR duplicates."),
            ("Perform two sequential captures with biotinylated viewpoint oligos", adapted,
             "Streptavidin enrichment retains fragments containing the selected viewpoint."),
        ))
    if key == "umi-4c":
        contact = _three_c("DpnII", ("viewpoint fragment", "reporter fragment"))
        one_end = molecule("single-end-adapter 3C fragment",
                           ("viewpoint", "contact junction", "reporter", "random sonication endpoint",
                            "single-end Illumina adapter"),
                           junctions=(("viewpoint", "contact junction", "proximity ligation"),
                                      ("random sonication endpoint", "single-end Illumina adapter", "adapter ligation")))
        return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "proximal chromatin")), (
            ("Digest with DpnII, proximity-ligate and reverse crosslinks", contact.scene(),
             "This is a linear 3C library; UMI-4C does not perform 4C recircularization."),
            ("Sonicate, end-repair and add one 3-prime dA to the ligatable end",
             da_tailed("dA-tailed UMI-4C fragment", "viewpoint-contact fragment with random sonication endpoint"),
             "The random shearing coordinate is already fixed before adapter ligation."),
            ("Ligate the single-end adapter", one_end,
             "The random shearing coordinate acts as the molecular identifier; it is not a synthesized random barcode."),
            ("Amplify with two nested viewpoint PCRs", one_end,
             "The nested viewpoint primer and universal adapter primer create the sequencing library."),
        ))
    if key == "tri-c":
        contact = multiway_contact("NlaIII", "viewpoint fragment",
                                   ("reporter fragment 1", "reporter fragment 2"))
        adapted = molecule("adapter-ligated multiway 3C fragment",
                           ("left TruSeq adapter", "random sonication endpoint", "viewpoint",
                            "reporter 1", "reporter 2", "right TruSeq adapter"),
                           junctions=(("left TruSeq adapter", "random sonication endpoint", "adapter ligation"),
                                      ("viewpoint", "reporter 1", "proximity ligation"),
                                      ("reporter 1", "reporter 2", "proximity ligation"),
                                      ("reporter 2", "right TruSeq adapter", "adapter ligation")))
        return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "two proximal regions")), (
            ("Digest with NlaIII and proximity-ligate in nuclei", contact.scene("multiway 3C concatemer"),
             "At least two ligation junctions in one molecule encode a concurrent contact."),
            ("Reverse crosslinks and sonicate gently to about 450 bp", adapted,
             "Longer fragments preserve multiple short NlaIII products and their junctions."),
            ("Add TruSeq adapters, index, then perform double viewpoint capture", adapted,
             "A biotinylated oligo targets the middle of the short viewpoint restriction fragment."),
        ))
    if key == "el-repli-seq":
        pulse = NascentPulse("BrdU", ("early S phase", "late S phase"))
        adapted = molecule("adapter-ligated nascent DNA",
                           ("left Illumina adapter", "BrdU-containing nascent DNA", "right Illumina adapter"),
                           junctions=(("left Illumina adapter", "BrdU-containing nascent DNA", "adapter ligation"),
                                      ("BrdU-containing nascent DNA", "right Illumina adapter", "adapter ligation")))
        return chain(molecule("asynchronously cycling cells", ("replicating chromosomes",)), (
            ("Pulse-label nascent DNA with BrdU", pulse.labelled_scene(),
             "Only DNA synthesized during the two-hour pulse contains BrdU."),
            ("Fix and sort early- and late-S fractions by DNA content", pulse.sorted_scene(),
             "Early and late fractions remain separate libraries."),
            ("Purify, sonicate, end-repair and ligate universal Illumina adapters", adapted,
             "Adapters are installed before immunoprecipitation."),
            ("Denature, immunoprecipitate BrdU DNA, then index and PCR-amplify", adapted,
             "Sequencing coverage is compared as an early/late ratio."),
        ))
    if key == "screpli-seq":
        dosage = CopyNumberReplication("mid-S-phase", "SeqPlex DOP-PCR whole-genome amplification")
        amplified = molecule("single-cell WGA product",
                             ("left SEQXE amplification handle", "cell genome copy",
                              "right SEQXE amplification handle"))
        adapted = molecule("single-cell Illumina library",
                           ("P5 / sample index", "WGA genomic fragment", "sample index / P7"),
                           junctions=(("P5 / sample index", "WGA genomic fragment", "adapter ligation"),
                                      ("WGA genomic fragment", "sample index / P7", "adapter ligation")))
        return chain(molecule("one flow-sorted S-phase cell", ("partly replicated genome",)), (
            ("Use copy number as the replication mark", dosage.states(),
             "Replicated domains have twice the DNA dosage of unreplicated domains in the same cell."),
            ("Whole-genome amplify each cell by SeqPlex DOP-PCR", amplified,
             "G1 cells processed in parallel provide the amplification and mappability reference."),
            ("Fragment WGA DNA, end-repair and add one 3-prime dA to each end",
             da_tailed("dA-tailed single-cell WGA fragment", "WGA genomic fragment"),
             "The amplified genomic fragment is shown before adapter ligation."),
            ("Ligate Illumina adapters", adapted,
             "The protocol permits KAPA LTP or NEBNext Ultra II library construction."),
            ("PCR-index each single-cell library and sequence", adapted,
             "The index assigns reads to one cell; it is not an inline cell barcode."),
        ))
    raise KeyError(key)


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    summary: str
    final_parts: tuple[Segment, ...]
    sequencing: str
    standard_truseq: bool
    caveat: str = ""


ENDPOINT = Segment("random sonication endpoint identifier",
                   "[random shearing endpoint used as UMI]", "umi", placeholder=True,
                   feature=feature("sonication_endpoint", "umi", "random",
                                   note="physical shearing coordinate, not a synthesized barcode"))
SAMPLE_INDEX = Segment("PCR sample index", "I" * 8, "cbc", placeholder=True,
                       feature=feature("sample_index", "sample_index", "unknown"))

DEFS = {
    "4c-seq": Definition(
        "4C-seq — circular chromosome conformation capture sequencing",
        "10.1016/j.ymeth.2019.07.014", "Krijger et al., Methods (2020)",
        "Two restriction digests and two ligations convert a viewpoint-centred 3C contact into a circle that can be opened by outward-facing inverse PCR.",
        (role("viewpoint PCR / Illumina end"), SAMPLE_INDEX, role("viewpoint fragment"),
         role("captured reporter"), role("second viewpoint PCR / Illumina end")),
        "A viewpoint-specific reading primer enters the captured reporter beside the primary restriction junction; sample and flow-cell sequences are added by the two PCR stages.", False,
        "Viewpoint primers are locus-specific. Their exact bases cannot be generalized, so the page shows their molecular roles rather than an invented universal oligo."),
    "ng-capture-c": Definition(
        "NG Capture-C — next-generation Capture-C", "10.1038/nmeth.3664",
        "Davies et al., Nature Methods (2016)",
        "A sonicated, indexed 3C library undergoes two rounds of capture with biotinylated DNA oligos directed to selected restriction fragments.",
        (role("captured viewpoint"), role("proximity-ligation junction"), role("reporter fragment")),
        "Paired-end Illumina reads recover the viewpoint-containing end and its ligated reporter; random sonication coordinates identify PCR duplicates.", True),
    "umi-4c": Definition(
        "UMI-4C — quantitative targeted chromosome contacts", "10.1038/nmeth.3922",
        "Schwartzman et al., Nature Methods (2016)",
        "Random shearing supplies molecule-specific endpoints before single-end adapter ligation and nested viewpoint PCR, enabling quantitative targeted 3C profiling.",
        (role("nested viewpoint / Illumina end"), SAMPLE_INDEX, role("viewpoint"),
         role("contact reporter"), ENDPOINT, role("universal adapter end")),
        "Paired-end sequencing recovers the viewpoint/contact molecule and the random shearing endpoint used to collapse PCR duplicates.", False,
        "The published ‘UMI’ is the genomic sonication endpoint, not a fixed-length random oligonucleotide; the exact locus-specific bait primers vary by experiment."),
    "tri-c": Definition(
        "Tri-C — targeted multiway chromatin contacts", "10.1038/s41588-018-0253-2",
        "Oudelaar et al., Nature Genetics (2018)",
        "Small NlaIII fragments, gentle sonication and viewpoint capture preserve multiple proximity-ligation junctions in one sequenced molecule.",
        (role("viewpoint"), role("reporter 1"), role("second proximity junction"),
         role("reporter 2"), ENDPOINT),
        "150-bp paired-end Illumina reads identify the viewpoint and two or more reporters within one roughly 450-bp captured 3C fragment.", True),
    "el-repli-seq": Definition(
        "E/L Repli-seq — early/late replication timing", "10.1038/nprot.2017.148",
        "Marchal et al., Nature Protocols (2018)",
        "BrdU-labelled nascent DNA is separated into early- and late-S populations, adapter-ligated, immunoprecipitated and sequenced as two quantitative fractions.",
        (role("BrdU-enriched genomic fragment"),),
        "Ordinary Illumina genomic reads measure coverage in separate early- and late-S libraries; replication timing is their normalized log2 ratio.", True),
    "screpli-seq": Definition(
        "scRepli-seq — single-cell DNA replication sequencing", "10.1038/s41596-020-0378-5",
        "Miura et al., Nature Protocols (2020)",
        "A single S-phase cell is whole-genome amplified and sequenced; twofold copy-number differences distinguish replicated from unreplicated domains.",
        (role("SeqPlex-amplified genomic fragment"),),
        "Each PCR index identifies one single-cell library. Genomic Read 1/Read 2 coverage, normalized to G1 controls, reports replicated and unreplicated domains.", True),
}


def get(key: str):
    d = DEFS[key]
    if d.standard_truseq:
        # These defining protocols predate routine unique-dual indexing and describe one
        # PCR sample index.  Do not silently modernize their read structure.
        final, primers = truseq_library(list(d.final_parts), d.title + " final library",
                                        dual_index=False)
    else:
        final = Construct([role("Illumina-compatible left end"), *d.final_parts,
                           role("Illumina-compatible right end")], name=d.title + " final library")
        primers = ()
    return SimpleNamespace(
        TITLE=d.title, NOTES="01_" + key.replace("-", "_") + ".html",
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final, SEQ_PRIMERS=primers,
        FINAL_CAPTION=d.sequencing, SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.sequencing,
        SEQUENCING_UNAVAILABLE=("The defining protocol uses locus-specific or custom primer sequence that cannot be represented by one universal binding sequence."
                                if not primers else ""),
        workflow=lambda: workflow_for(key), oligos=lambda: (),
    )
