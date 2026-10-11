"""Direct double-strand-break capture chemistries.

The shared model keeps the two molecular ends distinct: the proximal end is the
pre-existing cellular break labelled in situ, whereas the distal end is created later
by fragmentation.  A constructor cannot silently exchange those origins.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import illumina as il
import seqprimers as sp
from batch_ngs import seg, truseq_library
from chemdraw import Construct, Segment, feature
from protocol_model import chain, da_tailed, molecule, role


@dataclass(frozen=True)
class BreakCaptureDesign:
    name: str
    proximal_adapter: str
    distal_adapter: str
    immobilization: str
    direct_break_read: bool

    def __post_init__(self):
        if not all((self.name, self.proximal_adapter, self.distal_adapter,
                    self.immobilization)):
            raise ValueError("break capture requires named proximal/distal ends and immobilization")
        if self.proximal_adapter == self.distal_adapter:
            raise ValueError("the native-break and fragmentation ends must remain distinguishable")


# Crosetto et al. Supplementary Table 5; the hairpins are written continuously 5'->3'.
BLESS_PROXIMAL = "TACTACCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGGTAGTA"
BLESS_DISTAL = "CGTCGTCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGACGACG"
BLESS_PCR_PROXIMAL = "CCCTAGCGTAACTCTCGAGGTAGTA"
BLESS_PCR_DISTAL = "CTAGCGTAACTCTCGAGACGACG"
BLESS_PROXIMAL_READ_TAG = "TCGAGGTAGTA"
BLESS_DISTAL_READ_TAG = "TCGAGACGACG"

# Lensing et al. Fig. 1 / Supplementary Fig. 1c.
DSB_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
DSB_P7 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
DSB_PCR_F = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA"
DSB_PCR_R_AD06 = "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCAGACGTGT"

# Cameron et al. SITE-Seq supplementary protocol / Supplementary Table 26.
# Modifications are kept in the ordering strings returned by ``get``; these constants
# are the molecular bases used to derive the final library.
SITE_A1_F = "GTTGACATGCTGGATTGAGACTTCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
SITE_A1_R = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGAAGTCTCAATCCAGCATGTCAAC"
SITE_A2_R = "ACTACTAATACGAC"
SITE_A2_F7 = "GTCGTATTAGTAGTANNNNNNNAGATCGGAAGAGCACACGTCTGAACTCC"
SITE_A2_F6 = "GTCGTATTAGTAGTANNNNNNAGATCGGAAGAGCACACGTCTGAACTCC"
SITE_A2_F5 = "GTCGTATTAGTAGTANNNNNAGATCGGAAGAGCACACGTCTGAACTCC"
SITE_RECOVERY_F = "GGAGTTCAGACGTGTGCTC"
SITE_RECOVERY_R = "GTTGACATGCTGGATTGAGACTTC"
SITE_INDEX_F = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACG"
SITE_INDEX_R = "CAAGCAGAAGACGGCATACGAGATTANNNNNNNNGACTGGAGTTCAGACGTGTGCTC"


DESIGNS = {
    "bless": BreakCaptureDesign("BLESS", "biotinylated proximal hairpin",
                                 "distal hairpin", "intact fixed nuclei", False),
    "i-bless": BreakCaptureDesign("i-BLESS", "biotinylated proximal hairpin",
                                   "distal hairpin", "cells embedded in agarose beads", False),
    "dsbcapture": BreakCaptureDesign("DSBCapture", "biotinylated modified P5 adapter",
                                      "modified P7 adapter", "intact fixed nuclei", True),
    "site-seq": BreakCaptureDesign("SITE-Seq", "biotinylated Adapter 1",
                                    "degenerate Adapter 2", "purified HMW genomic DNA", True),
}


def _bless_workflow(key: str):
    d = DESIGNS[key]
    initial = molecule("yeast cells" if key == "i-bless" else d.immobilization,
                       (role("genomic DNA with a native DSB"),))
    prefix = ((
        "Embed, spheroplast and deproteinate cells in agarose beads",
        molecule("agarose-protected genomic DNA", (role("native DSB in deproteinated DNA"),)),
        "Agarose-bead immobilization is the defining i-BLESS change and protects yeast DNA from handling-induced shearing.",
    ),) if key == "i-bless" else ()
    return chain(initial, prefix + (
        ("Blunt and 5′-phosphorylate native break ends",
         molecule("ligation-ready native breaks", (role("blunt 5-prime-phosphorylated DSB"),)),
         "End repair makes 5′ overhangs, 3′ overhangs and blunt breaks enter the same ligation state."),
        ("Ligate the biotinylated proximal hairpin in situ",
         molecule("proximal-labelled break", (role("biotinylated proximal hairpin"), role("break-adjacent genomic DNA")),
                  junctions=(("biotinylated proximal hairpin", "break-adjacent genomic DNA", "T4 DNA ligase"),)),
         "The proximal barcode marks the end that existed in the cell before extraction."),
        ("Fragment DNA and capture proximal-labelled molecules on streptavidin",
         molecule("captured break fragment", (role("biotinylated proximal hairpin"), role("break-to-random genomic fragment"), role("fragmentation-created free end"))),
         "Fragmentation creates the distal end only after the native break has already been labelled."),
        ("Blunt the free end and ligate the distal hairpin",
         molecule("dual-hairpin break fragment", (role("proximal hairpin"), role("break-to-random genomic fragment"), role("distal hairpin")),
                  junctions=(("proximal hairpin", "break-to-random genomic fragment", "native-break ligation"),
                             ("break-to-random genomic fragment", "distal hairpin", "distal ligation"))),
         "Different proximal and distal barcodes preserve which end is the native DSB."),
        ("Linearize with I-SceI, PCR-amplify, then cut with XhoI",
         molecule("barcoded BLESS insert", (role("proximal 11-nt tag"), role("break-to-random genomic DNA"), role("distal 11-nt tag"))),
         "I-SceI releases the captured molecule from its hairpins; XhoI leaves the short end-identifying tags."),
        ("INFERRED — install the sequencing-kit ends",
         molecule("Illumina-compatible break library", (role("Illumina left end", inferred=True), role("proximal or distal tag"), role("genomic insert"), role("Illumina right end", inferred=True))),
         "The paper names the final library kit but does not print its complete adapter construct."),
    ))


def _dsbcapture_workflow():
    return chain(molecule("fixed nuclei", (role("genomic DNA with a native DSB"),)), (
        ("Blunt native break ends and add one 3′ dA",
         molecule("A-tailed native break", (role("repaired DSB"), role("3-prime dA"))),
         "End repair normalizes the break; Klenow exo− then creates a cohesive ligation end."),
        ("Ligate the biotinylated, T-overhang modified P5 adapter in situ",
         molecule("P5-labelled break", (role("biotinylated modified P5 adapter"), role("break-adjacent genomic DNA")),
                  junctions=(("biotinylated modified P5 adapter", "break-adjacent genomic DNA", "T/A ligation"),)),
         "Because the Read 1 site is already at the native break, R1 starts with the first genomic base."),
        ("Digest free adapter with lambda exonuclease, extract and sonicate",
         molecule("labelled sonication fragment", (role("modified P5 at native break"), role("break-to-random genomic fragment"), role("sonication-created free end"))),
         "Covalently break-labelled fragments retain biotin; free adapters are removed before extraction."),
        ("Streptavidin-capture, blunt and add one 3-prime dA to the distal end",
         da_tailed("distal-end dA-tailed DSBCapture fragment", "modified-P5-labelled break-to-random genomic DNA"),
         "The native-break P5 end is already labelled; only the sonication-created distal end is prepared."),
        ("Ligate modified P7 to the dA-tailed distal end",
         molecule("dual-adapter DSBCapture fragment", (role("modified P5"), role("break-to-random genomic DNA"), role("modified P7")),
                  junctions=(("modified P5", "break-to-random genomic DNA", "native-break ligation"),
                             ("break-to-random genomic DNA", "modified P7", "distal T/A ligation"))),
         "The second adapter is added only to the random sonication end."),
        ("PCR-add P5, i7 and P7; size-select the library",
         molecule("DSBCapture sequencing library", (role("P5 and Read 1 site"), role("DNA beginning at the DSB"), role("Read 2 and Index 1 site"), role("six-base sample index"), role("P7"))),
         "The published PCR primers complete a single-index Illumina library."),
    ))


def _site_seq_workflow():
    return chain(molecule("purified HMW genomic DNA", (role("intact genomic DNA"),)), (
        ("Cleave naked genomic DNA with Cas9–sgRNA in vitro",
         molecule("Cas9-digested genome", (role("genomic DNA"), role("programmable Cas9 cleavage ends"))),
         "A concentration series of each sgRNP exposes biochemically accessible on- and off-target sites; this is not a measurement of cleavage inside cells."),
        ("Add one 3-prime dA to each Cas9 cleavage end",
         da_tailed("dA-tailed Cas9-cleaved genome", "Cas9-cleavage-adjacent genomic DNA"),
         "The Cas9-created ends are shown in their cohesive ligation-ready state."),
        ("Ligate biotinylated Adapter 1 to the dA-tailed cleavage ends",
         molecule("cut-end-labelled genome", (role("5-prime-biotin Adapter 1"), role("Cas9-cleavage-adjacent genomic DNA")),
                  junctions=(("5-prime-biotin Adapter 1", "Cas9-cleavage-adjacent genomic DNA", "T/A ligation"),)),
         "Adapter 1 is installed before fragmentation, so its biotin identifies the Cas9-created end."),
        ("Fragment with dsDNA Fragmentase to about 200–800 bp",
         molecule("labelled genomic fragments", (role("biotinylated Cas9 end"), role("cut-to-random genomic fragment"), role("Fragmentase-created distal end"))),
         "The new distal end is physically distinct from the already labelled cleavage end."),
        ("End-repair and ligate the pooled N5/N6/N7 Adapter 2",
         molecule("dual-adapter SITE-Seq fragment", (role("biotinylated Adapter 1"), role("cut-to-random genomic DNA"), role("5/6/7-nt diversity spacer + Adapter 2")),
                  junctions=(("biotinylated Adapter 1", "cut-to-random genomic DNA", "cleavage-end ligation"),
                             ("cut-to-random genomic DNA", "5/6/7-nt diversity spacer + Adapter 2", "distal ligation"))),
         "Pooling three adapter lengths creates base diversity; the paper does not define this tract as a molecule-counting UMI."),
        ("Capture biotinylated fragments on streptavidin beads",
         molecule("enriched SITE-Seq fragments", (role("Adapter 1-labelled Cas9 end"), role("adjacent genomic DNA"), role("Adapter 2 distal end"))),
         "Only molecules retaining the Adapter 1-labelled cleavage end survive the affinity selection."),
        ("Recovery PCR, then dual-index PCR",
         molecule("SITE-Seq sequencing library", (role("P5"), role("i5 sample index"), role("Read 1 site"), role("DNA beginning at a Cas9 cut"), role("5/6/7-nt diversity spacer"), role("Index 1 site"), role("i7 sample index"), role("P7"))),
         "The two published PCR rounds first recover captured molecules and then install both flow-cell ends and eight-base sample indexes."),
    ))


def _digenome_workflow():
    return chain(molecule("purified genomic DNA", (role("intact genomic DNA"),)), (
        ("Cleave naked genomic DNA with Cas9–sgRNA in vitro",
         molecule("Cas9-digested genome", (role("uncut genomic DNA"), role("programmable Cas9 cleavage ends"))),
         "Cleavage is performed on cell-free DNA. All resulting DNA—not an enriched subset—continues into whole-genome library preparation."),
        ("Remove sgRNA and protein; purify the digested DNA",
         molecule("clean digenome", (role("genomic DNA retaining Cas9-created ends"),)),
         "RNase A and protein removal stop the reaction without adding a molecular tag to the cleavage ends."),
        ("Covaris-shear the whole digenome",
         molecule("mixed-origin WGS fragments", (role("fragments carrying one Cas9 end and one shear end"), role("fragments carrying two random shear ends"))),
         "Cas9 ends survive among a much larger background of ordinary random WGS fragment ends."),
        ("End-repair, size-select and add one 3-prime dA to each end",
         da_tailed("dA-tailed Digenome-seq WGS fragment", "unenriched genomic fragment"),
         "Both the Cas9-derived and randomly sheared endpoint classes pass through this end-prepared state."),
        ("Ligate TruSeq adapters without PCR",
         molecule("PCR-free Digenome-seq library", (role("TruSeq P5 arm"), role("unenriched genomic fragment"), role("TruSeq indexed P7 arm"))),
         "PCR-free preparation preserves read-start counts; there is no cut-end affinity capture or locus-specific amplification."),
        ("Paired-end whole-genome sequencing and cleavage-score analysis",
         molecule("aligned digenome reads", (role("forward and reverse 5-prime read-start pileups"), role("genome-wide depth background"))),
         "A true nuclease cut produces coincident read starts on both strands; ordinary shear endpoints remain dispersed."),
    ))


def _digenome_library() -> tuple[Construct, tuple]:
    return truseq_library(
        [seg("unenriched genomic DNA; Cas9 or random-shear boundary at either end",
             "X" * 42, placeholder=True)],
        "Digenome-seq TruSeq PCR-free WGS library", dual_index=False)


def _site_seq_library() -> tuple[Construct, tuple]:
    """Indexed product, with the N7 Adapter 2 branch shown explicitly."""
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5 sample index", "I" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i5", "sample_index", "unknown")),
        Segment("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        role("genomic DNA beginning at the in-vitro Cas9 cleavage end"),
        Segment("Adapter 2 length-diversity spacer (N7 branch)", "N" * 7, "cbc", placeholder=True),
        Segment("Index 1 + Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="SITE-Seq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("SITE-Seq final library: " + "; ".join(problems))
    return lib, primers


def _illumina_break_library(name: str, *, inferred: bool, tagged: bool) -> tuple[Construct, tuple]:
    left = Segment("P5 + Read 1 arm", il.TRUSEQ_P5_FULL, "p5", inferred=inferred)
    payload = ([Segment("proximal read tag", BLESS_PROXIMAL_READ_TAG, "cbc",
                        feature=feature("proximal_break_end", "feature_barcode", "fixed")),
                role("break-to-random genomic DNA"),
                Segment("distal read tag", BLESS_DISTAL_READ_TAG, "cbc",
                        feature=feature("distal_fragment_end", "feature_barcode", "fixed"))]
               if tagged else [role("DNA beginning at repaired DSB")])
    lib = Construct([left, *payload, Segment("dA junction", "A", inferred=inferred),
                     Segment("Index 1 + Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=inferred),
                     Segment("six-base i7 reverse complement", "I" * 6, "cbc", placeholder=True,
                             feature=feature("sample_index_i7", "sample_index", "unknown"), inferred=inferred),
                     Segment("P7 reverse complement", il.P7_RC, "p7", inferred=inferred)], name=name)
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def get(key: str):
    if key not in DESIGNS and key != "digenome-seq":
        raise KeyError(key)
    if key == "digenome-seq":
        lib, primers = _digenome_library()
        return SimpleNamespace(
            TITLE="Digenome-seq — digested-genome whole-genome sequencing",
            NOTES="01_digenome_seq.html",
            SOURCE='<a href="https://doi.org/10.1038/nmeth.3284">Kim et al., Nature Methods (2015)</a>; library architecture follows the specified Illumina TruSeq DNA PCR-free WGS chemistry.',
            SUMMARY="Cas9–sgRNA cleaves purified genomic DNA in vitro, after which the entire unenriched genome receives an ordinary PCR-free whole-genome library preparation. Coincident forward- and reverse-strand read starts reveal cleavage sites.",
            CAVEAT="Digenome-seq has no cut-site adapter, UMI or enrichment barcode. The library shell is the commercial TruSeq PCR-free chemistry used for WGS; the assay-specific information is the genomic fragment boundary itself.",
            oligos=lambda: (), workflow=_digenome_workflow,
            FINAL_LIBRARY=lib,
            FINAL_CAPTION="One ordinary TruSeq PCR-free WGS molecule. Its genomic insert may carry a Cas9-created boundary or only random shearing boundaries; the molecule itself has no label distinguishing them.",
            SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 8, "Read 2": 150},
            SEQUENCING_INTRO="The defining study used Illumina HiSeq X Ten whole-genome sequencing. Digenome analysis depends on paired genomic alignments and the concentration of identical 5′ endpoints, not an inline molecular identifier.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
    if key == "site-seq":
        lib, primers = _site_seq_library()
        return SimpleNamespace(
            TITLE="SITE-Seq — in-vitro Cas9 cleavage-site sequencing", NOTES="01_site_seq.html",
            SOURCE='<a href="https://doi.org/10.1038/nmeth.4284">Cameron et al., Nature Methods (2017)</a>; exact chemistry from the published supplementary protocol and Supplementary Table 26.',
            SUMMARY="Cas9–sgRNA cleaves purified genomic DNA in vitro; a biotinylated adapter labels those ends before enzymatic fragmentation, allowing streptavidin enrichment and Read 1 entry directly beside each candidate cut site.", CAVEAT="",
            oligos=lambda: (
                f"Adapter 1 Forward — 5′-/5Biosg/{SITE_A1_F[:-1]}*{SITE_A1_F[-1]}-3′",
                f"Adapter 1 Reverse — 5′-{SITE_A1_R}-3′",
                f"Adapter 2 N7 Forward — 5′-/5Phos/{SITE_A2_F7}-3′",
                f"Adapter 2 N6 Forward — 5′-/5Phos/{SITE_A2_F6}-3′",
                f"Adapter 2 N5 Forward — 5′-/5Phos/{SITE_A2_F5}-3′",
                f"Adapter 2 Reverse — 5′-{SITE_A2_R[:-1]}*{SITE_A2_R[-1]}-3′",
                f"Recovery PCR Forward — 5′-{SITE_RECOVERY_F}-3′",
                f"Recovery PCR Reverse — 5′-{SITE_RECOVERY_R}-3′",
                f"Index Forward — 5′-{SITE_INDEX_F}-3′",
                f"Index Reverse — 5′-{SITE_INDEX_R}-3′"),
            workflow=_site_seq_workflow, FINAL_LIBRARY=lib,
            FINAL_CAPTION="N7 Adapter 2 branch; N5 and N6 differ only in diversity-spacer length. Adapter 1 defines the Cas9-cleavage end, while Adapter 2 occupies the later Fragmentase end.",
            SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 8, "Index 2 (i5)": 8},
            SEQUENCING_INTRO="The supplementary protocol recommends 2–3 million 150-base single-end reads per sample. Both published index primers carry eight variable bases; the standard index primers land on the completed dual-index library.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")

    if key == "dsbcapture":
        lib, primers = _illumina_break_library("DSBCapture library", inferred=False, tagged=False)
        return SimpleNamespace(
            TITLE="DSBCapture — direct in-situ sequencing of DNA breaks", NOTES="01_dsbcapture.html",
            SOURCE='Lensing et al., Nature Methods (2016); <a href="https://doi.org/10.1038/nmeth.3960">doi:10.1038/nmeth.3960</a>.',
            SUMMARY="A biotinylated modified P5 adapter is ligated to repaired breaks inside fixed nuclei; after sonication and capture, modified P7 is installed only at the random distal end, so Read 1 begins at the DSB.", CAVEAT="",
            oligos=lambda: (
                "modified P5 top — 5′-/5Biosg/A*ATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′",
                "modified P5 bottom — 5′-/5Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT-3′",
                "modified P7 top — 5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC-3′",
                "modified P7 bottom — 5′-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3′",
                f"PCR F — 5′-{DSB_PCR_F}-3′",
                f"PCR R (AD06 example) — 5′-{DSB_PCR_R_AD06}-3′"),
            workflow=_dsbcapture_workflow, FINAL_LIBRARY=lib,
            FINAL_CAPTION="Published modified-P5/modified-P7 library; the proximal orientation is shown. Read 1 starts at the repaired native break.",
            SEQ_PRIMERS=primers, READ_LENGTHS={"Read 1": 75, "Index 1 (i7)": 6, "Read 2": 75},
            SEQUENCING_INTRO="The paper used 2×75-bp NextSeq runs (and 2×50-bp MiSeq for the EcoRV pilot); downstream break calling retains the proximal Read 1.", SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")

    lib, primers = _illumina_break_library(f"{DESIGNS[key].name} library", inferred=True, tagged=True)
    is_i = key == "i-bless"
    return SimpleNamespace(
        TITLE=("i-BLESS — immobilized Breaks Labeling, Enrichment on Streptavidin and Sequencing" if is_i else
               "BLESS — Breaks Labeling, Enrichment on Streptavidin and Sequencing"),
        NOTES=("01_i_bless.html" if is_i else "01_bless.html"),
        SOURCE=(
            'Biernacka et al., Communications Biology (2018); <a href="https://doi.org/10.1038/s42003-018-0165-9">doi:10.1038/s42003-018-0165-9</a>.' if is_i else
            'Crosetto et al., Nature Methods (2013); <a href="https://doi.org/10.1038/nmeth.2408">doi:10.1038/nmeth.2408</a>.'),
        SUMMARY=(("Agarose beads immobilize yeast cells through lysis and direct break labelling; the otherwise BLESS-like proximal/distal hairpins retain which sequenced end is the native DSB.") if is_i else
                 "A barcoded biotin hairpin labels repaired DSBs in fixed nuclei before fragmentation; streptavidin capture and a differently barcoded distal hairpin preserve the native-break end."),
        CAVEAT=("The proximal/distal tags are exact. The final ThruPLEX adapter geometry is not public and is shown as an inferred Illumina-compatible shell." if is_i else
                "The proximal/distal tags are exact. The paper names the final Illumina library-preparation reagents but does not print that kit's complete adapter geometry; those outer regions are marked inferred."),
        oligos=lambda: (
            f"proximal hairpin — 5′-/5Phos/{BLESS_PROXIMAL[:45]}/iBiodT/{BLESS_PROXIMAL[46:]}-3′",
            f"distal hairpin — 5′-/5Phos/{BLESS_DISTAL}-3′",
            f"proximal PCR primer — 5′-{BLESS_PCR_PROXIMAL}-3′",
            f"distal PCR primer — 5′-{BLESS_PCR_DISTAL}-3′"),
        workflow=lambda: _bless_workflow(key), FINAL_LIBRARY=lib,
        FINAL_CAPTION="One insert orientation is shown. Proximal and distal 11-nt tags identify which genomic end was the native break; inferred outer sequencing-kit regions are shaded.",
        SEQ_PRIMERS=primers, READ_LENGTHS={"Read 1": 75, "Index 1 (i7)": 6, "Read 2": 75},
        SEQUENCING_INTRO=("The defining i-BLESS study used 2×75-bp HiSeq 2500/4000 sequencing." if is_i else
                          "BLESS used paired-end Illumina sequencing; either end can carry the proximal tag, which is selected computationally."),
        SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
