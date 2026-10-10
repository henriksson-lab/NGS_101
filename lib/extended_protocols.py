"""Shared models for named protocols with unusual partitioning or molecular readouts."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import feature
from protocol_model import chain, molecule, role, role_library


def ident(name, identifier, kind, *, length=None, tag="cbc"):
    return role(name, tag=tag, length_bp=length,
                molecular_feature=feature(identifier, kind, "unknown"))


CELL = ident("cell barcode", "cell_barcode", "cell_barcode")
UMI = ident("UMI", "umi", "umi", tag="umi")
SAMPLE = ident("sample index", "sample_index", "sample_index", length=8)
TARGET = ident("target barcode", "target_barcode", "feature_barcode")
GUIDE = ident("guide identity", "guide_barcode", "guide_barcode")
RT_WELL = ident("RT-well barcode", "rt_well", "cell_barcode")
LIG_WELL = ident("ligation-well barcode", "ligation_well", "cell_barcode")
PCR_WELL = ident("PCR-well barcode", "pcr_well", "cell_barcode")


def _wellda():
    return chain(molecule("single nucleus in a well", (role("genomic DNA"), role("accessible chromatin"))), (
        ("Tagment accessible chromatin with indexed Tn5",
         molecule("ATAC-marked genome", (CELL, role("accessible-DNA fragments"), role("remaining genome"))),
         "Tagmentation records accessibility and supplies the first molecular handle."),
        ("Displace the tagged products and amplify the residual genome",
         molecule("paired products from one nucleus", (CELL, role("ATAC branch"), role("whole-genome branch"))),
         "The branches preserve their common well identity while measuring accessibility and copy number."),
        ("Add branch-compatible Illumina indexes and pool wells",
         molecule("wellDA-seq libraries", (CELL, SAMPLE, role("ATAC or genomic insert"))),
         "The same single-cell well barcode links the epigenome and genome."),
    ))


def _easysci():
    return chain(molecule("fixed nuclei distributed into RT wells", (role("nuclear RNA"),)), (
        ("Reverse-transcribe with matched barcoded oligo(dT) and random-hexamer primers",
         molecule("first-index cDNA in nuclei", (RT_WELL, UMI, role("cDNA"))),
         "The paired priming modes recover polyadenylated and broader nuclear RNA from the same RT well."),
        ("Pool, redistribute and ligate the second well barcode",
         molecule("two-round indexed nuclei", (RT_WELL, LIG_WELL, UMI, role("cDNA")),
                  junctions=(("RT-well barcode", "ligation-well barcode", "ligation"),)),
         "The ordered RT-plus-ligation pair narrows cell identity before PCR."),
        ("Pool, make second strand, tagment and add the PCR-well index",
         molecule("EasySci-RNA library", (RT_WELL, LIG_WELL, PCR_WELL, UMI, role("transcript insert"), SAMPLE)),
         "Three combinatorial rounds identify a nucleus; the P5 barcode is read separately from the P7 sample index."),
    ))


def _indrops2():
    return chain(molecule("cell and barcoded hydrogel bead in a droplet", (role("poly(A) RNA"), CELL, UMI)), (
        ("Lyse the cell and reverse-transcribe RNA on the bead oligo",
         molecule("barcoded first-strand cDNA", (CELL, UMI, role("3-prime transcript cDNA"))),
         "Cell barcode and UMI are copied into every captured molecule."),
        ("Choose linear IVT amplification or template-switching PCR amplification",
         molecule("amplified inDrops-2 cDNA", (CELL, UMI, role("transcript cDNA"), role("IVT or TSO handle"))),
         "The defining protocol validates two distinct amplification branches without treating them as identical chemistry."),
        ("Fragment or tagment the amplified product and complete Illumina adapters",
         molecule("inDrops-2 library", (CELL, UMI, role("3-prime transcript insert"), SAMPLE)),
         "Both branches retain the same barcode and UMI readout."),
    ))


def _multitag():
    return chain(molecule("permeabilized cell with chromatin targets", (role("target 1"), role("target 2"), role("target 3"))), (
        ("Bind the first primary antibody and its barcoded pA–Tn5 reagent",
         molecule("first target tagged", (TARGET, role("target-1 chromatin fragment"))),
         "The transposon barcode identifies the antibody target."),
        ("Repeat antibody-directed tagmentation sequentially for additional targets",
         molecule("multiplex tagged chromatin", (TARGET, role("target-proximal genomic DNA"))),
         "Sequential binding and stringent washes retain target-specific adapter identities in one cell."),
        ("Combinatorially index single cells and amplify",
         molecule("MulTI-Tag library", (CELL, TARGET, role("target-proximal genomic insert"), SAMPLE)),
         "Target barcode and cell barcode are independent structured identifiers."),
    ))


def _tipseq():
    return chain(molecule("antibody-bound chromatin", (role("chromatin protein target"),)), (
        ("Tagment target-proximal DNA with pA–Tn5 carrying a T7 promoter",
         molecule("promoter-tagged chromatin", (role("T7 promoter"), role("target-proximal DNA"),
                                                role("mosaic end"))),
         "Targeted insertion turns each bound site into a transcription template."),
        ("Transcribe tagged fragments with T7 RNA polymerase",
         molecule("linearly amplified RNA", (role("T7-derived RNA"), role("target-proximal sequence"))),
         "In-vitro transcription provides linear amplification before cDNA conversion."),
        ("Reverse-transcribe, amplify and index",
         molecule("TIP-seq library", (CELL, role("target-proximal genomic insert"), SAMPLE)),
         "The bulk and sciTIP implementations share promoter insertion; the single-cell form adds indexed transposomes."),
    ))


def _phage_atac():
    return chain(molecule("cell stained with nanobody-displaying M13 phage", (role("surface or intracellular epitope"), TARGET)), (
        ("Fix, lyse and tagment accessible chromatin in bulk",
         molecule("phage-bound tagmented cell", (TARGET, role("PAC-tag phagemid"), role("ATAC fragments"))),
         "The engineered PAC/RD1 sequence makes the phagemid compatible with 10x ATAC capture."),
        ("Co-encapsulate cells and 10x ATAC gel beads",
         molecule("droplet-barcoded products", (CELL, TARGET, role("phagemid tag"), role("ATAC fragment"))),
         "Linear amplification applies one droplet barcode to chromatin and phage-derived tags."),
        ("Split protein-derived-tag and ATAC libraries",
         molecule("PHAGE-ATAC readouts", (CELL, TARGET, role("ATAC or phage insert"), SAMPLE)),
         "The common cell barcode links protein abundance, accessibility and captured mitochondrial variants."),
    ))


def _sctem():
    return chain(molecule("one sorted cell", (role("genomic DNA"), role("mRNA"))), (
        ("Separate genomic DNA and mRNA for parallel processing",
         molecule("matched DNA and RNA branches", (CELL, role("genomic DNA"), role("full-length cDNA"))),
         "The split permits Smart-seq-type transcriptomics alongside targeted methylation."),
        ("Bisulfite-convert genomic DNA",
         molecule("converted DNA", (role("converted SINE Alu / LINE-1 templates"),)),
         "Unmethylated cytosines convert while methylated CpGs remain cytosine."),
        ("Amplify transposable elements with staggered, internally indexed primer pairs",
         molecule("scTEM amplicons", (CELL, role("internal dual index"), role("SINE Alu insert"), SAMPLE)),
         "A 0–5 N stagger and reversed primer orientation diversify the low-complexity amplicon pool."),
    ))


def _end5():
    return chain(molecule("single-cell RNA in a 10x droplet", (role("poly(A) and non-poly(A) RNA"),)), (
        ("Prime first-strand synthesis with oligo(dT) or random hexamers",
         molecule("nascent first strand", (role("transcript cDNA"),)),
         "sc-end5-dT and sc-end5-rand are explicit priming variants of the named method."),
        ("Template-switch at the RNA 5′ end onto the gel-bead oligo",
         molecule("5-prime-end cDNA", (CELL, UMI, role("template-switch junction"), role("transcript cDNA"))),
         "Barcode and UMI become associated with a nucleotide-resolution transcript 5′ end."),
        ("Amplify and make the sequencing library",
         molecule("sc-end5-seq library", (CELL, UMI, role("RNA 5-prime insert"), SAMPLE)),
         "The genomic end of the insert reports a candidate transcription start site."),
    ))


def _speac():
    return chain(molecule("pooled CRISPR-perturbed effector cells", (GUIDE, role("reporter cells"))), (
        ("Pair one perturbed cell with one reporter cell in a picolitre droplet",
         molecule("co-culture droplet", (GUIDE, role("perturbed cell"), role("fluorescent reporter cell"))),
         "Different dyes gate one-to-one pairing before the interaction develops."),
        ("Sort droplets by the functional fluorescent response",
         molecule("phenotype-enriched droplet fraction", (GUIDE, role("reporter state"))),
         "The screen's phenotype is a sorted population, not a single-cell transcriptome."),
        ("Break droplets, recover genomic DNA and PCR-amplify the integrated guide cassette",
         molecule("SPEAC-seq library", (GUIDE, SAMPLE, role("guide-flanking amplicon"))),
         "Guide counts in response-positive and control fractions identify interaction regulators."),
    ))


def _disco():
    return chain(molecule("low-input cell suspension and mRNA-capture beads", (role("single cells"), CELL, UMI)), (
        ("Machine-vision position exactly one cell and one bead at separate valves",
         molecule("deterministically paired cell and bead", (role("one cell"), CELL, UMI)),
         "Active positioning replaces Poisson co-loading for samples below about 500 cells."),
        ("Co-eject the pair into one droplet and sort it to the sample outlet",
         molecule("occupied DisCo droplet", (role("lysed cell RNA"), CELL, UMI)),
         "The device rejects empty or mis-paired events before molecular capture."),
        ("Capture poly(A) RNA and prepare a Drop-seq-derived library",
         molecule("DisCo library", (CELL, UMI, role("3-prime transcript insert"), SAMPLE)),
         "The molecular construct is bead-based 3′ RNA-seq; the innovation is deterministic partitioning."),
    ))


def _spindrop():
    return chain(molecule("fluorescently labelled cell and inDrop-v3 microgel", (role("cell or nucleus"), CELL, UMI)), (
        ("Co-encapsulate and fluorescence-sort occupied, viable or target droplets",
         molecule("selected droplet", (role("selected cell"), CELL, UMI)),
         "Inline FADS removes empty, damaged or unwanted droplets."),
        ("Heat-treat and picoinject an improved reverse-transcription mixture",
         molecule("barcoded first-strand cDNA", (CELL, UMI, role("3-prime transcript cDNA"))),
         "A second device adds RT reagent after sorting, enabling a multistep droplet reaction."),
        ("Recover, amplify and index cDNA",
         molecule("spinDrop library", (CELL, UMI, role("3-prime transcript insert"), SAMPLE)),
         "The final read structure follows the inDrop-v3 barcode scheme."),
    ))


def _limca():
    return chain(molecule("one intact cell", (role("cytoplasmic mRNA"), role("nuclear chromatin"))), (
        ("Physically separate cytoplasm and nucleus",
         molecule("matched compartments", (CELL, role("cytoplasm"), role("nucleus"))),
         "Compartment separation lets each modality retain its high-sensitivity chemistry."),
        ("Amplify cytoplasmic RNA by Smart-seq2",
         molecule("full-length transcriptome branch", (CELL, role("full-length cDNA"))),
         "The RNA branch reports gene expression from the same cell."),
        ("Crosslink, restriction-digest and proximity-ligate the nucleus",
         molecule("single-cell contact junctions", (role("restriction fragment A"), role("restriction fragment B")),
                  junctions=(("restriction fragment A", "restriction fragment B", "proximity ligation"),)),
         "The nuclear branch uses conventional 3C chemistry."),
        ("Apply META whole-genome amplification and build the Hi-C library",
         molecule("paired LiMCA libraries", (CELL, role("Smart-seq2 insert or Hi-C junction"), SAMPLE)),
         "Transcript and contact libraries are integrated by their shared cell-of-origin record."),
    ))


def _buttseq():
    return chain(molecule("purified chromatin-associated nascent RNA", (role("nascent RNA 5-prime body"), role("native RNA 3-prime end"))), (
        ("TGIRT template-switch from an RNA/DNA starter duplex onto the RNA 3′ end",
         molecule("first-strand cDNA", (UMI, role("nascent RNA 3-prime-end cDNA"))),
         "The first copied base reports the polymerase-proximal RNA end without RNA fragmentation."),
        ("Switch salt concentration and synthesize the second strand",
         molecule("double-stranded Butt-seq product", (UMI, role("nascent-RNA insert"))),
         "The salt switch favors the intended first template switch and suppresses later spurious switching."),
        ("PCR-amplify, size-select and index",
         molecule("Butt-seq library", (UMI, role("nascent-RNA terminal insert"), SAMPLE)),
         "A separate UMI read supports molecule-aware paired-end analysis."),
    ))


@dataclass(frozen=True)
class Definition:
    title: str; doi: str; citation: str; summary: str; endpoint: str; notes: str


DEFS = {
    "wellda": Definition("wellDA-seq", "10.1101/2024.09.06.611519", "Wang et al., bioRxiv (2024)", "Indexed accessibility tagging and a whole-genome branch co-profile chromatin accessibility and copy-number variation from the same single cells.", "Matched ATAC and genomic reads share a single-cell well identity.", "01_wellda_seq.html"),
    "easysci": Definition("EasySci-RNA", "10.1038/s41588-023-01572-y", "Sziraki et al., Nature Genetics (2023)", "Barcoded oligo(dT) and random-hexamer RT, barcode ligation and indexed PCR create a scalable three-round combinatorial single-nucleus RNA library.", "RT-well, ligation-well and PCR-well indexes jointly identify each nucleus; a UMI identifies each cDNA molecule.", "01_easysci.html"),
    "indrops2": Definition("inDrops-2", "10.1093/nar/gkae1312", "Juzėnas et al., Nucleic Acids Research (2025)", "Open-source droplet capture supports either IVT linear amplification or template-switching PCR amplification of bead-barcoded 3′ cDNA.", "The barcode read contains cell barcode and UMI; the insert read reports the transcript 3′ end.", "01_indrops_2.html"),
    "multitag": Definition("MulTI-Tag", "10.1038/s41587-022-01522-9", "Meers et al., Nature Biotechnology (2022)", "Sequential antibody-targeted pA–Tn5 reactions insert target-specific adapter barcodes before single-cell combinatorial indexing.", "Each genomic fragment carries both protein-target identity and cell identity.", "01_multi_tag.html"),
    "tipseq": Definition("TIP-seq", "10.1083/jcb.202103078", "Bartlett et al., Journal of Cell Biology (2021)", "Antibody-tethered Tn5 inserts a T7 promoter at chromatin targets so in-vitro transcription linearly amplifies adjacent genomic DNA.", "The genomic read maps target-proximal DNA; indexed transposomes provide cell identity in sciTIP-seq.", "01_tip_seq.html"),
    "phage-atac": Definition("PHAGE-ATAC", "10.1038/s41587-021-01065-5", "Fiskin et al., Nature Biotechnology (2022)", "Nanobody-displaying phagemids carry 10x-compatible PAC tags so one droplet barcode joins protein-derived tags to ATAC and mitochondrial DNA.", "Separate protein-tag and ATAC libraries share the 10x cell barcode.", "01_phage_atac.html"),
    "sctem": Definition("scTEM-seq", "10.1038/s41598-022-09765-x", "Hunt et al., Scientific Reports (2022)", "Bisulfite-converted single-cell DNA is amplified with indexed SINE Alu or LINE-1 primers while RNA from the same cell is processed in parallel.", "Paired reads span targeted transposable-element CpGs; internal and adapter indexes identify the cell library.", "01_sctem_seq.html"),
    "end5": Definition("sc-end5-seq", "10.1101/2021.04.04.438388", "Moody et al., bioRxiv (2021)", "Oligo(dT) or random priming followed by 10x template switching retains nucleotide-resolution RNA 5′ ends with cell barcodes and UMIs.", "The insert's genomic 5′ boundary reports a candidate transcription start site; cell barcode and UMI support single-cell counting.", "01_sc_end5_seq.html"),
    "speac": Definition("SPEAC-seq", "10.1038/s41596-024-01056-1", "Faust Akl et al., Nature Protocols (2025)", "CRISPR-perturbed and reporter cells are paired and cultured in droplets; functional fluorescence sorting precedes guide-cassette PCR and sequencing.", "Guide counts are compared between functionally sorted and control droplet fractions; this is a pooled screen, not a transcriptome library.", "01_speac_seq.html"),
    "disco": Definition("DisCo", "10.1038/s41592-021-01391-1", "Bues et al., Nature Methods (2022)", "Machine vision and multilayer valves deterministically co-encapsulate one low-input cell and one mRNA-capture bead for Drop-seq-derived 3′ RNA sequencing.", "The bead barcode identifies the captured cell and the UMI identifies each transcript molecule.", "01_disco.html"),
    "spindrop": Definition("spinDrop", "10.1038/s41467-023-40322-w", "Hochgerner et al., Nature Communications (2023)", "Fluorescence sorting enriches occupied inDrop-v3 droplets before a second device picoinjects reverse-transcription reagent.", "The inDrop-v3 barcode read identifies cell and molecule; the insert read reports the transcript 3′ end.", "01_spindrop.html"),
    "limca": Definition("LiMCA", "10.1038/s41592-024-02239-0", "Wu et al., Nature Methods (2024)", "Cytoplasm and nucleus are separated: Smart-seq2 reads full-length RNA while restriction digestion, proximity ligation and META amplification recover high-depth chromatin contacts.", "Matched full-length transcript and Hi-C libraries are joined by the recorded cell of origin.", "01_limca.html"),
    "buttseq": Definition("Butt-seq", "10.1101/gad.350434.123", "Yu and Rosbash, Genes & Development (2023)", "TGIRT template switching captures native nascent-RNA 3′ termini without RNA fragmentation, followed by second-strand synthesis and indexed PCR.", "The terminal insert base reports the nascent RNA 3′ end; a separate UMI read supports molecule deduplication.", "01_butt_seq.html"),
}


WORKFLOWS = {"wellda": _wellda, "easysci": _easysci, "indrops2": _indrops2,
             "multitag": _multitag, "tipseq": _tipseq, "phage-atac": _phage_atac,
             "sctem": _sctem, "end5": _end5, "speac": _speac, "disco": _disco,
             "spindrop": _spindrop, "limca": _limca, "buttseq": _buttseq}


def get(key: str):
    d = DEFS[key]
    payloads = {
        "wellda": (CELL, role("ATAC or genomic insert")),
        "easysci": (RT_WELL, LIG_WELL, PCR_WELL, UMI, role("transcript insert")),
        "indrops2": (CELL, UMI, role("3-prime transcript insert")),
        "multitag": (CELL, TARGET, role("target-proximal genomic insert")),
        "tipseq": (CELL, role("target-proximal genomic insert")),
        "phage-atac": (CELL, TARGET, role("ATAC or phage insert")),
        "sctem": (CELL, role("internal dual index"), role("TE amplicon")),
        "end5": (CELL, UMI, role("RNA 5-prime insert")),
        "speac": (GUIDE, role("guide-flanking amplicon")),
        "disco": (CELL, UMI, role("3-prime transcript insert")),
        "spindrop": (CELL, UMI, role("3-prime transcript insert")),
        "limca": (role("Smart-seq2 insert or Hi-C junction"),),
        "buttseq": (UMI, role("nascent-RNA terminal insert")),
    }
    final = role_library(d.title + " final library roles", payloads[key], index=SAMPLE)
    oligos = ()
    if key == "tipseq":
        oligos = ("ME-reverse — 5′-CTGTCTCTTATACACATCT-3′",)
    elif key == "speac":
        oligos = ("lentiCRISPRv2 P7 index primer — 5′-CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCAATTCCCACTCCTTTCAAGACCT-3′",)
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE="The source defines the read roles, but the primer sequence is panel-, kit- or implementation-dependent; no universal binding oligo is asserted.",
        workflow=WORKFLOWS[key], oligos=lambda: oligos,
    )
