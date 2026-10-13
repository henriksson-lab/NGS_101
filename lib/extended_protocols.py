"""Shared models for named protocols with unusual partitioning or molecular readouts."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from adapter_ligation import AdapterLigationGeometry
from chemdraw import feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
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


def _finish_pcr(graph: ProtocolGraph, prefix: str, template, final,
                endpoint: str, reviewed: list[str], *, label: str = "library"):
    forward = graph.input(prefix + "_forward", label + " PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " indexed PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Index and amplify the " + label, StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The source defines read roles, but the complete primer geometry is panel-, kit- or implementation-dependent."))
    graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library


def _rt_then_pcr(graph: ProtocolGraph, prefix: str, template, final,
                 endpoint: str, reviewed: list[str], *, rt_primer_label: str,
                 label: str = "library"):
    primer = graph.input(prefix + "_rt_primer", rt_primer_label, MaterialKind.OLIGO)
    rt = graph.input(prefix + "_rt_enzyme", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", label + " first-strand cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule(label + " first-strand cDNA", (role("cDNA"),)))
    graph.react(prefix + "_rt", "Reverse-transcribe the " + label + " RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    reviewed.extend((prefix + "_rt_primer", prefix + "_rt_enzyme", prefix + "_rt"))
    return _finish_pcr(graph, prefix, cdna, final, endpoint, reviewed, label=label)


def _end5_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "single_cell_rna", "single-cell poly(A) and non-poly(A) RNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("single-cell RNA in a 10x droplet", (
            role("poly(A) and non-poly(A) RNA"),)), primary=True)
    dt_rna = graph.product("dt_input", "RNA assigned to sc-end5-dT",
                           MaterialKind.MOLECULAR, state=rna.state)
    random_rna = graph.product("random_input", "RNA assigned to sc-end5-rand",
                               MaterialKind.MOLECULAR, state=rna.state)
    graph.react("split_priming_variants", "Separate oligo(dT) and random-hexamer variants",
                StepKind.SPLIT, ((rna, InputRole.TEMPLATE),), (dt_rna, random_rna),
                changes_molecule=False)
    reviewed = ["single_cell_rna", "split_priming_variants"]

    def branch(prefix: str, template, primer_label: str, label: str):
        primer = graph.input(prefix + "_primer", primer_label, MaterialKind.OLIGO)
        rt = graph.input(prefix + "_rt_enzyme", "10x template-switching reverse transcriptase",
                         MaterialKind.CATALYST)
        first = graph.product(prefix + "_first_strand", label + " nascent first strand",
                              MaterialKind.MOLECULAR,
                              state=molecule("nascent first strand", (
                                  role("transcript cDNA"),)))
        graph.react(prefix + "_rt", "Prime first-strand synthesis with " + primer_label,
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((template, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (first,))
        bead_tso = graph.input(prefix + "_bead_tso", "cell-barcode/UMI gel-bead template-switch oligo",
                               MaterialKind.OLIGO)
        switched = graph.product(
            prefix + "_switched", label + " 5-prime-end cDNA",
            MaterialKind.MOLECULAR,
            state=molecule("5-prime-end cDNA", (
                CELL, UMI, role("template-switch junction"), role("transcript cDNA"))))
        graph.react(prefix + "_template_switch", "Template-switch at the RNA 5-prime end",
                    StepKind.TEMPLATE_SWITCHING,
                    ((first, InputRole.TEMPLATE), (bead_tso, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (switched,))
        reviewed.extend((prefix + "_primer", prefix + "_rt_enzyme", prefix + "_rt",
                         prefix + "_bead_tso", prefix + "_template_switch"))
        _finish_pcr(graph, prefix, switched, final, label + ": " + d.endpoint,
                    reviewed, label=label)

    branch("dt", dt_rna, "oligo(dT) RT primer", "sc-end5-dT")
    branch("random", random_rna, "random-hexamer RT primer", "sc-end5-rand")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _speac_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "cell_pools", "CRISPR-perturbed effector cells and fluorescent reporter cells",
        MaterialKind.BIOLOGICAL,
        state=molecule("pooled screen cells", (GUIDE, role("effector cells"),
                                                role("reporter cells"))), primary=True)
    droplets = graph.input("pairing_condition", "picolitre one-to-one droplet pairing condition",
                           MaterialKind.CONDITION)
    paired = graph.product("paired_droplets", "one-effector/one-reporter droplets",
                           MaterialKind.MOLECULAR,
                           state=molecule("co-culture droplet", (
                               GUIDE, role("perturbed cell"), role("fluorescent reporter cell"))))
    pairing_failures = graph.product("pairing_failures", "empty or incorrectly paired droplets",
                                     MaterialKind.WASTE)
    graph.react("pair_cells", "Pair one perturbed and one reporter cell per droplet",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (paired, pairing_failures))
    graph.discard(pairing_failures, "does not contain the required one-to-one cell pair")
    response_gate = graph.input("response_gate", "functional fluorescence sort gates",
                                MaterialKind.SELECTION)
    positive = graph.product("positive_droplets", "response-positive droplet fraction",
                             MaterialKind.MOLECULAR, state=paired.state)
    control = graph.product("control_droplets", "response-control droplet fraction",
                            MaterialKind.MOLECULAR, state=paired.state)
    graph.react("sort_response", "Sort droplets by functional response",
                StepKind.SELECTION,
                ((paired, InputRole.TEMPLATE), (response_gate, InputRole.SELECTION)),
                (positive, control))
    reviewed = ["cell_pools", "pairing_condition", "pair_cells", "response_gate",
                "sort_response"]

    def recover(prefix: str, template, label: str):
        extraction = graph.input(prefix + "_extraction", label + " droplet-break and DNA-extraction condition",
                                 MaterialKind.CONDITION)
        dna = graph.product(prefix + "_dna", label + " integrated guide-cassette DNA",
                            MaterialKind.MOLECULAR,
                            state=molecule(label + " guide DNA", (
                                GUIDE, role("guide-flanking genomic DNA"))))
        graph.react(prefix + "_recover", "Break droplets and recover genomic DNA",
                    StepKind.EXTRACTION,
                    ((template, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                    (dna,))
        reviewed.extend((prefix + "_extraction", prefix + "_recover"))
        _finish_pcr(graph, prefix, dna, final, label + ": " + d.endpoint,
                    reviewed, label=label + " guide-cassette")

    recover("positive", positive, "response-positive")
    recover("control", control, "control")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _disco_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    inputs = graph.input(
        "cells_and_beads", "low-input cells and mRNA-capture beads",
        MaterialKind.BIOLOGICAL,
        state=molecule("low-input suspension and capture beads", (
            role("single cells"), CELL, UMI)), primary=True)
    positioning = graph.input("vision_positioning", "machine-vision valve positioning condition",
                              MaterialKind.CONDITION)
    paired = graph.product("paired_cell_bead", "deterministically paired cell and bead",
                           MaterialKind.MOLECULAR,
                           state=molecule("deterministically paired cell and bead", (
                               role("one cell"), CELL, UMI)))
    rejected = graph.product("positioning_rejects", "empty or mis-paired valve events",
                             MaterialKind.WASTE)
    graph.react("position_pair", "Position exactly one cell and one bead",
                StepKind.PARTITION,
                ((inputs, InputRole.TEMPLATE), (positioning, InputRole.CONDITION)),
                (paired, rejected))
    graph.discard(rejected, "fails deterministic one-cell/one-bead positioning")
    droplet_condition = graph.input("droplet_condition", "co-ejection and droplet-sorting condition",
                                    MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplet", "occupied DisCo droplet",
                             MaterialKind.MOLECULAR,
                             state=molecule("occupied DisCo droplet", (
                                 role("lysed cell RNA"), CELL, UMI)))
    failed_droplets = graph.product("failed_droplets", "unsorted or failed droplets",
                                    MaterialKind.WASTE)
    graph.react("make_droplet", "Co-eject the pair and sort the occupied droplet",
                StepKind.PARTITION,
                ((paired, InputRole.TEMPLATE), (droplet_condition, InputRole.CONDITION)),
                (occupied, failed_droplets))
    graph.discard(failed_droplets, "does not reach the selected sample outlet")
    reviewed = ["cells_and_beads", "vision_positioning", "position_pair",
                "droplet_condition", "make_droplet"]
    _rt_then_pcr(graph, "library", occupied, final, d.endpoint, reviewed,
                 rt_primer_label="Drop-seq bead oligo(dT) carrying cell barcode and UMI",
                 label="DisCo")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _spindrop_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    inputs = graph.input(
        "cells_and_microgels", "fluorescently labelled cells and inDrop-v3 microgels",
        MaterialKind.BIOLOGICAL,
        state=molecule("cell and inDrop-v3 microgel", (
            role("cell or nucleus"), CELL, UMI)), primary=True)
    sorting = graph.input("fads_gate", "occupied/viable/target fluorescence gates",
                          MaterialKind.SELECTION)
    selected = graph.product("selected_droplets", "selected occupied droplets",
                             MaterialKind.MOLECULAR,
                             state=molecule("selected droplet", (
                                 role("selected cell"), CELL, UMI)))
    rejected = graph.product("rejected_droplets", "empty, damaged or unwanted droplets",
                             MaterialKind.WASTE)
    graph.react("sort_droplets", "Co-encapsulate and fluorescence-sort droplets",
                StepKind.SELECTION,
                ((inputs, InputRole.TEMPLATE), (sorting, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "fails the inline fluorescence selection gates")
    heat = graph.input("heat_treatment", "droplet heat-treatment condition",
                       MaterialKind.CONDITION)
    treated = graph.product("treated_droplets", "heat-treated selected droplets",
                            MaterialKind.MOLECULAR, state=selected.state)
    graph.react("heat_treat", "Heat-treat selected droplets", StepKind.THERMAL_TREATMENT,
                ((selected, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                (treated,), changes_molecule=False)
    reviewed = ["cells_and_microgels", "fads_gate", "sort_droplets",
                "heat_treatment", "heat_treat"]
    _rt_then_pcr(graph, "library", treated, final, d.endpoint, reviewed,
                 rt_primer_label="picoinjected improved inDrop-v3 RT/bead-primer mixture",
                 label="spinDrop")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _butt_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "nascent_rna", "purified chromatin-associated nascent RNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("purified chromatin-associated nascent RNA", (
            role("nascent RNA 5-prime body"), role("native RNA 3-prime end"))),
        primary=True)
    starter = graph.input("starter_duplex", "UMI-bearing TGIRT RNA/DNA starter duplex",
                          MaterialKind.OLIGO)
    tgirt = graph.input("tgirt", "TGIRT reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product(
        "first_strand", "nascent-RNA 3-prime-end first-strand cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("first-strand cDNA", (
            UMI, role("nascent RNA 3-prime-end cDNA"))))
    graph.react("template_switch", "Template-switch TGIRT onto the native RNA 3-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((rna, InputRole.TEMPLATE), (starter, InputRole.STRUCTURAL),
                 (tgirt, InputRole.CATALYST)), (first,))
    second_primer = graph.input("second_strand_primer", "second-strand synthesis primer",
                                MaterialKind.OLIGO)
    dna_polymerase = graph.input("second_strand_polymerase", "second-strand DNA polymerase",
                                 MaterialKind.CATALYST)
    double = graph.product(
        "double_stranded", "double-stranded Butt-seq product",
        MaterialKind.MOLECULAR,
        state=molecule("double-stranded Butt-seq product", (
            UMI, role("nascent-RNA insert"))))
    graph.react("second_strand", "Switch salt and synthesize the second strand",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (dna_polymerase, InputRole.CATALYST)), (double,))
    reviewed = ["nascent_rna", "starter_duplex", "tgirt", "template_switch",
                "second_strand_primer", "second_strand_polymerase", "second_strand"]
    forward = graph.input("pcr_forward", "Butt-seq PCR forward primer", MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "indexed Butt-seq PCR reverse primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    amplified = graph.product("amplified_library", "amplified Butt-seq library",
                              MaterialKind.MOLECULAR, state=final)
    graph.react("library_pcr", "PCR-amplify and index", StepKind.PCR,
                ((double, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The source defines the UMI and read roles without one universal complete primer geometry."))
    size_rule = graph.input("size_selection", "Butt-seq library-size selection rule",
                            MaterialKind.SELECTION)
    library = graph.product("sequencing_library", "size-selected Butt-seq library",
                            MaterialKind.LIBRARY, state=final)
    rejected = graph.product("rejected_sizes", "off-size Butt-seq products",
                             MaterialKind.WASTE)
    graph.react("select_size", "Size-select the indexed library", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (library, rejected))
    graph.discard(rejected, "falls outside the selected nascent-RNA library size")
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    reviewed.extend(("pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
                     "size_selection", "select_size", "sequencing_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _wellda_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nucleus = graph.input(
        "nucleus", "single nucleus in a well",
        MaterialKind.BIOLOGICAL,
        state=molecule("single nucleus in a well", (
            role("genomic DNA"), role("accessible chromatin"))), primary=True)
    indexed_adapter = graph.input("indexed_transposome", "well-indexed loaded Tn5 adapter",
                                  MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    marked = graph.product(
        "marked_genome", "ATAC-marked genome retaining residual genomic DNA",
        MaterialKind.MOLECULAR,
        state=molecule("ATAC-marked genome", (
            CELL, role("accessible-DNA fragments"), role("remaining genome"))))
    losses = graph.product("tn5_losses", "non-library Tn5 products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment accessible chromatin with indexed Tn5",
                StepKind.TAGMENTATION,
                ((nucleus, InputRole.TEMPLATE), (indexed_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (marked, losses))
    graph.discard(losses, "does not retain an amplifiable indexed accessibility fragment")
    atac = graph.product("atac_branch", "displaced ATAC-tagged products",
                         MaterialKind.MOLECULAR, state=marked.state)
    genome = graph.product("genome_branch", "residual genomic DNA",
                           MaterialKind.MOLECULAR, state=marked.state)
    graph.react("separate_branches", "Separate tagged products from the residual genome",
                StepKind.SPLIT, ((marked, InputRole.TEMPLATE),), (atac, genome),
                changes_molecule=False)
    reviewed = ["nucleus", "indexed_transposome", "tn5", "tagment",
                "separate_branches"]
    _finish_pcr(graph, "atac", atac, final, "ATAC branch: " + d.endpoint,
                reviewed, label="wellDA ATAC branch")
    random_primer = graph.input("genome_random_primer", "whole-genome amplification primer",
                                MaterialKind.OLIGO)
    wga_polymerase = graph.input("genome_wga_polymerase", "strand-displacing DNA polymerase",
                                 MaterialKind.CATALYST)
    amplified = graph.product("amplified_genome", "whole-genome-amplified residual genome",
                              MaterialKind.MOLECULAR, state=genome.state)
    graph.react("whole_genome_amplification", "Amplify the residual genome",
                StepKind.EXTENSION,
                ((genome, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (wga_polymerase, InputRole.CATALYST)), (amplified,))
    reviewed.extend(("genome_random_primer", "genome_wga_polymerase",
                     "whole_genome_amplification"))
    _finish_pcr(graph, "genome", amplified, final, "whole-genome branch: " + d.endpoint,
                reviewed, label="wellDA genomic branch")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _easysci_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input(
        "rt_well_nuclei", "fixed nuclei distributed into reverse-transcription wells",
        MaterialKind.BIOLOGICAL,
        state=molecule("fixed nuclei in RT wells", (role("nuclear RNA"),)), primary=True)
    dt = graph.input("barcoded_dt", "barcoded oligo(dT) RT primer", MaterialKind.OLIGO)
    random = graph.input("barcoded_random", "matched barcoded random-hexamer RT primer",
                         MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product(
        "first_index_cdna", "RT-well-indexed cDNA retained in nuclei",
        MaterialKind.MOLECULAR,
        state=molecule("first-index cDNA in nuclei", (RT_WELL, UMI, role("cDNA"))))
    graph.react("reverse_transcription", "Reverse-transcribe with barcoded dT and random primers",
                StepKind.REVERSE_TRANSCRIPTION,
                ((nuclei, InputRole.TEMPLATE), (dt, InputRole.STRUCTURAL),
                 (random, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (first,))
    repartition = graph.input("repartition_condition", "pool-and-redistribute ligation-well condition",
                              MaterialKind.CONDITION)
    redistributed = graph.product("redistributed_nuclei", "nuclei redistributed into ligation wells",
                                  MaterialKind.MOLECULAR, state=first.state)
    collisions = graph.product("repartition_losses", "empty or collision-prone redistributed wells",
                               MaterialKind.WASTE)
    graph.react("repartition", "Pool and redistribute nuclei into ligation wells",
                StepKind.PARTITION,
                ((first, InputRole.TEMPLATE), (repartition, InputRole.CONDITION)),
                (redistributed, collisions), changes_molecule=False)
    graph.discard(collisions, "does not provide an interpretable combinatorial cell identity")
    lig_barcode = graph.input("ligation_barcode", "second-round ligation-well barcode adapter",
                              MaterialKind.ADAPTER)
    ligase = graph.input("barcode_ligase", "barcode ligase", MaterialKind.CATALYST)
    indexed = graph.product(
        "two_round_nuclei", "two-round-indexed nuclei",
        MaterialKind.MOLECULAR,
        state=molecule("two-round indexed nuclei", (RT_WELL, LIG_WELL, UMI, role("cDNA")),
                       junctions=(("RT-well barcode", "ligation-well barcode", "ligation"),)))
    graph.react("barcode_ligation", "Ligate the second well barcode",
                StepKind.LIGATION,
                ((redistributed, InputRole.TEMPLATE), (lig_barcode, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (indexed,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The combinatorial barcode role is defined, but the complete plate-specific ligation-end family is not a universal sequence."))
    second_primer = graph.input("second_strand_primer", "second-strand synthesis primer",
                                MaterialKind.OLIGO)
    second_polymerase = graph.input("second_strand_polymerase", "second-strand DNA polymerase",
                                    MaterialKind.CATALYST)
    double = graph.product("double_stranded_cdna", "double-stranded two-index cDNA",
                           MaterialKind.MOLECULAR, state=indexed.state)
    graph.react("second_strand", "Synthesize the second cDNA strand",
                StepKind.EXTENSION,
                ((indexed, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (second_polymerase, InputRole.CATALYST)), (double,))
    transposome = graph.input("transposome", "loaded cDNA tagmentation adapter",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product("tagmented_cdna", "tagmented combinatorially indexed cDNA",
                              MaterialKind.MOLECULAR, state=double.state)
    losses = graph.product("tn5_losses", "non-library Tn5 products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment the double-stranded cDNA", StepKind.TAGMENTATION,
                ((double, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, losses))
    graph.discard(losses, "does not retain an amplifiable two-ended tagmentation product")
    reviewed = [
        "rt_well_nuclei", "barcoded_dt", "barcoded_random", "reverse_transcriptase",
        "reverse_transcription", "repartition_condition", "repartition",
        "ligation_barcode", "barcode_ligase", "barcode_ligation",
        "second_strand_primer", "second_strand_polymerase", "second_strand",
        "transposome", "tn5", "tagment",
    ]
    _finish_pcr(graph, "library", tagmented, final, d.endpoint, reviewed,
                label="EasySci-RNA PCR-well-indexed")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _indrops2_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    droplet = graph.input(
        "droplet", "cell and barcoded hydrogel bead in a droplet",
        MaterialKind.BIOLOGICAL,
        state=molecule("cell and barcoded hydrogel bead", (
            role("poly(A) RNA"), CELL, UMI)), primary=True)
    bead_primer = graph.input("bead_primer", "barcoded bead oligo(dT) RT primer",
                              MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase",
                     MaterialKind.CATALYST)
    first = graph.product(
        "barcoded_cdna", "barcoded first-strand cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("barcoded first-strand cDNA", (
            CELL, UMI, role("3-prime transcript cDNA"), role("amplification handle"))))
    graph.react("reverse_transcription", "Lyse the cell and reverse-transcribe on the bead oligo",
                StepKind.REVERSE_TRANSCRIPTION,
                ((droplet, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    ivt_input = graph.product("ivt_input", "cDNA assigned to linear IVT amplification",
                              MaterialKind.MOLECULAR, state=first.state)
    tso_input = graph.product("tso_input", "cDNA assigned to template-switch PCR amplification",
                              MaterialKind.MOLECULAR, state=first.state)
    graph.react("choose_branch", "Retain the two validated amplification alternatives",
                StepKind.SPLIT, ((first, InputRole.TEMPLATE),), (ivt_input, tso_input),
                changes_molecule=False)
    reviewed = ["droplet", "bead_primer", "reverse_transcriptase",
                "reverse_transcription", "choose_branch"]

    t7 = graph.input("ivt_t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    amplified_rna = graph.product(
        "ivt_rna", "linearly amplified inDrops-2 RNA",
        MaterialKind.MOLECULAR,
        state=molecule("IVT-amplified RNA", (CELL, UMI, role("transcript RNA copy"))))
    graph.react("ivt", "Linearly amplify the T7-bearing cDNA by IVT",
                StepKind.TRANSCRIPTION,
                ((ivt_input, InputRole.TEMPLATE), (t7, InputRole.CATALYST)),
                (amplified_rna,))
    ivt_rt_primer = graph.input("ivt_rt_primer", "IVT-product reverse-transcription primer",
                                MaterialKind.OLIGO)
    ivt_rt = graph.input("ivt_reverse_transcriptase", "IVT-product reverse transcriptase",
                         MaterialKind.CATALYST)
    ivt_cdna = graph.product("ivt_cdna", "IVT-amplified inDrops-2 cDNA",
                             MaterialKind.MOLECULAR, state=first.state)
    graph.react("ivt_back_rt", "Convert amplified RNA back to cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((amplified_rna, InputRole.TEMPLATE), (ivt_rt_primer, InputRole.STRUCTURAL),
                 (ivt_rt, InputRole.CATALYST)), (ivt_cdna,))
    reviewed.extend(("ivt_t7_polymerase", "ivt", "ivt_rt_primer",
                     "ivt_reverse_transcriptase", "ivt_back_rt"))

    tso = graph.input("tso", "inDrops-2 template-switch oligo", MaterialKind.OLIGO)
    tso_rt = graph.input("tso_reverse_transcriptase", "template-switching reverse transcriptase",
                         MaterialKind.CATALYST)
    switched = graph.product("tso_cdna", "template-switched inDrops-2 cDNA",
                             MaterialKind.MOLECULAR, state=first.state)
    graph.react("template_switch", "Install the template-switch PCR handle",
                StepKind.TEMPLATE_SWITCHING,
                ((tso_input, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (tso_rt, InputRole.CATALYST)), (switched,))
    amp_forward = graph.input("tso_amp_forward", "template-switch amplification forward primer",
                              MaterialKind.OLIGO)
    amp_reverse = graph.input("tso_amp_reverse", "template-switch amplification reverse primer",
                              MaterialKind.OLIGO)
    amp_polymerase = graph.input("tso_amp_polymerase", "cDNA amplification polymerase",
                                 MaterialKind.CATALYST)
    tso_cdna = graph.product("tso_amplified_cdna", "PCR-amplified inDrops-2 cDNA",
                             MaterialKind.MOLECULAR, state=switched.state)
    graph.react("tso_amplification", "PCR-amplify the template-switched branch",
                StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (amp_forward, InputRole.STRUCTURAL),
                 (amp_reverse, InputRole.STRUCTURAL), (amp_polymerase, InputRole.CATALYST)),
                (tso_cdna,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The open-source implementation permits primer/index choices; no one immutable amplification pair is universal."))
    reviewed.extend(("tso", "tso_reverse_transcriptase", "template_switch",
                     "tso_amp_forward", "tso_amp_reverse", "tso_amp_polymerase",
                     "tso_amplification"))

    def finish(prefix: str, template, label: str):
        transposome = graph.input(prefix + "_transposome", label + " loaded library transposome",
                                  MaterialKind.ADAPTER)
        tn5 = graph.input(prefix + "_tn5", label + " Tn5 transposase",
                          MaterialKind.CATALYST)
        tagmented = graph.product(prefix + "_tagmented", label + " tagmented cDNA",
                                  MaterialKind.MOLECULAR, state=template.state)
        losses = graph.product(prefix + "_losses", label + " non-library Tn5 products",
                               MaterialKind.WASTE)
        graph.react(prefix + "_tagment", "Tagment the " + label + " product",
                    StepKind.TAGMENTATION,
                    ((template, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                     (tn5, InputRole.CATALYST)), (tagmented, losses))
        graph.discard(losses, "does not retain an amplifiable barcode-bearing insert")
        reviewed.extend((prefix + "_transposome", prefix + "_tn5", prefix + "_tagment"))
        _finish_pcr(graph, prefix, tagmented, final, label + ": " + d.endpoint,
                    reviewed, label=label)

    finish("ivt_library", ivt_cdna, "inDrops-2 IVT branch")
    finish("tso_library", tso_cdna, "inDrops-2 TSO-PCR branch")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _multitag_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "cells", "permeabilized cells containing several chromatin targets",
        MaterialKind.BIOLOGICAL,
        state=molecule("permeabilized cell", (
            role("target 1"), role("target 2"), role("target 3"))), primary=True)
    reviewed = ["cells"]

    def targeting_round(prefix: str, template, label: str):
        antibody = graph.input(prefix + "_antibody", label + " primary antibody",
                               MaterialKind.SELECTION)
        bound = graph.product(prefix + "_bound", label + " antibody-bound chromatin",
                              MaterialKind.MOLECULAR, state=template.state)
        unbound = graph.product(prefix + "_unbound", label + " unbound chromatin",
                                MaterialKind.WASTE)
        graph.react(prefix + "_bind", "Bind " + label + " primary antibody",
                    StepKind.BINDING,
                    ((template, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                    (bound, unbound))
        graph.discard(unbound, "does not retain the selected antibody-target complex")
        transposome = graph.input(prefix + "_transposome", label + " barcode-loaded pA–Tn5",
                                  MaterialKind.ADAPTER)
        tn5 = graph.input(prefix + "_tn5", label + " pA–Tn5 transposase",
                          MaterialKind.CATALYST)
        tagged = graph.product(prefix + "_tagged", label + " target-barcoded chromatin",
                               MaterialKind.MOLECULAR,
                               state=molecule(label + " tagged chromatin", (
                                   TARGET, role("target-proximal genomic DNA"))))
        losses = graph.product(prefix + "_losses", label + " non-library Tn5 products",
                               MaterialKind.WASTE)
        graph.react(prefix + "_tagment", "Tagment " + label + " target-proximal DNA",
                    StepKind.TAGMENTATION,
                    ((bound, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                     (tn5, InputRole.CATALYST)), (tagged, losses))
        graph.discard(losses, "does not retain the target-specific adapter identity")
        reviewed.extend((prefix + "_antibody", prefix + "_bind",
                         prefix + "_transposome", prefix + "_tn5", prefix + "_tagment"))
        return tagged

    first = targeting_round("first", cells, "first target")
    multiplex = targeting_round("additional", first, "additional-target pool")
    partition = graph.input("cell_index_partition", "single-cell combinatorial-index partitioning condition",
                            MaterialKind.CONDITION)
    indexed = graph.product("indexed_cells", "combinatorially indexed target-tagged cells",
                            MaterialKind.MOLECULAR, state=multiplex.state)
    collisions = graph.product("index_collisions", "ambiguous cell-index partitions",
                               MaterialKind.WASTE)
    graph.react("index_cells", "Combinatorially index single cells",
                StepKind.PARTITION,
                ((multiplex, InputRole.TEMPLATE), (partition, InputRole.CONDITION)),
                (indexed, collisions), changes_molecule=False)
    graph.discard(collisions, "does not yield an unambiguous cell-barcode combination")
    reviewed.extend(("cell_index_partition", "index_cells"))
    _finish_pcr(graph, "library", indexed, final, d.endpoint, reviewed,
                label="MulTI-Tag")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _phage_atac_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "stained_cells", "cells stained with nanobody-displaying M13 phage",
        MaterialKind.BIOLOGICAL,
        state=molecule("phage-stained cell", (
            role("surface or intracellular epitope"), TARGET, role("PAC-tag phagemid"))),
        primary=True)
    transposome = graph.input("atac_transposome", "loaded ATAC transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "ATAC Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product(
        "tagmented_cells", "phage-bound cells with tagmented accessible chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("phage-bound tagmented cell", (
            TARGET, role("PAC-tag phagemid"), role("ATAC fragments"))))
    losses = graph.product("tn5_losses", "non-library ATAC products", MaterialKind.WASTE)
    graph.react("tagment", "Fix, lyse and tagment accessible chromatin",
                StepKind.TAGMENTATION,
                ((cells, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, losses))
    graph.discard(losses, "does not retain an amplifiable accessible-chromatin fragment")
    droplets = graph.input("droplet_condition", "10x ATAC gel-bead co-encapsulation condition",
                           MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplets", "cells co-encapsulated with ATAC gel beads",
                             MaterialKind.MOLECULAR, state=tagmented.state)
    empty = graph.product("empty_droplets", "empty or collision droplets", MaterialKind.WASTE)
    graph.react("partition", "Co-encapsulate cells and ATAC gel beads",
                StepKind.PARTITION,
                ((tagmented, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (occupied, empty), changes_molecule=False)
    graph.discard(empty, "does not contain one interpretable cell/bead combination")
    bead_primer = graph.input("bead_primer", "10x ATAC cell-barcode bead primer",
                              MaterialKind.OLIGO)
    polymerase = graph.input("linear_polymerase", "droplet linear-amplification polymerase",
                             MaterialKind.CATALYST)
    barcoded = graph.product(
        "barcoded_products", "cell-barcoded ATAC and phagemid products",
        MaterialKind.MOLECULAR,
        state=molecule("droplet-barcoded products", (
            CELL, TARGET, role("phagemid tag"), role("ATAC fragment"))))
    graph.react("linear_amplification", "Apply the droplet barcode to chromatin and phage tags",
                StepKind.EXTENSION,
                ((occupied, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (barcoded,))
    atac = graph.product("atac_branch", "cell-barcoded ATAC products",
                         MaterialKind.MOLECULAR, state=barcoded.state)
    protein = graph.product("protein_branch", "cell-barcoded protein-derived phage tags",
                            MaterialKind.MOLECULAR, state=barcoded.state)
    graph.react("split_readouts", "Split ATAC and protein-derived-tag libraries",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (atac, protein),
                changes_molecule=False)
    reviewed = [
        "stained_cells", "atac_transposome", "tn5", "tagment",
        "droplet_condition", "partition", "bead_primer", "linear_polymerase",
        "linear_amplification", "split_readouts",
    ]
    _finish_pcr(graph, "atac", atac, final, "ATAC branch: " + d.endpoint,
                reviewed, label="PHAGE-ATAC chromatin")
    _finish_pcr(graph, "protein", protein, final, "protein-tag branch: " + d.endpoint,
                reviewed, label="PHAGE-ATAC protein tag")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _sctem_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input(
        "cell", "one sorted cell containing genomic DNA and mRNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("one sorted cell", (role("genomic DNA"), role("mRNA"))),
        primary=True)
    dna = graph.product("dna_branch", "single-cell genomic DNA branch",
                        MaterialKind.MOLECULAR, state=cell.state)
    rna = graph.product("rna_branch", "matched single-cell mRNA branch",
                        MaterialKind.MOLECULAR, state=cell.state)
    graph.react("separate_modalities", "Separate genomic DNA and mRNA",
                StepKind.SPLIT, ((cell, InputRole.TEMPLATE),), (dna, rna),
                changes_molecule=False)
    bisulfite = graph.input("bisulfite", "bisulfite conversion condition",
                            MaterialKind.CONDITION)
    converted = graph.product(
        "converted_dna", "bisulfite-converted transposable-element templates",
        MaterialKind.MOLECULAR,
        state=molecule("converted DNA", (
            role("converted SINE Alu / LINE-1 templates"),)))
    graph.react("bisulfite_conversion", "Bisulfite-convert genomic DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((dna, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    te_forward = graph.input("te_forward", "staggered internally indexed TE forward-primer pool",
                             MaterialKind.OLIGO)
    te_reverse = graph.input("te_reverse", "staggered internally indexed TE reverse-primer pool",
                             MaterialKind.OLIGO)
    te_polymerase = graph.input("te_polymerase", "transposable-element PCR polymerase",
                                MaterialKind.CATALYST)
    dna_library = graph.product("dna_library", "scTEM transposable-element amplicon library",
                                MaterialKind.LIBRARY, state=final)
    graph.react("te_pcr", "Amplify converted transposable elements",
                StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (te_forward, InputRole.STRUCTURAL),
                 (te_reverse, InputRole.STRUCTURAL), (te_polymerase, InputRole.CATALYST)),
                (dna_library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The panel uses staggered internally indexed primer families rather than one universal amplicon."))
    graph.sequence(dna_library, SequencingPlan(limitation="DNA branch: " + d.endpoint))
    rna_final = role_library(
        "scTEM matched transcriptome library roles",
        (CELL, role("full-length transcript cDNA")), index=SAMPLE)
    reviewed = [
        "cell", "separate_modalities", "bisulfite", "bisulfite_conversion",
        "te_forward", "te_reverse", "te_polymerase", "te_pcr", "dna_library",
    ]
    _rt_then_pcr(graph, "rna", rna, rna_final,
                 "matched RNA branch: Smart-seq-type full-length transcript reads",
                 reviewed, rt_primer_label="Smart-seq-type oligo(dT) RT primer",
                 label="scTEM transcriptome")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _limca_graph(d: Definition, final) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input(
        "cell", "one intact cell",
        MaterialKind.BIOLOGICAL,
        state=molecule("one intact cell", (
            role("cytoplasmic mRNA"), role("nuclear chromatin"))), primary=True)
    cytoplasm = graph.product("cytoplasm", "cytoplasmic RNA compartment",
                              MaterialKind.MOLECULAR, state=cell.state)
    nucleus = graph.product("nucleus", "matched nuclear chromatin compartment",
                            MaterialKind.MOLECULAR, state=cell.state)
    graph.react("separate_compartments", "Physically separate cytoplasm and nucleus",
                StepKind.SPLIT, ((cell, InputRole.TEMPLATE),), (cytoplasm, nucleus),
                changes_molecule=False)
    reviewed = ["cell", "separate_compartments"]
    rna_final = role_library("LiMCA Smart-seq2 transcriptome library roles",
                             (CELL, role("full-length transcript cDNA")), index=SAMPLE)
    _rt_then_pcr(graph, "rna", cytoplasm, rna_final,
                 "RNA branch: " + d.endpoint, reviewed,
                 rt_primer_label="Smart-seq2 oligo(dT) primer",
                 label="LiMCA Smart-seq2")

    formaldehyde = graph.input("formaldehyde", "formaldehyde cross-linking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed_nucleus", "crosslinked nuclear chromatin",
                          MaterialKind.MOLECULAR, state=nucleus.state)
    graph.react("crosslink", "Crosslink nuclear chromatin", StepKind.CROSS_LINKING,
                ((nucleus, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)),
                (fixed,))
    restriction = graph.input("restriction_enzyme", "Hi-C restriction endonuclease",
                              MaterialKind.CATALYST)
    digested = graph.product("restriction_fragments", "restriction-digested fixed chromatin",
                             MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("digest", "Restriction-digest the nucleus", StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)),
                (digested,))
    ligase = graph.input("proximity_ligase", "in-nucleus proximity DNA ligase",
                         MaterialKind.CATALYST)
    contacts = graph.product(
        "contacts", "single-cell proximity-ligation junctions",
        MaterialKind.MOLECULAR,
        state=molecule("single-cell contact junctions", (
            role("restriction fragment A"), role("restriction fragment B")),
            junctions=(("restriction fragment A", "restriction fragment B",
                        "proximity ligation"),)))
    graph.react("proximity_ligation", "Proximity-ligate restriction fragments",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
                (contacts,))
    reversal = graph.input("crosslink_reversal", "cross-link reversal condition",
                           MaterialKind.CONDITION)
    purified = graph.product("purified_contacts", "purified single-cell contact DNA",
                             MaterialKind.MOLECULAR, state=contacts.state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover contact DNA",
                StepKind.DEPROTEINIZATION,
                ((contacts, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (purified,))
    meta_primer = graph.input("meta_primer", "META whole-genome amplification primer",
                              MaterialKind.OLIGO)
    meta_polymerase = graph.input("meta_polymerase", "META strand-displacing polymerase",
                                  MaterialKind.CATALYST)
    amplified = graph.product("amplified_contacts", "META-amplified contact DNA",
                              MaterialKind.MOLECULAR, state=purified.state)
    graph.react("meta_wga", "Apply META whole-genome amplification",
                StepKind.EXTENSION,
                ((purified, InputRole.TEMPLATE), (meta_primer, InputRole.STRUCTURAL),
                 (meta_polymerase, InputRole.CATALYST)), (amplified,))
    reviewed.extend((
        "formaldehyde", "crosslink", "restriction_enzyme", "digest",
        "proximity_ligase", "proximity_ligation", "crosslink_reversal",
        "reverse_crosslinks", "meta_primer", "meta_polymerase", "meta_wga",
    ))
    _finish_pcr(graph, "hic", amplified, final, "Hi-C branch: " + d.endpoint,
                reviewed, label="LiMCA Hi-C")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


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
    spec = SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE="The source defines the read roles, but the primer sequence is panel-, kit- or implementation-dependent; no universal binding oligo is asserted.",
        workflow=WORKFLOWS[key], oligos=lambda: oligos,
    )
    factories = {
        "wellda": _wellda_graph,
        "easysci": _easysci_graph,
        "indrops2": _indrops2_graph,
        "multitag": _multitag_graph,
        "phage-atac": _phage_atac_graph,
        "sctem": _sctem_graph,
        "end5": _end5_graph,
        "speac": _speac_graph,
        "disco": _disco_graph,
        "spindrop": _spindrop_graph,
        "limca": _limca_graph,
        "buttseq": _butt_graph,
    }
    if key in factories:
        spec.protocol_graph = lambda: factories[key](d, final)
    return spec
