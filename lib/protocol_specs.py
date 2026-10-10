"""Source-backed specifications for the named-protocol expansion batch.

This module is deliberately data-oriented.  Protocol directories remain separate public
units, while shared constructors make the recurring molecular operations consistent.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Callable

import seqprimers as sp
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature
from duplex_error import codec_linked_duplex, nanoseq_restrictive_end_repair
from interaction_capture import ProteinMediatedRNAContact, ReversibleDuplexCrosslink


def role(name: str, width: int = 18, tag: str | None = None, *, length_bp=None,
         inferred: bool = False, molecular_feature=None) -> Segment:
    token = f"[{length_bp}-bp {name}]" if length_bp else f"[{name}]"
    return Segment(name, token, tag, placeholder=True, length_bp=length_bp,
                   inferred=inferred, feature=molecular_feature)


def molecule(label: str, names: tuple[str, ...], *, joined: bool = False,
             note: str = "", junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    parts = [Segment(name, chr(88 + (i % 2)) * 18, placeholder=True)
             for i, name in enumerate(names)]
    sc = Scene(); sc.strand("molecule", parts, label=label)
    # ``joined`` means only that the labels describe one covalent molecule.  A barcode
    # copied from a primer is just as joined as a ligated adapter, but it is not a
    # ligation site.  Only callers naming an exact boundary may draw the ** marker.
    for left, right, label_text in junctions:
        sc.junction("molecule", left, right, label_text)
    if note:
        sc.note("molecule", note)
    sc.labels("molecule")
    return sc


def two_branch(label: str, a: str, b: str, note: str = "") -> Scene:
    sc = Scene(); sc.strand("branch A", [role(a)], label=label)
    sc.strand("branch B", [role(b)], label=label)
    if note:
        sc.note("branch A", note)
    return sc


def workflow(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        w.react(action, state.rows(), name=action, note=note)
    return w


def _psoralen(key: str) -> Workflow:
    reagent, selection = {
        "paris": ("AMT", "cross-linked duplexes selected above the diagonal on a 2D gel"),
        "splash": ("biotinylated psoralen", "streptavidin selects biotin-psoralen RNA duplexes"),
        "ligr-seq": ("AMT", "RNase R depletes uncross-linked linear RNA"),
    }[key]
    x = ReversibleDuplexCrosslink(reagent, selection)
    return workflow(molecule("cellular RNA", ("RNA molecules",)), (
        (f"Cross-link paired RNA with {reagent} and long-wave UV", x.duplex(), selection),
        ("Fragment RNA and enrich cross-linked duplexes", x.duplex(), selection),
        ("Proximity-ligate the two RNA partners", x.chimera(),
         f"The covalent ligation is made while {reagent} holds the duplex together."),
        (f"Photoreverse {reagent}, then install sequencing adapters",
         molecule("cDNA library precursor", ("RNA A-derived cDNA", "RNA B-derived cDNA"), joined=True),
         "The chimeric junction, rather than the crosslink, remains in the sequenced molecule."),
    ))


def workflow_for(key: str) -> Workflow:
    if key in {"paris", "splash", "ligr-seq"}:
        return _psoralen(key)
    if key == "mario":
        x = ProteinMediatedRNAContact("biotinylated RNA linker", True)
        return workflow(molecule("cellular RNPs", ("protein-bound RNAs",)), (
            ("UV-cross-link RNA to its bound proteins", x.complex(),
             "Protein, rather than RNA base pairing, preserves the contact."),
            ("Attach the biotinylated linker and proximity-ligate co-bound RNAs", x.chimera(),
             "The linker creates an explicit RNA1–linker–RNA2 boundary."),
            ("Select linker-containing chimeras and convert them to cDNA", x.chimera(),
             "Streptavidin enrichment precedes paired-end library preparation."),
        ))
    if key == "codec":
        return workflow(molecule("original genomic duplex", ("Watson strand", "Crick strand")), (
            ("End-repair and ligate the preassembled CODEC adapter quadruplex",
             molecule("adapter-bound duplex", ("left CODEC adapter", "original duplex", "right CODEC adapter"), joined=True,
                      junctions=(("left CODEC adapter", "original duplex", "adapter ligation"),
                                 ("original duplex", "right CODEC adapter", "adapter ligation"))),
             "INFERRED — the public article establishes the roles; exact oligo sequence is supplementary."),
            ("Strand-displacing extension copies and physically links both complements",
             codec_linked_duplex(), "Both original strands become recoverable in one cluster molecule."),
            ("PCR and outward Illumina sequencing", codec_linked_duplex(),
             "Reads recover insert, UMI and index information for within-pair consensus."),
        ))
    if key == "nanoseq":
        repaired = nanoseq_restrictive_end_repair()
        return workflow(molecule("high-molecular-weight genomic DNA", ("genomic DNA",)), (
            ("Digest with a restriction enzyme instead of random shearing", repaired,
             "Defined intact ends limit lesion-generating end repair."),
            ("Perform restrictive end repair with dideoxynucleotides", repaired,
             "Nick extension is terminated before it can copy damage into both strands."),
            ("Ligate duplex adapters, amplify and recover both strand families",
             molecule("duplex-consensus library", ("duplex tag A", "insert", "duplex tag B"), joined=True,
                      junctions=(("duplex tag A", "insert", "adapter ligation"),
                                 ("insert", "duplex tag B", "adapter ligation"))),
             "Complementary tags associate Watson and Crick read families."),
        ))
    if key == "stlfr":
        return workflow(molecule("20–300-kb DNA", ("long genomic molecule",)), (
            ("Insert transposons carrying capture sequences along each long molecule",
             molecule("tagged long DNA", ("subfragment", "capture transposon", "subfragment")),
             "Transposition creates many addressable sites without losing long-range provenance."),
            ("Capture one long molecule on one barcoded bead",
             molecule("bead-bound molecule", ("bead barcode", "captured subfragments")),
             "Millions of bead barcodes provide practically nonredundant co-barcoding."),
            ("Ligate the bead barcode to captured subfragments and make DNBs",
             molecule("co-barcoded library", ("linked-read barcode", "short genomic subfragment"), joined=True,
                      junctions=(("linked-read barcode", "short genomic subfragment", "barcode ligation"),)),
             "Reads sharing a barcode are assigned back to the same long molecule."),
        ))
    if key == "nano3p-seq":
        return workflow(molecule("native RNA", ("RNA body", "native 3-prime end")), (
            ("Prime at the native RNA 3′ end", molecule("RNA–primer complex", ("RNA", "3-prime RT primer")),
             "Capture does not require a poly(A) tail or RNA-adapter ligation."),
            ("Template-switch during reverse transcription",
             molecule("full cDNA", ("template-switch handle", "RNA-derived cDNA"), joined=True),
             "The RT places a sequencing-compatible handle at the opposite cDNA end."),
            ("Attach the Nanopore sequencing adapter without PCR",
             molecule("Nanopore cDNA molecule", ("motor adapter", "full-length cDNA"), joined=True,
                      junctions=(("motor adapter", "full-length cDNA", "adapter ligation"),)),
             "One molecule reports RNA abundance and native 3′-tail composition."),
        ))
    if key == "maester":
        return workflow(molecule("barcoded 3′ scRNA-seq cDNA", ("cell barcode", "UMI", "cDNA")), (
            ("Multiplex-PCR mitochondrial transcript targets",
             molecule("targeted amplicon", ("cell barcode", "UMI", "mitochondrial target"), joined=True),
             "Target enrichment preserves the original cell barcode and UMI."),
            ("Pool target reactions and add Illumina indices",
             molecule("MAESTER library", ("sample index", "cell barcode", "UMI", "mtRNA target"), joined=True),
             "The enriched library is sequenced separately from the parent expression library."),
        ))
    if key == "redeem":
        return workflow(two_branch("fixed permeabilized whole cell", "nuclear chromatin", "retained mtDNA"), (
            ("Run 10x Multiome while retaining mitochondrial DNA",
             two_branch("matched cell barcode", "RNA + ATAC libraries", "mtDNA fragments"),
             "All branches retain a common droplet cell identity."),
            ("Hybrid-capture mtDNA with tiled probes",
             molecule("captured mtDNA branch", ("cell barcode", "mtDNA fragment", "capture probe")),
             "Fragment endpoints act as endogenous molecular identifiers."),
            ("Sequence mtDNA, ATAC and RNA as three matched libraries",
             two_branch("same cell barcode", "mtDNA library", "ATAC + RNA libraries"),
             "The mtDNA branch receives deep error-corrected sequencing."),
        ))
    if key == "rage-seq":
        return workflow(molecule("droplet scRNA-seq cDNA", ("cell barcode", "UMI", "full-length cDNA")), (
            ("Hybrid-capture immunoglobulin and T-cell-receptor cDNA",
             molecule("captured receptor cDNA", ("cell barcode", "UMI", "full-length V(D)J transcript")),
             "Biotinylated probes enrich receptor molecules while retaining cellular identifiers."),
            ("Attach long-read adapters and sequence the captured branch",
             molecule("RAGE-Seq long-read molecule", ("Nanopore adapter", "cell barcode", "full-length receptor cDNA"), joined=True,
                      junctions=(("Nanopore adapter", "cell barcode", "adapter ligation"),)),
             "Short-read gene expression and long-read receptor data join through the barcode."),
        ))
    if key == "malbac-dt":
        return workflow(molecule("single-cell poly(A) RNA", ("transcripts",)), (
            ("Reverse-transcribe with a cell barcode, UMI and poly(dT)",
             molecule("barcoded first-strand cDNA", ("cell barcode", "UMI", "cDNA"), joined=True),
             "The well barcode and UMI are installed before pooling."),
            ("Perform MALBAC quasi-linear random-primer amplification",
             molecule("MALBAC-DT amplicons", ("common handle", "barcoded cDNA", "common handle")),
             "Looping suppresses repeated copying of early amplicons."),
            ("Pool, tagment and add sequencing indices",
             molecule("indexed library", ("left Nextera end", "cell barcode", "UMI", "cDNA", "right Nextera end"), joined=True),
             "Read 1 reports barcode and UMI; Read 2 reports transcript sequence."),
        ))
    if key == "seq-scope":
        return workflow(molecule("Illumina flow cell", ("random HDMI oligos",)), (
            ("Cluster and sequence the HDMI array in situ",
             molecule("decoded capture surface", ("mapped HDMI cluster", "poly(dT) capture end")),
             "The first sequencing run maps each spatial barcode to an imaging coordinate."),
            ("Place tissue, capture mRNA and reverse-transcribe on the flow cell",
             molecule("surface cDNA", ("spatial HDMI", "UMI", "transcript cDNA"), joined=True),
             "Each cDNA inherits the barcode of the cluster beneath it."),
            ("Release, amplify and index the spatial transcript library",
             molecule("Seq-Scope library", ("sample index", "HDMI", "UMI", "transcript"), joined=True),
             "Sequencing reconnects transcript identity with the decoded surface coordinate."),
        ))
    if key == "chil-seq":
        return workflow(molecule("fixed chromatin", ("antibody target", "nearby DNA")), (
            ("Bind an antibody carrying a ChIL oligonucleotide",
             molecule("targeted chromatin", ("antibody–ChIL probe", "nearby genomic DNA")),
             "The probe contains a T7 promoter and sequence used to recover nearby DNA."),
            ("Tagment nearby genomic DNA in situ",
             molecule("tagged target DNA", ("ChIL probe", "target-proximal genomic DNA", "Tn5 end")),
             "The antibody position and the transposed DNA are joined during subsequent copying."),
            ("Extend, transcribe linearly with T7, then make a cDNA library",
             molecule("ChIL-seq cDNA", ("ChIL-probe end", "target-proximal DNA copy"), joined=True),
             "T7 transcription amplifies the selected genomic neighbourhood before PCR."),
        ))
    if key == "bliss":
        return workflow(molecule("fixed DNA", ("native double-strand break",)), (
            ("Blunt the break ends in situ",
             molecule("blunted break", ("genomic DNA", "blunt DSB end")),
             "Only ends already present in fixed material are exposed for tagging."),
            ("Ligate a UMI, sample barcode, RA5 and T7-promoter adapter",
             molecule("tagged break", ("genomic break", "UMI", "sample barcode", "RA5", "T7 promoter"), joined=True,
                      junctions=(("genomic break", "UMI", "break-adapter ligation"),)),
             "The ligation boundary marks the nucleotide adjacent to the break."),
            ("Linearly amplify by T7 transcription and complete the RNA library",
             molecule("BLISS library", ("RA5", "UMI", "sample barcode", "break-adjacent DNA"), joined=True),
             "IVT reduces dependence on exponential preamplification."),
        ))
    if key == "scsplat":
        return workflow(molecule("one cell", ("genomic DNA",)), (
            ("Bisulfite-convert DNA separately in each well",
             molecule("single-stranded converted DNA", ("bisulfite-converted fragments",)),
             "Conversion precedes adapter installation."),
            ("Random-prime second-strand synthesis with P5 and an inline cell barcode",
             molecule("barcoded strand", ("P5", "cell barcode", "random-primer copy"), joined=True),
             "Up to 32 individually barcoded wells can now be pooled."),
            ("Pool cells and install the other adapter by splinted ligation",
             molecule("scSPLAT library", ("P5", "cell barcode", "converted insert", "splint-ligated P7"), joined=True,
                      junctions=(("converted insert", "splint-ligated P7", "splinted ligation"),)),
             "The splint makes adapter tagging a bulk reaction after cell indexing."),
        ))
    if key == "cfmedip-seq":
        return workflow(molecule("plasma cfDNA", ("naturally fragmented cfDNA",)), (
            ("End-repair, dA-tail and ligate Illumina adapters",
             molecule("adapter-ligated cfDNA", ("P5-side adapter", "cfDNA", "P7-side adapter"), joined=True,
                      junctions=(("P5-side adapter", "cfDNA", "adapter ligation"),
                                 ("cfDNA", "P7-side adapter", "adapter ligation"))),
             "Adapters are installed before methylated-fragment enrichment."),
            ("Add methylated and unmethylated spike-ins plus λ filler DNA",
             molecule("low-input library mixture", ("cfDNA library", "controls", "filler DNA")),
             "Filler supports low-input immunoprecipitation; controls measure recovery."),
            ("Denature and immunoprecipitate 5-methylcytosine-containing fragments",
             molecule("5mC-enriched library", ("adapter", "methylated cfDNA fragment", "adapter"), joined=True),
             "The assay enriches methylated fragments without bisulfite conversion."),
        ))
    if key == "dart-seq":
        return workflow(molecule("cellular RNA", ("RNA containing m6A",)), (
            ("Express APOBEC1 fused to the YTH m6A-binding domain",
             molecule("targeted RNA", ("APOBEC1–YTH", "m6A-adjacent cytidine")),
             "YTH recruits cytidine deamination near m6A."),
            ("Convert nearby cytidine to uridine in the cell",
             molecule("marked RNA", ("RNA", "C-to-U mark near m6A")),
             "The modification is read as a C-to-T mismatch after reverse transcription."),
            ("Prepare an ordinary RNA-seq library",
             molecule("DART-seq library", ("left adapter", "cDNA carrying C-to-T mark", "right adapter"), joined=True,
                      junctions=(("left adapter", "cDNA carrying C-to-T mark", "adapter ligation"),
                                 ("cDNA carrying C-to-T mark", "right adapter", "adapter ligation"))),
             "The library reports an enzymatic mark, not direct sequencing of m6A."),
        ))
    if key == "snubar-atac":
        return workflow(molecule("sample-specific nuclei", ("chromatin",)), (
            ("Tagment each sample with the standard 10x ATAC transposome",
             molecule("tagmented chromatin", ("genomic fragment", "transposome universal tail")),
             "The sample barcode is not carried by Tn5 itself."),
            ("Anneal a SNuBar oligo to the transposome tail",
             molecule("sample-tagged fragment", ("genomic fragment", "sample barcode", "UMI", "PCR handle"), joined=True),
             "The SBO transfers sample identity before samples are pooled."),
            ("Pool samples and run the 10x single-cell ATAC workflow",
             molecule("SNuBar-ATAC library", ("10x cell barcode", "sample barcode", "UMI", "ATAC fragment"), joined=True),
             "Cell and sample identities occupy distinct molecular fields."),
        ))
    if key == "ric-seq":
        return workflow(molecule("fixed cells", ("spatially proximal RNAs",)), (
            ("Cross-link RNA–protein assemblies in situ and fragment RNA",
             two_branch("fixed RNP neighbourhood", "RNA fragment A", "RNA fragment B"),
             "Formaldehyde preserves spatial contacts, including protein-mediated ones."),
            ("Label RNA ends and proximity-ligate them in situ",
             molecule("RIC-seq chimera", ("RNA A fragment", "RNA B fragment"), joined=True,
                      junctions=(("RNA A fragment", "RNA B fragment", "RNA proximity ligation"),)),
             "The ligation junction records a spatial RNA–RNA contact."),
            ("Reverse crosslinks, enrich chimeras and make a paired-end library",
             molecule("RIC-seq library", ("left adapter", "RNA A cDNA", "RNA B cDNA", "right adapter"), joined=True,
                      junctions=(("left adapter", "RNA A cDNA", "adapter ligation"),
                                 ("RNA A cDNA", "RNA B cDNA", "copied proximity-ligation junction"),
                                 ("RNA B cDNA", "right adapter", "adapter ligation"))),
             "Paired sequence identifies both sides of the RNA contact."),
        ))
    if key == "plac-seq":
        return workflow(molecule("cross-linked nuclei", ("chromatin contacts",)), (
            ("Restriction-digest chromatin and fill ends with biotinylated nucleotides",
             two_branch("one chromatin contact", "filled locus A end", "filled locus B end"),
             "Biotin marks the ends that can become contact junctions."),
            ("Proximity-ligate inside intact nuclei",
             molecule("contact product", ("locus A", "locus B"), joined=True,
                      junctions=(("locus A", "locus B", "biotinylated proximity ligation"),)),
             "The junction preserves the same in-nucleus contact geometry as in situ Hi-C."),
            ("Shear, immunoprecipitate the chosen chromatin target, and enrich biotin",
             molecule("protein-anchored contact library", ("left Illumina adapter", "locus A", "locus B", "right Illumina adapter"), joined=True,
                      junctions=(("left Illumina adapter", "locus A", "adapter ligation"),
                                 ("locus A", "locus B", "contact ligation"),
                                 ("locus B", "right Illumina adapter", "adapter ligation"))),
             "ChIP occurs after proximity ligation, defining the distinction from HiChIP."),
        ))
    raise KeyError(key)


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    source_label: str
    summary: str
    library: str
    payload: tuple[str, ...]
    sequencing: str
    modality: str
    caveat: str = ""


DEFS = {
    "codec": Definition("CODEC — single-duplex error correction", "10.1038/s41588-023-01376-0", "Ding et al., Nature Genetics (2023)", "An adapter quadruplex and strand-displacing extension place readable copies of both strands of one original DNA duplex into a single Illumina cluster molecule.", "custom", ("Watson copy", "CODEC linkage", "Crick copy"), "Custom outward Illumina reads recover insert, UMI and index from the physically linked duplex.", "DNA", "Exact CODEC oligo sequences are supplementary; role-preserving regions are shown as inferred rather than as published bases."),
    "nanoseq": Definition("NanoSeq — nanorate duplex sequencing", "10.1038/s41586-021-03477-4", "Abascal et al., Nature (2021)", "Restriction fragmentation and dideoxynucleotide-based restrictive end repair prevent lesion copying before duplex consensus sequencing.", "truseq", ("duplex tag A", "genomic insert", "duplex tag B"), "Paired-end Illumina reads are grouped into complementary strand families and then a duplex consensus.", "DNA"),
    "stlfr": Definition("stLFR — single-tube long-fragment reads", "10.1101/gr.245126.118", "Wang et al., Genome Research (2019)", "Transposon-tagged subfragments from one long DNA molecule are captured on one bead and receive a shared linked-read barcode before DNBSEQ.", "dnb", ("linked-read barcode", "genomic subfragment"), "DNBSEQ reads include the bead barcode needed to group short reads back to a long source molecule.", "DNA"),
    "nano3p-seq": Definition("Nano3P-seq — Nanopore 3′-end capture", "10.1038/s41592-022-01714-w", "Begik et al., Nature Methods (2023)", "Template-switching reverse transcription captures native RNA 3′ ends for PCR-free Nanopore cDNA sequencing without RNA-adapter ligation.", "nanopore", ("motor adapter", "full-length cDNA from native 3-prime end"), "Nanopore traverses the cDNA from the captured native RNA 3′ end, retaining tail length and composition information.", "RNA"),
    "paris": Definition("PARIS — psoralen analysis of RNA interactions and structures", "10.1016/j.cell.2016.04.028", "Lu et al., Cell (2016)", "AMT photocrosslinking, two-dimensional gel enrichment and proximity ligation preserve the two partners of an RNA duplex in one chimera.", "custom", ("RNA A-derived cDNA", "RNA B-derived cDNA"), "Paired-end sequencing reads across or beside the chimeric RNA ligation junction.", "RNA"),
    "splash": Definition("SPLASH — sequencing psoralen-crosslinked RNA hybrids", "10.1016/j.molcel.2016.04.028", "Aw et al., Molecular Cell (2016)", "Biotinylated psoralen photocrosslinks RNA duplexes for streptavidin selection, proximity ligation and sequencing.", "custom", ("RNA A-derived cDNA", "RNA B-derived cDNA"), "Paired-end reads identify both partners retained in each proximity-ligated RNA chimera.", "RNA"),
    "ligr-seq": Definition("LIGR-seq — ligation of interacting RNA", "10.1016/j.molcel.2016.04.030", "Sharma et al., Molecular Cell (2016)", "AMT-crosslinked RNA duplexes survive RNase R enrichment and are proximity-ligated before cDNA library construction.", "custom", ("RNA A-derived cDNA", "RNA B-derived cDNA"), "Paired-end reads identify the two RNA regions joined at the proximity-ligation boundary.", "RNA"),
    "mario": Definition("MARIO — mapping RNA interactome in vivo", "10.1038/ncomms12023", "Nguyen et al., Nature Communications (2016)", "UV-crosslinked RNPs receive a selectable biotinylated RNA linker that joins two RNAs co-bound by one protein.", "custom", ("RNA A cDNA", "linker", "RNA B cDNA"), "Paired-end reads retain the linker boundary used to split and map the two interacting RNAs.", "RNA"),
    "maester": Definition("MAESTER — mitochondrial variant enrichment from scRNA-seq", "10.1038/s41587-022-01210-8", "Miller et al., Nature Biotechnology (2022)", "Multiplex PCR enriches mitochondrial transcripts from barcoded 3′ scRNA-seq cDNA while preserving cell barcodes and UMIs.", "truseq", ("cell barcode", "UMI", "mitochondrial target"), "Read 1 retains the parent scRNA-seq cell barcode and UMI; the opposite read reports the enriched mitochondrial target.", "RNA"),
    "redeem": Definition("ReDeeM — deep mitochondrial mutation profiling", "10.1038/s41586-024-07066-z", "Weng et al., Nature (2024)", "A whole-cell 10x Multiome modification retains mtDNA, then tiled hybrid capture produces a deeply sequenced mtDNA branch matched to RNA and ATAC.", "truseq", ("cell barcode", "mtDNA fragment"), "The mtDNA, ATAC and RNA libraries carry matched cell barcodes; mtDNA fragment endpoints support consensus error correction.", "multi"),
    "rage-seq": Definition("RAGE-Seq — repertoire and gene expression by sequencing", "10.1038/s41467-019-11049-4", "Singh et al., Nature Communications (2019)", "Hybrid capture enriches full-length immune-receptor cDNAs from barcoded droplet libraries for long-read sequencing alongside short-read expression.", "nanopore", ("motor adapter", "cell barcode", "full-length receptor cDNA"), "The long read retains the cell barcode and full receptor isoform; the matched short-read library supplies gene expression.", "RNA"),
    "malbac-dt": Definition("MALBAC-DT — high-detectability single-cell transcriptomics", "10.1101/2019.12.31.892190", "Chen et al., bioRxiv (2020; later PNAS 2022)", "Barcoded oligo(dT) reverse transcription is followed by MALBAC quasi-linear random amplification, pooling, tagmentation and indexed sequencing.", "nextera", ("cell barcode", "UMI", "transcript cDNA"), "Read 1 reports the cell barcode and UMI; Read 2 reports transcript sequence after TruePrep tagmentation.", "RNA"),
    "seq-scope": Definition("Seq-Scope — flow-cell spatial transcriptomics", "10.1016/j.cell.2021.05.010", "Cho et al., Cell (2021)", "An Illumina flow cell is first decoded as a dense HDMI barcode array and then reused to capture spatially resolved tissue transcripts.", "custom", ("spatial HDMI", "UMI", "transcript cDNA"), "The HDMI read maps a molecule to a decoded flow-cell coordinate; the transcript read identifies captured RNA.", "RNA"),
    "chil-seq": Definition("ChIL-seq — chromatin integration labeling", "10.1038/s41596-020-0375-8", "Handa et al., Nature Protocols (2020)", "An oligo-conjugated antibody targets chromatin; nearby Tn5 integration and T7 transcription linearly amplify target-proximal DNA before sequencing.", "custom", ("ChIL-probe end", "target-proximal DNA"), "Illumina sequencing starts from the ChIL-probe-derived end and reports nearby genomic sequence.", "DNA"),
    "bliss": Definition("BLISS — breaks labeling in situ and sequencing", "10.1038/ncomms15058", "Yan et al., Nature Communications (2017)", "Native DSBs are blunted and ligated in situ to a T7/RA5 adapter carrying a sample barcode and UMI, then linearly amplified by IVT.", "custom", ("UMI", "sample barcode", "break-adjacent DNA"), "The genomic read begins immediately beside the captured DSB; the inline UMI distinguishes molecular break events.", "DNA"),
    "scsplat": Definition("scSPLAT — scalable single-cell WGBS", "10.1101/2021.10.14.464375", "Raine et al., bioRxiv (2021; Scientific Reports 2022)", "Bisulfite-converted single-cell DNA receives a P5/cell-barcode random-primer copy before pooling and bulk splinted ligation of the opposite adapter.", "truseq", ("cell barcode", "bisulfite-converted insert"), "The first six bases of Read 1 are the inline cell barcode; paired genomic reads report the converted insert.", "DNA"),
    "cfmedip-seq": Definition("cfMeDIP-seq — cell-free DNA methylome enrichment", "10.1038/s41596-019-0202-2", "Shen et al., Nature Protocols (2019)", "Adapter-ligated low-input cfDNA is supplemented with controls and filler, denatured, and enriched with a 5mC antibody without bisulfite conversion.", "truseq", ("methylated cfDNA fragment",), "Ordinary paired-end Illumina reads report the endpoints and sequence of immunoprecipitated cfDNA fragments.", "DNA"),
    "dart-seq": Definition("DART-seq — antibody-free m6A mapping", "10.1038/s41592-019-0570-0", "Meyer, Nature Methods (2019)", "APOBEC1 fused to the YTH m6A-binding domain writes nearby C-to-U marks that become identifiable mismatches in an RNA-seq library.", "truseq", ("cDNA carrying m6A-adjacent C-to-T mark",), "Standard RNA-seq reads report C-to-T conversions adjacent to YTH-bound m6A sites.", "RNA"),
    "snubar-atac": Definition("SNuBar-ATAC — oligo sample multiplexing for scATAC", "10.1016/j.molcel.2021.09.026", "Wang et al., Molecular Cell (2021)", "A sample-barcode oligo anneals to the universal transposome tail during tagmentation before samples are pooled into 10x scATAC droplets.", "nextera", ("cell barcode", "sample barcode", "UMI", "accessible DNA"), "The standard 10x ATAC reads recover genomic insert and cell barcode; the SBO-derived read reports sample barcode and UMI.", "DNA"),
    "ric-seq": Definition("RIC-seq — RNA in situ conformation sequencing", "10.1038/s41586-020-2249-1", "Cai et al., Nature (2020)", "Crosslinked RNA–protein assemblies are fragmented and proximity-ligated in situ to capture global spatial RNA–RNA contacts.", "custom", ("RNA A cDNA", "RNA B cDNA"), "Paired-end reads identify the two RNA fragments joined at each in-situ ligation junction.", "RNA"),
    "plac-seq": Definition("PLAC-seq — proximity ligation-assisted ChIP-seq", "10.1038/cr.2016.137", "Fang et al., Cell Research (2016)", "In-nucleus restriction digestion and proximity ligation precede shearing and chromatin immunoprecipitation, enriching contacts anchored at a chosen protein or mark.", "truseq", ("locus A", "biotinylated contact junction", "locus B"), "Paired-end Illumina reads identify both genomic sides of a protein-anchored proximity-ligation product.", "DNA"),
}


def _payload(d: Definition) -> list[Segment]:
    out = []
    for name in d.payload:
        kwargs = {}
        tag = None
        low = name.lower()
        if "umi" in low:
            tag = "umi"; kwargs["feature"] = feature("umi", "umi", "random")
        elif "cell barcode" in low:
            tag = "cbc"; kwargs["feature"] = feature("cell_barcode", "cell_barcode", "unknown")
        elif "sample barcode" in low:
            tag = "cbc"; kwargs["feature"] = feature("sample_index", "sample_index", "unknown")
        elif "spatial hdmi" in low:
            tag = "cbc"; kwargs["feature"] = feature("spatial_hdmi", "spatial_barcode", "random")
        elif "linked-read barcode" in low:
            tag = "cbc"; kwargs["feature"] = feature("linked_read_barcode", "linked_read_barcode", "whitelist", whitelist="stLFR bead barcode set")
        out.append(Segment(name, "X" * (10 if tag else 28), tag, placeholder=True, **kwargs))
    return out


def get(key: str):
    d = DEFS[key]
    primers = ()
    if d.library == "truseq":
        final, primers = truseq_library(_payload(d), d.title + " final library")
    elif d.library == "nextera":
        final, primers = nextera_library(_payload(d), d.title + " final library")
    elif d.library == "nanopore":
        final = Construct([role(name) for name in d.payload], name=d.title + " final molecule")
    elif d.library == "dnb":
        final = Construct([role("DNBSEQ adapter"), *_payload(d), role("DNBSEQ adapter reverse complement")],
                          name=d.title + " final library")
    else:
        final = Construct([role("Illumina-compatible left adapter"), *_payload(d),
                           role("Illumina-compatible right adapter")],
                          name=d.title + " final library")
    return SimpleNamespace(
        TITLE=d.title, NOTES="01_" + key.replace("-", "_") + ".html",
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.source_label}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
        SEQ_PRIMERS=primers, FINAL_CAPTION=d.sequencing,
        SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.sequencing,
        SEQUENCING_UNAVAILABLE=("The defining source does not expose an exact sequencing-primer sequence; the molecular read entry is shown by role."
                                if not primers else ""),
        workflow=lambda: workflow_for(key), oligos=lambda: (),
    )
