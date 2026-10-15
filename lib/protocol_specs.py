"""Source-backed specifications for the named-protocol expansion batch.

This module is deliberately data-oriented.  Protocol directories remain separate public
units, while shared constructors make the recurring molecular operations consistent.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Callable

import seqprimers as sp
from rna_library_operations import finish_pcr, finish_rt_pcr, molecule
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import nextera_library, seg, truseq_library
from bliss_chem import adapter_a1 as bliss_adapter_a1, oligos as bliss_oligos
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature
from chil_chem import (AD1_NOMX as CHIL_AD1_NOMX, AD2_INDEX as CHIL_AD2_INDEX,
                       READ2_N11 as CHIL_READ2_N11, oligos as chil_oligos)
from circular import CircularizationGeometry
from duplex_error import (CODEC_ADAP5_1, CODEC_ADAP7_1, codec_linked_duplex,
                          codec_final_library, codec_oligos,
                          nanoseq_restrictive_end_repair)
from interaction_capture import (FixedRNPContact, ProteinMediatedRNAContact,
                                 ReversibleDuplexCrosslink)
from mario_chem import oligos as mario_oligos
from nano3p_chem import (COMPA_DNA as NANO3P_COMPA, D_DNA as NANO3P_D_DNA,
                        R_RNA as NANO3P_R_RNA, final_library as nano3p_final_library,
                        oligos as nano3p_oligos)
from pcr import PCRGeometry, PCRPrimer
from paris_chem import oligos as paris_oligos
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rage_seq_chem import oligos as rage_seq_oligos
from seqscope_chem import oligos as seqscope_oligos
from stlfr_chem import (barcode_transfer_geometry as stlfr_barcode_transfer_geometry,
                        oligos as stlfr_oligos)


def role(name: str, width: int = 18, tag: str | None = None, *, length_bp=None,
         inferred: bool = False, molecular_feature=None) -> Segment:
    token = f"[{length_bp}-bp {name}]" if length_bp else f"[{name}]"
    return Segment(name, token, tag, placeholder=True, length_bp=length_bp,
                   inferred=inferred, feature=molecular_feature)



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


def psoralen_workflow(reagent: str, selection: str) -> Workflow:
    """Construct the shared crosslink, selection and proximity-ligation workflow."""
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
    if key in {"paris", "ligr-seq"}:
        reagent, selection = {
        "paris": ("AMT", "cross-linked duplexes selected above the diagonal on a 2D gel"),
        "ligr-seq": ("AMT", "RNase R depletes uncross-linked linear RNA"),
        }[key]
        return psoralen_workflow(reagent, selection)
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
             "The four adapter oligos, modifications and matched index sets are printed in Supplementary Table 5."),
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
            ("Hybrid-capture mtDNA in four staggered probe pools",
             molecule("captured mtDNA branch", ("cell barcode", "mtDNA fragment", "capture probe")),
             "Published Supplementary Data 1 lists the probe sequences; fragment endpoints act as endogenous molecular identifiers."),
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
            ("End-repair and add one 3-prime dA to each end",
             molecule("dA-tailed cfDNA", ("3-prime dA", "cfDNA", "3-prime dA")),
             "The cohesive dA ends are the substrates for the following adapter ligation."),
            ("Ligate Illumina adapters to the dA-tailed cfDNA",
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
        contact = FixedRNPContact("formaldehyde crosslinking")
        return workflow(molecule("fixed cells", ("spatially proximal RNAs",)), (
            ("Cross-link RNA–protein assemblies in situ and fragment RNA",
             contact.complex(),
             "Formaldehyde preserves spatial contacts, including protein-mediated ones."),
            ("Label RNA ends and proximity-ligate them in situ",
             contact.chimera(),
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
    "codec": Definition("CODEC — single-duplex error correction", "10.1038/s41588-023-01376-0", "Bae et al., Nature Genetics (2023)", "An adapter quadruplex and strand-displacing extension place readable copies of both strands of one original DNA duplex into a single Illumina cluster molecule.", "custom", ("Watson copy", "CODEC linkage", "Crick copy"), "Standard TruSeq Read 1 and Read 2 primers bind within the two published bridge oligos and read outward through a 3-nt UMI, an 18-nt inline sample index and one original strand copy.", "DNA", "Supplementary Table 5 publishes all four adapter oligos in eight matched sample-index sets. Set 1 is shown as a representative; the linked post-extension product remains role-scaled where the source figure does not print one contiguous strand."),
    "nanoseq": Definition("NanoSeq — nanorate duplex sequencing", "10.1038/s41586-021-03477-4", "Abascal et al., Nature (2021)", "Restriction fragmentation and dideoxynucleotide-based restrictive end repair prevent lesion copying before duplex consensus sequencing.", "truseq", ("duplex tag A", "genomic insert", "duplex tag B"), "Paired-end Illumina reads are grouped into complementary strand families and then a duplex consensus.", "DNA"),
    "stlfr": Definition("stLFR — single-tube long-fragment reads", "10.1101/gr.245126.118", "Wang et al., Genome Research (2019)", "Transposon-tagged subfragments from one long DNA molecule are captured on one bead and receive a shared linked-read barcode before DNBSEQ.", "dnb", ("linked-read barcode", "genomic subfragment"), "The cited detailed protocol publishes the BGI Read 1, barcode and Read 2 sequencing primers and states paired-end 100-base genomic reads; barcode demultiplexing consumes a 42-base barcode/common-adapter string. The exact DNB circle-closing adapter or splint geometry is not published.", "DNA", "The paper-cited Protocol Exchange workbook publishes the common bead, bridge, transposon, PCR, branch-ligation, BGI and Illumina oligos, plus all 1,536 members of each of the three 10-nt barcode sets. The page shows the barcode-family architecture; the external workbook remains the authoritative full whitelist."),
    "nano3p-seq": Definition("Nano3P-seq — Nanopore 3′-end capture", "10.1038/s41592-022-01714-w", "Begik et al., Nature Methods (2023)", "Template-switching reverse transcription captures native RNA 3′ ends for PCR-free Nanopore cDNA sequencing without RNA-adapter ligation.", "nanopore", ("motor adapter", "full-length cDNA from native 3-prime end"), "Nanopore traverses the cDNA from the captured native RNA 3′ end, retaining tail length and composition information.", "RNA", "The defining paper prints all three custom Nano3P-seq oligos and their modifications. The SQK-DCS109 Adapter Mix motor-adapter sequence remains proprietary to Oxford Nanopore and is marked inferred."),
    "paris": Definition("PARIS — psoralen analysis of RNA interactions and structures", "10.1016/j.cell.2016.04.028", "Lu et al., Cell (2016)", "AMT photocrosslinking, two-dimensional gel enrichment and proximity ligation preserve the two partners of an RNA duplex in one chimera.", "custom", ("RNA A-derived cDNA", "RNA B-derived cDNA"), "The detailed PARIS protocol publishes the complete custom P6 sequencing primer and all adapter, RT and PCR oligos. Its circular-cDNA junction still needs to be represented before primer placement can be computed rather than stated by role.", "RNA", "The exact oligos were not absent: Lu et al.'s detailed PARIS protocol prints them, including modifications and the 24-member RT-index list. The current diagram retains a role-level final construct until the circular-cDNA library geometry is reconstructed."),
    "ligr-seq": Definition("LIGR-seq — ligation of interacting RNA", "10.1016/j.molcel.2016.04.030", "Sharma et al., Molecular Cell (2016)", "AMT-crosslinked RNA duplexes survive RNase R enrichment and are proximity-ligated before cDNA library construction.", "truseq", ("5-prime random barcode", "RNA A-derived cDNA", "RNA B-derived cDNA"), "The study analyzed single-end, strand-specific reads after trimming a 5-prime random barcode. The SMARTer kit supplies standard Illumina/TruSeq-compatible endpoints, so the standard sequencing-primer sites are derived centrally; the proprietary internal priming sequences remain unknown.", "RNA", "The defining paper and its complete supplement were checked. They identify the SMARTer Stranded RNA-Seq kit and reaction roles but do not print its N6, SMART, or PCR-primer sequences. The contemporary Clontech manual explicitly marks those component sequences proprietary, while identifying its Illumina-compatible endpoint and index families."),
    "mario": Definition("MARIO — mapping RNA interactome in vivo", "10.1038/ncomms12023", "Nguyen et al., Nature Communications (2016)", "UV-crosslinked RNPs receive a selectable biotinylated RNA linker that joins two RNAs co-bound by one protein.", "truseq", ("UMI", "RNA A cDNA", "MARIO linker", "RNA B cDNA"), "Standard paired-end Illumina reads traverse the two RNA-derived arms and linker. The paper publishes the full PE1/PE2 primers, three inline sample-barcode/UMI RT primers, and every custom linker/adapter oligo.", "RNA", "The original Methods section has a complete source-exact oligo list. Three RT-primer variants differ in their four-base sample barcode and random-barcode length; all three are shown."),
    "maester": Definition("MAESTER — mitochondrial variant enrichment from scRNA-seq", "10.1038/s41587-022-01210-8", "Miller et al., Nature Biotechnology (2022)", "Multiplex PCR enriches mitochondrial transcripts from barcoded 3′ scRNA-seq cDNA while preserving cell barcodes and UMIs.", "truseq", ("cell barcode", "UMI", "mitochondrial target"), "Read 1 retains the parent scRNA-seq cell barcode and UMI; the opposite read reports the enriched mitochondrial target.", "RNA"),
    "redeem": Definition("ReDeeM — deep mitochondrial mutation profiling", "10.1038/s41586-024-07066-z", "Weng et al., Nature (2024)", "A whole-cell 10x Multiome modification retains mtDNA, then tiled hybrid capture produces a deeply sequenced mtDNA branch matched to RNA and ATAC.", "truseq", ("cell barcode", "mtDNA fragment"), "The mtDNA, ATAC and RNA libraries carry matched cell barcodes; mtDNA fragment endpoints support consensus error correction.", "multi"),
    "rage-seq": Definition("RAGE-Seq — repertoire and gene expression by sequencing", "10.1038/s41467-019-11049-4", "Singh et al., Nature Communications (2019)", "Hybrid capture enriches full-length immune-receptor cDNAs from barcoded droplet libraries for long-read sequencing alongside short-read expression.", "nanopore", ("motor adapter", "cell barcode", "full-length receptor cDNA"), "The long read retains the cell barcode and full receptor isoform; the matched short-read library supplies gene expression. The paper prints the TSO and R1 primers used before and after both capture rounds; Oxford Nanopore kit adapters remain proprietary.", "RNA", "The original Methods and supplement were rechecked. Exact TSO and R1 amplification primers are published. The 644-exon Roche–NimbleGen capture panel is specified by target genes rather than probe sequences, and the SQK-LSK108/SQK-LSK308 adapter sequences are not disclosed."),
    "malbac-dt": Definition("MALBAC-DT — high-detectability single-cell transcriptomics", "10.1101/2019.12.31.892190", "Chen et al., bioRxiv (2020; later PNAS 2022)", "Barcoded oligo(dT) reverse transcription is followed by MALBAC quasi-linear random amplification, pooling, tagmentation and indexed sequencing.", "nextera", ("cell barcode", "UMI", "transcript cDNA"), "Read 1 reports the cell barcode and UMI; Read 2 reports transcript sequence after TruePrep tagmentation.", "RNA"),
    "seq-scope": Definition("Seq-Scope — flow-cell spatial transcriptomics", "10.1016/j.cell.2021.05.010", "Cho et al., Cell (2021)", "An Illumina flow cell is first decoded as a dense HDMI barcode array and then reused to capture spatially resolved tissue transcripts.", "truseq", ("spatial HDMI", "UMI", "transcript cDNA"), "The first run uses a published custom Read1-DraI primer to decode HDMI coordinates. The second paired-end run uses the completed TruSeq shell: one read reports HDMI/UMI and the other reports captured cDNA.", "RNA", "The open Methods print all primary and alternative seed oligos, both custom first-run primers, the random second-strand primer and both PCR rounds with phosphorothioates. These were previously missed."),
    "chil-seq": Definition("ChIL-seq — chromatin integration labeling", "10.1038/s41596-020-0375-8", "Handa et al., Nature Protocols (2020)", "An oligo-conjugated antibody targets chromatin; nearby Tn5 integration and T7 transcription linearly amplify target-proximal DNA before sequencing.", "custom", ("ChIL-probe end", "target-proximal DNA"), "Illumina sequencing starts from the ChIL-probe-derived end and reports nearby genomic sequence.", "DNA"),
    "bliss": Definition("BLISS — breaks labeling in situ and sequencing", "10.1038/ncomms15058", "Yan et al., Nature Communications (2017)", "Native DSBs are blunted and ligated in situ to a T7/RA5 adapter carrying a sample barcode and UMI, then linearly amplified by IVT.", "custom", ("UMI", "sample barcode", "break-adjacent DNA"), "The genomic read begins immediately beside the captured DSB; the inline UMI distinguishes molecular break events.", "DNA"),
    "scsplat": Definition("scSPLAT — scalable single-cell WGBS", "10.1101/2021.10.14.464375", "Raine et al., bioRxiv (2021; Scientific Reports 2022)", "Bisulfite-converted single-cell DNA receives a P5/cell-barcode random-primer copy before pooling and bulk splinted ligation of the opposite adapter.", "truseq", ("cell barcode", "bisulfite-converted insert"), "The first six bases of Read 1 are the inline cell barcode; paired genomic reads report the converted insert.", "DNA"),
    "cfmedip-seq": Definition("cfMeDIP-seq — cell-free DNA methylome enrichment", "10.1038/s41596-019-0202-2", "Shen et al., Nature Protocols (2019)", "Adapter-ligated low-input cfDNA is supplemented with controls and filler, denatured, and enriched with a 5mC antibody without bisulfite conversion.", "truseq", ("methylated cfDNA fragment",), "Ordinary paired-end Illumina reads report the endpoints and sequence of immunoprecipitated cfDNA fragments.", "DNA"),
    "dart-seq": Definition("DART-seq — antibody-free m6A mapping", "10.1038/s41592-019-0570-0", "Meyer, Nature Methods (2019)", "APOBEC1 fused to the YTH m6A-binding domain writes nearby C-to-U marks that become identifiable mismatches in an RNA-seq library.", "truseq", ("cDNA carrying m6A-adjacent C-to-T mark",), "Standard RNA-seq reads report C-to-T conversions adjacent to YTH-bound m6A sites.", "RNA"),
    "snubar-atac": Definition("SNuBar-ATAC — oligo sample multiplexing for scATAC", "10.1016/j.molcel.2021.09.026", "Wang et al., Molecular Cell (2021)", "A sample-barcode oligo anneals to the universal transposome tail during tagmentation before samples are pooled into 10x scATAC droplets.", "nextera", ("cell barcode", "sample barcode", "UMI", "accessible DNA"), "The standard 10x ATAC reads recover genomic insert and cell barcode; the SBO-derived read reports sample barcode and UMI.", "DNA"),
    "ric-seq": Definition("RIC-seq — RNA in situ conformation sequencing", "10.1038/s41586-020-2249-1", "Cai et al., Nature (2020)", "Crosslinked RNA–protein assemblies are fragmented and proximity-ligated in situ to capture global spatial RNA–RNA contacts.", "custom", ("RNA A cDNA", "RNA B cDNA"), "Paired-end reads identify the two RNA fragments joined at each in-situ ligation junction.", "RNA"),
    "plac-seq": Definition("PLAC-seq — proximity ligation-assisted ChIP-seq", "10.1038/cr.2016.137", "Fang et al., Cell Research (2016)", "In-nucleus restriction digestion and proximity ligation precede shearing and chromatin immunoprecipitation, enriching contacts anchored at a chosen protein or mark.", "truseq", ("locus A", "biotinylated contact junction", "locus B"), "Paired-end Illumina reads identify both genomic sides of a protein-anchored proximity-ligation product.", "DNA"),
}



def psoralen_graph(d, final, primers, *, reagent: str,
                   selection_label: str):
    """Build a psoralen contact graph from source-specific reagent choices."""
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    topology = ReversibleDuplexCrosslink(reagent, selection_label)
    rna = graph.input("rna", "cellular RNA containing paired regions",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("cellular RNA", ("RNA molecules",)), primary=True)
    crosslinker = graph.input("crosslinker", reagent + " and long-wave UV condition",
                              MaterialKind.CONDITION)
    crosslinked = graph.product("crosslinked", reagent + "-crosslinked RNA duplexes",
                                MaterialKind.MOLECULAR,
                                state=topology.duplex())
    graph.react("crosslink", "Photocrosslink paired RNA with " + reagent,
                StepKind.CROSS_LINKING,
                ((rna, InputRole.TEMPLATE), (crosslinker, InputRole.CONDITION)),
                (crosslinked,))
    fragmentation = graph.input("fragmentation_condition", "RNA fragmentation condition",
                                MaterialKind.CONDITION)
    fragments = graph.product("fragments", "fragmented crosslinked RNA duplexes",
                              MaterialKind.MOLECULAR, state=crosslinked.state)
    graph.react("fragment", "Fragment the crosslinked RNA", StepKind.FRAGMENTATION,
                ((crosslinked, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragments,))
    selector = graph.input("selector", selection_label, MaterialKind.SELECTION)
    selected = graph.product("selected", "enriched crosslinked RNA duplexes",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("rejected", "nonselected RNA fragments", MaterialKind.WASTE)
    graph.react("select_duplexes", "Enrich crosslinked RNA duplexes",
                StepKind.SELECTION,
                ((fragments, InputRole.TEMPLATE), (selector, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "does not retain the selected crosslinked duplex population")
    ligase = graph.input("proximity_ligase", "RNA proximity ligase", MaterialKind.CATALYST)
    chimera = graph.product("chimera", "proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=topology.chimera())
    graph.react("proximity_ligation", "Proximity-ligate the paired RNA fragments",
                StepKind.PROXIMITY_LIGATION,
                ((selected, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera,))
    reversal = graph.input("photoreversal", reagent + " photoreversal condition",
                           MaterialKind.CONDITION)
    reversed_chimera = graph.product("reversed_chimera", "photoreversed RNA chimera",
                                     MaterialKind.MOLECULAR, state=chimera.state)
    graph.react("photoreverse", "Photoreverse the psoralen crosslink",
                StepKind.CHEMICAL_CONVERSION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (reversed_chimera,))
    reviewed = ["rna", "crosslinker", "crosslink", "fragmentation_condition", "fragment",
                "selector", "select_duplexes", "proximity_ligase", "proximity_ligation",
                "photoreversal", "photoreverse"]
    finish_rt_pcr(graph, "library", reversed_chimera, final, primers, d.sequencing,
                   reviewed, label=d.title)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _mario_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rnps = graph.input("rnps", "cellular protein-bound RNA assemblies",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("cellular RNPs", ("protein-bound RNAs",)), primary=True)
    uv = graph.input("uv", "UV crosslinking condition", MaterialKind.CONDITION)
    fixed = graph.product("fixed", "UV-crosslinked RNP assemblies", MaterialKind.MOLECULAR,
                          state=rnps.state)
    graph.react("crosslink", "UV-crosslink RNA to bound proteins", StepKind.CROSS_LINKING,
                ((rnps, InputRole.TEMPLATE), (uv, InputRole.CONDITION)), (fixed,))
    linker = graph.input("linker", "biotinylated RNA linker", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "RNA proximity ligase", MaterialKind.CATALYST)
    chimera = graph.product("chimera", "RNA1-linker-RNA2 chimera", MaterialKind.MOLECULAR,
                            state=molecule("MARIO chimera", ("RNA A", "biotinylated linker", "RNA B"),
                                           junctions=(("RNA A", "biotinylated linker", "linker ligation"),
                                                      ("biotinylated linker", "RNA B", "linker ligation"))))
    graph.react("ligate_linker", "Attach the linker and ligate co-bound RNAs",
                StepKind.LIGATION,
                ((fixed, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (chimera,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The public schematic establishes the selectable linker boundary but not one universal end geometry."))
    beads = graph.input("streptavidin", "streptavidin capture beads", MaterialKind.SELECTION)
    captured = graph.product("captured", "linker-containing RNA chimeras", MaterialKind.MOLECULAR,
                             state=chimera.state)
    rejected = graph.product("rejected", "linker-negative RNA", MaterialKind.WASTE)
    graph.react("capture", "Select linker-containing chimeras", StepKind.AFFINITY_CAPTURE,
                ((chimera, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not contain the biotinylated interaction linker")
    reviewed = ["rnps", "uv", "crosslink", "linker", "ligase", "ligate_linker",
                "streptavidin", "capture"]
    finish_rt_pcr(graph, "library", captured, final, primers, d.sequencing, reviewed,
                   label="MARIO chimera")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _codec_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    duplex = graph.input("duplex", "original genomic Watson/Crick duplex",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("original genomic duplex", ("Watson strand", "Crick strand")),
                         primary=True)
    repair = graph.input("end_repair_enzyme", "end-repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired genomic duplex", MaterialKind.MOLECULAR,
                             state=duplex.state)
    graph.react("end_repair", "End-repair the original duplex", StepKind.END_REPAIR,
                ((duplex, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    adapter = graph.input("codec_adapter", "preassembled four-oligo CODEC adapter quadruplex",
                          MaterialKind.ADAPTER,
                          source=source + "; Supplementary Table 5")
    ligase = graph.input("ligase", "CODEC adapter ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "CODEC-adapter-bound duplex", MaterialKind.MOLECULAR,
                            state=molecule("adapter-bound duplex", (
                                "left CODEC adapter", "original duplex", "right CODEC adapter"),
                                junctions=(("left CODEC adapter", "original duplex", "adapter ligation"),
                                           ("original duplex", "right CODEC adapter", "adapter ligation"))))
    graph.react("ligate", "Ligate the preassembled CODEC adapter quadruplex",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive",
                    *reported_adapter_end(
                        "LD4-adap5-1 arm", kind="cohesive", polarity="3-prime",
                        terminal_sequence="T", body_sequence=CODEC_ADAP5_1[:-1],
                        statement="Supplementary Table 5 prints the phosphorothioate-linked terminal T."
                    ).scenes,
                    *reported_adapter_end(
                        "LD4-adap7-1 arm", kind="cohesive", polarity="3-prime",
                        terminal_sequence="T", body_sequence=CODEC_ADAP7_1[:-1],
                        statement="Supplementary Table 5 prints the phosphorothioate-linked terminal T."
                    ).scenes,
                    statement=("The paper dA-tails the insert and Supplementary Table 5 prints a "
                               "terminal phosphorothioate-linked T on each ligating adapter arm.")))
    strand_polymerase = graph.input("strand_polymerase", "strand-displacing polymerase",
                                    MaterialKind.CATALYST)
    linked = graph.product("linked", "physically linked copies of both original complements",
                           MaterialKind.MOLECULAR, state=codec_linked_duplex())
    graph.react("link_duplex", "Copy and physically link both complements by strand displacement",
                StepKind.EXTENSION,
                ((adapted, InputRole.TEMPLATE), (strand_polymerase, InputRole.CATALYST)),
                (linked,))
    reviewed = ["duplex", "end_repair_enzyme", "end_repair", "codec_adapter", "ligase",
                "ligate", "strand_polymerase", "link_duplex"]
    finish_pcr(graph, "library", linked, final, primers, d.sequencing, reviewed,
                label="CODEC linked duplex")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _nanoseq_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "high-molecular-weight genomic DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("high-molecular-weight genomic DNA", ("genomic DNA",)),
                      primary=True)
    restriction = graph.input("restriction_enzyme", "restriction endonuclease",
                              MaterialKind.CATALYST)
    fragments = graph.product("fragments", "restriction-defined genomic fragments",
                              MaterialKind.MOLECULAR, state=nanoseq_restrictive_end_repair())
    graph.react("digest", "Restriction-fragment genomic DNA", StepKind.DIGESTION,
                ((dna, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (fragments,))
    repair = graph.input("repair_mix", "dideoxynucleotide restrictive end-repair mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "restrictively end-repaired fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    graph.react("restrictive_end_repair", "Perform restrictive end repair",
                StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    adapters = graph.input("duplex_adapters", "complementary duplex-family adapters",
                           MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "duplex-adapter DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "duplex-tagged genomic fragments", MaterialKind.MOLECULAR,
                            state=molecule("duplex-consensus precursor", (
                                "duplex tag A", "genomic insert", "duplex tag B"),
                                junctions=(("duplex tag A", "genomic insert", "adapter ligation"),
                                           ("genomic insert", "duplex tag B", "adapter ligation"))))
    graph.react("ligate", "Ligate duplex-family adapters", StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source establishes complementary duplex tags; the role-level page does not assert one adapter-end sequence."))
    reviewed = ["dna", "restriction_enzyme", "digest", "repair_mix",
                "restrictive_end_repair", "duplex_adapters", "ligase", "ligate"]
    finish_pcr(graph, "library", adapted, final, primers, d.sequencing, reviewed,
                label="NanoSeq duplex")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _stlfr_graph(d, final):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "20–300-kb genomic DNA molecules", MaterialKind.BIOLOGICAL,
                      state=molecule("long genomic molecule", ("long genomic molecule",)),
                      primary=True)
    transposome = graph.input("transposome", "capture-sequence-loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "long molecules carrying addressable transposons",
                           MaterialKind.MOLECULAR,
                           state=molecule("tagged long DNA", (
                               "subfragment", "capture transposon", "subfragment")))
    waste = graph.product("tn5_waste", "non-library transposition products", MaterialKind.WASTE)
    graph.react("tagment", "Insert capture transposons along each long molecule",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an addressable genomic subfragment")
    bead = graph.input("barcoded_bead", "single barcoded capture bead", MaterialKind.SELECTION)
    captured = graph.product("captured", "one long molecule captured per barcoded bead",
                             MaterialKind.MOLECULAR,
                             state=molecule("bead-bound molecule", (
                                 "bead barcode", "captured subfragments")))
    rejected = graph.product("uncaptured", "uncaptured or multiply assigned fragments",
                             MaterialKind.WASTE)
    graph.react("capture", "Capture one long molecule on one barcoded bead",
                StepKind.AFFINITY_CAPTURE,
                ((tagged, InputRole.TEMPLATE), (bead, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not have an unambiguous long-molecule bead assignment")
    barcode_adapter = graph.input("barcode_adapter", "bead-barcode ligating oligo",
                                  MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "barcode DNA ligase", MaterialKind.CATALYST)
    barcoded = graph.product("barcoded", "co-barcoded genomic subfragments",
                             MaterialKind.MOLECULAR,
                             state=molecule("co-barcoded library precursor", (
                                 "linked-read barcode", "short genomic subfragment"),
                                 junctions=(("linked-read barcode", "short genomic subfragment",
                                             "barcode ligation"),)))
    graph.react("ligate_barcode", "Ligate the bead barcode to captured subfragments",
                StepKind.LIGATION,
                ((captured, InputRole.TEMPLATE), (barcode_adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (barcoded,),
                geometry=stlfr_barcode_transfer_geometry())
    circle_ligase = graph.input("circle_ligase", "DNB library circularization ligase",
                                MaterialKind.CATALYST)
    library = graph.product("library", "stLFR DNBSEQ library", MaterialKind.LIBRARY,
                            state=final)
    graph.react("circularize", "Circularize the co-barcoded products for DNB formation",
                StepKind.CIRCULARIZATION,
                ((barcoded, InputRole.TEMPLATE), (circle_ligase, InputRole.CATALYST)),
                (library,), geometry=CircularizationGeometry.unavailable(
                    "The stLFR source identifies DNB library circularization but does not disclose the circle-closing adapter/splint geometry."))
    graph.sequence(library, SequencingPlan(limitation=d.sequencing))
    reviewed = ["dna", "transposome", "tn5", "tagment", "barcoded_bead", "capture",
                "barcode_adapter", "ligase", "ligate_barcode", "circle_ligase",
                "circularize", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _nano3p_graph(d, final):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input("rna", "native RNA with an unmodified 3-prime end",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("native RNA", ("RNA body", "native 3-prime end")),
                      primary=True)
    starter = graph.input("rt_starter", "preannealed blocked R_RNA / 5-prime-phosphorylated D_DNA starter duplex",
                          MaterialKind.OLIGO,
                          source=source + "; Methods print both oligos and modifications")
    rt = graph.input("reverse_transcriptase", "template-switching reverse transcriptase",
                     MaterialKind.CATALYST)
    full = graph.product("full_cdna", "TGIRT template-switched cDNA copied from the native RNA 3-prime end",
                         MaterialKind.MOLECULAR,
                         state=molecule("full cDNA", ("D_DNA handle", "RNA-derived cDNA")))
    graph.react("template_switch_rt", "Template-switch TGIRT from the starter duplex onto the native RNA 3-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((rna, InputRole.TEMPLATE), (starter, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (full,))
    compa = graph.input("compa", "published CompA_DNA ligation complement",
                        MaterialKind.OLIGO, source=source + "; Methods")
    ready = graph.product("ligation_ready", "CompA-annealed cDNA end",
                          MaterialKind.MOLECULAR, state=full.state)
    graph.react("anneal_compa", "Anneal CompA_DNA to make the adapter-ligation end",
                StepKind.HYBRIDIZATION,
                ((full, InputRole.TEMPLATE), (compa, InputRole.STRUCTURAL)), (ready,))
    motor = graph.input("motor_adapter", "Nanopore motor sequencing adapter",
                        MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "Nanopore adapter ligase", MaterialKind.CATALYST)
    library = graph.product("library", "PCR-free Nano3P-seq molecule", MaterialKind.LIBRARY,
                            state=final)
    graph.react("ligate_motor", "Attach the Nanopore motor adapter without PCR",
                StepKind.LIGATION,
                ((ready, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source defines the motor-adapter role; kit-revision-specific end chemistry is not asserted."))
    graph.sequence(library, SequencingPlan(limitation=d.sequencing))
    reviewed = ["rna", "rt_starter", "reverse_transcriptase", "template_switch_rt",
                "compa", "anneal_compa", "motor_adapter", "adapter_ligase",
                "ligate_motor", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _maester_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cdna = graph.input("cdna", "barcoded 3-prime single-cell RNA-seq cDNA",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("barcoded scRNA-seq cDNA", ("cell barcode", "UMI", "cDNA")),
                       primary=True)
    target_forward = graph.input("target_forward", "multiplex mitochondrial target primers",
                                 MaterialKind.OLIGO)
    target_reverse = graph.input("target_reverse", "parent-library reverse primer",
                                 MaterialKind.OLIGO)
    polymerase = graph.input("target_polymerase", "target PCR polymerase",
                             MaterialKind.CATALYST)
    enriched = graph.product("enriched", "cell-barcode/UMI-preserving mitochondrial amplicons",
                             MaterialKind.MOLECULAR,
                             state=molecule("targeted amplicon", (
                                 "cell barcode", "UMI", "mitochondrial target")))
    graph.react("target_pcr", "Multiplex-PCR mitochondrial transcript targets",
                StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (target_forward, InputRole.STRUCTURAL),
                 (target_reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (enriched,), geometry=PCRGeometry.unavailable(
                    enriched.state, reason="The multiplex target panel, rather than one universal amplicon, defines the PCR geometry."))
    reviewed = ["cdna", "target_forward", "target_reverse", "target_polymerase", "target_pcr"]
    finish_pcr(graph, "library", enriched, final, primers, d.sequencing, reviewed,
                label="MAESTER mitochondrial-target")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _malbac_dt_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input("rna", "single-cell poly(A) RNA", MaterialKind.BIOLOGICAL,
                      state=molecule("single-cell poly(A) RNA", ("transcripts",)), primary=True)
    rt_primer = graph.input("rt_primer", "cell-barcode/UMI oligo(dT) RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "cell-barcoded first-strand cDNA", MaterialKind.MOLECULAR,
                         state=molecule("barcoded first-strand cDNA", (
                             "cell barcode", "UMI", "transcript cDNA")))
    graph.react("reverse_transcribe", "Reverse-transcribe with barcode, UMI and oligo(dT)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    random_primer = graph.input("malbac_primer", "MALBAC common-handle random primer",
                                MaterialKind.OLIGO)
    strand_polymerase = graph.input("malbac_polymerase", "strand-displacing MALBAC polymerase",
                                    MaterialKind.CATALYST)
    amplified = graph.product("amplified", "quasi-linearly amplified MALBAC-DT cDNA",
                              MaterialKind.MOLECULAR,
                              state=molecule("MALBAC-DT amplicons", (
                                  "common handle", "barcoded cDNA", "common handle")))
    graph.react("malbac", "Perform MALBAC quasi-linear random-primer amplification",
                StepKind.EXTENSION,
                ((cdna, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (strand_polymerase, InputRole.CATALYST)), (amplified,))
    transposome = graph.input("transposome", "TruePrep loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "TruePrep transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented MALBAC-DT products", MaterialKind.MOLECULAR,
                           state=amplified.state)
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Pool and tagment the amplified cDNA", StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an amplifiable barcoded transcript insert")
    reviewed = ["rna", "rt_primer", "reverse_transcriptase", "reverse_transcribe",
                "malbac_primer", "malbac_polymerase", "malbac", "transposome", "tn5", "tagment"]
    finish_pcr(graph, "library", tagged, final, primers, d.sequencing, reviewed,
                label="MALBAC-DT")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _seq_scope_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    surface = graph.input("surface", "Illumina flow cell bearing random HDMI capture oligos",
                          MaterialKind.BIOLOGICAL,
                          state=molecule("flow-cell array", ("random HDMI oligos",)), primary=True)
    decoder = graph.input("decoder", "in-situ HDMI sequencing primers",
                          MaterialKind.OLIGO)
    decoded = graph.product("decoded", "coordinate-decoded HDMI capture surface",
                            MaterialKind.MOLECULAR,
                            state=molecule("decoded capture surface", (
                                "mapped HDMI cluster", "poly(dT) capture end")))
    graph.react("decode", "Cluster and sequence the HDMI array in situ",
                StepKind.HYBRIDIZATION,
                ((surface, InputRole.TEMPLATE), (decoder, InputRole.STRUCTURAL)), (decoded,))
    tissue_rna = graph.input("tissue_rna", "poly(A) RNA released from overlaid tissue",
                             MaterialKind.MOLECULAR)
    captured = graph.product("captured", "RNA hybridized to decoded surface clusters",
                             MaterialKind.MOLECULAR, state=decoded.state)
    graph.react("capture", "Capture tissue mRNA on the decoded flow cell",
                StepKind.HYBRIDIZATION,
                ((decoded, InputRole.TEMPLATE), (tissue_rna, InputRole.STRUCTURAL)),
                (captured,))
    rt_primer = graph.input("rt_primer", "surface-bound poly(dT) RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "surface cDNA carrying HDMI and UMI",
                         MaterialKind.MOLECULAR,
                         state=molecule("surface cDNA", ("spatial HDMI", "UMI", "transcript cDNA")))
    graph.react("reverse_transcribe", "Reverse-transcribe captured RNA on the surface",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    reviewed = ["surface", "decoder", "decode", "tissue_rna", "capture", "rt_primer",
                "reverse_transcriptase", "reverse_transcribe"]
    finish_pcr(graph, "library", cdna, final, primers, d.sequencing, reviewed,
                label="Seq-Scope spatial transcript")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _bliss_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "fixed DNA containing native double-strand breaks",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("fixed DNA", ("native double-strand break",)), primary=True)
    repair = graph.input("blunting_enzyme", "break-end blunting enzyme mix",
                         MaterialKind.CATALYST)
    blunted = graph.product("blunted", "blunted native DNA break ends", MaterialKind.MOLECULAR,
                            state=molecule("blunted break", ("genomic DNA", "blunt DSB end")))
    graph.react("blunt", "Blunt the native break ends in situ", StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (blunted,))
    adapter = graph.input("break_adapter", "published UMI/sample-barcode/RA5/T7 break adapter",
                          MaterialKind.ADAPTER, state=bliss_adapter_a1(),
                          source=source + "; Supplementary Data 1 and Supplementary Fig. 1a")
    ligase = graph.input("ligase", "break-adapter DNA ligase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "adapter-tagged break ends", MaterialKind.MOLECULAR,
                           state=molecule("tagged break", (
                               "genomic break", "UMI", "sample barcode", "RA5", "T7 promoter"),
                               junctions=(("genomic break", "UMI", "break-adapter ligation"),)))
    graph.react("ligate", "Ligate the barcoded T7 adapter to break ends", StepKind.LIGATION,
                ((blunted, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (tagged,),
                geometry=AdapterLigationGeometry.known(
                    "blunt", bliss_adapter_a1(),
                    statement=("Supplementary Fig. 1a marks ligation at the blunt, "
                               "5-prime-phosphorylated sample-barcode end; Supplementary "
                               "Data 1 prints both oligos.")))
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    rna = graph.product("rna", "linearly amplified break-adjacent RNA", MaterialKind.MOLECULAR,
                        state=tagged.state)
    graph.react("transcribe", "Linearly amplify tagged breaks by T7 transcription",
                StepKind.TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (rna,))
    reviewed = ["dna", "blunting_enzyme", "blunt", "break_adapter", "ligase", "ligate",
                "t7_polymerase", "transcribe"]
    finish_rt_pcr(graph, "library", rna, final, primers, d.sequencing, reviewed,
                   label="BLISS break")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _scsplat_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "single-cell genomic DNA", MaterialKind.BIOLOGICAL,
                      state=molecule("one cell", ("genomic DNA",)), primary=True)
    bisulfite = graph.input("bisulfite", "bisulfite conversion condition",
                            MaterialKind.CONDITION)
    converted = graph.product("converted", "single-stranded bisulfite-converted DNA",
                              MaterialKind.MOLECULAR,
                              state=molecule("converted DNA", ("bisulfite-converted fragments",)))
    graph.react("convert", "Bisulfite-convert DNA in each well",
                StepKind.CHEMICAL_CONVERSION,
                ((dna, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)), (converted,))
    random_primer = graph.input("random_primer", "P5/cell-barcode random primer",
                                MaterialKind.OLIGO)
    polymerase = graph.input("strand_polymerase", "second-strand polymerase",
                             MaterialKind.CATALYST)
    barcoded = graph.product("barcoded", "P5/cell-barcode second-strand copy",
                             MaterialKind.MOLECULAR,
                             state=molecule("barcoded strand", (
                                 "P5", "cell barcode", "random-primer copy")))
    graph.react("second_strand", "Random-prime second-strand synthesis with P5 and cell barcode",
                StepKind.EXTENSION,
                ((converted, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (barcoded,))
    adapter = graph.input("p7_adapter", "splinted P7 adapter", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "splinted-adapter DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "pooled splint-ligated single-cell templates",
                            MaterialKind.MOLECULAR, state=barcoded.state)
    graph.react("ligate", "Pool wells and install P7 by splinted ligation",
                StepKind.LIGATION,
                ((barcoded, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source identifies splinted P7 installation; this role-level page does not assert the panel-specific splint bases."))
    reviewed = ["dna", "bisulfite", "convert", "random_primer", "strand_polymerase",
                "second_strand", "p7_adapter", "ligase", "ligate"]
    finish_pcr(graph, "library", adapted, final, primers, d.sequencing, reviewed,
                label="scSPLAT")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _dart_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input("rna", "cellular RNA containing m6A sites", MaterialKind.BIOLOGICAL,
                      state=molecule("cellular RNA", ("RNA containing m6A",)), primary=True)
    fusion = graph.input("fusion", "APOBEC1–YTH m6A-reader fusion", MaterialKind.SELECTION)
    targeted = graph.product("targeted", "RNA bound by APOBEC1–YTH near m6A",
                             MaterialKind.MOLECULAR, state=rna.state)
    unbound = graph.product("unbound", "RNA not selected by the reader fusion", MaterialKind.WASTE)
    graph.react("bind", "Recruit APOBEC1 to m6A-proximal RNA",
                StepKind.BINDING,
                ((rna, InputRole.TEMPLATE), (fusion, InputRole.SELECTION)),
                (targeted, unbound))
    graph.discard(unbound, "does not receive the targeted deamination mark")
    deamination = graph.input("deamination_condition", "APOBEC cytidine-deamination condition",
                              MaterialKind.CONDITION)
    marked = graph.product("marked", "RNA carrying m6A-adjacent C-to-U marks",
                           MaterialKind.MOLECULAR,
                           state=molecule("marked RNA", ("RNA", "C-to-U mark near m6A")))
    graph.react("deaminate", "Convert nearby cytidine to uridine",
                StepKind.CHEMICAL_CONVERSION,
                ((targeted, InputRole.TEMPLATE), (deamination, InputRole.CONDITION)),
                (marked,))
    reviewed = ["rna", "fusion", "bind", "deamination_condition", "deaminate"]
    finish_rt_pcr(graph, "library", marked, final, primers, d.sequencing, reviewed,
                   label="DART-seq marked RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _snubar_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "sample-separated nuclei containing accessible chromatin",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("sample-specific nuclei", ("accessible chromatin",)),
                         primary=True)
    transposome = graph.input("transposome", "standard 10x ATAC loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented chromatin carrying the universal tail",
                           MaterialKind.MOLECULAR,
                           state=molecule("tagmented chromatin", (
                               "genomic fragment", "transposome universal tail")))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment each sample with the 10x ATAC transposome",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an accessible-DNA insert with the universal tail")
    snubar = graph.input("snubar", "sample-barcode/UMI SNuBar oligo",
                         MaterialKind.OLIGO)
    sample_tagged = graph.product("sample_tagged", "SNuBar-tagged ATAC fragments",
                                  MaterialKind.MOLECULAR,
                                  state=molecule("sample-tagged fragment", (
                                      "genomic fragment", "sample barcode", "UMI", "PCR handle")))
    graph.react("anneal_snubar", "Anneal the SNuBar oligo to the transposome tail",
                StepKind.HYBRIDIZATION,
                ((tagged, InputRole.TEMPLATE), (snubar, InputRole.STRUCTURAL)),
                (sample_tagged,))
    droplets = graph.input("droplets", "10x single-cell ATAC droplet condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "cell- and sample-barcoded ATAC molecules",
                             MaterialKind.MOLECULAR,
                             state=molecule("SNuBar-ATAC products", (
                                 "10x cell barcode", "sample barcode", "UMI", "ATAC fragment")))
    collisions = graph.product("collisions", "empty or ambiguous droplets", MaterialKind.WASTE)
    graph.react("partition", "Pool samples and partition nuclei into 10x droplets",
                StepKind.PARTITION,
                ((sample_tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, collisions), changes_molecule=False)
    graph.discard(collisions, "does not retain an interpretable cell barcode")
    reviewed = ["nuclei", "transposome", "tn5", "tagment", "snubar", "anneal_snubar",
                "droplets", "partition"]
    finish_pcr(graph, "library", barcoded, final, primers, d.sequencing, reviewed,
                label="SNuBar-ATAC")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _redeem_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "whole cells retaining nuclear chromatin, RNA and mtDNA",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("whole cell", ("nuclear chromatin", "RNA", "mtDNA")),
                        primary=True)
    droplets = graph.input("multiome_droplets", "10x Multiome droplet/barcoding condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "matched cell-barcoded RNA, ATAC and mtDNA products",
                             MaterialKind.MOLECULAR,
                             state=molecule("matched cell products", (
                                 "cell barcode", "RNA products", "ATAC fragments", "mtDNA fragments")))
    failed = graph.product("failed_partitions", "empty or ambiguous droplets", MaterialKind.WASTE)
    graph.react("multiome", "Run 10x Multiome while retaining mitochondrial DNA",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, failed), changes_molecule=False)
    graph.discard(failed, "does not retain an interpretable shared cell barcode")
    mt = graph.product("mt_branch", "cell-barcoded mitochondrial-DNA branch",
                       MaterialKind.MOLECULAR, state=barcoded.state)
    atac = graph.product("atac_branch", "matched ATAC branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    rna = graph.product("rna_branch", "matched RNA branch", MaterialKind.MOLECULAR,
                        state=barcoded.state)
    graph.react("split_branches", "Separate mtDNA, ATAC and RNA branches",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (mt, atac, rna),
                changes_molecule=False)
    probes = graph.input("capture_probes", "four staggered mtDNA capture-probe pools (sequences in Supplementary Data 1)",
                         MaterialKind.SELECTION)
    captured = graph.product("captured_mt", "hybrid-captured cell-barcoded mtDNA",
                             MaterialKind.MOLECULAR, state=mt.state)
    rejected = graph.product("rejected_mt", "non-mitochondrial capture material",
                             MaterialKind.WASTE)
    graph.react("capture_mt", "Hybrid-capture mtDNA in four staggered probe pools",
                StepKind.AFFINITY_CAPTURE,
                ((mt, InputRole.TEMPLATE), (probes, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not belong to the selected mitochondrial branch")
    reviewed = ["cells", "multiome_droplets", "multiome", "split_branches",
                "capture_probes", "capture_mt"]
    finish_pcr(graph, "mt", captured, final, primers, "mtDNA branch: " + d.sequencing,
                reviewed, label="ReDeeM mtDNA")
    atac_final, atac_primers = nextera_library(
        [role("cell barcode"), role("accessible-DNA insert")], "ReDeeM ATAC library")
    rna_final, rna_primers = truseq_library(
        [role("cell barcode"), role("UMI"), role("transcript cDNA")], "ReDeeM RNA library")
    finish_pcr(graph, "atac", atac, atac_final, atac_primers,
                "matched ATAC branch: " + d.sequencing, reviewed, label="ReDeeM ATAC")
    finish_pcr(graph, "rna", rna, rna_final, rna_primers,
                "matched RNA branch: " + d.sequencing, reviewed, label="ReDeeM RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _rage_graph(d, final):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cdna = graph.input("cdna", "droplet scRNA-seq cDNA carrying cell barcode and UMI",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("droplet scRNA-seq cDNA", (
                           "cell barcode", "UMI", "full-length cDNA")), primary=True)
    expression = graph.product("expression_branch", "short-read expression branch",
                               MaterialKind.MOLECULAR, state=cdna.state)
    receptor = graph.product("receptor_branch", "immune-receptor enrichment branch",
                             MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("split", "Split expression and immune-receptor branches",
                StepKind.SPLIT, ((cdna, InputRole.TEMPLATE),), (expression, receptor),
                changes_molecule=False)
    probes = graph.input("capture_probes", "biotinylated immunoglobulin/TCR probes",
                         MaterialKind.SELECTION)
    captured = graph.product("captured", "captured full-length receptor cDNA",
                             MaterialKind.MOLECULAR,
                             state=molecule("captured receptor cDNA", (
                                 "cell barcode", "UMI", "full-length receptor transcript")))
    rejected = graph.product("rejected", "non-receptor cDNA in the capture branch",
                             MaterialKind.WASTE)
    graph.react("capture", "Hybrid-capture immunoglobulin and TCR cDNA",
                StepKind.AFFINITY_CAPTURE,
                ((receptor, InputRole.TEMPLATE), (probes, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not contain a selected immune-receptor transcript")
    motor = graph.input("motor_adapter", "Nanopore motor sequencing adapter",
                        MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "Nanopore adapter ligase", MaterialKind.CATALYST)
    long_library = graph.product("long_library", "RAGE-Seq long-read receptor library",
                                 MaterialKind.LIBRARY, state=final)
    graph.react("ligate_motor", "Attach long-read adapters to captured receptor cDNA",
                StepKind.LIGATION,
                ((captured, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (long_library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "Nanopore adapter end chemistry is kit-revision-specific; the source establishes its molecular role."))
    graph.sequence(long_library, SequencingPlan(limitation=d.sequencing))
    expression_final, expression_primers = truseq_library(
        [role("cell barcode"), role("UMI"), role("transcript cDNA")],
        "RAGE-Seq matched expression library")
    reviewed = ["cdna", "split", "capture_probes", "capture", "motor_adapter",
                "adapter_ligase", "ligate_motor", "long_library"]
    finish_pcr(graph, "expression", expression, expression_final, expression_primers,
                "matched short-read expression branch", reviewed,
                label="RAGE-Seq expression")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _chil_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    chromatin = graph.input("chromatin", "fixed chromatin containing an antibody target",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("fixed chromatin", ("antibody target", "nearby DNA")),
                            primary=True)
    antibody = graph.input("antibody", "antibody carrying the ChIL oligonucleotide",
                           MaterialKind.SELECTION)
    bound = graph.product("bound", "ChIL-probe-targeted chromatin", MaterialKind.MOLECULAR,
                          state=molecule("targeted chromatin", (
                              "antibody–ChIL probe", "nearby genomic DNA")))
    unbound = graph.product("unbound", "untargeted chromatin", MaterialKind.WASTE)
    graph.react("bind", "Bind the ChIL-oligo-conjugated antibody",
                StepKind.BINDING,
                ((chromatin, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (bound, unbound))
    graph.discard(unbound, "does not retain the selected chromatin target")
    transposome = graph.input("transposome", "loaded ChIL-compatible transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented target-proximal genomic DNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("tagged target DNA", (
                               "ChIL probe", "target-proximal genomic DNA", "Tn5 end")))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment nearby genomic DNA in situ", StepKind.TAGMENTATION,
                ((bound, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain target-proximal tagged DNA")
    extension_primer = graph.input("extension_primer", "published ChIL-probe extension architecture",
                                   MaterialKind.OLIGO,
                                   source=source + "; Supplementary Table 4")
    polymerase = graph.input("extension_polymerase", "DNA polymerase", MaterialKind.CATALYST)
    extended = graph.product("extended", "ChIL probe joined to a nearby DNA copy",
                             MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("extend", "Extend from the ChIL probe into target-proximal DNA",
                StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE), (extension_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (extended,))
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    rna = graph.product("rna", "linearly amplified target-proximal RNA copy",
                        MaterialKind.MOLECULAR, state=extended.state)
    graph.react("transcribe", "Linearly amplify with T7 transcription",
                StepKind.TRANSCRIPTION,
                ((extended, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (rna,))
    rt_primer = graph.input("library_rt_primer", "Read2+N11 first-strand synthesis primer",
                            MaterialKind.OLIGO,
                            source=source + "; Supplementary Table 4")
    rt = graph.input("library_rt_enzyme", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "ChIL target-proximal cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("ChIL cDNA", ("target-proximal cDNA",)))
    graph.react("library_rt", "Reverse-transcribe ChIL RNA with Read2+N11",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    forward = graph.input("library_forward", "Ad1_noMX PCR primer", MaterialKind.OLIGO,
                          source=source + "; Supplementary Table 4")
    reverse = graph.input("library_reverse", "Ad2_index PCR primer", MaterialKind.OLIGO,
                          source=source + "; Supplementary Table 4")
    pcr_polymerase = graph.input("library_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)
    library = graph.product("library_library", "ChIL target-proximal sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify and index the ChIL-seq library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    final, PCRPrimer("Ad1_noMX", CHIL_AD1_NOMX),
                    PCRPrimer("Ad2_index", CHIL_AD2_INDEX),
                    statement="Supplementary Table 4 prints both complete PCR primers."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    reviewed = ["chromatin", "antibody", "bind", "transposome", "tn5", "tagment",
                "extension_primer", "extension_polymerase", "extend", "t7_polymerase",
                "transcribe", "library_rt_primer", "library_rt_enzyme", "library_rt",
                "library_forward", "library_reverse", "library_polymerase",
                "library_pcr", "library_library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _cfmedip_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("cfdna", "naturally fragmented plasma cell-free DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("plasma cfDNA", ("naturally fragmented cfDNA",)), primary=True)
    repair = graph.input("repair_mix", "cfDNA end-repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired cfDNA", MaterialKind.MOLECULAR,
                             state=dna.state)
    graph.react("end_repair", "End-repair cfDNA fragments", StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "dA-tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "3-prime-dA-tailed cfDNA", MaterialKind.MOLECULAR,
                           state=molecule("dA-tailed cfDNA", ("3-prime dA", "cfDNA", "3-prime dA")))
    graph.react("tail", "Add one 3-prime dA to each fragment end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("adapters", "Illumina dT-overhang adapters", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "adapter DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "adapter-ligated cfDNA", MaterialKind.MOLECULAR,
                            state=molecule("adapter-ligated cfDNA", (
                                "P5-side adapter", "cfDNA", "P7-side adapter"),
                                junctions=(("P5-side adapter", "cfDNA", "adapter ligation"),
                                           ("cfDNA", "P7-side adapter", "adapter ligation"))))
    graph.react("ligate", "Ligate Illumina adapters to dA-tailed cfDNA", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The role-level source summary specifies dA/dT ligation without fixing one vendor adapter revision."))
    controls = graph.input("controls", "methylated/unmethylated controls and lambda filler DNA",
                           MaterialKind.MOLECULAR)
    mixture = graph.product("mixture", "low-input cfDNA library mixture",
                            MaterialKind.MOLECULAR,
                            state=molecule("library mixture", ("cfDNA library", "controls", "filler DNA")))
    graph.react("pool_controls", "Add controls and filler DNA",
                StepKind.POOL,
                ((adapted, InputRole.TEMPLATE), (controls, InputRole.TEMPLATE)), (mixture,))
    denaturing = graph.input("denaturing_condition", "single-strand denaturing condition",
                             MaterialKind.CONDITION)
    single = graph.product("single_strand", "denatured cfDNA library strands",
                           MaterialKind.MOLECULAR, state=mixture.state)
    complement = graph.product("complement", "opposite denatured strands", MaterialKind.WASTE)
    graph.react("denature", "Denature the low-input library mixture",
                StepKind.STRAND_SEPARATION,
                ((mixture, InputRole.TEMPLATE), (denaturing, InputRole.CONDITION)),
                (single, complement))
    graph.discard(complement, "represented by the retained single-stranded library population")
    antibody = graph.input("antibody", "5-methylcytosine immunoprecipitation antibody",
                           MaterialKind.SELECTION)
    enriched = graph.product("enriched", "5mC-enriched adapter-ligated cfDNA",
                             MaterialKind.MOLECULAR, state=single.state)
    depleted = graph.product("depleted", "unmethylated antibody-depleted material",
                             MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate 5mC-containing fragments",
                StepKind.AFFINITY_CAPTURE,
                ((single, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (enriched, depleted))
    graph.discard(depleted, "does not belong to the antibody-enriched methylated fraction")
    reviewed = ["cfdna", "repair_mix", "end_repair", "tailing_enzyme", "tail", "adapters",
                "ligase", "ligate", "controls", "pool_controls", "denaturing_condition",
                "denature", "antibody", "immunoprecipitate"]
    finish_pcr(graph, "library", enriched, final, primers, d.sequencing, reviewed,
                label="cfMeDIP-seq")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _ric_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    topology = FixedRNPContact("formaldehyde crosslinking")
    cells = graph.input("cells", "cells containing spatially proximal RNA–protein assemblies",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("fixed-cell input", ("spatially proximal RNAs",)), primary=True)
    formaldehyde = graph.input("crosslinker", "formaldehyde crosslinking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed", "crosslinked RNA–protein neighbourhoods",
                          MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("crosslink", "Crosslink RNA–protein assemblies in situ",
                StepKind.CROSS_LINKING,
                ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,))
    fragmentation = graph.input("fragmentation_condition", "RNA fragmentation condition",
                                MaterialKind.CONDITION)
    fragments = graph.product("fragments", "fragmented fixed RNA neighbourhoods",
                              MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("fragment", "Fragment RNA in the fixed assemblies", StepKind.FRAGMENTATION,
                ((fixed, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)), (fragments,))
    ligase = graph.input("proximity_ligase", "RNA proximity ligase", MaterialKind.CATALYST)
    chimera = graph.product("chimera", "in-situ proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=topology.chimera())
    graph.react("ligate", "Proximity-ligate labelled RNA ends in situ",
                StepKind.PROXIMITY_LIGATION,
                ((fragments, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera,))
    reversal = graph.input("reversal", "crosslink-reversal and protein-removal condition",
                           MaterialKind.CONDITION)
    purified = graph.product("purified", "de-crosslinked RNA chimeras",
                             MaterialKind.MOLECULAR, state=chimera.state)
    graph.react("reverse_crosslinks", "Reverse crosslinks and recover RNA chimeras",
                StepKind.DEPROTEINIZATION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (purified,))
    reviewed = ["cells", "crosslinker", "crosslink", "fragmentation_condition", "fragment",
                "proximity_ligase", "ligate", "reversal", "reverse_crosslinks"]
    finish_rt_pcr(graph, "library", purified, final, primers, d.sequencing, reviewed,
                   label="RIC-seq chimera")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _plac_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "crosslinked nuclei containing chromatin contacts",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("crosslinked nuclei", ("chromatin contacts",)), primary=True)
    restriction = graph.input("restriction_enzyme", "Hi-C restriction endonuclease",
                              MaterialKind.CATALYST)
    digested = graph.product("digested", "restriction-digested chromatin contacts",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    graph.react("digest", "Restriction-digest chromatin in nuclei", StepKind.DIGESTION,
                ((nuclei, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (digested,))
    fill_primer = graph.input("biotin_nucleotides", "biotinylated fill-in nucleotides",
                              MaterialKind.OLIGO)
    polymerase = graph.input("fill_polymerase", "restriction-end fill-in polymerase",
                             MaterialKind.CATALYST)
    filled = graph.product("filled", "biotin-filled restriction fragment ends",
                           MaterialKind.MOLECULAR,
                           state=molecule("filled contact ends", ("filled locus A end", "filled locus B end")))
    graph.react("fill", "Fill restriction ends with biotinylated nucleotides",
                StepKind.EXTENSION,
                ((digested, InputRole.TEMPLATE), (fill_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (filled,))
    ligase = graph.input("proximity_ligase", "in-nucleus DNA ligase", MaterialKind.CATALYST)
    contacts = graph.product("contacts", "biotinylated proximity-ligation products",
                             MaterialKind.MOLECULAR,
                             state=molecule("contact product", ("locus A", "locus B"),
                                            junctions=(("locus A", "locus B",
                                                        "biotinylated proximity ligation"),)))
    graph.react("proximity_ligation", "Proximity-ligate inside intact nuclei",
                StepKind.PROXIMITY_LIGATION,
                ((filled, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    shear = graph.input("shear_condition", "chromatin shearing condition", MaterialKind.CONDITION)
    sheared = graph.product("sheared", "sheared contact products", MaterialKind.MOLECULAR,
                            state=contacts.state)
    graph.react("shear", "Shear proximity-ligated chromatin", StepKind.FRAGMENTATION,
                ((contacts, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (sheared,))
    antibody = graph.input("antibody", "chosen chromatin-target ChIP antibody",
                           MaterialKind.SELECTION)
    chip = graph.product("chip", "protein-anchored contact products",
                         MaterialKind.MOLECULAR, state=sheared.state)
    chip_waste = graph.product("chip_waste", "non-target chromatin contacts", MaterialKind.WASTE)
    graph.react("chip_capture", "Immunoprecipitate the selected chromatin target",
                StepKind.AFFINITY_CAPTURE,
                ((sheared, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (chip, chip_waste))
    graph.discard(chip_waste, "does not contain the selected protein-anchored contact")
    streptavidin = graph.input("streptavidin", "streptavidin junction-capture beads",
                               MaterialKind.SELECTION)
    enriched = graph.product("enriched", "biotin-junction-enriched PLAC-seq contacts",
                             MaterialKind.MOLECULAR, state=chip.state)
    nonjunction = graph.product("nonjunction", "non-biotinylated ChIP material", MaterialKind.WASTE)
    graph.react("junction_capture", "Enrich biotinylated contact junctions",
                StepKind.AFFINITY_CAPTURE,
                ((chip, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (enriched, nonjunction))
    graph.discard(nonjunction, "does not retain a biotinylated proximity-ligation junction")
    reviewed = ["nuclei", "restriction_enzyme", "digest", "biotin_nucleotides",
                "fill_polymerase", "fill", "proximity_ligase", "proximity_ligation",
                "shear_condition", "shear", "antibody", "chip_capture", "streptavidin",
                "junction_capture"]
    finish_pcr(graph, "library", enriched, final, primers, d.sequencing, reviewed,
                label="PLAC-seq contact")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


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


def get(key: str, *, definition: Definition | None = None,
        psoralen: tuple[str, str, str, str] | None = None):
    """Create a page spec; caller-owned definitions supply named protocol claims."""
    d = definition if definition is not None else DEFS[key]
    primers = ()
    if key == "codec":
        final = codec_final_library()
        primers = sp.standard_primers(final, ("Read 1", "Read 2"))
    elif key == "chil-seq":
        final, primers = nextera_library(_payload(d), d.title + " final library",
                                          dual_index=False)
    elif key == "nano3p-seq":
        final = nano3p_final_library()
    elif d.library == "truseq":
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
    spec = SimpleNamespace(
        TITLE=d.title, NOTES="01_" + key.replace("-", "_") + ".html",
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.source_label}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
        SEQ_PRIMERS=primers, FINAL_CAPTION=d.sequencing,
        SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.sequencing,
        SEQUENCING_UNAVAILABLE=("The defining source does not expose an exact sequencing-primer sequence; the molecular read entry is shown by role."
                                if not primers else ""),
        workflow=(lambda: psoralen_workflow(*psoralen[:2])) if psoralen else
                 (lambda: workflow_for(key)),
        oligos=lambda: (),
    )
    if key == "codec":
        spec.oligos = codec_oligos
    elif key == "paris":
        spec.oligos = paris_oligos
        spec.SEQUENCING_UNAVAILABLE = (
            "The custom sequencing-primer bases are published; exact placement remains "
            "unresolved only because the circular-cDNA junction is not yet represented.")
    elif key == "mario":
        spec.oligos = mario_oligos
    elif key == "rage-seq":
        spec.oligos = rage_seq_oligos
        spec.SEQUENCING_MECHANISM = "nanopore_motor"
        spec.MOTOR_ADAPTER_NAME = "SQK-LSK108 / SQK-LSK308 motor adapter"
    elif key == "seq-scope":
        spec.oligos = seqscope_oligos
    elif key == "bliss":
        spec.oligos = bliss_oligos
    elif key == "chil-seq":
        spec.oligos = chil_oligos
    elif key == "nano3p-seq":
        spec.oligos = nano3p_oligos
        spec.SEQUENCING_MECHANISM = "nanopore_motor"
        spec.MOTOR_ADAPTER_NAME = "SQK-DCS109 AMX motor adapter"
    elif key == "stlfr":
        spec.oligos = stlfr_oligos
        # Primer bases are now source-exact and listed above.  Molecular placement in
        # the DNB circle remains role-level until its circle-closing geometry is known.
        spec.SEQUENCING_UNAVAILABLE = (
            "The BGI sequencing-primer sequences are published; only the DNB "
            "circle-closing adapter or splint geometry remains unavailable.")
    factories = {
        "paris": lambda: psoralen_graph(d, final, primers,
                                         reagent="AMT psoralen",
                                         selection_label="two-dimensional-gel duplex selection"),
        "ligr-seq": lambda: psoralen_graph(d, final, primers,
                                            reagent="AMT psoralen",
                                            selection_label="RNase-R-resistant duplex selection"),
        "mario": lambda: _mario_graph(d, final, primers),
        "codec": lambda: _codec_graph(d, final, primers),
        "nanoseq": lambda: _nanoseq_graph(d, final, primers),
        "stlfr": lambda: _stlfr_graph(d, final),
        "nano3p-seq": lambda: _nano3p_graph(d, final),
        "maester": lambda: _maester_graph(d, final, primers),
        "malbac-dt": lambda: _malbac_dt_graph(d, final, primers),
        "seq-scope": lambda: _seq_scope_graph(d, final, primers),
        "bliss": lambda: _bliss_graph(d, final, primers),
        "scsplat": lambda: _scsplat_graph(d, final, primers),
        "dart-seq": lambda: _dart_graph(d, final, primers),
        "snubar-atac": lambda: _snubar_graph(d, final, primers),
        "redeem": lambda: _redeem_graph(d, final, primers),
        "rage-seq": lambda: _rage_graph(d, final),
        "chil-seq": lambda: _chil_graph(d, final, primers),
        "cfmedip-seq": lambda: _cfmedip_graph(d, final, primers),
        "ric-seq": lambda: _ric_graph(d, final, primers),
        "plac-seq": lambda: _plac_graph(d, final, primers),
    }
    if key in factories:
        spec.protocol_graph = factories[key]
    return spec
