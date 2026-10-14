"""Reusable constructions for RNA end-capture and modification-stop libraries."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
import illumina as il
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role
import seqprimers as sp


@dataclass(frozen=True)
class StopReadout:
    reagent: str
    target: str
    control: str

    def __post_init__(self):
        if not self.reagent or not self.target or not self.control:
            raise ValueError("a modification-stop readout requires reagent, target and control")


@dataclass(frozen=True)
class BaseConversionReadout:
    reagent: str
    converted_base: str
    product: str
    read_as: str
    protected_base: str

    def __post_init__(self):
        if self.converted_base == self.protected_base:
            raise ValueError("converted and protected bases must be distinct chemical states")
        if not all((self.reagent, self.product, self.read_as)):
            raise ValueError("base-conversion chemistry must declare reagent, product and readout")


@dataclass(frozen=True)
class SelectiveCleavageReadout:
    target_modifications: tuple[str, ...]
    cleavage_reagent: str
    selected_end: str
    signal_offset: int

    def __post_init__(self):
        if self.selected_end != "5-prime phosphate":
            raise ValueError("adapter-positive cleavage readout must select a newly made 5-prime phosphate")
        if self.signal_offset != 1 or not self.target_modifications:
            raise ValueError("cleavage signal must be the N+1 fragment start for declared targets")


OTTR_PRIMER = Segment("universal +1T DNA primer", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT", "r2")
OTTR_DUPLEX_RNA = Segment("blocked RNA complement", "GATCGGAAGAGCACACGTCTGAACTCCAGTC", "r2")
OTTR_ACCEPTOR = Segment("universal adaptor template", "ACACTCTTTCCCTACACGACGCTCTTCCGATCTC", "r1")
PSEUDO_ADAPTER = Segment("preadenylated 3-prime adapter", "TGGAATTCTCGGGTGCCAAGG", "r2")
PSEUDO_RT = Segment("Pseudo-seq RT primer", "GATCGTCGGACTGTAGAACTCTGAACCTGTCGGTGGTCGCCGTATCATTCACTCAGCCTTGGCACCCGAGAATTCCA", "r2")
PSEUDO_FWD = Segment("PCR forward", "AATGATACGGCGACCACCGA", "p5")
PSEUDO_REV = Segment("indexed PCR reverse", "CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", "p7", placeholder=True,
                     feature=feature("sample_index", "sample_index", "unknown"))
DMS_LINKER = Segment("microRNA cloning linker-1", "CTGTAGGCACCATCAAT", "r2")
SS2_HAIRPIN = Segment("Structure-seq2 hairpin donor", "TGAAGAGCCTAGTCGCTGTTCANNNNNNCTGCCCATAGAG", "r1", placeholder=True)
SS2_RT_SEQUENCE = "CAGACGTGTGCTCTTCCGATC" + "N" * 6
SS2_FORWARD_SEQUENCE = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGAACAGCGACTAGGCTCTTCA"
SS2_REVERSE_SEQUENCE = "CAAGCAGAAGACGGCATACGAGAT" + "N" * 6 + "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"
SS2_FORWARD = Segment("Structure-seq2 TruSeq forward", SS2_FORWARD_SEQUENCE, "p5")
SS2_REVERSE = Segment("Structure-seq2 indexed reverse", SS2_REVERSE_SEQUENCE, "p7", placeholder=True,
                      feature=feature("sample_index", "sample_index", "unknown"))
SS2_CUSTOM = Segment("Structure-seq2 custom sequencing primer", "TCTTCCGATCTTGAACAGCGACTAGGCTCTTCA", "r1")
SS2_PCR_TEMPLATE = Construct((
    Segment("hairpin forward-primer site", SS2_FORWARD_SEQUENCE[-22:], "r1"),
    Segment("stop-defined cDNA", "X" * 42, placeholder=True),
    Segment("RT-primer-derived reverse site", revcomp(SS2_REVERSE_SEQUENCE[-21:]), "r2"),
), name="Structure-seq2 PCR-ready template")
SS2_PCR_GEOMETRY = PCRGeometry.exact(
    SS2_PCR_TEMPLATE,
    PCRPrimer("Structure-seq2 TruSeq forward", SS2_FORWARD_SEQUENCE),
    PCRPrimer("Structure-seq2 indexed reverse", SS2_REVERSE_SEQUENCE, anneal_nt=21),
    statement=("The main Methods print both PCR primers, and Supplementary Figure S1 "
               "draws their hairpin- and RT-primer-derived landing regions."))
SS2_FINAL_LIBRARY = Construct((
    Segment("forward PCR primer", SS2_FORWARD_SEQUENCE, "p5"),
    Segment("stop-defined cDNA", "X" * 42, placeholder=True),
    Segment("Index 1 primer site", revcomp(SS2_REVERSE_SEQUENCE[-33:]), "r2"),
    Segment("sample index complement", "N" * 6, "cbc", placeholder=True,
            feature=feature("sample_index", "sample_index", "unknown")),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="Structure-seq2 six-base-index library")
SS2_SEQ_PRIMERS = (
    sp.custom("Read 1", "Structure-seq2 custom Read 1", SS2_CUSTOM.top,
              "Ritchey et al. Supplementary Figure S1"),
    sp.TRUSEQ["I1"],
)

HRF_RT_RANDOM = "AGACGTGTGCTCTTCCGATCT" + "N" * 8 + "S"
HRF_RT_CASSETTE = "AGACGTGTGCTCTTCCGATCTGAACCGGACCGAAGCCCG"
HRF_LIGATION_ADAPTER = "N" * 7 + "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
HRF_PCR_FORWARD = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT"
HRF_INDEXES = ("GGAACT", "GGACGG", "CGTACG", "GCTACC")


def _hrf_reverse(index: str) -> str:
    return il.P7 + index + il.TRUSEQ_READ2


HRF_PCR_TEMPLATE = Construct((
    Segment("ligation-adapter forward site", HRF_PCR_FORWARD[-19:], "r1"),
    Segment("accessibility-defined cDNA", "X" * 42, placeholder=True),
    Segment("random-primer coordinate", "N" * 8 + "S", placeholder=True),
    Segment("RT-handle reverse site", revcomp(_hrf_reverse(HRF_INDEXES[0])[-21:]), "r2"),
), name="HRF-Seq PCR-ready template")
HRF_PCR_GEOMETRY = PCRGeometry.exact(
    HRF_PCR_TEMPLATE,
    PCRPrimer("PCR_forward", HRF_PCR_FORWARD, anneal_nt=19),
    PCRPrimer("PCR_REVERSE_INDEX.14_AGTTCC", _hrf_reverse(HRF_INDEXES[0]),
              anneal_nt=21),
    statement=("Supplementary Table 1 prints the random RT primer, ligation adapter, "
               "forward primer and all four indexed reverse primers."))
HRF_FINAL_LIBRARY = Construct((
    Segment("P5 / partial Read 1 PCR primer", HRF_PCR_FORWARD, "p5"),
    Segment("remaining Read 1 site", il.TRUSEQ_READ1[-10:], "r1"),
    Segment("7-nt molecular barcode", "N" * 7, "umi", placeholder=True,
            feature=feature("hrf_molecular_barcode", "umi", "random")),
    Segment("accessibility-defined cDNA", "X" * 42, placeholder=True),
    Segment("random-primer coordinate", "N" * 8 + "S", placeholder=True),
    Segment("Index 1 primer site", revcomp(il.TRUSEQ_READ2), "r2"),
    Segment("sample index", revcomp(HRF_INDEXES[0]), "cbc",
            feature=feature("hrf_sample_index", "sample_index", "whitelist",
                            whitelist="Kielpinski and Vinther Supplementary Table 1")),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="HRF-Seq representative indexed library")
HRF_SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])

SHAPE_A_T = il.P5 + il.TRUSEQ_READ1
SHAPE_A_B = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
SHAPE_HANDLE_FEATURE = feature("shape_treatment_handle", "feature_barcode", "whitelist",
                               whitelist="RRRY treated / YYYR control handle families")
SHAPE_RNA_FEATURE = feature("shape_rna_identity", "feature_barcode", "whitelist",
                            whitelist="designed RNA identity barcodes")
SHAPE_FINAL_LIBRARY = Construct((
    Segment("P5", il.P5, "p5"),
    Segment("A_adapter_t / Read 1 site", il.TRUSEQ_READ1, "r1"),
    Segment("plus/minus handle", "RRRY", "cbc", placeholder=True,
            feature=SHAPE_HANDLE_FEATURE),
    role("target-specific RT-primer site"),
    role("SHAPE-stop-defined cDNA"),
    role("RNA identity barcode", tag="cbc", length_bp=4,
         molecular_feature=SHAPE_RNA_FEATURE),
    Segment("A_adapter_b sequencing arm", SHAPE_A_B[:-len(il.P7_RC)], "r2"),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="SHAPE-Seq v1 library")


def _ottr() -> Workflow:
    return chain(molecule("input small RNA", (role("RNA of any end sequence"),)), (
        ("Cap each input 3′ end with one ddA and/or ddG",
         molecule("3′-blocked input", (role("input RNA"), role("single nonextendable ddR"))),
         "BoMoC terminal-transferase activity runs in Mn²⁺; phosphatase then destroys free ddRTP."),
        ("Switch BoMoC from the +1Y primer duplex onto input RNA",
         molecule("first template-jump product", (OTTR_PRIMER, role("input-derived cDNA"))),
         "A +1T/+1C DNA end pairs with the standardized ddA/ddG RNA end."),
        ("Add one dG and relay onto the 3′C adapter template",
         molecule("dual-adapter cDNA", (OTTR_PRIMER, role("input-derived cDNA"), OTTR_ACCEPTOR)),
         "The second, ordered template jump copies one terminating adapter template."),
        ("PCR-add flow-cell ends and indexes",
         molecule("OTTR-seq library", (role("P5 and i5"), role("R1 adapter"), role("input cDNA"),
                  role("R2 adapter"), role("i7 and P7"))),
         "Universal-adapter libraries use four to eight PCR cycles; full-length adapters can omit PCR."),
    ))


def _pseudo() -> Workflow:
    readout = StopReadout("CMC", "pseudouridine", "parallel mock −CMC")
    return chain(molecule("fragmented RNA", (role("RNA containing U and Ψ"),)), (
        ("Derivatize RNA with CMC and reverse U/G adducts",
         molecule("CMC-treated RNA", (role("RNA"), role("persistent Ψ–CMC adduct"))),
         f"{readout.control} distinguishes Ψ-dependent stops from background."),
        ("Repair ends, size-select, and ligate the preadenylated 3′ adapter",
         molecule("adapter-ligated RNA", (role("CMC-treated RNA fragment"), PSEUDO_ADAPTER),
                  junctions=(("CMC-treated RNA fragment", "preadenylated 3-prime adapter", "RNA ligation"),)),
         "CIP/PNK treatment precedes selection of 80–140-nt RNA fragments."),
        ("Reverse-transcribe until AMV RT stops one nucleotide 3′ of Ψ–CMC",
         molecule("truncated cDNA", (PSEUDO_RT, role("cDNA ending at Ψ-dependent stop"))),
         "The 5′ end of each selected cDNA reports the position of the block."),
        ("Circularize cDNA and amplify with indexed primers",
         molecule("Pseudo-seq library", (PSEUDO_FWD, role("stop-defined cDNA"), PSEUDO_REV)),
         "CircLigase closes the single-stranded cDNA before Illumina PCR."),
    ))


def _dms_seq() -> Workflow:
    return chain(molecule("poly(A)-selected RNA", (role("native cellular RNA structure"),)), (
        ("Probe accessible A and C with DMS",
         molecule("DMS-modified RNA", (role("RNA"), role("methylated exposed A/C"))),
         "Parallel in-vivo, refolded in-vitro and denatured samples separate cellular structure from intrinsic reactivity."),
        ("Fragment with Zn²⁺, select 60–80 nt, and repair 3′ phosphates",
         molecule("end-repaired DMS RNA fragment", (role("60–80-nt modified RNA"), role("3-prime hydroxyl"))),
         "T4 PNK resolves the 3′ phosphates left by chemical fragmentation."),
        ("Ligate preadenylated cloning linker-1",
         molecule("linker-ligated DMS RNA", (role("modified RNA fragment"), DMS_LINKER),
                  junctions=(("modified RNA fragment", "microRNA cloning linker-1", "Rnl2tr K227Q ligation"),)),
         "The linker is 5′ preadenylated and 3′ ddC-blocked, so ligation proceeds without ATP."),
        ("Reverse-transcribe and isolate prematurely terminated cDNA",
         molecule("DMS-stop cDNA", (role("adapter-primed cDNA ending before modified A/C"),)),
         "A 25–50-nt cDNA size selection enriches Superscript III termination products."),
        ("Circularize cDNA and add Illumina adapters by PCR",
         molecule("DMS-seq library", (role("P5"), role("DMS-stop cDNA"), role("sample index"), role("P7"))),
         "CircLigase precedes eight to ten PCR cycles and HiSeq 2000 sequencing."),
    ))


def _structure_seq2() -> Workflow:
    return chain(molecule("poly(A)-selected RNA", (role("RNA folded in vivo"),)), (
        ("Probe exposed A and C in vivo with DMS",
         molecule("DMS-probed RNA", (role("RNA with structure-dependent DMS adducts"),)),
         "A matched −DMS sample measures background reverse-transcription stops."),
        ("Random-prime reverse transcription to each DMS stop",
         molecule("stop-defined cDNA", (role("TruSeq handle plus random hexamer"), role("cDNA ending before modified base"))),
         "The optional branch incorporates biotin-dCTP for streptavidin purification instead of two PAGE steps."),
        ("Ligate the degenerate hairpin donor with T4 DNA ligase",
         molecule("hairpin-adapted cDNA", (SS2_HAIRPIN, role("stop-defined cDNA")),
                  junctions=(("Structure-seq2 hairpin donor", "stop-defined cDNA", "T4 DNA ligase"),)),
         "Six random bases let the donor base-pair at diverse cDNA ends and reduce ligation bias."),
        ("PCR with TruSeq forward and indexed reverse primers",
         molecule("Structure-seq2 library", (SS2_FORWARD, role("stop-defined cDNA"), SS2_REVERSE)),
         "A final 220–600-nt PAGE selection removes primer dimers and residual by-product."),
    ))


def _ribomethseq() -> Workflow:
    return chain(molecule("total or enriched RNA", (role("RNA containing 2-prime-O-methylated riboses"),)), (
        ("Partially hydrolyse RNA under mild alkaline conditions",
         molecule("alkaline RNA fragments", (role("random fragments with gaps beside protected Nm sites"),)),
         "A 2′-O-methyl ribose protects its 3′-adjacent phosphodiester bond from cleavage."),
        ("Convert cyclic-phosphate/5′-OH ends to ligatable termini",
         molecule("end-repaired RNA fragments", (role("5-prime phosphate"), role("RNA fragment"), role("3-prime hydroxyl"))),
         "Antarctic phosphatase followed by T4 PNK creates the ends required by the standard small-RNA kit."),
        ("Ligate 3′ and 5′ small-RNA adapters",
         molecule("dual-adapter RNA fragments", (role("NEBNext 5-prime small-RNA adapter"), role("hydrolysis fragment"), role("NEBNext 3-prime small-RNA adapter"))),
         "Direct ligation preserves the exact hydrolysis endpoints used for methylation scoring."),
        ("Reverse-transcribe, index by PCR and sequence fragment ends",
         molecule("RiboMethSeq library", (role("P5"), role("hydrolysis-end cDNA"), role("sample index"), role("P7"))),
         "Paired-end MiSeq can recover both ends; SR50 HiSeq analysis can use the less biased 5′-end profile."),
    ))


def _pars() -> Workflow:
    return chain(molecule("poly(A)-selected RNA", (role("purified transcriptome"),)), (
        ("Refold RNA in vitro and split it into matched aliquots",
         molecule("refolded transcriptome", (role("RNA secondary structures"),)),
         "The two branches start from the same refolded RNA pool."),
        ("Digest double-stranded or single-stranded regions",
         molecule("paired V1/S1 digests", (role("RNase V1: double-strand cleavage"), role("RNase S1: single-strand cleavage"))),
         "RNase V1 and S1 generate separate structure-selective 5′-phosphate endpoints."),
        ("Randomly fragment, then ligate the 5′ adapter selectively",
         molecule("cleavage-selected RNA", (role("5-prime adapter"), role("RNA downstream of nuclease cut"))),
         "Random fragmentation yields 5′-OH ends; only the nuclease-created 5′-P ends accept the first adapter."),
        ("Dephosphorylate 3′ ends and ligate the second adapter",
         molecule("dual-adapter PARS RNA", (role("5-prime adapter"), role("structure-defined RNA fragment"), role("3-prime adapter"))),
         "Alkaline phosphatase makes the second ligation compatible with nuclease and fragmentation end chemistries."),
        ("Reverse-transcribe, PCR-amplify and sequence both branches",
         molecule("PARS library pair", (role("SOLiD P1 end"), role("V1 or S1 cleavage-defined cDNA"), role("SOLiD P2 end"))),
         "The log2 normalized V1/S1 read-start ratio is the PARS score for the upstream nucleotide."),
    ))


def _fragseq() -> Workflow:
    return chain(molecule("purified nuclear RNA", (role("native folded RNA"),)), (
        ("Refold RNA and split treated and undigested-control aliquots",
         molecule("refolded RNA", (role("RNA secondary structure"),)),
         "The control measures fragmentation and library-preparation background."),
        ("Partially digest exposed single-stranded RNA with nuclease P1",
         molecule("P1 accessibility fragments", (role("5-prime phosphate at accessible site"),
                  role("downstream RNA fragment"))),
         "Nuclease P1 preferentially cleaves single-stranded nucleic acid; fragment starts encode accessibility."),
        ("Size-select and ligate small-RNA sequencing adapters",
         molecule("adapter-ligated FragSeq RNA", (role("5-prime sequencing adapter"),
                  role("P1-cleavage-defined RNA"), role("3-prime sequencing adapter"))),
         "Treated and control branches are converted independently so their endpoint enrichments can be compared."),
        ("Reverse-transcribe, amplify and sequence fragment ends",
         molecule("FragSeq library", (role("SOLiD P1 end"),
                  role("P1-cleavage-defined cDNA"), role("SOLiD P2 end"))),
         "Mapped fragment termini provide nucleotide-resolution nuclease-accessibility scores."),
    ))


def _hrf_seq() -> Workflow:
    return chain(molecule("folded RNA pool", (role("RNA tertiary structure"),)), (
        ("Expose RNA to hydroxyl radicals alongside an untreated control",
         molecule("hydroxyl-radical-cleaved RNA", (role("intact RNA regions"),
                  role("backbone break at solvent-accessible ribose"))),
         "Backbone cleavage reports solvent accessibility rather than base identity."),
        ("Random-prime reverse transcription from primers carrying a 5-prime Illumina tail",
         molecule("break-terminated cDNA", (role("Illumina-tailed random primer"),
                  role("cDNA ending at radical-induced break"))),
         "The cDNA endpoint records the probing position; random priming records the opposite end."),
        ("Ligate a seven-base molecular-barcode adapter to each cDNA 3-prime end",
         molecule("barcoded HRF cDNA", (role("Illumina-tailed cDNA"),
                  role("7-nt molecular barcode", tag="umi", length_bp=7,
                       molecular_feature=feature("umi", "umi", "random")),
                  role("ligation adapter"))),
         "Unique barcode counting corrects PCR and fragment-length sampling bias."),
        ("PCR-add indexes, size-select and paired-end sequence",
         molecule("HRF-Seq library", (role("P5"), role("random-primer coordinate"),
                  role("accessibility-defined cDNA"), role("7-nt molecular barcode", tag="umi",
                       length_bp=7, molecular_feature=feature("umi", "umi", "random")),
                  role("sample index"), role("P7"))),
         "The paired reads recover the priming position and the hydroxyl-radical termination position."),
    ))


def _shape_seq() -> Workflow:
    return chain(molecule("barcoded in-vitro RNA pool", (role("RNA-specific 3-prime DNA barcode"),
                  role("folded RNA"))), (
        ("Modify flexible riboses with a SHAPE reagent and retain a minus-reagent control",
         molecule("SHAPE-modified RNA", (role("RNA"),
                  role("2-prime-O-adduct at flexible nucleotide"))),
         "The pre-attached RNA barcode lets many defined RNAs share the chemical and library workflow."),
        ("Primer-extend until reverse transcriptase reaches a SHAPE adduct",
         molecule("SHAPE-stop cDNA", (role("common primer handle"),
                  role("cDNA ending before modified nucleotide"), role("RNA identity barcode"))),
         "Treated and control stop distributions are interpreted jointly rather than as presence/absence calls."),
        ("Ligate the second sequencing adapter to cDNA",
         molecule("dual-adapter SHAPE cDNA", (role("paired-end adapter"),
                  role("stop-defined cDNA"), role("RNA identity barcode"),
                  role("paired-end adapter"))),
         "The two ends separately identify the assayed RNA and the primer-extension stop."),
        ("PCR-amplify and paired-end sequence",
         molecule("SHAPE-Seq library", (role("P5"), role("SHAPE-stop coordinate"),
                  role("RNA identity barcode"), role("P7"))),
         "Maximum-likelihood comparison of plus- and minus-reagent channels yields nucleotide reactivities."),
    ))


def _protein_pip_seq() -> Workflow:
    return chain(molecule("cross-linked cellular ribonucleoproteins", (role("RNA"),
                  role("bound proteins"))), (
        ("Split material into footprint and structure-control branches",
         molecule("four PIP-seq branches", (role("protein-bound ssRNase branch"),
                  role("protein-bound dsRNase branch"), role("protein-removed ssRNase branch"),
                  role("protein-removed dsRNase branch"))),
         "Proteinase K precedes nuclease only in controls; footprint branches retain protein during digestion."),
        ("Digest exposed single- or double-stranded RNA",
         molecule("RNase-protected fragments", (role("ssRNase-resistant RNA"),
                  role("dsRNase-resistant RNA"))),
         "The paired nuclease specificities distinguish RNA structure while the branch order reveals protein protection."),
        ("Reverse cross-links, purify RNA and prepare strand-specific libraries",
         molecule("strand-specific PIP-seq libraries", (role("5-prime sequencing adapter"),
                  role("RNase-protected RNA-derived cDNA"), role("3-prime sequencing adapter"))),
         "Rehybridization and duplex-specific nuclease normalization reduce abundant-transcript domination."),
        ("Sequence and compare footprint with protein-depleted controls",
         molecule("PIP-seq library quartet", (role("P5"),
                  role("branch-specific protected-fragment cDNA"), role("sample index"), role("P7"))),
         "Enrichment in footprint libraries marks protein-protected sites; control branches report RNA accessibility."),
    ))


def _glori() -> Workflow:
    readout = BaseConversionReadout("glyoxal plus nitrite", "adenosine", "inosine", "G", "m6A")
    return chain(molecule("fragmented poly(A) RNA", (role("RNA containing A and m6A"),)), (
        ("Protect guanosines with glyoxal",
         molecule("glyoxal-protected RNA", (role("RNA"), role("protected G"), role("A and m6A exposed"))),
         "G protection suppresses unwanted guanosine deamination during the nitrite reaction."),
        ("Deaminate unmethylated A with nitrite",
         molecule("GLORI-converted RNA", (role(f"{readout.converted_base} converted to {readout.product}"), role(f"{readout.protected_base} retained as A"))),
         "The defining reaction uses 750 mM nitrite at 16 °C for 8 h followed by 4 °C overnight."),
        ("Heat-deprotect guanosines and repair RNA ends",
         molecule("deprotected converted RNA", (role("A-to-I converted RNA with intact m6A"), role("ligatable termini"))),
         "Deprotection at 95 °C restores G; T4 PNK prepares fragmented RNA for adapter ligation."),
        ("Ligate the modified-eCLIP 3-prime RNA adapter",
         molecule("GLORI adapter-ligated RNA", (role("converted RNA"),
                  Segment("published 3-prime RNA adapter", "AGATCGGAAGAGCGTCGTG", "r2"))),
         "The 2024 detailed protocol prints the pre-adenylated RNA adapter sequence in Supplementary Table 1."),
        ("Reverse-transcribe with the published RT primer",
         molecule("GLORI cDNA", (role("converted cDNA"),
                  Segment("published RT primer", "ACACGACGCTCTTCCGATCT", "r2"))),
         "The 2024 detailed protocol prints this RT primer in Supplementary Table 1."),
        ("Ligate the ten-base-UMI DNA adapter to cDNA",
         molecule("UMI-tagged GLORI cDNA", (role("converted cDNA"), role("10-nt UMI"),
                  Segment("published 3-prime cDNA adapter arm", "AGATCGGAAGAGCACACGTCTG", "r2"))),
         "The 2024 detailed protocol prints a 5-prime-phosphorylated DNA adapter with ten degenerate UMI bases."),
        ("Index-amplify and paired-end sequence",
         molecule("GLORI library", (role("P5"), role("10-nt UMI"), role("converted cDNA"), role("sample index"), role("P7"))),
         f"Unmodified A is read as {readout.read_as}; m6A remains A, enabling site-level stoichiometry."),
    ))


def _alkaniline() -> Workflow:
    readout = SelectiveCleavageReadout(("m7G", "m3C", "D"), "aniline", "5-prime phosphate", 1)
    targets = ", ".join(readout.target_modifications)
    return chain(molecule("total or enriched RNA", (role(f"RNA containing {targets}"),)), (
        ("Partially hydrolyse RNA in bicarbonate at 96 °C",
         molecule("alkali-treated RNA", (role("size-compatible RNA fragments"), role("abasic sites at alkali-labile modifications"))),
         "Ring opening and depurination create abasic intermediates at susceptible modified bases."),
        ("Remove every pre-existing terminal phosphate",
         molecule("dephosphorylated RNA", (role("5-prime hydroxyl"), role("RNA with abasic sites"), role("3-prime hydroxyl"))),
         "Antarctic phosphatase erases background 5′-P and fragmentation-derived 3′ phosphates."),
        ("Cleave abasic sites with 1 M aniline",
         molecule("aniline-cleaved RNA", (role("upstream cleavage fragment"), role("downstream fragment with new 5-prime phosphate"))),
         "Fifteen minutes at 60 °C and pH 4.5 creates the sole ligation-competent 5′ ends."),
        ("Ligate NEBNext small-RNA adapters",
         molecule("selected AlkAniline RNA", (
             Segment("NEBNext 5-prime SR adapter", "GTTCAGAGTTCTACAGTCCGACGATC", "r1"),
             role("RNA beginning N+1 after modification"),
             Segment("NEBNext 3-prime SR adapter", "AGATCGGAAGAGCACACGTCT", "r2"))),
         "The 5′ adapter positively selects the aniline-generated phosphate; the read start is N+1. The current NEB E7300/E7330 manual identifies both adapter sequences."),
        ("Reverse-transcribe, barcode by PCR and single-end sequence",
         molecule("AlkAniline-Seq library", (role("P5"), role("cleavage-defined cDNA"), role("sample index"), role("P7"))),
         "Only read 5′ ends are required for normalized-cleavage scoring."),
    ))


def _finish_pcr(graph, prefix, template, final, endpoint, reviewed, *, label,
                geometry=None, primers=(), acquired_roles=(), read_lengths=None,
                read_uses=None):
    forward = graph.input(prefix + "_forward", label + " forward PCR primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " indexed reverse PCR primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Amplify and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=(geometry or PCRGeometry.unavailable(
                    final, reason="The source establishes the readout roles but not one universal source-exact amplicon geometry.")))
    if primers:
        graph.sequence(library, SequencingPlan(
            library=final, primers=tuple(primers),
            acquired_roles=tuple(acquired_roles),
            read_lengths=read_lengths or {}, read_uses=read_uses or {}))
    else:
        graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library


def _apply_linear_step(graph, current, identifier, label, kind, reviewed):
    inputs = [(current, InputRole.TEMPLATE)]
    geometry = None
    if kind is StepKind.LIGATION:
        reagent = graph.input(identifier + "_adapter", label + " adapter or ligating oligo",
                              MaterialKind.ADAPTER)
        enzyme = graph.input(identifier + "_ligase", label + " ligase",
                             MaterialKind.CATALYST)
        inputs.extend(((reagent, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)))
        geometry = AdapterLigationGeometry.unavailable(
            "The source defines the ligated molecular role; the role-level shared page does not assert one universal end geometry.")
        reviewed.extend((identifier + "_adapter", identifier + "_ligase"))
    elif kind is StepKind.REVERSE_TRANSCRIPTION:
        primer = graph.input(identifier + "_primer", label + " RT primer", MaterialKind.OLIGO)
        enzyme = graph.input(identifier + "_enzyme", label + " reverse transcriptase",
                             MaterialKind.CATALYST)
        inputs.extend(((primer, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)))
        reviewed.extend((identifier + "_primer", identifier + "_enzyme"))
    elif kind is StepKind.TEMPLATE_SWITCHING:
        oligo = graph.input(identifier + "_oligo", label + " template oligo",
                            MaterialKind.OLIGO)
        enzyme = graph.input(identifier + "_enzyme", label + " template-switching polymerase",
                             MaterialKind.CATALYST)
        inputs.extend(((oligo, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)))
        reviewed.extend((identifier + "_oligo", identifier + "_enzyme"))
    elif kind in {StepKind.DIGESTION, StepKind.END_REPAIR, StepKind.PHOSPHORYLATION,
                  StepKind.DEPHOSPHORYLATION, StepKind.CIRCULARIZATION, StepKind.TAILING}:
        enzyme = graph.input(identifier + "_enzyme", label + " enzyme", MaterialKind.CATALYST)
        inputs.append((enzyme, InputRole.CATALYST))
        reviewed.append(identifier + "_enzyme")
        if kind is StepKind.CIRCULARIZATION:
            geometry = CircularizationGeometry.unavailable(
                "The source establishes cDNA circularization; this shared role-level model does not assert one universal closing-end geometry.")
    elif kind is StepKind.FRAGMENTATION:
        condition = graph.input(identifier + "_condition", label + " fragmentation condition",
                                MaterialKind.CONDITION)
        inputs.append((condition, InputRole.CONDITION))
        reviewed.append(identifier + "_condition")
    elif kind in {StepKind.CHEMICAL_CONVERSION, StepKind.THERMAL_TREATMENT,
                  StepKind.DEPROTEINIZATION}:
        condition = graph.input(identifier + "_condition", label + " condition",
                                MaterialKind.CONDITION)
        inputs.append((condition, InputRole.CONDITION))
        reviewed.append(identifier + "_condition")
    else:
        raise ValueError(f"unsupported linear RNA-end step {kind}")
    product = graph.product(identifier + "_product", label + " product",
                            MaterialKind.MOLECULAR,
                            state=molecule(label + " product", (role("retained molecular product"),)))
    graph.react(identifier, label, kind, inputs, (product,), geometry=geometry)
    reviewed.append(identifier)
    return product


_LINEAR_RECIPES = {
    "ottr-seq": (
        ("cap_end", "Cap each input 3-prime end with one ddA and/or ddG", StepKind.TAILING),
        ("first_jump", "Template-switch from the standardized primer duplex onto input RNA", StepKind.TEMPLATE_SWITCHING),
        ("second_jump", "Relay onto the terminating adapter template", StepKind.TEMPLATE_SWITCHING),
    ),
    "pseudo-seq": (
        ("cmc", "Derivatize pseudouridine with persistent CMC", StepKind.CHEMICAL_CONVERSION),
        ("repair", "Repair RNA fragment ends", StepKind.END_REPAIR),
        ("ligate", "Ligate the preadenylated 3-prime adapter", StepKind.LIGATION),
        ("rt", "Reverse-transcribe to the pseudouridine-dependent stop", StepKind.REVERSE_TRANSCRIPTION),
        ("circle", "Circularize the stop-defined cDNA", StepKind.CIRCULARIZATION),
    ),
    "dms-seq": (
        ("probe", "Probe accessible A and C with DMS", StepKind.CHEMICAL_CONVERSION),
        ("fragment", "Fragment DMS-modified RNA", StepKind.FRAGMENTATION),
        ("repair", "Repair fragmentation-derived RNA ends", StepKind.PHOSPHORYLATION),
        ("ligate", "Ligate the preadenylated cloning linker", StepKind.LIGATION),
        ("rt", "Reverse-transcribe to DMS-dependent stops", StepKind.REVERSE_TRANSCRIPTION),
        ("circle", "Circularize the stop-defined cDNA", StepKind.CIRCULARIZATION),
    ),
    "structure-seq2": (
        ("probe", "Probe exposed A and C in vivo with DMS", StepKind.CHEMICAL_CONVERSION),
        ("rt", "Random-prime reverse transcription to DMS stops", StepKind.REVERSE_TRANSCRIPTION),
        ("ligate", "Ligate the degenerate hairpin donor", StepKind.LIGATION),
    ),
    "ribomethseq": (
        ("hydrolyse", "Partially hydrolyse RNA under mild alkaline conditions", StepKind.CHEMICAL_CONVERSION),
        ("dephosphorylate", "Remove incompatible cyclic-phosphate termini", StepKind.DEPHOSPHORYLATION),
        ("phosphorylate", "Create ligatable 5-prime-phosphate termini", StepKind.PHOSPHORYLATION),
        ("ligate_3", "Ligate the 3-prime small-RNA adapter", StepKind.LIGATION),
        ("ligate_5", "Ligate the 5-prime small-RNA adapter", StepKind.LIGATION),
        ("rt", "Reverse-transcribe hydrolysis fragments", StepKind.REVERSE_TRANSCRIPTION),
    ),
    "hrf-seq": (
        ("cleave", "Expose RNA to hydroxyl radicals", StepKind.FRAGMENTATION),
        ("rt", "Random-prime reverse transcription to radical-induced breaks", StepKind.REVERSE_TRANSCRIPTION),
        ("ligate", "Ligate the molecular-barcode adapter to cDNA", StepKind.LIGATION),
    ),
    "shape-seq": (
        ("modify", "Modify flexible riboses with a SHAPE reagent", StepKind.CHEMICAL_CONVERSION),
        ("rt", "Primer-extend to SHAPE adducts", StepKind.REVERSE_TRANSCRIPTION),
        ("ligate", "Ligate the second sequencing adapter to cDNA", StepKind.LIGATION),
    ),
    "glori": (
        ("protect_g", "Protect guanosines with glyoxal", StepKind.CHEMICAL_CONVERSION),
        ("deaminate", "Deaminate unmethylated adenosine with nitrite", StepKind.CHEMICAL_CONVERSION),
        ("deprotect", "Heat-deprotect guanosines", StepKind.THERMAL_TREATMENT),
        ("repair", "Repair converted RNA ends", StepKind.PHOSPHORYLATION),
        ("ligate", "Ligate the modified-eCLIP 3-prime RNA adapter", StepKind.LIGATION),
        ("rt", "Reverse-transcribe converted RNA", StepKind.REVERSE_TRANSCRIPTION),
        ("ligate_cdna", "Ligate the UMI-bearing DNA adapter to cDNA", StepKind.LIGATION),
    ),
    "alkaniline-seq": (
        ("hydrolyse", "Partially hydrolyse RNA under alkaline heat", StepKind.CHEMICAL_CONVERSION),
        ("dephosphorylate", "Remove pre-existing terminal phosphates", StepKind.DEPHOSPHORYLATION),
        ("aniline", "Cleave abasic sites with aniline", StepKind.CHEMICAL_CONVERSION),
        ("ligate", "Ligate small-RNA adapters to newly phosphorylated fragments", StepKind.LIGATION),
        ("rt", "Reverse-transcribe cleavage-selected RNA", StepKind.REVERSE_TRANSCRIPTION),
    ),
}


def _linear_graph(key, spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    start = graph.input("input", spec.TITLE + " source RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("source RNA", (role("input RNA"),)), primary=True)
    reviewed = ["input"]
    current = start
    for identifier, label, kind in _LINEAR_RECIPES[key]:
        current = _apply_linear_step(graph, current, identifier, label, kind, reviewed)
    _finish_pcr(
        graph, "library", current, spec.FINAL_LIBRARY,
        spec.SEQUENCING_ENDING or spec.FINAL_CAPTION, reviewed, label=spec.TITLE,
        geometry=getattr(spec, "PCR_GEOMETRY", None),
        primers=getattr(spec, "SEQ_PRIMERS", ()),
        acquired_roles=getattr(spec, "ACQUIRED_ROLES", ()),
        read_lengths=getattr(spec, "READ_LENGTHS", None),
        read_uses=getattr(spec, "READ_USES", None))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _pars_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    rna = graph.input("input", "refolded poly(A)-selected transcriptome",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("refolded transcriptome", (role("RNA secondary structures"),)),
                      primary=True)
    v1 = graph.product("v1_branch", "RNase-V1 double-strand branch",
                       MaterialKind.MOLECULAR, state=rna.state)
    s1 = graph.product("s1_branch", "RNase-S1 single-strand branch",
                       MaterialKind.MOLECULAR, state=rna.state)
    graph.react("split", "Split refolded RNA into matched V1 and S1 aliquots",
                StepKind.SPLIT, ((rna, InputRole.TEMPLATE),), (v1, s1),
                changes_molecule=False)
    reviewed = ["input", "split"]
    for prefix, template, label in (("v1", v1, "RNase V1"), ("s1", s1, "RNase S1")):
        digested = _apply_linear_step(graph, template, prefix + "_digest",
                                      "Digest with " + label, StepKind.DIGESTION, reviewed)
        first = _apply_linear_step(graph, digested, prefix + "_ligate_5",
                                   "Select nuclease-created 5-prime phosphates by adapter ligation",
                                   StepKind.LIGATION, reviewed)
        dephos = _apply_linear_step(graph, first, prefix + "_dephosphorylate",
                                    "Dephosphorylate remaining RNA ends",
                                    StepKind.DEPHOSPHORYLATION, reviewed)
        second = _apply_linear_step(graph, dephos, prefix + "_ligate_3",
                                    "Ligate the second adapter", StepKind.LIGATION, reviewed)
        cdna = _apply_linear_step(graph, second, prefix + "_rt",
                                  "Reverse-transcribe the " + label + " branch",
                                  StepKind.REVERSE_TRANSCRIPTION, reviewed)
        _finish_pcr(graph, prefix, cdna, spec.FINAL_LIBRARY,
                    label + " branch: " + spec.FINAL_CAPTION, reviewed,
                    label="PARS " + label)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _fragseq_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    rna = graph.input("input", "purified refolded nuclear RNA", MaterialKind.BIOLOGICAL,
                      state=molecule("refolded nuclear RNA", (role("native folded RNA"),)),
                      primary=True)
    treated = graph.product("treated_branch", "nuclease-P1-treated aliquot",
                            MaterialKind.MOLECULAR, state=rna.state)
    control = graph.product("control_branch", "matched undigested control aliquot",
                            MaterialKind.MOLECULAR, state=rna.state)
    graph.react("split", "Split treated and undigested-control aliquots",
                StepKind.SPLIT, ((rna, InputRole.TEMPLATE),), (treated, control),
                changes_molecule=False)
    reviewed = ["input", "split"]
    treated = _apply_linear_step(graph, treated, "p1_digest",
                                 "Partially digest accessible RNA with nuclease P1",
                                 StepKind.DIGESTION, reviewed)
    for prefix, template, label in (("treated", treated, "P1-treated"),
                                    ("control", control, "undigested control")):
        ligated = _apply_linear_step(graph, template, prefix + "_ligate",
                                     "Ligate small-RNA sequencing adapters to the " + label + " branch",
                                     StepKind.LIGATION, reviewed)
        cdna = _apply_linear_step(graph, ligated, prefix + "_rt",
                                  "Reverse-transcribe the " + label + " branch",
                                  StepKind.REVERSE_TRANSCRIPTION, reviewed)
        _finish_pcr(graph, prefix, cdna, spec.FINAL_LIBRARY,
                    label + ": " + spec.FINAL_CAPTION, reviewed,
                    label="FragSeq " + label)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _pip_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    rnps = graph.input("input", "crosslinked cellular ribonucleoproteins",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("crosslinked RNPs", (role("RNA"), role("bound proteins"))),
                       primary=True)
    branches = []
    for identifier, label in (
            ("bound_ss", "protein-bound single-strand-RNase branch"),
            ("bound_ds", "protein-bound double-strand-RNase branch"),
            ("removed_ss", "protein-removed single-strand-RNase control"),
            ("removed_ds", "protein-removed double-strand-RNase control")):
        branches.append(graph.product(identifier, label, MaterialKind.MOLECULAR,
                                      state=rnps.state))
    graph.react("split", "Split into the four protein-state by nuclease branches",
                StepKind.SPLIT, ((rnps, InputRole.TEMPLATE),), tuple(branches),
                changes_molecule=False)
    reviewed = ["input", "split"]
    for branch in branches:
        prefix = branch.id
        current = branch
        if prefix.startswith("removed"):
            current = _apply_linear_step(graph, current, prefix + "_deproteinize",
                                         "Remove protein before nuclease digestion",
                                         StepKind.DEPROTEINIZATION, reviewed)
        current = _apply_linear_step(graph, current, prefix + "_digest",
                                     "Digest with the branch-specific RNA nuclease",
                                     StepKind.DIGESTION, reviewed)
        if prefix.startswith("bound"):
            current = _apply_linear_step(graph, current, prefix + "_deproteinize",
                                         "Reverse crosslinks after footprint digestion",
                                         StepKind.DEPROTEINIZATION, reviewed)
        current = _apply_linear_step(graph, current, prefix + "_rt",
                                     "Reverse-transcribe protected RNA fragments",
                                     StepKind.REVERSE_TRANSCRIPTION, reviewed)
        _finish_pcr(graph, prefix, current, spec.FINAL_LIBRARY,
                    branch.label + ": " + spec.FINAL_CAPTION, reviewed,
                    label="PIP-seq " + branch.label)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _attach_graph(key, spec):
    if key in _LINEAR_RECIPES:
        spec.protocol_graph = lambda: _linear_graph(key, spec)
    elif key == "pars":
        spec.protocol_graph = lambda: _pars_graph(spec)
    elif key == "fragseq":
        spec.protocol_graph = lambda: _fragseq_graph(spec)
    elif key == "protein-pip-seq":
        spec.protocol_graph = lambda: _pip_graph(spec)
    return spec


def _spec(key, **fields):
    return _attach_graph(key, SimpleNamespace(**fields))


def get(key: str):
    if key == "ottr-seq":
        return _spec(key,
            TITLE="OTTR-seq — Ordered Two-Template Relay", NOTES="01_ottr_seq.html",
            SOURCE='Upton et al., PNAS (2021); <a href="https://doi.org/10.1073/pnas.2107900118">doi:10.1073/pnas.2107900118</a>.',
            SUMMARY="BoMoC standardizes RNA 3′ ends with one ddR, then performs two ordered template jumps to copy the insert and the second adapter in one reverse-transcription reaction.", CAVEAT="",
            oligos=lambda: ("universal +1T DNA primer — 5′-/5Cy5/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT-3′", "blocked RNA complement — 5′-rGrArUrCrGrGrArArGrArGrCrArCrArCrGrUrCrUrGrArArCrUrCrCrArGrU/3SpC3/-3′", "universal adapter template — 5′-/5Cy5/ACACTCTTTCCCTACACGACGCTCTTCCrGrArUrCrUrC-3′"),
            workflow=_ottr, FINAL_LIBRARY=Construct((role("P5"), role("i5", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "random")), role("R1 adapter"), role("input cDNA"), role("R2 adapter"), role("i7", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "random")), role("P7")), name="OTTR-seq library"),
            FINAL_CAPTION="Universal-adapter branch after indexed PCR; the paper also validates a full-length, PCR-free branch.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The paper sequenced OTTR libraries on an Illumina MiniSeq 75-cycle kit.", SEQUENCING_UNAVAILABLE="The paper specifies Illumina adapter families but does not publish a distinct custom sequencing primer.")
    if key == "pseudo-seq":
        return _spec(key,
            TITLE="Pseudo-seq — transcriptome-wide pseudouridine mapping", NOTES="01_pseudo_seq.html",
            SOURCE='Carlile et al., Nature (2014); <a href="https://doi.org/10.1038/nature13802">doi:10.1038/nature13802</a>.',
            SUMMARY="Persistent Ψ–CMC adducts stop reverse transcriptase; enriched cDNA termini relative to a mock-treated control locate pseudouridines.", CAVEAT="",
            oligos=lambda: ("3′ adapter — 5′-/5Phos/TGGAATTCTCGGGTGCCAAGG/3ddC/-3′", "RT primer — 5′-/5Phos/GATCGTCGGACTGTAGAACTCTGAACCTGTCGGTGGTCGCCGTATCATT/iSp18/CACTCA/iSp18/GCCTTGGCACCCGAGAATTCCA-3′", "PCR forward — 5′-AATGATACGGCGACCACCGA-3′", "indexed PCR reverse — 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′"),
            workflow=_pseudo, FINAL_LIBRARY=Construct((role("P5"), role("stop-defined RNA cDNA"), role("sample index", tag="cbc", length_bp=6, molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="Pseudo-seq library"),
            FINAL_CAPTION="PCR product from circularized, stop-defined cDNA.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining experiment sequenced the amplified library on an Illumina HiSeq 2000.", SEQUENCING_UNAVAILABLE="The paper publishes library oligos but does not identify a separate custom sequencing-primer oligo.")
    if key == "dms-seq":
        return _spec(key,
            TITLE="DMS-seq — in-vivo RNA structure probing", NOTES="01_dms_seq.html",
            SOURCE='Rouskin et al., Nature (2014); <a href="https://doi.org/10.1038/nature12894">doi:10.1038/nature12894</a>.',
            SUMMARY="DMS marks solvent-accessible A and C in living cells; reverse-transcriptase termination sites report nucleotide accessibility after adapter ligation and cDNA circularization.", CAVEAT="",
            oligos=lambda: ("microRNA cloning linker-1 — 5′-/5rApp/CTGTAGGCACCATCAAT/3ddC/-3′",),
            workflow=_dms_seq, FINAL_LIBRARY=Construct((role("P5"), role("DMS-stop cDNA"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="DMS-seq library"),
            FINAL_CAPTION="Role-level PCR product; the paper publishes linker-1 but not the adapter PCR or sequencing-primer oligos.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining experiment sequenced RT-termination libraries on an Illumina HiSeq 2000.", SEQUENCING_UNAVAILABLE="Exact PCR and sequencing-primer oligos are not reproduced in the defining paper.")
    if key == "structure-seq2":
        return _spec(key,
            TITLE="Structure-seq2 — low-bias in-vivo RNA structure probing", NOTES="01_structure_seq2.html",
            SOURCE='Ritchey et al., Nucleic Acids Research (2017); <a href="https://doi.org/10.1093/nar/gkx533">doi:10.1093/nar/gkx533</a>.',
            SUMMARY="Random-primed DMS-stop cDNA receives a degenerate hairpin adapter by T4 DNA ligase; PAGE or biotin–streptavidin purification removes the primer–adapter by-product.", CAVEAT="",
            oligos=lambda: ("RT primer — 5′-CAGACGTGTGCTCTTCCGATCNNNNNN-3′", "hairpin donor — 5′-/5Phos/TGAAGAGCCTAGTCGCTGTTCANNNNNNCTGCCCATAGAG/3SpC3/-3′", "TruSeq forward — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGAACAGCGACTAGGCTCTTCA-3′", "indexed reverse — 5′-CAAGCAGAAGACGGCATACGAGAT[6–8-nt barcode]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-3′", "custom sequencing primer — 5′-TCTTCCGATCTTGAACAGCGACTAGGCTCTTCA-3′"),
            workflow=_structure_seq2, FINAL_LIBRARY=SS2_FINAL_LIBRARY,
            PCR_GEOMETRY=SS2_PCR_GEOMETRY,
            FINAL_CAPTION="Representative six-base-index PCR product. The custom primer abuts the unique insert, bypassing the 22-nt low-diversity prefix.",
            SEQ_PRIMERS=SS2_SEQ_PRIMERS,
            ACQUIRED_ROLES=("Read 1", "Index 1 (i7)"),
            READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 6},
            READ_USES={"Read 1": "include", "Index 1 (i7)": "demultiplex"},
            SEQUENCING_INTRO="The custom primer and standard multiplexing-index primer both land on source-exact PCR-derived sites.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
    if key == "ribomethseq":
        return _spec(key,
            TITLE="RiboMethSeq — 2′-O-methylation by protected cleavage", NOTES="01_ribomethseq.html",
            SOURCE='Birkedal et al., Angewandte Chemie (2015), <a href="https://doi.org/10.1002/anie.201408362">doi:10.1002/anie.201408362</a>; Illumina implementation: Marchand et al., Nucleic Acids Research (2016), <a href="https://doi.org/10.1093/nar/gkw547">doi:10.1093/nar/gkw547</a>.',
            SUMMARY="Partial alkaline hydrolysis avoids phosphodiester bonds immediately 3′ of 2′-O-methylated riboses; depletion of mapped fragment ends identifies and quantifies Nm sites.", CAVEAT="",
            oligos=lambda: ("NEBNext 3′ small-RNA adapter — 5′-/rApp/AGATCGGAAGAGCACACGTCT/3AmMO/-3′",),
            workflow=_ribomethseq, FINAL_LIBRARY=Construct((role("P5"), role("hydrolysis-end cDNA"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="RiboMethSeq library"),
            FINAL_CAPTION="Standard small-RNA library retaining the exact alkaline-cleavage endpoints.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The Illumina implementation used paired-end MiSeq or single-read 50-cycle HiSeq 1000 sequencing; standard kit primers enter from the ligated adapters.", SEQUENCING_UNAVAILABLE="The method paper names the commercial NEBNext Small RNA kit but does not reproduce a distinct sequencing-primer oligo.")
    if key == "pars":
        return _spec(key,
            TITLE="PARS — parallel analysis of RNA structure", NOTES="01_pars.html",
            SOURCE='Kertesz et al., Nature (2010), <a href="https://doi.org/10.1038/nature09322">doi:10.1038/nature09322</a>; detailed protocol: Wan et al., Nature Protocols (2013), <a href="https://doi.org/10.1038/nprot.2013.045">doi:10.1038/nprot.2013.045</a>.',
            SUMMARY="Matched RNase V1 and S1 digests expose double- and single-stranded positions; selective capture of nuclease-created 5′ phosphates turns branch-specific read starts into a per-base structure score.", CAVEAT="",
            oligos=lambda: (
                "SOLiD Small RNA kit common P1 / 5′ PCR primer — 5′-CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT-3′",
                "SOLiD Small RNA kit indexed 3′ PCR primer 1 (kit example; experiment's primer set unstated) — 5′-CTGCCCCGGGTTCCTCATTCTCTAAGCCCCTGCTGTACGGCCAAGGCG-3′",
            ), workflow=_pars,
            FINAL_LIBRARY=Construct((role("SOLiD P1 adapter"), role("V1- or S1-cleavage-defined cDNA"), role("SOLiD P2 adapter")), name="PARS library"),
            FINAL_CAPTION="One library is made for each nuclease branch; the defining experiment used the SOLiD platform.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining libraries were sequenced on Applied Biosystems SOLiD; the vendor guide publishes library PCR primers but not the complete adaptor mix or sequencing primer.", SEQUENCING_UNAVAILABLE="Adaptor Mix A remains a degenerate kit oligo set without complete published bases; the exact sequencing primer is also unresolved.")
    if key == "fragseq":
        return _spec(key,
            TITLE="FragSeq — nuclease accessibility profiling of RNA", NOTES="01_fragseq.html",
            SOURCE='Underwood et al., Nature Methods (2010), <a href="https://doi.org/10.1038/nmeth.1529">doi:10.1038/nmeth.1529</a>.',
            SUMMARY="Partial nuclease P1 digestion converts exposed single-stranded RNA into mapped fragment termini, directly profiling RNA accessibility across many transcripts.", CAVEAT="FragSeq probes purified, refolded RNA. It reports nuclease accessibility, not native intracellular accessibility or abundance alone.",
            oligos=lambda: (), workflow=_fragseq,
            FINAL_LIBRARY=Construct((role("SOLiD P1 adapter"), role("P1-cleavage-defined cDNA"), role("SOLiD P2 adapter")), name="FragSeq library"),
            FINAL_CAPTION="The informative coordinate is a nuclease-P1-defined RNA fragment end.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining experiment converted treated and control RNA fragments with the SOLiD small-RNA workflow and sequenced their ends.", SEQUENCING_UNAVAILABLE="Exact proprietary SOLiD adapter and sequencing-primer sequences are not represented.")
    if key == "hrf-seq":
        return _spec(key,
            TITLE="HRF-Seq — hydroxyl-radical RNA accessibility", NOTES="01_hrf_seq.html",
            SOURCE='Kielpinski and Vinther, Nucleic Acids Research (2014), <a href="https://doi.org/10.1093/nar/gku167">doi:10.1093/nar/gku167</a>.',
            SUMMARY="Hydroxyl radicals cleave solvent-exposed ribose backbones; random-primed RT, a seven-base molecular barcode and paired-end sequencing recover both break and priming coordinates.", CAVEAT="The defining demonstration probes purified RNAs. Hydroxyl-radical accessibility reflects backbone solvent exposure and tertiary packing, not only base pairing.",
            oligos=lambda: (
                "RT_random_primer — 5′-AGACGTGTGCTCTTCCGATCTNNNNNNNNS-3′",
                "RT_structure_cassette — 5′-AGACGTGTGCTCTTCCGATCTGAACCGGACCGAAGCCCG-3′",
                "LIGATION_ADAPTER_RB — 5′-PHO-NNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-3NHC3-3′",
                "PCR_forward — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT-3′",
                *tuple(f"PCR reverse {index} — 5′-{_hrf_reverse(index)}-3′" for index in HRF_INDEXES)),
            workflow=_hrf_seq, FINAL_LIBRARY=HRF_FINAL_LIBRARY,
            PCR_GEOMETRY=HRF_PCR_GEOMETRY,
            FINAL_CAPTION="Representative indexed member of the published four-primer family. Read 1 begins with the seven-base molecular barcode; Read 2 reports the random-priming end.",
            SEQ_PRIMERS=HRF_SEQ_PRIMERS,
            ACQUIRED_ROLES=("Read 1", "Index 1 (i7)", "Read 2"),
            READ_LENGTHS={"Read 1": 100, "Index 1 (i7)": 6, "Read 2": 100},
            READ_USES={"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex", "Read 2": "include"},
            SEQUENCING_INTRO="The source uses standard paired-end Illumina primers. Their sites, the first-read molecular barcode, sample index and second-read random-primer coordinate are derived from the published oligos.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
    if key == "shape-seq":
        return _spec(key,
            TITLE="SHAPE-Seq — multiplexed RNA flexibility probing", NOTES="01_shape_seq.html",
            SOURCE='Lucks et al., PNAS (2011), <a href="https://doi.org/10.1073/pnas.1106501108">doi:10.1073/pnas.1106501108</a>.',
            SUMMARY="SHAPE chemistry acylates flexible riboses, reverse transcriptase stops at the adducts, and paired-end sequencing links each stop to a pre-attached RNA identity barcode.", CAVEAT="The original SHAPE-Seq protocol assays designed RNA pools in vitro; later in-cell and transcriptome-wide SHAPE methods use different capture architectures.",
            oligos=lambda: (
                f"A_adapter_t / PE_forward — 5′-{SHAPE_A_T}-3′",
                f"A_adapter_b — 5′-/5Phos/{SHAPE_A_B}/3SpC3/-3′",
                "treated RT primer pattern — 5′-[A_adapter_t]-RRRY-[target-specific RT primer]-3′",
                "control RT primer pattern — 5′-[A_adapter_t]-YYYR-[target-specific RT primer]-3′"),
            workflow=_shape_seq, FINAL_LIBRARY=SHAPE_FINAL_LIBRARY,
            FINAL_CAPTION="Source-exact v1 adapter arms flank a target-specific RT product. The first paired read begins with the four-base treatment handle; the opposite read locates the SHAPE-dependent RT stop.",
            SEQ_PRIMERS=(), READ_LENGTHS={"Read 1": 50, "Read 2": 50},
            SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining method used 2×50-cycle paired-end sequencing on a Genome Analyzer IIx.",
            SEQUENCING_UNAVAILABLE="The original supplement prints both adapter molecules and PCR-primer construction, but does not name a separate immutable oligo for the historical Genome Analyzer IIx paired-end sequencing-primer reagent.")
    if key == "protein-pip-seq":
        return _spec(key,
            TITLE="PIP-seq — protein interaction profile sequencing", NOTES="01_protein_pip_seq.html",
            SOURCE='Silverman et al., Genome Biology (2014), <a href="https://doi.org/10.1186/gb-2014-15-1-r3">doi:10.1186/gb-2014-15-1-r3</a>.',
            SUMMARY="Four matched libraries combine single- and double-strand-specific RNases with protein-present and protein-removed branches, revealing both RNA accessibility and protein-protected sites.", CAVEAT="This is protein-interaction-profile PIP-seq, not the unrelated particle-templated single-cell RNA-seq protocol with the same acronym.",
            oligos=lambda: (), workflow=_protein_pip_seq,
            FINAL_LIBRARY=Construct((role("P5"), role("RNase-protected strand-specific cDNA"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="protein-interaction-profile PIP-seq library"),
            FINAL_CAPTION="One strand-specific sequencing library is made for each protein-state × nuclease-specificity branch.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining study sequenced strand-specific libraries from all four branches; branch identity must be retained through sample indexing.", SEQUENCING_UNAVAILABLE="Exact commercial library-kit primer sequences are not represented in the defining schematic.")
    if key == "glori":
        return _spec(key,
            TITLE="GLORI — absolute m6A quantification by A-to-I conversion", NOTES="01_glori.html",
            SOURCE='Liu et al., Nature Biotechnology (2023), <a href="https://doi.org/10.1038/s41587-022-01487-9">doi:10.1038/s41587-022-01487-9</a>; detailed protocol: Shen et al., Nature Protocols (2024), <a href="https://doi.org/10.1038/s41596-023-00937-1">doi:10.1038/s41596-023-00937-1</a>.',
            SUMMARY="Glyoxal-assisted nitrite treatment converts unmethylated adenosine to inosine while m6A remains adenosine; A retention therefore reports absolute methylation stoichiometry.", CAVEAT="This page follows the modified-eCLIP library branch. Exact ligation adapters and RT primer come from the authors' 2024 detailed protocol; equivalence to every 2023 experimental lot is not established. The authors also validated NEBNext Small RNA and SMARTer stranded kit alternatives.",
            oligos=lambda: (
                "2024 detailed protocol, 3′ RNA adaptor — 5′-rAppAGATCGGAAGAGCGTCGTG-3′",
                "2024 detailed protocol, RT primer — 5′-ACACGACGCTCTTCCGATCT-3′",
                "2024 detailed protocol, 3′ DNA adaptor — 5′-Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG-3′"),
            workflow=_glori,
            FINAL_LIBRARY=Construct((role("P5"), role("10-nt UMI", tag="umi", length_bp=10, molecular_feature=feature("umi", "umi", "random")), role("A-to-I-converted cDNA"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="GLORI modified-eCLIP library"),
            FINAL_CAPTION="Modified-eCLIP branch. The 10-nt UMI is inline; converted A is read as G while m6A is retained as A.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="Defining modified-eCLIP libraries were sequenced paired-end on an Illumina HiSeq X Ten; reads begin with the 10-nt UMI before converted insert sequence.", SEQUENCING_UNAVAILABLE="The 2024 protocol's Supplementary Table 1 gives both ligation adapters and the RT primer. Complete library amplification and sequencing-primer oligos remain unreported for this branch; its final primer landing sites are not reconstructed.")
    if key == "alkaniline-seq":
        return _spec(key,
            TITLE="AlkAniline-Seq — positive selection of modification-derived RNA cleavage", NOTES="01_alkaniline_seq.html",
            SOURCE='Marchand et al., Angewandte Chemie (2018), <a href="https://doi.org/10.1002/anie.201810946">doi:10.1002/anie.201810946</a>.',
            SUMMARY="Alkali-labile m7G, m3C and related residues form abasic sites; after background dephosphorylation, aniline cleavage creates the only 5′ phosphates accepted by the first sequencing adapter.", CAVEAT="AlkAniline-Seq also responds to dihydrouridine and some other abasic-site-forming lesions. A signal identifies the chemistry, not uniquely one modification class.",
            oligos=lambda: (
                "NEBNext E7300/E7330 3′ SR adapter — 5′-rAppAGATCGGAAGAGCACACGTCT-NH2-3′",
                "NEBNext E7300/E7330 5′ SR adapter — 5′-rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC-3′",
                "NEBNext E7300/E7330 SR RT primer — 5′-AGACGTGTGCTCTTCCGATCT-3′",
                "NEBNext E7300/E7330 SR PCR primer — 5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG*A-3′",
                "NEBNext E7300/E7330 index 1 PCR primer (kit example; experiment's index choice unstated) — 5′-CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3′",
            ),
            workflow=_alkaniline,
            FINAL_LIBRARY=Construct((role("P5"), role("cDNA beginning N+1 after aniline-sensitive modification"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="AlkAniline-Seq library"),
            FINAL_CAPTION="The informative coordinate is the cDNA/read start one nucleotide downstream of the modified base.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining workflow used single-read 50-nt Illumina sequencing; a standard small-RNA read primer enters from the selectively ligated 5′ adapter.", SEQUENCING_UNAVAILABLE="NEB confirms standard Illumina sequencing primers for this kit family; the role-level final molecule has not yet reconstructed their exact landing site.")
    raise KeyError(key)
