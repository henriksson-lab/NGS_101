"""Reusable constructions for RNA end-capture and modification-stop libraries."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import Construct, Segment, Workflow, feature
from protocol_model import chain, molecule, role


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
SS2_FORWARD = Segment("Structure-seq2 TruSeq forward", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGAACAGCGACTAGGCTCTTCA", "p5")
SS2_REVERSE = Segment("Structure-seq2 indexed reverse", "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "p7", placeholder=True,
                      feature=feature("sample_index", "sample_index", "unknown"))
SS2_CUSTOM = Segment("Structure-seq2 custom sequencing primer", "TCTTCCGATCTTGAACAGCGACTAGGCTCTTCA", "r1")


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
        ("Ligate modified-eCLIP adapters carrying a 10-nt UMI",
         molecule("UMI-tagged GLORI RNA", (role("5-prime adapter"), role("10-nt UMI"), role("converted RNA"), role("3-prime adapter"))),
         "The UMI is removed computationally after PCR-duplicate handling."),
        ("Reverse-transcribe, index-amplify and paired-end sequence",
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
         molecule("selected AlkAniline RNA", (role("NEBNext 5-prime adapter"), role("RNA beginning N+1 after modification"), role("NEBNext 3-prime adapter"))),
         "The 5′ adapter positively selects the aniline-generated phosphate; the read start is N+1."),
        ("Reverse-transcribe, barcode by PCR and single-end sequence",
         molecule("AlkAniline-Seq library", (role("P5"), role("cleavage-defined cDNA"), role("sample index"), role("P7"))),
         "Only read 5′ ends are required for normalized-cleavage scoring."),
    ))


def get(key: str):
    if key == "ottr-seq":
        return SimpleNamespace(
            TITLE="OTTR-seq — Ordered Two-Template Relay", NOTES="01_ottr_seq.html",
            SOURCE='Upton et al., PNAS (2021); <a href="https://doi.org/10.1073/pnas.2107900118">doi:10.1073/pnas.2107900118</a>.',
            SUMMARY="BoMoC standardizes RNA 3′ ends with one ddR, then performs two ordered template jumps to copy the insert and the second adapter in one reverse-transcription reaction.", CAVEAT="",
            oligos=lambda: ("universal +1T DNA primer — 5′-/5Cy5/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT-3′", "blocked RNA complement — 5′-rGrArUrCrGrGrArArGrArGrCrArCrArCrGrUrCrUrGrArArCrUrCrCrArGrU/3SpC3/-3′", "universal adapter template — 5′-/5Cy5/ACACTCTTTCCCTACACGACGCTCTTCCrGrArUrCrUrC-3′"),
            workflow=_ottr, FINAL_LIBRARY=Construct((role("P5"), role("i5", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "random")), role("R1 adapter"), role("input cDNA"), role("R2 adapter"), role("i7", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "random")), role("P7")), name="OTTR-seq library"),
            FINAL_CAPTION="Universal-adapter branch after indexed PCR; the paper also validates a full-length, PCR-free branch.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The paper sequenced OTTR libraries on an Illumina MiniSeq 75-cycle kit.", SEQUENCING_UNAVAILABLE="The paper specifies Illumina adapter families but does not publish a distinct custom sequencing primer.")
    if key == "pseudo-seq":
        return SimpleNamespace(
            TITLE="Pseudo-seq — transcriptome-wide pseudouridine mapping", NOTES="01_pseudo_seq.html",
            SOURCE='Carlile et al., Nature (2014); <a href="https://doi.org/10.1038/nature13802">doi:10.1038/nature13802</a>.',
            SUMMARY="Persistent Ψ–CMC adducts stop reverse transcriptase; enriched cDNA termini relative to a mock-treated control locate pseudouridines.", CAVEAT="",
            oligos=lambda: ("3′ adapter — 5′-/5Phos/TGGAATTCTCGGGTGCCAAGG/3ddC/-3′", "RT primer — 5′-/5Phos/GATCGTCGGACTGTAGAACTCTGAACCTGTCGGTGGTCGCCGTATCATT/iSp18/CACTCA/iSp18/GCCTTGGCACCCGAGAATTCCA-3′", "PCR forward — 5′-AATGATACGGCGACCACCGA-3′", "indexed PCR reverse — 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′"),
            workflow=_pseudo, FINAL_LIBRARY=Construct((role("P5"), role("stop-defined RNA cDNA"), role("sample index", tag="cbc", length_bp=6, molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="Pseudo-seq library"),
            FINAL_CAPTION="PCR product from circularized, stop-defined cDNA.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining experiment sequenced the amplified library on an Illumina HiSeq 2000.", SEQUENCING_UNAVAILABLE="The paper publishes library oligos but does not identify a separate custom sequencing-primer oligo.")
    if key == "dms-seq":
        return SimpleNamespace(
            TITLE="DMS-seq — in-vivo RNA structure probing", NOTES="01_dms_seq.html",
            SOURCE='Rouskin et al., Nature (2014); <a href="https://doi.org/10.1038/nature12894">doi:10.1038/nature12894</a>.',
            SUMMARY="DMS marks solvent-accessible A and C in living cells; reverse-transcriptase termination sites report nucleotide accessibility after adapter ligation and cDNA circularization.", CAVEAT="",
            oligos=lambda: ("microRNA cloning linker-1 — 5′-/5rApp/CTGTAGGCACCATCAAT/3ddC/-3′",),
            workflow=_dms_seq, FINAL_LIBRARY=Construct((role("P5"), role("DMS-stop cDNA"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="DMS-seq library"),
            FINAL_CAPTION="Role-level PCR product; the paper publishes linker-1 but not the adapter PCR or sequencing-primer oligos.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining experiment sequenced RT-termination libraries on an Illumina HiSeq 2000.", SEQUENCING_UNAVAILABLE="Exact PCR and sequencing-primer oligos are not reproduced in the defining paper.")
    if key == "structure-seq2":
        return SimpleNamespace(
            TITLE="Structure-seq2 — low-bias in-vivo RNA structure probing", NOTES="01_structure_seq2.html",
            SOURCE='Ritchey et al., Nucleic Acids Research (2017); <a href="https://doi.org/10.1093/nar/gkx533">doi:10.1093/nar/gkx533</a>.',
            SUMMARY="Random-primed DMS-stop cDNA receives a degenerate hairpin adapter by T4 DNA ligase; PAGE or biotin–streptavidin purification removes the primer–adapter by-product.", CAVEAT="",
            oligos=lambda: ("hairpin donor — 5′-/5Phos/TGAAGAGCCTAGTCGCTGTTCANNNNNNCTGCCCATAGAG/3SpC3/-3′", "TruSeq forward — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGAACAGCGACTAGGCTCTTCA-3′", "indexed reverse — 5′-CAAGCAGAAGACGGCATACGAGAT[6–8-nt barcode]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-3′", "custom sequencing primer — 5′-TCTTCCGATCTTGAACAGCGACTAGGCTCTTCA-3′"),
            workflow=_structure_seq2, FINAL_LIBRARY=Construct((role("P5"), SS2_CUSTOM, role("stop-defined cDNA"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="Structure-seq2 library"),
            FINAL_CAPTION="The custom primer abuts the unique insert, bypassing the 22-nt low-diversity prefix.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The custom primer 5′-TCTTCCGATCTTGAACAGCGACTAGGCTCTTCA-3′ reads directly into the stop-defined cDNA in a 150-cycle MiSeq run.", SEQUENCING_UNAVAILABLE="")
    if key == "ribomethseq":
        return SimpleNamespace(
            TITLE="RiboMethSeq — 2′-O-methylation by protected cleavage", NOTES="01_ribomethseq.html",
            SOURCE='Birkedal et al., Angewandte Chemie (2015), <a href="https://doi.org/10.1002/anie.201408362">doi:10.1002/anie.201408362</a>; Illumina implementation: Marchand et al., Nucleic Acids Research (2016), <a href="https://doi.org/10.1093/nar/gkw547">doi:10.1093/nar/gkw547</a>.',
            SUMMARY="Partial alkaline hydrolysis avoids phosphodiester bonds immediately 3′ of 2′-O-methylated riboses; depletion of mapped fragment ends identifies and quantifies Nm sites.", CAVEAT="",
            oligos=lambda: ("NEBNext 3′ small-RNA adapter — 5′-/rApp/AGATCGGAAGAGCACACGTCT/3AmMO/-3′",),
            workflow=_ribomethseq, FINAL_LIBRARY=Construct((role("P5"), role("hydrolysis-end cDNA"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="RiboMethSeq library"),
            FINAL_CAPTION="Standard small-RNA library retaining the exact alkaline-cleavage endpoints.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The Illumina implementation used paired-end MiSeq or single-read 50-cycle HiSeq 1000 sequencing; standard kit primers enter from the ligated adapters.", SEQUENCING_UNAVAILABLE="The method paper names the commercial NEBNext Small RNA kit but does not reproduce a distinct sequencing-primer oligo.")
    if key == "pars":
        return SimpleNamespace(
            TITLE="PARS — parallel analysis of RNA structure", NOTES="01_pars.html",
            SOURCE='Kertesz et al., Nature (2010), <a href="https://doi.org/10.1038/nature09322">doi:10.1038/nature09322</a>; detailed protocol: Wan et al., Nature Protocols (2013), <a href="https://doi.org/10.1038/nprot.2013.045">doi:10.1038/nprot.2013.045</a>.',
            SUMMARY="Matched RNase V1 and S1 digests expose double- and single-stranded positions; selective capture of nuclease-created 5′ phosphates turns branch-specific read starts into a per-base structure score.", CAVEAT="",
            oligos=lambda: (), workflow=_pars,
            FINAL_LIBRARY=Construct((role("SOLiD P1 adapter"), role("V1- or S1-cleavage-defined cDNA"), role("SOLiD P2 adapter")), name="PARS library"),
            FINAL_CAPTION="One library is made for each nuclease branch; the defining experiment used the SOLiD platform.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining libraries were sequenced on Applied Biosystems SOLiD; their ligation-sequencing primer oligos are not printed in the defining paper.", SEQUENCING_UNAVAILABLE="The public source describes the SOLiD adapter roles but does not reproduce exact adapter or sequencing-primer sequences.")
    if key == "glori":
        return SimpleNamespace(
            TITLE="GLORI — absolute m6A quantification by A-to-I conversion", NOTES="01_glori.html",
            SOURCE='Liu et al., Nature Biotechnology (2023), <a href="https://doi.org/10.1038/s41587-022-01487-9">doi:10.1038/s41587-022-01487-9</a>; detailed protocol: Shen et al., Nature Protocols (2024), <a href="https://doi.org/10.1038/s41596-023-00937-1">doi:10.1038/s41596-023-00937-1</a>.',
            SUMMARY="Glyoxal-assisted nitrite treatment converts unmethylated adenosine to inosine while m6A remains adenosine; A retention therefore reports absolute methylation stoichiometry.", CAVEAT="This page follows the defining modified-eCLIP library branch. The authors also validated NEBNext Small RNA and SMARTer stranded library kits after the same GLORI conversion.",
            oligos=lambda: ("3′ adapter sequence used by the ligation branch — 5′-/5rApp/AGATCGGAAGAGCACACGTC/3AmMO/-3′", "RT primer — 5′-GACGTGTGCTCTTCCGATCT-3′"),
            workflow=_glori,
            FINAL_LIBRARY=Construct((role("P5"), role("10-nt UMI", tag="umi", length_bp=10, molecular_feature=feature("umi", "umi", "random")), role("A-to-I-converted cDNA"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="GLORI modified-eCLIP library"),
            FINAL_CAPTION="Modified-eCLIP branch. The 10-nt UMI is inline; converted A is read as G while m6A is retained as A.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="Defining modified-eCLIP libraries were sequenced paired-end on an Illumina HiSeq X Ten; reads begin with the 10-nt UMI before converted insert sequence.", SEQUENCING_UNAVAILABLE="The protocol identifies standard Illumina entry sites but does not publish a distinct custom sequencing-primer oligo.")
    if key == "alkaniline-seq":
        return SimpleNamespace(
            TITLE="AlkAniline-Seq — positive selection of modification-derived RNA cleavage", NOTES="01_alkaniline_seq.html",
            SOURCE='Marchand et al., Angewandte Chemie (2018), <a href="https://doi.org/10.1002/anie.201810946">doi:10.1002/anie.201810946</a>.',
            SUMMARY="Alkali-labile m7G, m3C and related residues form abasic sites; after background dephosphorylation, aniline cleavage creates the only 5′ phosphates accepted by the first sequencing adapter.", CAVEAT="AlkAniline-Seq also responds to dihydrouridine and some other abasic-site-forming lesions. A signal identifies the chemistry, not uniquely one modification class.",
            oligos=lambda: ("NEBNext 3′ small-RNA adapter — 5′-/rApp/AGATCGGAAGAGCACACGTCT/3AmMO/-3′",),
            workflow=_alkaniline,
            FINAL_LIBRARY=Construct((role("P5"), role("cDNA beginning N+1 after aniline-sensitive modification"), role("sample index", tag="cbc", molecular_feature=feature("sample_index", "sample_index", "unknown")), role("P7")), name="AlkAniline-Seq library"),
            FINAL_CAPTION="The informative coordinate is the cDNA/read start one nucleotide downstream of the modified base.", SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining workflow used single-read 50-nt Illumina sequencing; a standard small-RNA read primer enters from the selectively ligated 5′ adapter.", SEQUENCING_UNAVAILABLE="The paper names the commercial NEBNext Small RNA kit but does not reproduce a distinct sequencing-primer oligo.")
    raise KeyError(key)
