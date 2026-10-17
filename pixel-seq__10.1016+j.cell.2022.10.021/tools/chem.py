"""Pixel-seq polony-gel stamping and library chemistry."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, feature, inline_sequence, revcomp
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Disposition, Highlight, InputRole,
                            MaterialKind, ProtocolGraph, SequencingPlan, StepKind,
                            Terminal)
from protocol_model import chain, molecule, role


SOURCE_LOCATION = ("author manuscript PMC9691594 (JATS), STAR Methods 'Polony template "
                   "construction' to 'Library construction and sequencing' and Key Resources "
                   "Table oligonucleotides")
UMI_PRIMER = "CTACACGACGCTCTTCCGATCTNNNNNNNCCGCTAATACGACTCACTATAGGG"
TRUSEQ_SEQUENCING_PRIMER = "CCTACACGACGCTCTTCCGATC"
TSO_PRIMER = "/5’BIOSG/AAGCAGTGGTATCAACGCAGAGT"
TSO_PRIMER_BASES = "AAGCAGTGGTATCAACGCAGAGT"
TSO = "AAGCTGGTATCAACGCAGAGTGAATrG+GrG"
TRUSEQ_LIB_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
POLONY_SEQUENCING_PRIMER = "TTCCATCATTGTTCCCGGGTTCCTCATTCTC"
BA_PLUS = "AATGATACGGCGACCACCGAGAUCTACAC"
BA_MINUS = "CAAGCAGAAGACGGCATACGAGAT"
PS_BA_PLUS = "/5PS/TTTTTTTTTTAATGATACGGCGACCACCGAGAUCTACAC"
PS_BA_MINUS = "/5PS/TTTTTTTTTTCAAGCAGAAGACGGCATACGAGAT"
CY5_PROBE = "/5Cy5/CTTGCGCTAGTCCGGCAACCTCGGTGG"
NEXTERA_INDEX_PRIMER = "CAAGCAGAAGACGGCATACGAGAT[INDEX]GTCTCGTGGGCTCGG"
# Seed-template components exactly as labelled in STAR Methods "Polony template
# construction" (author-manuscript JATS text lines 186-200).
TEMPLATE_P5 = "AATGATACGGCGACCACCGAGATCTACAC"
TEMPLATE_T7 = "TGCGGCCGCTAATACGACTCACTATAGGGATCT"
TEMPLATE_BARCODE = "N" * 24
POLONY_SEQUENCING_SITE = "AGAGAATGAGGAACCCGGGAACAATGATGGAA"
POLY_T_PROBE = "T" * 19
TAQI_SITE = "TCGA"
CY5_PROBE_SITE = "ACCACCGAGGTTGCCGGACTAGCGCAAG"
INTERNAL_SEGMENT = (
    "TACTTGTCCATTCCTGAAGAAATATTATATTTATACAACTTACCCATAGAATCCTATTTACTAGGAAAGGAAAAGCCTCCTATTTATAAAAATTGGATAGAGCTTTCTCAACAACAGTGGAATATCAATGATAGAACAATTGCCGATTTATTAGATGGGGTCTTAATAATACCATCGA")
TEMPLATE_P7 = "TCTCGTATGCCGTCTTCTGCTTG"
POLONY_TEMPLATE = (TEMPLATE_P5 + TEMPLATE_T7 + TEMPLATE_BARCODE + POLONY_SEQUENCING_SITE
                   + POLY_T_PROBE + TAQI_SITE + CY5_PROBE_SITE + INTERNAL_SEGMENT
                   + TEMPLATE_P7)

# The UMI primer is TruSeq Read 1 sequence, seven N bases and a 3' anchor that
# matches the template T7-promoter segment; derive the junctions rather than retype.
UMI_PRIMER_TRUSEQ = UMI_PRIMER[:UMI_PRIMER.index("N")]
UMI_PRIMER_ANCHOR = UMI_PRIMER[UMI_PRIMER.rindex("N") + 1:]
if not il.TRUSEQ_READ1.endswith(UMI_PRIMER_TRUSEQ):
    raise AssertionError("UMI primer 5' arm is not TruSeq Read 1 sequence")
if UMI_PRIMER_ANCHOR not in TEMPLATE_T7:
    raise AssertionError("UMI primer anchor is absent from the template T7 segment")
T7_AFTER_ANCHOR = TEMPLATE_T7[TEMPLATE_T7.index(UMI_PRIMER_ANCHOR) + len(UMI_PRIMER_ANCHOR):]
if il.P5[:-4] + il.TRUSEQ_READ1 != TRUSEQ_LIB_P5:
    raise AssertionError("published TruSeq libP5 primer differs from the P5 + Read 1 arm")
if revcomp(NEXTERA_INDEX_PRIMER.split("]")[1]) != nx.S7_RC:
    raise AssertionError("Nextera index primer 3' end does not prime the S7 site")


SPATIAL = feature("pixel_seq.spatial", "spatial_barcode", "random",
                  note="24-base polony identity mapped once per replicated gel family")
UMI = feature("pixel_seq.umi", "umi", "random",
              note="seven random bases introduced by the published UMI primer")


def pixel_library():
    """Informative product of the published TruSeq libP5 + Nextera i7 index PCR.

    STAR Methods: Nextera XT tagmentation, then PCR with the TruSeq libP5 primer and
    Nextera index primers. The P5 end therefore comes from the UMI primer's TruSeq
    Read 1 arm (no i5 window) and only the Tn5 end carries a sample index.
    """
    lib = Construct([
        seg("P5", il.P5[:-4], "p5"),
        seg("TruSeq Read 1 arm (UMI primer)", il.TRUSEQ_READ1, "r1"),
        seg("7-nt UMI", "U" * 7, "umi", placeholder=True, feature=UMI),
        seg("T7-promoter anchor (UMI primer 3' end)", UMI_PRIMER_ANCHOR, "t7"),
        seg("T7-promoter end", T7_AFTER_ANCHOR, "t7"),
        seg("24-nt polony spatial barcode", "B" * 24, "cbc", placeholder=True,
            feature=SPATIAL),
        seg("polony sequencing-primer site", POLONY_SEQUENCING_SITE, "r3"),
        seg("19-nt poly(T) capture probe", POLY_T_PROBE),
        role("captured transcript cDNA (antisense)"),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 sample index reverse complement (length not printed; 8 shown)", "I" * 8,
            "cbc", placeholder=True, inferred=True,
            feature=feature("sample_index_i7", "sample_index", "unknown")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="informative Pixel-seq library")
    primers = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("Pixel-seq library: " + "; ".join(problems))
    return lib, primers


def workflow():
    return chain(molecule("crosslinked polyacrylamide stamp gel", (
        role("covalently grafted plus primer"), role("covalently grafted minus primer"))), (
        ("Seed random 24-base templates and bridge-amplify surface polonies",
         molecule("stamp gel with clonal barcode polonies", (
             role("gel-bound primer"), role("one random 24-base spatial barcode"),
             role("USER-linearizable opposite strand"))),
         "Each approximately one-micrometre polony is a clonal coordinate feature. USER treatment leaves a single-stranded stamp template."),
        ("Press stamp and copy gels together; extend copied strands across the interface",
         molecule("copied barcode pattern on a second gel", (
             role("same polony geometry"), role("complementary copied barcode strands"))),
         "Taq polymerase copies the whole array without consuming the stamp. Bridge amplification on the copy gel restores high local DNA density, so many gels inherit one barcode map."),
        ("TaqI-cut copied polonies to expose their 3-prime poly(dT) capture probes",
         molecule("RNA-capture polony array", (
             role("gel-bound spatial barcode"), role("TaqI cleavage boundary"),
             role("exposed poly(dT)"))),
         "A few sibling gels are sequenced in situ to map the 24-base polony barcodes; other gels with the replicated pattern receive tissue."),
        ("Hybridize tissue RNA and template-switch first-strand synthesis on the gel",
         molecule("spatially barcoded cDNA", (
             role("24-base polony barcode"), role("poly(dT)"),
             role("captured transcript cDNA"), role("template-switch handle"))),
         "Maxima H-minus RT copies poly(A) RNA and appends the TSO-derived universal end. The barcode is inherited from the polony beneath the RNA."),
        ("TSO-prime second strands, elute with KOH, extend the UMI primer and PCR-amplify",
         molecule("amplified barcode–UMI–transcript cDNA", (
             role("24-base spatial barcode"), role("seven-base UMI"),
             role("captured transcript cDNA"), role("universal PCR handles"))),
         "Bst 2.0 WarmStart extends the TSO primer; a single Q5 extension of the UMI primer installs the UMI before thermolabile Exonuclease I cleanup; PCR then uses the TSO and TruSeq sequencing primers."),
        ("Tagment 1 ng amplified cDNA with Nextera XT; PCR with TruSeq libP5 and a Nextera i7 index primer",
         molecule("paired-end Pixel-seq library", (
             role("P5 and TruSeq Read 1 arm"), role("UMI, T7 anchor and spatial barcode"),
             role("transcript insert"), role("mosaic end, S7, i7 and P7"))),
         "Only fragments retaining the TruSeq barcode end, a UMI and a transcript alignment contribute spatial molecules."),
    ))


def _protocol_graph(final, primers) -> ProtocolGraph:
    """Typed Pixel-seq process, including reusable stamp and uninformative fragments."""
    source = "Fu et al., Cell (2022); doi:10.1016/j.cell.2022.10.021"
    graph = ProtocolGraph("Pixel-seq", authority=source)
    stamp = graph.input(
        "stamp_gel", "crosslinked polyacrylamide stamp gel",
        MaterialKind.MOLECULAR,
        state=molecule("crosslinked polyacrylamide stamp gel", (
            role("covalently grafted plus primer"), role("covalently grafted minus primer"))),
        primary=True)
    random_templates = graph.input("random_barcode_templates",
                                   "370-bp seed templates (24 random bases), BA(±)-PCR 15 cycles, gel-purified, 8–12 pM",
                                   MaterialKind.OLIGO)
    plus_primer = graph.input("stamp_plus_primer", "gel-grafted PS-BA(+) (contains U)",
                              MaterialKind.OLIGO)
    minus_primer = graph.input("stamp_minus_primer", "gel-grafted PS-BA(−)",
                               MaterialKind.OLIGO)
    bridge_polymerase = graph.input("stamp_polymerase",
                                    "Taq first-strand extension (74 °C), then Bst bridge amplification (22 formamide cycles, 60 °C)",
                                    MaterialKind.CATALYST)
    polonies = graph.product(
        "stamp_polonies", "stamp gel with clonal barcode polonies",
        MaterialKind.MOLECULAR,
        state=molecule("stamp gel with clonal barcode polonies", (
            role("gel-bound primer"), role("one random 24-base spatial barcode"),
            role("USER-linearizable opposite strand"))))
    double_polonies = graph.product(
        "double_stranded_polonies", "double-stranded bridge-amplified polonies",
        MaterialKind.MOLECULAR,
        state=molecule("double-stranded stamp polonies", (
            role("PS-BA(+)-anchored strand with U"), role("PS-BA(−)-anchored strand"),
            role("one random 24-base spatial barcode"))))
    graph.react(
        "seed_polonies", "Seed and bridge-amplify random barcode polonies",
        StepKind.PCR,
        ((stamp, InputRole.TEMPLATE), (random_templates, InputRole.STRUCTURAL),
         (plus_primer, InputRole.STRUCTURAL), (minus_primer, InputRole.STRUCTURAL),
         (bridge_polymerase, InputRole.CATALYST)), (double_polonies,),
        geometry=PCRGeometry.unavailable(
            double_polonies.state,
            reason="A surface polony family is not one linear amplicon with one fixed random barcode."))
    user = graph.input("user_enzyme", "USER enzyme in CutSmart (37 °C, 1 h)", MaterialKind.CATALYST)
    released = graph.product("user_released_strands", "USER-released PS-BA(+)-anchored strands",
                             MaterialKind.WASTE)
    graph.react(
        "user_linearize", "USER-linearize the polonies at the BA(+) uracil",
        StepKind.DIGESTION,
        ((double_polonies, InputRole.TEMPLATE), (user, InputRole.CATALYST)),
        (polonies, released))
    graph.discard(released, "is washed away (formamide/elution buffer) before stamping")

    copy_gel = graph.input("copy_gel", "primer-functionalized copy gel",
                           MaterialKind.MOLECULAR, primary=True)
    transfer_polymerase = graph.input("transfer_polymerase",
                                     "Taq stamp mix (95 °C contact 1 min, 60 °C 1 min, formamide 95 °C release)",
                                     MaterialKind.CATALYST)
    copied = graph.product(
        "copied_pattern", "copied barcode pattern on a second gel",
        MaterialKind.MOLECULAR,
        state=molecule("copied barcode pattern on a second gel", (
            role("same polony geometry"), role("complementary copied barcode strands"))))
    recovered_stamp = graph.product("recovered_stamp", "reusable original stamp gel",
                                    MaterialKind.MOLECULAR, state=polonies.state)
    graph.react(
        "stamp_transfer", "Press stamp and copy gels together and extend across the interface",
        StepKind.EXTENSION,
        ((polonies, InputRole.TEMPLATE), (copy_gel, InputRole.TEMPLATE),
         (transfer_polymerase, InputRole.CATALYST)), (copied, recovered_stamp))
    graph.mark_terminal(Terminal(
        recovered_stamp.id, Disposition.STORED,
        "the nonconsumptive stamp can generate additional sibling copy gels"))

    copy_plus = graph.input("copy_plus_primer", "copy-gel bridge-PCR plus primer",
                            MaterialKind.OLIGO)
    copy_minus = graph.input("copy_minus_primer", "copy-gel bridge-PCR minus primer",
                             MaterialKind.OLIGO)
    copy_polymerase = graph.input("copy_polymerase", "Bst bridge amplification, 20 cycles",
                                  MaterialKind.CATALYST)
    dense_copy = graph.product("amplified_copy", "bridge-amplified copy-gel polonies",
                               MaterialKind.MOLECULAR, state=copied.state)
    graph.react(
        "amplify_copy", "Bridge-amplify the copied barcode pattern",
        StepKind.PCR,
        ((copied, InputRole.TEMPLATE), (copy_plus, InputRole.STRUCTURAL),
         (copy_minus, InputRole.STRUCTURAL), (copy_polymerase, InputRole.CATALYST)),
        (dense_copy,), geometry=PCRGeometry.unavailable(
            dense_copy.state,
            reason="The surface contains a spatial family of random-barcode polonies, not one fixed amplicon."))

    taqi = graph.input("taqi", "TaqI restriction endonuclease (60 °C, 1.5 h)", MaterialKind.CATALYST)
    capture_array = graph.product(
        "capture_array", "RNA-capture polony array with exposed poly(dT)",
        MaterialKind.MOLECULAR,
        state=molecule("RNA-capture polony array", (
            role("gel-bound spatial barcode"), role("TaqI cleavage boundary"),
            role("exposed poly(dT)"))))
    cleaved = graph.product("taqi_cleavage_products", "released TaqI cleavage products",
                            MaterialKind.WASTE)
    graph.react("taqi_digest", "TaqI-cut copied polonies to expose poly(dT)",
                StepKind.DIGESTION,
                ((dense_copy, InputRole.TEMPLATE), (taqi, InputRole.CATALYST)),
                (capture_array, cleaved))
    graph.discard(cleaved, "is released while the gel-bound poly(dT) capture array is retained")

    tissue_rna = graph.input("tissue_rna", "poly(A) RNA from a fresh 10-µm frozen section",
                             MaterialKind.BIOLOGICAL, primary=True)
    tso = graph.input("template_switch_oligo", "template-switching oligo AAGCTGGTATCAACGCAGAGTGAATrG+GrG",
                      MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "Maxima H-minus reverse transcriptase (42 °C, 1 h)",
                                       MaterialKind.CATALYST)
    barcoded_cdna = graph.product(
        "barcoded_cdna", "gel-bound spatially barcoded first-strand cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("gel-bound first-strand cDNA", (
            role("gel anchor and P5/T7 segment"), role("24-base polony barcode"),
            role("polony sequencing-primer site"), role("poly(T) probe"),
            role("first-strand cDNA (antisense)"), role("TSO complement"))))
    graph.react(
        "capture_rt", "Hybridize tissue RNA for 15 min, then template-switch first-strand synthesis",
        StepKind.REVERSE_TRANSCRIPTION,
        ((capture_array, InputRole.STRUCTURAL), (tissue_rna, InputRole.TEMPLATE),
         (tso, InputRole.STRUCTURAL), (reverse_transcriptase, InputRole.CATALYST)),
        (barcoded_cdna,))

    proteinase = graph.input("proteinase_k", "proteinase K in PKD buffer (55 °C, 30 min)",
                             MaterialKind.CATALYST)
    cleaned = graph.product("cleaned_cdna_gel", "washed gel carrying first-strand cDNA",
                            MaterialKind.MOLECULAR, state=barcoded_cdna.state)
    tissue_debris = graph.product("tissue_debris", "digested proteins, genomic DNA and tissue debris",
                                  MaterialKind.WASTE)
    graph.react(
        "tissue_cleanup", "Proteinase-K digest the tissue and wash with 2× SSC/0.1% SDS and 0.1× SSC",
        StepKind.DEPROTEINIZATION,
        ((barcoded_cdna, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST)),
        (cleaned, tissue_debris))
    graph.discard(tissue_debris, "is washed out of the reaction chamber while cDNA stays gel-bound")

    tso_primer = graph.input("tso_primer", "biotinylated TSO primer /5’BIOSG/AAGCAGTGGTATCAACGCAGAGT",
                             MaterialKind.OLIGO)
    second_polymerase = graph.input("second_strand_polymerase",
                                    "Bst 2.0 WarmStart DNA polymerase (65 °C, 15 min)",
                                    MaterialKind.CATALYST)
    duplex_on_gel = graph.product(
        "gel_duplex", "gel-bound first strand paired with its second strand",
        MaterialKind.MOLECULAR,
        state=molecule("second-strand duplex on gel", (
            role("TSO primer"), role("transcript-sense cDNA"), role("poly(A)"),
            role("polony-site, barcode and T7/P5 complements"))))
    graph.react(
        "second_strand", "Prime second-strand synthesis with the TSO primer",
        StepKind.EXTENSION,
        ((cleaned, InputRole.TEMPLATE), (tso_primer, InputRole.STRUCTURAL),
         (second_polymerase, InputRole.CATALYST)), (duplex_on_gel,))

    elution = graph.input("strand_elution", "0.08 M KOH, 10 min RT, neutralized with Tris pH 7.0",
                          MaterialKind.CONDITION)
    eluted = graph.product(
        "eluted_second_strands", "eluted spatially barcoded second strands",
        MaterialKind.MOLECULAR,
        state=molecule("eluted second strand", (
            role("TSO primer"), role("transcript-sense cDNA"), role("poly(A)"),
            role("polony-site complement"), role("24-base barcode complement"),
            role("T7-promoter and P5 complement"))))
    spent_gel = graph.product("spent_capture_gel", "gel retaining covalent first strands",
                              MaterialKind.WASTE)
    graph.react("elute", "Denature with KOH and recover the second strands",
                StepKind.STRAND_SEPARATION,
                ((duplex_on_gel, InputRole.TEMPLATE), (elution, InputRole.CONDITION)),
                (eluted, spent_gel))
    graph.discard(spent_gel, "stays in the chamber; about 35 µL eluate is transferred to a tube")

    umi_primer = graph.input("umi_primer", "published seven-random-base UMI primer",
                             MaterialKind.OLIGO,
                             state=Construct([seg("UMI primer", UMI_PRIMER)]))
    extension_polymerase = graph.input("umi_extension_polymerase",
                                       "Q5 Ultra II master mix (95 °C 30 s, 65 °C 30 s, 72 °C 5 min)",
                                       MaterialKind.CATALYST)
    umi_cdna = graph.product(
        "umi_cdna", "barcode–UMI–transcript cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("UMI-extended cDNA", (
            role("TruSeq Read 1 arm"), role("seven-base UMI"), role("T7-promoter anchor"),
            role("24-base spatial barcode"), role("polony sequencing-primer site"),
            role("poly(T)"), role("cDNA (antisense)"), role("TSO complement"))))
    graph.react(
        "umi_extension", "Anneal the UMI primer to the T7-promoter complement and extend once",
        StepKind.EXTENSION,
        ((eluted, InputRole.TEMPLATE), (umi_primer, InputRole.STRUCTURAL),
         (extension_polymerase, InputRole.CATALYST)), (umi_cdna,))

    exo = graph.input("thermolabile_exo_i", "thermolabile Exonuclease I (37 °C 4 min, 80 °C 1 min)",
                      MaterialKind.CATALYST)
    exo_cleaned = graph.product("exo_cleaned_cdna", "UMI-extended cDNA after primer removal",
                                MaterialKind.MOLECULAR, state=umi_cdna.state)
    digested_primer = graph.product("digested_umi_primer", "Exonuclease I-digested unincorporated UMI primer",
                                    MaterialKind.WASTE)
    graph.react(
        "exo_cleanup", "Digest non-incorporated UMI primer with thermolabile Exonuclease I",
        StepKind.DIGESTION,
        ((umi_cdna, InputRole.TEMPLATE), (exo, InputRole.CATALYST)),
        (exo_cleaned, digested_primer))
    graph.mark_terminal(Terminal(
        digested_primer.id, Disposition.UNRESOLVED,
        "degraded in the same tube and carried into PCR; no physical separation is described"))

    universal_forward = graph.input("universal_forward", "TSO primer (Methods: 'TSO' primer for cDNA PCR)",
                                    MaterialKind.OLIGO)
    universal_reverse = graph.input("universal_reverse",
                                    "TruSeq sequencing primer CCTACACGACGCTCTTCCGATC",
                                    MaterialKind.OLIGO,
                                    state=Construct([seg("TruSeq sequencing primer",
                                                         TRUSEQ_SEQUENCING_PRIMER)]))
    universal_polymerase = graph.input("universal_polymerase",
                                       "Q5 HotStart polymerase (4 + 8 cycles for an OB section)",
                                       MaterialKind.CATALYST)
    amplified = graph.product("amplified_cdna", "amplified barcode–UMI–transcript cDNA (5–10 ng)",
                              MaterialKind.MOLECULAR, state=umi_cdna.state)
    graph.react(
        "universal_pcr", "PCR-amplify cDNA with the TSO and TruSeq sequencing primers", StepKind.PCR,
        ((exo_cleaned, InputRole.TEMPLATE), (universal_forward, InputRole.STRUCTURAL),
         (universal_reverse, InputRole.STRUCTURAL),
         (universal_polymerase, InputRole.CATALYST)), (amplified,),
        geometry=PCRGeometry.unavailable(
            amplified.state,
            reason="The amplified pool contains variable tissue transcripts and random spatial/UMI sequences."))

    transposome = graph.input("nextera_transposome", "Nextera XT (FC-131-1024) tagmentation mix",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Nextera XT Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product("tagmented_cdna", "Nextera-tagmented Pixel-seq cDNA pool (~1 ng input)",
                              MaterialKind.MOLECULAR, state=amplified.state)
    tn5_losses = graph.product("tagmentation_losses", "non-library tagmentation products",
                               MaterialKind.WASTE)
    graph.react(
        "tagment", "Tagment about 1 ng amplified cDNA with Nextera XT", StepKind.TAGMENTATION,
        ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
         (tn5, InputRole.CATALYST)), (tagmented, tn5_losses))
    graph.discard(tn5_losses, "does not form an amplifiable TruSeq-to-Nextera molecule")

    index_forward = graph.input("index_forward", "TruSeq libP5 primer",
                                MaterialKind.OLIGO,
                                state=Construct([seg("TruSeq libP5 primer", TRUSEQ_LIB_P5)]))
    index_reverse = graph.input("index_reverse", "Nextera index primer CAAGCAGAAGACGGCATACGAGAT[INDEX]GTCTCGTGGGCTCGG",
                                MaterialKind.OLIGO)
    index_polymerase = graph.input("index_polymerase", "Nextera XT index-PCR polymerase",
                                  MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "informative Pixel-seq library",
                            MaterialKind.LIBRARY, state=final)
    uninformative = graph.product(
        "uninformative_fragments", "fragments lacking the TruSeq barcode end or a Nextera i7 end",
        MaterialKind.WASTE)
    graph.react(
        "index_pcr", "PCR with the TruSeq libP5 primer and a Nextera i7 index primer",
        StepKind.PCR,
        ((tagmented, InputRole.TEMPLATE), (index_forward, InputRole.STRUCTURAL),
         (index_reverse, InputRole.STRUCTURAL), (index_polymerase, InputRole.CATALYST)),
        (library, uninformative), geometry=PCRGeometry.unavailable(
            final,
            reason="Random tagmentation makes a population; the final construct depicts the TruSeq-barcode-to-Nextera subset that both primers amplify exponentially."))
    graph.discard(uninformative, "is not exponentially amplified or lacks a spatial barcode, UMI or transcript")
    graph.sequence(library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers)))
    graph.annotate(Highlight(
        "index_pcr", "Single-indexed TruSeq/Nextera hybrid (model correction)",
        "STAR Methods use the TruSeq libP5 primer and Nextera index primers after Nextera XT "
        "tagmentation. The P5 end is the UMI primer's TruSeq Read 1 arm without an i5 "
        "window; only the Tn5 end carries an i7 index. An earlier repository model drew a "
        "dual-indexed Nextera library with the UMI before the barcode. The inspected author "
        "manuscript prints neither the i7 length nor run read lengths.",
        AnnotationSlot.SEQUENCING))
    graph.declare_reviewed_scope((
        "stamp_gel", "random_barcode_templates", "stamp_plus_primer", "stamp_minus_primer",
        "stamp_polymerase", "seed_polonies", "user_enzyme", "user_linearize", "copy_gel", "transfer_polymerase",
        "stamp_transfer", "copy_plus_primer", "copy_minus_primer", "copy_polymerase",
        "amplify_copy", "taqi", "taqi_digest", "tissue_rna", "template_switch_oligo",
        "reverse_transcriptase", "capture_rt", "proteinase_k", "tissue_cleanup",
        "tso_primer", "second_strand_polymerase", "second_strand", "strand_elution",
        "elute", "umi_primer", "umi_extension_polymerase", "umi_extension",
        "thermolabile_exo_i", "exo_cleanup", "universal_forward",
        "universal_reverse", "universal_polymerase", "universal_pcr",
        "nextera_transposome", "tn5", "tagment", "index_forward", "index_reverse",
        "index_polymerase", "index_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def oligos():
    return (
        f"Published polony seed template (24 random barcode bases) — {inline_sequence(POLONY_TEMPLATE, 'Pixel-seq polony template')}",
        f"Bridge-amplification BA(+), as written — {BA_PLUS}",
        f"Bridge-amplification BA(−), as written — {BA_MINUS}",
        f"Gel-grafted PS-BA(+), as written — {PS_BA_PLUS}",
        f"Gel-grafted PS-BA(−), as written — {PS_BA_MINUS}",
        f"Polony Cy5 probe, as written — {CY5_PROBE}",
        f"UMI primer — {inline_sequence(UMI_PRIMER, 'Pixel-seq UMI primer')}",
        f"Template-switching oligo, as written — {TSO}",
        f"Biotinylated TSO primer, as written — {TSO_PRIMER}",
        f"TruSeq sequencing primer (cDNA PCR) — {inline_sequence(TRUSEQ_SEQUENCING_PRIMER, 'TruSeq sequencing primer')}",
        f"TruSeq library P5 primer — {inline_sequence(TRUSEQ_LIB_P5, 'TruSeq libP5')}",
        f"Nextera index primer, as written — {NEXTERA_INDEX_PRIMER}",
        f"Polony sequencing primer — {inline_sequence(POLONY_SEQUENCING_PRIMER, 'polony sequencing primer')}",
    )


def get():
    final, primers = pixel_library()
    spec = SimpleNamespace(
        TITLE="Pixel-seq — stamped polony-gel spatial transcriptomics",
        NOTES="01_pixel_seq.html",
        SOURCE='<a href="https://doi.org/10.1016/j.cell.2022.10.021">Fu et al., <i>Cell</i> (2022)</a>.',
        SUMMARY="A gel bearing clonal random-barcode polonies enzymatically stamps its complete coordinate pattern onto many copy gels. TaqI exposes poly(dT), tissue RNA is reverse-transcribed on the replicated array, and a TSO-primed second strand is eluted, a UMI primer is extended across the barcode, and Nextera XT tagmentation with TruSeq libP5 and Nextera i7 index primers completes the library.",
        CAVEAT="The displayed final molecule is the selected sequencing product that retains spatial barcode, UMI and transcript information after Nextera XT tagmentation; random tagmentation also makes products lacking one of those analytical regions. INFERRED — the Nextera index primer is printed as [INDEX] without a length, so the eight-position i7 window is a placeholder. Run read lengths are not printed in the inspected author manuscript and are not assigned.",
        oligos=oligos, workflow=workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="Representative informative Pixel-seq fragment: TruSeq libP5 end, UMI-primer UMI and T7-promoter anchor, the template's 24-base polony barcode, polony sequencing-primer site and poly(T), then antisense cDNA ending at a Nextera i7 end. The i7 window length is a placeholder (INFERRED).",
        SEQ_PRIMERS=primers,
        READ_LENGTHS=None,
        SEQUENCING_INTRO="TruSeq Read 1 reads the UMI, the T7-promoter anchor and then the spatial barcode; Read 2 starts at the Nextera mosaic end inside the transcript. The inspected author manuscript names the primers but not the instrument run lengths (it cites NovaSeq S4 only for a cost estimate). Barcode and UMI are extracted with Flexbar, barcodes mapped to the polony map with up to two mismatches, and paired reads aligned with STAR.",
        SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="",
    )
    spec.protocol_graph = lambda: _protocol_graph(final, primers)
    return spec


SPEC = get()
globals().update(vars(SPEC))


def protocol_graph():
    return _protocol_graph(SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)
