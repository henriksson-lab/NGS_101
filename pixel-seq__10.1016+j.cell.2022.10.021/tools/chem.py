"""Pixel-seq polony-gel stamping and library chemistry."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))

from batch_ngs import nextera_library, seg
from chemdraw import feature, inline_sequence
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from protocol_model import chain, molecule, role


UMI_PRIMER = "CTACACGACGCTCTTCCGATCTNNNNNNNCCGCTAATACGACTCACTATAGGG"
TRUSEQ_SEQUENCING_PRIMER = "CCTACACGACGCTCTTCCGATC"
TSO_PRIMER = "/5’BIOSG/AAGCAGTGGTATCAACGCAGAGT"
TSO = "AAGCTGGTATCAACGCAGAGTGAATrG+GrG"
TRUSEQ_LIB_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
POLONY_SEQUENCING_PRIMER = "TTCCATCATTGTTCCCGGGTTCCTCATTCTC"
BA_PLUS = "AATGATACGGCGACCACCGAGAUCTACAC"
BA_MINUS = "CAAGCAGAAGACGGCATACGAGAT"
PS_BA_PLUS = "/5PS/TTTTTTTTTTAATGATACGGCGACCACCGAGAUCTACAC"
PS_BA_MINUS = "/5PS/TTTTTTTTTTCAAGCAGAAGACGGCATACGAGAT"
CY5_PROBE = "/5Cy5/CTTGCGCTAGTCCGGCAACCTCGGTGG"
NEXTERA_INDEX_PRIMER = "CAAGCAGAAGACGGCATACGAGAT[INDEX]GTCTCGTGGGCTCGG"
POLONY_TEMPLATE = (
    "AATGATACGGCGACCACCGAGATCTACAC"
    "TGCGGCCGCTAATACGACTCACTATAGGGATCT"
    "NNNNNNNNNNNNNNNNNNNNNNNN"
    "AGAGAATGAGGAACCCGGGAACAATGATGGAA"
    "TTTTTTTTTTTTTTTTTTT"
    "TCGA"
    "ACCACCGAGGTTGCCGGACTAGCGCAAG"
    "TACTTGTCCATTCCTGAAGAAATATTATATTTATACAACTTACCCATAGAATCCTATTTACTAGGAAAGGAAAAGCCTCCTATTTATAAAAATTGGATAGAGCTTTCTCAACAACAGTGGAATATCAATGATAGAACAATTGCCGATTTATTAGATGGGGTCTTAATAATACCATCGA"
    "TCTCGTATGCCGTCTTCTGCTTG"
)


SPATIAL = feature("pixel_seq.spatial", "spatial_barcode", "random",
                  note="24-base polony identity mapped once per replicated gel family")
UMI = feature("pixel_seq.umi", "umi", "random",
              note="seven random bases introduced by the published UMI primer")


def pixel_library():
    return nextera_library([
        seg("24-nt polony spatial barcode", "B" * 24, "cbc", placeholder=True,
            feature=SPATIAL),
        seg("7-nt UMI", "U" * 7, "umi", placeholder=True, feature=UMI),
        role("captured transcript cDNA"),
    ], "representative Pixel-seq Nextera XT library", dual_index=True)


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
        ("Elute second strands, extend the seven-random-base UMI primer and PCR-amplify",
         molecule("amplified barcode–UMI–transcript cDNA", (
             role("24-base spatial barcode"), role("seven-base UMI"),
             role("captured transcript cDNA"), role("universal PCR handles"))),
         "A single Q5 extension installs the UMI before Exonuclease I cleanup; PCR then uses the TSO and TruSeq primers."),
        ("Tagment 1 ng amplified cDNA with Nextera XT and add dual sample indexes",
         molecule("paired-end Pixel-seq library", (
             role("P5, i5, S5 and mosaic end"), role("spatial barcode and UMI"),
             role("transcript insert"), role("mosaic end, S7, i7 and P7"))),
         "The defining libraries were sequenced PE150 on NovaSeq. Only fragments retaining a mapped spatial barcode, UMI and transcript alignment contribute spatial molecules."),
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
    random_templates = graph.input("random_barcode_templates", "random 24-base barcode templates",
                                   MaterialKind.OLIGO)
    plus_primer = graph.input("stamp_plus_primer", "gel-grafted plus primer",
                              MaterialKind.OLIGO)
    minus_primer = graph.input("stamp_minus_primer", "gel-grafted minus primer",
                               MaterialKind.OLIGO)
    bridge_polymerase = graph.input("stamp_polymerase", "stamp-gel bridge-PCR polymerase",
                                    MaterialKind.CATALYST)
    polonies = graph.product(
        "stamp_polonies", "stamp gel with clonal barcode polonies",
        MaterialKind.MOLECULAR,
        state=molecule("stamp gel with clonal barcode polonies", (
            role("gel-bound primer"), role("one random 24-base spatial barcode"),
            role("USER-linearizable opposite strand"))))
    graph.react(
        "seed_polonies", "Seed and bridge-amplify random barcode polonies",
        StepKind.PCR,
        ((stamp, InputRole.TEMPLATE), (random_templates, InputRole.STRUCTURAL),
         (plus_primer, InputRole.STRUCTURAL), (minus_primer, InputRole.STRUCTURAL),
         (bridge_polymerase, InputRole.CATALYST)), (polonies,),
        geometry=PCRGeometry.unavailable(
            polonies.state,
            reason="A surface polony family is not one linear amplicon with one fixed random barcode."))

    copy_gel = graph.input("copy_gel", "primer-functionalized copy gel",
                           MaterialKind.MOLECULAR, primary=True)
    transfer_polymerase = graph.input("transfer_polymerase", "gel-stamping extension polymerase",
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
    copy_polymerase = graph.input("copy_polymerase", "copy-gel bridge-PCR polymerase",
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

    taqi = graph.input("taqi", "TaqI restriction endonuclease", MaterialKind.CATALYST)
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

    tissue_rna = graph.input("tissue_rna", "poly(A) RNA from overlaid tissue",
                             MaterialKind.BIOLOGICAL, primary=True)
    tso = graph.input("template_switch_oligo", "published template-switch oligo",
                      MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "Maxima H-minus reverse transcriptase",
                                       MaterialKind.CATALYST)
    barcoded_cdna = graph.product(
        "barcoded_cdna", "gel-bound spatially barcoded cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("spatially barcoded cDNA", (
            role("24-base polony barcode"), role("poly(dT)"),
            role("captured transcript cDNA"), role("template-switch handle"))))
    graph.react(
        "capture_rt", "Hybridize tissue RNA and template-switch first-strand synthesis",
        StepKind.REVERSE_TRANSCRIPTION,
        ((capture_array, InputRole.STRUCTURAL), (tissue_rna, InputRole.TEMPLATE),
         (tso, InputRole.STRUCTURAL), (reverse_transcriptase, InputRole.CATALYST)),
        (barcoded_cdna,))

    elution = graph.input("strand_elution", "second-strand elution condition",
                          MaterialKind.CONDITION)
    eluted = graph.product("eluted_second_strands", "eluted spatially barcoded second strands",
                           MaterialKind.MOLECULAR, state=barcoded_cdna.state)
    spent_gel = graph.product("spent_capture_gel", "spent RNA-capture gel",
                              MaterialKind.WASTE)
    graph.react("elute", "Elute spatially barcoded second strands",
                StepKind.STRAND_SEPARATION,
                ((barcoded_cdna, InputRole.TEMPLATE), (elution, InputRole.CONDITION)),
                (eluted, spent_gel))
    graph.discard(spent_gel, "the soluble barcode–transcript strand is the library substrate")

    umi_primer = graph.input("umi_primer", "published seven-random-base UMI primer",
                             MaterialKind.OLIGO)
    extension_polymerase = graph.input("umi_extension_polymerase", "Q5 UMI-extension polymerase",
                                      MaterialKind.CATALYST)
    umi_cdna = graph.product(
        "umi_cdna", "barcode–UMI–transcript cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("amplified barcode–UMI–transcript cDNA", (
            role("24-base spatial barcode"), role("seven-base UMI"),
            role("captured transcript cDNA"), role("universal PCR handles"))))
    graph.react(
        "umi_extension", "Extend the seven-random-base UMI primer",
        StepKind.EXTENSION,
        ((eluted, InputRole.TEMPLATE), (umi_primer, InputRole.STRUCTURAL),
         (extension_polymerase, InputRole.CATALYST)), (umi_cdna,))

    universal_forward = graph.input("universal_forward", "TSO-derived universal PCR primer",
                                    MaterialKind.OLIGO)
    universal_reverse = graph.input("universal_reverse", "TruSeq-derived universal PCR primer",
                                    MaterialKind.OLIGO)
    universal_polymerase = graph.input("universal_polymerase", "cDNA amplification polymerase",
                                      MaterialKind.CATALYST)
    amplified = graph.product("amplified_cdna", "amplified barcode–UMI–transcript cDNA",
                              MaterialKind.MOLECULAR, state=umi_cdna.state)
    graph.react(
        "universal_pcr", "PCR-amplify barcode–UMI–transcript cDNA", StepKind.PCR,
        ((umi_cdna, InputRole.TEMPLATE), (universal_forward, InputRole.STRUCTURAL),
         (universal_reverse, InputRole.STRUCTURAL),
         (universal_polymerase, InputRole.CATALYST)), (amplified,),
        geometry=PCRGeometry.unavailable(
            amplified.state,
            reason="The amplified pool contains variable tissue transcripts and random spatial/UMI sequences."))

    transposome = graph.input("nextera_transposome", "loaded Nextera XT adapter transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Nextera XT Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product("tagmented_cdna", "Nextera-tagmented Pixel-seq cDNA pool",
                              MaterialKind.MOLECULAR, state=amplified.state)
    tn5_losses = graph.product("tagmentation_losses", "non-library tagmentation products",
                               MaterialKind.WASTE)
    graph.react(
        "tagment", "Tagment amplified cDNA with Nextera XT", StepKind.TAGMENTATION,
        ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
         (tn5, InputRole.CATALYST)), (tagmented, tn5_losses))
    graph.discard(tn5_losses, "does not form an amplifiable two-ended Nextera molecule")

    index_forward = graph.input("index_forward", "Nextera XT i5 indexing primer",
                                MaterialKind.OLIGO)
    index_reverse = graph.input("index_reverse", "Nextera XT i7 indexing primer",
                                MaterialKind.OLIGO)
    index_polymerase = graph.input("index_polymerase", "Nextera XT index-PCR polymerase",
                                  MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "informative Pixel-seq library",
                            MaterialKind.LIBRARY, state=final)
    uninformative = graph.product(
        "uninformative_fragments", "amplified fragments lacking barcode, UMI or transcript",
        MaterialKind.WASTE)
    graph.react(
        "index_pcr", "Add dual sample indexes and select informative molecules",
        StepKind.PCR,
        ((tagmented, InputRole.TEMPLATE), (index_forward, InputRole.STRUCTURAL),
         (index_reverse, InputRole.STRUCTURAL), (index_polymerase, InputRole.CATALYST)),
        (library, uninformative), geometry=PCRGeometry.unavailable(
            final,
            reason="Random tagmentation makes a population; the final construct depicts the informative barcode–UMI–transcript subset."))
    graph.discard(uninformative, "lacks one or more analytical regions needed for a spatial molecule")
    graph.sequence(library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers),
        {"Read 1": 150, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 150}))
    graph.declare_reviewed_scope((
        "stamp_gel", "random_barcode_templates", "stamp_plus_primer", "stamp_minus_primer",
        "stamp_polymerase", "seed_polonies", "copy_gel", "transfer_polymerase",
        "stamp_transfer", "copy_plus_primer", "copy_minus_primer", "copy_polymerase",
        "amplify_copy", "taqi", "taqi_digest", "tissue_rna", "template_switch_oligo",
        "reverse_transcriptase", "capture_rt", "strand_elution", "elute", "umi_primer",
        "umi_extension_polymerase", "umi_extension", "universal_forward",
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
        SUMMARY="A gel bearing clonal random-barcode polonies enzymatically stamps its complete coordinate pattern onto many copy gels. TaqI exposes poly(dT), tissue RNA is reverse-transcribed on the replicated array, and barcode–UMI–transcript cDNA is completed as a Nextera XT library.",
        CAVEAT="The displayed final molecule is the selected sequencing product that retains spatial barcode, UMI and transcript information after Nextera XT tagmentation; random tagmentation also makes products lacking one of those analytical regions.",
        oligos=oligos, workflow=workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="Representative informative Pixel-seq fragment after Nextera XT. The 24-base random barcode maps the molecule to a polony coordinate and the adjacent seven-base random sequence identifies the captured molecule.",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 8,
                      "Index 2 (i5)": 8, "Read 2": 150},
        SEQUENCING_INTRO="The defining cDNA libraries were run PE150 on NovaSeq. Paired genomic/transcript sequence is retained only when the same library molecule also carries a spatial barcode matching the in-situ polony map and a UMI.",
        SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="",
    )
    spec.protocol_graph = lambda: _protocol_graph(final, primers)
    return spec


SPEC = get()
globals().update(vars(SPEC))


def protocol_graph():
    return _protocol_graph(SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)
