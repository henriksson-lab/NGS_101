"""Reusable polony-gel stamping and Pixel-seq library chemistry."""
from __future__ import annotations

from types import SimpleNamespace

from batch_ngs import nextera_library, seg
from chemdraw import feature, inline_sequence
from protocol_model import chain, molecule, role


UMI_PRIMER = "CTACACGACGCTCTTCCGATCTNNNNNNNCCGCTAATACGACTCACTATAGGG"
TRUSEQ_SEQUENCING_PRIMER = "CCTACACGACGCTCTTCCGATC"
TSO_PRIMER = "/5’BIOSG/AAGCAGTGGTATCAACGCAGAGT"
TSO = "AAGCTGGTATCAACGCAGAGTGAATrG+GrG"
TRUSEQ_LIB_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
POLONY_SEQUENCING_PRIMER = "TTCCATCATTGTTCCCGGGTTCCTCATTCTC"


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


def oligos():
    return (
        f"UMI primer — {inline_sequence(UMI_PRIMER, 'Pixel-seq UMI primer')}",
        f"Template-switching oligo, as written — {TSO}",
        f"Biotinylated TSO primer, as written — {TSO_PRIMER}",
        f"TruSeq library P5 primer — {inline_sequence(TRUSEQ_LIB_P5, 'TruSeq libP5')}",
        f"Polony sequencing primer — {inline_sequence(POLONY_SEQUENCING_PRIMER, 'polony sequencing primer')}",
    )


def get():
    final, primers = pixel_library()
    return SimpleNamespace(
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
