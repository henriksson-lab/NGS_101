"""DR-seq co-amplification and source-backed branch intermediates."""
from __future__ import annotations

from batch_ngs import seg
from chemdraw import (Construct, MolecularState, Row, Scene, Workflow,
                      complement_segments, duplex_rows, feature)
from multibranch import illumina_branch

TITLE = "DR-seq — genome and transcriptome without physical separation"
NOTES = "01_dr-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nbt.3129">'
          'Dey et al., <i>Nature Biotechnology</i> (2015)</a>.')
SUMMARY = ("A cell-barcoded oligo(dT) first makes cDNA in the same lysate as genomic DNA. "
           "Ad-2 random primers then quasilinearly copy both analytes before the pool is "
           "split into genome and transcriptome finishing reactions.")
CAVEAT = ("The paper prints Ad-1x, Ad-2 and the branch-selection primers, but refers the "
          "RNA branch to the original CEL-Seq small-RNA finishing workflow without printing "
          "those outer library oligos. The final shells therefore show inferred standard "
          "Illumina adapters; the co-amplification intermediates are source-exact.")

AD1_SPACER = "CGATTGAGGCCGG"
T7_PROMOTER = "TAATACGACTCACTATAGGG"
RA5 = "GTTCAGAGTTCTACAGTCCGACGATC"
CELL_BARCODE_NT = 8
POLYT_NT = 24
AD2_HANDLE = "GTGAGTGATGGTTGAGGTAGTGTGGAG"
AD2_RANDOM_NT = 8
P1 = "CGATTGAGGCCGGTAATAC"
P2 = AD2_HANDLE
P3 = "GTGAGCTGGAGTTGAGGTAGTGTGGAG"


def ad1_segments():
    return [seg("Ad-1x spacer", AD1_SPACER), seg("T7 promoter", T7_PROMOTER, "t7"),
            seg("RA5", RA5, "r1"),
            seg("cell barcode", "B" * CELL_BARCODE_NT, "cbc", placeholder=True,
                feature=feature("cell_barcode", "cell_barcode", "fixed")),
            seg("oligo(dT)24", "T" * POLYT_NT),
            seg("V anchor", "V", placeholder=True)]


def ad2_segments():
    return [seg("Ad-2 common handle", AD2_HANDLE, "tso"),
            seg("Ad-2 random N8", "N" * AD2_RANDOM_NT, placeholder=True)]


def rt_scene() -> Scene:
    mrna = [seg("transcript", "X" * 34, placeholder=True),
            seg("poly(A)", "A" * POLYT_NT)]
    sc = Scene(); sc.strand("mRNA", mrna, label="poly(A) RNA")
    sc.anneal("Ad-1x", ad1_segments(), to="mRNA",
              pair=("oligo(dT)24", "poly(A)"), label="cell-barcoded Ad-1x",
              unpaired=("V anchor",))
    sc.arrow("Ad-1x", "reverse transcription")
    sc.mark("Ad-1x", "cell barcode", "8-nt cell identity")
    sc.labels("Ad-1x")
    return sc


def ad2_priming_scene(template_name: str) -> Scene:
    template = [seg(template_name, "X" * 40, placeholder=True),
                seg("random priming site", "NNNNNNNN", placeholder=True)]
    sc = Scene(); sc.strand("template", template, label=template_name)
    sc.anneal("Ad-2", ad2_segments(), to="template",
              pair=("Ad-2 random N8", "random priming site"), label="Ad-2 random primer")
    sc.arrow("Ad-2", "strand-displacing extension")
    sc.labels("Ad-2")
    return sc


def gdna_amplicon() -> Construct:
    left = [*ad2_segments(), seg("genomic copy", "X" * 40, placeholder=True)]
    return Construct([*left, *complement_segments(ad2_segments(), suffix=" right complement")],
                     name="Ad-2/Ad-2 genomic full amplicon")


def cdna_amplicon() -> Construct:
    # The productive transcript-derived species is Ad-2 at its random-primed end and
    # retains Ad-1x at the poly(A)-primed end.
    return Construct([*ad2_segments(), seg("cDNA", "X" * 38, placeholder=True),
                      seg("copied poly(A)", "A" * POLYT_NT),
                      seg("cell barcode complement", "b" * CELL_BARCODE_NT, "cbc",
                          placeholder=True,
                          feature=feature("cell_barcode", "cell_barcode", "fixed")),
                      *complement_segments(ad1_segments()[:3], suffix=" complement")],
                     name="Ad-2/Ad-1x transcript amplicon")


GDNA_AMPLICON = gdna_amplicon()
CDNA_AMPLICON = cdna_amplicon()


def arna_product() -> Construct:
    return Construct([
        seg("T7 initial G-rich end", "GGG", "t7"), seg("RA5", RA5, "r1"),
        seg("cell barcode", "B" * CELL_BARCODE_NT, "cbc", placeholder=True,
            feature=feature("cell_barcode", "cell_barcode", "fixed")),
        seg("poly(U) and antisense transcript", "UUU" + "X" * 38, placeholder=True),
        *complement_segments(ad2_segments(), suffix=" complement"),
    ], name="cell-barcoded antisense RNA")


ARNA = arna_product()


def arna_rows():
    sc = Scene(); sc.strand("aRNA", list(ARNA), label="cell-barcoded aRNA")
    sc.labels("aRNA")
    return sc.rows()

# The exact post-sonication DNA insert and exact post-IVT RNA insert are retained, while
# outer Illumina completion is marked inferred because those oligos are not printed here.
DNA, DP = illumina_branch("DR-seq genome", "sonicated genomic insert",
                          architecture="truseq", inferred_adapters=True)
RNA, RP = illumina_branch("DR-seq transcriptome", "3-prime RNA-derived insert",
                          cell=CELL_BARCODE_NT, inferred_adapters=True,
                          encoding="fixed")
FINAL_LIBRARIES = (
    ("Genome library", DNA, DP,
     "NEBNext Ultra library from sonicated, Ad-2-depleted genomic amplicons; outer shell is inferred because its lot-specific oligos are not printed.",
     "Canonical TruSeq-family primer landing on the inferred NEBNext shell."),
    ("Transcriptome library", RNA, RP,
     "Directional CEL-Seq-derived library retaining the 8-nt Ad-1x cell barcode; outer small-RNA completion is shown as inferred.",
     "Primer landing is shown on the inferred Illumina shell; see the evidence caveat."),
)


def oligos():
    return (
        f"<b>Ad-1x</b> — 5′-{AD1_SPACER}{T7_PROMOTER}{RA5}[8-nt cell barcode]{'T' * POLYT_NT}V-3′",
        f"<b>Ad-2</b> — 5′-{AD2_HANDLE}[8 random nt]-3′",
        f"<b>P1</b> — 5′-{P1}-3′",
        f"<b>P2</b> — 5′-{P2}-3′",
        f"<b>5′-biotin P3</b> — 5′-biotin-{P3}-3′",
    )


def workflow() -> Workflow:
    start = Scene.duplex([seg("genomic DNA", "X" * 42, placeholder=True)],
                         label="single-cell genomic DNA")
    wf = Workflow(MolecularState("One lysed cell: genomic DNA plus poly(A) RNA",
                                 tuple(start.rows())))
    wf.react("Reverse-transcribe RNA with cell-barcoded Ad-1x",
             rt_scene().rows(), name="Ad-1x-tagged cDNA beside genomic DNA",
             note="The 5′ tail installs the T7 promoter, RA5 and cell barcode; genomic DNA remains in the tube.")
    wf.react("Random-prime genomic DNA and cDNA with Ad-2",
             [*ad2_priming_scene("genomic DNA").rows(),
              *ad2_priming_scene("Ad-1x-tagged cDNA").rows()],
             name="First Ad-2 semi-amplicons")
    wf.react("Run seven MALBAC-like quasilinear cycles",
             [*duplex_rows(GDNA_AMPLICON, label="genome-derived Ad-2/Ad-2 species"),
              *duplex_rows(CDNA_AMPLICON, label="transcript-derived Ad-2/Ad-1x species"),
              Row(chunks=[("full Ad-2/Ad-2 amplicons self-anneal through complementary ends and leave the quasilinear pool", None, False)])],
             name="Co-amplified genome and transcript products")
    wf.react("Split the common amplification pool into DNA and RNA aliquots",
             [Row(chunks=[("DNA aliquot → P2 PCR → 5′-biotin P3 PCR → sonicate → streptavidin depletion", None, False)]),
              Row(chunks=[("RNA aliquot → P1 second strand → T7 IVT", None, False)])],
             name="Modality-specific finishing reactions")
    wf.react("Select the RNA branch with P1 and amplify it by T7 transcription",
             [*duplex_rows(CDNA_AMPLICON, label="P1-completed T7 template"),
              Row(chunks=[("T7 RNA polymerase →", "t7", False)]),
              *arna_rows()],
             name="Antisense RNA copies from transcript-derived templates",
             note="Genomic Ad-2/Ad-2 products lack the Ad-1x T7 end and do not make this aRNA.")
    wf.react("Finish the DNA and RNA aliquots as separate Illumina libraries",
             [*duplex_rows(DNA, label="genome library"),
              *duplex_rows(RNA, label="transcriptome library")],
             name="DR-seq genome and transcriptome libraries",
             note="The source names NEBNext Ultra for DNA and the CEL-Seq workflow for RNA; the unprinted outer oligos are marked inferred.")
    return wf
