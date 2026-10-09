"""Human NET-seq library model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, circle_rows, revcomp
from circular import circularize_ssdna

TITLE = "Human NET-seq — nascent RNA 3′ ends"
NOTES = "01_human-net-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.03.010">Mayer et al., <i>Cell</i> (2015)</a>; '
          '<a href="https://doi.org/10.1038/nprot.2016.086">detailed protocol</a>.')
SUMMARY = ("Ligate a pre-adenylated, six-base-barcoded DNA linker to nascent RNA 3′ ends, reverse-transcribe, circularize the cDNA, deplete abundant mature-RNA products, and amplify for single-read sequencing.")
CAVEAT = ("The protocol names a commercial Illumina Index forward primer without printing its sequence; that outer PCR arm is shown as a role token. The custom read-primer site and molecule-facing geometry are exact.")

LINKER = "N" * 6 + "CTGTAGGCACCATCAAT"
RT_PRIMER = "ATCTCGTATGCCGTCTTCTGCTTG" + "CACTCA" + "TCCGACGATCATTGATGGTGCCTACAG"
REVERSE_PCR = "CAAGCAGAAGACGGCATACGA"
SEQUENCING = "TCCGACGATCATTGATGGTGCCTACAG"

def circular_cdna():
    linear = Construct([
        Segment("RT-primer body", RT_PRIMER[:-len(SEQUENCING)]),
        Segment("custom sequencing-primer site", SEQUENCING, "r1"),
        Segment("six-base molecular barcode", "N" * 6, "umi", placeholder=True),
        Segment("nascent-RNA 3′-end cDNA", "X" * 36, placeholder=True),
    ], name="human NET-seq first-strand cDNA")
    return circularize_ssdna(linear, five_prime_phosphate=True, three_prime_oh=True)

CIRCLE = circular_cdna()

def final_library() -> tuple[Construct, tuple]:
    lib = Construct([
        Segment("Illumina P5/index PCR arm", "[Illumina P5 + index]", "p5+cbc", placeholder=True),
        Segment("custom sequencing-primer site", SEQUENCING, "r1"),
        Segment("six-base molecular barcode", "N" * 6, "umi", placeholder=True),
        Segment("nascent-RNA 3′-end cDNA", "X" * 36, placeholder=True),
        Segment("P7-side reverse-primer complement", revcomp(REVERSE_PCR), "p7"),
    ], name="human NET-seq library")
    primer = sp.custom("Read 1", "oLSC006", SEQUENCING, "Mayer and Churchman 2016, Table 1")
    problems = sp.verify(lib, (primer,), required_roles=("Read 1",))
    if problems: raise ValueError("invalid human NET-seq endpoint: " + "; ".join(problems))
    return lib, (primer,)

FINAL_LIBRARY, SEQ_PRIMERS = final_library()
FINAL_CAPTION = ("The custom primer reads the six-base molecular barcode first, followed immediately by cDNA copied from the nascent RNA 3′ end.")
SEQUENCING_INTRO = ("oLSC006 is the published custom sequencing primer. The run is single-read; its first genomic base marks the engaged polymerase.")

def sections():
    lig = Construct([
        Segment("nascent RNA", "X" * 28, placeholder=True),
        Segment("ligation junction", "N", placeholder=True),
        Segment("six-base molecular barcode", "N" * 6, "umi", placeholder=True),
        Segment("linker common region", LINKER[6:]),
    ], name="linker-ligated nascent RNA")
    sc = Scene.duplex(list(lig), label="schematic RNA/linker product", unpaired=tuple(s.name + "'" for s in lig))
    sc.junction("top", "nascent RNA", "ligation junction", "RNA–DNA linker ligation")
    sc.labels("top")
    return [
        ("Capture the nascent 3′ end", sc.rows(), "A 5′ pre-adenylated DNA linker with N6 and a blocked 3′ ddC is ligated to the RNA 3′-OH. ** marks the ligation."),
        ("Reverse-transcribe and circularize", [Row(chunks=[("linker-primed cDNA → gel select 85–160 nt → CircLigase closes 5′-phosphate to 3′-OH", None, False)]), *circle_rows(CIRCLE.linear, "CircLigase closure")], "Circularization juxtaposes the primer-derived arms so one PCR can recover the RNA-derived insert."),
        ("Deplete and amplify", [Row(chunks=[("biotinylated depletion oligos remove abundant mature-RNA cDNAs → 6–12 cycle indexed PCR", None, False)])], "The minimal productive PCR cycle count is used; the informative product is about 150 nt."),
    ]
