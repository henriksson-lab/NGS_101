"""Conventional ligation-first MethylC-seq / WGBS model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment

TITLE = "MethylC-seq — conventional whole-genome bisulfite sequencing"
NOTES = "01_methylc-seq.html"
SOURCE = 'Protocol source: <a href="https://doi.org/10.1038/nprot.2014.114">Urich et al., <i>Nature Protocols</i> (2015)</a>.'
SUMMARY = "Fragment genomic DNA, repair and dA-tail it, ligate fully methylated adapters, then denature and bisulfite-convert before limited uracil-tolerant PCR."
CAVEAT = ("The paper specifies NEXTflex Bisulfite-Seq Barcodes–12 (Bioo Scientific, "
          "cat. 511912), not an unspecified adapter. It is a single-index, 6-nt i7, "
          "TruSeq-compatible library: standard Read 1, Index 1 and Read 2 sequencing "
          "primers bind; there is no i5 index.")

# The 12 six-base indices supplied with the product family.  They are kept as data rather
# than collapsed to an 8-base generic barcode so the library constructor enforces the
# actual single-index geometry used by this protocol.
NEXTFLEX_I7 = (
    "CGATGT", "TGACCA", "ACAGTG", "GCCAAT", "CAGATC", "CTTGTA",
    "ATCACG", "TTAGGC", "ACTTGA", "GATCAG", "TAGCTT", "GGCTAC",
)


def nextflex_library(index: str = NEXTFLEX_I7[0], *,
                     insert_label: str = "C→T-converted genomic insert") -> tuple[Construct, tuple]:
    """Exact NEXTflex/TruSeq-compatible single-index endpoint used by MethylC-seq.

    P5 and the Read 1 site share their terminal/initial ACAC.  Modelling that overlap
    explicitly avoids inventing four bases, while the dA junction supplies the leading A
    of the reverse-complemented Read 2 site on the opposite adapter arm.
    """
    if index not in NEXTFLEX_I7:
        raise ValueError("index is not one of the 12 NEXTflex Bisulfite-Seq barcodes")
    lib = Construct([
        Segment("P5 before shared ACAC", il.P5[:-4], "p5"),
        Segment("Read 1 site (shares ACAC with P5)", il.TRUSEQ_READ1, "r1"),
        Segment(insert_label, "X" * 36, placeholder=True),
        Segment("dA junction / Read 2 site start", "A"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("6-nt i7", index, "cbc"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="NEXTflex MethylC-seq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid NEXTflex MethylC-seq endpoint: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = nextflex_library()
FINAL_CAPTION = ("PCR-restored NEXTflex library: a single 6-nt i7 and the standard "
                 "TruSeq-compatible sequencing-primer sites. Unmethylated genomic C is "
                 "read as T; protected adapter 5mC remains C through bisulfite treatment.")
SEQUENCING_INTRO = ("Standard TruSeq Read 1, Index 1 and Read 2 primers bind to the "
                    "NEXTflex library. There is no Index 2 read. Paired reads traverse "
                    "the converted insert, whose C/T outcome is compared with the "
                    "unconverted reference sequence.")


def ligated_nextflex_scene() -> Scene:
    ligated, _ = nextflex_library(insert_label="unconverted genomic insert")
    sc = Scene.duplex(list(ligated), label="NEXTflex adapter-ligated molecule")
    sc.junction("top", "Read 1 site (shares ACAC with P5)",
                "unconverted genomic insert", "adapter ligation")
    sc.junction("top", "dA junction / Read 2 site start",
                "Index 1 / Read 2 arm", "adapter ligation")
    sc.labels("top")
    return sc

def sections():
    return [
        ("Fragment, end-repair and dA-tail", [Row(chunks=[("genomic DNA → ~200-bp fragments → blunt 5′-phosphorylated ends → 3′ dA", None, False)])],
         "Mechanical fragmentation is followed by conventional end preparation."),
        ("Ligate protected NEXTflex adapters before conversion", ligated_nextflex_scene().rows(),
         "NEXTflex Bisulfite-Seq Barcodes–12 is a methylated, single-i7, "
         "TruSeq-compatible Y-adapter; ** marks both ligations."),
        ("Convert and amplify", [Row(chunks=[("unmethylated C → U → T in PCR     5mC → C", "w1", False)]),
                                  Row(chunks=[("adapter-ligated duplex → denatured non-complementary single strands → PCR duplex", None, False)])],
         "Bisulfite conversion occurs after adapter ligation and can break already tagged molecules."),
    ]
