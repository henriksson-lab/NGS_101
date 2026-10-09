"""Liu et al. 2019 TAPS chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from base_conversion import taps_path
from chemdraw import Construct, Row, Scene, Segment

TITLE = "TAPS — bisulfite-free methylation sequencing"
NOTES = "01_taps.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41587-019-0041-2">Liu et al., <i>Nature Biotechnology</i> (2019)</a>.'
SUMMARY = "TET oxidation and pyridine–borane reduction make 5mC and 5hmC read as T, while unmodified C remains C."
CAVEAT = ""
INDEX6 = "GCCAAT"  # TruSeq Index 6, printed in the defining supplement


def library():
    lib = Construct([
        Segment("P5 before shared ACAC", il.P5[:-4], "p5"),
        Segment("Read 1 site", il.TRUSEQ_READ1, "r1"),
        Segment("TAPS-treated genomic insert", "X" * 36, placeholder=True),
        Segment("dA junction / Read 2 site start", "A"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 reverse complement", INDEX6, "cbc"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="TAPS TruSeq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems: raise ValueError("invalid TAPS library: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = "Single-index TruSeq library after TAPS treatment and PCR."
SEQUENCING_INTRO = "Paired reads interrogate the treated insert; i7 identifies the library."


def sections():
    lig = Scene.duplex(list(FINAL_LIBRARY), label="adapter-ligated molecule")
    lig.junction("top", "Read 1 site", "TAPS-treated genomic insert", "adapter ligation")
    lig.junction("top", "dA junction / Read 2 site start", "Index 1 / Read 2 arm", "adapter ligation")
    return [
        ("Ligate TruSeq adapters", lig.rows(), "The supplement prints the complete single-index construct; ** marks ligation."),
        ("Oxidize and reduce modified cytosines", [
            Row(chunks=[(taps_path("5mC").text(), "w1", False)]),
            Row(chunks=[(taps_path("5hmC").text(), "w1", False)]),
            Row(chunks=[(taps_path("C").text(), None, False)]),
        ], "TET oxidation makes 5caC; pyridine–borane makes DHU, which PCR copies as T."),
        ("PCR", [Row(chunks=[("DHU → T in amplified library; C remains C", None, False)])],
         "The conversion polarity is opposite to bisulfite sequencing."),
    ]
