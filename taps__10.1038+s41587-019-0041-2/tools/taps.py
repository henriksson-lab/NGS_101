"""Liu et al. 2019 TAPS chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from base_conversion import taps_path
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature
from endprep import (ThreePrimeOverhangAdapter, cohesive_ligation_geometry,
                     dA_tailed_scene, ligation_product_scene, repair_and_dA_tail)

TITLE = "TAPS — bisulfite-free methylation sequencing"
NOTES = "01_taps.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41587-019-0041-2">Liu et al., <i>Nature Biotechnology</i> (2019)</a>.'
SUMMARY = "TET oxidation and pyridine–borane reduction make 5mC and 5hmC read as T, while unmodified C remains C."
CAVEAT = ""
INDEX6 = "GCCAAT"  # TruSeq Index 6, printed in the defining supplement
TRUSEQ_ADAPTER = ThreePrimeOverhangAdapter(
    "TAPS TruSeq adapter ligating stem", "GCTCTTCCGATCT", "GATCGGAAGAGC")


def library():
    lib = Construct([
        Segment("P5 before shared ACAC", il.P5[:-4], "p5"),
        Segment("Read 1 site", il.TRUSEQ_READ1, "r1"),
        Segment("TAPS-treated genomic insert", "X" * 36, placeholder=True),
        Segment("dA junction / Read 2 site start", "A"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 reverse complement", INDEX6, "cbc",
                feature=feature("sample_i7", "sample_index", "whitelist", whitelist="TruSeq index set")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="TAPS TruSeq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems: raise ValueError("invalid TAPS library: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = "Single-index TruSeq library after TAPS treatment and PCR."
SEQUENCING_INTRO = "Paired reads interrogate the treated insert; i7 identifies the library."


def workflow():
    insert = Construct([Segment("TAPS genomic insert", "X" * 36, placeholder=True)])
    tailed = repair_and_dA_tail(insert)
    w = Workflow(MolecularState("fragmented genomic DNA", tuple(Scene.duplex(list(insert)).rows())))
    w.react("End-repair and add 3′ dA", dA_tailed_scene(tailed).rows(), name="dA-tailed DNA",
            note="The single-base insert overhang is the substrate for the printed TruSeq adapter.")
    geometry = cohesive_ligation_geometry(tailed, TRUSEQ_ADAPTER,
        statement="The TAPS supplement prints the complete TruSeq construct; its adapter stem leaves a 3′-T complementary to insert 3′-dA.")
    lig = ligation_product_scene(FINAL_LIBRARY, "TAPS-treated genomic insert", label="adapter-ligated molecule")
    w.ligate_adapter("Ligate TruSeq adapters", geometry, lig.rows(), name="adapter-ligated TAPS molecule",
                     note="The supplement prints the complete single-index construct.")
    w.react("Oxidize and reduce modified cytosines", [
            Row(chunks=[(taps_path("5mC").text(), "w1", False)]),
            Row(chunks=[(taps_path("5hmC").text(), "w1", False)]),
            Row(chunks=[(taps_path("C").text(), None, False)]),
        ], name="TAPS-converted library", note="TET oxidation makes 5caC; pyridine–borane makes DHU, which PCR copies as T.")
    w.react("PCR", [Row(chunks=[("DHU → T in amplified library; C remains C", None, False)])],
            name="amplified TAPS library", note="The conversion polarity is opposite to bisulfite sequencing.")
    return w
