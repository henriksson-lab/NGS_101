"""Chen et al. single-cell CAPS+ selective 5hmC library.

Scope of the scTAPS import: the joint Genome Biology paper prints one set of
transposome oligos (Additional file 5 Table S3), one front end ("scTAPS and
scCAP + library construction and sequencing") and one set of custom sequencing
primers (Table S4) for both protocols.  Only those shared declarations and the
front-/back-end graph builders are imported; the conversion chemistry is local.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "sctaps__10.1101+2024.02.12.579864" / "tools"))
from base_conversion import BasePath
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
import sctaps as shared
from sctaps import (INITIAL_NAME, INITIAL_ROWS, READ_LENGTHS, READ_USES, SEQ_PRIMERS,  # noqa: F401
                    SHOW_READ_CONTENT, SOURCE, base_scene, final_library)


def caps_plus_path(base: str) -> BasePath:
    """CAPS+ selective 5hmC oxidation, borane reduction and PCR readout."""
    if base == "5hmC":
        return BasePath(("5hmC", "5fC", "5caC", "DHU", "T"))
    if base in {"5mC", "C"}:
        return BasePath((base, base))
    raise ValueError("CAPS+ input must be C, 5mC or 5hmC")

TITLE = "scCAPS+ — single-cell 5hmC sequencing"
NOTES = "01_sccaps-plus.html"
SUMMARY = "Barcoded Tn5 fragments each single-cell genome before pooling. ACT+BF4− and Pinnick oxidation selectively take 5hmC through 5fC to 5caC; pyridine–borane then makes only 5hmC read as T."
CAVEAT = shared.CAVEAT
SOURCE_OLIGO_EXCEPTIONS = shared.SOURCE_OLIGO_EXCEPTIONS
CONVERTS = ("5hmC",)
CAPS_PATHS = {base: caps_plus_path(base) for base in ("5hmC", "5mC", "C")}

FINAL_LIBRARY = final_library("CAPS+-converted genomic insert", "scCAPS+ library")
FINAL_CAPTION = ("Table S3 transposon arms flank the insert; Nextera XT PCR adds P5/P7 and sample indices (INFERRED, kit sequences not printed). "
                 "In the reads, 5hmC—but not 5mC or C—appears as T.")
SEQUENCING_INTRO = shared.SEQUENCING_INTRO


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    pooled = shared.front_end(g, "200 ng 2-kb filler DNA plus 0.1% spike-ins (methylated lambda, 2-kb unmodified, 5hmC spike-in)")
    act = g.input("act_bf4", "ACT+BF4− in sodium phosphate buffer, 25 µl, 25 °C 16 h; 1.8× AMPure XP", MaterialKind.CONDITION)
    formyl = g.product("formyl", "pooled DNA after selective 5hmC-to-5fC oxidation", MaterialKind.MOLECULAR,
                       state=base_scene(CAPS_PATHS, "ACT+BF4− oxidation", 1, "after ACT+BF4− oxidation"))
    g.react("oxidize_formyl", "Oxidize 5hmC to 5fC with ACT+BF4−", StepKind.CHEMICAL_CONVERSION,
            ((pooled, InputRole.TEMPLATE), (act, InputRole.CONDITION)), (formyl,))
    pinnick = g.input("pinnick", "Pinnick oxidation: NaClO2 and 2-methyl-2-butene in sodium acetate buffer, 30 µl, 25 °C 16 h",
                      MaterialKind.CONDITION)
    carboxyl = g.product("carboxyl", "pooled DNA after 5fC-to-5caC oxidation", MaterialKind.MOLECULAR,
                         state=base_scene(CAPS_PATHS, "Pinnick oxidation with NaClO2", 2, "after Pinnick oxidation"))
    g.react("oxidize_carboxyl", "Pinnick-oxidize 5fC to 5caC", StepKind.CHEMICAL_CONVERSION,
            ((formyl, InputRole.TEMPLATE), (pinnick, InputRole.CONDITION)), (carboxyl,))
    borane = g.input("pyridine_borane", "pyridine–borane reduction; Zymo column with Oligo Binding Buffer", MaterialKind.CONDITION)
    reduced = g.product("reduced", "DNA in which original 5hmC has become DHU", MaterialKind.MOLECULAR,
                        state=base_scene(CAPS_PATHS, "pyridine–borane reduction", 3, "after borane reduction"))
    g.react("reduce", "Reduce 5caC to DHU with pyridine–borane", StepKind.CHEMICAL_CONVERSION,
            ((carboxyl, InputRole.TEMPLATE), (borane, InputRole.CONDITION)), (reduced,))
    return shared.library_back_end(g, reduced, FINAL_LIBRARY,
                                   ["act_bf4", "oxidize_formyl", "formyl", "pinnick", "oxidize_carboxyl",
                                    "carboxyl", "pyridine_borane", "reduce", "reduced"])
