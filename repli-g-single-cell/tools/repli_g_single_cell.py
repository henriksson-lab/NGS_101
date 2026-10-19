"""Source-owned REPLI-g Single Cell graph: 2012 handbook MDA plus an unresolved library boundary.

Only the October 2012 QIAGEN REPLI-g Single Cell Handbook is on disk and identified.
The cited named-kit library protocol (QIAGEN RG35, August 2014) and the Illumina
Nextera DNA Sample Preparation Guide 15027987 Rev B (October 2012) it cites are not
available; the November 2018 RG35 names the REPLI-g Advanced DNA Single Cell kit and
is not substituted. No canonical Nextera arms, index widths, dual indexing, read
lengths or adapter ends are supplied: the downstream library is a role-only boundary.
"""
from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import seqprimers as sp
from chemdraw import Scene, Segment
from protocol_graph import (Disposition, InputRole as I, MaterialKind as M,
                            ProtocolGraph, SequencingPlan, StepKind as K, Terminal)

HANDBOOK = ("QIAGEN REPLI-g Single Cell Handbook, October 2012 (supplier mirror, md5 6411c5e2): "
            "Introduction/Principle pp. 6–7 (text lines 224–315); Protocol: Amplification of "
            "Genomic DNA from Single Cells pp. 12–14 (lines 477–655)")
AUTHORITY = HANDBOOK
UNAVAILABLE_LIBRARY = (
    "The 2012 handbook (step 13, p. 14) states only that amplified DNA can be used for "
    "next-generation sequencing. The cited named-kit library route, QIAGEN RG35 Nextera "
    "protocol (August 2014), and its cited Illumina Nextera DNA Sample Preparation Guide "
    "15027987 Rev B (October 2012) are unavailable; the November 2018 RG35 (recorded in the "
    "source manifest, not inspected, as naming the REPLI-g Advanced DNA Single Cell kit) "
    "is not a revision substitute. Tagmentation "
    "adapters, PCR primers, index scheme/widths, adapter ends, cleanup and input amount "
    "are therefore unassigned.")
SEQUENCING_ENDING = ("Sequencing of a library prepared from the amplified DNA; platform, "
                     "primer sites, index scheme and read lengths unassigned because the "
                     "cited library-preparation originals are unavailable.")

TITLE = "REPLI-g Single Cell — phi29 MDA; downstream library route unresolved"
NOTES = "01_repli-g-single-cell.html"
SOURCE = ('QIAGEN <a href="https://www.ys-bio.cn/app/wp-content/uploads/2021/12/'
          'EN-REPLI-g-Single-Cell-Handbook.pdf">REPLI-g Single Cell Handbook '
          '(October 2012; supplier mirror)</a>, pp. 6–7 and 12–14. The cited downstream '
          'QIAGEN RG35 Nextera protocol (August 2014) and Illumina Nextera DNA Sample '
          'Preparation Guide 15027987 Rev B (October 2012) were not available.')
SUMMARY = ("A single cell in 4 µl PBS is lysed and its DNA denatured with 3 µl Buffer D2 "
           "(Buffer DLB + DTT) for 10 min at 65 °C, then 3 µl Stop Solution neutralizes the "
           "denaturation. A 40 µl master mix (H2O sc, REPLI-g sc Reaction Buffer, REPLI-g sc "
           "DNA Polymerase) drives phi29 multiple-displacement amplification with "
           "exonuclease-resistant primers for 8 h at 30 °C, stopped at 65 °C for 3 min. An "
           "aliquot of the amplified DNA enters a sequencing-library preparation whose "
           "cited original is unavailable, so that stage is shown as an unresolved boundary.")
CAVEAT = ("The 2012 handbook does not print primer bases, length, protection positions or "
          "which kit component carries the primers, and does not call Buffer D2 alkaline. "
          "It reports primer-dimer DNA in no-template controls, shown as a discarded "
          "parallel control. The downstream library stage is an explicit unresolved "
          "boundary: " + UNAVAILABLE_LIBRARY + " No canonical Nextera sequence, adapter "
          "geometry or read plan is substituted. MDA yields a branched population; it is "
          "not drawn as one linear molecule.")
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
READ_SPECS = (sp.read_spec("Sequencing output", None,
                           "Genomic DNA from the MDA product; read/index acquisition "
                           "unassigned (cited library originals unavailable)", "include"),)


def role(name: str) -> Segment:
    return Segment(name, "[" + name + "]", placeholder=True)


def state(stage: str) -> Scene:
    """Independent molecular populations at each stage; role tokens carry no bases."""
    sc = Scene()
    if stage == "cell":
        sc.strand("genome", [role("cellular genomic DNA")], label="single-cell genome")
    elif stage == "denatured":
        sc.strand("genome", [role("denatured single-stranded genomic DNA")],
                  label="denatured genomic strands in Buffer D2/Stop Solution")
    elif stage in {"amplified", "inactivated", "aliquot", "stored"}:
        sc.strand("copies", [role("long phi29 copies of genomic tracts (>10 kb typical, 2–100 kb)"),
                             role("displaced strands and further primed branches")],
                  label="branched MDA daughter population")
        sc.strand("originals", [role("original genomic template strands")],
                  label="independent original DNA")
        sc.strand("primers", [role("residual exonuclease-resistant primers; identity unassigned")],
                  label="independent residual primer population")
        if stage == "inactivated":
            sc.note("copies", "Heat inactivation stops the polymerase; it does not change DNA identity.")
    elif stage == "control":
        sc.strand("dimer", [role("high-molecular-weight primer-dimer extension DNA")],
                  label="no-template-control product")
    elif stage == "dimer_seed":
        sc.strand("dimer", [role("intermolecular primer pairing; bases unassigned")],
                  label="primer-dimer seed in no-template control")
    elif stage == "library":
        sc.strand("library", [role("sequencing library from MDA DNA; adapters, ends, indexes unassigned")],
                  label="native library population; cited library originals unavailable")
        sc.strand("carried", [role("possible carried MDA DNA, primers and library reagents; recovery unassigned")],
                  label="independent possible carried populations")
    else:
        raise ValueError(stage)
    return sc


def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph("REPLI-g Single Cell", authority=AUTHORITY)
    cell = g.input("cell", "single cell in 4 µl PBS (step 2)", M.BIOLOGICAL,
                   state=state("cell"), primary=True)
    d2 = g.input("buffer_d2", "3 µl Buffer D2 (33 µl reconstituted Buffer DLB + 3 µl 1 M DTT "
                 "per 12 reactions); 10 min at 65 °C", M.CONDITION)
    denatured = g.product("denatured_dna", "lysed, denatured single-cell genomic DNA",
                          M.MOLECULAR, state=state("denatured"))
    g.react("lyse_denature", "Lyse the cell and denature genomic DNA in Buffer D2 "
            "(10 min, 65 °C; steps 3–4)", K.EXTRACTION,
            ((cell, I.TEMPLATE), (d2, I.CONDITION)), (denatured,))
    stop = g.input("stop_solution", "3 µl Stop Solution (neutralization buffer); ice",
                   M.CONDITION)
    neutral = g.product("neutralized_dna", "neutralized denatured genomic DNA (10 µl)",
                        M.MOLECULAR, state=state("denatured"))
    g.react("neutralize", "Stop denaturation with Stop Solution; hold on ice (step 5)",
            K.QUENCHING, ((denatured, I.TEMPLATE), (stop, I.CONDITION)), (neutral,),
            changes_molecule=False)
    primers = g.input("mda_primers", "kit-supplied exonuclease-resistant primers (p. 7); "
                      "carrying component, bases, length and protection pattern not printed",
                      M.OLIGO, inferred=True,
                      source=HANDBOOK + " (primer composition unassigned)")
    buffer = g.input("reaction_buffer", "REPLI-g sc Reaction Buffer 29 µl + H2O sc 9 µl (Table 5)",
                     M.CONDITION)
    phi29 = g.input("phi29", "REPLI-g sc DNA Polymerase, optimized Phi 29 formulation (2 µl)",
                    M.CATALYST)
    amplified = g.product("mda_product", "multiple-displacement amplification reaction (50 µl)",
                          M.MOLECULAR, state=state("amplified"))
    g.react("mda", "Add 40 µl master mix to 10 µl denatured DNA; phi29 multiple-displacement "
            "amplification, 8 h at 30 °C (steps 7–9)", K.EXTENSION,
            ((neutral, I.TEMPLATE), (primers, I.STRUCTURAL), (buffer, I.CONDITION),
             (phi29, I.CATALYST)), (amplified,))
    heat = g.input("heat_inactivation", "65 °C for 3 min (step 10)", M.CONDITION)
    inactivated = g.product("mda_inactivated", "amplified DNA after polymerase inactivation "
                            "(~40 µg per 50 µl typical)", M.MOLECULAR, state=state("inactivated"))
    g.react("inactivate", "Inactivate REPLI-g sc DNA Polymerase (65 °C, 3 min)",
            K.THERMAL_TREATMENT, ((amplified, I.TEMPLATE), (heat, I.CONDITION)), (inactivated,),
            changes_molecule=False)
    aliquot = g.product("library_input", "aliquot of amplified DNA diluted in water or TE; "
                        "amount per unavailable library protocol", M.MOLECULAR,
                        state=state("aliquot"))
    stored = g.product("stored_dna", "amplified DNA not used directly, stored at 4 °C "
                       "short-term or –20 °C long-term (step 11)", M.MOLECULAR, state=state("stored"))
    g.react("allocate", "Store amplified DNA not used directly; dilute the amount required "
            "downstream in water or TE (steps 11–12)", K.SPLIT, ((inactivated, I.TEMPLATE),), (aliquot, stored),
            changes_molecule=False)
    g.mark_terminal(Terminal(stored.id, Disposition.STORED,
                             "Physically separate portion not used directly (step 11); later use unassigned."))

    blank = g.input("no_template_control", "parallel no-template control tube (p. 12 note)",
                    M.CONDITION)
    seed = g.product("primer_dimer_seed", "annealed primer dimer in no-template control",
                     M.MOLECULAR, state=state("dimer_seed"))
    g.react("control_primer_pair", "Primers pair in the no-template control", K.HYBRIDIZATION,
            ((primers, I.TEMPLATE), (primers, I.STRUCTURAL), (blank, I.CONDITION)), (seed,))
    dimer = g.product("primer_dimer_dna", "~40 µg high-molecular-weight primer-dimer DNA in "
                      "no-template control", M.WASTE, state=state("control"))
    g.react("control_primer_dimer", "Random extension of primer dimers in the no-template "
            "control", K.EXTENSION,
            ((seed, I.TEMPLATE), (primers, I.STRUCTURAL), (buffer, I.CONDITION),
             (phi29, I.CATALYST)), (dimer,))
    g.discard(dimer, "Parallel negative-control artifact; not sampled for the sample library.")

    library_prep = g.input("library_preparation", "sequencing-library preparation; cited "
                           "route QIAGEN RG35 (August 2014) unavailable, reagents unassigned",
                           M.CONDITION)
    library = g.product("library", "sequencing library from MDA DNA; native ends unassigned",
                        M.LIBRARY, state=state("library"))
    g.react("library_preparation_boundary", "UNRESOLVED — prepare a sequencing library from "
            "the amplified-DNA aliquot", K.MAGIC,
            ((aliquot, I.TEMPLATE), (library_prep, I.CONDITION)), (library,),
            reason=UNAVAILABLE_LIBRARY)
    g.sequence(library, SequencingPlan(acquisition="unspecified", limitation=SEQUENCING_ENDING))
    g.declare_reviewed_scope((
        "cell", "buffer_d2", "lyse_denature", "denatured_dna", "stop_solution", "neutralize",
        "neutralized_dna", "mda_primers", "reaction_buffer", "phi29", "mda", "mda_product",
        "heat_inactivation", "inactivate", "mda_inactivated", "allocate", "library_input",
        "stored_dna", "no_template_control", "control_primer_pair", "primer_dimer_seed",
        "control_primer_dimer", "primer_dimer_dna", "library_preparation",
        "library_preparation_boundary", "library"), AUTHORITY)
    return g.require_publishable()


def oligos():
    """No kit oligo sequence is printed in the identified on-disk original."""
    return ()
