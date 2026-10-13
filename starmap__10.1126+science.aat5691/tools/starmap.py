"""Source-backed molecular model of Wang et al. STARmap."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, circle_rows
from circular import CircularizationGeometry, circularize_ssdna
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind

TITLE = "STARmap — SNAIL amplification and SEDAL in-situ sequencing"
SOURCE = "Wang et al., Science (2018), doi:10.1126/science.aat5691"

def seg(name, seq, tag=None, **kw): return Segment(name, seq, tag, **kw)

TARGET_LEFT = "N" * 22
TARGET_RIGHT = "N" * 22
PAIRING_ARM = "N" * 6
GENE_CODE = "B" * 5
DOCK = "N" * 11

def snail_scene() -> Scene:
    sc = Scene()
    sc.strand("RNA", [seg("left RNA target", TARGET_LEFT, placeholder=True),
                      seg("0–2-nt target gap", "N" * 2, placeholder=True),
                      seg("right RNA target", TARGET_RIGHT, placeholder=True)],
              label="target RNA")
    sc.anneal("primer probe", [seg("primer target arm", TARGET_LEFT, placeholder=True),
                               seg("primer 6-nt assembly arm", PAIRING_ARM, placeholder=True)],
              to="RNA", pair=("primer target arm", "left RNA target"),
              unpaired=("primer 6-nt assembly arm",), label="SNAIL primer probe")
    sc.anneal("padlock probe", [seg("padlock target arm", TARGET_RIGHT, placeholder=True),
                                seg("SEDAL docking region", DOCK, placeholder=True),
                                seg("5-base gene code", GENE_CODE, "cbc", placeholder=True),
                                seg("padlock 6-nt assembly arm", PAIRING_ARM, placeholder=True)],
              to="RNA", pair=("padlock target arm", "right RNA target"),
              unpaired=("SEDAL docking region", "5-base gene code",
                        "padlock 6-nt assembly arm"), label="SNAIL padlock probe",
              mod5="p")
    sc.arrow("primer probe", "adjacent probes assemble; T4 ligase closes one DNA circle")
    return sc

def circle_template() -> Construct:
    return Construct([seg("primer-probe copy", "N" * 28, placeholder=True),
                      seg("SEDAL docking region", DOCK, placeholder=True),
                      seg("5-base gene code", GENE_CODE, "cbc", placeholder=True),
                      seg("padlock-probe copy", "N" * 28, placeholder=True)],
                     "closed SNAIL circle")

CIRCLE = circularize_ssdna(circle_template(), five_prime_phosphate=True)

def rolony() -> Construct:
    return Construct([seg(f"SNAIL circle repeat {i}", CIRCLE.sequence, placeholder=True)
                      for i in range(1, 4)], "three repeats of one STARmap rolony")

def sedal_scene() -> Scene:
    sc = Scene(); sc.strand("rolony", [seg("reading-probe docking site", "N" * 11,
                                                     placeholder=True),
                                      seg("gene-code interrogation base", "B", "cbc",
                                          placeholder=True),
                                      seg("remaining gene-code bases", "B" * 4, "cbc",
                                          placeholder=True)], label="one rolony repeat")
    sc.anneal("reading probe", [seg("11-nt partly degenerate reading probe", "N" * 11,
                                    placeholder=True)], to="rolony",
              pair=("11-nt partly degenerate reading probe", "reading-probe docking site"),
              label="5′-phosphorylated SEDAL reading probe", mod5="p")
    sc.strand("decoder", [seg("fluorescent decoding probe", "N" * 8,
                              placeholder=True)], label="one of four fluorescent decoders")
    sc.arrow("reading probe", "T4 ligase retains the matched decoder; image, strip, shift")
    return sc

def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("fixed_rna", "RNA in fixed, permeabilized tissue", MaterialKind.BIOLOGICAL,
                  state=snail_scene(), primary=True)
    probes = g.input("snail_probes", "paired target-specific SNAIL probes",
                     MaterialKind.OLIGO)
    ligase = g.input("t4_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    circle = g.product("snail_circle", "target-coded DNA circle assembled on RNA",
                       MaterialKind.MOLECULAR,
                       state=tuple(circle_rows(circle_template(), "one T4-ligase junction")))
    g.react("snail_ligation", "Adjacent SNAIL probes assemble and close into a circle",
            StepKind.CIRCULARIZATION,
            ((rna, InputRole.TEMPLATE), (probes, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (circle,),
            geometry=CircularizationGeometry.known(
                "splinted single-strand ligation",
                statement="Adjacent RNA-bound primer and padlock probes create one T4-ligase junction."))
    phi29 = g.input("phi29", "phi29 DNA polymerase", MaterialKind.CATALYST)
    aadutp = g.input("aminoallyl_dutp", "5-(3-aminoallyl)-dUTP", MaterialKind.MOLECULAR)
    product = g.product("rolony", "aminoallyl-bearing SNAIL rolling-circle product",
                        MaterialKind.MOLECULAR, state=rolony())
    g.react("rca", "Rolling-circle amplification in place", StepKind.EXTENSION,
            ((circle, InputRole.TEMPLATE), (phi29, InputRole.CATALYST),
             (aadutp, InputRole.STRUCTURAL)), (product,))
    gel = g.input("hydrogel", "acrylamide/bis-acrylamide hydrogel and AA-NHS anchoring chemistry",
                  MaterialKind.CONDITION)
    anchored = g.product("anchored_rolony", "hydrogel-anchored STARmap rolony",
                         MaterialKind.MOLECULAR, state=rolony())
    g.react("anchor", "Covalently anchor rolonies in hydrogel and clear proteins",
            StepKind.CROSS_LINKING,
            ((product, InputRole.TEMPLATE), (gel, InputRole.CONDITION)), (anchored,))
    protease = g.input("proteinase_k", "proteinase K", MaterialKind.CATALYST)
    library = g.product("cleared_rolony", "cleared, hydrogel-anchored STARmap rolony",
                        MaterialKind.LIBRARY, state=rolony())
    waste = g.product("digested_proteins", "digested tissue proteins", MaterialKind.WASTE)
    g.react("clear_tissue", "Digest proteins while retaining covalently anchored rolonies",
            StepKind.DIGESTION,
            ((anchored, InputRole.TEMPLATE), (protease, InputRole.CATALYST)),
            (library, waste))
    g.discard(waste, "proteinase-K digestion products are washed from the hydrogel")
    g.sequence(library, SequencingPlan(
        acquired_roles=("Spatial gene code",),
        limitation="SEDAL is six cycles of ligation imaging in situ, not R1/I1/I2/R2 sequencing."))
    g.declare_reviewed_scope(("fixed_rna", "snail_probes", "t4_ligase", "snail_ligation",
                              "phi29", "aminoallyl_dutp", "rca", "hydrogel", "anchor",
                              "proteinase_k", "clear_tissue", "cleared_rolony"), SOURCE)
    return g.require_publishable()
