"""Source-backed molecular model of Wang et al. STARmap."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, circle_rows, feature
from circular import CircularizationGeometry, circularize_ssdna
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from topology import PrimerTemplatedPadlockTopology

TITLE = "STARmap — SNAIL amplification and SEDAL in-situ sequencing"
SOURCE = "Wang et al., Science (2018), doi:10.1126/science.aat5691"

def seg(name, seq, tag=None, **kw): return Segment(name, seq, tag, **kw)

TARGET_LEFT = "N" * 22
TARGET_RIGHT = "N" * 22
# The primer carries one 12-nt tail. The two 6-nt padlock ends hybridize to its two
# halves, making that tail both the ligation splint and the retained Phi29 primer.
# Exact panel-specific bases were not recovered from the defining-paper table.
PRIMER_SPLINT = "N" * 12
PADLOCK_ARM_5 = "N" * 6
PADLOCK_ARM_3 = "N" * 6
GENE_CODE = "B" * 5
DOCK = "N" * 11
# Main text (NSF manuscript text line 128): a 5-base gene-unique barcode from a 1,024-code
# library is built into each padlock; the per-gene codes are in original Table S1.
GENE_CODE_FEATURE = feature("snail_gene_code", "feature_barcode", "whitelist",
                            whitelist="original Table S1 SNAIL 5-base gene codes (not recovered)")

def _target_rna(sc: Scene) -> None:
    sc.strand("RNA", [seg("left RNA target", TARGET_LEFT, placeholder=True),
                      seg("0–2-nt target gap", "N" * 2, placeholder=True),
                      seg("right RNA target", TARGET_RIGHT, placeholder=True)],
              label="target RNA")


def target_rna_scene() -> Scene:
    sc = Scene(); _target_rna(sc); return sc


def snail_scene() -> Scene:
    sc = Scene(); _target_rna(sc)
    sc.anneal("primer probe", [seg("primer target arm", TARGET_LEFT, placeholder=True),
                               seg("primer 12-nt splint / RCA primer", PRIMER_SPLINT,
                                   placeholder=True)],
              to="RNA", pair=("primer target arm", "left RNA target"),
              unpaired=("primer 12-nt splint / RCA primer",), label="SNAIL primer probe")
    sc.anneal("padlock probe", [seg("padlock 5-prime 6-nt assembly arm", PADLOCK_ARM_5,
                                    placeholder=True),
                                seg("padlock target arm", TARGET_RIGHT, placeholder=True),
                                seg("SEDAL docking region", DOCK, placeholder=True),
                                seg("5-base gene code", GENE_CODE, "cbc", placeholder=True, feature=GENE_CODE_FEATURE),
                                seg("padlock 3-prime 6-nt assembly arm", PADLOCK_ARM_3,
                                    placeholder=True)],
              to="RNA", pair=("padlock target arm", "right RNA target"),
              unpaired=("padlock 5-prime 6-nt assembly arm", "SEDAL docking region",
                        "5-base gene code", "padlock 3-prime 6-nt assembly arm"),
              label="SNAIL padlock probe",
              mod5="p")
    sc.arrow("primer probe", "padlock ends dock on the primer tail; T4 ligase closes one junction")
    return sc

def circle_template() -> Construct:
    return Construct([seg("padlock 5-prime assembly arm", PADLOCK_ARM_5, placeholder=True),
                      seg("padlock target arm", TARGET_RIGHT, placeholder=True),
                      seg("SEDAL docking region", DOCK, placeholder=True),
                      seg("5-base gene code", GENE_CODE, "cbc", placeholder=True, feature=GENE_CODE_FEATURE),
                      seg("padlock 3-prime assembly arm", PADLOCK_ARM_3,
                          placeholder=True)],
                     "closed SNAIL circle")

CIRCLE = circularize_ssdna(circle_template(), five_prime_phosphate=True)
SNAIL_TOPOLOGY = PrimerTemplatedPadlockTopology(
    "STARmap SNAIL primer-templated padlock complex",
    TARGET_LEFT,
    TARGET_RIGHT,
    "N" * 2,
    PRIMER_SPLINT,
    PADLOCK_ARM_5,
    DOCK + "N" * len(GENE_CODE),
    PADLOCK_ARM_3,
)

def rolony() -> Construct:
    # Phi29 extends the RNA-bound primer probe around the circle, so each repeat is the
    # complement of the circle (the reading-probe docking site and gene code included).
    return Construct([seg("SNAIL primer probe (RCA primer)", "N" * 6, placeholder=True)] +
                     [seg(f"circle-complement repeat {i}", "N" * len(CIRCLE), placeholder=True)
                      for i in range(1, 4)],
                     "STARmap rolony: primer-probe-initiated tandem circle complements")

def sedal_scene() -> Scene:
    sc = Scene(); sc.strand("rolony", [seg("reading-probe docking site", "N" * 11,
                                                     placeholder=True),
                                      seg("5-base gene code (complement)", "B" * 5, "cbc",
                                          placeholder=True, feature=GENE_CODE_FEATURE)],
              label="one rolony repeat")
    sc.anneal("reading probe", [seg("11-nt partly degenerate reading probe", "N" * 11,
                                    placeholder=True)], to="rolony",
              pair=("11-nt partly degenerate reading probe", "reading-probe docking site"),
              label="SEDAL reading probe (5′-phosphate, 3′ inverted dT)", mod5="p", mod3="InvT")
    sc.strand("decoder", [seg("8-nt fluorescent decoding probe", "N" * 8,
                              placeholder=True)],
              label="one of sixteen 2-base fluorescent decoders (four 1-base decoders in the 28-gene large-volume variant)",
              mod3="OH")
    sc.arrow("reading probe", "T4 ligase joins the matched decoder 3′-OH to the reading-probe 5′-phosphate; image, strip with 60% formamide, shift one base")
    return sc

def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("fixed_rna", "RNA in fixed, permeabilized tissue", MaterialKind.BIOLOGICAL,
                  state=target_rna_scene(), primary=True)
    # NSF manuscript supplement, "STARmap procedure for cell culture and thin tissue
    # sections", Library construction (text lines 1147-1168).
    probes = g.input("snail_probes", "paired target-specific SNAIL probes (100 nM each; heated 90 °C, cooled)",
                     MaterialKind.OLIGO)
    hyb = g.input("hybridization_buffer", "2X SSC, 10% formamide, 1% Tween-20, RVC, salmon sperm DNA; 40 °C overnight",
                  MaterialKind.CONDITION)
    hybridized = g.product("hybridized_probes", "SNAIL primer and padlock probes hybridized to adjacent RNA",
                           MaterialKind.MOLECULAR, state=snail_scene())
    unbound = g.product("unbound_probes", "unhybridized SNAIL probes", MaterialKind.WASTE)
    g.react("hybridize", "Hybridize SNAIL probe pairs to cellular RNA", StepKind.HYBRIDIZATION,
            ((rna, InputRole.TEMPLATE), (probes, InputRole.STRUCTURAL), (hyb, InputRole.CONDITION)),
            (hybridized, unbound))
    g.discard(unbound, "PBSTR and 4X SSC washes remove unhybridized probes")
    ligase = g.input("t4_ligase", "T4 DNA ligase (Thermo EL0011), 1:50, RT 2 h", MaterialKind.CATALYST)
    circle = g.product("snail_circle", "target-coded DNA circle assembled on RNA",
                       MaterialKind.MOLECULAR,
                       state=tuple(circle_rows(circle_template(), "one T4-ligase junction")))
    g.react("snail_ligation", "Adjacent SNAIL probes assemble and close into a circle",
            StepKind.CIRCULARIZATION,
            ((hybridized, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (circle,),
            geometry=CircularizationGeometry.known(
                "splinted single-strand ligation",
                statement="Adjacent RNA-bound primer and padlock probes create one T4-ligase junction."))
    phi29 = g.input("phi29", "Phi29 DNA polymerase (Thermo EP0094), 1:50, 30 °C 2 h", MaterialKind.CATALYST)
    dntp = g.input("dntps", "250 µM dNTPs", MaterialKind.MOLECULAR)
    aadutp = g.input("aminoallyl_dutp", "20 µM 5-(3-aminoallyl)-dUTP", MaterialKind.MOLECULAR)
    product = g.product("rolony", "aminoallyl-bearing SNAIL rolling-circle product",
                        MaterialKind.MOLECULAR, state=rolony())
    g.react("rca", "Rolling-circle amplification primed by the primer probe", StepKind.EXTENSION,
            ((circle, InputRole.TEMPLATE), (phi29, InputRole.CATALYST),
             (dntp, InputRole.STRUCTURAL), (aadutp, InputRole.STRUCTURAL)), (product,))
    aanhs = g.input("aa_nhs", "20 mM acrylic acid NHS ester in PBST, RT 2 h", MaterialKind.CONDITION)
    acrylated = g.product("acrylated_rolony", "acrylamide-functionalized rolony",
                          MaterialKind.MOLECULAR, state=rolony())
    g.react("acrylate", "Convert aminoallyl groups to polymerizable acrylamide moieties",
            StepKind.CHEMICAL_CONVERSION,
            ((product, InputRole.TEMPLATE), (aanhs, InputRole.CONDITION)), (acrylated,))
    gel = g.input("hydrogel", "4% acrylamide, 0.2% bis-acrylamide, 2X SSC; APS/TEMED, RT 1 h",
                  MaterialKind.CONDITION)
    anchored = g.product("anchored_rolony", "hydrogel-anchored STARmap rolony",
                         MaterialKind.MOLECULAR, state=rolony())
    g.react("anchor", "Copolymerize rolonies into the tissue hydrogel",
            StepKind.CROSS_LINKING,
            ((acrylated, InputRole.TEMPLATE), (gel, InputRole.CONDITION)), (anchored,))
    protease = g.input("proteinase_k", "Proteinase K 0.2 mg/ml, 37 °C 1 h", MaterialKind.CATALYST)
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
    g.declare_reviewed_scope(("fixed_rna", "snail_probes", "hybridization_buffer", "hybridize",
                              "unbound_probes", "t4_ligase", "snail_ligation",
                              "phi29", "dntps", "aminoallyl_dutp", "rca", "aa_nhs", "acrylate",
                              "hydrogel", "anchor", "proteinase_k", "clear_tissue",
                              "cleared_rolony"), SOURCE)
    return g.require_publishable()
