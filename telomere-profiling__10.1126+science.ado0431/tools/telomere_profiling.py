"""Telomere Profiling enrichment workflow (Karimian et al. 2024).

Provenance (round 2026-10-17-deferred-followup-03): reaction details and oligos come
from the published Science Supplementary Materials (science.ado0431_sm.pdf, printed
doi:10.1126/science.ado0431): Materials and Methods pp. 3-4 (.pdf.txt lines 77-147),
fig. S1A (lines 262-281; NB50 TeloTag with its six splints) and Table S3 (p. 20, lines
2296-2357; fourteen TeloTags).  The bioRxiv v2 preprint used in earlier rounds is
superseded and no longer supplies any modelled value.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import (Construct, Row, Scene, Segment, complement_segments, feature,
                      feature_rows, oligo, revcomp, strand_row)
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
from telomere import AffinityCapture, HUMAN_C_REPEAT, TelomereEnd

TITLE = "Telomere Profiling"
SOURCE = 'Karimian et al. 2024, <a href="https://doi.org/10.1126/science.ado0431">doi:10.1126/science.ado0431</a>.'
PUBLISHED = ('Karimian et al. 2024 Science Supplementary Materials, Materials and Methods '
             '(science.ado0431_sm.pdf pp. 3–4)')
PUBLISHED_TABLE = ('Karimian et al. 2024 Science Supplementary Materials, Table S3 (p. 20) and fig. S1A')
SOURCE_OLIGO_EXCEPTIONS = {
    "telotag_splints": ("Fig. S1A of the published Supplementary Materials (.pdf.txt lines 275–281) prints the six "
                        "NB50 splints written 3′→5′ beneath the TeloTag, so they are not 5′→3′ text matches; "
                        "each reversed row equals one derived NB50 splint (checked by script this round)."),
}
SITE_NAMES = {"GAATTC": "EcoRI", "CGATCG": "PvuI-HF"}

TELOTAG_TOPS = {
    "NB50": "CCCTCCGATAATGGACTTTGGTAACTTCCTGCGTTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB65": "CCCTCCGATATTCTCAGTCTTCCTCCAGACAAGGTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB68": "CCCTCCGATAGAATCTAAGCAAACACGAAGGTGGTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB88": "CCCTCCGATATCTTCTACTACCGATCCGAAGCAGTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB01": "CCCTCCGATACACAAAGACACCGACAACTTTCTTTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB02": "CCCTCCGATAACAGACGACTACAAACGGAATCGATGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB15": "CCCTCCGATAAGGTCTACCTCGCTAACACCACTGTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB19": "CCCTCCGATAGTTCCTCGTGCAGTGTCAAGAGATTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB20": "CCCTCCGATATTGCGTCCTGTTACGAGAACTCATTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB66": "CCCTCCGATACCGATCCTTGTGGCTTCTAACTTCTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB67": "CCCTCCGATAGTTTGTCATACTCGTGTGCTCACCTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB69": "CCCTCCGATATACAGTCCGAGCCTCATGTGATCTTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB70": "CCCTCCGATAACCGAGATCCTACGAATGGAGTGTTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB72": "CCCTCCGATATAGCTGACTGTCTTCCATACCGACTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
}
# Table S3 prints the 14 TeloTags in this order: NB50, NB65, NB68, NB88, NB01, NB02, NB15,
# NB19, NB20, NB66, NB67, NB69, NB70, NB72 (top strands only; no splint table).
TELOMERE_ARMS = tuple(
    (HUMAN_C_REPEAT * 4)[offset:offset + 18]
    for offset in (5, 4, 3, 2, 1, 0)
)


def top_parts(name="NB50"):
    sequence = TELOTAG_TOPS[name]
    return (
        Segment("5-prime TeloTag ligating flank", sequence[:10], "r3"),
        Segment("24-nt TeloTag sample barcode", sequence[10:34], "cbc",
                feature=feature("telotag_sample_index", "sample_index", "whitelist",
                                whitelist="Karimian et al. Extended Data Table 2")),
        Segment("barcode flank", sequence[34:49]),
        Segment("AscI", sequence[49:57], "r3"),
        # NB50/65/68/88 carry EcoRI before PvuI; the other ten Table S3 tags carry
        # PvuI before EcoRI. Label each site by its bases.
        Segment(SITE_NAMES[sequence[57:63]], sequence[57:63], "r3"),
        Segment(SITE_NAMES[sequence[63:69]], sequence[63:69], "r3"),
        Segment("PacI", sequence[69:77], "r3"),
        Segment("paired post-PacI flank", sequence[77:87]),
        Segment("3-prime biotin tail", sequence[87:]),
    )


def splint_family(name="NB50"):
    """The six exact phase splints, derived from their paired TeloTag strand."""
    paired_top = TELOTAG_TOPS[name][:-12]
    return tuple(revcomp(paired_top) + arm for arm in TELOMERE_ARMS)


# Published discrepancies shown on the page (both locations are in the published
# Supplementary Materials or article; the model's choice follows each entry).
PUBLISHED_DISCREPANCIES = (
    "Order: the published Results summary (author manuscript JATS line 122) lists TeloTag ligation, then DNA digestion; "
    "the Supplementary Materials and Methods (p. 3) digest 40 µg gDNA with ClaI and heat-inactivate before the tagging "
    "ligation. The model follows the Methods, which give the step conditions",
    "Units: TeloTag annealing is printed as '5mM' of each splint and '30 mM' of biotinylated adapter, diluted 1:100 and then "
    "used as '0.3μM duplex TeloTag adapter' (Methods p. 3). The printed units are inconsistent and are shown unresolved",
    "Table label: the Methods say the TeloTag and splint sequences are 'listed in table S2', but Table S2 is telomere lengths; "
    "Table S3 lists the TeloTag top strands, and only the NB50 splints are printed (fig. S1A). The model uses Table S3 and "
    "fig. S1A",
)

TELOTAG_CONSTRUCT=Construct(top_parts(), name="NB50 barcoded TeloTag")
CAPTURE=AffinityCapture("biotin","adapter-encoded restriction site")


def oligos():
    top_rows = [
        oligo(name + " TeloTag", top_parts(name), mods="/5Phos/",
              three="/3BioTEG/-3′")
        for name in TELOTAG_TOPS
    ]
    splint_rows = [  # fig. S1A prints these NB50 splints 3′→5′
        oligo(f"NB50.C{i}", (Segment(f"NB50.C{i}", sequence),))
        for i, sequence in enumerate(splint_family(), 1)
    ]
    return (*top_rows, *splint_rows)


def _telotag_scene():
    paired = complement_segments(top_parts()[:-1])
    sc = Scene()
    sc.strand("splint", [*paired, Segment("telomere-complement arm", TELOMERE_ARMS[0], "r2")],
              label="NB50.C1 splint")
    sc.anneal("TeloTag", top_parts(), to="splint",
              pair=("5-prime TeloTag ligating flank", "5-prime TeloTag ligating flank'"),
              unpaired=("3-prime biotin tail",), label="NB50 TeloTag",
              above=True, mod5="p", mod3="BioTEG")
    end = TelomereEnd(overhang_copies=3)
    chromosome = [Segment("subtelomere", "X" * end.subtelomere_bases, placeholder=True),
                  Segment("duplex telomere", end.duplex, "r1"),
                  # The six splints cover the six cyclic telomere phases.  Show
                  # the chromosome phase that actually pairs to NB50.C1 rather
                  # than forcing the canonical repeat origin onto this example.
                  Segment("G-rich overhang", revcomp(TELOMERE_ARMS[0]), "r2")]
    sc.anneal("G-rich strand", chromosome, to="splint",
              pair=("G-rich overhang", "telomere-complement arm"),
              unpaired=("subtelomere", "duplex telomere"), label="chromosome")
    bottom = complement_segments(chromosome[:2])
    sc.anneal("C-rich strand", bottom, to="G-rich strand",
              pair=("duplex telomere'", "duplex telomere"), label="chromosome")
    sc.mark("TeloTag", "24-nt TeloTag sample barcode", "published 24-nt barcode")
    sc.ligation_end("TeloTag", "5-prime TeloTag ligating flank", kind="splinted",
                    polarity="5-prime", owner="adapter",
                    text="splinted terminal TeloTag ligation")
    sc.footer("** HiFi Taq DNA ligase seals the TeloTag 5′ phosphate at the chromosome end",
              "TeloTag", "5-prime TeloTag ligating flank")
    return sc

def _cleanup(g, substrate, ident, label, beads_label, kept_label, waste_label, why):
    beads = g.input(ident + "_beads", beads_label, MaterialKind.SELECTION, source=PUBLISHED)
    kept = g.product(ident + "_eluate", kept_label, MaterialKind.MOLECULAR,
                     state=substrate.state, source=PUBLISHED)
    waste = g.product(ident + "_waste", waste_label, MaterialKind.WASTE, source=PUBLISHED)
    g.react(ident, label, StepKind.CLEANUP,
            ((substrate, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, waste),
            source=PUBLISHED, changes_molecule=False)
    g.discard(waste, why)
    return kept


def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph(TITLE, authority=SOURCE)
    top = g.input("telotag_top", "NB50 TeloTag (/5Phos/ … /3BioTEG/), representative of the 14 Table S3 barcodes",
                  MaterialKind.OLIGO, state=TELOTAG_TOPS["NB50"], source=PUBLISHED_TABLE)
    splints = g.input("telotag_splints", "NB50.C1–C6 phase splints", MaterialKind.OLIGO,
                      state=splint_family("NB50"), source=PUBLISHED_TABLE)
    buffer = g.input("hifi_taq_buffer", "HiFi Taq DNA Ligase Reaction Buffer (NEB M0647S); printed 5mM each splint, 30 mM adapter",
                     MaterialKind.CONDITION, source=PUBLISHED)
    duplex = g.product("telotag_duplex", "annealed duplex TeloTag (diluted 1:100 in 1x Taq buffer)",
                       MaterialKind.ADAPTER, state=_telotag_scene(), source=PUBLISHED)
    g.react("anneal_telotag", "Anneal TeloTag to its six splints (99 °C, −1 °C/min)", StepKind.HYBRIDIZATION,
            ((splints, InputRole.TEMPLATE), (top, InputRole.STRUCTURAL), (buffer, InputRole.CONDITION)),
            (duplex,), source=PUBLISHED)
    gdna = g.input("genomic_dna", "40 µg high-molecular-weight genomic DNA", MaterialKind.BIOLOGICAL,
                   state=molecule("chromosome ends", (role("subtelomere"), role("native telomere terminus"))),
                   primary=True, source=PUBLISHED)
    enzyme = g.input("internal_enzyme", "3 µL ClaI (NEB R0197S)", MaterialKind.CATALYST, source=PUBLISHED)
    digested = g.product("digested", "restriction-shortened DNA with native telomere ends",
                         MaterialKind.MOLECULAR, state=gdna.state, source=PUBLISHED)
    g.react("digest_internal", "Digest genomic DNA (37 °C 2 h; 65 °C 20 min)", StepKind.DIGESTION,
            ((gdna, InputRole.TEMPLATE), (enzyme, InputRole.CATALYST)), (digested,), source=PUBLISHED)
    ligase = g.input("hifi_taq_ligase", "1 µL HiFi Taq DNA Ligase (NEB M0647S) per 50 µL with 4 µg DNA",
                     MaterialKind.CATALYST, source=PUBLISHED)
    tagged = g.product("telotagged", "biotin-TeloTag ligated to the native chromosome terminus",
                       MaterialKind.MOLECULAR, state=digested.state, source=PUBLISHED)
    g.react("ligate_telotag", "Ligate duplex TeloTag (4 µL 0.3 µM; 65 °C 5 min, then 15 × [65 °C 1 min, 45 °C 3 min])",
            StepKind.LIGATION,
            ((digested, InputRole.TEMPLATE), (duplex, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
            (tagged,), source=PUBLISHED,
            geometry=AdapterLigationGeometry.unavailable(
                "Splinted ligation of the TeloTag 5′ phosphate at the chromosome end; the native overhang phase is sample-dependent, so six phase splints are supplied rather than one end geometry."))
    pooled = _cleanup(g, tagged, "spri", "Pool tagging reactions; SPRI select 45 µL : 100 µL, 20 min; two 85% ethanol washes; elute in 65 °C 1x rCutSmart, 20 min",
                      "SPRI beads (Beckman Coulter B23318)", "TeloTagged gDNA",
                      "excess TeloTag adapters and bead supernatant", "removed with SPRI supernatant and 85% ethanol washes")
    beads = g.input("streptavidin", "Dynabeads MyOne Streptavidin C1 (65001), 1 µg : 250 ng DNA, kilobaseBINDER binding buffer",
                    MaterialKind.SELECTION, source=PUBLISHED)
    captured = g.product("captured", "bead-bound biotin-TeloTag telomere DNA", MaterialKind.MOLECULAR,
                         state=pooled.state, source=PUBLISHED)
    unbound = g.product("unbound", "background genomic DNA without biotin-TeloTag", MaterialKind.WASTE, source=PUBLISHED)
    g.react("capture", "Bind to streptavidin (RT 20 min); wash 2× kilobaseBINDER wash, 1× EB, 1× rCutSmart, 1× EB",
            StepKind.AFFINITY_CAPTURE, ((pooled, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
            (captured, unbound), source=PUBLISHED)
    g.discard(unbound, "removed in the binding supernatant and washes")
    release_enzyme = g.input("release_enzyme", "3 µL PvuI (NEB R3150S) in 72 µL 1X rCutSmart",
                             MaterialKind.CATALYST, source=PUBLISHED)
    released = g.product("released", "telomere molecules released at a TeloTag restriction site",
                         MaterialKind.MOLECULAR, state=captured.state, source=PUBLISHED)
    stub = g.product("bead_stub", "bead-bound biotinylated TeloTag remnant", MaterialKind.WASTE, source=PUBLISHED)
    g.react("release", "Release telomeres by PvuI digestion at the TeloTag (37 °C 30 min; 65 °C 20 min)",
            StepKind.DIGESTION, ((captured, InputRole.TEMPLATE), (release_enzyme, InputRole.CATALYST)),
            (released, stub), source=PUBLISHED)
    g.discard(stub, "remains on the magnetically separated streptavidin beads")
    motor = g.input("motor_adapter", "ONT Ligation Sequencing Kit V10 (SQK-LSK110) motor adapter", MaterialKind.ADAPTER,
                    inferred=True, source=PUBLISHED + " The kit's motor-bearing adapter assembly is not disclosed.")
    kit = g.input("motor_ligase", "SQK-LSK110 kit reagents (kit instructions not inspected); libraries eluted in 40 µl EB", MaterialKind.CATALYST,
                  inferred=True, source=PUBLISHED)
    library = g.product("library", "Telomere Profiling nanopore library", MaterialKind.LIBRARY, state=Construct([
        Segment("motor-loaded ONT adapter", "[motor adapter]", placeholder=True, inferred=True),
        Segment("enriched TeloTag–telomere fragment", "[released TeloTag–telomere–subtelomere]", placeholder=True)],
        name="Telomere Profiling nanopore library"))
    g.react("ligate_motor", "Standard ONT ligation library preparation", StepKind.LIGATION,
            ((released, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL), (kit, InputRole.CATALYST)),
            (library,), source=PUBLISHED, geometry=_ont_geometry())
    g.sequence(library, SequencingPlan(acquired_roles=("Sequencing output",),
        limitation="MinION, R9.4.1 FLO-MIN106D flow cells, each library split into three loads of ~16 h with EXP-WSH004 washes; MinKNOW 5.7.5 without live basecalling; Guppy v6 super-high-accuracy basecalling (published Supplementary Materials pp. 4 and 6). Motor entry is platform-driven; no sequencing primer binds the library."))
    g.declare_reviewed_scope(tuple(list(g.materials) + list(g.steps)), SOURCE)
    return g.require_publishable()


def _ont_geometry():
    return AdapterLigationGeometry.unavailable("The Telomere Profiling protocol invokes the historical ONT ligation workflow, but the cited source does not disclose its complete motor-bearing adapter assembly.")

def panels(): return [("Digest high-molecular-weight genomic DNA away from telomere ends",[Row(chunks=[("internal restriction cut | subtelomere -- telomere -- native 3′ overhang",None,False)])],"Restriction reduces non-telomeric length while leaving the native chromosome terminus available for tagging."),("Anneal six telomeric splints and ligate the biotinylated TeloTag",_telotag_scene(),"NB50 and its six phase splints are source-exact. The splint bridges the TeloTag to the native G-rich overhang; ** marks the terminal ligation."),("Select tagged ends on streptavidin and release at the adapter restriction site",[Row(chunks=[("biotin-TeloTag -> streptavidin capture -> wash -> restriction release",None,False)])],"Affinity selection removes most untagged genomic fragments; release retains the tagged telomere molecule."),("ONT adapter terminal geometry",_ont_geometry().rows(),_ont_geometry().statement),("Attach Oxford Nanopore sequencing adapters",[Row(chunks=[("[motor-loaded ONT adapter] ** [enriched TeloTag–telomere fragment]",None,False)])],"The selected material enters the historical ligation library workflow; ** marks adapter ligation.")]
