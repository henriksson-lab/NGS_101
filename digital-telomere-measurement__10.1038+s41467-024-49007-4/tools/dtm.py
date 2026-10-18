"""Digital Telomere Measurement workflow (Sanchez et al. 2024)."""
from __future__ import annotations
import sys
from dataclasses import dataclass
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row, Scene, Segment, feature, oligo
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from telomere import TelomereEnd

SEQ_TETHER="AACCTTGGAGATGCACGGAGCAAGCAAT"
CAPTURE="TGCTCCGTGCATCTCCAAGGTTCACAAAGACACCGACAACTTTCTTCCTAACC"
# Supplementary Data 3 prints the oligos unpartitioned.  The partition is derived
# within that table: all eight canonical capture oligos share their first 22 nt
# (the seqTether-complementary arm) and last 7 nt (CCTAACC); the seven TSQ variants
# share the first 46 nt with p-t3-nb01 and differ only in the last 7 nt.  So the
# barcode is 24 nt and the telomere-targeting arm 7 nt.
CAPTURE_TETHER_ARM=CAPTURE[:22]
CAPTURE_BARCODE=CAPTURE[22:46]
CAPTURE_TELOMERE_ARM=CAPTURE[46:]

def oligos():
    return (
        oligo("seqTether",(Segment("seqTether",SEQ_TETHER),)),
        oligo("p-t3-nb01 capture oligo",(
            Segment("seqTether-complementary arm",CAPTURE_TETHER_ARM),
            Segment("capture barcode",CAPTURE_BARCODE,"cbc", feature=feature("capture_sample", "sample_index", "whitelist", whitelist="Sanchez 2024 Supplementary Data 3")),
            Segment("telomere arm",CAPTURE_TELOMERE_ARM)),mods="/5Phos/"),
    )

def _duplex():
    sc=Scene()
    sc.strand("capture",[Segment("seqTether-complementary arm",CAPTURE_TETHER_ARM),
                          Segment("capture barcode",CAPTURE_BARCODE,"cbc", feature=feature("capture_sample", "sample_index", "whitelist", whitelist="Sanchez 2024 Supplementary Data 3")),
                          Segment("telomere arm",CAPTURE_TELOMERE_ARM)],
              label="p-t3-nb01 capture",mod5="p")
    sc.anneal("tether",[Segment("paired seqTether",SEQ_TETHER[:22]),
                         Segment("seqTether 3-prime tail",SEQ_TETHER[22:])],
              to="capture",pair=("paired seqTether","seqTether-complementary arm"),
              unpaired=("seqTether 3-prime tail",),label="seqTether")
    sc.mark("capture","capture barcode",f"{len(CAPTURE_BARCODE)}-nt capture barcode (variable region of Supplementary Data 3)")
    return sc

def _telomere_capture():
    end=TelomereEnd(repeat="TTAGGG", duplex_copies=3, overhang_copies=2, subtelomere_bases=16)
    sc=end.scene()
    sc.anneal("capture oligo",[Segment("capture tether arm",CAPTURE_TETHER_ARM),
                                Segment("capture barcode",CAPTURE_BARCODE,"cbc", feature=feature("capture_sample", "sample_index", "whitelist", whitelist="Sanchez 2024 Supplementary Data 3")),
                                Segment("CCTAACC telomere arm",CAPTURE_TELOMERE_ARM, inferred=True)],
              to="G-rich strand",pair=("CCTAACC telomere arm","G-rich 3-prime overhang"),
              shift=4,unpaired=("capture tether arm","capture barcode"),
              label="capture oligo",mod5="p")
    sc.anneal("seqTether",[Segment("paired seqTether",SEQ_TETHER[:22]),
                           Segment("seqTether 3-prime tail",SEQ_TETHER[22:])],
              to="capture oligo",pair=("paired seqTether","capture tether arm"),
              unpaired=("seqTether 3-prime tail",),label="seqTether")
    sc.footer("the 5′-phosphorylated capture end is distal (seqTether-paired); its 3′ end faces the C-rich 5′ terminus across a potential gap. Which nick T4 ligase seals is not stated in the Methods.",
              "C-rich strand")
    return sc

def panels():
    return [
        ("Anneal each barcoded capture oligo to seqTether",_duplex(),
         "The representative p-t3-nb01 sequence and seqTether are printed in Supplementary Data 3. The 22-nt tether-complementary arm, 24-nt barcode and 7-nt CCTAACC telomere arm are derived by comparing the printed canonical and variant capture oligos."),
        ("Anneal and ligate the capture duplex at native telomeres",_telomere_capture(),
         "CCTAACC is shown in one representative register on the terminal G-rich overhang; the native terminal phase is not drawn from the source. The Methods ligate the duplex overnight with T4 DNA ligase at 37 °C and heat-inactivate at 65 °C, but do not identify the sealed junction."),
        ("Fill the capture–telomere gap with Sulfolobus DNA polymerase IV",
         [Row(chunks=[("capture/tether ** gap-filled native telomere -- subtelomere",None,False)])],
         "Sulfolobus DNA polymerase IV (56 °C 2 min, 72 °C 15 min) fills potential gaps between the capture oligo and the double-strand/single-strand junction; ProNex beads (1.6×) then recover the post-capture library. The published core protocol does not mechanically fragment this product."),
        ("Ligate the Oxford Nanopore sequencing adapter",
         [Row(chunks=[("[motor-loaded ONT adapter] ** [gap-filled captured chromosome end]",None,False)])],
         "AMII (R9.4.1) or NA (R10.4) adapter with Quick T4 DNA ligase, 1 h at room temperature. The proprietary motor-adapter terminal bases are unavailable; the published capture and tether bases remain explicit."),
    ]

TITLE = "Digital Telomere Measurement (DTM)"
SOURCE = 'Sanchez et al. 2024, <a href="https://doi.org/10.1038/s41467-024-49007-4">doi:10.1038/s41467-024-49007-4</a>.'
METHODS = "Sanchez et al. 2024 PMC11189511, Methods: Telomere capture sequencing nanopore library preparation (lines 158–160)"
TABLE = "Sanchez et al. 2024 Supplementary Data 3, MOESM6, rows seqTether and p-t3-nb01"


def _single(name, sequence, *, phosphate=False):
    sc = Scene()
    parts = ([Segment("seqTether-complementary arm", CAPTURE_TETHER_ARM),
              Segment("capture barcode", CAPTURE_BARCODE, "cbc",
                      feature=feature("capture_sample", "sample_index", "whitelist",
                                      whitelist="Sanchez 2024 Supplementary Data 3")),
              Segment("telomere arm", CAPTURE_TELOMERE_ARM)]
             if name == "p-t3-nb01" else [Segment(name, sequence)])
    sc.strand(name, parts, mod5="p" if phosphate else "")
    return sc


@dataclass(frozen=True)
class CapturedPopulation:
    """An unresolved assembly with unplaced identities and independent residues.

    Unplaced segments are metadata, not extra free strands or an invented bond.
    Residual scenes are possible co-resident physical classes, not observed yields.
    """
    assembly: Scene
    unplaced_segments: tuple[Segment, ...]
    possible_residuals: tuple[Scene, ...]

    def rows(self):
        rows = self.assembly.rows()
        for residual in self.possible_residuals:
            rows.extend([Row(), *residual.rows()])
        return rows


def _possible_residual(name, *, capture=False):
    sc = Scene()
    if capture:
        parts = [Segment("possible free capture tether arm", CAPTURE_TETHER_ARM, inferred=True),
            Segment("possible free capture barcode", CAPTURE_BARCODE, "cbc", inferred=True,
                feature=feature("capture_sample", "sample_index", "whitelist", whitelist="Sanchez 2024 Supplementary Data 3")),
            Segment("possible free capture telomere arm", CAPTURE_TELOMERE_ARM, inferred=True)]
    else:
        parts = [Segment(name, "[" + name + "]", placeholder=True, inferred=True)]
    sc.strand(name, parts, label="INFERRED — possible " + name, mod5="p" if capture else "")
    sc.footer("Candidate residual class only; abundance and cleanup carryover are unassigned, not claimed absent or fully removed.", name)
    return sc


def _product_state(label, *, gap_filled=False, motor=False):
    """Keep unassigned covalent structure, retained identity and mixture scope."""
    sc = Scene()
    sc.strand("captured end", [Segment("captured telomere assembly",
        "[capture/tether–telomere assembly: sealed nick unresolved]", placeholder=True,
        inferred=True)], label=label)
    sc.footer("Source-limited assembly: native G-rich and C-rich DNA, capture and tether; covalent connectivity is unassigned.", "captured end")
    barcode = Segment("retained capture barcode (position unassigned)", CAPTURE_BARCODE, "cbc",
        feature=feature("capture_sample", "sample_index", "whitelist", whitelist="Sanchez 2024 Supplementary Data 3"))
    sc.footer("Retained unplaced identifier metadata: capture_sample (source-exact capture barcode); no junction or coordinate is assigned.", "captured end")
    if gap_filled:
        sc.footer("Pol IV fills potential gaps; tract length and participating strand remain unassigned in this model.", "captured end")
    if motor:
        sc.footer("AMII (R9) or NA (R10) adapter is ligated; historical motor assembly and sealed adapter boundary remain unresolved.", "captured end")
    residuals = [
        _possible_residual("unreacted native genomic DNA"),
        _possible_residual("partially capture-adapted genomic DNA"),
        _possible_residual("free capture oligo", capture=True),
        _possible_residual("free or incompletely incorporated seqTether"),
    ]
    if motor:
        residuals.extend((_possible_residual("captured DNA without motor adapter"),
            _possible_residual("free or incompletely attached motor adapter")))
    return CapturedPopulation(sc, (barcode,), tuple(residuals))


def _cleanup(g, current, ident, ratio, label):
    beads = g.input(ident + "_beads", "Promega ProNex beads (" + ratio + ")", MaterialKind.SELECTION, source=METHODS)
    kept = g.product(ident + "_eluate", label, current.kind, state=current.state, inferred=current.inferred, source=METHODS)
    waste = g.product(ident + "_waste", "supernatant and two 80% ethanol washes", MaterialKind.WASTE, source=METHODS)
    g.react(ident, "ProNex " + ratio + " cleanup; elute in ONT EB (37 °C, 15 min)", StepKind.CLEANUP,
            ((current, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, waste),
            source=METHODS, changes_molecule=False)
    g.discard(waste, "physically removed on the magnet in supernatant and ethanol washes; no barcode/read filtering is implied")
    return kept


def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph(TITLE, authority=SOURCE)
    capture = g.input("capture", "p-t3-nb01 (/5Phos/), representative canonical capture barcode",
        MaterialKind.OLIGO, state=_single("p-t3-nb01", CAPTURE, phosphate=True), source=TABLE)
    tether = g.input("tether", "seqTether", MaterialKind.OLIGO,
        state=_single("seqTether", SEQ_TETHER), source=TABLE)
    low_te = g.input("low_te", "equimolar capture and seqTether in low TE", MaterialKind.CONDITION, source=METHODS)
    duplex = g.product("capture_duplex", "freshly duplexed capture/seqTether", MaterialKind.ADAPTER,
        state=_duplex(), source=TABLE + "; partition derived by comparison of table rows")
    g.react("duplex_capture", "Anneal capture to seqTether (95 °C 2 min; cool at RT 1 h)", StepKind.HYBRIDIZATION,
        ((capture, InputRole.TEMPLATE), (tether, InputRole.STRUCTURAL), (low_te, InputRole.CONDITION)), (duplex,), source=METHODS)
    # Explicit local display dimensions are schematic, not measured telomere lengths.
    end = TelomereEnd(repeat="TTAGGG", duplex_copies=3, overhang_copies=2, subtelomere_bases=16)
    gdna = g.input("genomic_dna", "approximately 3 µg HMW genomic DNA, native chromosome termini",
        MaterialKind.BIOLOGICAL, state=end.scene(), primary=True, source=METHODS + "; representative TTAGGG register and drawn lengths are schematic")
    hybrid = g.product("capture_hybrid", "INFERRED — representative capture register at native telomere",
        MaterialKind.MOLECULAR, state=_telomere_capture(), source=METHODS + "; annealing register is schematic")
    g.react("anneal_capture", "INFERRED — capture telomere arm annealed in one representative register",
        StepKind.HYBRIDIZATION, ((gdna, InputRole.TEMPLATE), (duplex, InputRole.STRUCTURAL)), (hybrid,),
        source=METHODS, changes_molecule=False)
    ligase = g.input("capture_ligase", "4000 U T4 DNA ligase; rCutSmart + ATP, 100 µL", MaterialKind.CATALYST, source=METHODS)
    tagged = g.product("tagged", "capture-oligo-ligated chromosome ends; connectivity unresolved",
        MaterialKind.MOLECULAR, state=_product_state("post-ligation DNA"), inferred=True, source=METHODS)
    # The duplex has already entered hybrid. Do not consume the same aliquot again.
    # NICK_SEALING classifies source-stated ligation within this preassembled
    # complex only; it does not establish a particular sealed nick or substrate bond.
    g.react("ligate_capture", "INFERRED — unresolved product connectivity after capture ligation (37 °C overnight; 65 °C 10 min)",
        StepKind.NICK_SEALING, ((hybrid, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (tagged,),
        source=METHODS, geometry=AdapterLigationGeometry.unavailable(
            "The source specifies ligation within the preassembled capture/telomere complex. NICK_SEALING is an operational classification only: neither a particular nick nor its covalent product is assigned from the inspected Methods."))
    polymerase = g.input("pol_iv", "4 U Sulfolobus DNA polymerase IV; ThermoPol, dNTPs, ATP in same tube", MaterialKind.CATALYST, source=METHODS)
    filled = g.product("gap_filled", "INFERRED — potential-gap-filled capture/telomere assembly",
        MaterialKind.MOLECULAR, state=_product_state("post-Pol-IV DNA", gap_filled=True), inferred=True, source=METHODS)
    g.react("gap_fill", "INFERRED — fill potential gaps; tract unresolved (56 °C 2 min; 72 °C 15 min, 500 RPM)", StepKind.EXTENSION,
        ((tagged, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (filled,), source=METHODS)
    post_capture = _cleanup(g, filled, "capture_cleanup", "1.6×", "post-capture library, 60 µL ONT EB")
    motor = g.input("motor_adapter", "5 µL ONT AMII (R9) or NA (R10) sequencing adapter mix",
        MaterialKind.ADAPTER, inferred=True, source=METHODS + "; historical assembly is unassigned, no current-kit sequence substituted")
    motor_ligase = g.input("motor_ligase", "10 µL Quick T4 DNA ligase + 20 µL 5× Quick ligase buffer + 5 µL water",
        MaterialKind.CATALYST, source=METHODS)
    ligated = g.product("motor_ligated", "INFERRED — ONT adapter-ligated captured chromosome ends", MaterialKind.LIBRARY,
        state=_product_state("ONT adapter-ligated DNA", gap_filled=True, motor=True), inferred=True, source=METHODS)
    g.react("ligate_motor", "INFERRED — ligate ONT adapter; assembly unresolved (complete post-capture eluate, RT 1 h)", StepKind.LIGATION,
        ((post_capture, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL), (motor_ligase, InputRole.CATALYST)),
        (ligated,), source=METHODS, geometry=AdapterLigationGeometry.unavailable(
            "Historical AMII/NA motor adapter end geometry remains unassigned; the 2025 author revision is not substituted for the 2024 preparation."))
    library = _cleanup(g, ligated, "motor_cleanup", "1.1×", "DTM nanopore sequencing library")
    g.sequence(library, SequencingPlan(acquired_roles=("Sequencing output",), limitation=
        "R9.4.1 or R10.4 PromethION/P2Solo; 20–50 fmol, 24–72 h, per defining Methods. Platform motor entry; no sequencing primer or Illumina index cycles. This graph models one representative sample; the up-to-eight-sample 1000-ng pool is deferred."))
    # Explicit curator-written scope; newly generated nodes do not enrol themselves.
    g.declare_reviewed_scope(("capture", "tether"), TABLE)
    g.declare_reviewed_scope((
        "low_te", "genomic_dna", "capture_ligase", "pol_iv", "capture_cleanup_beads",
        "motor_adapter", "motor_ligase", "motor_cleanup_beads",
        "duplex_capture", "anneal_capture", "ligate_capture", "gap_fill",
        "capture_cleanup", "ligate_motor", "motor_cleanup", "motor_cleanup_eluate"), METHODS)
    return g.require_publishable()
