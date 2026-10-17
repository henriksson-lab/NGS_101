"""Oxford Nanopore SQK-RNA004 direct-RNA library architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import AdapterLigationGeometry
from nanopore_library import motor_library, motor_sequencing_plan, role
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role as molecular_role

TITLE = "Oxford Nanopore Direct RNA Sequencing — SQK-RNA004"
NOTES = "01_direct-rna.html"
SOURCE = 'Commercial protocol: Oxford Nanopore <a href="https://nanoporetech.com/document/direct-rna-sequencing-sqk-rna004">Direct RNA sequencing (SQK-RNA004)</a>, DRS_9195_v4_revL_29June2026; architecture paper: <a href="https://doi.org/10.1038/nmeth.4577">Garalde et al. (2018)</a>.'
SUMMARY = "The RT Adapter (RTA) is ligated to the RNA 3′ end with T4 DNA ligase, together with the optional RNA CS control strand. Induro reverse transcriptase then makes a stabilizing cDNA, and the hybrid is cleaned up on RNAClean XP beads. A second T4 DNA ligase step attaches the RNA Ligation Adapter (RLA) sequencing adapter, and the library is washed with WSB and eluted in REB. The nanopore reads the native RNA strand from its 3′ end toward 5′."
CAVEAT = ("The inspected protocol sections give no bases for the RTA or RLA, no motor identity and no ligating-end geometry. Bracketed components below are structural roles, not invented bases. "
          "The protocol accepts 300 ng poly(A)-tailed RNA or 1 µg total RNA. That RTA targets poly(A) is inferred from this input requirement and from the technical document; the protocol does not print an oligo(dT) sequence.")
FINAL_LIBRARY = motor_library([
    role("RLA/RTA boundary", "[RLA to RTA, junction undisclosed]", "me"),
    role("RT Adapter (RTA)", "[RTA]", "r2"),
    role("native RNA:cDNA hybrid", "[native RNA with stabilizing cDNA]"),
], adapter="sequencing adapter", name="SQK-RNA004 pore-entry complex")
SEQ_PRIMERS = ()
SEQUENCING_MECHANISM = "nanopore_motor"
MOTOR_ADAPTER_NAME = "RNA sequencing adapter"
SEQUENCING_ENDING = "No sequencing primer extends this library. The motor-loaded sequencing adapter captures the RNA–cDNA hybrid at the pore and controls passage of the native RNA strand, beginning at its poly(A)/3′ end and proceeding toward the RNA 5′ end. The cDNA stabilizes the RNA but is not the strand base-called in direct-RNA mode."
SEQUENCING_UNAVAILABLE = ("Direct RNA sequencing is motor-driven; current RTA and "
                          "sequencing-adapter sequences are proprietary.")
FINAL_CAPTION = "The motor presents native RNA from its poly(A)/3-prime end; cDNA is stabilizing material."
RTA_END = AdapterLigationGeometry.unavailable(
    "The vendor specifies RTA ligation to the RNA 3-prime end but withholds the current ligating-end sequence.")
RLA_END = AdapterLigationGeometry.unavailable(
    "The vendor ligates RLA to the RNA-cDNA hybrid with T4 DNA ligase but does not disclose the adapter ends or junction.")

def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("rna", "native poly(A) RNA", MaterialKind.MOLECULAR,
                  state=molecule("native RNA", (molecular_role("RNA body"),
                                                 molecular_role("poly(A) tail"))),
                  primary=True)
    rta = g.input("rta", "reverse-transcription adapter", MaterialKind.ADAPTER,
                  inferred=True, source=SOURCE + " (current RTA bases proprietary)")
    rcs = g.input("rcs", "RNA CS (RCS) control strand, optional (yeast enolase II)",
                  MaterialKind.MOLECULAR)
    ligase = g.input("ligase", "T4 DNA ligase (2M U/ml) in NEBNext Quick Ligation buffer, with murine RNase inhibitor",
                     MaterialKind.CATALYST)
    adapted = g.product("adapted", "RTA-ligated native RNA", MaterialKind.MOLECULAR,
                        inferred=True, source=SOURCE + " (RTA bases proprietary)")
    g.react("ligate_rta", "Ligate RTA (and optional RCS) to RNA, 10 min room temperature",
            StepKind.LIGATION,
            ((rna, InputRole.TEMPLATE), (rta, InputRole.STRUCTURAL),
             (rcs, InputRole.SPIKE_IN), (ligase, InputRole.CATALYST)), (adapted,),
            geometry=RTA_END)
    rt = g.input("rt", "Induro reverse transcriptase with dNTPs", MaterialKind.CATALYST)
    hybrid = g.product("hybrid", "native RNA with stabilizing cDNA",
                       MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Synthesize the stabilizing cDNA, 60 °C 30 min then 70 °C 10 min",
            StepKind.REVERSE_TRANSCRIPTION,
            ((adapted, InputRole.TEMPLATE), (rta, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (hybrid,))
    rt_beads = g.input("rt_beads", "RNAClean XP beads (72 µl), 70% ethanol washes",
                       MaterialKind.SELECTION)
    rt_clean = g.product("rt_clean", "bead-purified RNA–cDNA hybrid in 23 µl water",
                         MaterialKind.MOLECULAR)
    rt_waste = g.product("rt_waste", "RT supernatant and ethanol washes",
                         MaterialKind.WASTE)
    g.react("rt_cleanup", "Purify the RT reaction on RNAClean XP beads", StepKind.CLEANUP,
            ((hybrid, InputRole.TEMPLATE), (rt_beads, InputRole.SELECTION)),
            (rt_clean, rt_waste), changes_molecule=False)
    g.discard(rt_waste, "unbound RT reaction components removed with the bead supernatant")
    adapter = g.input("sequencing_adapter", "RNA Ligation Adapter (RLA) sequencing adapter",
                      MaterialKind.ADAPTER, inferred=True,
                      source=SOURCE + " (adapter bases, motor and ligating ends undisclosed)")
    ligase2 = g.input("rla_ligase", "T4 DNA ligase (2M U/ml) in NEBNext Quick Ligation buffer",
                      MaterialKind.CATALYST)
    ligated = g.product("ligated", "RLA-ligated RNA–cDNA hybrid", MaterialKind.MOLECULAR,
                        inferred=True, source=SOURCE + " (adapter bases undisclosed)")
    g.react("ligate_rla", "Ligate RLA to the RNA–cDNA hybrid, 10 min room temperature",
            StepKind.LIGATION,
            ((rt_clean, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (ligase2, InputRole.CATALYST)), (ligated,), geometry=RLA_END)
    beads = g.input("rla_beads", "RNAClean XP beads (16 µl), two 150 µl WSB washes, 13 µl REB elution",
                    MaterialKind.SELECTION)
    library = g.product("library", "motor-adapted native RNA library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (adapter bases proprietary)")
    waste = g.product("rla_waste", "supernatant and WSB washes with free adapter",
                      MaterialKind.WASTE)
    g.react("rla_cleanup", "Remove free adapter on RNAClean XP beads with WSB washes",
            StepKind.CLEANUP,
            ((ligated, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
            (library, waste), changes_molecule=False)
    g.discard(waste, "free RLA and ligation components removed in the WSB washes")
    g.sequence(library, motor_sequencing_plan(
        FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE,
        read_uses={"Nanopore strand": "include native RNA; stabilizing cDNA is not base-called"}))
    g.declare_reviewed_scope(("rna", "rta", "rcs", "ligase", "ligate_rta", "rt",
                              "reverse_transcribe", "rt_beads", "rt_cleanup",
                              "sequencing_adapter", "rla_ligase", "ligate_rla",
                              "rla_beads", "rla_cleanup", "library"), SOURCE)
    return g.require_publishable()
