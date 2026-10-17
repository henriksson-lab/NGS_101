"""Telo-seq targeted Oxford Nanopore workflow (Schmidt et al. 2024).

Source choices (round 2026-10-17-breadth-first-20): Methods "Telo-seq library
preparation" (PMC11189484 JATS ``.xml.txt`` line 169) and Supplementary Table 7
(Supplementary Information ``MOESM1_ESM.pdf.txt`` lines 1780-1796).  The graph is
assembled here from those passages; only the neutral phase algebra
(``PhasedTelomereAdapter``) is imported from ``lib/telomere``.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Row, Segment, oligo, revcomp
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
from telomere import PhasedTelomereAdapter
from restriction import ECORV

TITLE="Telo-seq"
SOURCE='Schmidt et al. 2024, <a href="https://doi.org/10.1038/s41467-024-48917-7">doi:10.1038/s41467-024-48917-7</a>.'
METHODS=SOURCE+" Methods, Telo-seq library preparation (PMC11189484 JATS line 169)."
TABLE7=SOURCE+" Supplementary Table 7 (Supplementary Information PDF lines 1780–1796)."

CORE="AGCAATACGTAACTGAACGAAGT"
# Supplementary Table 7 prints these six /5Phos/ telorettes as T1, T2, T3, T3, T5, T6.
TELORETTES=(
"AGCAATACGTAACTGAACGAAGTCCCTAACCCTAACCCTAA",
"AGCAATACGTAACTGAACGAAGTTAACCCTAACCCTAACCC",
"AGCAATACGTAACTGAACGAAGTCTAACCCTAACCCTAACC",
"AGCAATACGTAACTGAACGAAGTCCTAACCCTAACCCTAAC",
"AGCAATACGTAACTGAACGAAGTAACCCTAACCCTAACCCT",
"AGCAATACGTAACTGAACGAAGTACCCTAACCCTAACCCTA")
PRINTED_LABELS=("T1","T2","T3","T3","T5","T6")
S1="ACTTCGTTCAGTTACGTATTGCTAGCAAT"
ADAPTER=PhasedTelomereAdapter(CORE)
# Derived: the first 23 nt of S1 are the reverse complement of the telorette core;
# the remaining 3' extension is what S1 offers the AMII adapter (not verifiable,
# because the AMII sequence is not disclosed in the source).
S1_CORE_PAIRING=S1[:len(CORE)]==revcomp(CORE)
S1_EXTENSION=S1[len(CORE):]

DISCREPANCY=("Source discrepancy: Schmidt et al. Methods (JATS line 169) specify six phased "
    "telorette adapters, but Supplementary Table 7 (Supplementary Information PDF lines "
    "1780–1796) labels both the third and fourth distinct sequences T3 and has no T4 row. "
    "The model retains all six printed sequences in table order, treating the fourth "
    "sequence as the fourth phase; the intended fourth label has not been confirmed by the authors.")


def oligos():
    rows=[oligo(f"{label} (printed; row {i})",(Segment("telorette core",CORE),Segment("telomere-complement arm",seq[len(CORE):])),mods="/5Phos/")
          for i,(label,seq) in enumerate(zip(PRINTED_LABELS,TELORETTES),1)]
    rows.append(oligo("S1 splint",(Segment("telorette-core complement",S1[:len(CORE)]),Segment("3′ extension",S1_EXTENSION))))
    return tuple(rows)


def _ont_geometry():
    return AdapterLigationGeometry.unavailable("The Telo-seq source identifies the ONT AMII motor adapter and the disclosed S1 recruitment splint, but ONT does not publish the motor adapter's terminal ligating bases in the cited protocol.")


def _cleanup(g, substrate, ident, label, beads_label, kept_label):
    beads=g.input(ident+"_beads", beads_label, MaterialKind.SELECTION, source=METHODS)
    kept=g.product(ident+"_eluate", kept_label, MaterialKind.MOLECULAR, state=substrate.state, source=METHODS)
    waste=g.product(ident+"_waste", "AMPure XP supernatant and washes", MaterialKind.WASTE, source=METHODS)
    g.react(ident, label, StepKind.CLEANUP, ((substrate,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
            (kept,waste), source=METHODS, changes_molecule=False)
    g.discard(waste, "removed with the bead supernatant and washes")
    return kept


def protocol_graph() -> ProtocolGraph:
    g=ProtocolGraph(TITLE, authority=METHODS)
    gdna=g.input("genomic_dna","15 µg extracted genomic DNA with native telomeres",MaterialKind.BIOLOGICAL,
                 state=molecule("chromosome ends",(role("subtelomere"),role("native telomere terminus"))),primary=True)
    telorettes=g.input("telorettes","six /5Phos/ phased telorettes (Supplementary Table 7), pre-mixed equally to 1 µM; 20 µL",
                       MaterialKind.OLIGO, state=TELORETTES, source=TABLE7)
    ligase=g.input("t4_ligase","10,000 U T4 DNA ligase (NEB M0202) in 1x rCutSmart + 0.4 mM ATP, 200 µL",MaterialKind.CATALYST)
    tagged=g.product("telorette_ligated","telorette ligated to the native C-rich strand at the G-overhang",MaterialKind.MOLECULAR,state=gdna.state)
    g.react("ligate_telorettes","Anneal telorettes to the G-overhang and ligate to the C-strand (35 °C 16 h; 65 °C 10 min)",StepKind.LIGATION,
            ((gdna,InputRole.TEMPLATE),(telorettes,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(tagged,),
            geometry=AdapterLigationGeometry.unavailable("Telorette-to-C-strand junction: the native overhang phase is sample-dependent, so the source supplies six phased arms rather than one end geometry."))
    ecorv=g.input("ecorv","80 U EcoRV-HF (NEB R3195)",MaterialKind.CATALYST)
    cut=g.product("ecorv_cut","telorette-tagged DNA with distal blunt EcoRV ends",MaterialKind.MOLECULAR,state=tagged.state)
    g.react("digest_ecorv","Cut genomic DNA with EcoRV-HF (37 °C 30 min; 65 °C 20 min)",StepKind.DIGESTION,
            ((tagged,InputRole.TEMPLATE),(ecorv,InputRole.CATALYST)),(cut,))
    klenow=g.input("klenow","45 U Klenow Fragment (3′→5′ exo−) (NEB E6053) in 1x NEBNext dA-tailing buffer, 250 µL as 5 × 50 µL",MaterialKind.CATALYST)
    tailed=g.product("distal_tailed","dA-tailed blunt EcoRV ends (reduces concatemer ligation)",MaterialKind.MOLECULAR,state=cut.state)
    g.react("tail_distal","dA-tail the blunt EcoRV ends (37 °C 30 min)",StepKind.TAILING,
            ((cut,InputRole.TEMPLATE),(klenow,InputRole.CATALYST)),(tailed,))
    clean1=_cleanup(g,tailed,"ampure_1x","Purify with 1x AMPure XP; elute in 190 µL water (37 °C 15 min)",
                    "1x v/v AMPure XP (Beckman Coulter A63881)","purified dA-tailed telorette-tagged DNA")
    s1=g.input("s1_splint","S1 splint oligo (Supplementary Table 7), 100 nM, with 50 mM NaCl",MaterialKind.OLIGO,state=S1,source=TABLE7)
    annealed=g.product("s1_annealed","S1 annealed to the telorette core",MaterialKind.MOLECULAR,state=clean1.state)
    g.react("anneal_s1","Anneal the S1 splint (50 °C 1 h)",StepKind.HYBRIDIZATION,
            ((clean1,InputRole.TEMPLATE),(s1,InputRole.STRUCTURAL)),(annealed,))
    clean2=_cleanup(g,annealed,"ampure_05x","Purify with 0.5x AMPure XP; elute in 30 µL water (37 °C 15 min)",
                    "0.5x v/v AMPure XP","S1-annealed telorette-tagged DNA")
    motor=g.input("motor_adapter","AMII motor-loaded sequencing adapter (SQK-NBD111 kit instructions)",MaterialKind.ADAPTER,inferred=True,
                  source=METHODS+" The kit- and revision-specific AMII assembly is not disclosed.")
    kit_ligase=g.input("motor_ligase","SQK-NBD111 sequencing-adapter ligation reagents",MaterialKind.CATALYST,inferred=True,
                       source=METHODS+" Reagent identity follows the kit instructions, which were not inspected.")
    library=g.product("library","Telo-seq nanopore library",MaterialKind.LIBRARY,state=Construct([
        Segment("AMII motor adapter","[motor adapter]",placeholder=True,inferred=True),
        Segment("S1:telorette core","[S1:telorette]",placeholder=True),
        Segment("telomere to subtelomere","[telomere -- subtelomere -- EcoRV|A]",placeholder=True)],
        name="Telo-seq nanopore library"))
    g.react("ligate_motor","Ligate the AMII sequencing adapter via the S1 splint (SQK-NBD111 instructions)",StepKind.LIGATION,
            ((clean2,InputRole.TEMPLATE),(motor,InputRole.STRUCTURAL),(kit_ligase,InputRole.CATALYST)),(library,),
            geometry=_ont_geometry())
    g.sequence(library,SequencingPlan(acquired_roles=("Sequencing output",),
        limitation="GridION, R9.4.1 flow cells, one flow cell per library; bases called with a customized Bonito telomere model. Motor entry is platform-driven; no sequencing primer binds the library."))
    g.declare_reviewed_scope(tuple(list(g.materials)+list(g.steps)),METHODS)
    return g.require_publishable()


def panels(): return [("Anneal six phased telorettes and ligate the native telomere terminus",ADAPTER.annealed_scene(),"All oligos are 5′ phosphorylated. ** marks ligation to the native C-rich strand."),("EcoRV cuts distal genomic DNA; Klenow exo− dA-tails only those new blunt ends",[Row(chunks=[(f"telorette ** telomere -- subtelomere -- {ECORV.name}|A     (internal cut end only)",None,False)])],"The terminal telorette end is preserved while internal EcoRV ends receive dA, suppressing wrong-end adapter ligation."),("ONT adapter terminal geometry",_ont_geometry().rows(),_ont_geometry().statement),("S1 splint recruits the Oxford Nanopore sequencing adapter",[Row(chunks=[("[motor-loaded ONT adapter] ** [S1:telorette core] ** telomere -- subtelomere -- EcoRV|A",None,False)])],"S1 bridges the disclosed telorette core to the AMII adapter. ** marks the ligated boundaries.")]
