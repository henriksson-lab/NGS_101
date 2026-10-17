"""NEBNext Enzymatic Methyl-seq v2, E8015 manual v2.0 7/26."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment
from adapter_ligation import AdapterLigationGeometry
from rna_library_operations import molecule
from endprep import dA_tailed_scene, repair_and_dA_tail
from nebnext import STANDARD_PRIMERS, indexed_library
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
TITLE="NEBNext Enzymatic Methyl-seq v2"
NOTES="01_nebnext-em-seq-v2.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en/-/media/nebus/files/manuals/manuale8015.pdf">NEB E8015 manual v2.0 (July 2026)</a>.'
SUMMARY="The EM-seq Adaptor is ligated before enzymatic conversion. TET2 and T4-BGT protect 5mC/5hmC; after formamide (or NaOH) denaturation, APOBEC converts unmodified C to U; Q5U then amplifies the converted library with LV UDI primers."
OUTCOME={"C":"U","5mC":"C","5hmC":"C"}
DNA=Construct([Segment("genomic DNA", "X"*42, placeholder=True)], name="genomic DNA")
FRAGMENT=Construct([Segment("fragmented DNA (~350 bp)", "X"*42, placeholder=True)], name="fragmented DNA")
END_PREP=repair_and_dA_tail(FRAGMENT)
# Manual §1.3 / Fig. 1: end-repaired, dA-tailed DNA is ligated to the NEBNext EM-seq
# Adaptor. The adaptor's ligating end, strands and cytosine chemistry are not printed.
ADAPTOR_END=AdapterLigationGeometry.unavailable(
    "E8015 v2.0_7/26 names the NEBNext EM-seq Adaptor (E7165AVIAL) and ligates it to dA-tailed DNA; "
    "its strands, ligating end and how its cytosines survive APOBEC are not disclosed.")
def final_library(): return indexed_library(Construct([Segment("enzymatically converted insert", "X"*42, placeholder=True)]), "EM-seq v2 library")
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Paired reads sequence the converted insert; C calls represent protected 5mC/5hmC and T calls at genomic C positions represent unmodified cytosine. i5 and i7 are separate index reads."
FINAL_CAPTION="PCR-completed, dual-index EM-seq library. The manual says the libraries resemble TruSeq libraries and prints both adaptor trimming sequences (Overview p. 3), so primer binding is computed from the canonical arms."
TRIM_READ1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"  # E8015 v2.0_7/26 p. 3, Adaptor Read 1
TRIM_READ2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # E8015 v2.0_7/26 p. 3, Adaptor Read 2

def _beads(g, ident, label, template, beads, state):
    kept=g.product(ident,label,MaterialKind.MOLECULAR,state=state)
    waste=g.product(ident+"_supernatant",label+": bead supernatant and ethanol washes",MaterialKind.WASTE)
    g.react(ident+"_cleanup","Sample Purification Bead clean-up: "+label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(kept,waste))
    g.discard(waste,"unbound reaction components removed with the supernatant")
    return kept

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    dna=g.input("dna","0.1–200 ng double-stranded sample DNA (48 µl)",MaterialKind.BIOLOGICAL,state=DNA,primary=True)
    lam=g.input("lambda_control","diluted Control DNA Unmethylated Lambda (1 µl)",MaterialKind.BIOLOGICAL)
    puc=g.input("puc19_control","diluted Control DNA CpG methylated pUC19 (1 µl)",MaterialKind.BIOLOGICAL)
    shear=g.input("fragmentation","NEBNext UltraShear (M7634) or Covaris, to ~350 bp",MaterialKind.CONDITION)
    frag=g.product("fragmented","fragmented sample + control DNA (50 µl, no clean-up)",MaterialKind.MOLECULAR,state=FRAGMENT)
    g.react("fragment","Combine controls and fragment (§1.1)",StepKind.FRAGMENTATION,((dna,InputRole.TEMPLATE),(lam,InputRole.SPIKE_IN),(puc,InputRole.SPIKE_IN),(shear,InputRole.CONDITION)),(frag,))
    endmix=g.input("endprep_mix","Ultra II End Prep Buffer + Enzyme Mix; 20 °C 15 min, 65 °C 15 min",MaterialKind.CATALYST)
    repaired=g.product("repaired","end-repaired fragmented DNA",MaterialKind.MOLECULAR,state=END_PREP.insert)
    g.react("end_repair","Repair fragmented DNA ends (§1.2)",StepKind.END_REPAIR,((frag,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(repaired,))
    tailed=g.product("tailed","end-repaired, dA-tailed DNA",MaterialKind.MOLECULAR,state=dA_tailed_scene(END_PREP))
    g.react("da_tail","Add one 3-prime dA to each repaired end (same End Prep reaction)",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(tailed,))
    adapter=g.input("adapter","NEBNext EM-seq Adaptor (2.5 µl; composition not disclosed)",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (adaptor strands not printed)")
    ligase=g.input("ligase","Ultra II Ligation Master Mix + Ligation Enhancer; 20 °C 15 min",MaterialKind.CATALYST)
    adapted=g.product("adapted","EM-seq adaptor-ligated DNA",MaterialKind.MOLECULAR,
                      state=molecule("adaptor-ligated DNA",("EM-seq adaptor","dA-tailed insert","EM-seq adaptor"),joined=True),inferred=True,source=SOURCE+" (adaptor strands not printed)")
    g.react("ligate","Ligate the EM-seq Adaptor (§1.3)",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=ADAPTOR_END)
    beads=g.input("beads","NEBNext Sample Purification Beads (1X, 1X, 0.8X)",MaterialKind.SELECTION)
    ligated=_beads(g,"ligated_clean","1X (93 µl) clean-up of ligated DNA; ≤10 ng adds Carrier DNA",adapted,beads,adapted.state)
    protect=g.input("protection_enzymes","TET2 + T4-BGT (diluted 1:10 for ≤10 ng), UDP-glucose, DTT, TET2 buffer; Fe(II) starts oxidation, 37 °C 1 h",MaterialKind.CATALYST)
    protected=g.product("protected","5mC/5hmC → 5caC; 5hmC → 5gmC (protected)",MaterialKind.MOLECULAR,
                        state=molecule("protected adaptor-ligated DNA",("EM-seq adaptor","insert: 5mC/5hmC as 5caC/5gmC","EM-seq adaptor"),joined=True))
    g.react("protect","Protect 5mC and 5hmC (§1.5)",StepKind.CHEMICAL_CONVERSION,((ligated,InputRole.TEMPLATE),(protect,InputRole.CATALYST)),(protected,))
    stop=g.input("stop_reagent","Stop Reagent (1 µl), 37 °C 30 min",MaterialKind.CONDITION)
    stopped=g.product("stopped","oxidation-stopped protected DNA",MaterialKind.MOLECULAR,state=protected.state)
    g.react("stop","Stop the oxidation reaction",StepKind.QUENCHING,((protected,InputRole.TEMPLATE),(stop,InputRole.CONDITION)),(stopped,),changes_molecule=False)
    clean_protected=_beads(g,"protected_clean","1X (50 µl) clean-up of protected DNA",stopped,beads,protected.state)
    denat=g.input("denaturant","4 µl formamide (recommended) or 0.05 N NaOH, 85 °C 10 min, then cooling block on ice",MaterialKind.CONDITION)
    top=g.product("denatured_top","single-stranded protected top strand",MaterialKind.MOLECULAR,
                  state=molecule("protected top strand",("EM-seq adaptor","insert top: 5caC/5gmC, unmodified C","EM-seq adaptor")))
    bottom=g.product("denatured_bottom","single-stranded protected bottom strand",MaterialKind.MOLECULAR,
                     state=molecule("protected bottom strand",("EM-seq adaptor","insert bottom: 5caC/5gmC, unmodified C","EM-seq adaptor")))
    g.react("denature","Denature DNA (§1.7)",StepKind.STRAND_SEPARATION,((clean_protected,InputRole.TEMPLATE),(denat,InputRole.CONDITION)),(top,bottom))
    apobec=g.input("apobec","APOBEC + Deamination Buffer + Recombinant Albumin; 37/42/50/75 °C (≤100 ng) or 37 °C 3 h",MaterialKind.CATALYST)
    converted=g.product("converted","deaminated DNA: unmodified C → U; 5caC/5gmC unchanged",MaterialKind.MOLECULAR,
                        state=molecule("deaminated single strands (both original strands)",("EM-seq adaptor","insert: U at unmodified C","EM-seq adaptor")))
    g.react("deaminate","Deaminate unmodified cytosine (§1.8); no clean-up before PCR",StepKind.CHEMICAL_CONVERSION,((top,InputRole.TEMPLATE),(bottom,InputRole.TEMPLATE),(apobec,InputRole.CATALYST)),(converted,))
    fwd=g.input("pcr_forward","NEBNext LV UDI primer (i5 side), pair 5 µl",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (primer synthesis strands not printed)")
    rev=g.input("pcr_reverse","NEBNext LV UDI primer (i7 side)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (primer synthesis strands not printed)")
    q5u=g.input("q5u","NEBNext Q5U Master Mix; 98 °C, 4–14 cycles 98/62/65 °C",MaterialKind.CATALYST)
    amplified=g.product("amplified","PCR-amplified EM-seq library (U copied as T)",MaterialKind.LIBRARY,state=final_library())
    g.react("pcr","Amplify converted strands with Q5U (§1.9)",StepKind.PCR,((converted,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(q5u,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(final_library(),reason="The manuals do not disclose the complete adaptor or UDI primer synthesis sequences."))
    library=g.product("library","0.8X (72 µl) bead-purified EM-seq library",MaterialKind.LIBRARY,state=final_library())
    pcr_waste=g.product("pcr_supernatant","post-PCR bead supernatant and ethanol washes",MaterialKind.WASTE)
    g.react("pcr_cleanup","Sample Purification Bead clean-up 0.8X (§1.10)",StepKind.CLEANUP,((amplified,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(library,pcr_waste))
    g.discard(pcr_waste,"unbound PCR components and short fragments removed with the supernatant")
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),acquired_roles=REQUIRED_ROLES))
    g.annotate(SourceDiscrepancy(
        target="pcr", title="LV index kit identifiers in the E8015 manual",
        claims=(
            SourceClaim("E8015 v2.0_7/26, p. 2, Required Materials Not Included", "Names E3390 and E3392 among the compatible LV Unique Dual Index Primer Sets."),
            SourceClaim("E8015 v2.0_7/26, p. 2, Fig. 1 caption", "Prints E3900 and E3902 for the first two primer sets."),
            SourceClaim("NEB LV UDI manual, title and index tables", "Identifies the sold sets as E3390 and E3392.")),
        model_choice="The graph refers to the actual E3390/E3392 primer kits, with the chosen well left unspecified.",
        evidence="The E8015 required-materials list agrees with the dedicated LV UDI manual.",
        uncertainty="NEB has not confirmed a correction to the E8015 Fig. 1 caption. Earlier revisions (1.2 9/25, 1.3 4/26 in the revision history) have not been compared.",
        slot=AnnotationSlot.INTRO))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    return g.require_publishable()
