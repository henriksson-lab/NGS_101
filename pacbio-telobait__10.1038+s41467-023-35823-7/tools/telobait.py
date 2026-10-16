"""Telobait enrichment for PacBio HiFi (Tham et al. 2023)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,Row,Segment,feature,feature_rows,strand_row
from dumbbell import (Dumbbell,dumbbell_rows,smrtbell_circuit,
                      smrtbell_hairpin_geometry,smrtbell_sequencing_plan)
from endprep import dA_tailed_scene,repair_and_dA_tail
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,StepKind
from telomere import AffinityCapture,PhasedTelomereAdapter,TelomereEnd
from restriction import RSAI,HINFI,ECORI

END=TelomereEnd(); F01="ATCGACGGTTCAAGATGCCAGATGCACGGAGCAGAATTCTTTTAC"
R_PHASES=("CCCTAA","ACCCTA","AACCCT","TAACCC","CTAACC","CCTAAC")
SAMPLE_INDEX=feature("telobait_sample_index", "sample_index", "whitelist",
                     whitelist="Tham et al. Supplementary Data 1 telobait barcodes")
F01_PARTS=(Segment("telobait prefix",F01[:5],"r3"),
           Segment("8-nt telobait barcode",F01[5:13],"cbc",feature=SAMPLE_INDEX),
           Segment("telobait capture / EcoRI region",F01[13:],"r3"))
TELOBAIT=PhasedTelomereAdapter(F01,arm_copies=1,core_segments=F01_PARTS)
CAPTURE=AffinityCapture("3′ biotin",ECORI.name)
INSERT=Construct([Segment("subtelomere","X"*30,placeholder=True),
                  Segment("telomere","TTAGGG"*4,"r1"),
                  *F01_PARTS],name="captured telomere fragment")
SMRTBELL=Dumbbell(INSERT,"N"*16,"N"*16,name="telobait SMRTbell")
SOURCE='Tham et al., Nature Communications (2023), doi:10.1038/s41467-023-35823-7'

def panels(): return [("Restriction digest and terminal ligation in one incubation",[Row(chunks=[(f"subtelomere -- {RSAI.name}/{HINFI.name} cuts -- telomere -- native 3′ overhang",None,False)])],"RsaI, HinfI, phased telobaits and T4 ligase are co-incubated overnight; neither enzyme cuts canonical TTAGGG repeats."),("Anneal and ligate one of six phased, barcoded 3′-biotin telobaits",TELOBAIT.annealed_scene(END),"The six ends cover all repeat phases. ** marks the terminal ligation; displayed F01 is one of the published barcoded sets."),("Capture, fill the gap, and release with EcoRI",[Row(chunks=[("3′ biotin -> streptavidin capture -> T4 polymerase gap fill -> EcoRI release",None,False)])],"Gap filling happens on the captured molecules before EcoRI releases the enriched fragments."),("Build a PacBio HiFi SMRTbell",[*dumbbell_rows(SMRTBELL),strand_row(INSERT),*feature_rows(INSERT,prefix_width=4)],"Hairpins close the selected fragment; the eight-base telobait sample index remains inside the insert. ** marks each sealed insert–adapter junction.")]


def protocol_graph():
    graph=ProtocolGraph("PacBio telobait HiFi sequencing",authority=SOURCE)
    dna=graph.input("genomic_dna","genomic DNA with native telomere overhangs",
                    MaterialKind.BIOLOGICAL,state=END.scene(),primary=True)
    restriction_enzymes=graph.input("restriction_enzymes","RsaI and HinfI",
                                    MaterialKind.CATALYST)
    trimmed=graph.product("trimmed_telomeres","restriction-trimmed telomere-bearing DNA",
                          MaterialKind.MOLECULAR,state=END.scene())
    restriction_waste=graph.product("restriction_waste","distal restriction fragments",
                                    MaterialKind.WASTE)
    graph.react("trim_genome","Digest with RsaI/HinfI during terminal telobait ligation",
                StepKind.DIGESTION,
                ((dna,InputRole.TEMPLATE),(restriction_enzymes,InputRole.CATALYST)),
                (trimmed,restriction_waste))
    graph.discard(restriction_waste,"distal genomic fragments do not carry the native telomere end")
    bait=graph.input("telobait","six-phase barcoded 3-prime-biotin telobait pool",
                     MaterialKind.ADAPTER,state=TELOBAIT.annealed_scene(END))
    ligase=graph.input("terminal_ligase","telobait DNA ligase",MaterialKind.CATALYST)
    tagged=graph.product("tagged_telomeres","covalently telobait-tagged telomere DNA",
                         MaterialKind.MOLECULAR,state=TELOBAIT.annealed_scene(END))
    terminal_geometry=AdapterLigationGeometry.known(
        "cohesive",TELOBAIT.annealed_scene(END),
        statement=("The source prints six terminal repeat phases; the phased bait arm "
                   "anneals to the native G-rich overhang before ligase seals the end."))
    graph.react("ligate_telobait","Anneal and ligate a phased telobait at the native terminus",
                StepKind.LIGATION,
                ((trimmed,InputRole.TEMPLATE),(bait,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(tagged,),geometry=terminal_geometry)
    streptavidin=graph.input("streptavidin","streptavidin capture matrix",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","streptavidin-captured telobait molecules",
                           MaterialKind.MOLECULAR,state=tagged.state)
    uncaptured=graph.product("uncaptured","untagged genomic fragments",MaterialKind.WASTE)
    graph.react("capture","Capture 3-prime-biotin telobait molecules",StepKind.AFFINITY_CAPTURE,
                ((tagged,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,uncaptured))
    graph.discard(uncaptured,"does not carry the biotinylated telobait")
    polymerase=graph.input("fill_polymerase","T4 DNA polymerase",MaterialKind.CATALYST)
    filled=graph.product("gap_filled","gap-filled captured telobait molecules",
                         MaterialKind.MOLECULAR,state=captured.state)
    graph.react("gap_fill","Fill the captured telobait–chromosome gap",StepKind.EXTENSION,
                ((captured,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(filled,))
    ecori=graph.input("ecori","EcoRI",MaterialKind.CATALYST)
    released=graph.product("released","EcoRI-released enriched telomere fragments",
                           MaterialKind.MOLECULAR,state=INSERT)
    graph.react("release","Release captured telomere fragments with EcoRI",StepKind.DIGESTION,
                ((filled,InputRole.TEMPLATE),(ecori,InputRole.CATALYST)),(released,))
    end_mix=graph.input("end_mix","PacBio damage-repair, end-repair and dA-tailing mix",
                        MaterialKind.CATALYST)
    end_product=repair_and_dA_tail(INSERT,inferred=True)
    tailed=graph.product("tailed","repaired dA-tailed telomere inserts",
                         MaterialKind.MOLECULAR,state=dA_tailed_scene(end_product))
    graph.react("end_prep","Repair and dA-tail enriched telomere inserts",StepKind.END_REPAIR,
                ((released,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(tailed,))
    hairpins=graph.input("hairpins","T-overhang PacBio SMRTbell hairpin adapters",
                         MaterialKind.ADAPTER,inferred=True,
                         source=SOURCE+" (hairpin-body bases proprietary)")
    smrt_ligase=graph.input("smrt_ligase","SMRTbell ligase",MaterialKind.CATALYST)
    circuit=smrtbell_circuit(SMRTBELL)
    library=graph.product("library",SMRTBELL.name,MaterialKind.LIBRARY,state=circuit,
                          inferred=True,source=SOURCE+" (hairpin-body bases proprietary)")
    graph.react("ligate_hairpins","Ligate hairpins to both telomere-insert ends",
                StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(hairpins,InputRole.STRUCTURAL),
                 (smrt_ligase,InputRole.CATALYST)),(library,),
                geometry=smrtbell_hairpin_geometry(
                    SMRTBELL,statement="PacBio SMRTbell preparation ligates T-overhang hairpins to dA-tailed inserts."))
    graph.sequence(library,smrtbell_sequencing_plan(
        SMRTBELL,circuit,
        limitation="The source names PacBio primer/polymerase entry but does not print the proprietary hairpin primer-site bases."))
    graph.declare_reviewed_scope((
        "genomic_dna","restriction_enzymes","trim_genome","telobait","terminal_ligase",
        "ligate_telobait","fill_polymerase","gap_fill","streptavidin","capture",
        "ecori","release","end_mix","end_prep","hairpins","smrt_ligase",
        "ligate_hairpins","library"),SOURCE)
    return graph.require_publishable()
