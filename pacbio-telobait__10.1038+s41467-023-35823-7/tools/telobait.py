"""Telobait enrichment for PacBio HiFi (Tham et al. 2023)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct,Row,Scene,Segment,complement_segments,feature,
                      feature_rows,revcomp,strand_row)
from dumbbell import (Dumbbell,dumbbell_rows,smrtbell_circuit,
                      smrtbell_hairpin_geometry,smrtbell_sequencing_plan)
from endprep import dA_tailed_scene,repair_and_dA_tail
from protocol_graph import (InputRole,MaterialKind,ProtocolGraph,SourceClaim,
                            SourceDiscrepancy,StepKind)
from telomere import AffinityCapture,TelomereEnd
from restriction import RSAI,HINFI,ECORI

END=TelomereEnd(); F01="ATCGACGGTTCAAGATGCCAGATGCACGGAGCAGAATTCTTTTAC"
R_PHASES=("CCCTAA","ACCCTA","AACCCT","TAACCC","CTAACC","CCTAAC")
SAMPLE_INDEX=feature("telobait_sample_index", "sample_index", "whitelist",
                     whitelist="Tham et al. Supplementary Data 1 telobait barcodes")
F01_PARTS=(Segment("telobait prefix",F01[:5],"r3"),
           Segment("8-nt telobait barcode",F01[5:13],"cbc",feature=SAMPLE_INDEX),
           Segment("telobait capture / EcoRI region",F01[13:],"r3"))
# Supplementary Data 2 (Set 1-12, F01/R01-1..6) and Supplementary Fig. S2a: the
# telobait is a duplex.  F (/5Phos/ ... /3Bio/) continues the native G-rich strand;
# each R partner is revcomp(F) plus a 3-prime six-base C-rich phase arm that anneals to
# the terminal bases of the G-rich 3-prime overhang.  T4 ligase joins the G-rich 3-prime
# end to the F 5-prime phosphate (Methods, "Telomere enrichment using telobaits").
R01=tuple(revcomp(F01)+arm for arm in R_PHASES)   # R01-1..R01-6 as printed

def telobait_scene(end:TelomereEnd=END)->Scene:
    """F01/R01-1 annealed to the native overhang that ends in 5'-TTAGGG-3'."""
    arm=R_PHASES[0]
    top=[Segment("subtelomere","X"*end.subtelomere_bases,placeholder=True),
         Segment("duplex telomere",end.duplex,"r1"),
         Segment("G-rich overhang",end.overhang,"r2")]
    sc=Scene(); sc.strand("G-rich strand",top,label="chromosome")
    sc.anneal("C-rich strand",complement_segments(top[:2]),to="G-rich strand",
              pair=("duplex telomere'","duplex telomere"),label="chromosome")
    r_strand=[*complement_segments(F01_PARTS),Segment("C-rich phase arm",arm,"r2")]
    sc.anneal("telobait R01-1",r_strand,to="G-rich strand",
              pair=("C-rich phase arm","G-rich overhang"),
              shift=len(end.overhang)-len(arm),label="telobait R")
    sc.anneal("telobait F01",list(F01_PARTS),to="telobait R01-1",
              pair=(F01_PARTS[0].name,F01_PARTS[0].name+"'"),above=True,
              label="telobait F",mod5="p",mod3="Bio")
    sc.ligation_end("telobait R01-1","C-rich phase arm",kind="cohesive",
                    polarity="3-prime",owner="adapter",
                    text="R 3′ phase arm pairs with the G-rich overhang")
    sc.footer("** T4 ligase joins the native G-rich 3′ end to the F01 5′ phosphate; "
              "the R 3′ end is later extended by T4 polymerase across the gap",
              "telobait F01",F01_PARTS[0].name)
    return sc

def telobait_reagent_scene()->Scene:
    """The free F01/R01-1 duplex before it meets a chromosome end."""
    r_strand=[*complement_segments(F01_PARTS),Segment("C-rich phase arm",R_PHASES[0],"r2")]
    sc=Scene(); sc.strand("telobait F01",list(F01_PARTS),label="telobait F",mod5="p",mod3="Bio")
    sc.anneal("telobait R01-1",r_strand,to="telobait F01",
              pair=(F01_PARTS[0].name+"'",F01_PARTS[0].name),label="telobait R")
    sc.ligation_end("telobait R01-1","C-rich phase arm",kind="cohesive",
                    polarity="3-prime",owner="adapter",
                    text="3′ phase arm (R01-1..R01-6 vary only here)")
    return sc

def ligated_scene(end:TelomereEnd=END,filled:bool=False)->Scene:
    """After ligation: G-rich strand continues into F01; a gap remains on the R side."""
    arm=R_PHASES[0]
    top=[Segment("subtelomere","X"*end.subtelomere_bases,placeholder=True),
         Segment("duplex telomere",end.duplex,"r1"),
         Segment("G-rich overhang",end.overhang,"r2"),*F01_PARTS]
    sc=Scene(); sc.strand("G-rich strand",top,label="chromosome + F01",mod3="Bio")
    sc.anneal("C-rich strand",complement_segments(top[:2]),to="G-rich strand",
              pair=("duplex telomere'","duplex telomere"),label="chromosome")
    gap=end.overhang[:len(end.overhang)-len(arm)]
    fill=[Segment("T4 polymerase fill-in",revcomp(gap),"r2")] if filled and gap else []
    r_strand=[*complement_segments(F01_PARTS),Segment("C-rich phase arm",arm,"r2"),*fill]
    sc.anneal("telobait R01-1",r_strand,to="G-rich strand",
              pair=("C-rich phase arm","G-rich overhang"),
              shift=len(end.overhang)-len(arm),label="telobait R")
    sc.junction("G-rich strand","G-rich overhang",F01_PARTS[0].name,text="ligated")
    if fill:
        sc.mark("telobait R01-1","T4 polymerase fill-in","fill-in; nick to the C-rich 5′ end")
    return sc

TELOBAIT_SCENE=telobait_scene()
TELOBAIT_LIGATED=ligated_scene()
TELOBAIT_FILLED=ligated_scene(filled=True)
TELOBAIT_REAGENT=telobait_reagent_scene()
CAPTURE=AffinityCapture("3′ biotin",ECORI.name)
# EcoRI (G^AATTC) cuts inside F01, so the released G-rich strand ends at the G and
# the R strand keeps a 5'-AATT overhang; biotin stays on the beads.
_CUT=F01.index(ECORI.site)+ECORI.top_cut
_TELOMERE=Segment("telomere","TTAGGG"*4,"r1")
_SUBTEL=Segment("subtelomere","X"*30,placeholder=True)
RELEASED_PARTS=(F01_PARTS[0],F01_PARTS[1],
                Segment("telobait arm to EcoRI cut",F01[13:_CUT],"r3"))
# Released duplex drawn by its G-rich strand; the 5'-AATT overhang on the R strand is
# stated in prose because a Construct carries one blunt duplex.
RELEASED=Construct([_SUBTEL,_TELOMERE,*RELEASED_PARTS],name="EcoRI-released telomere fragment")
# INFERRED: SMRTbell end repair is assumed to fill the 5'-AATT overhang (kit guide unread).
INSERT=Construct([_SUBTEL,_TELOMERE,*RELEASED_PARTS,
                  Segment("EcoRI overhang fill-in",ECORI.overhang,"r3",inferred=True)],
                 name="end-repaired telomere fragment")
SMRTBELL=Dumbbell(INSERT,"N"*16,"N"*16,name="telobait SMRTbell")
SOURCE='Tham et al., Nature Communications (2023), doi:10.1038/s41467-023-35823-7'

def panels(): return [("Restriction digest and terminal ligation in one incubation",[Row(chunks=[(f"subtelomere -- {RSAI.name}/{HINFI.name} cuts -- telomere -- native 3′ overhang",None,False)])],"RsaI, HinfI, phased telobaits and T4 ligase are co-incubated overnight; neither enzyme cuts canonical TTAGGG repeats."),("Anneal and ligate one of six phased, barcoded 3′-biotin telobaits",TELOBAIT_SCENE,"F01 pairs with six R partners whose 3′ arms cover all repeat phases; R01-1 is drawn. ** marks the G-rich 3′ end to F01 5′-phosphate ligation."),("Capture, fill the gap, and release with EcoRI",[Row(chunks=[("3′ biotin -> streptavidin capture -> T4 polymerase gap fill -> EcoRI release",None,False)])],"Gap filling happens on the captured molecules (Methods order); EcoRI (G^AATTC) leaves the biotinylated stub on the beads, and 0.45× AMPure PB purifies the released fragments."),("Build a PacBio HiFi SMRTbell",[*dumbbell_rows(SMRTBELL),strand_row(INSERT),*feature_rows(INSERT,prefix_width=4)],"INFERRED — SMRTbell end repair is drawn as filling the EcoRI 5′-AATT overhang. Hairpins close the selected fragment; the eight-base telobait sample index remains inside the insert. ** marks each sealed insert–adapter junction.")]


def protocol_graph():
    graph=ProtocolGraph("PacBio telobait HiFi sequencing",authority=SOURCE)
    dna=graph.input("genomic_dna","genomic DNA with native telomere overhangs",
                    MaterialKind.BIOLOGICAL,state=END.scene(),primary=True)
    restriction_enzymes=graph.input("restriction_enzymes","FastDigest RsaI and HinfI (2 µl each, 1× FastDigest buffer)",
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
    bait=graph.input("telobait","F01 annealed 1:6 with R01-1..R01-6 (95 °C 10 min, slow cool); 5 µl of 10 µM",
                     MaterialKind.ADAPTER,state=TELOBAIT_REAGENT)
    ligase=graph.input("terminal_ligase","T4 DNA ligase (2 µl, 30 U/µl) + 0.5 mM ATP, 37 °C overnight",MaterialKind.CATALYST)
    tagged=graph.product("tagged_telomeres","covalently telobait-tagged telomere DNA",
                         MaterialKind.MOLECULAR,state=TELOBAIT_LIGATED)
    terminal_geometry=AdapterLigationGeometry.known(
        "cohesive",TELOBAIT_REAGENT,
        statement=("Supplementary Data 2 and Fig. S2a: the R oligo's 3-prime C-rich phase "
                   "arm anneals to the G-rich overhang, and ligase joins the G-rich 3-prime "
                   "end to the 5-prime phosphate of F01."))
    graph.react("ligate_telobait","Anneal and ligate a phased telobait at the native terminus",
                StepKind.LIGATION,
                ((trimmed,InputRole.TEMPLATE),(bait,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(tagged,),geometry=terminal_geometry)
    streptavidin=graph.input("streptavidin","Dynabeads MyOne Streptavidin T1 (150 µl), RT 2 h, three 1× B&W washes",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","streptavidin-captured telobait molecules",
                           MaterialKind.MOLECULAR,state=tagged.state)
    uncaptured=graph.product("uncaptured","untagged genomic fragments",MaterialKind.WASTE)
    graph.react("capture","Capture 3-prime-biotin telobait molecules",StepKind.AFFINITY_CAPTURE,
                ((tagged,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,uncaptured))
    graph.discard(uncaptured,"does not carry the biotinylated telobait")
    polymerase=graph.input("fill_polymerase","T4 DNA polymerase (5 U) + 1 mM dNTP on beads, 37 °C 30 min",MaterialKind.CATALYST)
    filled=graph.product("gap_filled","gap-filled captured telobait molecules",
                         MaterialKind.MOLECULAR,state=TELOBAIT_FILLED)
    graph.react("gap_fill","Fill the captured telobait–chromosome gap",StepKind.EXTENSION,
                ((captured,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(filled,))
    ecori=graph.input("ecori","EcoRI",MaterialKind.CATALYST)
    released=graph.product("released","EcoRI-released enriched telomere fragments",
                           MaterialKind.MOLECULAR,state=RELEASED)
    bead_stub=graph.product("bead_stub","biotinylated EcoRI-distal telobait stub",
                            MaterialKind.WASTE)
    graph.react("release","Release captured telomere fragments with EcoRI",StepKind.DIGESTION,
                ((filled,InputRole.TEMPLATE),(ecori,InputRole.CATALYST)),(released,bead_stub))
    graph.discard(bead_stub,"remains bound to the streptavidin beads after EcoRI cleavage")
    ampure=graph.input("ampure_pb","0.45× AMPure PB beads",MaterialKind.SELECTION)
    purified=graph.product("purified","AMPure PB-purified released fragments",
                           MaterialKind.MOLECULAR,state=RELEASED)
    ampure_waste=graph.product("ampure_waste","AMPure PB supernatant and washes",
                               MaterialKind.WASTE)
    graph.react("ampure_cleanup","Purify released fragments with 0.45× AMPure PB",
                StepKind.CLEANUP,((released,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),
                (purified,ampure_waste))
    graph.discard(ampure_waste,"short material and buffer not retained by 0.45× beads")
    end_mix=graph.input("end_mix","SMRTbell Express Template Prep Kit 2.0 damage/end repair and dA-tailing",
                        MaterialKind.CATALYST)
    end_product=repair_and_dA_tail(INSERT,inferred=True)
    tailed=graph.product("tailed","repaired dA-tailed telomere inserts",
                         MaterialKind.MOLECULAR,state=dA_tailed_scene(end_product))
    graph.react("end_prep","Repair and dA-tail enriched telomere inserts",StepKind.END_REPAIR,
                ((purified,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(tailed,))
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
    graph.annotate(SourceDiscrepancy(
        "gap_fill","Published order of gap filling and streptavidin capture differs",
        (SourceClaim("Results, Telomere enrichment using telobaits (PMC9845338 JATS line 128)",
                     "after ligation T4 DNA polymerase fills the gap; fragments are subsequently purified on MyOne T1 beads"),
         SourceClaim("Methods, Telomere enrichment using telobaits (PMC9845338 JATS line 224)",
                     "beads are added to the ligation mix, washed, and T4 DNA polymerase gap filling is done on the beads")),
        "The model follows the Methods order: capture, then on-bead gap fill, then EcoRI release.",
        "Methods give volumes, buffers and washes for each on-bead step; Results is a summary.",
        "Whether the Results sentence describes order or only the set of operations is not stated."))
    graph.declare_reviewed_scope((
        "genomic_dna","restriction_enzymes","trim_genome","telobait","terminal_ligase",
        "ligate_telobait","fill_polymerase","gap_fill","streptavidin","capture",
        "ecori","release","ampure_pb","ampure_cleanup","end_mix","end_prep","hairpins","smrt_ligase",
        "ligate_hairpins","library"),SOURCE)
    return graph.require_publishable()
