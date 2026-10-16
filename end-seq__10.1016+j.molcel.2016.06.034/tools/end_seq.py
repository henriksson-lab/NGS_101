"""Canela et al. END-seq model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import joined_scene, seg, truseq_library
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState, Row, Scene, Segment, Workflow
from end_capture import UserHairpin, captured_break_scene, user_opening_rows
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from topology import HairpinTopology
import seqprimers as sp

TITLE = "END-seq — direct capture of DNA double-strand-break ends"
NOTES = "01_end-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1016/j.molcel.2016.06.034">Canela et al., <i>Molecular Cell</i> (2016)</a>.'
SUMMARY = "Native DNA ends are blunted and trapped in agarose with a biotinylated P5 hairpin; after shearing and enrichment, a P7 hairpin captures the opposite library end and USER opens both loops for indexed PCR."
PCR_CYCLES = 16
READ_SPECS = (sp.read_spec("Read 1", "50 (HiSeq2500) or75 (NextSeq500)", "processed break-derived genomic end; source alternative runs", "include"),)
CAVEAT = "INFERRED reference endpoint — actual index members/widths are unassigned; the source specifies single-end sequencing, not acquired dual-index or paired-end reads. Blunting deliberately erases the original overhang chemistry. Read 1 starts at the first base retained at the processed break end."
A1 = UserHairpin("ENDseq adaptor 1", "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGUUTUTUUACACTCTTTCCCTACACGACGCTCTTCCGATCT", (33,35))
A2 = UserHairpin("ENDseq adaptor 2", "GATCGGAAGAGCACACGTCUUUUUUUUAGACGTGTGCTCTTCCGATCT")
A1_STEM_NT, A1_LOOP_NT = 31, 8
A2_STEM_NT, A2_LOOP_NT = 19, 9
A1_TOPOLOGY = HairpinTopology(A1.name, A1.sequence, 0,
                              A1_STEM_NT + A1_LOOP_NT, A1_STEM_NT,
                              colour="var(--c-r1)")
A2_TOPOLOGY = HairpinTopology(A2.name, A2.sequence, 0,
                              A2_STEM_NT + A2_LOOP_NT, A2_STEM_NT,
                              colour="var(--c-r2)")
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([seg("break-derived genomic insert", "X"*42, placeholder=True)], "END-seq library", inferred_adapters=True)
SEQ_PRIMERS = tuple(p for p in SEQ_PRIMERS if p.role == "Read 1")
FINAL_CAPTION = "INFERRED reference outer library — USER-opened hairpin arms are completed by TruSeq index PCR; the captured DSB boundary is adjacent to the Read 1 arm."
SEQUENCING_INTRO = "Single-end END-seq uses Read 1 to enter the captured genomic end; the complete construct also shows the standard TruSeq index and opposite-end primer sites created by PCR."

def oligos():
    """Source-written phosphate, biotin-dT, U and terminal phosphorothioate retained."""
    rows=[]
    for adapter in (A1,A2):
        written="".join("[Biotin-dT]" if i in adapter.biotin_dt else b
                        for i,b in enumerate(adapter.sequence))
        written=written[:-1]+"*"+written[-1]
        rows.append(adapter.name+" — 5′-Phos-"+written+"-3′; U denotes USER-sensitive deoxyuridine")
    return tuple(rows)

def hairpin_ligation_geometry(adapter: UserHairpin):
    sc=Scene();sc.strand("prepared DNA",[Segment("DNA duplex","X"*24,placeholder=True),Segment("DNA 3′ dA","A")],label="dA-tailed DNA")
    sc.overhang("prepared DNA","DNA 3′ dA",polarity="3-prime",owner="substrate",text="substrate 3′-A")
    sc.strand("hairpin adapter",[Segment("USER hairpin body",adapter.sequence[:-1],"r1"),Segment("hairpin 3′ dT",adapter.sequence[-1],"r1")],label=adapter.name,mod5="p")
    sc.overhang("hairpin adapter","hairpin 3′ dT",polarity="3-prime",owner="adapter",text="adapter 3′-T")
    return AdapterLigationGeometry.known("cohesive",sc,statement=f"The END-seq protocol prints {adapter.name}; its terminal 3′-T pairs with the prepared DNA 3′-A overhang.")

def sections():
    captured = [seg("P5 hairpin", "X"*20, "r1", placeholder=True),
                seg("break-derived fragment", "X"*34, placeholder=True),
                seg("P7 hairpin", "X"*20, "r2", placeholder=True)]
    return [
        ("Deproteinize DNA in agarose, blunt breaks and dA-tail", [Row(chunks=[("genomic DNA ---- processed DSB end-A", None, False)])],
         "Embedding before manipulation limits adventitious shearing; polymerases make the captured end blunt before one dA is added."),
        ("Ligate biotinylated P5 hairpin adaptor 1 in the plug", captured_break_scene().rows(),
         "The 3-prime T overhang accepts the break-end dA; two biotin-dU residues later support streptavidin capture and USER opening."),
        ("Shear, enrich captured molecules and ligate P7 hairpin adaptor 2", joined_scene(captured, (("P5 hairpin", "break-derived fragment", "adaptor 1 ligation"), ("break-derived fragment", "P7 hairpin", "adaptor 2 ligation")), label="two-hairpin captured fragment").rows(),
         "Sonication creates the distal end; its repair and dA-tail precede the second hairpin ligation on streptavidin beads."),
        ("Open both hairpins with USER and index-amplify", user_opening_rows(A1,A2),
         "USER cleavage at deoxyuridines exposes the two PCR-addressable Illumina arms."),
    ]

def workflow():
    s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])))
    w.ligate_adapter(s[1][0],hairpin_ligation_geometry(A1),s[1][1],note=s[1][2])
    w.ligate_adapter(s[2][0],hairpin_ligation_geometry(A2),s[2][1],note=s[2][2])
    w.react(s[3][0],s[3][1],note=s[3][2]);return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("genomic_dna","agarose-embedded genomic DNA with native DSBs",
                    MaterialKind.MOLECULAR,primary=True)
    deproteinize=graph.input("deproteinization_condition","in-plug deproteinization condition",
                             MaterialKind.CONDITION)
    clean=graph.product("clean","deproteinized genomic DNA",MaterialKind.MOLECULAR)
    graph.react("deproteinize","Deproteinize genomic DNA in agarose",
                StepKind.DEPROTEINIZATION,
                ((dna,InputRole.TEMPLATE),(deproteinize,InputRole.CONDITION)),(clean,))
    repair=graph.input("repair_enzymes","DSB-end blunting enzymes",MaterialKind.CATALYST)
    blunt=graph.product("blunted","blunted native DSB ends",MaterialKind.MOLECULAR)
    graph.react("blunt","Blunt native DSB ends",StepKind.END_REPAIR,
                ((clean,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(blunt,))
    tailing=graph.input("tailing_polymerase","dA-tailing polymerase",MaterialKind.CATALYST)
    tailed1=graph.product("tailed1","dA-tailed captured DSB ends",MaterialKind.MOLECULAR)
    graph.react("tail_break","Add 3-prime dA at processed break ends",StepKind.TAILING,
                ((blunt,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed1,))
    a1=graph.input("a1","biotinylated P5 USER hairpin adaptor 1",MaterialKind.ADAPTER)
    ligase1=graph.input("ligase1","P5 hairpin ligase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","P5-hairpin-tagged DSB ends",MaterialKind.MOLECULAR,
                         state=captured_break_scene())
    graph.react("ligate_a1","Ligate biotinylated P5 hairpin adaptor 1",
                StepKind.LIGATION,
                ((tailed1,InputRole.TEMPLATE),(a1,InputRole.STRUCTURAL),
                 (ligase1,InputRole.CATALYST)),(tagged,),geometry=hairpin_ligation_geometry(A1))
    gelase=graph.input("gelase","GELase for agarose digestion",MaterialKind.CATALYST)
    melting=graph.input("agarose_melting","agarose melting condition",MaterialKind.CONDITION)
    recovered=graph.product("recovered","DNA recovered from the agarose plug",MaterialKind.MOLECULAR)
    graph.react("recover_plug","Melt and digest agarose to recover tagged DNA",StepKind.EXTRACTION,
                ((tagged,InputRole.TEMPLATE),(gelase,InputRole.CATALYST),(melting,InputRole.CONDITION)),(recovered,))
    shear=graph.input("shear_condition","sonication condition",MaterialKind.CONDITION)
    fragments=graph.product("fragments","sheared P5-hairpin-tagged fragments",
                            MaterialKind.MOLECULAR)
    graph.react("shear","Shear distal genomic DNA",StepKind.FRAGMENTATION,
                ((recovered,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(fragments,))
    streptavidin=graph.input("streptavidin","streptavidin enrichment beads",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","captured P5-hairpin DSB fragments",
                           MaterialKind.MOLECULAR)
    untagged=graph.product("untagged","fragments lacking the biotinylated break end",
                           MaterialKind.WASTE)
    graph.react("capture","Enrich biotinylated break-end fragments",
                StepKind.AFFINITY_CAPTURE,
                ((fragments,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,untagged))
    graph.discard(untagged,"fragments lacking adaptor 1 are washed away")
    repaired2=graph.product("repaired2","distal-end-repaired captured fragments",
                            MaterialKind.MOLECULAR)
    graph.react("repair_distal","Repair the sonicated distal ends",StepKind.END_REPAIR,
                ((captured,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired2,))
    tailed2=graph.product("tailed2","distal-dA-tailed captured fragments",
                          MaterialKind.MOLECULAR)
    graph.react("tail_distal","Add 3-prime dA at distal fragment ends",StepKind.TAILING,
                ((repaired2,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed2,))
    a2=graph.input("a2","P7 USER hairpin adaptor 2",MaterialKind.ADAPTER)
    ligase2=graph.input("ligase2","P7 hairpin ligase",MaterialKind.CATALYST)
    double_hairpin=graph.product("double_hairpin","two-hairpin captured fragment",
                                 MaterialKind.MOLECULAR)
    graph.react("ligate_a2","Ligate P7 hairpin adaptor 2",StepKind.LIGATION,
                ((tailed2,InputRole.TEMPLATE),(a2,InputRole.STRUCTURAL),
                 (ligase2,InputRole.CATALYST)),(double_hairpin,),geometry=hairpin_ligation_geometry(A2))
    user=graph.input("user","USER enzyme",MaterialKind.CATALYST)
    opened=graph.product("opened","USER-opened PCR-addressable library precursor",
                         MaterialKind.MOLECULAR)
    loops=graph.product("loops","USER hairpin-cleavage products; recovery unassigned",MaterialKind.MOLECULAR)
    graph.react("open_hairpins","Open both uracil-containing hairpins",
                StepKind.DIGESTION,
                ((double_hairpin,InputRole.TEMPLATE),(user,InputRole.CATALYST)),
                (opened,loops))
    graph.mark_terminal(Terminal(loops.id,Disposition.UNRESOLVED,"USER opens the hairpins. Products are outside the retained adapter arms, but their subsequent physical removal is not established by cleavage itself."))
    forward=graph.input("pcr_forward","TruSeq PCR forward primer",MaterialKind.OLIGO)
    reverse=graph.input("pcr_reverse","indexed TruSeq PCR reverse primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("index_pcr",f"Amplify the opened END-seq library for {PCR_CYCLES} source cycles",StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The final TruSeq sites are exact; the opened two-hairpin precursor is not reconstructed as a complete duplex."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS), acquired_roles=("Read 1",), limitation="The original reports alternative50-base HiSeq2500 and75-base NextSeq500 single-end runs; exact index-cycle recipe is unassigned. Reference primer geometry does not assign actual loaded oligos."))
    graph.declare_reviewed_scope((
        "genomic_dna","deproteinization_condition","deproteinize","repair_enzymes",
        "blunt","tailing_polymerase","tail_break","a1","ligase1","ligate_a1",
        "gelase","agarose_melting","recover_plug","shear_condition","shear","streptavidin","capture","repair_distal",
        "tail_distal","a2","ligase2","ligate_a2","user","open_hairpins",
        "pcr_forward","pcr_reverse","pcr_polymerase","index_pcr","library"),SOURCE)
    return graph.require_publishable()
