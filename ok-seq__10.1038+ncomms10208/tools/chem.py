"""Human OK-seq strand-preserving Okazaki-fragment library."""
import illumina as il
import seqprimers as sp
from protocol_extensions import native_three_prime_end
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      oligo, revcomp)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
TITLE="OK-seq — Okazaki-fragment sequencing"
NOTES="01_ok-seq.html"
SOURCE='Defining human method: <a href="https://doi.org/10.1038/ncomms10208">Petryk et al., <i>Nature Communications</i> (2016)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-022-00793-5">Nat. Protocols (2023)</a>.'
SUMMARY="EdU-labelled, short single-stranded Okazaki fragments are size-purified, click-biotinylated and ligated directionally to random-overhang adapters so read strand reports replication-fork direction."
A1_WATSON="ACACTCTTTCCCTACACGACGCTCTTCC"
A1_CRICK="NNNNNNGGAAGAGCGTCGTGTAGGGAAAGAGTG"
A2_WATSON="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
A2_CRICK="TGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNN"
PCR_FORWARD="ATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCC"
PCR_REVERSE="CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
CUSTOM_READ1=A1_WATSON
PRE_PCR=Construct([Segment("A1 Watson",A1_WATSON,"r1"),
                   Segment("Okazaki fragment","X"*36,placeholder=True),
                   Segment("A2 Watson",A2_WATSON)],name="adapter-ligated OK-seq fragment")
FINAL_LIBRARY=Construct([
    Segment("truncated P5",il.P5[1:],"p5"),Segment("A1 remainder",A1_WATSON[4:],"r1"),
    Segment("Okazaki fragment","X"*36,placeholder=True),
    Segment("completed A2 / Index-1 site",revcomp(il.TRUSEQ_READ2),"r2"),
    Segment("six-base i7 reverse complement","N"*6,"cbc",placeholder=True),
    Segment("P7 reverse complement",il.P7_RC,"p7")],name="OK-seq library")
SEQ_PRIMERS=(sp.custom("Read 1","OK-seq custom Read 1",CUSTOM_READ1,
                       "Petryk et al. 2016 supplementary methods"),sp.TRUSEQ["I1"])
_problems=sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)"))
if _problems: raise ValueError("invalid OK-seq endpoint: "+"; ".join(_problems))
PCR_GEOMETRY=PCRGeometry.exact(PRE_PCR,
    PCRPrimer("forward PCR primer",PCR_FORWARD,anneal_nt=len(A1_WATSON)),
    PCRPrimer("indexed reverse PCR primer",PCR_REVERSE,anneal_nt=len(A2_WATSON)),
    statement="The detailed protocol prints both complete PCR primers; their 3′ arms are the A1 Watson and reverse-complementary A2 Watson adapter strands.")
FINAL_CAPTION="Adapter orientation preserves whether the captured Okazaki fragment came from the Watson or Crick strand."
SEQUENCING_INTRO="The source specifies the A1 Watson oligo as a custom sequencing primer; the six-base sample index is read with the standard TruSeq Index 1 primer."
ACQUIRED_ROLES=("Read 1","Index 1 (i7)")
READ_LENGTHS={"Read 1":51,"Index 1 (i7)":6}
READ_USES={"Read 1":"include","Index 1 (i7)":"demultiplex"}
def oligos():
    return [oligo("A1 Watson",(Segment("A1 Watson",A1_WATSON),)),
            oligo("A1 Crick",(Segment("random overhang","N"*6,placeholder=True),Segment("paired arm",A1_CRICK[6:])),),
            oligo("A2 Watson",(Segment("A2 Watson",A2_WATSON),),mods="/5Phos/",three="/3ddC/-3′"),
            oligo("A2 Crick",(Segment("paired arm",A2_CRICK[:-6]),Segment("random overhang","N"*6,placeholder=True)),three="/3ddC/-3′"),
            oligo("forward PCR primer",(Segment("forward PCR primer",PCR_FORWARD),)),
            oligo("indexed reverse PCR primer",(Segment("indexed reverse PCR primer",PCR_REVERSE),)),
            oligo("custom Read 1 sequencing primer",(Segment("custom Read 1",CUSTOM_READ1),))]
def adapter_geometry():
    a1=Scene(); a1.strand("A1 Watson",[Segment("A1 Watson ligating strand",A1_WATSON)],label="A1 adapter")
    a1.anneal("A1 Crick",[Segment("A1 random annealing overhang","N"*6,placeholder=True),Segment("A1 paired arm",A1_CRICK[6:])],to="A1 Watson",pair=("A1 paired arm","A1 Watson ligating strand"),shift=1,label="A1 adapter",unpaired=("A1 random annealing overhang",))
    a1.overhang("A1 Crick","A1 random annealing overhang",polarity="5-prime",owner="adapter",text="6-nt random annealing overhang")
    a1.ligation_end("A1 Watson","A1 Watson ligating strand",kind="splinted",polarity="3-prime",owner="adapter",text="3′-OH joins the fragment 5′ phosphate")
    a2=Scene(); a2.strand("A2 Watson",[Segment("A2 Watson ligating strand",A2_WATSON)],label="A2 adapter",mod5="p",mod3="ddC")
    a2.anneal("A2 Crick",[Segment("A2 paired arm",A2_CRICK[:-6]),Segment("A2 random annealing overhang","N"*6,placeholder=True)],to="A2 Watson",pair=("A2 paired arm","A2 Watson ligating strand"),label="A2 adapter",mod3="ddC",unpaired=("A2 random annealing overhang",))
    a2.overhang("A2 Crick","A2 random annealing overhang",polarity="3-prime",owner="adapter",text="6-nt random annealing overhang")
    a2.ligation_end("A2 Watson","A2 Watson ligating strand",kind="splinted",polarity="5-prime",owner="adapter",text="5′ phosphate joins the fragment 3′-OH")
    return AdapterLigationGeometry.known("splinted",a1,a2,statement="The detailed OK-seq protocol prints both duplex adapters. Their six random bases anneal to each ssDNA terminus, positioning the opposite adapter strand for directional nick sealing; the random strand is ddC-blocked where printed.")
def sections():
    geometry=adapter_geometry()
    wf=Workflow(MolecularState("Cells in S phase",(Row(chunks=[("replicating DNA",None,False)]),)))
    wf.react("Pulse-label nascent DNA with EdU",[Row(chunks=[("lagging strand: short Okazaki fragment — EdU","w1",False)])],note="EdU marks newly synthesized fragments without requiring ligase-deficient cells.")
    wf.react("Denature and size-purify short single strands",[Row(chunks=[("genomic DNA → denature → sucrose gradient → <200-nt ssDNA",None,False)])],note="The physical selection enriches Okazaki fragments.")
    wf.react("Click biotin onto EdU and phosphorylate 5-prime ends",native_three_prime_end(captured="Okazaki-fragment end").rows(),note="Azide–alkyne chemistry creates the affinity handle; PNK makes ligatable 5-prime ends.")
    wf.ligate_adapter("Ligate A1 and A2 and capture",geometry,[Row(chunks=[("A1 ** Okazaki ssDNA ** A2 → streptavidin capture → PCR","me",False)])],note="Adapter orientation and single-end sequencing preserve the genomic strand signal.")
    return [(r.action,r.after.rows,r.note) for r in wf.reactions]

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/ncomms10208; detailed protocol: https://doi.org/10.1038/s41596-022-00793-5"
    graph=ProtocolGraph("OK-seq",authority=source)
    cells=graph.input("cells","cells in S phase",MaterialKind.BIOLOGICAL,primary=True)
    edu=graph.input("edu_pulse","EdU pulse-label condition",MaterialKind.CONDITION)
    labelled=graph.product("labelled_dna","EdU-labelled nascent DNA",MaterialKind.MOLECULAR)
    graph.react("pulse_label","Pulse-label nascent DNA with EdU",StepKind.CHEMICAL_CONVERSION,
                ((cells,InputRole.TEMPLATE),(edu,InputRole.CONDITION)),(labelled,))
    denature=graph.input("denaturation","DNA denaturation condition",MaterialKind.CONDITION)
    single=graph.product("single_strands","denatured single-stranded DNA",MaterialKind.MOLECULAR)
    graph.react("denature","Denature genomic DNA",StepKind.THERMAL_TREATMENT,
                ((labelled,InputRole.TEMPLATE),(denature,InputRole.CONDITION)),(single,))
    gradient=graph.input("size_gradient","less-than-200-nt sucrose-gradient selection",
                         MaterialKind.SELECTION)
    short=graph.product("okazaki_fragments","size-selected short Okazaki fragments",
                        MaterialKind.MOLECULAR)
    long=graph.product("long_dna","rejected long genomic DNA",MaterialKind.WASTE)
    graph.react("size_select","Select short single-stranded Okazaki fragments",StepKind.SELECTION,
                ((single,InputRole.TEMPLATE),(gradient,InputRole.SELECTION)),(short,long))
    graph.discard(long,"long genomic DNA is excluded by the sucrose-gradient size selection")
    click=graph.input("click_reagents","azide-alkyne biotin click reagents",MaterialKind.CONDITION)
    biotin=graph.product("biotinylated_fragments","click-biotinylated Okazaki fragments",
                         MaterialKind.MOLECULAR)
    graph.react("click_biotin","Click biotin onto EdU-labelled fragments",
                StepKind.CHEMICAL_CONVERSION,
                ((short,InputRole.TEMPLATE),(click,InputRole.CONDITION)),(biotin,))
    pnk=graph.input("pnk","polynucleotide kinase",MaterialKind.CATALYST)
    phosphorylated=graph.product("phosphorylated_fragments","5-prime-phosphorylated Okazaki fragments",
                                 MaterialKind.MOLECULAR)
    graph.react("phosphorylate","Phosphorylate fragment 5-prime ends",StepKind.PHOSPHORYLATION,
                ((biotin,InputRole.TEMPLATE),(pnk,InputRole.CATALYST)),(phosphorylated,))
    adapters=graph.input("directional_adapters","A1 and A2 random-overhang adapters",
                         MaterialKind.ADAPTER)
    ligase=graph.input("ligase","DNA ligase",MaterialKind.CATALYST)
    ligated=graph.product("ligated_fragments","directionally adapter-ligated Okazaki fragments",
                          MaterialKind.MOLECULAR)
    graph.react("adapter_ligation","Directionally ligate A1 and A2",StepKind.LIGATION,
                ((phosphorylated,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=adapter_geometry())
    beads=graph.input("streptavidin_beads","streptavidin capture beads",MaterialKind.SELECTION)
    captured=graph.product("captured_fragments","biotin-captured adapter-ligated fragments",
                           MaterialKind.MOLECULAR)
    uncaptured=graph.product("uncaptured_fragments","uncaptured DNA",MaterialKind.WASTE)
    graph.react("capture","Capture EdU-derived biotinylated molecules",StepKind.AFFINITY_CAPTURE,
                ((ligated,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(captured,uncaptured))
    graph.discard(uncaptured,"unbound DNA is removed during streptavidin washes")
    forward=graph.input("pcr_forward","OK-seq forward PCR primer",MaterialKind.OLIGO,
                        inferred=True,source=source+" (PCR-primer bases not represented on this page)")
    reverse=graph.input("pcr_reverse","indexed OK-seq reverse PCR primer",MaterialKind.OLIGO,
                        inferred=True,source=source+" (PCR-primer bases not represented on this page)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Amplify captured molecules and complete Illumina adapters",
                StepKind.PCR,
                ((captured,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCR_GEOMETRY)
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=ACQUIRED_ROLES, read_lengths=READ_LENGTHS,
        read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
    graph.declare_reviewed_scope((
        "cells","edu_pulse","pulse_label","denaturation","denature","size_gradient",
        "size_select","click_reagents","click_biotin","pnk","phosphorylate",
        "directional_adapters","ligase","adapter_ligation","streptavidin_beads",
        "capture","pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
