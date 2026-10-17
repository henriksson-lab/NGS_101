"""Human OK-seq strand-preserving Okazaki-fragment library."""
import illumina as il
import seqprimers as sp
from protocol_extensions import native_three_prime_end
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      feature, oligo, revcomp)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
TITLE="OK-seq — Okazaki-fragment sequencing"
NOTES="01_ok-seq.html"
SOURCE='Defining human method: <a href="https://doi.org/10.1038/ncomms10208">Petryk et al., <i>Nature Communications</i> (2016)</a>, Methods. A later detailed protocol (<a href="https://doi.org/10.1038/s41596-022-00793-5">Nat. Protocols 2023</a>) has not been inspected for this page.'
CAVEAT=("Published-source note: Methods call the forward primer &quot;PE PCR Primer 1.0 truncated by 5 bases&quot; "
        "(consistent with the 3′ GATCT removed, as for A1watson), but the printed sequence also begins ATGATAC… "
        "without the leading A of P5. The model uses the printed 52-nt primer, so the drawn library carries a "
        "P5 lacking its first base; whether the printed sequence or the description is the intended reagent "
        "remains unresolved. The six-base i7 index (six asterisks in the printed reverse primer) is drawn, "
        "but the acquired index-read cycle count is not printed; only the 51-base single read is stated.")
SUMMARY="EdU-labelled, short single-stranded Okazaki fragments are size-purified, click-biotinylated and ligated directionally to random-overhang adapters so read strand reports replication-fork direction."
A1_WATSON="ACACTCTTTCCCTACACGACGCTCTTCC"
A1_CRICK="NNNNNNGGAAGAGCGTCGTGTAGGGAAAGAGTG"
A2_WATSON="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
A2_CRICK="TGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNN"
PCR_FORWARD="ATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCC"
PCR_REVERSE="CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
CUSTOM_READ1=A1_WATSON
I7_FEATURE=feature("sample_index_i7","sample_index","unknown")
# Printed as P7 arm, six asterisks ("Illumina index"), then the TruSeq index arm.
REVERSE_PRIMER=Construct([Segment("P7 arm",PCR_REVERSE[:24]),
                          Segment("six-base Illumina index","N"*6,"cbc",placeholder=True,feature=I7_FEATURE),
                          Segment("TruSeq index arm",PCR_REVERSE[30:])],name="TruSeq Index reverse PCR primer")
ADAPTER_OLIGOS=(A1_WATSON,A1_CRICK,A2_WATSON,A2_CRICK)
PRE_PCR=Construct([Segment("A1 Watson",A1_WATSON,"r1"),
                   Segment("Okazaki fragment","X"*36,placeholder=True),
                   Segment("A2 Watson",A2_WATSON)],name="adapter-ligated OK-seq fragment")
FINAL_LIBRARY=Construct([
    Segment("truncated P5",il.P5[1:],"p5"),Segment("A1 remainder",A1_WATSON[4:],"r1"),
    Segment("Okazaki fragment","X"*36,placeholder=True),
    Segment("completed A2 / Index-1 site",revcomp(il.TRUSEQ_READ2),"r2"),
    Segment("six-base i7 reverse complement","N"*6,"cbc",placeholder=True,feature=I7_FEATURE),
    Segment("P7 reverse complement",il.P7_RC,"p7")],name="OK-seq library")
SEQ_PRIMERS=(sp.custom("Read 1","OK-seq custom Read 1",CUSTOM_READ1,
                       "Petryk et al. 2016 Methods, Library amplification and sequencing"),sp.TRUSEQ["I1"])
_problems=sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)"))
if _problems: raise ValueError("invalid OK-seq endpoint: "+"; ".join(_problems))
PCR_GEOMETRY=PCRGeometry.exact(PRE_PCR,
    PCRPrimer("forward PCR primer",PCR_FORWARD,anneal_nt=len(A1_WATSON)),
    PCRPrimer("indexed reverse PCR primer",PCR_REVERSE,anneal_nt=len(A2_WATSON)),
    statement="Petryk et al. 2016 Methods print both PCR primers (12 cycles, Phusion); their 3′ arms are the A1 Watson and reverse-complementary A2 Watson adapter strands.")
FINAL_CAPTION="Adapter orientation preserves whether the captured Okazaki fragment came from the Watson or Crick strand."
SEQUENCING_INTRO="The source specifies the A1 Watson oligo as a custom sequencing primer; the six-base sample index is read with the standard TruSeq Index 1 primer."
ACQUIRED_ROLES=("Read 1","Index 1 (i7)")
READ_LENGTHS={"Read 1":51}
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
    return AdapterLigationGeometry.known("splinted",a1,a2,statement="Petryk et al. 2016 Methods print both duplex adapters (annealed 95 °C 5 s, cooled 0.1 °C/s to 4 °C). Their six random bases anneal to each ssDNA terminus, positioning the opposite adapter strand for directional nick sealing; the random strand is ddC-blocked where printed.")
def sections():
    geometry=adapter_geometry()
    wf=Workflow(MolecularState("Cells in S phase",(Row(chunks=[("replicating DNA",None,False)]),)))
    wf.react("Pulse-label nascent DNA with EdU",[Row(chunks=[("lagging strand: short Okazaki fragment — EdU","w1",False)])],note="EdU marks newly synthesized fragments without requiring ligase-deficient cells.")
    wf.react("Denature and size-purify short single strands",[Row(chunks=[("genomic DNA → denature → sucrose gradient → <200-nt ssDNA",None,False)])],note="The physical selection enriches Okazaki fragments.")
    wf.react("Click biotin onto EdU and phosphorylate 5-prime ends",native_three_prime_end(captured="Okazaki-fragment end").rows(),note="Azide–alkyne chemistry creates the affinity handle; PNK makes ligatable 5-prime ends.")
    wf.ligate_adapter("Ligate A1 and A2 and capture",geometry,[Row(chunks=[("A1 ** Okazaki ssDNA ** A2 → streptavidin capture → PCR","me",False)])],note="Adapter orientation and single-end sequencing preserve the genomic strand signal.")
    return [(r.action,r.after.rows,r.note) for r in wf.reactions]

def protocol_graph():
    source="Defining source: Petryk et al. 2016, https://doi.org/10.1038/ncomms10208, Methods"
    graph=ProtocolGraph("OK-seq",authority=source)
    def cleanup(prefix,template,reagent,kept_label,label,kind=StepKind.CLEANUP,waste_label=None,reason="removed by the cleanup"):
        agent=graph.input(prefix+"_reagent",reagent,MaterialKind.SELECTION)
        kept=graph.product(prefix+"_kept",kept_label,MaterialKind.MOLECULAR)
        waste=graph.product(prefix+"_waste",waste_label or reagent+" discarded fraction",MaterialKind.WASTE)
        graph.react(prefix,label,kind,((template,InputRole.TEMPLATE),(agent,InputRole.SELECTION)),(kept,waste))
        graph.discard(waste,reason)
        reviewed.extend((agent.id,prefix))
        return kept
    reviewed=[]
    cells=graph.input("cells","asynchronous HeLa or GM06990 cells",MaterialKind.BIOLOGICAL,primary=True)
    edu=graph.input("edu_pulse","20 µM EdU for 2 min",MaterialKind.CONDITION)
    labelled=graph.product("labelled_cells","cells with EdU-labelled nascent DNA",MaterialKind.MOLECULAR)
    graph.react("pulse_label","Pulse-label nascent DNA with EdU",StepKind.CHEMICAL_CONVERSION,
                ((cells,InputRole.TEMPLATE),(edu,InputRole.CONDITION)),(labelled,))
    extraction=graph.input("dna_extraction","proteinase K/phenol extraction and ethanol precipitation",MaterialKind.CONDITION)
    genomic=graph.product("genomic_dna","genomic DNA including EdU-labelled Okazaki fragments",MaterialKind.MOLECULAR)
    graph.react("extract_dna","Extract genomic DNA",StepKind.DEPROTEINIZATION,
                ((labelled,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(genomic,))
    denature=graph.input("denaturation","94 °C for 10 min, chill on ice",MaterialKind.CONDITION)
    single=graph.product("single_strands","denatured single-stranded DNA",MaterialKind.MOLECULAR)
    graph.react("denature","Denature genomic DNA",StepKind.THERMAL_TREATMENT,
                ((genomic,InputRole.TEMPLATE),(denature,InputRole.CONDITION)),(single,))
    reviewed.extend(("cells","edu_pulse","pulse_label","dna_extraction","extract_dna","denaturation","denature"))
    short=cleanup("size_select",single,"5–30% sucrose gradient (SW28, 17 h, 26,000 g); <200-nt fractions, Amicon 10k concentration",
                  "size-selected short single strands (<200 nt)","Select short single strands on a sucrose gradient",
                  kind=StepKind.SELECTION,waste_label="longer genomic DNA strands",
                  reason="long strands are excluded by sucrose-gradient fractionation")
    click=graph.input("click_reagents","2 mM CuSO4, 1 mM biotin-TEG-azide, 10 mM sodium ascorbate, 45 min RT in the dark; ethanol precipitation",MaterialKind.CONDITION)
    biotin=graph.product("biotinylated_fragments","click-biotinylated Okazaki fragments",MaterialKind.MOLECULAR)
    graph.react("click_biotin","Click biotin onto EdU-labelled fragments",StepKind.CHEMICAL_CONVERSION,
                ((short,InputRole.TEMPLATE),(click,InputRole.CONDITION)),(biotin,))
    naoh=graph.input("alkaline_hydrolysis","0.25 M NaOH-equivalent addition (10 µL 2.5 M NaOH), 30 min 37 °C; acetic-acid neutralization",MaterialKind.CONDITION)
    hydrolysed=graph.product("rna_hydrolysed","Okazaki-fragment DNA after RNA-primer hydrolysis",MaterialKind.MOLECULAR)
    graph.react("hydrolyse_rna","Hydrolyse RNA primers with alkali",StepKind.CHEMICAL_CONVERSION,
                ((biotin,InputRole.TEMPLATE),(naoh,InputRole.CONDITION)),(hydrolysed,))
    reviewed.extend(("click_reagents","click_biotin","alkaline_hydrolysis","hydrolyse_rna"))
    desalted=cleanup("p30_cleanup",hydrolysed,"Micro Bio-Spin P-30 column","column-purified DNA strands",
                     "Purify on Bio-Gel P-30 spin columns",waste_label="retained small molecules and hydrolysed ribonucleotides",
                     reason="retained by the gel-filtration column")
    pnk=graph.input("pnk","30 U T4 polynucleotide kinase (Fermentas), 30 min 37 °C; 75 °C 10 min inactivation",MaterialKind.CATALYST)
    phosphorylated=graph.product("phosphorylated_fragments","5-prime-phosphorylated Okazaki fragments",
                                 MaterialKind.MOLECULAR)
    graph.react("phosphorylate","Phosphorylate fragment 5-prime ends",StepKind.PHOSPHORYLATION,
                ((desalted,InputRole.TEMPLATE),(pnk,InputRole.CATALYST)),(phosphorylated,))
    adapters=graph.input("directional_adapters","A1 and A2 random-overhang adapters (80 pM each; 65 °C 10 min, 16 °C 5 min annealing)",
                         MaterialKind.ADAPTER,state=ADAPTER_OLIGOS)
    ligase=graph.input("ligase","4 Weiss U T4 DNA ligase (Fermentas), 500 mM betaine, 5% PEG8000, 16 °C overnight",MaterialKind.CATALYST)
    ligated=graph.product("ligated_fragments","directionally adapter-ligated Okazaki fragments",
                          MaterialKind.MOLECULAR)
    graph.react("adapter_ligation","Directionally ligate A1 and A2",StepKind.LIGATION,
                ((phosphorylated,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=adapter_geometry())
    reviewed.extend(("pnk","phosphorylate","directional_adapters","ligase","adapter_ligation"))
    captured=cleanup("capture",ligated,"200 µg Dynabeads MyOne Streptavidin T1; 1 M NaCl/Tween washes",
                     "bead-bound biotinylated fragments","Capture biotinylated fragments and wash away unligated adaptors",
                     kind=StepKind.AFFINITY_CAPTURE,waste_label="unbound DNA and unligated adaptors",
                     reason="removed in the streptavidin bead washes")
    adapters2=graph.input("directional_adapters_repeat","A1 and A2 adapters (40 pM each) for the repeated on-bead ligation",
                          MaterialKind.ADAPTER,state=ADAPTER_OLIGOS)
    ligase2=graph.input("ligase_repeat","T4 DNA ligase, repeated ligation conditions",MaterialKind.CATALYST)
    religated=graph.product("religated_fragments","bead-bound fragments after a second annealing/ligation",MaterialKind.MOLECULAR)
    graph.react("adapter_ligation_repeat","Repeat adapter annealing and ligation on 20 µL bead suspension",StepKind.LIGATION,
                ((captured,InputRole.TEMPLATE),(adapters2,InputRole.STRUCTURAL),(ligase2,InputRole.CATALYST)),
                (religated,),geometry=adapter_geometry())
    reviewed.extend(("directional_adapters_repeat","ligase_repeat","adapter_ligation_repeat"))
    washed=cleanup("rewash",religated,"repeat bead washes; resuspend in 20 µL 10 mM Tris-HCl",
                   "washed bead-bound adapter-ligated fragments","Wash beads as before",
                   waste_label="unligated adaptors from the repeated ligation",reason="removed in the bead washes")
    forward=graph.input("pcr_forward","OK-seq forward PCR primer (printed)",MaterialKind.OLIGO,state=PCR_FORWARD,
                        source=source+", Library amplification and sequencing")
    reverse=graph.input("pcr_reverse","TruSeq Index reverse PCR primer (six-base index printed as asterisks)",MaterialKind.OLIGO,
                        state=REVERSE_PRIMER,source=source+", Library amplification and sequencing")
    polymerase=graph.input("pcr_polymerase","Phusion polymerase (NEB), 12 cycles, 10 µL bead suspension",MaterialKind.CATALYST)
    amplified=graph.product("amplified","amplified OK-seq library before size selection",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    graph.react("library_pcr","Amplify captured molecules and complete Illumina adapters",
                StepKind.PCR,
                ((washed,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),
                geometry=PCR_GEOMETRY)
    reviewed.extend(("pcr_forward","pcr_reverse","pcr_polymerase","library_pcr"))
    columned=cleanup("minelute",amplified,"Qiagen MinElute column","column-purified PCR products",
                     "Purify PCR products on MinElute columns",waste_label="column flow-through",
                     reason="unbound material removed by the column")
    gel=graph.input("gel_selection","4% small-fragment agarose gel, 150–400 bp excision, Qiagen gel extraction",MaterialKind.SELECTION)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    off=graph.product("gel_excluded","products outside 150–400 bp",MaterialKind.WASTE)
    graph.react("gel_select","Gel-select 150–400 bp fragments",StepKind.SELECTION,
                ((columned,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(library,off))
    graph.discard(off,"left in the agarose gel")
    reviewed.extend(("gel_selection","gel_select","library"))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=ACQUIRED_ROLES, read_lengths=READ_LENGTHS,
        read_uses=READ_USES))
    graph.declare_reviewed_scope(tuple(reviewed),source)
    return graph.require_publishable()
