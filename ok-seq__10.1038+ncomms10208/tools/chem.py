"""Human OK-seq strand-preserving Okazaki-fragment library."""
from multibranch import illumina_branch
from protocol_extensions import native_three_prime_end
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Scene,Segment,Workflow
TITLE="OK-seq — Okazaki-fragment sequencing"
NOTES="01_ok-seq.html"
SOURCE='Defining human method: <a href="https://doi.org/10.1038/ncomms10208">Petryk et al., <i>Nature Communications</i> (2016)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-022-00793-5">Nat. Protocols (2023)</a>.'
SUMMARY="EdU-labelled, short single-stranded Okazaki fragments are size-purified, click-biotinylated and ligated directionally to random-overhang adapters so read strand reports replication-fork direction."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("OK-seq","Okazaki fragment")
FINAL_CAPTION="Adapter orientation preserves whether the captured Okazaki fragment came from the Watson or Crick strand."
SEQUENCING_INTRO="Declared Illumina primers bind the completed OK-seq library."
def oligos():
    return ["<b>A1 Watson</b> — 5′-ACACTCTTTCCCTACACGACGCTCTTCC-3′","<b>A1 Crick</b> — 5′-NNNNNNGGAAGAGCGTCGTGTAGGGAAAGAGTG-3′","<b>A2 Watson</b> — 5′-[Phos]-AGATCGGAAGAGCACACGTCTGAACTCCAGTCA[ddC]-3′","<b>A2 Crick</b> — 5′-TGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNN[ddC]-3′"]
def sections():
    a1=Scene(); a1.strand("A1 Watson",[Segment("A1 Watson ligating strand","ACACTCTTTCCCTACACGACGCTCTTCC")],label="A1 adapter")
    a1.anneal("A1 Crick",[Segment("A1 random annealing overhang","N"*6,placeholder=True),Segment("A1 paired arm","GGAAGAGCGTCGTGTAGGGAAAGAGTG")],to="A1 Watson",pair=("A1 paired arm","A1 Watson ligating strand"),shift=1,label="A1 adapter",unpaired=("A1 random annealing overhang",))
    a1.overhang("A1 Crick","A1 random annealing overhang",polarity="5-prime",owner="adapter",text="6-nt random annealing overhang")
    a1.ligation_end("A1 Watson","A1 Watson ligating strand",kind="splinted",polarity="3-prime",owner="adapter",text="3′-OH joins the fragment 5′ phosphate")
    a2=Scene(); a2.strand("A2 Watson",[Segment("A2 Watson ligating strand","AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")],label="A2 adapter",mod5="p")
    a2.anneal("A2 Crick",[Segment("A2 paired arm","TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"),Segment("A2 random annealing overhang","N"*6,placeholder=True)],to="A2 Watson",pair=("A2 paired arm","A2 Watson ligating strand"),label="A2 adapter",mod3="ddC",unpaired=("A2 random annealing overhang",))
    a2.overhang("A2 Crick","A2 random annealing overhang",polarity="3-prime",owner="adapter",text="6-nt random annealing overhang")
    a2.ligation_end("A2 Watson","A2 Watson ligating strand",kind="splinted",polarity="5-prime",owner="adapter",text="5′ phosphate joins the fragment 3′-OH")
    geometry=AdapterLigationGeometry.known("splinted",a1,a2,statement="The detailed OK-seq protocol prints both duplex adapters. Their six random bases anneal to each ssDNA terminus, positioning the opposite adapter strand for directional nick sealing; the random strand is ddC-blocked where printed.")
    wf=Workflow(MolecularState("Cells in S phase",(Row(chunks=[("replicating DNA",None,False)]),)))
    wf.react("Pulse-label nascent DNA with EdU",[Row(chunks=[("lagging strand: short Okazaki fragment — EdU","w1",False)])],note="EdU marks newly synthesized fragments without requiring ligase-deficient cells.")
    wf.react("Denature and size-purify short single strands",[Row(chunks=[("genomic DNA → denature → sucrose gradient → <200-nt ssDNA",None,False)])],note="The physical selection enriches Okazaki fragments.")
    wf.react("Click biotin onto EdU and phosphorylate 5-prime ends",native_three_prime_end(captured="Okazaki-fragment end").rows(),note="Azide–alkyne chemistry creates the affinity handle; PNK makes ligatable 5-prime ends.")
    wf.ligate_adapter("Ligate A1 and A2 and capture",geometry,[Row(chunks=[("A1 ** Okazaki ssDNA ** A2 → streptavidin capture → PCR","me",False)])],note="Adapter orientation and single-end sequencing preserve the genomic strand signal.")
    return [(r.action,r.after.rows,r.note) for r in wf.reactions]
