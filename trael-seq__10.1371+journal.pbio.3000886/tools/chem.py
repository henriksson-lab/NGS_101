"""TrAEL-seq native 3-prime-end capture."""
from multibranch import illumina_branch
from protocol_extensions import native_three_prime_end
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState, Row, Scene, Workflow
TITLE="TrAEL-seq — native DNA 3-prime-end sequencing"
NOTES="01_trael-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1371/journal.pbio.3000886">Kara et al., <i>PLOS Biology</i> (2021)</a>.'
SUMMARY="A pre-adenylated, biotinylated hairpin adapter ligates selectively to native DNA 3-prime OH ends; fill-in, fragmentation and affinity selection retain the captured end for directional Illumina sequencing."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("TrAEL-seq","native-end-adjacent genomic DNA")
FINAL_CAPTION="Read orientation starts at genomic sequence adjacent to the captured native 3-prime end."
SEQUENCING_INTRO="Declared Illumina primers bind the completed TrAEL-seq library."
def oligos():
    return ["<b>TrAEL adapter 1</b> — 5′-[Phos]NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTU…[BtndT]…-3′","<b>Adapter 2</b> — 5′-[Phos]GATCGGAAGAGCACACGTCTGAACTCCAGTCUUUUGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3′"]
def _adapter1_geometry():
    sc=Scene(); sc.strand("adapter 1",[seg("5′ pre-adenylated ligating end","N"*8,placeholder=True),seg("hairpin arm with internal biotin","X"*26,placeholder=True)],label="published TrAEL adapter 1")
    sc.ligation_end("adapter 1","5′ pre-adenylated ligating end",kind="single-strand",polarity="5-prime",owner="adapter",text="5′-App single-strand end")
    return AdapterLigationGeometry.known("single-strand",sc,statement="The protocol prints adapter 1 as a 5′-adenylated hairpin and joins that activated 5′ end to a native DNA 3′-OH without ATP.")

def _adapter2_geometry():
    sc=Scene.duplex([seg("adapter 2 paired ligating end","X"*20,placeholder=True)],label="published adapter 2 hairpin")
    sc.ligation_end("top","adapter 2 paired ligating end",kind="blunt",owner="adapter",text="blunt hairpin end")
    return AdapterLigationGeometry.known("blunt",sc,statement="The protocol prints phosphorylated adapter 2 and ligates it after sonication and blunt end repair.")

def workflow():
    start=native_three_prime_end()
    w=Workflow(MolecularState("native DNA 3′ end",tuple(start.rows())))
    first=Scene(); first.strand("captured molecule",[seg("native-end genomic DNA","X"*30,placeholder=True),seg("adapter 1","X"*28,placeholder=True)],label="adapter-1-ligated end"); first.junction("captured molecule","native-end genomic DNA","adapter 1","RNA-ligase junction")
    w.ligate_adapter("Ligate pre-adenylated adapter 1",_adapter1_geometry(),first.rows(),name="adapter-1-tagged end",note="Its random bases identify molecules and internal biotin enables capture.")
    w.react("Fill the adapter hairpin and displace its short arm",[Row(chunks=[("captured 3′ end—adapter hairpin → Bst fill/displacement → tagged duplex","r1",False)])],name="tagged duplex",note="Polymerase stabilizes the captured junction.")
    w.react("Sonicate and blunt end-repair",[Row(chunks=[("tagged duplex → sonicated fragment with repaired blunt end","me",False)])],name="blunt-ended tagged fragment",note="End repair creates the substrate for adapter 2.")
    second=Scene.duplex([seg("tagged genomic fragment","X"*34,placeholder=True),seg("adapter 2","X"*20,placeholder=True)],label="two-adapter fragment"); second.junction("top","tagged genomic fragment","adapter 2","adapter-2 ligation")
    w.ligate_adapter("Ligate adapter 2",_adapter2_geometry(),second.rows(),name="two-adapter TrAEL molecule",note="The second adapter supplies the opposite amplification end.")
    w.react("Capture biotin and open USER sites",[Row(chunks=[("streptavidin ● biotin-tagged fragment → USER release → indexed PCR","w1",False)])],name="indexed TrAEL library",note="Affinity selection retains molecules descending from a native 3′ end.")
    return w
