"""TrAEL-seq native 3-prime-end capture."""
import illumina as il
from multibranch import illumina_branch
from protocol_extensions import native_three_prime_end
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState, Row, Scene, Workflow, feature
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from topology import HairpinTopology
TITLE="TrAEL-seq — native DNA 3-prime-end sequencing"
NOTES="01_trael-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1371/journal.pbio.3000886">Kara et al., <i>PLOS Biology</i> (2021)</a>.'
SUMMARY="A pre-adenylated, biotinylated hairpin adapter ligates selectively to native DNA 3-prime OH ends; fill-in, fragmentation and affinity selection retain the captured end for directional Illumina sequencing."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("TrAEL-seq","native-end-adjacent genomic DNA")
FINAL_CAPTION="Read orientation starts at genomic sequence adjacent to the captured native 3-prime end."
SEQUENCING_INTRO="Declared Illumina primers bind the completed TrAEL-seq library."
ADAPTER2 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCUUUUGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
ADAPTER1 = ("NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTU"
            "GCGCAGGCCATTGGCCTGCGCUACACTCTTTCCCTACACGACGCT")
ADAPTER1_TOPOLOGY = HairpinTopology(
    "TrAEL adapter 1", ADAPTER1, first_stem_start=18,
    second_stem_start=64, stem_length=23, colour="var(--c-r1)")
ADAPTER2_TOPOLOGY = HairpinTopology(
    "TrAEL adapter 2", ADAPTER2, first_stem_start=0,
    second_stem_start=35, stem_length=31, colour="var(--c-r2)")
def oligos():
    return ["<b>TrAEL adapter 1</b> — 5′-[rApp]NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTUGCGCAGGCCATTGGCC[BtndT]GCGCUACACTCTTTCCCTACACGACGCT-3′",f"<b>Adapter 2</b> — 5′-[Phos]{ADAPTER2[:-1]}*{ADAPTER2[-1]}-3′"]
def _adapter1_geometry():
    sc=Scene(); sc.strand("adapter 1",[
        seg("UMI","N"*8,"umi",placeholder=True,feature=feature("umi","umi","random")),
        seg("5′ adapter prefix","AGATCGGAAG"),
        seg("hairpin stem arm 1","AGCGTCGTGTAGGGAAAGAGTGT","r1"),
        seg("uracil/biotin loop","UGCGCAGGCCATTGGCCTGCGCU","w1"),
        seg("hairpin stem arm 2","ACACTCTTTCCCTACACGACGCT","r1")],label="published TrAEL adapter 1")
    sc.ligation_end("adapter 1","UMI",kind="single-strand",polarity="5-prime",owner="adapter",text="5′-App single-strand end")
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

def protocol_graph():
    source="Defining source: https://doi.org/10.1371/journal.pbio.3000886"
    graph=ProtocolGraph("TrAEL-seq",authority=source)
    dna=graph.input("native_dna","DNA bearing native 3-prime-OH ends",MaterialKind.BIOLOGICAL,
                    state=native_three_prime_end(),primary=True)
    adapter1=graph.input("adapter1","pre-adenylated biotinylated hairpin adapter 1",
                         MaterialKind.ADAPTER)
    rna_ligase=graph.input("rna_ligase","RNA ligase",MaterialKind.CATALYST)
    tagged=graph.product("adapter1_tagged","adapter-1-tagged native DNA end",MaterialKind.MOLECULAR)
    graph.react("adapter1_ligation","Ligate pre-adenylated adapter 1 to native 3-prime OH ends",
                StepKind.LIGATION,
                ((dna,InputRole.TEMPLATE),(adapter1,InputRole.STRUCTURAL),
                 (rna_ligase,InputRole.CATALYST)),(tagged,),geometry=_adapter1_geometry())
    bst=graph.input("bst","Bst strand-displacing polymerase",MaterialKind.CATALYST)
    filled=graph.product("filled_hairpin","hairpin-filled tagged duplex",MaterialKind.MOLECULAR)
    graph.react("hairpin_fill","Fill the adapter hairpin and displace its short arm",
                StepKind.EXTENSION,
                ((tagged,InputRole.TEMPLATE),(bst,InputRole.CATALYST)),(filled,))
    sonication=graph.input("sonication","sonication condition",MaterialKind.CONDITION)
    fragmented=graph.product("fragmented","sonicated tagged duplex",MaterialKind.MOLECULAR)
    graph.react("sonicate","Sonicate the tagged duplex",StepKind.FRAGMENTATION,
                ((filled,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),(fragmented,))
    repair=graph.input("end_repair_mix","blunt-end repair enzymes",MaterialKind.CATALYST)
    blunt=graph.product("blunt_fragment","blunt-ended tagged fragment",MaterialKind.MOLECULAR)
    graph.react("end_repair","Blunt end-repair the sonicated fragment",StepKind.END_REPAIR,
                ((fragmented,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(blunt,))
    adapter2=graph.input("adapter2","phosphorylated hairpin adapter 2",MaterialKind.ADAPTER)
    ligase=graph.input("dna_ligase","DNA ligase",MaterialKind.CATALYST)
    dual=graph.product("dual_adapter_fragment","two-adapter TrAEL molecule",MaterialKind.MOLECULAR)
    graph.react("adapter2_ligation","Ligate adapter 2 to the repaired end",StepKind.LIGATION,
                ((blunt,InputRole.TEMPLATE),(adapter2,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(dual,),geometry=_adapter2_geometry())
    beads=graph.input("streptavidin_beads","streptavidin capture beads",MaterialKind.SELECTION)
    captured=graph.product("captured","biotin-captured TrAEL molecules",MaterialKind.MOLECULAR)
    waste=graph.product("uncaptured","uncaptured DNA fragments",MaterialKind.WASTE)
    graph.react("capture","Capture molecules descending from a native 3-prime end",
                StepKind.AFFINITY_CAPTURE,
                ((dual,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(captured,waste))
    graph.discard(waste,"unbound molecules are removed during streptavidin washes")
    user=graph.input("user_enzyme","USER enzyme",MaterialKind.CATALYST)
    opened=graph.product("user_opened","USER-opened captured fragments",MaterialKind.MOLECULAR)
    graph.react("user_open","Open uracil-containing adapter sites",StepKind.DIGESTION,
                ((captured,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(opened,))
    forward=graph.input("pcr_forward","NEBNext Universal PCR primer",MaterialKind.OLIGO,
                        source=source+"; named NEB E7335 oligo set")
    reverse=graph.input("pcr_reverse","NEBNext indexed PCR primer",MaterialKind.OLIGO,
                        source=source+"; named NEB E7335 oligo set")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("indexed_pcr","Amplify captured molecules and complete Illumina adapters",
                StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    FINAL_LIBRARY,
                    PCRPrimer("NEBNext Universal",il.NEBNEXT_UNIVERSAL_PRIMER),
                    PCRPrimer("NEBNext indexed primer",il.nebnext_i7_primer("N"*8)),
                    statement="The paper names the NEBNext Universal and Index primers from E7335; their canonical kit sequences are centralized in illumina.py."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
    graph.declare_reviewed_scope((
        "native_dna","adapter1","rna_ligase","adapter1_ligation","bst","hairpin_fill",
        "sonication","sonicate","end_repair_mix","end_repair","adapter2","dna_ligase",
        "adapter2_ligation","streptavidin_beads","capture","user_enzyme","user_open",
        "pcr_forward","pcr_reverse","pcr_polymerase","indexed_pcr","library"),source)
    return graph.require_publishable()
