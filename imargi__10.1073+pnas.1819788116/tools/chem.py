"""iMARGI in-situ RNA–DNA bridge ligation."""
from multibranch import illumina_branch
from interaction_capture import BivalentRNAChromatinBridge
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Scene,Segment,Workflow
from chemdraw import Construct
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
TITLE="iMARGI — in situ RNA–genome interaction mapping"
NOTES="01_imargi.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1819788116">Yan et al., <i>PNAS</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-019-0229-4">Nat. Protocols (2019)</a>.'
SUMMARY="A bivalent linker is ligated in situ first to fragmented chromatin-associated RNA and then to nearby genomic DNA; reverse transcription and Illumina preparation preserve the RNA–linker–DNA order."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("iMARGI","RNA–linker–DNA contact chimera",inferred_adapters=True)
FINAL_CAPTION="INFERRED — paired reads recover the RNA-derived and DNA-derived ends; the exact PCR landing geometry is not represented by this compact library view. The linker has two ligation-bias bases, not a UMI."
SEQUENCING_INTRO="Declared Illumina primers bind the completed iMARGI library."
CAVEAT=(
    "The author protocol prints the modified bivalent linker, cut oligo and PCR primers. "
    "This compact final-library view does not draw the strand-resolved BamHI split and PCR landing junctions. "
    "Published citation discrepancy: Yan et al. 2019 reference 10 cites "
    "<a href=\"https://doi.org/10.1016/j.cub.2017.01.068\">DOI 10.1016/j.cub.2017.01.068</a> "
    "for Sridhar et al. MARGI, but that DOI is the erratum; "
    "the <a href=\"https://doi.org/10.1016/j.cub.2017.01.011\">original article</a> "
    "has DOI 10.1016/j.cub.2017.01.011. The model uses the "
    "iMARGI detailed-protocol oligos. Original Sridhar 2017 MARGI Document S1 "
    "(Supplemental Experimental Procedures, lines 204–215) prints a top linker "
    "ending ...TCCGA/BiodT/CT without a C* mark; later Wu 2019 iMARGI Nature "
    "Protocols (JATS lines 275–285) prints ...TCC*GA/iBiodT/CT and defines C* "
    "as phosphorothioated C. The bases agree, but the printed modifications differ. "
    "This model uses the later iMARGI-specific C* oligo; the 2019 PNAS paper's "
    "statement that it reused the MARGI linker sequence does not establish "
    "which modification was in its experimental oligo."
    " Apparent reagent typo in Wu et al. 2019 detailed protocol: stock-preparation text "
    "(extracted JATS lines 395–397) says ATP for a reagent labeled 10 mM dATP, "
    "while its dA-tailing reaction table (lines 673–685) and Yan et al. 2019 "
    "Methods (line 203) specify dATP. This model uses a 3′-dA tail from the "
    "reaction descriptions; the reagent actually used in the experiment remains unverified."
)
READ_LENGTHS={"Read 1":100,"Read 2":100}
def oligos():
 return [
  "<b>bivalent linker top</b> — 5′-/Phos/NNAGATCGGAAGAGCGTCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCTTCC*GA/iBiodT/CT-3′ (C* is phosphorothioated)",
  "<b>bivalent linker bottom</b> — 5′-/Phos/GATCGGAAGAGCACACGTCTGAACGGATCCTCCCTACACGACGCTCT-3′",
  "<b>BamHI cut oligo</b> — 5′-TCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCT/3InvdT/-3′",
  "<b>NEBNext universal PCR</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′",
  "<b>NEBNext index PCR</b> — 5′-CAAGCAGAAGACGGCATACGAGAT[6 bp barcode]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-S-T-3′",
 ]
BRIDGE=BivalentRNAChromatinBridge("iMARGI bivalent linker")
def rna_ligation_geometry():
    rna=Scene(); rna.strand("RNA arm",[Segment("preadenylated RNA-ligation end","[RNA-arm end]",placeholder=True),Segment("annealed bivalent-linker body","[bivalent linker]",placeholder=True)],label="iMARGI linker")
    rna.ligation_end("RNA arm","preadenylated RNA-ligation end",kind="single-strand",polarity="5-prime",owner="adapter",text="preadenylated 5′ RNA-ligation end")
    return AdapterLigationGeometry.known("single-strand",rna,statement="The detailed protocol identifies a preadenylated RNA-compatible end on the annealed bivalent linker; this is a single-strand end ligation, not a duplex overhang.")
def dna_ligation_geometry():
    linker=Scene()
    linker.strand("linker top",[
        Segment("5-prime RNA ligation tail", "NNAGATCGGA"),
        Segment("paired linker core", "AGAGCGTCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCTTCCGATC"),
        Segment("3-prime T overhang", "T")],label="iMARGI bivalent linker")
    linker.anneal("linker bottom",[
        Segment("paired bottom core", "GATCGGAAGAGCACACGTCTGAACGGATCCTCCCTACACGACGCTCT")],
        to="linker top",pair=("paired bottom core","paired linker core"),
        label="5′-phosphorylated bottom strand")
    linker.ligation_end("linker top","3-prime T overhang",kind="cohesive",
                      polarity="3-prime",owner="adapter",text="3′-T linker overhang")
    substrate=dA_tailed_scene(repair_and_dA_tail(Construct([
        Segment("AluI-cut genomic fragment", "X"*24,placeholder=True)])),
        label="dA-tailed genomic DNA end")
    return AdapterLigationGeometry.known("cohesive",substrate,linker,
        statement="The detailed iMARGI protocol describes dA tailing and sticky-end DNA ligation (Wu et al. 2019, lines 101–103); the original MARGI Fig. 1 identifies a linker 3′-T overhang paired with the genomic 3′-A tail. The 2019 linker strands print the same paired core and terminal T (lines 273–285).")
def sections():
    wf=Workflow(MolecularState("Crosslinked and fragmented chromatin",tuple(BRIDGE.unjoined().rows())))
    wf.ligate_adapter("Ligate the RNA end of the bivalent linker in situ",rna_ligation_geometry(),BRIDGE.rna_joined().rows(),note="The linker end joins fragmented RNA; its DNA-reactive end remains free.")
    wf.ligate_adapter("Ligate the bridge DNA end to nearby genomic DNA",dna_ligation_geometry(),BRIDGE.chimera().rows(),note="The second ligation makes the diagnostic RNA–linker–DNA chimera in the nucleus.")
    wf.react("Reverse-transcribe, purify and make the paired-end library",[Row(chunks=[("RNA-derived cDNA—linker—genomic DNA → paired-end adapters","r1",False)])],note="The linker orientation distinguishes RNA and DNA read ends.")
    return [(r.action,r.after.rows,r.note) for r in wf.reactions]

def protocol_graph():
 source="Defining source: https://doi.org/10.1073/pnas.1819788116; detailed protocol: https://doi.org/10.1038/s41596-019-0229-4"
 graph=ProtocolGraph("iMARGI",authority=source)
 chromatin=graph.input("chromatin","cross-linked RNA–chromatin contacts",
                       MaterialKind.BIOLOGICAL,primary=True)
 nuclease=graph.input("restriction_enzyme","chromatin restriction enzyme",MaterialKind.CATALYST)
 fragmented=graph.product("fragmented_contacts","fragmented RNA and restriction-cut chromatin",
                          MaterialKind.MOLECULAR,state=BRIDGE.unjoined())
 graph.react("fragmentation","Fragment cross-linked chromatin",StepKind.DIGESTION,
             ((chromatin,InputRole.TEMPLATE),(nuclease,InputRole.CATALYST)),(fragmented,))
 linker=graph.input("bivalent_linker","annealed iMARGI bivalent linker",MaterialKind.ADAPTER)
 rna_ligase=graph.input("rna_ligase","RNA-end ligase",MaterialKind.CATALYST)
 rna_joined=graph.product("rna_joined","RNA-end-ligated bivalent linker contacts",
                         MaterialKind.MOLECULAR,state=BRIDGE.rna_joined())
 graph.react("rna_ligation","Ligate the RNA end of the bivalent linker in situ",
             StepKind.LIGATION,
             ((fragmented,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
              (rna_ligase,InputRole.CATALYST)),(rna_joined,),geometry=rna_ligation_geometry())
 dna_ligase=graph.input("dna_ligase","DNA-end ligase",MaterialKind.CATALYST)
 chimera=graph.product("contact_chimera","RNA–linker–DNA contact chimera",
                       MaterialKind.MOLECULAR,
                       state=BRIDGE.chimera())
 graph.react("dna_ligation","Ligate the bridge DNA end to nearby genomic DNA",
             StepKind.LIGATION,
             ((rna_joined,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
              (dna_ligase,InputRole.CATALYST)),(chimera,),geometry=dna_ligation_geometry())
 rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
 cdna=graph.product("cdna_chimera","cDNA–linker–genomic DNA chimera",MaterialKind.MOLECULAR)
 graph.react("reverse_transcription","Reverse-transcribe the RNA side of the chimera",
             StepKind.REVERSE_TRANSCRIPTION,
             ((chimera,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
              (rt,InputRole.CATALYST)),(cdna,))
 p5=graph.input("p5_primer","P5 library PCR primer",MaterialKind.OLIGO)
 p7=graph.input("p7_primer","indexed P7 library PCR primer",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","Purify and PCR-complete the paired-end library",StepKind.PCR,
             ((cdna,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The compact Illumina branch records the completed library rather than separate PCR landing objects."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
     read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
 graph.declare_reviewed_scope((
     "chromatin","restriction_enzyme","fragmentation","bivalent_linker","rna_ligase",
     "rna_ligation","dna_ligase","dna_ligation","reverse_transcriptase",
     "reverse_transcription","p5_primer","p7_primer","polymerase","library_pcr","library"),source)
 return graph.require_publishable()
