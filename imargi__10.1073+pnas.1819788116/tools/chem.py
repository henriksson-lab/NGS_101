"""iMARGI in-situ RNA–DNA bridge ligation."""
from multibranch import illumina_branch
from protocol_extensions import proximity_product
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Scene,Segment,Workflow
TITLE="iMARGI — in situ RNA–genome interaction mapping"
NOTES="01_imargi.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1819788116">Yan et al., <i>PNAS</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-019-0229-4">Nat. Protocols (2019)</a>.'
SUMMARY="A bivalent linker is ligated in situ first to fragmented chromatin-associated RNA and then to nearby genomic DNA; reverse transcription and Illumina preparation preserve the RNA–linker–DNA order."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("iMARGI","RNA–linker–DNA contact chimera",umi=8)
FINAL_CAPTION="Paired reads recover the RNA-derived and DNA-derived ends separated by the diagnostic bridge."
SEQUENCING_INTRO="Declared Illumina primers bind the completed iMARGI library."
def sections():
    rna=Scene(); rna.strand("RNA arm",[Segment("preadenylated RNA-ligation end","[RNA-arm end]",placeholder=True),Segment("annealed bivalent-linker body","[bivalent linker]",placeholder=True)],label="iMARGI linker")
    rna.ligation_end("RNA arm","preadenylated RNA-ligation end",kind="single-strand",polarity="5-prime",owner="adapter",text="preadenylated 5′ RNA-ligation end")
    rna_geometry=AdapterLigationGeometry.known("single-strand",rna,statement="The detailed protocol identifies a preadenylated RNA-compatible end on the annealed bivalent linker; this is a single-strand end ligation, not a duplex overhang.")
    dna_geometry=AdapterLigationGeometry.unavailable("The sources read establish the second linker-to-restriction-fragment ligation but the protocol note does not transcribe a source-supported adapter-overhang sequence or polarity.")
    wf=Workflow(MolecularState("Crosslinked and fragmented chromatin",(Row(chunks=[("RNA … crosslinked proximity … restriction-cut genomic DNA",None,False)]),)))
    wf.ligate_adapter("Ligate the RNA end of the bivalent linker in situ",rna_geometry,[Row(chunks=[("chromatin RNA ** RNA arm—annealed bridge","me",False)])],note="The linker end joins fragmented RNA.")
    wf.ligate_adapter("Ligate the bridge DNA end to nearby genomic DNA",dna_geometry,proximity_product(left="chromatin RNA",right="nearby genomic DNA",bridge="iMARGI bivalent linker").rows(),note="The second ligation makes the diagnostic RNA–linker–DNA chimera in the nucleus.")
    wf.react("Reverse-transcribe, purify and make the paired-end library",[Row(chunks=[("RNA-derived cDNA—linker—genomic DNA → paired-end adapters","r1",False)])],note="The linker orientation distinguishes RNA and DNA read ends.")
    return [(r.action,r.after.rows,r.note) for r in wf.reactions]
