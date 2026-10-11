"""ChAR-seq RNA–DNA proximity ligation."""
from multibranch import illumina_branch
from protocol_extensions import proximity_product
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Workflow
TITLE="ChAR-seq — chromatin-associated RNA sequencing"
NOTES="01_char-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.7554/eLife.27024">Bell et al., <i>eLife</i> (2018)</a>.'
SUMMARY="A biotinylated bridge joins chromatin-associated RNA to nearby genomic DNA in situ; cDNA conversion, fragmentation and junction enrichment yield paired RNA–DNA contact reads."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("ChAR-seq","RNA–bridge–DNA contact chimera",umi=8)
FINAL_CAPTION="The internal bridge orients the RNA-derived and genomic sides of each contact molecule."
SEQUENCING_INTRO="Declared Illumina primers bind the final ChAR-seq library."
def sections():
    return [("Crosslink and fragment chromatin-associated RNA",[Row(chunks=[("chromatin RNA — fragmented 3′ end … nearby DNA", None, False)])],"RNA remains spatially fixed relative to chromatin."),("Ligate the RNA side of the biotinylated bridge",[Row(chunks=[("chromatin RNA ** RNA bridge arm—biotin—DNA arm", "me", False)])],"The bridge records molecule orientation and provides affinity selection."),("Ligate the DNA side in situ",proximity_product(left="chromatin RNA",right="nearby genomic DNA",bridge="ChAR-seq bridge").rows(),"The resulting chimera contains one RNA-derived and one DNA-derived locus."),("Convert RNA to cDNA and enrich bridge-containing fragments",[Row(chunks=[("RNA cDNA—bridge(biotin)—DNA → streptavidin → paired-end library", "w1", False)])],"Junction enrichment removes molecules that lack the diagnostic bridge."),]

def workflow():
 s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])))
 rna=AdapterLigationGeometry.unavailable("ChAR-seq establishes RNA-side bridge ligation, but the cut-facing bridge strand and terminal modification have not yet been verified from the supplementary source.")
 dna=AdapterLigationGeometry.unavailable("ChAR-seq establishes DNA-side bridge ligation, but the strand-resolved DNA-facing bridge end has not yet been verified from the supplementary source.")
 w.ligate_adapter(s[1][0],rna,s[1][1],note=s[1][2]);w.ligate_adapter(s[2][0],dna,s[2][1],note=s[2][2]);w.react(s[3][0],s[3][1],note=s[3][2]);return w
