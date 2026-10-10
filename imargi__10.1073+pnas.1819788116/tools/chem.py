"""iMARGI in-situ RNA–DNA bridge ligation."""
from multibranch import illumina_branch
from protocol_extensions import proximity_product
from chemdraw import Row
TITLE="iMARGI — in situ RNA–genome interaction mapping"
NOTES="01_imargi.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1819788116">Yan et al., <i>PNAS</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-019-0229-4">Nat. Protocols (2019)</a>.'
SUMMARY="A bivalent linker is ligated in situ first to fragmented chromatin-associated RNA and then to nearby genomic DNA; reverse transcription and Illumina preparation preserve the RNA–linker–DNA order."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("iMARGI","RNA–linker–DNA contact chimera",umi=8)
FINAL_CAPTION="Paired reads recover the RNA-derived and DNA-derived ends separated by the diagnostic bridge."
SEQUENCING_INTRO="Declared Illumina primers bind the completed iMARGI library."
def sections():
    return [("Crosslink and fragment chromatin-associated RNA and DNA",[Row(chunks=[("RNA … crosslinked proximity … restriction-cut genomic DNA", None, False)])],"The physical interaction is fixed before either ligation."),("Ligate the RNA end of the bivalent linker in situ",[Row(chunks=[("chromatin RNA ** RNA arm—annealed bridge", "me", False)])],"A pre-adenylated RNA-compatible linker end joins fragmented RNA."),("Ligate the bridge DNA end to nearby genomic DNA",proximity_product(left="chromatin RNA",right="nearby genomic DNA",bridge="iMARGI bivalent linker").rows(),"The second ligation makes the diagnostic RNA–linker–DNA chimera in the nucleus."),("Reverse-transcribe, purify and make the paired-end library",[Row(chunks=[("RNA-derived cDNA—linker—genomic DNA → paired-end adapters", "r1", False)])],"The linker orientation distinguishes RNA and DNA read ends."),]
