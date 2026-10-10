"""TrAEL-seq native 3-prime-end capture."""
from multibranch import illumina_branch
from protocol_extensions import native_three_prime_end
from chemdraw import Row
TITLE="TrAEL-seq — native DNA 3-prime-end sequencing"
NOTES="01_trael-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1371/journal.pbio.3000886">Kara et al., <i>PLOS Biology</i> (2021)</a>.'
SUMMARY="A pre-adenylated, biotinylated hairpin adapter ligates selectively to native DNA 3-prime OH ends; fill-in, fragmentation and affinity selection retain the captured end for directional Illumina sequencing."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("TrAEL-seq","native-end-adjacent genomic DNA")
FINAL_CAPTION="Read orientation starts at genomic sequence adjacent to the captured native 3-prime end."
SEQUENCING_INTRO="Declared Illumina primers bind the completed TrAEL-seq library."
def oligos():
    return ["<b>TrAEL adapter 1</b> — 5′-[Phos]NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTU…[BtndT]…-3′","<b>Adapter 2</b> — 5′-[Phos]GATCGGAAGAGCACACGTCTGAACTCCAGTCUUUUGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3′"]
def sections():
    return [("Ligate pre-adenylated adapter 1 to native 3-prime OH ends",native_three_prime_end().rows(),"Truncated RNA ligase joins the adapter without ATP; its random bases identify molecules and internal biotin enables capture."),("Fill the adapter hairpin and displace its short arm",[Row(chunks=[("captured 3′ end—adapter hairpin → Bst fill/displacement → tagged duplex", "r1", False)])],"Polymerase converts the captured junction into a stable duplex."),("Fragment, end-repair and ligate adapter 2",[Row(chunks=[("tagged duplex — sonication — end repair ** adapter 2", "me", False)])],"The second adapter supplies the opposite amplification end."),("Capture biotin and open USER sites",[Row(chunks=[("streptavidin ● biotin-tagged fragment → USER release → indexed PCR", "w1", False)])],"Affinity selection retains molecules descending from a native 3-prime end."),]
