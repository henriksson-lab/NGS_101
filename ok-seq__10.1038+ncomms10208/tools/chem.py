"""Human OK-seq strand-preserving Okazaki-fragment library."""
from multibranch import illumina_branch
from protocol_extensions import native_three_prime_end
from chemdraw import Row
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
    return [("Pulse-label nascent DNA with EdU",[Row(chunks=[("lagging strand: short Okazaki fragment — EdU", "w1", False)])],"EdU marks newly synthesized fragments without requiring ligase-deficient cells."),("Denature and size-purify short single strands",[Row(chunks=[("genomic DNA → denature → sucrose gradient → <200-nt ssDNA", None, False)])],"The physical selection enriches Okazaki fragments."),("Click biotin onto EdU and phosphorylate 5-prime ends",native_three_prime_end(captured="Okazaki-fragment end").rows(),"Azide–alkyne chemistry creates the affinity handle; PNK makes ligatable 5-prime ends."),("Ligate random-overhang adapters and capture",[Row(chunks=[("A1 ** Okazaki ssDNA ** A2 → streptavidin capture → PCR", "me", False)])],"Adapter orientation and single-end sequencing preserve the genomic strand signal."),]
