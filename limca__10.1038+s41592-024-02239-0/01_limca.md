# LiMCA

**Defining source:** doi:10.1038/s41592-024-02239-0.

## Molecular boundary

🟢 Cytoplasm and nucleus are physically separated. Smart-seq2 profiles full-length mRNA; the nucleus undergoes restriction digestion, proximity ligation and META amplification.

## Evidence boundary

🟢 The contact ligation junction is explicit. The modalities are joined by cell provenance, not by a shared molecular barcode.

## Oligo source review

🟢 The defining article's open full text (PMC11166570, “LiMCA protocol”)
specifies Smart-seq2 `oligo-dT30VN`, a template-switch oligo and an
`ISPCR` primer for the cytoplasmic RNA branch. It describes NlaIII
digestion followed by T4 ligase proximity ligation of the nuclear branch,
then whole-genome amplification with **20 META tags**. Final cDNA uses a
Nextera Tn5 kit; the 3C and RNA libraries are sequenced separately at
2 × 150 bases. The article does not print the exact RNA oligo sequences
or the twenty META tag sequences in its methods.

🟡 The article points to a separate author protocol, Wu et al.,
doi:10.21203/rs.3.pex-2423/v1, for step-by-step LiMCA and its
“Supplementary Table 1” oligos. That Research Square document was not
retrievable through the source fetcher. The Nature Methods supplementary
files fetched for this paper are figures, peer-review correspondence and
source-data tables, not an ordering-level primer table. The current
role-only page therefore cannot distinguish verified LiMCA oligo bases
from inherited Smart-seq2/META chemistry; the author protocol and the
version-specific META tag table need review before exact primer-binding
geometry can be constructed.
