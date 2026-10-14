# SNU-Seq

**Defining paper:** Gerlevik et al., *Nucleic Acids Research* (2026), doi:10.1093/nar/gkag652.

## Molecular boundary

🟢 Pulse-labelled 4sU RNA is thiol-biotinylated and affinity-purified. Intact enriched RNA is polyadenylated with bacterial poly(A) polymerase before a 3′-anchored library is made.

🟢 The protocol selects labelled molecules without converting 4sU. Its informative coordinate is the native end retained by the intact-RNA preparation.

## Evidence boundary

🟡 SNU-Seq and sf4sU-Seq share a defining paper but are separate protocols: bacterial polyadenylation of intact RNA and size fractionation of short RNA do not produce equivalent endpoint information.

## Oligo source audit

🟢 The defining article's Supplementary Methods (`PMC13309788_supp_gkag652_supplemental_file.pdf.txt`, lines 380–388) specify Lexogen QuantSeq 3′ mRNA library preparation with either the Ion Torrent or Illumina FWD kit, 13 PCR cycles, and single-end 100-nt sequencing on Ion Proton or Illumina NextSeq2000. The authors' [GEO record for a SNU-Seq library](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5240711) confirms the Ion Torrent branch in one deposited run. The supplement prints no oligo bases. The oligo(dT) capture primer, opposite library end, amplification/index primers and platform-specific sequencing entry remain unresolved from the paper and these run metadata; the schematic now keeps both platform ends symbolic and does not assert an 8-base index.
