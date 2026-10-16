# SNU-Seq

**Defining paper:** Gerlevik et al., *Nucleic Acids Research* (2026), doi:10.1093/nar/gkag652.

## Molecular boundary

🟢 Pulse-labelled 4sU RNA is thiol-biotinylated and affinity-purified. Intact enriched RNA is polyadenylated with bacterial poly(A) polymerase before a 3′-anchored library is made.

🟢 The protocol selects labelled molecules without converting 4sU. Its informative coordinate is the native end retained by the intact-RNA preparation.

## Evidence boundary

🟡 SNU-Seq and sf4sU-Seq share a defining paper but are separate protocols: bacterial polyadenylation of intact RNA and size fractionation of short RNA do not produce equivalent endpoint information.

## Oligo source audit

🟢 The defining article's Supplementary Methods (`PMC13309788_supp_gkag652_supplemental_file.pdf.txt`, lines 380–388) specify Lexogen QuantSeq 3′ mRNA library preparation with either the Ion Torrent or Illumina FWD kit, 13 PCR cycles, and single-end 100-nt sequencing on Ion Proton or Illumina NextSeq2000. The authors' [GEO record for a SNU-Seq library](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5240711) confirms the Ion Torrent branch in one deposited run. The supplement prints no kit oligo bases.

🟢 The official [Ion Torrent guide 012UG036V0130 (2018)](https://www.lexogen.com/wp-content/uploads/2018/06/012UG036V0130_QuantSeq-IonTorrent.pdf), overview and Appendices D–E, gives the kit's oligo(dT)-primed first strand with a P1 arm, random-primed second strand with an A arm, 9–11-base inline barcode choices, and the finished library shell. Its Appendix E prints a `TCAG–barcode–CGAT–insert` read start; an index read is not needed for this branch. The official [Illumina FWD guide 015UG009V0260 (2021)](https://www.lexogen.com/wp-content/uploads/2021/08/015UG009V0260_QuantSeq_Illumina_2021-08-16.pdf), overview and Appendix G, prints the P5–Read 1–insert–Read 2–i7–P7 shell and sequencing-primer bases. It uses oligo(dT) first-strand and random second-strand priming, with i7 introduced by PCR. The paper and GEO do not identify each sample's exact kit revision or chosen index/barcode, so the single schematic keeps its platform-specific ends symbolic; this is a limit of the displayed combined branch, not an absence of published vendor library architecture.

🔴 Neither original paper nor the checked guides identify every proprietary primer reagent as an ordering oligo for the experiments. The actual revision and barcode well for each library remain unresolved; source-review status is blocked rather than confirmed missing.

🟢 The same Supplementary Methods (lines 375–388) place bacterial poly(A) polymerase before Lexogen library preparation. They do not establish a separate RNA-adapter ligation in the SNU-Seq branch, so the graph passes the tailed RNA directly to kit-based reverse transcription rather than drawing an invented ligation.
