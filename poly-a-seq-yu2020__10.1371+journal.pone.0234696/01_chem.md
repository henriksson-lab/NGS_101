# Poly(A)-seq — Yu2020 — capped source-qualified notes

Yu et al. (2020), Poly(A)-seq: A method for direct sequencing and analysis of the transcriptomic poly(A)-tails. DOI10.1371/journal.pone.0234696.

🟢 The published main Methods77–99/Results126–130,299 were read. S1/S2 workbook headers/initial rows were inspected to establish their actual scope: they contain own-study/comparison mapping statistics, respectively. These administrative identity checks count toward the three-original cap; the remaining tables/figures and historical vendor originals remain unread. Mandatory whole-bundle primer scraper completed. Spike-template hits are not automatically sequencing-library primers.

This entry denotes **Yu2020 Poly(A)-seq**, not the Derti2012 PolyA-seq namesake. It retains native poly(A) tails in RNA-derived fragments and combines direct3′adapter ligation with a different second-handle component. It is also distinct from PAT-seq's Klenow end-extension and biotin-DNA/RT-primer capture.

🟢 Initial NextSeq route: Trizol/RQ1DNase RNA →separate spectrophotometric purity/quantity and agarose integrity assays →partial RNaseT1 digestion37°C3min →5.1µg fragmentedRNA input →oligo-dT bead capture **after** fragmentation →GnomeGen sRNA-seq3′adapter ligation →adapter-primedRT →ScriptSeqv2 Terminal-Tagging DNA synthesis →cDNA purification →18cyclePCR →250–350bp or300–500bp product purification →quantitation/storage →NextSeq500 single-end300nt. Exact per-sample size-window assignments and historical kit sequences/conditions are unassembled.

🟢 Later HiSeqXTen/NovaSeq route is a separate source preparation:5µg HeLaRNA →RNaseT1 →oligo-dT capture →**then** independently prepared40A/120AspikeRNA addition →GnomeGenK02420-L3′adapter ligation →RT →SMARTer StrandedRNA-seq634837 second-strand component →source finishing interface →HiSeqXTenPE150 orNovaSeq6000PE250. The graph does not transfer ScriptSeqv2,18PCRcycles or initial size-window choices into that later route. The same method's later platform experiments are not extra protocol entries.

🟢 Main spike preparation separately describes DNA-template PCR →bead purification →BsrDI removal of downstream sequence →4%agarose/Qiagen recovery →T7 transcription →12%ureaPAGE. This upstream preparation remains unassembled, with the recovered RNA controls introduced explicitly as independent physical inputs after poly(A) capture. The source uses30pg/9pg additions in two repeats; the model leaves actual sample choice unresolved. Printed39/40nt upstream-template descriptions and actual transcript boundaries are not silently forced into a complete spike sequence.

🟢 RNaseT1 generates new phosphate-ended fragments while native3′OH ends remain candidates for ligation; no PNK/end-healing step is inserted before3′adapter ligation. Capture can bias toward longer tails, as the source itself discusses, so no perfect tail-independent recovery is claimed. Kit terminal-tagging/second-strand interfaces are kept opaque rather than treated as an established modern SMART/tagmentation operation.

🟡 A poly(T)-containing first DNA strand acts as sequencing template for poly(A) calls. This does not chemically convert the original RNA or DNA T into A. pA-finder, poly(G) trimming and sequence-quality filters occur after acquisition. Exact library oligos, primer landings, identifiers, kit revisions and full reactions remain unassembled; no index/UMI is inferred from generic kit names. Private incomplete graphs only; no public promotion or accepted absence.

Deferred source scope:

- HistoricalGnomeGen sRNA-seq/K02420-L, ScriptSeqv2 and SMARTerStranded634837 originals: identity/complete internals not inspected, no modern-version substitution.
- Source-cited terminal-tagging/strand-synthesis component originals and exact sample size-window assignments remain unreviewed/unassigned.
- PublishedS3/S4tables and sourceFigure1/S3plusremainingfigures remain unreviewed; spike-template endpoints and exactchemistry unresolved.
