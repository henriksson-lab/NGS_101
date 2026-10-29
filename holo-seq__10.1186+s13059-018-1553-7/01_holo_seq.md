# Holo-seq

Xiao et al.(2018), *Holo-Seq: single-cell sequencing of holo-transcriptome*. Defining original: https://doi.org/10.1186/s13059-018-1553-7 .

## Source scope

🟢 Published non-manuscript PMC6193298: Results75, carrier preparation165, library Methods177–197 and capture-library primers219–233. Published Additionalfile1/MOESM1 original DOCX: FigS1 original embedded image inspected, and FigS2/S8/S13 captions; Additionalfile4/MOESM4 original DOCX: step-by-step mRNA/total/small-RNA and dual branch recipes read. Mandatory full-bundle scraper run; its library PCR primers in main Methods are for the explicitly captured/re-amplified directional library branch, not automatically every Holo-seq branch. Supplementary expression workbooks and Additionalfile5 remain unreviewed. One capped source-reading pass only.

## Carrier and physical material

🟢 T7-transcribed carrier RNA supports low-input recovery. Source preparation digests template DNA with DNaseI before RNA purification. The source mixes poly(A)-negative/positive carrier species and supplies a separate25nt5′phosphorylated small-RNA carrier. FigS1 shows artificial DNA transcription templates with repeated NotI recognition motifs and a separate Cas9-targeted alternative. The read image establishes the architecture, not a full base transcription of every template. Poly(A) capture in the dual branch retains mRNA on beads and small RNAs in the supernatant; these are separate physical fractions from one cell.

NotI acts on carrier-derived library DNA, not by directly deleting carrier RNA. The source acknowledges that endogenous DNA can also contain NotI sites. The local selected-cleavage helper therefore preserves both DNA pieces and original outer-adapter associations for either declared origin; it predicts no site efficiency, no universal carrier specificity and no removal before physical purification. Exact staggered ends and full adapters remain unassembled. Source cleanup retains an unquantified mixture, not perfectly carrier-free native products.

## Branch order and boundaries

🟢 The mRNA E7490/E7530 route applies USER after ligation, purifies with1×AMPure (the detailed protocol adds100ng~2kb DNA carrier), then NotI37°C2h,1×AMPure, then16cycleQ5 enrichment and0.9×AMPure. The total-RNA E7420 route uses dUTP in second-strand synthesis, ligation,1×AMPure with the late DNA carrier, NotI,1×AMPure, then USER37°C15min in the pre-PCR mixture and16cycleQ5. These USER operations are source-reported; full historical adapter/second-strand substrate and USER cleavage geometry remain unreviewed rather than supplied from a newer kit.

🟢 Small-RNA E7300 places NotI after library PCR/1.8×AMPure, followed by1×AMPure with added100ng~2kb DNA carrier. Main Methods193 reports15PCRcycles, whereas original Additionalfile4 small-RNA step-by-step protocol prints22cycles. Both are retained; branch applicability/revision is unreconciled. The graph starts after source PCR and assigns no numerical count. A graph-owned source-value comparison preserves both locations. This is not silently repaired by taking a current small-RNA kit.

The partial graph constructors represent three source-owned post-adapter/PCR branches and the actual cleanup order. Late DNA carrier is a separate structural input at the specified cleanup and does not disappear by being called a carrier. Residual carrier-derived, native and late-carrier material remains unquantified through physical selection. Full upstream library chemistry, exact historical oligos, USER processing and acquired read layout remain explicit uncovered source items. No finished schematic or absent-information claim.
