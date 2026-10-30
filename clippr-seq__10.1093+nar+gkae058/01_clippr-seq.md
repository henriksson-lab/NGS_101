# CLiPPR-seq

Singh et al., *Global identification of mRNA-interacting circular RNAs by CLiPPR-Seq*, Nucleic Acids Research (2024), [10.1093/nar/gkae058](https://doi.org/10.1093/nar/gkae058).

## Selected evidence and limits

🟢 The PMC-hosted PDF is visibly the printed OUP/NAR 2024 published advance article, pages 1–14, with its original OUP download footer. This pass uses its Methods PDF-text lines 77–122, Results methodology page 4, and original Figure 1 page 5. The associated supplementary PDF page 1 is a table/figure role map. Original Table S8 workbook headers and whole entries were inspected in ZIP/XML: PCR/qPCR primers, LNA knockdown oligos, biotin-ASO validation reagents and reporter-construction oligos have separate roles. They are not assigned as Collibri sequencing adapters. The mandatory bundle scraper reported 675 sequence-like hits in four representations and skipped three large result tables by its size limit; those result tables were not read or used as chemistry evidence.

The original investigation is closed. A version bookkeeping mistake exposed the associated XML's first identity header after three selected originals: it identifies EMS194846/UKPMC author-manuscript metadata, whereas the PDF is visibly the published layout. This fourth conservative administrative identity exposure is recorded, not hidden as cap compliance. No XML chemistry was adopted, no full published/manuscript equivalence is asserted, and no further source was opened. The whole atom maps, kit adapter bodies and finishing geometry remain unassembled in this partial declaration. This is not an absence audit.

## What is selected

🟢 AMT crosslinking secures RNA-duplex associations before poly(A) capture. Photoreversal releases associated RNAs, and RNase R enriches the captured circular-RNA population. Untreated CLiPP is the paired control branch, not another method. βTC6 is the selected biological branch; HeLa and differentiated C2C12 are applications of the same preparation.

The abstract and Methods heading call this “proximity ligation.” Figure 1 and the executed overview describe AMT/UV crosslinking and photoreversal; this model does not turn that wording into a phosphodiester RNA–RNA junction or a sequenced chimeric partner join. Predicted partners are assigned by BLAST comparison with the control mRNA library and validation, not by a molecular partner barcode or joined paired read.

## Physical trace

🟢 AMT stock 1 mg/ml is diluted in DMEM/2× PBS to irradiation medium 20 µg/ml; cells incubate 37°C for 30 minutes, then receive two 7-minute 365-nm exposures on ice, 5 cm from the source. The selected conventional TRIzol isolation route includes DNase treatment and actual RNA quality/quantity allocations; magnetic isolation is an alternative, not silently merged.

🟢 Select 10 µg total RNA for NEBNext poly(A) magnetic isolation. The source Results identifies oligo-dT bead capture. Capture support and its attached oligo inventory are represented without assuming a disclosed oligo length or kit recipe. Nonselected RNA and residual support are physical fractions. RNA QC consumes actual allocations with amounts unassigned.

🟢 Irradiate captured RNA for 5 minutes for photoreversal. Methods prints “254 nM”; Figure 1 and Results print 254 nm UV. The model uses Figure 1's wavelength interpretation and preserves the literal Methods notation without calling it an author's error. The exact damage/crosslink reversal efficiency is unassigned.

🟢 Divide the de-crosslinked preparation into actual halves: RNase R treatment for CLiPPR and untreated CLiPP. RNA degradation products remain copresent until later workup; no instantaneous physical separation or perfectly pure circRNA is assumed. The RNase R recipe is delegated to reference 20 and remains unreviewed.

🟢 Select 10 ng treated RNA and 25 ng untreated RNA for separate Collibri Stranded RNA libraries. The kit's actual supply inventory is represented with source-positive oligo, reverse-transcription/nucleotide classes; exact stock packaging, bodies and finishing reactions are unassembled. DNA ScreenTape/TapeStation QC and separate loading allocations precede NextSeq 550 single-end 150-base sequencing. No physical pooling is asserted.

## Open scope

Kit originals, reference 20 treatment details, original read-primer landing and full library molecules are unreviewed. Table S8 validation reagents are not merged into production. Crosslinking and RNase R enrichment do not prove exclusive direct interactions or complete removal of linear RNA. No barcode/UMI meaning is invented.
