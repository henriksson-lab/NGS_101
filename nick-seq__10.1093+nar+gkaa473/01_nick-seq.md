# Nick-seq: complementary nick-translation and terminal-tailing libraries

Cao et al. (2020), [Nick-seq for single-nucleotide resolution genomic maps of DNA modifications and damage](https://doi.org/10.1093/nar/gkaa473).

## Capped sources and scope

One reader,28October2026: published PMC7337925 main article, complete Materials and Methods and mechanism/discussion; published Supporting Information DOCX from the original supplementary ZIP (captions and source-file inventory). Mandatory complete-bundle scraper run; its high-scoring sequence hits are mapped genomic neighborhoods, not ordering oligos. The DOCX contains EMF images for Supplementary Figures1/2 and Table1; text extraction omits those graphics. Rendering was attempted with installed LibreOffice but failed in this sandbox, so those image details remain unread. No absence claim follows from the extraction. Publisher supplement originals are identified separately from the fetched superseded bioRxivv2; no preprint chemistry is adopted.

The graph represents the **published E.coli nicking-enzyme validation preparation**, followed by the positive NT and TdT libraries. The source's negative-control processing and lesion-specific applications are described below; their exact kit reactions are not silently borrowed or declared complete.

## Common preparation

1. Fragment1µg genomic DNA in **three separate restriction preparations**: NciI; HindIII+XhoI; SalI+XbaI+NdeI. Include RNaseA. Qiagen PCR-purify each and pool. These are not six enzymes mixed in one reaction.
2. Shrimp alkaline phosphatase37°C30min removes3′ phosphates; heat70°C10min. PolI with ddNTPs37°C40min blocks pre-existing extendable ends. Add phosphatase37°C30min to process excess ddNTPs, heat75°C10min and DyeEx desalt. Nucleotide-processing products are not physically removed merely by phosphatase catalysis.
3. Create new nicks with Nb.BsmI and Nb.BsrDI65°C1h, heat80°C20min, cool slowly. The source proceeds without further purification. The new3′OH class is separate from the previously dd-blocked background. Split material between NT and TdT; do not label both halves as the same library chemistry.

## NT branch

PolI incorporates α-thio-dNTPs at15°C90min, generating nuclease-resistant phosphorothioate-containing copied DNA. The source negative control substitutes water for PolI. Heat75°C20min and DyeEx purify. ExoIII37°C1h,95°C3min/ice, RecJf37°C1h then80°C10min digest original/background DNA. Some complex samples receive additional RecJf. **Digestion products remain in the tube until Zymo Oligo Clean & Concentrator recovery.** Resistance is not a guarantee that every background molecule is destroyed.

The recovered material is converted to a library with the source-named **Clontech SMART ChIP-seq kit**, following its instructions, with12 final PCR cycles. Exact historical kit adapter/primer/template-switch polymers remain unassembled in this pass. AMPure size selection follows pooling with the corresponding source negative control. NextSeq500 paired-end75bp.

## TdT branch

Nick-converted DNA is denatured95°C3min, then **polyT-tailed**,37°C20min using TdT and the source DNA SMART T-Tailing Mix; terminate70°C10min. Primer annealing/template switching uses the source SMART ChIP-seq kit; Illumina primers supplied with the kit amplify for12cycles, followed by barcode-specific pooling/AMPure size selection. NextSeq500 paired-end75bp. Do not substitute polyA tails or invent a standard Nextera endpoint. T-tail length and exact kit oligos are unassigned here.

The combined analytical workflow uses complementary NT and TdT data, not a chemical ligation of their products. The source trims the first3bases of Read1 during processing; this is not evidence of a3-base UMI. Tail-based site calling is ambiguous near genomic thymines, and NT helps resolve that limitation. A readout intersection does not validate the chemical identity of every detected nick.

## Other source-qualified front ends

For natural phosphorothioate mapping, iodine treatment of blocked Salmonella DNA creates breaks, followed by phosphatase removal of3′ phosphate. These natural PT sites differ from PT backbones deliberately installed during NT. For oxidative AP-site mapping, EndoIV treats DNA from H2O2-exposed E.coli. Main Methods lists exposure concentrations0.1/0.5/1/2mM; Results describes a0.2mM application. This published concentration difference is retained rather than harmonized by inference; the represented nickase-validation graph depends on neither value. The article's suggested applicability to more than25 damage/modification types is not a claim that all were experimentally demonstrated.

## Deferred sources and ancestry

Historical SMART ChIP-seq kit originals and exact oligos; rendered DOCX graphics; source code version and complete index assignment; negative-control molecular library processing, alternative lesion-conversion branches, and independent reader. These are open source/model tasks, not confirmed omissions. The article cites Zhu et al.2001 on reverse-transcriptase template switching as conceptual SMART background; exact study kit implementation requires its own source. Comparisons with SSiNGLe, SSB-seq, Excision-seq and Click-code-seq do not establish recipe inheritance.
