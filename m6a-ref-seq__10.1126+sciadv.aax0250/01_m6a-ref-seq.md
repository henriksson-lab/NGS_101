# m6A-REF-seq

Zhang et al., *Single-base mapping of m6A by an antibody-independent method*, Science Advances (2019), [defining article](https://doi.org/10.1126/sciadv.aax0250).

## Sources inspected

🟢 Published original PMC6609220: the manifest identifies open-access publisher content, not an author manuscript. Read Results, recognition and control design (JATS49–69); Methods “RNA materials”, MazF validation, FTO demethylation and “m6A-REF-seq library construction” (139–169); read analysis (171–175); and the separate ligation/qPCR validation context (177–183).

🟢 Published supplementary-material PDF `aax0250_SM.pdf`: printed title/DOI/date and contents; Figure S1–S3/S15 captions; Table S1 data headings and treatment groups; Table S2/S3 headings identifying site-validation probes and qPCR primers. These are not the library kit's adapter table. The scraper ran on the cached article and supplement text. Relevant validation-oligo hits were read in their source context. No rendered image-only sequence table is claimed inspected.

These two originals are the bounded pass. The cached article PDF/XML represent the same original, and the supplemental index is metadata, not an additional independent chemistry original. No superseded-preprint content is used. See [retrieval manifest](ref/MANIFEST.md).

## Boundary relative to MAZTER-seq

🟢 Both independently published methods use MazF sensitivity to the first A in ACA. The represented m6A-REF route explicitly prepares matched FTO-demethylated and minus-FTO-buffer samples before MazF, then uses Vazyme T4PNK and the VAHTS Small RNA kit NR801. MAZTER's inspected no-IP route uses FastAP/PNK, source 3iLL and 5iLL-22 ligations and SuperScript III. Those are different published preparation declarations, not an inferred equality or recipe inheritance. No directed descent between the two methods is asserted from motif or name similarity.

The ChpBK UAC experiments and MazF-K56A experiments are separate validation/screening branches in this article. They are not treated as additional genome-wide sequencing protocols.

## Matched sample prefixes

1. 🟢 Start with source poly(A) RNA: HEK293T total RNA was extracted with TRIzol then purified with Dynabeads, whereas tissue mRNA was purchased. The graph begins at the purified mRNA boundary and keeps these upstream preparations distinct.
2. 🟢 First heat and chill mRNA, before preparing matched FTO-positive and minus-enzyme samples. The FTO reaction uses Fe(II), alpha-ketoglutarate, ascorbate, RNase inhibitor and Tris; incubate at room temperature, then stop with EDTA. The source says the MazF control receives the corresponding reaction mixture without FTO protein. The model records treatment history without declaring complete m6A erasure or assigning residual/intermediate fractions.
3. 🟢 Heat both samples again, chill and digest with MazF. The represented event cuts immediately before ACA; methylated and other uncleaved RNA remains in the reaction population. This is not a complete digest simulator. Source secondary-structure and control filtering further limit interpretation.
4. 🟢 Purify with RNA Clean & Concentrator-5. Treat with source Vazyme T4PNK, then purify again. Exact newly created RNA end chemistry and its mechanistic repair remain unassembled. The graph records these operations without assuming an ATP formulation from another supplier's PNK buffer.
5. 🟢 The source continues with historical VAHTS NR801 small-RNA library preparation, {{= ref_seq.PCR_CYCLES }} PCR cycles, urea-PAGE library purification and HiSeq X10 sequencing. The current molecular prefix ends before that unreviewed kit; an explicit source item keeps its missing model coverage visible. No adapter, reverse transcriptase or PCR primer is supplied from MAZTER or a newer kit.

## Oligos and readout

The main Methods' unmodified RNA-ACA validation substrate is `{{= ref_seq.CONTROL_RNA }}` (RNA, written 5′→3′). The corresponding RNA-m6ACA control replaces the first A of that ACA with m6A. These are validation substrates, not assumed spike-ins or adapters in the NGS preparation. The unmodified control has a source-transcription check. Exact NR801 library oligos/modifications and run-primer landing sites remain unassigned.

🟢 Published analysis explicitly processes paired-end reads and counts internal versus terminal ACA/TGT motifs. Candidate coverage, read-length, RNAfold-structure and FTO-associated change filters are computational. They do not physically remove RNA or change a methylated base into a sequencing mutation. Complete run cycles, index roles/orientation and physical read geometry remain unassigned.

## Deferred questions

- Identify and read the applicable historical NR801 original and construct its physical end, ligation, RT, PCR and run states. An unassigned historical revision is an open question, not an invented inaccessible document.
- Read the directly cited Imanishi et al.2017 MazF methyl-sensitivity original, source reagent documentation and remaining original panels/tables; reconstruct exact cleavage-end/PNK chemistry only where established.
- Further inspect the distinct ChpBK, FTO-validation and ligation/qPCR controls, plus source-linked supporting methods, without merging them into this NGS route.
- An independent source/linked-information pass remains pending. No missing-information claim is confirmed.

## Bounded publication review, 30 October 2026

The published non-manuscript original PMC6609220 Methods139–175 was re-read, especially library-construction167–169 and its separate analysis context. Mandatory bundle oligo scraping was repeated; the printed validation substrates and actual reaction/library paragraph were checked directly. This pass adds one defining original, not a full supplemental or linked-source review. The source-reported 10 ng NR801 input, 15-cycle PCR, 8% urea PAGE and HiSeqX10 remain positive facts; they do not establish the historical kit's complete adapter/RT/PCR/read polymers.

The model keeps matched FTO/MazF and minus-FTO control branches, then stops at purified repaired RNA. Required `historical_nr801` remains unaccounted for. Actual source Qubit analysis, 10 ng kit-input allocation/remaining stock, source-positive DNA-library production, physical gel purification and loading/remainder interfaces remain to be assembled locally with the historical kit/end/read geometry. The current source names NR801 and Vazyme T4PNK; it does not license replacing those declarations with MAZTER-seq, contemporary NR801 or another supplier's PNK recipe.

The applicable historical NR801 documentation identity/revision remains an open target, not an invented unavailable document. The source main and published supplement are available; kit-specific full assembly and linked-source disposition remain unresolved. Validation oligos and the separate dot-blot/ChpBK/K56A/ligation-qPCR experiments do not become sequencing adapters or automatic controls in this selected NGS graph. No missing sequence, reaction detail or read layout is declared confirmed absent.
