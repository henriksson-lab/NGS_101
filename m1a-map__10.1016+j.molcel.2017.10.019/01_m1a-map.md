# m1A-MAP

Li et al., *Base-Resolution Mapping Reveals Distinct m1A Methylome in Nuclear- and Mitochondrial-Encoded Transcripts*, Molecular Cell (2017), [defining article](https://doi.org/10.1016/j.molcel.2017.10.019).

## Sources inspected

🟢 Associated accepted manuscript PMC5722686, canonical PMC Cloud XML with assigned URL/MD5. It is not established as procedurally identical to the final version of record. Read RNA isolation (Methods144–146), the separate general demethylation assay (172), the actual m1A-MAP section (176–180), distinct targeted validations (184–188), processing/deduplication (264–272), Methods S1 identity (324–334), and Figure1 caption (600–608). The mandatory primer scraper ran; main-library oligos and modifications were checked in Methods180 rather than borrowed from the shRNA or targeted-PCR hits.

The explicitly linked Methods S1 original, `NIHMS915125-supplement-2.pdf`, reached a PMC download gate. Other article supplements are enumerated in the external manifest but not read. Superseded preprint material is excluded. This capped pass reads one chemistry original; its inaccessible companion does not establish absent information. See [retrieval manifest](ref/MANIFEST.md).

## Source preparation and library branches

🟢 The represented route starts at poly(A)+ RNA. The main RNA-isolation section describes TRIzol, DNase treatment, small-RNA depletion and two oligo(dT) selections upstream; that full preparative recipe is outside this graph.

1. 🟢 Fragment the source poly(A)+ RNA using magnesium fragmentation. Keep the fragmented-input preparation separate from anti-m1A-enriched RNA. The cited prior IP original is unreviewed; neither perfect enrichment nor an exact physical partition is predicted.
2. 🟢 Prepare IP-minus-demethylase and IP-plus-AlkB samples. The source expresses an N-terminally truncated ΔN11 AlkB, not an implicitly interchangeable wild-type preparation. Its main-MAP mixture/incubation, EDTA quench, phenol/chloroform extraction and ethanol recovery are represented. Residual m1A, intermediates and nonspecific material remain unquantified. The separate general demethylation assay uses different input/mixture conditions; its heating/program is not silently imposed on the main-MAP route. The separately described Dimroth procedure is outside this graph and is not assigned an A product here. The source does not establish a matched-buffer treatment for every control.
3. 🟢 All three source RNA preparations enter the same locally declared library operations: PNK3′end treatment and heat inactivation; ethanol recovery; 3′linker ligation with truncated KQ RNA ligase2; sequential5′deadenylase and RecJf treatment/inactivation; ethanol recovery; heat RNA with the printed RT primer; TGIRT copying; ExonucleaseI treatment and silane recovery; source cDNA-adapter ligation with RNA ligase1; silane recovery; PCR; and8%TBE-gel recovery.

Digested primers/adapters and other molecular products remain in their reaction mixtures until the specified physical recovery. Complete cleanup purity is not inferred. The cDNA recovery does not acquire an invented RNA-hydrolysis step. Historical buffer/end chemistry, full strand assembly and detailed Methods S1 remain explicit model boundaries.

## Printed library oligos

The following are source strings with the main modifications retained. Canonical sequence parts are imported from `illumina`; local source checks establish their agreement with this article, without importing an eCLIP recipe.

| Reagent | Source ordering string, 5′→3′ |
|---|---|
| RNA3′linker | rAPP-{{= m1amap.RIL19 }}-3SpC3 |
| RT primer | {{= m1amap.RT_PRIMER }} |
| cDNA adaptor | 5Phos-{{= m1amap.cdna_adapter('example').top() }}-3SpC3 |
| Forward PCR | {{= m1amap.PCR_F }} |
| Reverse PCR | {{= illumina.P7 }}XXXXXX{{= m1amap.PCR_R_TAIL }} |

🟢 The source identifies XXXXXX as the sample index. Its ten random cDNA-adapter bases are used for duplicate handling in Methods272, supporting a structured UMI role. They are not native RNA sequence. C3 is a spacer, not an additional cytidine. The model represents the unassigned sample-index pool as N6 while preserving the source's X notation above.

🟡 Scene pairing derives the printed RT binding subrelation, retaining its nonpairing5′A. That does not establish the full RNA/linker substrate geometry. The printed PCR reverse suffix is shorter than the canonical sequencing-primer string; the model uses the printed PCR sequence, without treating a primer-length difference as a publication error or silently extending it.

## Chemical product versus readout

🟢 The source detects a change in mismatch rates between demethylase preparations. RT does not erase an RNA modification simply because its DNA copy mismatches. The local selected-outcome class therefore retains the RNA class separately from a reference-oriented reported base and its complementary physical cDNA base. Arrest yields no incorporated nucleotide at that site. Caller-selected events are not mutation-rate or complete-copy predictions; unmodified controls are not forced to copy perfectly.

🟢 The associated manuscript reports HiSeqX10 paired-end {{= m1amap.READ_LENGTH }}. Methods268–272 assign the first10bases of Read2 to the random barcode and remove that prefix computationally. Exact physical run-primer landing and index acquisition remain unassembled. Computational trimming, mismatch calling and duplicate collapsing do not discard physical RNA/DNA during preparation. The model makes no genomic-site, modification-stoichiometry or detection-specificity conclusion.

## Lineage and remaining sources

The source explicitly adapts eCLIP library construction with local modifications (Methods180). This permits a scoped methodological link, not inheritance of eCLIP's full recipe, RT enzyme or product geometry. m6A-SAC later cites this method for a library subroute; the detection chemistry is different.

Pending: final-VOR applicability; Methods S1 and remaining actual supplements; the directly cited prior IP/eCLIP originals and the remainder of this article's AlkB preparation; precise RNA/linker/cDNA/PCR/run geometry; and an independent original/linked-material pass. No missing-information claim is accepted.
