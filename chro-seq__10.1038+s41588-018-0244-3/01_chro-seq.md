# ChRO-seq — research notes

Chu *et al.*, *Nature Genetics* (2018),
[10.1038/s41588-018-0244-3](https://doi.org/10.1038/s41588-018-0244-3),
[original author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6204104/).

## Sources and scope

- 🟢 Original article Methods: chromatin isolation, ChRO library preparation and mapping
  (`s41588-018-0244-3_PMC6204104.xml.txt`, lines 150–158), plus the actual external-link
  inventory and references 27–28.
- 🟢 Original Supplementary Table 1, first sheet of the MOESM3 workbook (text lines 1–40):
  sample/replicate-specific run-on method and “Unique molecular index barcode” column.
- 🟢 Mahat *et al.*, *Nature Protocols* (2016),
  [10.1038/nprot.2016.086](https://doi.org/10.1038/nprot.2016.086): original PRO-seq
  overview, Table 1 and the adapter/end-repair steps. This is the explicitly cited
  predecessor, not proof that every ChRO library used unchanged oligos.
- 🟢 Supplementary Notes 1–7 and the table legends were read for referring scope: assay
  comparisons, biological analysis and the dREG-HD computational method. This does not
  constitute a complete visual audit of the 30 supplementary figures or all result tables.

The oligo scraper was run on the complete bundles before source reading. Its output
is a reading aid, not proof that source information is absent.

## Source-owned ChRO route

1. 🟢 NUN isolation retains a chromatin pellet, followed by washes and sonication in
   chromatin storage buffer to suspend the pellet. The source retains engaged polymerases.
2. 🟢 Run-on supplies ATP and GTP with biotin-11-CTP and biotin-11-UTP in Sarkosyl-containing
   buffer. The nascent RNA acquires affinity labels during extension; those labels are
   not UMIs and have no assigned eight-base length.
3. 🟢 Trizol LS stops the reaction; GlycoBlue-assisted recovery yields RNA. ChRO then uses
   base hydrolysis, intended to yield 40–100-base RNA fragments.
4. 🟢 First streptavidin capture and Trizol release recover biotin-containing RNA.
5. 🟢 The 3′ RNA adapter is ligated using the source M0204L reagent.
6. 🟢 Second streptavidin capture is followed by RppH decapping, PNK 5′ phosphorylation
   and Trizol purification. The cited predecessor makes the intervening bead release
   explicit; these steps prepare RNA for 5′-adapter ligation.
7. 🟢 Ligate the 5′ RNA adapter, then perform the third streptavidin capture/recovery.
8. 🟢 Reverse transcription (source 18080-044) and PCR (M0491L) generate an Illumina-compatible
   library. The defining article names NextSeq 500, but no universal cycle mask is selected
   by this model.

The local graph declares all three captures, both adapter ligations and their source
order. Extraction and binding/wash removed fractions have physical discard fates; these
are not computational deduplication or a guarantee that every unbound molecule lacks label.
The graph does not infer an exact adapter-ligation or PCR geometry from a finished endpoint.

## ChRO and leChRO differ

🟢 **leChRO is length-extension ChRO.** RNase cocktail replaces the RNase inhibitor during
initial isolation; post-run-on hydrolysis is omitted. It is not an extra linear-amplification
branch. The previous note's description was a repository error and has been corrected.
This page's reaction graph shows the ChRO branch and annotates the leChRO differences.

🟢 Source mapping interprets the nascent RNA **3′ end from the read 5′ end** in ChRO, and
RNA **5′ end from the read 3′ end** in leChRO. These are read-orientation statements, not
an assignment of an arbitrary modern sequencing-primer family.

## Oligos and identifiers

🟢 Mahat 2016 Table 1 prints these **predecessor** oligos:

- VRA3: `5′-pGAUCGUCGGACUGUAGAACUCUGAAC-/Inverted dT/-3′`.
- VRA5: `5′-CCUUGGCACCCGAGAAUUCCA-3′`, unphosphorylated at its 5′ end.
- RP1, for RT and PCR: `5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA-3′`.
- RPI-n: `5′-CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′`.
  Its six Ns are a **sample index**, not the optional ChRO molecule identifier.

🟢 The ChRO article says some libraries have a molecule identifier in the first six
sequenced bases. Original Table 1 actually marks samples, including `UMU94-13_leChr`,
`GBM-88-04_Chr` and `PDX-88-04_P57`, among others. Marks vary between replicate/subsample
rows. Blank cells and combined rows do not define an immutable ordered adapter identity.
The previous statement that sample membership was undocumented was incorrect.

The public page lists predecessor sequences without assigning them unchanged to every
ChRO sample. The endpoint retains molecular roles while exact study-specific adapter arms,
UMI-bearing variants, primer sites and run masks remain open model questions. This does
not assert that those details were never published.

## Remaining source obligations

Reference 27 is Wuarin and Schibler's 1994 chromatin-isolation method; reference 28 is the
2016 PRO-seq protocol. The article also directly links GSE117832 and the Danko-Lab PRO-seq
processing pipeline. The review ledger records actual original-source retrieval and reading
scope. An unspecified manufacturer protocol is an unresolved source identity, not an
invented unavailable document. Full supplementary visual/result-table scope and relevant
linked material still require bounded dispositions and independent original readings.
No missing-information claim is accepted at this checkpoint.
