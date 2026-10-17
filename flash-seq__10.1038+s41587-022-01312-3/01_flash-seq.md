# FLASH-seq: source-owned standard branch

Hahaut et al., “Fast and highly sensitive full-length single-cell RNA sequencing using
FLASH-seq,” *Nature Biotechnology* (2022),
[doi:10.1038/s41587-022-01312-3](https://doi.org/10.1038/s41587-022-01312-3).

Evidence: 🟢 source-positive · 🟡 model-derived/inferred · 🔴 not yet represented or checked.
These model limits do not assert that information is absent from the publication.

## Originals inspected

- Defining JATS Methods, lines 133–187: FS, FS-LA and FS-UMI reagent recipes;
  commercial XT, homemade Tn5 and plexWell choices; pooling and acquired read modes.
- Original Supplementary Fig. 1, p3, visually inspected by two readers: combined
  RT-PCR, amplifying TSO and FS-dT, crossed-out ISPCR, standard cleanup/QC versus
  the FS-LA shortcut, indexing, pooling and final cleanup.
- Supplementary Note 2, p79, both extracted text and original table image:
  the FLASH-seq table contains dT30VN and FS TSO. ISPCR belongs to the separate
  Smart-seq2 comparator table.
- Recovered article-linked author standard protocol **v3, 23 March 2022**, printed
  integer ID59800: oligos pp4–6, RT-PCR pp14–15, commercial library preparation
  pp21–22, sequencing/data processing pp26–29. Both readers checked its originals.
  The article's materials link (`protocols.io/view/flash-seq-protocol-b6myrc7w`, JATS
  line 167) is the URI that the integer-59800 public metadata assigns to the
  `kxygxzkrwv8j/v3` DOI, so this v3 is the article-cited index-adapter source.
  The later **v4, 12 January 2023**, remains a separate comparator. Matching selected
  components does not establish complete equivalence between revisions.
- The separate FS-UMI protocol and article identify its distinct oligos, primers and
  read modes. Previous supplement/branch checks remain in the source ledger; they
  do not establish a fully reconstructed FS-UMI graph.

All originals are obtained through `tools/get_sources.py`; the bundle and mandatory
oligo-scraper output are reading aids, not public page content.

## Standard branch and repository correction

🟢 FS uses biotinylated H25-dT30VN and H25-rGrGrG, with SuperScript IV and KAPA
HiFi HotStart in one combined RT-PCR mix. It does **not** add the separate ISPCR
reagent. The old repository page inherited Smart-seq2's ISPCR amplification and
LNA-ending TSO through the SMART-family page builder. That was a modelling error,
not an error in the publication. The original figure crosses out ISPCR; text
extraction alone loses that strike-through.

```text
FS dT30VN  /5Biosg/AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN
FS TSO    /5Biosg/AAGCAGTGGTATCAACGCAGAGTACrGrGrG
```

🟡 The amplified cDNA reference is derived from these two primers. Its inward
sites can be checked on that reference, but this does not prove RNA-primer
extension, every first-cycle intermediate or the suppressive folding mechanism.
The RT/terminal-addition/template-switch intermediates remain explicitly unresolved.

The local graph selects **standard FS with commercial Nextera XT, NT neutralization
and NPM enrichment**. It retains pre-tagmentation cDNA cleanup (0.8:1 SeraMag), PicoGreen QC and
normalization to 100–200 pg/µl, and the pooled-library 0.8:1 cleanup (JATS lines
155–161, 182). Index primers are already present in the initial NPM extension
mixture, and the same mixture continues through index PCR. The in-house Tn5/KAPA,
plexWell, FS-LA and FS-UMI routes are not mixed into this branch.

## Endpoint and sequencing

🟢 The article-linked author v3 ordering table gives complete biotinylated, terminal-phosphorothioate
N714 and S502 primers. Their index members are respectively `TCATGAGC` and
`CTCTCTAT`. The page uses them **illustratively**, without assigning them to an
experimental well or inferring a plate map.

🟡 Final primer arms are derived locally from these complete oligos. Canonical
Nextera ME segments in the commercial endpoint are explicitly inferred by their
constructor and marked on the page; the commercial loaded transposome and its exact
covalent intermediates have not been independently source-checked here. The separate
in-house ME table is not treated as proof of commercial reagent identity. Sequencing
sites on those inferred segments remain conditional on that assignment.

🟢 The defining Methods report standard NextSeq550 **75-8-8-0** acquisition: Read1
75 cycles, i7 and i5 eight cycles each, no Read2. The standard branch has no UMI.
FS-UMI's 100-8-8-50 mode and asymmetric enrichment primers belong to that distinct branch.

Tagmentation follows the represented library class without inventing discarded
products or asserting complete recovery. Other commercial product structures/fates
remain unassigned. Declared cleanup waste remains visible in the graph.

## Published recipe discrepancy

Both readers inspected the actual author v3 commercial tagmentation table on pp21–22.
Its rows give 0.250 µl ATM, 2.000 µl TD and 0.500 µl water: **2.750 µl**, while the
continuation declares **3.000 µl** and says to dispense 3 µl before adding 1 µl cDNA.
The v4 table retains these entries. The defining article instead gives 0.5 µl ATM,
2 µl TD and 1 µl cDNA, with water to 4 µl total. These are distinct printed recipes;
the inconsistent table does not identify which row supplies its missing 0.250 µl.
The public graph represents the commercial route without assigning volumes, and
shows the discrepancy rather than silently repairing a component.

## Remaining source and model scope

- Finish the separate FS-LA/FS-UMI v3 chemical reconstructions. Their original
  revisions are now verified, but downloading or recognizing identical source bytes
  does not perform an independent chemical reading.
- Follow the article's actual cleanup-method chapter, code repository, study
  metadata PRJNA816486 and chemistry-relevant linked documents. Comparator studies
  and processing links need explicit relevance dispositions, not title-only exclusion.
- Reconstruct the separately published FS-UMI/FS-LA variants and any further
  source-supported commercial transposome/primer geometry before claiming that scope.
- Per-well index assignment is unassigned model data. No identified unavailable
  “plate map document” or confirmed publication absence is asserted.

The source review remains in progress. No missing-information claim is accepted.
