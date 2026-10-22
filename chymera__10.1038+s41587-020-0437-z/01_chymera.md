# CHyMErA: two-round pooled cloning and sequencing

Defining publication: Gonatopoulos-Pournatzis T et al. (2020),
[Genetic interaction mapping and exon-resolution functional genomics with a hybrid Cas9–Cas12a platform](https://doi.org/10.1038/s41587-020-0437-z).
Nature Biotechnology 38:638–648.

Detailed recipe modelled here: Aregger M, Xing K, Gonatopoulos-Pournatzis T (2021),
[Application of CHyMErA Cas9-Cas12a combinatorial genome-editing platform for genetic interaction mapping and gene fragment deletion screening](https://doi.org/10.1038/s41596-021-00595-1).

## Scope and source identity

The graph implements the 2021 pooled two-round cloning and nested-PCR readout.
It does not relabel that detailed recipe as a verified copy of the 2020 Methods.
The 2020 main article remains unretrieved: the author-hosted PDF returned 404,
and the publisher PDF route returned subscription-preview HTML. Those route
failures do not establish source absence. Source-owned 2020 Supplementary Table14
independently matches the six named primer rows used by the model. Its A502-Cpf1-F1
name becomes A502-Cas12a-F1 in the 2021 table; their sequences agree.

Read the 2021 published-version author manuscript PMC12314557, own DOI/title
verified: Experimental Design on pooled cloning and NGS primer design; pooled
oligo template and amplification (text584–690); intermediate digestion and
second assembly (text790–879); nested PCR and acquisition (text1330–1436);
Table1 primer rows (text2321–2448); Figures4/5 captions (text2315–2319).
The actual Figure4/5 graphics have not been inspected: PMC's ordinary page is
gated, and publisher figure routes did not return a readable original.
The source original is the XML, not an upstream protocol drawing.

Read the 2020 publisher preview for article identity, supplementary inventory
and its explicitly named GEO series. The published Supplementary Table14 is
MOESM16 XLSX: inspected the six relevant amplification/readout primer rows.
Other guide or assay rows are not assigned to this model. An earlier pass
inspected only the 2020 supplementary PDF's title/DOI and extracted captions;
its original molecular graphics remain uninspected.

Retrieved the official GSE144281 series text after the publisher identified it.
Read its title, publication ID, summary and overall-design identity fields only;
its sample-level preparation methods remain a follow-up target. No raw data
were downloaded. Ran `tools/scrape_primers.py` on the complete bundle. Its
large-workbook exclusions remain explicit; they are not a negative primer audit.
No complete or independent source review has been claimed.

## Molecular construction

🟢 The 2021 source's synthetic oligo preserves one 20-base Cas9 spacer and one
23-base Cas12a spacer, separated by the source-written BsmBI-containing stuffer.
These remain distinct structured guide-identity features. No gene target or
editing design is supplied. The two guide windows are not deduplicated even
if a study-specific pair were to share sequence.

🟡 TKO_F1/TKO_R1 sites are located on the source oligo and primer tails are
added computationally. An independent `plasmid.amplify()` calculation agrees
with the declared amplicon and the source's 157-bp expected band.
The starting pool is drawn as single-stranded DNA. Complementary-strand synthesis
within the reported PCR is unexpanded, so the graph does not assert an already
duplex input or a complete first-cycle binding state. The independent calculation
checks primer ends and tails on the reference axis only.

🟢 A separately opened, rSAP-treated recipient receives the amplified insert in
assembly round one. Physical cleanup and clone recovery precede BsmBI opening
of that intermediate. A separate TOPO donor supplies the tracrRNA/direct-repeat
interval in assembly round two. The donor is a separate molecule until assembly;
it is not covalently appended merely because it shares the reaction tube.

🟡 `opened_intermediate()` derives the known recipient-side BsmBI ends from the
source-written internal stuffer. An independent digest checks those two sites,
not all possible sites in a complete recipient map. Donor-side compatibility
and complete vector maps are unassigned, so matching donor junctions are marked
inferred. No complete cohesive-assembly product relation or quantitative yield
is claimed. BfuAI release geometry remains explicitly unassigned rather than
being reconstructed from a similar method.

The downstream assay and genomic recovery form an explicit unmodelled boundary.
The model provides no delivery, integration, cellular editing or phenotype recipe.
Outer-PCR primer sequences are exact, but genomic landing coordinates and
outer-to-inner intervals remain unassigned. The intended inner-PCR sites permit
a local primer/site consistency check without establishing a genomic map.

## Indexes and acquisition

🟢 The example uses source-listed A502-Cas12a-F1 and A701-Cas9-R1, each without
an added stagger. Their combination is illustrative, not an assigned experimental
sample. Sample i5/i7 features remain distinct from both guide-identity features.
The canonical P5/P7 and TruSeq arms are imported from `illumina`; the named full
primers are checked against the actual 2020 Table14 rows.

Unresolved source read-side assignment: 2021 Step139 names Read1 as Cas12a
and Read2 as Cas9, with 29 dark cycles before 27 Read1 cycles, two 8-cycle
index reads, and 20 dark cycles before 27 Read2 cycles. Step34 orders the
encoded oligo Cas9 then Cas12a, while Table1 and Experimental Design place
i5 on the U6 side and i7 on the termination side. These source statements
have not been reconciled against the original Figure5 graphic and complete
maps. Both source statements and oligos are retained; the model assigns no
corrective recipe and no guide-cycle coverage. Potential canonical primer
sites are displayed without applying run lengths or dark-cycle offsets.
This is an unresolved assignment, not a confirmed publication-error conclusion.

Physical band/clone selection remains separate from the source's later
computational trimming, anchor filtering, mapping and counting. Sequence
consistency does not establish library representation, abundance balance,
guide processing, target specificity, editing, functional effect or coverage.

## Deferred source work

- Defining 2020 main Methods and source-specific applicability of the 2021 recipe.
- Original Figures4/5 graphics and reconciliation of read-side/orientation statements.
- Original pLCHKO and TOPO-tracrRNA-LbDR maps, source-specific ends, complete cuts
  and insert junctions; the named 2021 reagent records are Addgene155048/155049.
- Full genomic primer coordinates and context; GSE144281 sample-level methods.
- Focused single-construct and higher-order variants require their own source scope.
- Full linked-source inventory and independent source reader. No confirmed-missing
  information claim has been made.

## Supplied 2020 article, 2026-10-22

The user supplied the original publisher PDF `s41587-020-0437-z.pdf`; its printed own
DOI/title were verified and its chemistry-relevant Methods subsections inspected:
“Lentiviral hgRNA library construction” (construction portion only) and
“Preparation of CRISPR sequencing libraries and Illumina sequencing”. Ran the primer
scraper on that original and reviewed the relevant construction/readout pointers.
The downstream delivery, target optimization and functional-selection procedures are
outside the current model and were not incorporated.

🟢 The 2020 Methods independently describe the113-nt paired-spacer pool and two
successive construction rounds, with32-nt stuffer, BfuAI-flanked pool, BveI insert
processing, BfuAI/BspMI recipient processing and later Esp3I/BsmBI opening. This confirms
the broad architecture, not identity of every2021 reagent choice, condition or map.
The page continues to identify the exact2021 detailed protocol as its implemented recipe.

🟢 The2020 readout explicitly differs:29 dark cycles then31 imaged bases for Cas12a,
20 dark cycles then30 imaged bases for Cas9, with8-base indexes on each side. Step139
of the2021 protocol instead gives27 imaged bases for each side. The page now presents
this as version-specific acquisition, retaining the2021 values and the unresolved
read-side orientation. A later recipe must not silently overwrite the defining article.

The previously requested2020 main is now available and partially reviewed; remaining
work is exact component/map applicability, original figure/read-side reconciliation and
independent full review. No complete review or missing-information conclusion is implied.
