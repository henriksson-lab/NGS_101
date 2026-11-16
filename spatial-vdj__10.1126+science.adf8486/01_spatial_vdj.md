# Spatial VDJ: published SR-TRB preparation

Engblom, Thrane, Lin et al., *Spatial transcriptomics of B cell and T cell receptors reveals lymphocyte clonal dynamics*, [Science382,eadf8486](https://doi.org/10.1126/science.adf8486) (2023).

## Source identity and resumed scope

🟢 The primary source is now the **final published article**, [Science, DOI 10.1126/science.adf8486](https://www.science.org/doi/full/10.1126/science.adf8486). An authorized visible browser displayed full institutional access through Umeå University; its full article HTML was saved, source-manifested, validated and converted in the external store on 30 October 2026. Inspected: Results1580 (constant-region capture); spatial GEX Methods1623–1624, capture1625–1642, LR boundary1643–1644, selected SR preparation1645–1667 and supplementary index1767–1782 of `science_adf8486_final_browser.html.txt`. The mandatory bundle oligo scraper was rerun; sequence-like analysis barcodes and bulk-validation oligos are not reassigned to the spatial preparation.

The final published Methods directly support the selected short-read TRB workflow below. This establishes the inspected preparation against the publication; it does not assert equivalence of every accepted-manuscript or supplementary detail. The earlier [Zenodo13148082 accepted author version](https://zenodo.org/records/13148082) and its independent review remain historical evidence. No preprint or biological data archive supplies chemistry.

The user supplied both published downloads on 30 October 2026. The supplementary PDF and intact Data S1–S5 archive were moved from `/tmp/camilla` to the external source store, checksummed, manifested and revalidated. The ZIP contains five named XLSX workbooks. Only **Data S1** was selected for chemistry inspection; Data S2–S5 were extracted and converted but not read as chemistry evidence.

Inspected supplementary PDF: published title/DOI; Methods1–72 distinguishing the barcode-losing validation PCR and capture-enrichment evaluation from Spatial VDJ. The selected spatial oligos come from the actual Data S1 worksheets, not the separate bulk SS3 or validation recipes. The graph still begins with prepared full-length spatial cDNA, without inheriting a Visium recipe.

## Selected physical preparation

🟢 Choose the **short-read TRB** route, with the second experiment's single14-cycle postcapture PCR. TRA's14-cycle outer-PCR alternative, LR SMRTbell preparation, bulkSS3 validation and the FigS1F V/C screen are separate scopes. In particular, the FigS1F screen loses the spatial barcode; it is not substituted for SR-SpatialVDJ.

1. Start with full-length spatially barcoded Visium cDNA **before fragmentation**. Select10µl from the library, divide into10ng PCR subsamples whose number depends on concentration, then amplify with25µl KAPAHiFiHotStartReadyMix2X/7.5µl10xcDNA primers/17.5µl sample and water **per reaction**. The graph preserves a variable array of separate reactions through this opaque interface; it does not clone a single PCR or pool templates before amplification. Pool completed products, clean with0.6xSPRIselect and elute40µlQiagenEB.
2. The source delegates concentration to xGenv4AppendixA, then PlateProtocolstep13 onward. Coenrich IG and TRA/TRB constant-region transcripts with custom adapter blockers and IDTxGenHybridizationandWashKit1080584,16h overnight. ChooseIG:TR1:3. This is achieved by **triplicating TR probe sequences in the order and adding equal volumes of the produced IG/TR pools**, not by adding three volumes of TR pool. The graph does not invent bait biotin or a streptavidin formulation from the kit name.
3. Use the selected14-cycle postcapture amplification and purification. Four actual TRB outer-PCR replicates each receive2.5–5ng captured cDNA. OuterTRBV primers and partRead1 amplify for12cycles; source primer stock-dilution wording is100µM→40x→2.5µM, with final-reaction concentration convention unassigned. PreserveqPCR calibration/Qubit/Bioanalyzer analytical allocations separately. Clean each outer product with0.6xAMPure/two80%EtOH washes and eluteEB after15min37°C.
4. Select5ng from each outer stock for its own innerTRBV/partRead1 PCR. Run8cycles; retain the printed72°C30s followed by72°C90s in the cycle programme. Purify and quantify each replicate separately. Select10ng from each of the four inner stocks; **pool at the indexing stage**, assigning the same index to the four replicates.
5. Index-PCR8cycles, bead-purify/quantify and select sequencing allocations. The source locally names KAPAReadyMix for common/outer/inner reactions; its exact indexing formulation/packaging is not assigned from that earlier name. Aggregate index-PCR enzyme/nucleotide/buffer supply classes remain unassembled. Samples are pooled for NovaSeq6000SP withR1=35nt, one6-cycle index read andR2=259nt. Primer landing, index orientation and exact amplicons remain unresolved.

Known ReadyMix components are exposed from each actual named ReadyMix allocation; they are not fresh extra enzyme/nucleotide additions. Stock remainders and analytical aliquots are kept apart. Amounts not specified by the source remain unassigned.

## Sequence and validation boundaries

The final published article, “LR Spatial VDJ read processing” (text1685), gives two short computational search patterns (`CGACGCTCTTCCGATCT` and `TCTGCGTTGATACCACT`) in **LR analysis**, followed by16spatial-barcode and12UMI bases. These are source analysis locators, not complete spatial primers or permission to import Visium/TruSeq constructs. No exact source oligo state is manufactured from a platform name or matching anchor. Spatial barcodes identify tissue spots, not single cells or physically paired receptor chains; receptor pairing is computational.

## Published Data S1 oligo inventories

🟢 Data S1 now supplies the complete source sequence fields for the selected preparation:

- `SR-Spatial VDJ PartRead1 primer`: the whole22nt PartRead1 order, source line4.
- `SR-Spatial VDJ Inner TRBV`:62 named whole orders and their explicit target/handle columns, lines55–119. All whole orders equal the published34nt Read2 handle followed by that row’s target sequence. The constructor checks this relation rather than guessing a handle from a kit name.
- `SR-Spatial VDJ outer TRBV`:62 named rows, lines171–235, with58 distinct sequence bodies. The worksheet’s interior caption says “Outer TRAV” while the worksheet name and individual IDs say TRBV; the model follows the worksheet/row identities and retains this label discrepancy. It does not silently relabel the original.
- `Hyb. capture - Probes`:174 Ig rows and20 TCR rows, all120nt, lines236–431. Ig has169 distinct bodies. Every row is preserved, including repeated sequences. Source row multiplicity does not independently establish pool molarity or the experiment’s separate triplicate TR-order instruction.
- `Hyb. capture - Blocking Oligos`:both whole patterns, lines432–435. Read1+polyT is82nt:PartRead1(22)+N28+T30+VN. The N/V positions are a degenerate blocking-reagent pattern, not assigned sample barcodes or UMIs. The second blocker, TSO, is30nt.

The source sequence cells do not assign terminal chemistry or occupied binding sites; no5′phosphate, blocked3′end, probe biotin or concrete N/V value is inferred from plain sequence text. These are oligo-pool inventory views, not annealed multimolecular complexes. The four actual replicate allocations retain their respective whole outer/inner pools and PartRead1.

The selected whole index-primer orders, historical capture-kit geometry, captured transcript sequences and complete downstream molecular products remain unassembled. Preamplification/postcapture PN2000089 primer sequences and the input cDNA are now supplied by the exact cited Visium RevA declaration, as described in the current follow-up below. This is an open modelling/source question, not a conclusion that those details are unpublished. The28N blocker window alone does not prove the16/12 barcode/UMI partition; the final article’s LR parsing independently specifies those lengths. This is an analysis-coordinate statement, not a fully assembled SR read map. No bulk SS3 or barcode-losing validation primer is substituted.

The final published preamplification programme (line1632) literally says “Repeat steps 2 to 5 six times for a total of five cycles.” Both printed values are retained, and executed cycle count is unassigned. The same ambiguity was previously observed in the accepted manuscript. The model does not choose a corrected cycle count.

The public source-qualified workflow records the final article’s positive source order and stock/fraction boundaries. It remains an incomplete library model: selected cDNA/index primer preparations, occupied transcript geometry, molecular products and primer landing remain unresolved. No complete source audit or accepted missing-information conclusion is claimed.


## Reader-facing page revision (31 October 2026)

🟢 Re-read the final published HTML Results1574–1580 and receptor pairing1601; capture/SR preparation1623–1668; LR parsing1685 and SR clonotype processing1693–1700. Re-read Data S1 PartRead1 line4, innerTRBV55–62, outerTRBV171–177 and blocker432–435. The mandatory bundle scraper was rerun. No new original was downloaded and no preprint recipe was adopted.

The public page now explains the biological question, why the cDNA must be taken before fragmentation, what constant-region capture and seminested V-region PCR select, the four-replicate pooling boundary, oligo roles, run lengths, identifiers and the clone-by-spot output. The additional views use existing source-order constructors for one outer/inner TRBV example and both blockers; they assert no occupied transcript geometry. Complete source tables and the full graph remain accessible in disclosures.

🟢 The final published LR-processing paragraph independently verifies the16nt spatial barcode and subsequent12nt UMI; the previous accepted-version-only qualification for these analysis lengths is superseded by this specific final-article reading. No complete primer or exact SR cycle map follows from that paragraph. Final SR analysis explicitly combines barcode and UMI for collapsing, requires three supporting reads in tonsil, and describes index-hopping/overamplification filtering. Receptor pairing remains computational; the source explicitly says individual captured chains do not establish receptor pairs.

This was a bounded presentation/source-excerpt pass, not an additional independent missing-information audit. Historical Visium/xGen kit geometry, index-primer orders and complete amplicons remain deferred model questions. The existing publication discrepancy notices are retained on the schematic.


## Primer-design and upstream-input follow-up — 1 November 2026

Historical checkpoint: the Ståhl source blocker below was resolved on 16 November; see the source-linked indexing update at the end.

🟢 Re-read final published Methods1623–1668, DataS1 PartRead1/innerTRBV10-1/outerTRBV10-1
and explicit handle/target columns. Source primer scraper rerun. The requested exact
VisiumCG000239RevA original was retrieved to the external **Visium** source bundle,
its printed revision verified, and Steps4.2–4.3 (printedpp49–50) and original Appendixp75
inspected. The copied TSO end includes the source `ACAT` DNA extension beyond the
generic SMART handle; independent source-transcription checks protect that detail.

The SpatialVDJ graph now imports the actual source-owned full-length, pre-fragmentation
cDNA export from the Visium page. The requirement pins **RevA**, not that page's default
RevF sequencing route. It rejects the wrong revision, a fragmented/indexed input,
mutated molecular declarations or loss of the fixed PartRead1/cDNA-reverse primer sites.
This establishes a checked handoff, not validation of unknown native reaction kinetics,
upstream released-strand identity or all source-specific products.

🟢 The RevA appendix identifies cDNA PrimersPN2000089: forward equals source PartRead1;
reverse is the partialTSO primer. Those whole sequences are now declared for both common
pre/postcapture amplifications. The public page walks through a source TRBV10-1 outer
order and the inner order as source Read2 handle plus source3′target sequence, with
whole ordering orientation, model-derived lengths, roles and a computed universal
primer-binding view. No occupied V-gene coordinate or invented amplicon is supplied.

Remaining: source-specific indexing primer orders/orientation, occupied receptor
transcripts, complete captured products and exact capture-reagent/terminal chemistry.
The final indexing paragraph directly cites Ståhl2016 (reference6,
DOI10.1126/science.aaf2403). That original's existing external manifest records a manual,
unavailable source with no cached original. Its indexing chemistry is deferred; a
modern Visium/TruSeq index recipe is not substituted. This is not an absence claim.

## Step-by-step presentation follow-up — 16 November 2026

🟢 Re-read the final published target-enrichment and SR Spatial VDJ Methods,
text1625–1668, together with the existing DataS1 primer declarations. The source
primer scraper was rerun over the external bundle. No additional molecular source
or historical kit version is presumed read in this presentation pass.

The public page now leads with one connected molecular Workflow rather than a
six-row summary and a hidden audit graph. It shows the source-owned Visium input,
preamplification, post-PCR pooling/concentration, capture, postcapture PCR, outer
TRBV amplification, per-replicate cleanup, inner amplification, replicate pooling,
index PCR and sequencing. Reaction notes supply the published inputs, conditions,
replicate boundaries and the unresolved preamplification-cycle discrepancy.

🟡 The receptor-specific output drawings are **expected region-level classes**,
not new validated source-instance amplicons. They preserve the exact upstream
spatial-barcode/UMI feature identities. The unknown receptor interval is labelled
constant / J / CDR3 / retained V in the direction of the displayed tag-side DNA
strand, which is opposite to the receptor RNA. It has no assigned nucleotide length.
Outer/inner primer-derived ends and the copied Read2 handle are computed from the
published source primer orders with the appropriate reverse complement. A literal
V allele, occupied V-priming complex or exact CDR3 coordinate is not fabricated.
The paired strands come from Scene rather than manually aligned sequence columns.

The inner-product duplex represents the bounded class reached after copying both
primer ends. The first inner V primer initially carries an unpaired5′ handle;
its complementary copy cannot appear in the same initial binding event. This
walkthrough explains that distinction without claiming a full reaction simulation.
At this checkpoint, Index/P5/P7 extensions were unassembled role tokens; the source-linked update below supersedes that boundary. The process graph
and its unresolved transitions remain available separately; these explanatory
drawings do not increase the protocol's chemical-validation completeness score.


## Can sequencing reads resolve the index-primer ends? — 16 November 2026

🟢 **Public versus controlled data.** The final published article, Data and materials availability (text1764), deposits tonsil raw FASTQs in [EGA EGAD00001011062](https://ega-archive.org/datasets/EGAD00001011062), under DAC EGAC00001003294. Its public metadata requires an access application. No controlled files were accessed. The accompanying [Zenodo7961605](https://zenodo.org/records/7961605) record supplies one 25,609,313,924-byte archive, `SpatialVDJ_forZenodo_resubmission.zip` (record MD5 `e5fb109059ea51d46a48dcf174300902`).

🟢 **A concrete public read-derived lead remains.** A bounded HTTP Range inspection of the ZIP64 central directory (2,965 entries;587,380 directory bytes, plus131,072 trailer bytes) identifies these members below `SpatialVDJ_forZenodo/data/tonsil/2_SR-SpatialVDJ/`. No FASTQ/BAM/SAM/CRAM-named member was listed in that outer directory inventory; this is a bounded archive-listing observation, not a source-wide absence finding. The two nested ZIP entries are analysis-code archives, and their contents were not inspected in this pass.

| Candidate member | Compressed bytes | Uncompressed bytes | Scope |
|---|---|---|---|
| `TRA/Tonsil_merged.vdjca` | 2434678 | 2734593 | Public MiXCR alignment candidate; read-bearing content not downloaded |
| `TRB/Tonsil_merged.vdjca` | 11008212 | 12248978 | Public MiXCR alignment candidate; read-bearing content not downloaded |
| `TRA/Tonsil_merged.matrix.tsv` | 456551 | 4047600 | Public exported alignment-table candidate; rows not inspected |
| `TRB/Tonsil_merged.matrix.tsv` | 1807200 | 19167636 | Public exported alignment-table candidate; rows not inspected |

🟢 **What “original reads” means here.** The final published SR analysis (text1693–1698) first filters/pRESTO-processes the short reads, checks spatial barcodes/inner primers and collapses UMI groups; then MiXCR3.0.13 is invoked with `-OsaveOriginalReads=true`. Its reported table export requests `-targetSequences`. The small public `TRB/Tonsil_merged.report` was read by exact HTTP ranges and validated against its ZIP CRC32 `fe9b1c96`: it confirms MiXCR3.0.13, inputs `Tonsil.R1.merged.fastq` / `Tonsil.R2.merged.fastq`, `-OsaveOriginalReads=true`,43,999 input reads and42,863 aligned reads. Such saved sequences would be original **MiXCR-input** sequences, not automatically the original untrimmed sequencer FASTQs; the actual read-bearing binary members have not yet been inspected. The public SR readme and downstream analysis script were also inspected; they describe MiXCR outputs and count analysis, without establishing a complete upstream trimming recipe. Whether a discriminating adapter tail survives that preprocessing remains open.

🟡 **Actionable next check.** First inspect a bounded, manifested sample of the public `matrix.tsv` target-sequence fields (or a version-compatible bounded `.vdjca` read export), without downloading the25.6-GB archive. Verify export headers, paired-read identity, adapter/primer trimming and UMI-consensus provenance before counting short-fragment read-through. Seek reads that retain both a genuine receptor/inner-primer context and an opposite-end adapter flank; compare both orientations, support across independent spatial-barcode/UMI groups and base quality where retained. Adapter dimers and off-target short products cannot establish the intended receptor-library product. If preprocessing removed those bases, the necessary source is an authorized small sample of pre-trim R1/R2 with sample-sheet/index information or original indexing-oligo records; controlled EGA access is an actual access requirement, not a route to bypass.

🟡 **Identifiability limit.** A supported read-through consensus could establish copied terminal sequence of the final library. It cannot by itself uniquely identify which original primer synthesis order, overhang arrangement, terminal modification or chemically equivalent PCR design produced that sequence. Ordinary R1/R2 also do not automatically retain the six-cycle index read or both outer flow-cell ends. Any future sequence inferred from reads needs a read-derived provenance label and its exact scope; a matching TruSeq endpoint alone does not resolve the source-specific whole index-primer order. At that checkpoint the P5/index/P7 extension roles remained unassembled; the subsequent original-source update below resolves their ordering patterns.

This was a bounded feasibility pass over the final article, its directly linked public data record and archive metadata/small processing-text members, and the directly linked EGA metadata. It is not a complete source audit or a confirmed absence review. No biological read dataset, third-party archive or source material was placed in the checkout.


### Separate bounded public-target-sequence follow-up

🟢 A separate follow-up inspected the first 100 **processed TRB alignment rows** from the exact public Zenodo member `SpatialVDJ_forZenodo/data/tonsil/2_SR-SpatialVDJ/TRB/Tonsil_merged.matrix.tsv`. The protocol-local standard-library sampler `tools/sample_public_alignments.py` supports `--max-reads` (1–100 processed rows) and a strict compressed-prefix cap (maximum131072bytes), checks exact HTTP206 intervals and ZIP member identity, and writes the sample plus manifest outside the checkout. No full ZIP or `.vdjca` was downloaded. The sampled member prefix is131072bytes; its whole member CRC from the central directory is retained as provenance but was **not** presented as verified from a partial download.

The external sample, manifest and aggregate summary are under `$CHEM_DATA/read_samples/spatial-vdj__10.1126+science.adf8486/index-tail-probe-2026-11-16/`, with basename `TRB_first_100_processed_alignments`. The sample SHA256 is `4dab7e97e6e6e17dea09d44ad66931c0d456b969e4ab46d37034c5e9c003dc62`. Only the TSV `targetSequences` fields were inspected, not saved-input VDJCA contents or raw EGA reads. The presence of a `cloneId` column does not by itself identify the target sequences as clone consensus. Its actual columns are `cloneId`, `descrsR2`, `targetSequences`, `nSeqCDR3` and `bestJGene`; base qualities and a separate index read are not present in this export. The100 rows contain200 target strings:100 of35nt and100 of259nt. Mate acquisition length is not an insert-size measurement, and MiXCR target ordering/orientation must not be assumed to be raw R1/R2 order.

🟡 The bounded sample has no exact full match, or exact12/16/20nt window, to canonical P5, P7, TruSeq Read1/Read2 or source PartRead1 in either orientation. This means **no informative technical read-through was found in these first 100 rows**; the sample is not random or exhaustive and cannot establish absence in the study or full public archive. It does not reject the strategy for genuine shorter inserts. No source-specific terminal primer sequence is inferred from this negative sample, and no index-primer/P5/P7 chemistry was assigned from this data sample. The original-source update below now supplies the primer patterns.

To reproduce the bounded retrieval, run `python3 spatial-vdj__10.1126+science.adf8486/tools/sample_public_alignments.py --chain TRB --max-reads 100`. The script keeps raw/processed data and derived summaries in the external store and refuses a source-identity mismatch when reusing the output filename. Further evidence needs informative short-insert pre-trim reads, a version-qualified retained-read export with qualities/processing provenance, the actual sample sheet/index information, or the original indexing-oligo records. Equivalent primer designs yielding identical copied sequence remain indistinguishable by final-read sequence alone.


Sampler robustness follow-up: reuse now checks the declared member CRC, compressed/uncompressed sizes and layout metadata as well as URL/member identity, verifies the existing sample against its manifest SHA256, and refuses a different incoming sample. Suffix responses must match their declared archive total. An incomplete final byte-line is removed before UTF-8 decoding. Run `python3 spatial-vdj__10.1126+science.adf8486/tools/selftest_public_alignments.py` for the offline regression suite; its mocked responses and temporary fixtures require no network or downloaded biological data. The separate CI sampler step runs these checks without adding them to the chemistry-suite count.


## Source-linked Ståhl indexing resolved — 16 November 2026

🟢 User-supplied `science.aaf2403.pdf` and `aaf2403_stahl_sm.pdf` were moved into the canonical external `sources/doi_10.1126+science.aaf2403/` bundle, checksummed, source-manifested and revalidated. Both originals print DOI **10.1126/science.aaf2403**, *Visualization and analysis of gene expression in tissue sections by spatial transcriptomics*, Science353,78–82 (2016). The main PDF begins with the preceding article’s tail; its own title, text and DOI establish identity. The supplement cover identifies Materials and Methods, FigsS1–S11 and TablesS1–S2, with original locator `www.sciencemag.org/content/353/6294/78/suppl/DC1`. User download endpoints were not supplied; the manifest records the printed/DOI locators rather than inventing an endpoint.

Inspected: main article’s array-capture → IVT library explanation and supplementary index; original supplement PDF pp7–9, library preparation and full PCR oligo list. The supplement text layer has a shifted/garbled font mapping. Mandatory primer scraping first returned zero hits, **not evidence of absence**. Clearly labelled OCR of these pages was appended to the PDF’s text twin; the primer orders, six Xs and precursor qPCR composition were checked against rendered original pp8–9, then the scraper was rerun (four hits). No OCR-only sequence was accepted. No recursively linked predecessor recipe or complete missing-information audit is claimed.

🟢 The final Engblom indexing Methods explicitly says TruSeq indexes “as described previously(6)”; reference6 is this Ståhl paper. Its supplement p9 lists, verbatim 5′→3′:

- **PCR primer InPE1.0:** `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`.
- **PCR primer InPE2.0:** `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`.
- **PCR Index primer:** `CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTC`.

These are **three** oligos, not a guessed two-primer kit. The model imports canonical constant sequences and declares the selected six-base index as an unknown sample-index feature. XXXXXX denotes the chosen index member; it is not an N-randomized synthesis, a spatial barcode or a UMI. The selected index sequence for each 2023 sample remains unassigned. No additional i5 index is introduced.

🟡 InPE1.0 contains the retained22-base PartRead1 at its 3′ end and contributes the upstream36-base P5/Read1 extension. InPE2.0 is the complete34-base handle already supplied by the inner V primer; the PCR Index primer has P7(24) + selected index(6) + the first13 bases of that handle, for43 bases total. Scene computes the forward and index-foot binding against the conditional inner-product class. InPE2.0 and PCR Index share the same template end; they are not shown simultaneously occupying overlapping sites. The expected later bounded duplex now carries explicit source-derived arms; the displayed tag-side strand has the reverse complement of the P7/index contribution. This is a conditional copied-product class, not a full occupied receptor reconstruction or validated physical PCR population.

🟢 Precursor supplement p8 specifies a10µl qPCR containing1×KAPA HiFi HotStart ReadyMix **KK2601**,1×EVA green, InPE1.0 **0.5µM**, InPE2.0 **0.01µM**, PCR Index **0.5µM** and2µl purified cDNA. Its programme is98°C3min, then98°C20s/60°C30s/72°C30s; qPCR selects cycling, followed by25µl amplification. These are **2016 precursor conditions**, not independently verified 2023 reagent doses. Engblom provides its own8-cycle programme with72°C90s extension and72°C5min final; that local programme remains authoritative. The 2023 graph retains unassigned index-PCR stock identities, doses and formulation rather than silently assigning the2016 ReadyMix/EVA green recipe.

Only the explicitly delegated indexing design is reused. Ståhl’s older surface-capture/IVT/RNA-ligation preparation does not replace the source-owned **Visium CG000239 RevA pre-fragmentation cDNA** input. Existing preamplification-cycle and DataS1 worksheet-heading discrepancies remain visible. Historical read-sampling results above are retained as limited negative results; reads are no longer needed to recover these three source-published primer patterns. Historical retrieval request171 in `geturl.txt` is resolved. Remaining local work concerns selected sample indexes, receptor occupancy/full amplicons, terminal chemistry and the exact historical xGen capture original; no missing-information conclusion follows.
