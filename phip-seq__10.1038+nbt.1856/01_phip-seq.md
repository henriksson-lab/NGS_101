# PhIP-Seq: enrich displayed peptides, sequence their phage DNA

Larman et al. (2011), *Autoantigen discovery with a synthetic human peptidome*, [defining article](https://doi.org/10.1038/nbt.1856).

This is a capped first pass and a private provisional model, not a verified reconstruction of the final published procedure. The genotype/display association and assay branch are represented; the complete synthetic-library cloning, packaging and propagation remain unmodelled. Exact PCR geometry and final procedural applicability remain open.

## Originals inspected and version limits

Three originals, 28 October 2026:

- Associated accepted-manuscript JATS `nbt.1856_PMC4169279.xml.txt`: source identity and manuscript flag; Introduction/Results60–93; assay/background101–103; Discussion127–135; Fig1 caption353. Its Methods139–141 points to online methods. It is not the version-of-record procedural text.
- Author supplement `nbt.1856_supp_NIHMS284271-supplement-1.pdf`: Supplementary Methods1055–1088, detailed protocol1107–1178. Original PDFpp14–15 were rendered and inspected for colored/underlined primer bases and the image-based PCR table. These are author-supplement details; applicability to the final publication has not been established.
- Publisher supplement `nbt.1856_41587_2011_BFnbt1856_MOESM57_ESM.pdf`: document identity, page count13, and rendered p13 Table3/reference inspected. Remaining pages were not fully read. Its contents are not declared equivalent to the longer author supplement. An empty oligo scrape of this file does not establish absence of information from the publication.

The mandatory bundle scraper was run, including `--all` after the publisher supplement arrived. Ten oligo hits in the author supplement were checked against its text and the rendered p14 original. No missing-information claim is accepted. Source investigation stopped at the three-original cap; no cited references were followed.

## What the assay links

🟢 The associated article describes T7-Pep: a synthetic peptide-display library, with each displayed peptide associated with its encoding phage DNA. Immunoprecipitation enriches a population of particles; PCR and sequencing count their DNA. This does not mean every retained particle bound a specific patient antibody. The source explicitly measures substantial bead-associated background and compares against bead-only controls.

🟡 The local `DisplayParticle` state retains that qualitative DNA/display association through incubation and selection. It neither translates an unassigned insert nor predicts affinity, display copy number, genotype integrity or capture yield. Capsid protein is not joined covalently to the drawn DNA.

## Provisional author-supplement assay order

1. 🟢 Preblock tubes with BSA. Incubate T7-Pep library with patient antibody, including the named SAPK4 positive-control antibody, at4°C for24h.
2. 🟢 Add a1:1 mixture of ProteinA and ProteinG magnetic bead slurry; incubate4h at4°C. Wash six times with150mM NaCl,50mM Tris-HCl,0.1% NP-40 pH7.5. Only the actual removed supernatants/washes become waste.
3. 🟢 Resuspend beads in30µL water and heat90°C10min to disrupt phage and release DNA. Transfer the aqueous DNA template for PCR. The model does not assert complete protein removal or invent a DNA extraction kit.
4. 🟢 TaKaRa HS Ex RR006A PCR adds the printed forward barcode family and reverse primer. The rendered p15 table states98°C10s,56°C15s,72°C25s, repeated for40cycles, then72°C7min; text permits increasing to45cycles. Exact unamplified phage flanks and introduced mismatches remain unassembled, so the PCR geometry is explicitly unavailable.
5. 🟢 Gel-purify each product separately, quantify and combine500ng from each barcoded sample. Use the custom forward sequencing primer below. The source does not establish a total acquired read length in the inspected section; none is assigned here.

## Author-supplement oligos

All sequences are5′→3′ with no chemical modifications printed for these oligos. They remain provisional author-supplement transcriptions. The forward family contains a20nt P5-related prefix, the custom sequencing-primer site, three variable bases and fixed GCGC. Do not replace its historical prefix or the printed reverse oligo with a full modern canonical P5/P7 primer.

| Oligo | Sequence |
|---|---|
| HsORF-FL-mmBC1-F | {{= phip.forward_primer('ACT') }} |
| HsORF-FL-mmBC2-F | {{= phip.forward_primer('GCC') }} |
| HsORF-FL-mmBC3-F | {{= phip.forward_primer('CCT') }} |
| HsORF-FL-mmBC4-F | {{= phip.forward_primer('TCT') }} |
| HsORF-FL-mmBC5-F | {{= phip.forward_primer('GAT') }} |
| HsORF-FL-mmBC6-F | {{= phip.forward_primer('GGT') }} |
| HsORF-FL-mmBC7-F | {{= phip.forward_primer('GTT') }} |
| HsORF-FL-mmBC8-F | {{= phip.forward_primer('GCG') }} |
| P7-T7Down | {{= phip.REVERSE_PRIMER }} |
| mmBC-FL_seq_prim | {{= phip.READ_PRIMER }} |

🟢 Original p14 colors three bases blue as the barcode and underlines the custom-primer site. Textp15 says to use the first seven calls to parse samples. 🟡 The represented primer-derived order predicts three variable bases followed by the fixed four-base GCGC tail. Thus a seven-call sample-parsing window need not mean seven independently variable bases. This is a conditional layout derivation, not an assertion that a particular study sample used any selected member. The final model uses an unassigned sample-index family, not a UMI.

## Upstream cloning and interpretation remain scoped

The associated article Fig1 describes peptide-coding oligo pools cloned into a FLAG-modified T7Select10-3b derivative, T7FNS2. The author supplementp13 names EcoRI/SalI insertion and NotI/XhoI vector modification. Exact vector sequence, oligo amplification inputs and packaging/propagation are not reconstructed. Its in-silico subsection describes72nt tiles whereas associated Results describe36aa peptides. Those source/version statements are retained as unresolved; no publication error or whole-document equivalence is established.

No registered predecessor is adopted as a molecular constructor. T7Select10-3b is a named platform lineage question, not a new authoritative document invented for a desired component. Similar Illumina handles do not prove ancestry from another NGS protocol.

## Actual deferred targets

Version-of-record article/online Methods at the defining DOI; final procedural applicability of the author supplement and the remaining publisher-supplement pages; original T7FNS2/vector and oligo-pool sequence identities; exact mismatch-introducing PCR geometry; acquired instrument/read configuration; final sample-to-barcode assignment and independent second reader. These are unperformed checks or unassigned values, not confirmed omissions.

## Round04 final-publication retrieval follow-up

The exact publisher article HTML was retrieved and inspected: it exposes identity and abstract but explicitly says it is subscription-preview content. Its `citation_pdf_url`, `https://www.nature.com/articles/nbt.1856.pdf`, was then requested through the standard fetcher; the response was HTML, not a PDF, and was rejected. The final procedural original remains a real access blocker. The bundle scraper was rerun; no new chemistry evidence or author/publisher equivalence was established. The provisional model and all existing limits remain unchanged.

## Unlinked-page adjudication: required preparation scope retained

The bounded read-only pass on 30 October 2026 revisited only the existing
accepted-manuscript main (identity, Methods pointer and Figure 1 caption) and
the associated author supplement, lines 1055–1088 and 1107–1178. No publisher
original, vendor source or reference was added. The final-published procedure
and author-supplement applicability remain qualified by the existing actual
rejected publisher-PDF response; this pass makes no new availability claim.

The current normal graph gate still reports the unaccounted authoritative
source item `upstream_display_library`. That is a real preparation obligation,
not coverage to exclude merely to publish a coarse assay endpoint. The author
supplement specifies 23 independent oligo pools for amplification and cloning
into T7FNS2, with subsequent packaging/propagation still unassembled locally.
It also explicitly requires two independent immunoprecipitation reactions per
sample after the initial phage/antibody setup. Those preparations must exist
before selection; a single completed IP or PCR cannot be relabeled or split
into supposed independent preparations.

The detailed assay transfers the entire 30 µl heat-released aqueous DNA into
a 50 µl PCR. Individually gel-purified samples are quantified and contribute
500 ng each to the sequencing pool. Actual independent preparations, retained
stock, analytical allocations and pool contributions therefore need a faithful
model before publication. The current external other-library input does not
prove those source preparations were produced by the modeled lineage.

The model and private page remain unchanged. This is a source-qualified
preparation/model hold, not an accepted absence claim or a completed review of
the remaining published supplement, vector originals or whole library.
