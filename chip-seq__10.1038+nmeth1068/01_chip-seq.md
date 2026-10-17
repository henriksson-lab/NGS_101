# Conventional ChIP-seq — research notes

## 1. What it is

ChIP-seq enriches crosslinked chromatin fragments associated with an antibody target and
turns the recovered DNA into a conventional end-repaired, dA-tailed, adapter-ligated
sequencing library.

## 2. Sources

- Robertson G. et al. "Genome-wide profiles of STAT1 DNA association using chromatin
  immunoprecipitation and massively parallel sequencing." *Nature Methods* 4:651–657
  (2007), DOI [10.1038/nmeth1068](https://doi.org/10.1038/nmeth1068), PMID 17558387.
  🟢 Crossref resolves `10.1038/nmeth1068` (no dot) to this article; the dotted
  `10.1038/nmeth.1068` formerly used by notes and the catalogue does not resolve (both
  are now corrected). Not in PMC. 🟢 Read: the publisher PDF (printed DOI 10.1038/NMETH1068,
  hand-downloaded by the maintainer on 2026-10-17), Methods "Illumina library
  construction and sequencing" and "Processing 1G data" (`nmeth1068.pdf.txt` lines
  750–815), discussion of read length (lines 727–729) and results read counts (lines
  144–157).
- 🟢 Read: Robertson et al. Supplementary Methods (Springer MOESM305 PDF), section
  "STAT1 ChIP DNA" (text lines 545–566), plus its contents list (lines 1–33).
- 🟢 Read: Illumina, *Preparing Samples for ChIP Sequencing of DNA*, part 11257047
  Rev. A (© 2007), whole guide (text lines 1–494). The article follows Illumina's
  "manufacturer's instructions" for end repair, A-tailing, ligation and PCR without
  naming a document; identity of those instructions with this guide is assumed.
- 🟢 Read: Robertson et al. Supplementary Data 1 (Springer MOESM306 ZIP, "Data files for
  significant peaks"), `README.txt` lines 1–60 (read length only; peak tables not used).
- Schmidt D. et al. *Methods* 48:240–248 (2009), DOI `10.1016/j.ymeth.2009.03.001`,
  PMC4052679 (author manuscript JATS XML). A later protocol from a different laboratory,
  not a source cited by Robertson et al. Read only sections 2.2.7–2.2.11 and the
  reamplification primer list (text lines 192–340, 456–460).

## 3. Chromatin preparation (Robertson Supplementary Methods)

🟢 HeLa S3 cells ± 5 ng/ml IFN-γ for 30 min were fixed in 1% formaldehyde for 10 min at
room temperature, quenched with 125 mM glycine, swelled in hypotonic buffer, Dounce
homogenized, and nuclei lysed in RIPA. Chromatin was sonicated (Branson 250, setting 2,
seven 30-s pulses) to an average of 700 bp. Lysates were incubated overnight at 4 °C with
anti-STAT1 rabbit polyclonal (Santa Cruz sc-345); Protein A agarose beads bound for
60 min. Immunoprecipitates were washed three times in RIPA and once in PBS, eluted in
TE/1% SDS, and crosslinks reversed overnight at 65 °C. DNA was purified with RNase A
(1 h, 37 °C), proteinase K (200 µg/ml, 2 h, 45 °C), phenol:chloroform:isoamyl alcohol
extraction and ethanol precipitation with pellet paint.

## 4. Library construction

🟢 (article Methods, lines 751–786) 5–50 ng of immunoprecipitate size-fractionated by
12% PAGE, excising 100–300 bp; overnight elution at 4 °C in 300 µl LoTE/7.5 M ammonium
acetate (5:1); QIAquick → end repair with a 1:5 mixture of T4 and Klenow DNA polymerases
(Illumina, manufacturer's instructions), 30 min at 20 °C; phenol–chloroform–isoamyl
alcohol in phase-lock gel tubes; ethanol precipitation with mussel glycogen and NH4OAc →
single-A addition with Klenow exo− (3′ and 5′ exo minus; Illumina), 30 min at 37 °C;
extraction and precipitation; two 70% ethanol washes; 10 µl LoTE → ligation of Illumina
adapters to the A-tailed DNA per manufacturer's instructions, except the adaptor
oligonucleotide mix diluted 1/10; QIAquick with 30 µl elution buffer at 50 °C → PCR
enrichment with Phusion (Finnzymes) and PCR primers 1.1 and 2.1 (Illumina), "separate
10- and 15-cycle reactions were run"; MinElute, 15 µl elution buffer at 50 °C → Agilent
DNA 1000 and Nanodrop QC, 10 nM dilution → cluster generation and 27 cycles of
sequencing on the cluster station and 1G analyzer. Three biological replicates.

🟢 (guide; differs from the article as shown on the page) ~10 ng ChIP DNA: end repair with T4 DNA polymerase, Klenow and T4 PNK with
dNTPs in T4 ligase buffer with ATP, 30 min at 20 °C, QIAquick → Klenow fragment
(3′→5′ exo−) with 1 mM dATP, 30 min at 37 °C, MinElute. The guide states the adapters
"have a single 'T' base overhang at their 3' end" → Adapter Oligo Mix (part 1000521)
diluted 1:10, DNA ligase, 15 min at room temperature, MinElute → 2% agarose gel, excise
200 ± 25 bp (also removes unligated adapter), QIAquick gel extraction → Phusion PCR with
PCR primers 1.1 (part 1000537) and 2.1 (part 1000538), 18 cycles of 98/65/72 °C, MinElute
→ Bioanalyzer validation.

🟡 Schmidt et al. 2009 (later, different laboratory) substitutes commercial enzymes but
keeps the Illumina "Genomic Adapter Oligo mix" and "Genomic PCR primer 1.1/2.1", uses a
40-fold adapter dilution, and runs PCR **before** gel size selection (section 2.2.10).
This is a different protocol, not a published discrepancy for Robertson et al.

🟢/🟡 Article versus guide: the article size-selects the ChIP DNA by 12% PAGE
(100–300 bp) *before* end repair and describes no post-ligation gel, and runs 10- and
15-cycle PCRs; the guide gels after ligation (2% agarose, 200 ± 25 bp) and uses 18
cycles. The model follows the article order; the guide supplies only reagent details
the article leaves to Illumina's instructions (T4 PNK, dNTP/ATP buffer, 1 mM dATP,
ligation time). The page shows this as a source discrepancy. Whether the authors also
ran the guide's post-ligation gel, and which PCR product was sequenced, is not stated.

## 5. Oligo bases

🔴 (in the inspected originals) Neither the 2007 Illumina guide, the Robertson article
Methods nor the Robertson supplement prints the Adapter Oligo Mix, PCR primer 1.1/2.1 or sequencing-primer bases.
Schmidt et al. 2009 print reamplification primers, citing a UConn core-facility wiki:

- Primer 1.2: `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`
- Primer 2.2: `5′-CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT`

🟡 The model reconstructs the library from these: primer 1.2 equals P5 followed by
nt 5–33 of the Read 1 sequence used as the adapter's T-overhang strand, and primer 2.2
supplies the PCR primer 2.1 bases from which the phosphorylated adapter oligo is derived
(its reverse complement without the 3′ T junction base). Equality of 1.2/2.2 with
1.1/2.1, and of the reconstructed oligos with part 1000521, is assumed; all these bases
are drawn as inferred on the page.

## 6. Read layout

🟢 One single-end genomic read of {{= chip_seq.FINAL_READ_LENGTH }} cycles: article Methods line 785 ("27 cycles
of sequencing on the Illumina cluster station and 1G analyzer"); lines 804–805 ("a 27-bp
read directionally represents one end of a DNA fragment (SET)"); also the Supplementary
Data 1 README (MOESM306, lines 18–21, "number of 27 bp Illumina reads").
🟡 No sample index is described; the unindexed layout is inferred.

## 7. Open questions

- Whether a post-ligation gel was also run, and which of the 10- and 15-cycle PCR
  products was sequenced (article silent).
- Historical Illumina adapter/oligonucleotide sequence document giving part 1000521 and
  1000537/1000538 bases: not obtained.
- No missing-information claim is accepted; this is a capped single-reader pass.

## Bounded curation round 2 (2026-10-16)

Inspected: Local source bundle MANIFEST.tsv inspected; existing get_sources --slugs chip-seq__10.1038+nmeth1068 attempted:0 originals returned.

No correction made: defining original and historical vendor originals were not available in this bundle. No source claim was re-verified.

Deferred: Defining Robertson2007 article, Schmidt2009 protocol and archived Illumina original11257047/adapter guide; retrieval yielded no source and remains deferred. This capped pass does not establish a complete source review or publication-wide absence.

## Deferred defining-source pass (2026-10-17)

Retrieval: `get_sources --slugs` returned nothing because the catalogue DOI is wrong;
`get_sources --doi 10.1038/nmeth1068` fetched the Springer supplements; `--doi
10.1016/j.ymeth.2009.03.001` fetched the PMC Cloud author-manuscript XML; the Illumina
guide PDF was fetched by hand from its support URL. All three are now in the protocol's
external source bundle with manifest rows. Corrections: DOI, chromatin conditions,
guide reagents and conditions, gel waste including unligated adapter, adapter/primer
bases and read length downgraded to inferred/unverified, constant-restating selftest
check removed. Deferred: article Methods and the historical adapter sequence document.

## Deferred follow-up pass 02 (2026-10-17)

Retrieval routes for the main article: OpenAlex and Unpaywall report it closed (no
repository copy); Europe PMC lists it as subscription-only (no PMCID); nature.com
`/articles/nmeth1068` and `.pdf` return the abstract page with institutional-access
options only. The nature.com page lists a second supplement, MOESM306 ("Supplementary
Data 1, data files for significant peaks"), which `get_sources --slugs` now fetched. Its
README states 27 bp Illumina reads, so the read length is now sourced. Main-article
Methods remain unread (manual URL `https://doi.org/10.1038/nmeth1068`).

## Deferred follow-up pass 03 (2026-10-17)

Source: the hand-placed `nmeth1068.pdf` (17 pp.; printed DOI 10.1038/NMETH1068, title
"Genome-wide profiles of STAT1 DNA association using chromatin immunoprecipitation and
massively parallel sequencing", Robertson … Jones — identity matches). Read lines 1–75,
144–157, 727–729 and Methods 750–815, 849–855. The scraper finds no oligo in the article.
Corrections: the graph now follows the article order (12% PAGE 100–300 bp size selection
of the ChIP DNA before end repair; phenol–chloroform/ethanol cleanups after end repair
and A-tailing; 1/10 adapter mix; QIAquick after ligation; no post-ligation agarose gel;
Phusion with primers 1.1/2.1 in separate 10- and 15-cycle reactions; MinElute; Agilent
QC; 27 sequencing cycles). The article/guide difference is a graph-owned discrepancy on
the page. Single-end reading is now sourced; no index remains inferred. Deferred: the
historical Illumina sequence document for parts 1000521/1000537/1000538, the cited
predecessor ChIP method (ref. 30) beyond the supplement. (The catalogue and directory
already use the corrected DOI.)
