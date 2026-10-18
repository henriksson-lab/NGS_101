# eccDNA Circle-Seq — research notes

## Sources read

- 🟢 Møller et al. 2015, DOI [10.1073/pnas.1508825112](https://doi.org/10.1073/pnas.1508825112), PMC4475933, full text and Methods.
- 2026-10-17 bounded pass: PMC article page text (Results "Circle-Seq, A Method for
  Genome-Scale Detection of eccDNA"; Methods "Purification and enrichment", "Digestion of
  remaining linear DNA", "Amplification", "PCR", "Sequencing", "Mapping of Reads").
  Not read: SI Appendix PDF and Datasets S1–S3 (PMC download gate; URLs in the
  source manifest), Plasmid Mini AX, REPLI-g Mini and Apollo 324 vendor documents, and
  the cited plasmid-purification predecessors (refs 26, 27).
- 2026-10-17 follow-up retrieval (round deferred-followup-02), SI still not obtained:
  PMC `bin/pnas.201508825SI.pdf` returns the JavaScript "Preparing to download" gate;
  PMC Cloud has no `PMC4475933` prefix; Europe PMC `supplementaryFiles` is empty,
  `fullTextXML` returns HTTP 500 and its `bin/` route HTTP 403; the PNAS
  `suppl_file/pnas.201508825si.pdf` route is behind Cloudflare (HTTP 403);
  OpenAlex/Unpaywall list only PMC and a University of Copenhagen research-portal
  record with no file. ENA PRJEB6368 run metadata give `library_construction_protocol`
  "eccDNA purification" only (HiSeq 2000, SINGLE). Takara's current PrepX DNA Library
  Protocol-At-A-Glance documents (2018/2019) exist, but the paper names no Apollo kit
  and those revisions postdate the work, so they were not used. Manual URL:
  https://www.pnas.org/doi/suppl/10.1073/pnas.1508825112/suppl_file/pnas.201508825si.pdf

- 2026-10-17 round deferred-followup-03: the maintainer hand-downloaded the published PNAS
  article PDF (`pnas.1508825112.pdf`, printed doi:10.1073/pnas.1508825112, Møller et al.)
  and the SI Appendix (`pnas.1508825112_pnas.201508825si.pdf`, header "Supporting
  Information Møller et al. 10.1073/pnas.1508825112", 10 pages). Read: article Results
  "Circle-Seq, A Method…" (`.txt` lines 280–327), Fig. 1 legend (270–271), Materials and
  Methods "Strains and Plasmids" through "Mapping of Reads" (700–788) and PCR (771–777,
  894–913); SI Appendix in full (Figs. S1–S10 legends, Dataset S1–S3 legends, lines
  1–147). The SI contains figure legends and dataset descriptions only, no additional
  library chemistry; Datasets S1–S3 are described there as eccDNA region, gene-copy and
  repeat analyses (still not downloaded; not chemistry by their legends). The scraper finds
  only the ACT1/qPCR and inverse-PCR validation primers (PDF lines 774, 906–909); none in
  the SI.

## Molecular order

🟢 Yeast cells (10<sup>10</sup>) are disrupted with glass beads in Plasmid Mini AX L1 solution;
the kit's alkaline treatment, protein precipitation and ion-exchange column separate
circular from chromosomal DNA. Before the column, plasmid spike-ins are added to the lysate
at 1:1 (pBR322), 1:50 (pUC19) and 1:2,500 (pUG72) per cell, assuming 8 × 10<sup>9</sup>
genomes; they serve as internal detection-limit standards (Fig. S5) and are now an input
of the purification step. The eluate (6–16 µg) is precipitated into TE. NotI (1 FastDigest unit, 37 °C, 17 h; 80 °C 5 min) cuts
chromosomal DNA to expose more ends, then Plasmid-Safe ATP-dependent DNase (37 °C, ATP
and 16 U enzyme added every 24 h) runs for 130 hours and is heat-inactivated (70 °C,
30 min). One percent (5 µL) of the enriched fraction is amplified with phi29 polymerase
by the REPLI-g Mini protocol (30 °C, 16 h). The product is sonicated to about 300 nt
(Covaris LE220); 300 ng is loaded on an Apollo 324 robot, which adds adapters and 6-base
barcode index labels, and libraries are sequenced single-end on HiSeq 2000.
🟡 Model correction (2026-10-17): the circle is no longer built by the shared ssDNA
circularisation helper, which implied a ligation that never happens; it is declared as a
pre-existing closed circle, one strand drawn. The former representative phi29 sequence/site has since been removed; REPLI-g
priming sites and sequence geometry remain unassigned in the current model.
🟡 The page shows the Apollo step as a warned boundary: the paper does not state the
joining chemistry or whether any library PCR follows. The earlier model's explicit
adapter ligation and library-PCR steps were model insertions and were removed.

🟢 Read length: Results give 141-nt single-end reads; Mapping of Reads says the longest
reads (215 nt; samples S2, Z4, R1) were adapter-trimmed with Cutadapt before the
141-nt reads were mapped. Reads were demultiplexed by index barcode (0 mismatches in the
first lane, 1 in later lanes); how the index was read is not stated.

Published discrepancy (shown on the page): Methods "Sequencing" says samples were
"sequenced up to 100 Mb", whereas Results report 40–100 million reads per population;
100 Mb of 141-nt reads would be about 0.7 million reads. The model sets no depth.

## Boundary and uncertainty

The circles are biological inputs, not library-prep intermediates. The source names
Apollo324 and six-base barcode oligos; actual adapter/index identities, geometry and
molecular placement remain unassigned in this model. The public page now uses an
explicit library boundary instead of a conventional single-index Illumina shell.
This is separate from CRISPR CIRCLE-seq.

🟢 The Methods also print ACT1 PCR/qPCR and inverse-PCR pairs for several confirmed
circles and spiked plasmids. Those oligos test exonuclease removal or validate
individual circles; none enters the sequenced library. No Apollo ordering table or
adapter lot is identified, so raw reads could at most support inferred adapter
sequence and index content, not chemical modification or reagent identity.

## Source-owned validation migration (2026-10-18)

This capped pass read two originals: the defining PMC4475933 article Results and
Methods (Circle-Seq purification/digestion/amplification, Sequencing and Mapping of
Reads), plus the published SI original identity and FigS1/FigS5 captions (linear-DNA
hydrolysis assay and independent plasmid spike-ins). Mandatory bundle primer-scraper
pointers were checked in these originals. Detailed predecessor plasmid methods,
REPLI-g primer chemistry and Apollo324 kit documents remain deferred.

`batch_ngs.truseq_library` and its assumed dA junction/standard arm placement are
removed. Apollo324 supplies source-named adapter reagents and six-base barcode
oligos, represented as a reagent inventory with unknown molecular placement.
The final molecule’s adapter/index arrangement and library PCR status remain
unassigned, rather than choosing a canonical TruSeq endpoint.

The arbitrary16-base native-circle sequence and8-base phi29 primer were removed.
Native eccDNA and independent plasmid standards retain source-qualified circular
roles; REPLI-g priming sequences, repeat count and branching geometry remain
unassigned. Different cellular DNA populations occupy separate Scene rows.
Plasmid-Safe hydrolysis produces an explicitly unresolved molecular product class,
not physically discarded waste. Source-described1%(5µl) sampling is explicit, with
remaining99% disposition unassigned. Subsequent sonication produces fragment roles
rather than retaining an unfragmented RCA drawing.

The source141base single-end reads and215base S2/Z4/R1 reads remain distinct sample
scopes. Six-base barcode oligos do not prove a six-cycle index acquisition. The
published100Mb versus40–100million-read depth discrepancy remains visible on the
page; it does not select a molecular endpoint. This is not a full source audit or
an accepted publication-absence claim.
