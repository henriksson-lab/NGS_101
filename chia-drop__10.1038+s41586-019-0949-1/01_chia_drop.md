# ChIA-Drop — research notes

Zheng et al., *Nature* (2019),
[doi:10.1038/s41586-019-0949-1](https://doi.org/10.1038/s41586-019-0949-1).

## Sources read

- 🟢 Defining article and Figure 1 (PMCID: PMC7001875).
- 🟢 [Nature Supplementary Methods](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-019-0949-1/MediaObjects/41586_2019_949_MOESM1_ESM.pdf), "ChIA-Drop and RNAPII ChIA-Drop library construction and sequencing" (p. 6).
- 🟢 Supplementary Table 2 prints many **observed** 16-nt GEMcodes alongside called complexes; these observations are not the complete manufactured gel-bead barcode whitelist.
- 🟢 Primary read metadata in GEO series GSE109355, including GSM3347523 and
  GSM3347525.

## Protocol boundary

ChIA-Drop is a ligation-free chromatin-conformation assay. It must not be drawn as Hi-C:
the remote genomic fragments never become one covalent contact molecule.

1. 🟢 About ten million cells are fixed with 1% formaldehyde for 10 minutes, quenched
   with 0.125 M glycine and lysed while crosslinked complexes remain intact.
2. 🟢 Permeabilized nuclei are restriction-digested in situ. The study compared the
   roughly 300-bp MboI products with roughly 3-kb HindIII products and approximately
   6-kb sonicated chromatin; longer templates produced more useful mapped reads.
3. 🟢 Sonication opens the nuclei and releases fragmented but still crosslinked chromatin.
   There is no end fill-in, biotin marking or proximity ligation, and the DNA is not
   purified before droplet partitioning.
4. 🟢 The chromatin sample is adjusted to 0.5 ng/µl with 50 µg/ml BSA; 0.5 ng is loaded
   with Chromium Genome v2 Library and Gel Bead Kit PN-120258 into the Chromium Genome
   Chip PN-120257.
5. 🟢 Each GEM contains a gel bead bearing many oligos with P5, a Read 1 priming site, a
   16-nt GEM barcode and a random 8-mer. Overnight 30°C isothermal linear amplification
   random-primes accessible stretches of several DNA fragments in the same complex.
6. 🟢 Products from one GEM are separate molecules bearing the same barcode. After GEM
   breakup they are pooled, end-repaired, A-tailed and adapter-ligated using the **Chromium
   Multiplex Kit PN-120262**; Sample Index PCR and 600-bp size selection finish the library.

## Barcode and read meaning

- 🟢 The 16-nt GemCode is a **molecular-complex/partition barcode**, not a cell barcode:
  it groups all fragment copies made in one droplet.
- 🟢 Figure 1 calls the adjacent sequence a random 8-mer. GEO calls those same eight
  bases a UMI. These statements are compatible: the random priming sequence is retained
  on the product and can distinguish copies.
- 🟢 The defining libraries were sequenced 2 × 150 bp. Read 1 starts with the 16-bp
  barcode and 8-bp UMI before genomic sequence; Read 2 is biological/genomic sequence.
- 🟢 ChIA-DropBox groups reads by GemCode, merges overlapping reads into fragments and
  groups fragments into putative complexes. A barcode records a GEM, not an infallible
  one-complex measurement: collisions can put unrelated complexes in one droplet.

## Sequence boundary

🟡 The exact Chromium Genome v2 bead oligo, manufactured barcode whitelist,
PN-120262 adapter and Sample Index PCR primers remain unresolved in this bounded
review. Historical vendor originals are deferred. This is not an accepted claim
that these sequences are absent from all relevant sources. The compatible Illumina
arms and eight-base sample-index window in the schematic are inferred; acquired
index cycles have not been established here. Published FASTQ may help characterize
read-visible sequence, but does not alone establish ordered oligo modifications.

## Bounded source check — 2026-10-17

Read the defining published article's preparation/analysis overview (PMC7001875
JATS text lines 64–66 and Figure 1 caption line 491) and Supplementary Methods
pp. 4–6 (text 137–226), plus Supplementary Note 1 pp. 15–17 (920–1002).
The oligo scraper was run across the bundle. Its example GEMcode at Methods
lines 952–958 was verified there; observed barcodes are not ordered oligo sequences.
Other output-table hits were not promoted to reagent facts.

The graph selects the published HindIII branch followed by sonication, with
MboI and sonication-only preparations retained as alternatives in the caveat.
Removed illustrative 1/2/3-kb per-fragment assignments: no exact molecular lengths
are assigned to those symbolic example molecules. Added explicit pooled DNA
recovery, end repair, A-tailing, adapter ligation, Sample Index PCR and post-PCR
selection around 600 bp. Exact commercial intermediate ends remain unresolved;
no adapter end or junction is invented. Library-adapter ligation does not join
remote contact fragments.

Corrected the former partition waste branch: colliding complexes can be copied
and sequenced. Their inter-chromosomal reads are subsequently split into
intra-chromosomal groups; singleton removal is also computational (article line
66; Supplementary Note pp. 15–17). Neither is a physical separation of GEMs.
The source acquired 2 × 150 bp but analyzed only R1; R2 was reserved for future
analysis (Methods lines 963–968). The index-cycle count remains open.

Deferred: historical PN-120258/PN-120262 vendor originals and cited linked
predecessor methods, exact bead/adapter/PCR oligos and modifications, full
Figure 1 rendered geometry, other supplements and complete per-product recovery
mechanisms. Existing GEO-based UMI assignments above are retained as earlier
provenance, not newly reverified: the generic GSM3347523 endpoint is a manual
challenge entry in the local manifest. This capped check is not a full protocol
review or confirmation of missing information.
