# ZipSeq — source notes

Defining paper: Hu et al., *Nature Methods* (2020), [doi:10.1038/s41592-020-0880-2](https://doi.org/10.1038/s41592-020-0880-2).

- 🟢 PMC manuscript, Methods and oligonucleotide list were read.
- 🟢 Antibody- or lipid-linked anchors are protected by NPOM-caged strands. Patterned **365-nm** light uncages a selected region (published Methods, microscopy and on-scope validation; manuscript lines 244–264). The former 405-nm page value was a repository error, not a published discrepancy.
- 🟢 A complementary zipcode with a poly(A) tail hybridizes after each illumination round; rounds can be repeated with different cages/zipcodes. The printed zipcode 1 has 28 A bases.
- 🟢 Dissociated cells undergo 10x single-cell RNA-seq, which captures zipcode and endogenous transcript molecules in the same GEM.
- 🟢 The paper's Methods oligonucleotide list gives the antibody anchor as
  `Amine-C6–ATCGTTTTTTTTTTTTGATGTTGACGG`, the O1 blocking strand as
  `CGATCTGTGGTTGCTAC`, and zipcode 1 as
  `GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGC` followed by 28 A bases.
  The final eight bases before poly(A), `CATGATGC`, are the region barcode;
  the structured feature is spatial identity, distinct from the droplet cell barcode.
- 🟢 The same list gives the caged O1 strand, second-layer O1′ strand,
  second-layer caged O2 strand, zipcodes 2–4, fluorescent readout oligo,
  lipid-anchor adaptors and the alternative design-2 anchor/zipcode oligos.
- 🟢 The caged O1 strand is printed as
  `CGATCTGTGGTTGCTACCCGTCAACATCAAAAAAAAAAAA`; the page now lists these
  bases alongside the anchor and first zipcode.
- 🟢 The rendered bioRxiv v1 Methods p. 17 marks O1 caged-strand thymidines
  **4, 8, 11 and 15** in red; the publisher's JATS transcription drops the
  colour. The preprint's O1 bases agree with the published manuscript.
- 🟢 The publisher's Reporting Summary, supplementary figure 1 methods
  (extracted lines 233–244), states that caged O1 and the anti-CD45
  antibody–anchor conjugate were prehybridized 1:1 at 37 °C for 15 minutes
  before cooling to room temperature and illuminating with 365-nm light.
- 🟡 The published manuscript's library-construction paragraph (lines 280–286)
  says that a partial Read 2 small-RNA primer was added at cDNA amplification,
  short zipcode products were retained from the SPRI supernatant, and that
  CITE-seq primers amplified the zipcode fraction. The CITE-seq supplementary
  protocol (p. 5) gives the primer *classes* and PCR recipe, but this does not
  identify the exact indexed primer or 10x kit version for every ZipSeq run.
- 🟢 The schematic's 10x endpoint is the **v3 tumor example**, rather than a
  universal ZipSeq library. Hu et al. specify v2 for the in-tube, wound-healing,
  and lymph-node experiments and v3 for tumor (published Methods lines 238,
  256–280). Official CG00052 Rev D, Sequencing p. 44, specifies v2
  26/8/0/98 read cycles and the v2 bead has a 10-base UMI; CG000183 Rev A,
  Step 4 and oligo appendix pp. 45 and 53–55, specifies v3 28/8/0/91
  cycles and a 12-base UMI. The page now labels its v3 choice.
- 🟢 McGinnis et al.'s original MULTI-seq author manuscript, Online Methods
  p. 8 (extracted lines 465–494), gives the 5′ lignoceric-acid anchor
  `GTAACGATCCAGCTGTCACTTGGAATTCTCGGGTGCCAAGG` and 3′ palmitic-acid
  co-anchor `AGTGACAGCTGGATCGTTAC`, matching Hu et al.'s lipid-anchor
  list. Its sample barcode has a 30-base poly(A) tail; this does not replace
  the 28-base tail printed for Hu et al.'s zipcode 1.
- 🟢 The synthesis reference cited by MULTI-seq is Weber et al.,
  [doi:10.1021/bm501467h](https://doi.org/10.1021/bm501467h). Its original
  article, Experimental Section (extracted lines 95–136), couples a fatty-acid
  carboxyl group to an oligo's 3′ or 5′ free amine while on solid support,
  making the terminal amide linkage. Its original supplement, p. S3
  (extracted lines 25–41), lists the 20-base co-anchor and C16–C24 saturated
  fatty-acid choices. Weber's construct and synthesis establish the inherited
  attachment chemistry; Hu et al.'s own Methods remains the authority for
  ZipSeq-specific oligo identities.


## Local endpoint reconstruction (2026-10-16)

The shared endpoint recipe had inserted a Nextera Read 1 arm into both
libraries. Hu's represented v3 gene-expression branch instead uses the
TruSeq arm printed in CG000183. Its assembly and 28/8/91 source run now live
in the ZipSeq module.

The zipcode branch has its own source boundary. Hu's Library Construction
Methods says that the additive primer and zipcode amplification follow
CITE-seq. The original CITE-seq supplemental protocol (2017-07-03, p. 6)
prints the 10x PE 1.0 P5 oligo and the RPI1 primer, whose 33-base Small RNA
Read 2 tail is different from ordinary TruSeq Read 2. The reconstructed
zipcode endpoint is derived from Hu's exact Zipcode 1 and this predecessor
PCR architecture, with a six-position RPI-family index window. The actual
indexed RPI member and loaded sequencing-primer oligos for Hu's experiment
remain unverified, so the page marks the extension inferred and uses a
source-declared read table rather than asserting exact I1/R2 binding.

The inherited RPI 3′ footprint spans the native partial Small RNA handle
plus the initial `CA` of the named eight-base Zipcode 1. That computed PCR
footprint does not by itself establish the experiment's sequencing entry
or the number of zipcode bases acquired by its loaded Read 2 primer.
Constructor checks compare the final endpoint to the original predecessor
PCR pair; incompatible Nextera Read 1 and ordinary TruSeq Read 2
substitutions are rejected. A separate local-source check reads the original
CITE-seq protocol text to guard the RPI1 transcription.
