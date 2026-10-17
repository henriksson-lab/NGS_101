# sn-m3C-seq

## 1. What it is

A single-nucleus multi-omic protocol that reads chromosome contacts and cytosine methylation from the same bisulfite-converted molecules.

## 2. Sources and evidence

- 🟢 Lee et al., *Nature Methods* 2019, DOI [10.1038/s41592-019-0547-z](https://doi.org/10.1038/s41592-019-0547-z); full text [PMC6765423](https://pmc.ncbi.nlm.nih.gov/articles/PMC6765423/).
- 🟢 Proximity-ligation products are carried into bisulfite conversion without a biotin-junction pull-down.
- 🟢 Library construction follows the snmC-seq2 strategy: random priming with a
  partial Read 1 tail, inline index and 3′ H9 arm, followed by adapter tagging and
  PCR. The inherited oligo is printed in Luo et al. Supplementary Methods under
  “Random primers” (for example, P5L_AD002_H); the P5/i5 and P7/i7 arms are added
  by the later indexed PCR primers. The sn-m3C article identifies this predecessor
  procedure in Methods under “Bulk and single-cell m3C-seq library preparation.”
- 🟢 The cited snmC-seq2 source (Luo et al., DOI
  [10.1038/s41467-018-06355-2](https://doi.org/10.1038/s41467-018-06355-2))
  prints eight random-primer designs in its Supplementary Methods and 768 paired
  P5/i5 and P7/i7 amplification-primer pairs in Supplementary Table 2. Those are
  available in the separate [snmC-seq2 note](../snmc-seq2__10.1038+s41467-018-06355-2/01_snmc-seq2.md).
- 🟡 The page draws a representative inherited snmC-seq2 endpoint using the
  source-printed P5L_AD002_H inline barcode and dual-index PCR pair 1. Per-library
  member assignments remain unresolved from the original source chain inspected here. The representative must not be read as that sample's
  actual barcode or index identity. 🔴 The Adaptase-derived bases
  remain unresolved here; converted contact insert bases are sample-dependent.
- 🟢 The authors trim the first 25 bases and final 3 bases of each read to remove
  random-primer sequence and the low-complexity Adaptase tail. These are read
  processing lengths; the drawn random H9 and schematic Adaptase span must not
  be read as the full lengths of those trimmed read regions.

## 3. Distinguishing chemistry

Bisulfite conversion happens after the contact molecule has been created. Consequently each paired read supports both a ligation contact and local methylation inference, unlike a conventional Hi-C read pair.

## 4. Bounded source pass — 2026-10-17

Read defining published article archived manuscript PMC6765423 JATS68,
80 and Supplemental Methods136–152: Hi-C versus in situ 3C, FANS, methylome
library preparation and Data Processing. Manifest identifies the article's
published DOI and PMC original. The cached supplement1 DOCX was inspected
at its initial AOP/Issue capsule only; it does not provide the detailed
chemistry in this pass. Ran the bundle oligo scraper to completion, including
the large extracted workbooks; no sequence hit was returned. This no-hit
result is not evidence of absence, and workbook bodies were not read as
chemistry originals during this pass.

🟢 The selected in situ 3C route uses DpnII digestion, 62 C for 10 min heat
inactivation and T4 DNA ligation, omitting the biotin fill-in used for Hi-C.
The model now represents heat inactivation and names T4 ligase. Intact
singlet nuclei are selected by FANS/Hoechst/DNA-content gates and deposited
into 384-well plates. Proteinase treatment precedes conversion (article68)
and is now a separate graph operation, with detailed conditions unresolved.
Other nuclei remain outside this lineage with unresolved fate.

🟢 Article146 states PE150 sequencing for both cultured-cell and cortical
samples, with instrument choices varying by sample. The graph and public
readout now declare two 150-base acquired reads. Trimming 25 bases upstream
and 3 bases downstream from both reads occurs computationally after acquisition;
it is not experimental removal of those molecules.

🟡 The old graph retained its unconverted contact sequence unchanged after
bisulfite conversion, copying and Adaptase treatment. These sequence-changing
intermediates now remain unrepresented rather than claiming that unchanged
molecule. The inherited snmC-seq2 endpoint is retained provisionally, not newly
verified here. Its generic schematic index widths do not establish real
index lengths or acquisition. The decoupling inventory already flags the
shared endpoint for comparison with its owner's source; that follow-up is
not an architecture migration in this capped pass.

Deferred: full linked Luo 2018 snmC-seq2 protocol/primer originals; original
workbooks, especially supplement7's sample metadata; predecessor Hi-C/vendor
originals and Adaptase evidence; full proteinase/conversion/cleanup/priming
states and study-specific primer/index assignments. Current missing-information
claims remain blocked. This is one bounded defining-article pass, not a full
chemistry validation or accepted absence review.
