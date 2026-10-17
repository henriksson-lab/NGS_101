# TUC-seq — original source boundary

[The original Riml 2017 publication](https://doi.org/10.1002/anie.201707465) describes osmium-mediated transformation of 4-thiouridine into **native cytidine**. This is distinct from the substituted cytidine analogue C* produced by TimeLapse chemistry; the two reactions do not share a molecular-product declaration merely because both yield apparent U-to-C sequence changes.

## Sources read (deferred follow-up 03, 2026-10-17)

- 🟢 **Supporting Information** of the defining article, `anie201707465-sup-0001-misc_information.pdf`
  (hand-downloaded from Wiley on the institutional network; internal name
  `anie_201707465_sm_miscellaneous_information.pdf`; printed title "Osmium-Mediated
  Transformation of 4-Thiouridine to Cytidine as Key To Study RNA Dynamics by Sequencing",
  Riml, Amort, Rieder, Gasser, Lusser, Micura — matches PubMed 28817234). Text twin lines:
  metabolic labelling 674–680, RNA isolation 681–684, "4sU-to-C conversion. See
  Experimental Section in the main text" 694, tRNA-Val RT/cloning 703–714, amplicon
  sequencing 716–724, data analysis 726–748, Supporting Fig. 1 legend 801–810, Supporting
  Table 4 27434–27490. Supporting Tables 1–3 are mutation-frequency data.
- 🟢 PubMed abstract (`pubmed_28817234_abstract.txt`).
- 🔴 The **main article** (Angew. Chem. Int. Ed. 2017, 56, 13479–13483) is still not on disk;
  its Experimental Section holds the conversion protocol the SI defers to.

## What the SI establishes

- 🟢 HEK293 cells (293-HA-StrepIII-eGFP, Dox-inducible eGFP), 5 × 10^5 per well, labelled
  with 0.05 or 0.1 mM 4sU for 1 h; doxycycline (1 µg/ml) added after 30 min.
- 🟢 Total RNA: TRIzol, DNase I (NEB), phenol–chloroform, 0.3 M NaAc precipitation.
- 🟢 Conversion: synthetic 17-nt and 42-nt 4sU RNAs "treated with OsO4 (450 µM) in
  ammonium chloride buffer (180 mM)" (Supporting Fig. 1 legend); cellular RNA "treated with
  OsO4 as described above for 3 h at room temperature" (amplicon section). 🔴 pH,
  the OsO4 concentration used on cellular RNA and the purification after conversion are in
  the unread main text.
- 🟢 Amplicon route: purification, GoScript RT (Promega) with random hexamers,
  amplification of selected targets "with specific primers containing barcode overhangs
  using standard PCR conditions", 2% agarose gel, gel purification, equimolar pooling;
  "Library preparation from the amplicon pool and sequencing using the Illumina HighSeq
  platform was performed by GATC Biotech". Reads demultiplexed by sample-specific forward
  and reverse barcodes (flexbar); U-to-C positions called within 146 bases of amplicon ends.
- 🟢 Supporting Table 4 "Amplicon Sequencing" prints 5′→3′ gene-specific cores for CcnA2,
  CcnB1, CcnD1, CcnE1, p21/Cdkn1A, eGFP and PCNA (forward and reverse; shown on the page
  from the model). The barcode overhangs are not printed 🔴 (not confirmed absent).
- 🟢 **Published cross-reference discrepancy**: the amplicon section says the primers
  "are listed in Supporting Table S1", but Supporting Table 1 is "T-to-C mutation
  frequencies in amplicon sequences" and the primers are in Supporting Table 4. The model
  uses Table 4; shown on the page.
- 🟢 tRNA-Val branch: GoScript RT with a tRNA-Val primer, a stem-loop primer annealed to
  the cDNA 5′ end (65 °C 5 min, 25 °C 10 min), Taq second strand, PCR with stem-loop and
  3′-end primers (Table 4), pGEM-T cloning and sequencing of 40 clones. A non-NGS endpoint,
  not drawn.

Model changes (2026-10-17): the RNA-handling magic step became a sourced TRIzol/DNase
extraction; conversion now carries the SI conditions (partial, with the unassigned parts
named); GoScript random-hexamer RT, barcoded gene-specific PCR (Table 4 cores, overhangs as
role tokens, primer landing unavailable), gel purification/equimolar pooling with
gel-excluded waste were added; the GATC Biotech library preparation remains one warned
boundary. The earlier statement that the SI routes return HTTP 403 is superseded by the
hand download.

The local model distinguishes RNA carrying a pulse-incorporated label from RNA without a pulse-installed label. This classifies incorporation, not the age of every base or the entire molecule. A representative converted site is native C; a possible unconverted labelled class remains explicit, without claiming conversion yields or a physical separation.

## Earlier source-route history

The actual Micura institutional publication entry links the original Wiley abstract and a named 2018 commentary, without an author PDF. In round 21 the NCBI E-utilities abstract (PMID 28817234) was read. The earlier model took the "osmium tetroxide in ammonium chloride buffer" pairing from the same authors' TUC-seq DUAL main text (PMC7186826 line 140); the original SI now states it directly (Supporting Fig. 1), so DUAL is no longer the source for that reagent.

## Later positive comparisons

The same authors' [TUC-seq DUAL publication](https://doi.org/10.1002/anie.201916272) and its original supplement describe a later implementation. Supplement section 1.12 explicitly uses OsO4-hydrazine-treated RNA, GoScript with random hexamers, gene-specific PCR with barcode overhangs, gel purification and pooled amplicons sent to Eurofins for paired-end 150 sequencing. Supporting Table 2 prints target-specific primer cores without the barcodes. The fluorescent tRNA primer belongs to a different primer-extension experiment. These are positive later facts, not a reconstruction of the original 2017 library.

The [later detailed TUC-seq chapter](https://doi.org/10.1007/978-1-4939-9822-7_10) and its [identified 2024 correction](https://doi.org/10.1007/978-1-4939-9822-7_27) are also separate sources. The 2024 correction PDF (`978-1-4939-9822-7_27.pdf`, printed DOI verified) was read on 2026-10-17: it corrects only an author name ("Jeffrey M. Bhasin"), no chemistry. The revised chapter PDF (`978-1-4939-9822-7_10.pdf`) is now on disk but its procedure was not read in this capped pass.

## Remaining source gate

No missing-information claim is accepted. The original main article (Experimental Section), the barcode overhang sequences, actual later chapter/correction contents and named chemistry-relevant linked sources remain explicit in the review. Library selection and historical reagent identities are model questions until supported by an actual original; a later common readout cannot supply an unobserved reaction route.
