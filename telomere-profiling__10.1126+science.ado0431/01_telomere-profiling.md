# Telomere Profiling

## 1. What it is

Telomere Profiling combines restriction digestion, a phased biotinylated TeloTag, streptavidin enrichment and Oxford Nanopore sequencing to concentrate intact human chromosome ends.

## 2. Sources read

- 🟢 Karimian et al. (2024), *Science*, [doi:10.1126/science.ado0431](https://doi.org/10.1126/science.ado0431): author manuscript PMC12425158 JATS (`science.ado0431_PMC12425158.xml.txt`), Results lines 122–138 and the Supplementary Material listing (lines 182–186).
- 🟢 Published Supplementary Materials (`science.ado0431_sm.pdf`), hand-downloaded by the maintainer on 2026-10-17 from the Science `suppl_file` URL; the printed title, DOI and "Karimian et al." were checked. Read in round deferred-followup-03: Materials and Methods "Annealing of TeloTags", "Telomere Tagging" and "Telomere Enrichment and Nanopore Sequencing" (`.pdf.txt` lines 77–147), basecalling (lines 205–213), fig. S1A (lines 262–281, legend 454–456), fig. S5 legend (lines 977–981) and Table S3 "TeloTag Sequences" (p. 20, lines 2296–2357).
- 🟡 Karimian et al. bioRxiv v2, [doi:10.1101/2023.12.21.572870](https://doi.org/10.1101/2023.12.21.572870): superseded preprint. It was used in earlier rounds only while the published supplement was unavailable, and it no longer supplies any modelled value (see review trail).
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), workflow comparison (not used for values).
- Unread leads, not needed for this pass: the inventors' patent PCT/US2023/073375 (JATS line 207) and a 2026 preprint protocol (doi:10.64898/2026.08.15.745035).

## 3. Construct and reaction order

🟢 Published Supplementary Materials and Methods (pp. 3–4):

1. TeloTag duplex assembly in 100 µL: "5mM of each of the 6 permutations of telomere splint (fig. S1A) and 30 mM of biotinylated adapter" in 1x HiFi Taq DNA Ligase Reaction Buffer (NEB M0647S); 99 °C, then −1 °C/min; diluted 1:100 in 1x Taq buffer, kept at 4 °C (lines 77–83).
2. 40 µg HMW gDNA with 3 µL ClaI (NEB R0197S), 2 h at 37 °C; heat-inactivated at 65 °C for 20 min (lines 86–89).
3. Ligation with 4 µg DNA per 50 µL: 4 µL of "0.3μM duplex TeloTag adapter", 5 µL 10X HiFi Taq ligase buffer and 1 µL HiFi Taq DNA Ligase. 65 °C for 5 min, then 15 cycles of 65 °C for 1 min / 45 °C for 3 min (lines 92–100).
4. Pooling; SPRI select (Beckman B23318), 45 µL : 100 µL, 20 min on a Hula mixer; two 85% ethanol washes; elution in 65 °C 1X rCutSmart for 20 min (lines 106–117).
5. Dynabeads MyOne Streptavidin C1 (65001), 1 µg : 250 ng DNA in kilobaseBINDER binding buffer, room temperature for 20 min. Washes: 2× kilobaseBINDER wash, 1× EB, 1× rCutSmart, 1× EB (lines 118–130).
6. Release in 72 µL 1X rCutSmart with 3 µL PvuI (NEB R3150S), 37 °C for 30 min, then 65 °C for 20 min; expected recovery about 0.1–0.01% of input (lines 131–139).
7. Standard ONT library preparation with "Ligation Sequencing Kit V10 (SQK-LSK110)" on R9.4.1 FLO-MIN106D flow cells; elution in 40 µl EB. Each library is split into three loads of about 16 h with EXP-WSH004 washes; MinKNOW 5.7.5 without live basecalling (lines 140–147). Basecalling uses Guppy v6 super-high-accuracy (line 211).

🟢 Fig. S5 (lines 977–981) also compares AsiSI/PvuI and BamHI/EcoRI digests. These alternative branches are not drawn.

🟢 Table S3 prints 14 TeloTag top strands with `/5Phos/` and `/3BioTEG/` (NB50, NB65, NB68, NB88, NB01, NB02, NB15, NB19, NB20, NB66, NB67, NB69, NB70, NB72); all are now in the model. NB50/65/68/88 carry EcoRI before PvuI; the other ten carry PvuI before EcoRI. Restriction-site segments are labelled from their bases. 🟢 Fig. S1A prints the six NB50 splints, written 3′→5′. Each one, reversed, equals a splint derived from the NB50 strand (reverse complement of the 87-nt tag body plus an 18-nt telomere-repeat phase). The selftest checks the 14 Table S3 rows and the fig. S1A splints against the source text. 🟡 Splints for the other 13 TeloTags are not printed; the model derives them so that the duplex pairs, and draws only NB50.

Published discrepancies (shown on the page):

- Order: the Results summary (JATS line 122) reads "We ligated telomeric ends with a biotinylated oligonucleotide (TeloTag), digested the DNA, and isolated the tagged telomeres". The Supplementary Methods digest with ClaI and heat-inactivate before the tagging ligation. The model follows the Methods, which give conditions.
- Units: "5mM" splints / "30 mM" adapter, then a 1:100 dilution used as "0.3μM". The printed units are inconsistent and are left unresolved.
- Table label: Methods line 83 says the sequences are "listed in table S2"; Table S2 is telomere lengths, and the sequences are in Table S3 (top strands) and fig. S1A (NB50 splints).

## 4. Final library and sequencing

The motor-loaded Oxford Nanopore adapter, rather than a sequencing primer, controls strand entry through the pore. 🟢 MinION sequencing with Guppy basecalling (published Results lines 122, 126; SM lines 140–147, 211).

🔴 (not printed in the inspected texts; not confirmed absent) The defining workflow does not print the complete historical motor-bearing ONT adapter assembly. Current ONT adapter oligonucleotides are not substituted into this earlier library without evidence that the kit revision is identical.

## 5. Open questions

- The SQK-LSK110 motor-adapter assembly and kit-internal reaction steps (ONT kit document not read).
- The unit inconsistency in TeloTag annealing ("5mM"/"30 mM" vs "0.3μM").
- The splint families of the 13 non-NB50 TeloTags (not printed; derived).
- The fig. S5 alternative digests (AsiSI/PvuI, BamHI/EcoRI) are not drawn.

## 6. Review trail

- Rounds breadth-first-20 and deferred-followup-02: the published supplement was unavailable. Reaction conditions and oligos were taken from bioRxiv v2 (Methods lines 1857–2149, Extended Data Table 2 lines 3275–3446). That table's splint rows for NB68 (one barcode base), NB70/N01 (a forward-written `CTTAAGGCTAGC` window) and the N01 labels (C1, C2, C3, C3, C4, C5) did not equal the derived splints. The preprint also printed SQK-LSK114 with R9.4.1 flow cells, internal enzymes ClaI/AsiSI/PmeI/BamHI and release enzymes PvuI/PacI/EcoRI. These preprint observations are history only. The published supplement prints no splint rows for those barcodes, and its Methods specify ClaI, PvuI and SQK-LSK110.
- Round deferred-followup-03 (2026-10-17): read the published Supplementary Materials (section 2). Provenance moved to the published source; the model gained the eight additional Table S3 TeloTags; enzymes, washes and the ONT kit were corrected; and three published discrepancies are now shown on the page.
