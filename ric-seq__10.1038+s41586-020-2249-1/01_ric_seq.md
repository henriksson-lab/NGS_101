# RIC-seq — research notes

Defining source: Cai *et al.*, [10.1038/s41586-020-2249-1](https://doi.org/10.1038/s41586-020-2249-1).

## Evidence used

- 🟢 RIC-seq globally profiles RNA–RNA spatial interactions in situ.
- 🟢 Crosslinking preserves RNA–protein assemblies, fragmented RNA ends are labeled and spatially proximal RNAs are ligated.
- 🟢 Crosslinks are reversed and chimeric molecules are enriched and converted to a paired-end sequencing library.

## Steps after crosslink reversal (Extended Data Fig. 1)

🟢 The publisher-page Extended Data Fig. 1 legend (accessed 2026-10-17; Methods behind the paywall) states: Bioanalyzer traces of "total RNA after rRNA depletion", "RNA after fragmentation" and "RNA after C1 beads pull-down"; libraries "purified from the gel (the bracketed region)"; "−pCp denotes samples without pCp–biotin labelling". The graph adds rRNA depletion, fragmentation of recovered RNA, C1 streptavidin pull-down (unbound fraction discarded) and library gel purification in that order. 🟢 Reporting Summary (MOESM2) line 36: reads acquired on HiSeq X Ten. In-situ fragmentation, RT primer and the adapter-ligation workflow are now taken from GEO (next section). The labelling enzyme and proximity-ligation conditions are now taken from the defining Methods (next section).

## Defining Methods (Cai 2020, read 2026-10-17)

🟢 Source: publisher PDF `s41586-020-2249-1.pdf` (printed doi:10.1038/s41586-020-2249-1, hand-downloaded by the maintainer; MANIFEST row), Methods "Construction of RIC-seq libraries", `.pdf.txt` lines 2067–2195, and "RIC-seq data processing", lines 2372–2380.

- 🟢 Fixation in 1% formaldehyde, 10 min, "quenched by 0.125 M glycine" (the GEO text gives the 2.5 M stock added at 1/20 volume).
- 🟢 After MNase: 10 U FastAP, 37 °C 15 min; washes. "To deposit pCp–biotin": 4 µl of 1 mM pCp–biotin (Thermo Fisher 20160), 100 U T4 RNA ligase, RNase inhibitor and 50 µl of 30% PEG 20000 in 100 µl, overnight at 16 °C.
- 🟢 A second 10 U FastAP treatment (15 min), then T4 PNK (EK0032) with ATP, 37 °C 45 min. 🟡 The Methods give no purpose for these two steps. The 2021 Nature Protocols overview (lines 355–361) states that the second phosphatase removes the 3′ phosphate of the added Cp–biotin and that PNK phosphorylates RNA 5′ ends. The graph uses those purposes.
- 🟢 In-situ proximity ligation: 10 U/µl T4 RNA ligase (EL0021, 10 µl), RNase inhibitor and BSA in a 200-µl mix, overnight at 16 °C. This resolves the earlier inference that the proximity ligase is T4 RNA ligase.
- 🟢 Library: after A-tailing, "Illumina Y-shaped adaptors (2 pmol)" are ligated with Quick T4 DNA ligase. PCR uses "Illumina PE1.0" plus an index primer with Pfx and USER. 200–450-bp products are purified with the MinElute Gel Extraction Kit, and sequencing is on Illumina HiSeq X Ten. 🔴 (not printed in the Methods; not confirmed absent) the Y-shaped adaptor bases and ends.
- 🟢 Libraries were "paired-end-sequenced at Novogene"; adapters were trimmed with Trimmomatic 0.36. 🔴 (not printed in the inspected Methods; not confirmed absent) read lengths.
- 🟢 The Methods state that dUTP marking with Illumina adaptor ligation "performed better than the sequential ligation of RNA adaptors". The model therefore draws only the dUTP route.

The graph now replaces the earlier pCp–biotin MAGIC boundary with explicit FastAP → pCp–biotin ligation → FastAP → PNK → proximity-ligation steps. The pCp ligation carries an `unavailable` end geometry, because the MNase-generated RNA ends have no defined bases.

## Reaction text from GEO GSE127188 (author-deposited sample protocols)

🟢 Source: GEO GSE127188 sample records (SOFT text, `targ=gsm`, fetched 2026-10-17; cached as `geo_GSE127188_acc-targ-gsm_brief.soft.txt`), GSM3629915 `treatment_protocol` (line 13) and `extract_protocol` (lines 16–17); identical text repeats for the other RIC-seq samples. The series links PubMed 32499643 (the defining article). This is the authors' own deposited protocol text, not the Methods; it was used because the Methods are paywalled.

- 🟢 Fixation: 1% formaldehyde in PBS, 10 min at room temperature, quenched with 1/20 volume 2.5 M glycine.
- 🟢 Permeabilization buffer (10 mM Tris-HCl pH 7.5, 10 mM NaCl, 0.5% NP-40, 0.3% Triton X-100, 0.1% Tween 20, protease inhibitors, SUPERase In), 15 min on ice; PNK-buffer washes.
- 🟢 In-situ fragmentation: 1× MN mixture with 0.03 U/µl micrococcal nuclease, 10 min at 37 °C; stopped with PNK + EGTA washes.
- 🟢 "After pCp-biotin labeling and proximity ligation": proteinase K (37 °C 60 min, 56 °C 15 min) and TRIzol LS extraction. 🔴 (not in GEO text; Methods unread) the pCp–biotin attaching enzyme and proximity-ligation conditions. 🟡 The proximity ligase is probably T4 RNA ligase: the next step removes "RNA-DNA ligation byproducts mediated by T4 RNA ligase".
- 🟢 RQ1 DNase I (37 °C, 20 min), acid phenol:chloroform and ethanol precipitation.
- 🟢 rRNA depletion is optional (Ribo-off, Vazyme N406-01). GSE127188 contains both `rRNA_depleted` and `Total` libraries.
- 🟢 RNA fragmentation: 94 °C for 5 min in 1× first-strand buffer (15 mM MgCl2 in the 5× buffer), then ice.
- 🟢 Enrichment on MyOne Streptavidin C1 beads in TWB, 30 min at room temperature; three elutions at 95 °C in PK buffer, then acid phenol:chloroform and precipitation.
- 🟢 First strand: N6 primer and SuperScript II (25 °C 10 min, 42 °C 40 min, 70 °C 15 min). Second strand: RNase H, E. coli DNA polymerase I and dNTPs with 80% of dTTP replaced by dUTP, 16 °C for 2 h.
- 🟢 1.8× AMPure; end repair with T4 DNA polymerase, Klenow and T4 PNK; A-tailing with Klenow exo−; ligation of "adaptors" with Quick T4 DNA ligase (20 °C, 15 min); two AMPure cleanups.
- 🟢 PCR mix: Illumina PE1.0, an index primer, Platinum Pfx and USER enzyme. USER digests the dUTP strand (37 °C, 15 min), then 94 °C for 2 min and 11–13 cycles. Products of 200–450 bp are gel-purified. Sequencing is Illumina paired-end on HiSeq X Ten (line 35).
- 🔴 (unread, not confirmed absent) The GEO text does not give the 2020 adapter oligo sequences, their ligating ends or the read lengths.

The graph now models each of these steps. Adapter ligation uses a structured `unavailable` geometry. pCp–biotin labelling now follows the defining Methods (above). The 2021 PEI adapters (table below) fit the A-tail/T-overhang workflow but are not assumed for 2020.

## Protocol boundary

RIC-seq captures spatial proximity including protein-mediated contacts; it is not represented as a psoralen-selected base-paired duplex assay.

## Oligo audit

The [defining Nature article](https://doi.org/10.1038/s41586-020-2249-1) is the authority for this catalogue entry. Its fetched **Supplementary Table 9**, `MOESM12_ESM.xlsx.txt`, lines 248–262, directly lists the library PCR primers: 🟢 PE1.0 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` and twelve P7 index primers. Each P7 primer has `CAAGCAGAAGACGGCATACGAGAT` followed by a source-listed sample index and `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`. Splitting each printed primer at the shared P7 and `GTGACTGGAG…` arms gives the twelve index windows `AACGTGAT`, `ACAGGAGC`, `ACCGAGCT`, `ACCTTGAA`, `ACACATCG`, `ACGCCTAA`, `ACTGGTCA` (8 nt, primers 1–7) and `CGAAAC`, `GATCTG`, `TTGACT`, `GGAACT`, `TGACAT` (6 nt, primers 8–12). (Model correction 2026-10-17: an earlier version of this note listed primers 9–11 as 7-nt `GATCTGG`, `TTGACTG`, `GGAACTG` by including the first base of the downstream arm.) These have unequal lengths as printed, so they must not be flattened to a six-base placeholder without selecting a concrete primer.

The authors' [2021 Nature Protocols RIC-seq protocol](https://doi.org/10.1038/s41596-021-00524-2), Cao *et al.*, was read on 2026-10-17 from the hand-downloaded publisher PDF `s41596-021-00524-2.pdf` (printed DOI checked). Its Table 1, "Required oligos for RIC-seq library construction" (`.pdf.txt` lines 510–534), prints the oligos below. The authors' later [RIC-seq protocol for infected cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC11245969/) repeats the same four sequences. They are 🟢 for the 2021 protocol. Step 174 there also specifies HiSeq X Ten "paired-end mode" without a read length. The 2020 Methods say only "Illumina Y-shaped adaptors". The 2020 Supplementary Table 9 index primers carry 8- or 6-nt windows, unlike the 2021 6-N window, so these adapters are not assumed for the 2020 libraries.

| Reagent in the 2021 protocol (Table 1 name) | Sequence as printed, 5′→3′ |
|---|---|
| PEI adapter oligo A | `/5Phos/GATCGGAAGAGCACACGTCT` |
| PEI adapter oligo B | `ACACTCTTTCCCTACACGACGCTCTTCCGATCT` |
| Illumina PE 1.0 primer | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` |
| Index primer | `CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` |

Here `/5Phos/` is a 5′ phosphate and the six `N` bases are the index window. The later protocol also specifies an N6 random primer for first-strand synthesis, pCp-biotin for RNA-end labeling, dUTP in second-strand synthesis followed by USER digestion, and an A-tailed DNA intermediate before PEI adapter ligation. Its adapter and primer table is in the *Oligonucleotides* section; molecular steps are in *Adapter preparation* and *Strand-specific library construction*.

The 2020 Methods and GEO text describe the same N6 / dUTP / A-tail / USER workflow, and the schematic models these steps. 🟡 Adapter A and the adapter-B 3′ end pair as a Y adapter with a single 3′-T overhang on adapter B. Read from the printed sequences, that geometry fits the A-tailed cDNA. The exact 2020 adapter oligos and PCR arms on the final construct are still not drawn.

## Review trail

2026-10-17 bounded pass (round breadth-first-16): publisher page abstract and Extended Data Fig. 1 legend; Supplementary Table 9 (MOESM12) lines 250–262 re-verified by computation; Reporting Summary (MOESM2) line 36. Defining Methods (paywalled PDF) and Cao et al. 2021 Nature Protocols remain unread.

2026-10-17 deferred follow-up (round deferred-followup-02): no open copy of the defining article or of Cao 2021 found (OpenAlex/Unpaywall closed, no PMC/Europe PMC full text, no preprint found in a Europe PMC preprint search). Read the GEO GSE127188 sample protocols, GSM3629915 lines 13 and 16–17, and modelled the steps above.

2026-10-17 deferred follow-up (round deferred-followup-03): the maintainer hand-placed the defining Nature PDF and the Cao 2021 Nature Protocols PDF; printed DOIs, titles and authors were checked. Read the Methods library construction (lines 2067–2195) and data processing (2372–2380), plus Nat Protoc overview lines 355–365, Table 1 (lines 510–534), adapter preparation (1547–1580) and step 174. The pCp–biotin boundary was replaced by explicit steps, and reagent conditions were updated. The final library stays role-level, and read lengths remain unstated in the inspected texts.
