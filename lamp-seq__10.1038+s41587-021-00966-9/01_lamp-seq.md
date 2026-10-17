# LAMP-Seq — source notes

Defining source: Ludwig *et al.*, *Nature Biotechnology* (2021),
[doi:10.1038/s41587-021-00966-9](https://doi.org/10.1038/s41587-021-00966-9).

## Sources read

- 🟢 Main article and Methods, PMC8678193.
- 🟢 Supplementary Data 1, including the primer list, 192 validated Tier-1 FIP
  barcodes, the larger Tier-2 set and pool-index primer table.
- 🟢 Supplementary Notes and application note.

## Molecular interpretation

🟢 Six target primers drive RT-LAMP: F3, B3, FIP, BIP, LF and LB. The 10-nt
sample identifier is physically embedded between the F1c and F2 portions of FIP.
The published open-source implementation uses LNA substitutions in F3 and B3 and
a strand-displacing Bst polymerase at 65 °C for 80 minutes.

🟢 The paper gives two PCR branches. The two-stage branch shown on the page pools
and dilutes LAMP products 1:100,000 before an 18-cycle target PCR, then uses
another 18-cycle pool-index PCR. Adding 5 µl of diluted pool to the 50 µl first
PCR makes the effective dilution 1:1,000,000. The open-source one-step branch
instead dilutes 1:40,000 before a 25-cycle indexed PCR (effective 1:400,000).
The dilution limits barcode swapping caused by relabelling across pooled wells.

🟢 The PCR primers install two additional 10-nt pool identifiers together with P5,
P7 and standard TruSeq sequencing arms. Supplementary Data 1 publishes every
primer; the schematic uses the first validated FIP barcode (`AGCAAACTAC`) and
Pool_fwd_1/Pool_rev_1 as representatives.

🟢 The pool-barcoding primers carry only partial sequencing-primer tails
(Supplementary Data 1 `Pool barcoding primers` row 1; Methods two-step PCR, PMC8678193 text line 277–279):
Pool_fwd_1 is P5 + `TCGTGGAGCG` + `{{= lamp_seq.LAMP_POOL_FWD_TAIL }}` ({{= len(lamp_seq.LAMP_POOL_FWD) }} nt) and Pool_rev_1 is
P7 + `CGCTCAGTTC` + `{{= lamp_seq.LAMP_POOL_REV_TAIL }}` ({{= len(lamp_seq.LAMP_POOL_REV) }} nt). The first PCR supplies the rest of
both arms. PCR-N-rev-11 prints its tail as `TGACTGG…` in both the Methods and
Supplementary Data 1; the 5' G of the pool reverse tail completes canonical `GTGACT`,
so this is a design, not a source error. *Model correction (2026-10-17):* earlier
versions appended full 33/34-nt TruSeq arms to the pool primers (72/68 nt).
The full-arm forms belong to the open-source one-step Ultramer primers (Methods, line 353–355).

🟢 The two-step branch drawn is the commercial-reagent protocol (Methods lines 205–289):
QuickExtract lysis with pX330 carrier DNA, ≥10 min RT then 95 °C 5 min; 30 min
activated-carbon/ion-exchange bead treatment and HCl pH adjustment; NEB LAMP master mix
E1700L plus Bst 3.0; C-F3-LNA, C-B3-LNA, C-BIP, C-LF, C-LB and barcoded C-FIP; 65 °C
80 min; pooling on ice; 1:100,000 dilution; two 18-cycle NEBNext 2× PCRs (98 °C 20 s,
65 °C 20 s, 72 °C 30 s); pooling, 2% E-gel band excision, two silica spin columns,
NanoDrop. The Methods call the target primers PCR-C-fwd/PCR-C-rev; Supplementary
Data 1 prints the identical sequences as PCR-N-fwd-10/PCR-N-rev-11. The commercial
master mix also contains an ACTB control primer set and barcoded ACTB-FIP, not drawn.

🟡 LAMP does not yield a unique linear intermediate. It yields stem-looped,
branched concatemers containing repeated and inverted target segments. The final
PCR selects compatible units from that population. Consequently, the page preserves
the exact primer landing sequences and embedded barcode but does not assert one
universal order or length for the intervening LAMP-derived sequence.

## Sequencing

🟢 The MiSeq Nano and iSeq runs use Read 1 100 cycles, Index 1 8 cycles and Index 2
8 cycles (Methods lines 365–385, open-source section); Read 2 is not acquired. Only
the first 8 nt of each 10-nt pool index are therefore read. 🟡 The commercial branch
names MiSeq or iSeq without its own cycle numbers; the page applies these values. Read 1 is retained after extracting the embedded
FIP barcode and confirming non-primer viral sequence. The two index reads identify
the post-pooling PCR group.
