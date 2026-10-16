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

🟡 LAMP does not yield a unique linear intermediate. It yields stem-looped,
branched concatemers containing repeated and inverted target segments. The final
PCR selects compatible units from that population. Consequently, the page preserves
the exact primer landing sequences and embedded barcode but does not assert one
universal order or length for the intervening LAMP-derived sequence.

## Sequencing

🟢 The MiSeq Nano and iSeq runs use Read 1 100 cycles, Index 1 8 cycles and Index 2
8 cycles; Read 2 is not acquired. Read 1 is retained after extracting the embedded
FIP barcode and confirming non-primer viral sequence. The two index reads identify
the post-pooling PCR group.
