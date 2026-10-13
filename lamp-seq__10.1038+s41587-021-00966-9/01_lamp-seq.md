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

🟢 Saturated reactions are pooled and diluted 1:40,000 before a one-step 25-cycle
PCR. The dilution is part of the molecular protocol: insufficient dilution allowed
new barcoded primers to relabel products from other wells.

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
