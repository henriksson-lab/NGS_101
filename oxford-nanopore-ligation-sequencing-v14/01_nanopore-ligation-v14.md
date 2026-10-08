# Oxford Nanopore ligation sequencing V14 (SQK-LSK114)

## 1. What it is

Oxford Nanopore's standard ligation workflow repairs and dA-tails long double-stranded
DNA, ligates motor-loaded sequencing adapters to both ends, and feeds one strand through a
nanopore without PCR or a sequencing primer.

Commercial defining source: Oxford Nanopore, “Ligation sequencing gDNA V14 —
SQK-LSK114,” **GDE_9161_v114_revAC_24Sep2025**.

## 2. Sources read

- 🟢 Oxford Nanopore SQK-LSK114 genomic-DNA protocol, GDE_9161 v114 Rev AC.
- 🟢 Oxford Nanopore Chemistry Technical Document, CHTD_500 v1 Rev AV, 24 July 2026.
- 🟢 Oxford Nanopore “Nanopore sensing” technical explanation.

## 3. End preparation

🟢 NEBNext FFPE DNA Repair v2 and Ultra II End Prep act at 20 °C for 5 min and 65 °C
for 5 min. The chemistry document states that repair/end prep produces uniform 5′
phosphates and single 3′ dA overhangs. `repair_and_dA_tail()` records these as properties
of both insert ends.

## 4. Ligation Adapter

🟢 The current chemistry document prints the two disclosed LA oligos:

```text
top     CCTGTACTTCGTTCAGTTACGTATTGCT
bottom  GCAATACGTAACTGAACGAAGTACAGG
```

🟡 The bottom is the reverse complement of the first 27 top-strand bases. The result is a
27-bp duplex with one top-strand 3′ dT overhang. The shared adapter constructor verifies
that relationship and refuses an incompatible dA/dT ligation.

🟢 Ligation uses LA, Ligation Buffer and Salt-T4 DNA ligase for 10 min at room
temperature. Optimal output requires adapted DNA ends on both sides. Long Fragment Buffer
enriches fragments above approximately 3 kb; Short Fragment Buffer retains all sizes.

## 5. Motor, tether and sequencing

🟢 LA is supplied with a motor protein. The motor associates with a nanopore, unzips the
duplex and controls translocation of one intact strand. This is direct sequencing: there
is no conventional sequencing primer.

🟢 Flow Cell Tether is mixed into the flow-cell priming solution and concentrates library
at the membrane. 🔴 The vendor does not establish a covalent tether-to-library bond.

🔴 The complete Y-shaped leader sequence, current motor identity, attachment chemistry,
full tether structure and definitive strand polarity are not disclosed. The page therefore
shows a named proprietary leader/motor feature and does not invent bases or a covalent
tether segment.
