# Oxford Nanopore ligation sequencing V14 (SQK-LSK114)

## 1. What it is

Oxford Nanopore's standard ligation workflow repairs and dA-tails long double-stranded
DNA, ligates motor-loaded sequencing adapters to both ends, and feeds one strand through a
nanopore without PCR or a sequencing primer.

Commercial defining source: Oxford Nanopore, “Ligation sequencing gDNA V14 —
SQK-LSK114,” **GDE_9161_v114_revAC_24Sep2025**.

## 2. Sources read

- 🟢 Oxford Nanopore SQK-LSK114 genomic-DNA protocol, GDE_9161 v114 Rev AC (text lines
  166/249 revision; 490–611 end-prep and cleanup; 615–800 adapter ligation, cleanup and
  loading amounts; 817–931 priming with Flow Cell Tether; 1207–1212 troubleshooting).
- 🟢 Oxford Nanopore Chemistry Technical Document, CHTD_500 v1 Rev AV, 24 July 2026
  (lines 248, 345, 770, 813, 2067–2081, 2122–2134).
- 🟢 Oxford Nanopore “Nanopore sensing” technical explanation.

## 3. End preparation

🟢 NEBNext FFPE DNA Repair v2 and Ultra II End Prep act at 20 °C for 5 min and 65 °C
for 5 min. The chemistry document states that repair/end prep produces uniform 5′
phosphates and single 3′ dA overhangs. The local explicit `EndPreparedDuplex` declaration records these as properties
of both insert ends.

🟢 Input is 1 µg (or 100–200 fmol) DNA in 47 µl, with optional DNA Control Sample
(a 3.6 kb lambda amplicon). After end prep, 60 µl AMPure XP beads bind the DNA; two
200 µl 80% ethanol washes and a 61 µl water elution follow.

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
temperature (60 µl DNA, 5 µl LA, 25 µl LNB, 10 µl Salt-T4). The chemistry document says
molecules with adapters (and tethers) on both ends are more potent than singly adapted
ones. Cleanup uses 40 µl AMPure XP, two 250 µl washes with Long Fragment Buffer (enriches
fragments of about 3 kb and longer) or Short Fragment Buffer (retains all sizes), and
15 µl Elution Buffer (37 °C can improve recovery of high-molecular-weight DNA). Ethanol
must not be used here because it can denature the adapter motor protein.

🟡 Model correction (2026-10-17): the graph now includes both bead cleanups and their
discarded supernatants. Earlier it went directly from end prep to ligation to sequencing.

## 5. Motor, tether and sequencing

🟢 LA is supplied with a motor protein. The motor associates with a nanopore, unzips the
duplex and controls translocation of one intact strand. This is direct sequencing: there
is no conventional sequencing primer.

🟢 Flow Cell Tether is mixed into the flow-cell priming solution and concentrates library
at the membrane. The inspected SQK-LSK114 protocol, chemistry document and sensing
explanation do not specify a covalent tether-to-library bond, so none is drawn.

🔴 The complete Y-shaped leader sequence, current motor identity, attachment chemistry,
full tether structure and definitive strand polarity remain unverified in this source
audit. The older sensing explanation describes motor-on-leader strand capture but is
not a version-specific molecular inventory for SQK-LSK114. The page therefore shows
a named unresolved leader/motor feature and does not invent bases or a covalent
tether segment. An independent source check of that explanation and its links remains
open before any absence claim can be accepted.

## Native graph scope, 2026-10-18

Native graph migration (2026-10-18): source-qualified LA strands and the selected
no-DCS genomic-DNA route are owned locally. Explicit prepared ends and paired
adapter strands feed the ligation geometry; the terminal T is represented once.
The 61 µl eluate is split into a 1 µl Qubit aliquot and 60 µl ligation input;
measurement-aliquot fate is unresolved. The disclosed adapter DNA and sealed
junctions are separate from the associated, unplaced motor/leader component.
No protein or speculative leader is appended as a covalent DNA segment. The
September 2025 workflow and July 2026 oligo reference retain distinct identities;
they do not establish full historical reagent equivalence. Optional DCS is outside
this selected no-control route. Source-known operation order and both physical
cleanups remain explicit. This is a native graph migration within the stated source
scope, not complete leader/motor validation or an accepted absence audit.

Independent implementation review retained unligated/partly adapted molecular
classes and unresolved later bead recovery; ligation does not establish complete
adaptation or physical purity.
This code review read no additional originals and does not count as an independent
second source-reading pass.
