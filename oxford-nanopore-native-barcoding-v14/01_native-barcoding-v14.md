# Oxford Nanopore Native Barcoding V14

## 1. What it is

The SQK-NBD114.24 workflow repairs and dA-tails genomic DNA, ligates a different
dT-tailed native barcode adapter to each sample, pools the samples, and ligates the
motor-bearing Native Adapter through the barcode adapter's cohesive end. It is PCR-free
and retains native base modifications.

Commercial defining source: Oxford Nanopore, “Ligation sequencing gDNA — Native
Barcoding Kit 24 V14 (SQK-NBD114.24),” **NBE_9169_v114_revX_15Jul2026** (the
revision read on 2026-10-17).

## 2. Sources read

- 🟢 [Oxford Nanopore protocol NBE_9169_v114](https://nanoporetech.com/document/ligation-sequencing-gdna-native-barcoding-v14-sqk-nbd114-24), revision X (15 Jul 2026): end-prep and cleanup (text lines 642–777), native barcode ligation, EDTA stop, pooling and 0.4X cleanup (lines 781–998), Native Adapter ligation and cleanup (lines 1002–1214), barcode table (lines 452–549; all 24 forward/reverse pairs match the model).
- 🟢 [Oxford Nanopore kit description](https://store.nanoporetech.com/native-barcoding-kit-24-v14.html), for the dT-tailed barcode adapter and cohesive-end hook architecture.
- 🟢 Oxford Nanopore, *Chemistry Technical Document*, CHTD_500 v1 rev AV, section 14 (native-barcode flanking sequences, lines 1139–1145) and section 15 (both Native Adapter strands, lines 2074–2081). This bundle's copy was fetched from `oxfordnanoporedx.com/document/chemistry-technical-document`; the copies fetched from `nanoporetech.com/document/chemistry-technical-document` (other ONT bundles) have a different HTML md5, but both print `CHTD_500_v1_revAV_24July2026` and "Last updated: 7/24/2026", and their extracted text is identical line for line (checked 2026-10-17); the HTML differs only in page-build markup. Section 14 lists the NB01–96 V14 native barcodes; it does not say the flanks apply to older kits such as EXP-NBD103, so the ARTIC page (EXP-NBD002/EXP-NBD103) does not inherit them.

## 3. Molecular steps

🟢 FFPE repair and Ultra II end prep run together for 5 min at 20 °C and 5 min at
65 °C, yielding dA-tailed duplexes (400 ng per sample for more than four barcodes,
1000 ng for four or fewer; optional DNA Control Sample). Each reaction is cleaned with
15 µl AMPure XP and 80% ethanol washes, eluted in water and quantified, and an
**equimolar mass of each sample** is taken into barcode ligation. A unique Native Barcode
is ligated to each sample with NEB Blunt/TA Ligase Master Mix for 20 min at room
temperature; EDTA stops each reaction; the **entire volume of every barcoded sample is
pooled** and cleaned with 0.4X AMPure XP, Short Fragment Buffer washes and water elution at
37 °C. Native Adapter is then ligated with NEBNext Quick Ligation Reaction Buffer and
Quick T4 DNA Ligase for 20 min at room temperature, followed by AMPure XP cleanup with Long
or Short Fragment Buffer washes (ethanol is prohibited because it denatures the motor) and
Elution Buffer at 37 °C.

🟡 Model correction (2026-10-17): the earlier model placed equimolar normalisation after
barcode ligation and omitted the EDTA stop and all three bead cleanups. The graph now
follows revision X: cleanup → equimolar sampling → barcode ligation → EDTA → pool →
cleanup → Native Adapter ligation → cleanup.

🟢 The kit description states that each native barcode adapter has a dT overhang for the
dA-tailed template and a cohesive end that hooks the sequencing adapter. The public
24-nt NB01 sequence is `CACAAAGACACCGACAACTTTCTT`; all 24 sequences are published in
the protocol.

🟢 The current Native Adapter top strand is `CCTGTACTTCGTTCAGTTACGTATTGCT`; its
published bottom strand is `ACGTAACTGAACGAAGTACAGG`.

🟢 The Chemistry Technical Document prints barcode flanking sequences “which add an
extra level of context during analysis”: forward `5′-AAGGTTAA-barcode-CAGCACCT-3′`,
reverse `5′-GGTGCTG-barcode-TTAACCTTAGCAAT-3′`.

🟡 Derived (computed in `native_barcoding.py`, which raises if the relations fail): the
reverse 5′ flank is the reverse complement of the forward 3′ flank without its final T,
so that T is drawn as the single dT tail ligated to the insert dA. The published Native
Adapter bottom strand pairs with the first 22 bases of the top strand, leaving a 6-nt
3′ overhang `ATTGCT`; the reverse barcode flank extends six bases beyond the forward
5′ flank and ends in its complement `AGCAAT`. This is drawn as the cohesive hook. ONT
presents the flanks as analysis context and does not publish a barcode-adapter oligo
inventory, so the strand placement, dT split and hook pairing remain model derivations;
the hook geometry is still declared unavailable.

🔴 The complete motor assembly, motor attachment chemistry and the exact strand carried
into the pore are not published in the inspected sources. The schematic keeps a role
token for the motor assembly.

## 4. Read structure

🟢 Nanopore sequencing is motor-controlled translocation rather than primer extension.
The barcode lies adjacent to the genomic insert and is recognized in the raw read for
demultiplexing. The drawing shows one productive pore-entry orientation; molecules can
carry barcode adapters at both insert ends.
