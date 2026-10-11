# Oxford Nanopore Native Barcoding V14

## 1. What it is

The SQK-NBD114.24 workflow repairs and dA-tails genomic DNA, ligates a different
dT-tailed native barcode adapter to each sample, pools the samples, and ligates the
motor-bearing Native Adapter through the barcode adapter's cohesive end. It is PCR-free
and retains native base modifications.

Commercial defining source: Oxford Nanopore, “Ligation sequencing gDNA — Native
Barcoding Kit 24 V14 (SQK-NBD114.24),” NBE_9169_v114.

## 2. Sources read

- 🟢 [Oxford Nanopore protocol NBE_9169_v114](https://nanoporetech.com/document/ligation-sequencing-gdna-native-barcoding-v14-sqk-nbd114-24), including the public NB01–NB24 barcode table.
- 🟢 [Oxford Nanopore kit description](https://store.nanoporetech.com/native-barcoding-kit-24-v14.html), for the dT-tailed barcode adapter and cohesive-end hook architecture.

## 3. Molecular steps

🟢 FFPE repair and Ultra II end prep run together for 5 min at 20 °C and 5 min at
65 °C, yielding dA-tailed duplexes. A unique Native Barcode is ligated to each sample
with Blunt/TA Ligase Master Mix before bead cleanup and equimolar pooling. Native Adapter
is then ligated for 20 min at room temperature.

🟢 The kit description states that each native barcode adapter has a dT overhang for the
dA-tailed template and a cohesive end that hooks the sequencing adapter. The public
24-nt NB01 sequence is `CACAAAGACACCGACAACTTTCTT`; all 24 sequences are published in
the protocol.

🔴 Oxford Nanopore does not publish the barcode-adapter flanks, cohesive-end sequence,
Native Adapter sequence, motor attachment chemistry, or the exact strand carried into
the pore. The schematic uses bracketed role tokens for those parts and does not substitute
sequences from older kits.

## 4. Read structure

🟢 Nanopore sequencing is motor-controlled translocation rather than primer extension.
The barcode lies adjacent to the genomic insert and is recognized in the raw read for
demultiplexing. The drawing shows one productive pore-entry orientation; molecules can
carry barcode adapters at both insert ends.
