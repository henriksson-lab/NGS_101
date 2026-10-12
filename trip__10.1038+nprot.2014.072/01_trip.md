# TRIP — Thousands of Reporters Integrated in Parallel

## Citation

Akhtar W, de Jong J, Pindyurin AV, Pagie L, Meuleman W, de Ridder J, Berns A,
Wessels LFA, van Lohuizen M, van Steensel B. *Using TRIP for genome-wide position
effect analysis in cultured cells.* Nature Protocols 9, 1760–1769 (2014).
[doi:10.1038/nprot.2014.072](https://doi.org/10.1038/nprot.2014.072), PMID 24810036.

## Sources read

- 🟢 Main protocol text and reagent list, including the barcode, amplification and
  inverse-PCR oligos.
- 🟢 Supplementary Table 1, which prints all 48 ten-base-indexed
  `PB-cDNA-forward-1-BC` primers.
- 🟢 Supplementary Figure 3, which diagrams the normalization/expression and mapping
  library structures and the paired-end mapping orientation.

## What the protocol measures

Each piggyBac reporter receives a random 16-nt barcode before integration. Three
libraries retain that same identifier:

1. reporter cDNA barcode counts measure expression;
2. genomic reporter barcode counts normalize clone abundance;
3. an inverse-PCR library joins the barcode to the genomic sequence beside the
   piggyBac 3′ terminal repeat, identifying the integration site.

## Published oligos

🟢 `PB-barcode-long-7`:

`5′-GTGACACCTGCAGGATCA(N)16CTCGAGTTGTGGCCGGCCCTTGTGACTG-3′`

🟢 `PB-cDNA-forward-1-BC-` consists of the TruSeq Read 1 sequence, a 10-nt inline
sample index, and the reporter annealing sequence. Supplementary Table 1 index 1 is:

`5′-ACACTCTTTCCCTACACGACGCTCTTCCGATCTACGAGTGCGTGTCACAAGGGCCGGCCACAA-3′`

🟢 `PB-cDNA-reverse`:

`5′-CAAGCAGAAGACGGCATACGAGATGATCACATGGTCCTGCTGGAGTT-3′`

🟢 `Solexa-Adapter-SE1`:

`5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′`

🟢 `InvPCR-F-TruSeqPE2`:

`5′-GATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGTACGTCACAATATGATTATCTTTCTAG-3′`

🟢 `Solexa-R-Ad-TruSeqPE2`:

`5′-CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGT-3′`

## Library logic

🟢 Expression RNA is oligo(dT)-primed and reverse-transcribed. The cDNA and genomic
normalization aliquots are amplified across the same reporter barcode cassette. The
forward primer introduces a 10-nt inline sample index immediately after the Read 1
site; a second PCR with `Solexa-Adapter-SE1` completes P5.

🟢 For mapping, genomic DNA is digested with DpnII and self-ligated at low DNA
concentration. The published inverse-PCR primer faces the piggyBac 3′ terminal repeat;
after circularization it amplifies through the adjacent genomic flank toward the
barcoded reporter. A second PCR completes the P5 and P7 ends.

## Read layout

- Expression and normalization: Read 1 reports the 10-nt sample index followed by
  reporter-constant sequence and the 16-nt reporter barcode. No index read is required.
- Mapping: Read 1 identifies the sample and reporter barcode. Read 2 starts from the
  non-barcode end, traverses the piggyBac terminal repeat and then reports the unknown
  genomic integration flank.

The protocol permits 65–100 sequencing cycles; the page uses 100 cycles to show the
maximum published read content. Genomic-flank length is biological and therefore not
assigned a fixed value.
