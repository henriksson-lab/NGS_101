# TRIP — Thousands of Reporters Integrated in Parallel

## Citation

Akhtar W, de Jong J, Pindyurin AV, Pagie L, Meuleman W, de Ridder J, Berns A,
Wessels LFA, van Lohuizen M, van Steensel B. *Using TRIP for genome-wide position
effect analysis in cultured cells.* Nature Protocols 9, 1255–1281 (2014).
[doi:10.1038/nprot.2014.072](https://doi.org/10.1038/nprot.2014.072), PMID 24810036.

## Sources read

- 🔴 Main protocol text and reagent list: **not available in this checkout**. The
  publisher full text is paywalled (`MANIFEST.tsv` keeps the manual URL
  <https://doi.org/10.1038/nprot.2014.072>), so the barcode-cloning, amplification and
  inverse-PCR oligos other than the indexed forward primer are inherited transcriptions
  that this pass could not re-read at the source. Cycling conditions and the permitted
  read lengths come from the same unread text.
- 🟢 Supplementary Table 1 (`…MOESM451_ESM.pdf.txt`, names lines 40–87, sequences lines
  89–136): all 48 ten-base-indexed `PB-cDNA-forward-1-BC` primers; index 1 is line 89
  and the 48 rows differ only at nt 34–43.
- 🟢 Supplementary Figure 3 (`…MOESM450_ESM.pdf.txt` lines 1–63): panel a, the
  expression/normalization library (lines 9–28); panel b, the mapping library
  (lines 35–59), including the Read 2 primer position.
- 🟡 Supplementary Figures 1 and 2 (`…MOESM448`, `…MOESM449`): iPCR artefact scenarios
  and the IR copy-number amplicons; the extracted text is figure labels only.

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

Only `PB-cDNA-forward-1-BC` is verified verbatim here (Supplementary Table 1 line 89).
The other five are consistent with the ends drawn in Supplementary Figure 3 —
`PB-cDNA-reverse` and `Solexa-R-Ad-TruSeqPE2` appear as the reverse complement of the
library 3′ end, `Solexa-Adapter-SE1` as the 5′ end, and the `InvPCR-F-TruSeqPE2`
annealing region as the piggyBac 3′-TR stretch — but their exact ordering sequences
come from the unread reagent list.

## Library structures from Supplementary Figure 3

🟢 Panel a prints the counting amplicon continuously. Between the reporter right
constant (`…GAATTCTTACTT`) and the `PB-cDNA-reverse` annealing site the amplicon runs
through the 3′ end of eGFP: `GTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACG`
({{= len(trip.EGFP_3END) }} nt). The model previously joined the reporter constant
straight to the reverse-primer site; that repository omission is corrected, and the
drawn library now reproduces the figure base for base.

🟢 Panel b prints `TTAACC` between the variable genomic flank and the inverse-PCR primer
site: the duplicated `TTAA` piggyBac target site followed by the `CC` that begins the
3′ terminal repeat. Both were missing from the model and are now drawn; the inverse-PCR
primer still anneals to the 28-nt repeat stretch behind them.

🟢 Panel b also confirms the mapping library ends `…AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`
+ P7 with no `GAT` spacer: PCR2 with `Solexa-R-Ad-TruSeqPE2` fuses P7 directly to the
Read 2 site, dropping the 5′ `GAT` of `InvPCR-F-TruSeqPE2`.

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

## Open questions

- The defining article text is paywalled and absent from the local bundle: reagent-list
  oligo transcriptions other than Supplementary Table 1, the cycling conditions, the
  65–100-cycle read allowance and the actual acquired read lengths all await it.
- Supplementary Figure 3 draws the eGFP body and the reporter cassette as continuous
  sequence only at the two ends shown; the full reporter plasmid map is not published
  with the protocol.
- `PB-barcode-long-7`, `BssH2-A` and the `eGFP-int3` sequencing primer appear in no
  available supplement and remain unverified transcriptions.
