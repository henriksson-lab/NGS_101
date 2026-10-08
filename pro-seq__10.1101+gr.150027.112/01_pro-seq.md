# PRO-seq — precision nuclear run-on sequencing

## 1. What it is

PRO-seq maps the active sites and directions of engaged RNA polymerases by allowing a short nuclear run-on with biotinylated nucleotides and sequencing the marked nascent RNA ends.

## 2. Sources and evidence

- 🟢 Kwak H et al. *Precise maps of RNA polymerase reveal how promoters direct initiation and pausing.* Genome Research (2013). [doi:10.1101/gr.150027.112](https://doi.org/10.1101/gr.150027.112).
- 🟢 Mahat DB et al. *Base-pair-resolution genome-wide mapping of active RNA polymerases using precision nuclear run-on (PRO-seq).* Nature Protocols (2016). [doi:10.1038/nprot.2016.086](https://doi.org/10.1038/nprot.2016.086).
- 🟢 Mahat et al. Table 1 gives VRA3 (`GAUCGUCGGACUGUAGAACUCUGAAC-/Inverted dT/`), VRA5 (`CCUUGGCACCCGAGAAUUCCA`), RP1 (`AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA`) and RPI-n (`CAAGCAGAAGACGGCATACGAGAT NNNNNN GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`) verbatim. VRA3 is ordered 5′ phosphorylated.
- 🟢 The detailed protocol specifies TruSeq Small RNA-compatible single-end sequencing. RPI-n carries one six-base i7 barcode; no i5/Index 2 read is present.
- 🟢 Biotin-NTP run-on marks the nascent RNA 3′ end and permits stringent streptavidin enrichment; RNA end processing and sequential adapter ligations precede RT and PCR.

## 3. Molecular path

The information-bearing boundary is the run-on RNA 3′ end. Its genomic coordinate and strand identify an engaged polymerase active site.

The final construct is therefore modeled from the published PCR primers and adapter-derived landing sites, not as generic dual-index TruSeq. Read 1 sequences the reverse complement of the nascent RNA beginning at its informative 3′ end; the separate index primer reads the six-base RPI-n barcode.
