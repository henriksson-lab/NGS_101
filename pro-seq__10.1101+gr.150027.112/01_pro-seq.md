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

The molecular workflow now retains the two RNA ligations separately. VRA3 is drawn as a
5′-phosphorylated single-strand adapter whose 5′ end joins the run-on RNA 3′-OH; its
3′ inverted dT remains visible as the ligation blocker. After biotin selection, decapping
and 5′-end repair, VRA5's 3′-OH joins the RNA 5′ phosphate. Both covalent junctions are
marked on the dual-adapter RNA rather than collapsing them into a generic library icon.

RP1 anneals through its 21-nt 3′ suffix within VRA3 and contributes P5 as a non-annealing
tail. The RPI-n primer family anneals through the VRA5 sequence and contributes its
six-base i7 index plus P7. Primer landing, extension direction and amplicon boundaries are
derived from those exact Table 1 sequences. The final construct is therefore modeled from
the published PCR primers and adapter-derived landing sites, not as generic dual-index
TruSeq. Read 1 sequences the reverse complement of nascent RNA beginning at its
informative 3′ end; the separate index primer reads the six-base RPI-n barcode.
