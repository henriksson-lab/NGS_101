# ProBac-seq — source notes

Defining protocol: Samanta et al., *Nature Protocols* (2024), [doi:10.1038/s41596-024-01002-1](https://doi.org/10.1038/s41596-024-01002-1).

## Sources read

- 🟢 The PMC author manuscript, including the complete stepwise protocol and figure legends.
- 🟢 Supplementary Table 1, containing the ordered *B. subtilis* probe pool.
- 🟢 Supplementary Table 2, containing the RCA, blocked-extension and library-PCR primers.

## Probe preparation

- 🟢 Each of the 29,765 ordered protoprobes in Supplementary Table 1 is 5′-`AAGCTTTGCTCGAAGAAACTATG`–target-specific region–`TGTATGGACCGGGTCATC`-3′: the 23-nt 5′ flank and 18-nt 3′ flank are the exact common prefix and suffix of the whole pool (computed over the table). The 18-nt 3′ flank is the exact reverse complement of the InDrop reverse primer; the model enforces that relation. (An earlier model used a 17-nt `GTATGG…` handle, leaving the InDrop primer's 3′-terminal base unpaired; corrected 2026-10-17.)
- 🟢 Published discrepancy: the article's experimental-design text says a common 23 bp PCR handle is added to the probe 5′ end and a 17 bp HindIII-containing extension handle to the 3′ end. The table instead orders a 23-nt HindIII-containing 5′ flank and an 18-nt InDrop-complementary 3′ flank; in the final mRNA-complementary orientation the positions agree but the stated lengths are swapped. The model and page use the table sequences.
- 🟡 The table's variable regions between those flanks are 55–60 nt, whereas the text designs 50 bp binding sequences. Whether the extra bases are spacers or longer binding regions is not stated in the inspected passages; the drawing keeps a schematic 50-nt placeholder and declares no known molecular length for it.
- 🟢 Step 1 phosphorylates the crude Twist pool with T4 polynucleotide kinase (37 °C 30 min, 65 °C 20 min) before the first splinted T4 DNA ligation.
- 🟢 RCA-1 and RCA-3 use 5′-`TTCGAGCAAAGCTTGATGACCCGGTCC`-3′. RCA-2 uses its reverse complement, 5′-`GGACCGGGTCATCAAGCTTTGCTCGAA`-3′. Each primer bridges the two probe ends for circularization, primes phi29 and creates the HindIII-containing duplex needed to release unit-length products.
- 🟢 Three successive ligation–RCA–HindIII cycles alternate product orientation. The third releases the strand complementary to bacterial mRNA.
- 🟢 After the final HindIII digestion (directed by the RCA-2 "nicking" primer) the protoprobes are ethanol-precipitated and the ~100–120 nt band is excised from a 15% TBE–urea gel, removing primers and short digestion products (Steps 15–24).
- 🟢 The blocked extension template is 5′-T30–V–N11–`AAGCTTGCTCGAAGAAACTATG`-/3Phos/-3′. Its twelve degenerate positions and T30 tract are copied onto the protoprobe as a 12-nt probe UMI followed by poly(A)30. The copy is made by DNA polymerase I large (Klenow) fragment at room temperature for 30 min; the product is purified again by 15% TBE–urea gel (Steps 20–24), separating the unextended blocked primer.
- 🟢 The two supplementary tables disagree by one T at this arm: every ordered probe begins `AAGCTTTGCT…`, whereas the extension template ends `AAGCTTGCT…`. The public model retains both transcriptions and does not draw them as an exact fully paired duplex.

## Cellular barcoding and final library

- 🟢 Fixed bacteria bearing hybridized probes are encapsulated with 10x 3′ beads (overview: “V3 chemistry”; Materials: 10x Next GEM chip G kit PN-1000128 and 10x V3.1 single cell reagent kit 2000164). The ProBac-seq sources print only the 10X forward primer that matches the 5′ end of the released gel-bead oligo. 🟢 The gel-bead oligo is taken from the vendor guide for that kit, 10x CG000204 Rev D (Chromium Next GEM Single Cell 3′ Reagent Kits v3.1; in this protocol's source bundle as `10x_CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf.txt`), whose kit list names PN-1000128 and Gel Beads PN 2000164 (lines 11, 308, 330, 422) and whose Appendix “Oligonucleotide Sequences”, Protocol Step 1.5 (line 5139), prints 5′-`CTACACGACGCTCTTCCGATCT`-N16-N12-T30VN-3′: partial Read 1, a 16-nt cell barcode, a 12-nt GEM UMI and poly(dT)30VN. ProBac-seq names no guide revision, so Rev D is an explicit choice. 🟡 The VN anchor is not drawn as a separate segment; it pairs within the first two copied probe-UMI positions. Six-cycle in-GEM PCR copies bound DNA probes without reverse transcription.
- 🟢 PCR-1 and PCR-2 use the 18-nt InDrop reverse primer 5′-`GATGACCCGGTCCATACA`-3′; PCR-2 uses the source-listed 10x forward primer 5′-`CTACACGACGCTCTTCCGATCT`-3′.
- 🟢 PCR-3 adds the full P5/Read 1 arm and a P7/index/Read 2 arm. Supplementary Table 2 provides sixteen distinct eight-base i7 primers.
- 🟢 The stated NextSeq 2000 run acquires Index 1 for 8 cycles and Read 1 for 119 cycles with 30% PhiX. Read 1 contains, in order, the 10x cell barcode, GEM UMI, poly(dT), probe UMI, extension handle and the beginning of the probe identity. The analysis splits this long read and substitutes the probe UMI for the GEM UMI when counting molecules.
- 🟢 The commercial bead oligo is not reprinted in the ProBac-seq Supplementary Table 2; the schematic draws it from CG000204 Rev D as above (Gel Beads PN 2000164 as named in the ProBac-seq Materials).
