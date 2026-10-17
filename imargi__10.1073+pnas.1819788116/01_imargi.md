# iMARGI — research notes

## Sources read

- 🟢 Yan et al. 2019, DOI `10.1073/pnas.1819788116`, PMC full text and supplement.
- 🟢 Wu et al. 2019, DOI `10.1038/s41596-019-0229-4`, PMC detailed protocol, especially the “iMARGI linkers and library PCR primers” section and Fig. 1.
- 🟢 Sridhar et al. 2017, DOI `10.1016/j.cub.2017.01.011`, original MARGI author manuscript and publisher Document S1 Supplemental Experimental Procedures (`1-s2.0-S0960982217300118-mmc1.pdf.txt`, lines 178–325). Two independent readers checked the original and linked supplement; the checklist is in `catalogue/source_reviews/`.

## Published citation discrepancy

🟢 Yan et al. 2019 reference 10 points to DOI `10.1016/j.cub.2017.01.068` for Sridhar et al.'s MARGI method. [PubMed's record for that DOI](https://pubmed.ncbi.nlm.nih.gov/28222283/) identifies it as an erratum to the [original method article](https://pmc.ncbi.nlm.nih.gov/articles/PMC5319903/), DOI `10.1016/j.cub.2017.01.011`. The model uses the iMARGI oligos printed by Wu et al. 2019. The recovered MARGI supplement prints the inherited linker bases and processing steps.

🟢 **Linker modification varies across source versions.** Sridhar et al. 2017 original
Document S1, lines 204–215, prints the top oligo ending
`...TCCGA/BiodT/CT`, with no phosphorothioate `*` mark. Wu et al. 2019
detailed iMARGI protocol, JATS lines 275–285, prints
`...TCC*GA/iBiodT/CT` and explicitly defines `C*` as phosphorothioated C.
The unmodified base order agrees. The model follows the later iMARGI-specific
ordering string, including `C*`; the 2019 PNAS paper's “same linker sequence”
statement does not resolve whether its actual linker included this modification.
The two biotin spellings both denote internal biotin-dT, but direct identity of
the ordered oligos is not established by the printed sources.

🟢 Wu et al. 2019 detailed protocol, extracted JATS lines 395–397, labels a stock
“10 mM dATP” but says to dilute **ATP** to prepare it. Its dA-tailing reaction
table at lines 673–685 instead calls for **10 mM dATP**, and Yan et al. 2019
Methods line 203 also specifies dATP for Klenow dA-tailing. The ATP wording in
the stock recipe appears to be a publication error. The model uses a genomic
3′-dA tail, supported by the reaction table and original paper; no numeric
stock preparation is modeled. The actual reagent placed in the published
experiment cannot be independently verified from these texts.

## Construct path

🟢 iMARGI improves MARGI by performing the linker ligations in situ. Fragmented chromatin RNA is joined to one end of an annealed bivalent linker, and the linker’s DNA end is joined to nearby restriction-cut genomic DNA. Reverse transcription converts the RNA side, and paired-end sequencing reports RNA and DNA partners on opposite sides of the bridge.

## Exact published structural oligos

🟢 The detailed protocol prints both linker strands, the BamHI cut oligo and the PCR pair verbatim:

- Linker top: `5′-/Phos/NNAGATCGGAAGAGCGTCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCTTCC*GA/iBiodT/CT-3′`. `NN` reduces RNA-ligation sequence bias; the source defines `C*` as a phosphorothioated C and `/iBiodT/` as internal biotin-dT.
- Linker bottom: `5′-/Phos/GATCGGAAGAGCACACGTCTGAACGGATCCTCCCTACACGACGCTCT-3′`.
- BamHI cut oligo: `5′-TCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCT/3InvdT/-3′`.
- NEBNext universal PCR: `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′`.
- NEBNext index PCR: `5′-CAAGCAGAAGACGGCATACGAGAT[6 bp barcode]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-S-T-3′`, retaining the protocol's printed terminal `-S-T` notation.

🟢 The top strand is 5′ adenylated after synthesis and annealed to the bottom strand. The bottom linker strand itself primes first-strand cDNA synthesis; the prior model's separate RT primer was unsupported and is removed. The cDNA is circularized, annealed to the cut oligo, and digested at the linker BamHI site. The split linker halves then serve as PCR landing sites. The run is paired-end 100 cycles. 🔴 The exact six-base sample-index member used in each library is not specified here; the source prints only a degenerate index slot. 🟡 The compact schematic still leaves the strand-resolved post-BamHI junction as a declared geometry limitation.

🟢 The DNA-side ligation is a cohesive-end reaction, not an undisclosed end:
the detailed iMARGI protocol lines 101–103 specifies dA-tailing and sticky-end
ligation, and its printed linker strands at lines 273–285 pair for 47 bases,
leaving a terminal 3′-T on the top strand. The original MARGI Fig. 1 text
likewise says the linker 3′-T joins the DNA 3′-A tail. The model now draws
these complementary one-base ends; it does not infer the downstream BamHI
product junction from this ligation view.

## Reaction order and read layout (bounded pass 2026-10-17, round 11)

🟢 Yan et al. 2019 Methods (PMC JATS lines 197–207) and Wu et al. 2019 Procedure
(JATS lines 900–1240, Steps 66–114) give this order, now declared step by step in
the protocol graph:

1. AluI chromatin digestion overnight in SDS-permeabilized, Triton-quenched nuclei; RNase I (1:10, 3 min, 37 °C) RNA fragmentation.
2. T4 PNK in pH 6.5 phosphatase buffer converts RNA 3′ phosphates to 3′-OH; Klenow exo− with 0.1 mM dATP adds a genomic 3′-dA.
3. In situ RNA–linker ligation with T4 RNA ligase 2 truncated KQ and 20% PEG 8000 (22 °C 6 h, 16 °C overnight); five washes remove excess free linker.
4. In situ DNA proximity ligation with T4 DNA ligase in 2 mL.
5. Proteinase K/SDS cross-link reversal at 65 °C, phenol:chloroform extraction and ethanol precipitation.
6. Exonuclease I then T4 DNA polymerase (dATP + dGTP, 12 °C) remove biotin from linkers not ligated to DNA.
7. MyOne Streptavidin C1 capture with seven 4 M NaCl high-salt washes.
8. On-bead SuperScript III RT (50 °C) primed by the linker bottom strand; T4 PNK + ATP 5′ phosphorylation on beads.
9. 0.1 M NaOH releases the cDNA-containing strand; 98 °C 20 min hydrolyses RNA; neutralization and SILANE purification.
10. CircLigase circularization (60 °C 4 h), Cut_oligo annealing, BamHI relinearization, SILANE purification.
11. NEBNext High-Fidelity 2× PCR with the Universal and Index primers (12 cycles in the authors' experiments); MinElute and 250–1000 bp 2% E-gel size selection.
12. Paired-end 100-cycle sequencing on HiSeq 4000.

🟢 Read 1 is the RNA end and Read 2 the DNA end (Yan et al. Methods line 207).
The first two Read 1 bases are the linker `NN` and are trimmed; read pairs whose
Read 2 does not start with AluI-derived `CT` are filtered (line 213). 🟡 The
compact final-library view remains inferred; the strand-resolved post-BamHI
junction is not drawn.

🟢 **Within-source reagent-name discrepancy.** Wu et al. 2019 Step 104 (JATS
line 1159) names “2× NEBNext HC PCR Master Mix” for the pilot PCR, while
Step 106 and the reagent list (line 259) name NEBNext High-Fidelity 2× PCR
Master Mix (M0541S). The model uses the High-Fidelity mix; “HC” appears to be
a typo, but no other source resolves it.

Not inspected in this capped pass: the MARGI Document S1 library steps were not
re-read against the round-11 graph (prior passes record lines 178–325); the
iMARGI_methods GitHub pipeline; the selected NEBNext index member (unassigned).
