# miCLIP for m6A and m6Am

## Sources read

- 🟢 Linder et al., *Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome*, Nature Methods (2015), [doi:10.1038/nmeth.3453](https://doi.org/10.1038/nmeth.3453). Author manuscript PMC4487409 (JATS text): Online Methods “Crosslinking of cellular RNA”, “3′-dephosphorylation, linker ligation & labeling”, “Library preparation” and “Read pre-processing” (text lines 153–177); in vitro mutagenesis assay (line 143); Supplementary Data README (MOESM79). The article prints **no oligonucleotide** and delegates the 3′ adapter, circularization and re-linearization to the iCLIP protocol (ref. 10).
- 🟢 Grozhik et al., *Mapping m6A at individual-nucleotide resolution using crosslinking and immunoprecipitation*, Methods in Molecular Biology (2017), [doi:10.1007/978-1-4939-6807-7_5](https://doi.org/10.1007/978-1-4939-6807-7_5), PMC5562447: section 2.3 “Oligonucleotides”, sections 3.4–3.12 and Notes 23–27. Same laboratory, detailed miCLIP protocol.
- 🟢 König et al., [iCLIP (2010)](https://doi.org/10.1038/nsmb.1838), reference 10 of Linder et al.: Supplementary Methods (publisher MOESM21 PDF, now on disk; identical bytes to the PMC3000544 supplement) pp. 1–3, L3 ligation, RT primers, circularization/BamHI linearization and PCR.
- 🟢 Huppertz et al., [the same-group detailed iCLIP protocol (2014)](https://doi.org/10.1016/j.ymeth.2013.10.011), original PDF pp. 2–3 (read in an earlier pass), prints L3-App and the Rtclip barcoded primer family.

## Oligos from the detailed 2017 protocol

🟢 L3 linker (HPLC purified): `/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/`. Cut oligo (PAGE purified): `GTTCAGGATCCACGACGCTCTTC/3ddC/`.

🟢 P5 Solexa PCR primer: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`. P3 Solexa PCR primer: `CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT`.

🟢 Sixteen RT primers share `/5phos/NN[4-base experimental barcode]NNNAGATCGGAAGAGCGTCGTGgatcCTGAACCGC`. The protocol prints barcodes `aacc`, `acaa`, `attg`, `aggt`, `cgcc`, `ccgg`, `ctaa`, `catt`, `gcca`, `gacc`, `ggtt`, `gtgg`, `tccg`, `tgcc`, `tatt`, `ttaa` in RT1–RT16 respectively. RT1 writes its first random position as lowercase `n` (`NnaaccNNN`), while the other rows use uppercase `N`. The 5′ phosphate and lowercase `gatc` are source modifications/case notation, not guessed bases. The local self-test checks these strings against the PMC5562447 text.

🟢 The 2010 iCLIP generation differs: L3 was an RNA/DNA linker `5′phosphate-UGAGAUCGGAAGAGCGGTTCAG-3′-Puromycin`, RT primers carried a 3-nt random plus 2-nt experimental barcode (`5′-phosphate-NNNCAAGATCGGAAGAGCGTCGTGGATCCTGAACCGCTC-3′` and two siblings), and the cut primer ended `…CGCTCTTCAAAA-3′`. The PCR primers are the same Illumina paired-end pair. 🟡 Linder et al.'s read processing removes a 4-nt experimental barcode and keeps a 2-nt second random part (Methods line 173), which matches the 2014/2017 `NN-BBBB-NNN` primers, not the 2010 `NNN-BB` primers. The model therefore uses the 2017 oligos, marked inferred.

## Molecular path (Linder 2015 conditions)

🟢 RNA fragmented to 30–130 nt (20 µg) is bound by 1–5 µg anti-m6A antibody and crosslinked twice with 0.15 J cm⁻² 254-nm UV. Protein A/G beads capture the complexes, followed by high-salt, IP-buffer and PNK-buffer washes.

🟢 On beads, PNK in pH 6.5 buffer dephosphorylates RNA 3′ ends (30 min). After washes the 3′ adaptor is ligated with T4 RNA ligase (NEB) overnight. Complexes are eluted in NuPAGE buffer, run on SDS-PAGE, transferred to nitrocellulose, and the heavy- and light-chain regions are excised. Proteinase K releases the RNA, leaving a peptide at the crosslink, followed by phenol–chloroform extraction.

🟢 SuperScript III either reads through the adduct with a characteristic substitution or stops at it. These are alternative molecular outcomes. First-strand cDNA is size-selected on a 6% TBE-urea gel (80–120 nt for total RNA, 70–100 nt for mRNA). It is then circularized and re-linearized as in iCLIP, PCR-amplified for 18–21 cycles and sequenced on a HiSeq 2500. Short-insert (~30–40 nt) libraries were sequenced paired-end (Methods line 177).

🟡 With the 2017 oligos, the cut oligo pairs with the RT primer's `{{= chemdraw.revcomp(miclip.CUT_OLIGO) }}` region and restores a BamHI site. A single G^GATCC cut opens the circle without excising cDNA, which moves the barcode and random bases next to the P5-side arm. P5 Solexa then primes from the RT-primer arm and P3 Solexa from the L3 copy. These relations are asserted in `miclip.py`, not hand-aligned.

## Final library and reads

{{= ' '.join(s.name for s in miclip.FINAL_LIBRARY) }}

🟡 Read 1 is primed on the TruSeq Read 1 arm that the P5 Solexa primer carries. It reads NNN–BBBB–NN, where BBBB is the reverse complement of the RT-primer barcode (2017 Note 27: “barcode 1 is NNNGGTT”), and then RNA starting at the base after an RT stop. Linder et al. state neither read lengths nor a Read 2 primer. The model uses the 2017 recommendation of 50-cycle paired-end reads (Note 25), declares only Read 1, and leaves Read 2 undrawn.

## Model corrections (2026-10-17 bounded pass)

- The earlier model placed the Illumina small-RNA RA5 sequence (`GTTCAGAGTTCTACAGTCCGACGATC`) as the Read 1 site, with 51 cycles. No inspected miCLIP or iCLIP original uses that site. It was replaced by the P5 Solexa/TruSeq Read 1 arm computed from the printed primers.
- The ligase was relabelled as T4 RNA ligase, as printed, instead of “truncated” ligase. The model also adds on-bead PNK dephosphorylation, SDS-PAGE/membrane selection, proteinase K release and the TBE-urea cDNA size selection.
- An “excised restriction fragment” waste product was removed. A single BamHI cut linearizes a circle without releasing any cDNA.

## Open questions

- Which RT barcode (and, strictly, which oligo generation) each 2015 library used. No inspected source maps the 2015 samples to RT1–RT16.
- The read lengths actually acquired in 2015 and the Read 2 primer. The GEO/SRA run records were not inspected in this pass.
- Exact CircLigase and BamHI conditions in 2015 (“as described in the iCLIP protocol”). The 2010 and 2017 conditions differ, and neither is asserted for 2015.
