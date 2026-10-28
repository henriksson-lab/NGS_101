# GLOE-Seq: capture native 3′-OH termini before fragmentation

Sriramachandran et al. (2020), *Genome-wide Nucleotide-Resolution Mapping of DNA Replication Patterns, Single-Strand Breaks, and Lesions by GLOE-Seq*. [Published article](https://doi.org/10.1016/j.molcel.2020.03.027).

## Sources read in the capped pass

One reader, 28 October 2026; three originals: published PMC7276987 JATS (Introduction; Results “Development of the GLOE-Seq Method”, EndoSeq comparison; STAR Methods “GLOE-Seq Library Preparation”, sequencing and GLOE-Pipe); published Methods S1 (12 pages, oligos p4, steps y23–50/m22–28, notes7–17); published Supplemental Information Figure S1A, p2, inspected as a rendered original as well as extracted text. The merged mmc3 PDF is not a fourth independent source. Manifest confirms published article, not an author manuscript. Mandatory oligo scraper run; sequence/modification hits checked against Methods S1 and Figure S1. This is a capped pass, not a completed source/absence audit.

## What is distinctive

🟢 Capture the original free 3′-OH with a phosphorylated, biotinylated proximal adapter while the DNA remains otherwise intact. Sonication follows capture ligation, so new sonication ends do not define the captured break position. Mammalian lysis and ligation take place in agarose plugs; yeast uses gently extracted DNA. Both converge on steps29–50. SSBs and 3′ ends of DSBs can contribute: a read is not proof of an SSB. Lesion-specific nuclease treatments are separate optional front ends, not an assumed property of native-end libraries.

## Oligos, verbatim 5′→3′

Methods S1 p4. `*` is a phosphorothioate bond; NH2 is the 3′ amino modification; Bio is terminal 3′ biotin-dT. The six random splint bases are hybridization degeneracy, not a UMI or sequenced sample barcode.

| Name | Source sequence |
|---|---|
| 3898 splinter | CTACACGACGCTCTTCCGATCTNNNNN*N-NH2 |
| 3899 proximal | PO4-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGTTTT-Bio |
| 3790 second-strand primer | CGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT |
| 3791 distal | GACTGGAGTTCAGACGTGTGCTCTTCCGATCT |
| 3792 distal | GATCGGAAGAGCACACGTCTGAACTCCAGTC |
| P5 | AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT |
| P7 | CAAGCAGAAGACGGCATACGAGAT(X)6GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT |

## Source-defined reaction order

1. 🟢 Optional lesion conversion occurs before denaturation: T4 endonuclease V plus APE1 for pyrimidine dimers, APE1 for AP sites, RNase H for ribonucleotides. These create a different analytical substrate from native-end capture. Yeast Okazaki-fragment experiment additionally dephosphorylates genomic DNA with Antarctic Phosphatase (STAR Methods), not a universal core step.
2. 🟢 Denature DNA and quench on ice. T4 DNA ligase seals the native 3′-OH to adapter3899 5′ phosphate;3898 bridges its end with six random terminal bases. Its 3′ amine prevents extension/ligation at that end; source specifies phosphorothioate near that end. Keep adapter annealing separate from covalent incorporation.
3. 🟢 Sonicate to about200nt, perform AMPure cleanups, capture biotin on streptavidin beads. Wash twice in SSC and once in20mM NaOH to remove the splinter. Heat elute at95°C. FigureS1 depicts loss of the biotin linkage upon release, but this pass does not assign a specific covalent cleavage product or chemically alter the DNA sequence on that basis.
4. 🟢 Primer3790 and Phusion FLASH synthesize the complementary strand; AMPure purification follows. NEBNext UltraII EndPrep performs end polishing/A-tailing. Distal adapter3791/3792 is ligated with the kit ligation mix; salt/PEG and two rounds of AMPure cleanup precede PCR. FigureS1 shows the terminal T/A junction; it is not a blunt adapter.
5. 🟢 Q5 with P5 and indexed P7 amplifies for eight cycles. Two AMPure cleanups and Qubit quantification finish the library. The six-base index is a selected P7 member, not a fixed study-wide value.

## Read layout and analytical coordinate

🟢 STAR Methods: NextSeq500, single-end, read lengths75/84/92/150/160 depending on run; seven index cycles for single-index libraries or8+8 for dual-index libraries. MethodsS1 prints a six-base P7-index family. Do not replace that family with an invented eight-base index or silently equate molecular index length and acquired cycles. Exact dual-index molecular architecture/run assignments remain open in this pass.

🟢 GLOE-Pipe explicitly states reads are reverse complements of captured fragments. It assigns the read5′ end to the opposite genomic strand. Direct mode reports the original fragment3′ terminal nucleotide; indirect mode shifts to the base immediately3′ of the break for treatments cleaving5′ of a lesion. The physical captured end and the computationally shifted lesion coordinate are different claims.

## Ancestry and boundaries

🟢 Results “Development”: Gansauge et al.2017 random-splinter/T4DNA-ligase approach is an explicitly cited molecular inspiration, applied here to intact heat-denatured DNA. It does not establish inheritance of all their library reagents or adapter sequences. EndoSeq/emRiboSeq use a related splinter approach, but Results “Comparison” explicitly distinguishes order: GLOE proximal ligation precedes fragmentation; EndoSeq fragmentation/distal ligation precede lesion conversion and denaturation. This is a comparison, not wholesale recipe inheritance. OK-seq is a biological readout comparison, not a parent supplying GLOE chemistry.

## Deferred, without absence claims

Full exact reaction validation, including native-end splint geometry across random sequence classes, bead-release chemistry and distal cohesive junction; individual run/index configurations; linked Gansauge2017 and emRiboSeq originals; second reader. These are model/inspection tasks, not established omissions from the publication. A source-complete universal dual-index model is not justified by this capped pass.
