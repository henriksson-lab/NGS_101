# BL-Hi-C: source-qualified standard HaeIII / S2 route

Liang et al., *BL-Hi-C is an efficient and sensitive approach for capturing structural and regulatory chromatin interactions*, Nature Communications (2017), DOI [10.1038/s41467-017-01754-3](https://doi.org/10.1038/s41467-017-01754-3).

Private incomplete graph: `tools/chem.py`. This bounded source pass is not a full source review or an accepted missing-information claim. The selected published route uses HaeIII and the two-step S2 bridge; the S1 one-step and overhang-enzyme alternatives are not merged into it.

## Sources inspected

🟢 Published PMC5696378 article, Methods90–106 and Fig1; published MOESM1 chemistry supplement, lines734–957. Original SI pages10–11 oligos and page12 fixation were rendered and inspected. Mandatory primer scrape returned30candidates/3files; relevant hits were checked in the original. The peer-review file was not read as chemistry authority.

## Exact oligos and modifications

🟢 S2-F `5′-/5P/CGCGATATC/iBIOdT/TATCTGACT-3′`; S2-R `5′-/5P/GTCAGATAAGATATCGCGT-3′`. The internal biotin-dT is an actual thymidine, not a zero-length modifier: each strand has19nucleotides. 🟡 The first18 bases pair and each strand retains one terminal3′T overhang; the duplex spans20columns. The local typed bridge retains the actual biotin-T position and both5′phosphates. This fixed bridge has no cell-barcode or UMI feature.

🟢 Y-adapter-F `5′-/5P/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC-3′`; Y-adapter-R `5′-TACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′`. A phosphate on Y-R is not asserted. The source annealing operation is retained; the complete genomic/Y junction is unassembled.

🟢 Universal1.0 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGAC`; index1.0 `CAAGCAGAAGACGGCATACGAGAT/CGTGATGT/GACTGGAGTTCAGACGTGT`; NTR-F `AATGATACGGCGACCACCGA`; NTR-R `CAAGCAGAAGACGGCATACGA`. Original p11 explicitly brackets the eight-base `CGTGATGT` window with `/index/`. The literal window has a structured sample-index feature. It is not silently replaced by the familiar six-base index interpretation, and validated index-read/sample-sheet geometry is not claimed.

## Source order and fractions

🟢 Fixation/PBS washes → sequential SDS lysis/pellets → HaeIII → direct Klenow-exo−/dATP A-tailing → S2/T4 ligation16°C4h → pellet → lambda exonuclease/ExoI → PK/SDS and NaCl heating → PCI/MaXtract aqueous recovery and isopropanol precipitation → Covaris400bp target. The one continuous ligation reaction can form two genomic joins; no second independent reaction is invented. Exonuclease debris remains physical until extraction.

🟢 M280 beads are washed/I-Block treated and washed, then blocked with1µg Salmon Sperm DNA and washed before contact-DNA capture. Blocking DNA is an independent physical input and possible carryover through all later work; proteolysis cannot erase DNA. Source capture and actual washes precede bead-bound end repair, A-tailing, Y-adapter ligation and further washes. No perfect capture specificity or wash partition is assumed.

🟢 The60µl bead suspension is divided into20µl aliquots stored at−20°C; one aliquot feeds four-primer Q5 PCR9–12cycles and the independent40µl remainder is stored. PCR product is transferred without Dynabeads. The first0.6xAMPure stage retains the **supernatant**; adding7.5µl beads to cumulative0.75x retains the **beads**. Following80%EtOH washes, DNA is eluted in20µl water. Qubit/Bioanalyzer consume separate analytical aliquots, while the production remainder feeds HiSeq X Ten PE150 acquisition.

## Explicit source boundary

🟢 Original SI p12 prints1ml medium, then278µl37%formaldehyde “final1%,” and894µl2.5M glycine “final0.2M.” These numerical formulation statements are incompatible as written. The graph retains the known timing and source wording but leaves numeric fixation composition unassembled; it does not change1ml to10ml or assert an erratum.

Full genomic contact, end/junction/PCR and sequencing-primer geometry remain unassembled. Native contact-DNA/biotin and PCR copies are physical molecules; contact calls are downstream computation. No source completeness or missing-sequence claim follows from this pass. Source-known alternatives and exact historical component originals remain deferred, rather than filled from Hi-C family recipes.
