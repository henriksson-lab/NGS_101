# Single-stranded library preparation for ancient or damaged DNA

## 1. What it is

This protocol recovers individual highly degraded DNA strands instead of requiring an intact duplex with two usable ends.

## 2. Sources and evidence

- 🟢 Gansauge MT, Meyer M. *Single-stranded DNA library preparation for the sequencing of ancient or damaged DNA.* Nature Protocols 8:737–748 (2013). [doi:10.1038/nprot.2013.038](https://doi.org/10.1038/nprot.2013.038).
- 🟢 A phosphorylated, biotinylated adapter is ligated to the 3′ ends of heat-denatured DNA, captured on streptavidin beads, copied, and joined to a second adapter by blunt ligation.
- 🟢 Published oligos include CL78 with internal C3 spacers and terminal TEG-biotin, the CL53/CL73 second adapter, CL9 extension primer and CL72 custom Read 1 primer.
- 🟢 CL72 is `ACACTCTTTCCCTACACGACGCTCTTCC`; the custom primer is required because the P5-side adapter is truncated relative to ordinary TruSeq libraries.
- 🟢 The publisher's [15 May 2013 correction](https://www.nature.com/articles/nprot.2013.038) explicitly restores the missing 5′ phosphate on CL78. The corrected oligo is `5′-[Phosphate]AGATCGGAAG[C3Spacer]10[TEG-biotin]-3′`, not an unphosphorylated adapter. The author group's [2017 supplementary Table 2](https://pure.mpg.de/pubman/item/item_2459133_8/component/file_2459156/Gansauge_Single-stranded_NAR_2017_Suppl.pdf) republishes the original CL oligos and identifies which were already in the 2013 protocol: CL9 `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`, CL53 `CGACGCTCTTC-ddC`, and CL73 `Pho-GGAAGAGCGTCGTGTAGGGAAAGAG*T*G*T*A` (asterisks mark phosphorothioate bonds). These are not interchangeable with the later paper's new splinter oligos.
- 🟢 The 2013 Supplementary Table 2 supplies indexed P5 and P7 primer families. Its members differ at the 7-base index window; there is no single universal index-primer pair. The diagram's sample-index placeholders represent a selected family member, not an unpublished barcode.

## 3. Molecular path

Every denatured template strand is independently captured through its 3′ end. Copying converts it to a duplex only after capture, avoiding the need for complementary ancient strands to survive together.
