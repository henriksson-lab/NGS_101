# PacBio telobait HiFi sequencing

## 1. What it is

PacBio telobait sequencing ligates phased, barcoded and biotinylated duplex baits to native human telomere G-rich overhangs, affinity-purifies intact chromosome ends, releases them with EcoRI, then makes SMRTbell HiFi libraries.

## 2. Sources read

- 🟢 Tham et al. (2023), Nat Commun 14:281, [doi:10.1038/s41467-023-35823-7](https://doi.org/10.1038/s41467-023-35823-7), PMC9845338 JATS (`s41467-023-35823-7_PMC9845338.xml.txt`): Results "Telomere enrichment using telobaits" (line 128); Methods "Synthesis of barcoded Telobaits" (lines 218–220), "Telomere enrichment using telobaits" (line 224), "Analysis of Telobait-enriched PacBio HiFi sequencing data" (lines 230–232).
- 🟢 Supplementary Information (MOESM3), Figure S2a drawing and legend (lines 155–330).
- 🟢 Supplementary Data 2 (MOESM5; same md5 as the earlier hand-added `supp_data_2.xlsx`), Set 1–12 F01 and R01-1…R01-6 (lines 1–9).
- Supplementary Data 1 (MOESM4) is per-run ZMW/demultiplexing data, not chemistry.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), a review; comparison context only, not used for any base or reaction.
- Not read: the PacBio SMRTbell Express Template Prep Kit 2.0 guide named in Methods line 224.

## 3. Telobait oligos

🟢 Each telobait is a duplex. The forward oligo (one of 48 barcoded sets) is `/5Phos/ATCGACGGTTCAAGATGCCAGATGCACGGAGCAGAATTCTTTTAC/3Bio/` (F01). Its six reverse partners are the reverse complement of F01 followed by a 3′ six-base C-rich phase arm, e.g. R01-1 `5'-GTAAAAGAATTCTGCTCCGTGCATCTGGCATCTTGAACCGTCGATCCCTAA-3'`; the six arms are `CCCTAA`, `ACCCTA`, `AACCCT`, `TAACCC`, `CTAACC` and `CCTAAC`. Forward and reverse oligos are annealed 1:6 (95 °C 10 min, slow cool; Methods line 220). F carries an 8-nt barcode after `ATCGA`, the EcoRI site `GAATTC` and the 3′ biotin (Fig. S2a).

🟡 Derived geometry (Fig. S2a, computed in `telobait.py`): the R phase arm pairs with the terminal six bases of the G-rich overhang, so ligase joins the native G-rich 3′ end to the F 5′ phosphate. The tagged G-rich strand therefore reads telomere repeats followed by the F oligo, which matches the read analysis (capture oligo at the 3′ end after the telomeric region, line 232).

Model correction (2026-10-17): the earlier model drew F01 and the phase arm as one strand ligated to a "C-rich 3′ end". That strand is not printed in the source and failed the source guard; it was replaced by the F01/R01-1 duplex.

## 4. Reaction order

1. 🟢 About 20 µg gDNA, 5 µl 10 µM telobaits, ATP, T4 DNA ligase and FastDigest RsaI/HinfI incubate together overnight at 37 °C (Methods line 224).
2. 🟢 MyOne Streptavidin T1 beads bind the ligation mix 2 h at room temperature, then three washes.
3. 🟢 On the beads, T4 DNA polymerase with dNTP fills the gap (37 °C, 30 min); the model extends the R 3′ end to the C-rich 5′ end and leaves a nick (🟡 the article does not mention nick sealing).
4. 🟢 EcoRI releases the fragments; 🟡 the cut (G^AATTC) leaves the biotinylated stub on the beads and a 5′-AATT overhang on the released fragment.
5. 🟢 0.45× AMPure PB purification.
6. 🟢 SMRTbell Express Template Prep Kit 2.0 "as per manufacturer's recommendation"; 🟡 end repair (drawn as filling the AATT overhang), dA-tailing and hairpin ligation are inferred from the kit name.

Published discrepancy: Results line 128 says T4 polymerase gap filling follows ligation and fragments are "subsequently purified" on streptavidin beads; Methods line 224 adds beads to the ligation mix and gap-fills on the beads. The model follows the Methods order because it gives the volumes, buffers and washes for each on-bead step. The page shows this notice.

## 5. Final library and sequencing

A PacBio primer/polymerase complex enters at the hairpin primer site and repeatedly traverses the selected telomere fragment and its reverse complement. Reads are demultiplexed by the telobait barcode within the last 150 bases (lines 230–232). The article does not print the hairpin or sequencing-primer bases.

## 6. Open questions

- SMRTbell kit end-repair/hairpin geometry and primer site (kit guide not read).
- Sequencer and chemistry version are not given in the inspected Methods passages.
- Whether the gap-fill nick is sealed before SMRTbell preparation.
