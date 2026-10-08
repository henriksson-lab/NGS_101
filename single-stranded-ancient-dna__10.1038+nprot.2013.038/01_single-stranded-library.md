# Single-stranded library preparation for ancient or damaged DNA

## 1. What it is

This protocol recovers individual highly degraded DNA strands instead of requiring an intact duplex with two usable ends.

## 2. Sources and evidence

- 🟢 Gansauge MT, Meyer M. *Single-stranded DNA library preparation for the sequencing of ancient or damaged DNA.* Nature Protocols 8:737–748 (2013). [doi:10.1038/nprot.2013.038](https://doi.org/10.1038/nprot.2013.038).
- 🟢 A phosphorylated, biotinylated adapter is ligated to the 3′ ends of heat-denatured DNA, captured on streptavidin beads, copied, and joined to a second adapter by blunt ligation.
- 🟢 Published oligos include CL78 with internal C3 spacers and terminal TEG-biotin, the CL53/CL73 second adapter, CL9 extension primer and CL72 custom Read 1 primer.
- 🟢 CL72 is `ACACTCTTTCCCTACACGACGCTCTTCC`; the custom primer is required because the P5-side adapter is truncated relative to ordinary TruSeq libraries.

## 3. Molecular path

Every denatured template strand is independently captured through its 3′ end. Copying converts it to a duplex only after capture, avoiding the need for complementary ancient strands to survive together.
