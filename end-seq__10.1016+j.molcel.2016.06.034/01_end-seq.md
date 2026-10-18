# END-seq — research notes

## Sources read

- 🟢 Canela et al. 2016, DOI [10.1016/j.molcel.2016.06.034](https://doi.org/10.1016/j.molcel.2016.06.034), PMC6299834, main text and Experimental Procedures.

## Molecular order

🟢 Cells are embedded in agarose before deproteinization. Native DNA ends are blunted,
dA-tailed and ligated in-plug to a biotinylated, USER-openable hairpin bearing the P5
side. DNA is recovered, sonicated to 150–200 bp and captured on streptavidin. The new
sonication end is repaired, dA-tailed and ligated to the P7 hairpin. USER opens both
hairpins; TruSeq index PCR makes the sequencing library. Read 1 begins at the first base
of the blunted DSB end.

## Oligos

🟢 The Methods print both ordered hairpins, including 5-prime phosphate, the two
biotin-dT residues in adaptor 1 and terminal phosphorothioate bonds. The Python model
stores the base strings and biotin-dT coordinates; the page keeps the chemical
modifications in the research boundary rather than pretending they are ordinary bases.

The source-printed adaptor 1 is
`5′-Phos-GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGUU[Biotin-dT]U[Biotin-dT]UUACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′`;
adaptor 2 is
`5′-Phos-GATCGGAAGAGCACACGTCUUUUUUUUAGACGTGTGCTCTTCCGATC*T-3′`.
Here `U` is the USER-cleavable deoxyuridine and `*` marks the terminal
phosphorothioate bond. Both have a 3′-T overhang for ligation to dA-tailed DNA.
The source names TruSeq index adapters for amplification but does not identify
the particular index oligo used in each library. The paper's control input
libraries use universal adapters in a separate branch; these do not replace the
END-seq hairpins.

## Interpretation boundary

🟢 End repair means END-seq maps the processed boundary but does not retain the original
overhang sequence or chemistry. The final model now retains source-written USER arms and leaves PCR-added outer
ends and index architecture unassigned; it does not impose a canonical shell.

## Bounded primary-source check (2026-10-16)

Inspected: Original PMC6299834 JATS152–164, especially159,161.

Added source-written modified hairpins, GELase plug recovery, 16-cycle END-seq PCR and alternative 50/75-base single-end runs. Canonical outer library is inferred; index acquisition unassigned. USER fragments remain unresolved rather than presumed discarded.

Deferred: Original NIHMS998981 Supplement1 PDF and Supplement2–6 tables; cited BLESS method; exact index assignment. This is a capped check, not a full source audit or confirmation of missing information.

## Source-owned validation migration (2026-10-18)

One original was re-read: defining PMC6299834 JATS Results61 and Experimental
Procedures159–161 (source hairpins/modifications, recovery, distal enzymes,
USER-before-PCR,16cycles and alternative50/75base single-end runs). Mandatory
bundle-scraper hits were verified directly. The adapter1 bracketed `[Biotin-dT]`
notation splits the scraper’s whole-strand extraction; a documented exemption is
paired with direct full modified-sequence and biotin-position transcription tests.
No supplement, Figure1A original, exact index-primer sheet or further vendor
revision was followed in this bounded pass.

The local endpoint now derives the retained adaptor1 Read1 arm and adaptor2 arm
from the actual source hairpins after USER, with PCR-added outer/index ends
unassigned. No generic dual-index TruSeq shell, eight-position index defaults or
arbitrary42-position insert remains. The Read1 reference exactly matches the
retained source adaptor1 arm. Actual physical sequencing reagent, index architecture,
member/width and acquired index cycles remain unassigned;50base HiSeq2500 and75base
NextSeq500 are source alternative runs.

Both modified hairpin inputs and their known folded stem relations are explicit.
Each ligation geometry uses the actual dA-prepared preceding material. PCR uses
actual USER-opened retained-arm DNA, with its complete complementary ends and
physical PCR primers unresolved, rather than the expected final product. Distal
T4DNApolymerase/Klenow/PNK, Klenow-exo-minus tailing, QuickLigation,150–200bp Covaris
and MyOneC1 enrichment labels follow the original. The source-first in-plug
blunting/tailing enzyme identities remain deferred to the original Figure1A.
USER loop products retain an unresolved physical fate; cleavage is not a discard.

Legacy fallback workflow and source-selecting `end_capture` drawing shortcuts were
removed; the explicit-parameter `UserHairpin` record remains a neutral container.
This pass does not accept absence claims or constitute a full chemistry review.
