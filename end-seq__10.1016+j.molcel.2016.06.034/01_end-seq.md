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
biotin-dU residues in adaptor 1 and terminal phosphorothioate bonds. The Python model
stores the base strings and biotin-dU coordinates; the page keeps the chemical
modifications in the research boundary rather than pretending they are ordinary bases.

## Interpretation boundary

🟢 End repair means END-seq maps the processed boundary but does not retain the original
overhang sequence or chemistry. 🟡 The final PCR-completed library uses the canonical
TruSeq arm model; variable index bases are structured placeholders.
