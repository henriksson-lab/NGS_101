# Telomere-end tagging

Telomere sequencing methods must distinguish the physical chromosome terminus from an
internal fragment containing telomeric repeats. Most targeted long-read workflows solve
this by reacting with the native G-rich 3′ overhang before fragmentation or enrichment.

## Repeat phase

A human telomere can terminate at any of the six offsets of `TTAGGG`. An adapter with a
telomere-complement arm therefore needs six cyclic phases. `lib/telomere.py` derives that
family from one repeat; protocol modules select published arms rather than maintaining six
independent hand-written identities.

## Recurring workflow states

- A native end is duplex telomere followed by a G-rich 3′ overhang.
- A phased terminal adapter anneals across the overhang; terminal ligation is marked at
  the native C-rich 3′ end.
- Polymerase fill may complete the adapter–chromosome duplex.
- Biotinylated tags can select chromosome ends on streptavidin. Some workflows then cut a
  restriction site encoded in the tag to release the selected molecule.
- Oxford Nanopore libraries finish at a motor-loaded adapter and have no sequencing
  primer. PacBio libraries finish as closed SMRTbells with primer sites in the hairpins.

The six protocol pages keep their paper-specific order and oligos; this concept supplies
only the molecular invariants shared between them.
