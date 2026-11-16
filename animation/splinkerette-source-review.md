# Splinkerette PCR: source and animation review

The protocol is Potter and Luo, PLoS ONE5:e10168 (2010),
[doi:10.1371/journal.pone.0010168](https://doi.org/10.1371/journal.pone.0010168).
The defining published article and Protocol S1 are the chemistry authorities.
See the protocol's [source notes](../splinkerette-pcr__10.1371+journal.pone.0010168/01_splinkerette.md)
and external manifest for exact inspected sections and original identities.

## Scope

The adapter oligos and the three P-element3′ primers are source transcriptions.
Their short locus, inter-site distances and genomic flank are explicitly synthetic.
The selected path demonstrates a restriction-derived5′ GATC end, a folded TOP
oligo with a23-bp BOT core and7-bp terminal stem, one-strand ligation, selective
first copying, S1 priming, two nested template-copying cycles and endpoint
reassociation. The final molecule is a flanking-sequence amplicon for Sanger
sequencing, not an NGS library. Cleanup and Sanger chemistry are described on
the page but are outside the movie.

The hairpin is a checked complementary geometry, not a folding-energy result.
It retains a3′ OH; the model does not invent a blocked terminus or promise zero
background. Spectator restriction fragments, the original TOP, parental templates
and earlier PCR products remain in the molecular inventory when moved offscreen.
No yield, PCR population, melting competition or kinetic prediction is claimed.

## Shared molecular trajectory

`tools/splinkerette.py` owns the source-qualified reaction trajectory. Both the
graph's sequence maps and the movie consume its snapshots and nucleotide IDs.
The previously animation-local standard-library engine now lives in
`lib/nucleotide_simulator.py`; the historical import remains an identity-preserving
compatibility module. Declared same-strand folding is checked without changing
any covalent bond. Source-selected ligation requires the actual3′ OH/5′ P nick
and ATP. The opposite genomic/TOP junction remains unsealed.

Extension walks actual physical templates. It creates the S1 complement rather
than inserting a finished PCR product. Both nested primer complexes extend in
each demonstrated cycle. The first nested products retain one outer boundary;
the following cycle establishes the bounded complementary pair.

## Presentation and verification

The reviewed export is
`/tmp/chemistry-animation-splinkerette-complete/splinkerette.mp4`:
1920×1080,30fps,H.264/yuv420p,128.3s,default1.5× playback. The timeline footer
is off. Its adjacent manifest records source scope, exact scene declarations and
the full molecular reaction trace. Generated videos and decoded images remain
outside the repository.

Visual review used decoded MP4 action and transition frames, plus full-size
native frames. Revision removed duplicate overlapping BstYI actors: one selected
restriction event cuts both backbones. Persistent colours identify the known
insertion end and unknown flank. Primers enter from outside the viewport;
polymerases appear after binding. Both nested complexes have their own polymerase.
The hairpin and core are visible, the new ligation bond stays highlighted, the
opposite nick stays open, and the final anatomy identifies both boundaries,
the T3 direction, spacers, genomic flank, restriction junction and copied core.

Automated checks cover unchanged molecular state and conserved nucleotide
coordinates at every shot transition,24-pixel separation of final paired bases,
explicit reagent births, correct polymerase timing, the persistent opposite nick
and native text overflow. The shared simulator has standard-library tests;
protocol/source transcription checks and the full site/link build also pass.
The audit tool compares the exported manifest with current definitions, verifies
the encoded stream and hashes the shared engine along with the presentation code.
These checks do not replace source review or establish biochemical efficiency.
