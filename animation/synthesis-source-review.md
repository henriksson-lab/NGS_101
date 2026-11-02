# Source review for the synthesis animation examples

Selected mechanisms, reviewed 1 November 2026. These are not complete assay movies
or accepted missing-information audits. Molecular inputs and design choices are
explicit; native yields, whole-transcript sequences and reaction kinetics are not
predicted. Sources are external under `/husky/henriksson/NGS101_papers/sources/`.
The complete bundle primer scraper was run separately for each of the three slugs;
its hits were used as pointers, then the original text was read. No upstream drawn
page or superseded preprint supplied a chemistry claim.

## CEL-Seq

Original: published supplement `j.celrep.2012.08.003_mmc4.pdf`, printed CEL-Seq
Protocol revision 2012-07-16. Text lines 38–46 describe the anchored poly(T),
eight-base barcode, small-RNA adapter and T7 promoter, and give primer #1 verbatim.
The movie selects that source barcode **CATCACGC**, selects V=A, and appends a short
explicitly illustrative cDNA insert. The prepared complementary strand is supplied
as an input, not fabricated by a claimed second-strand enzyme. The source primer
components come from the protocol owner's `native.py` / `celseq.py` definitions.

Inspected protocol first-/second-strand and IVT sections: the source proceeds from
second-strand cDNA preparation to T7 RNA amplification. The movie begins at the
prepared duplex, shows de novo RNA synthesis, and retains both DNA strands and the
first RNA when a second transcript is synthesized. Initiation includes the printed
promoter terminal GGG; it is not shifted to the first downstream adapter base.
RNA identity features are copied from the explicitly annotated DNA templates.
The icon and RNA departure are a schematic elongation view, not a structural model
of the T7 transcription bubble. Fragmentation and downstream adapter ligation are
outside the selected mechanism.

## Drop-seq

Originals: published Cell 2015 Supplemental Information
`j.cell.2015.05.002_1-s2.0-S0092867415005498-mmc1.pdf`, text lines 119–132
(reverse transcription and exonuclease treatment), and published Table S6 in
`j.cell.2015.05.002_1-s2.0-S0092867415005498-mmc7.pdf`, text lines 2891–2914.
Table S6 prints the bead-linked SeqB oligo, seven-base T linker, SMART handle,
shared AC, twelve random barcode bases, eight random UMI bases, and dT30.
The original Methods supplies Maxima H− RT, dNTPs and the template-switch oligo.

The movie uses the owner module's **SeqB** scaffold and explicitly synthetic
barcode/UMI members. An illustrative RNA body ending C and a 30-base poly(A) window
make the selected binding alignment unambiguous. Only the dT30 portion binds RNA;
the DNA identifiers and bead linkage are retained during actual RNA-templated
extension. The source supplies the template-switch oligo during RT, but this
focused clip stops before terminal addition/switching; it does not depict a complete
RT vessel or imply that the source omits the TSO. Cleanup and later PCR are outside
this selected view. No RNA is silently deleted.

## scRRBS

Original: published defining article, `PMC3847781.html`, text lines 1225–1227,
“Construction of single-cell RRBS sequencing libraries”. Read its complete methods
paragraph: MspI C^CGG digestion; Klenow exo− fill plus 3′A; premethylated indexed
adapter ligation; bisulfite denaturation/conversion; on-column desulfonation/recovery;
first-round PfuTurbo Cx then second-round Phusion PCR.

The animation begins after the preparative steps with one explicitly synthetic,
denatured genomic window. Susceptible C becomes physical DNA U; selected 5mC
remains. PfuTurbo Cx copies U to A, and another complementary copying event makes
T. These are three separately retained molecular identities. Both illustrated
cycles use PfuTurbo Cx within the first PCR round; the second drawn copy is not
labelled as the source's later Phusion round. Original DNA and previous daughter
remain after denaturation. Short synthetic genomic primers are introduced intact
only after conversion/recovery, then arrive from outside the viewport and anneal.
They are visibly scoped as illustrative priming windows, **not** claimed source PCR
oligos. The movie does not invent undisclosed adapter geometry or full native
primer assignments. These preparative boundaries are stated in the initial scope
and provenance rather than disguised as simulated steps.

## Checks and artifacts

`python3 -m animation.selftest_synthesis_examples` exercises causal dependencies:
paired T7 promoter, actual insert/identifier copying, distinct transcript identities,
RT template-sugar capability and physical primer contacts, retained bead primer
atoms, protected methyl-C, U/A/T atom identity, atomic stock-primer arrival and
identical molecular states across every incoming transition. Nine tests pass, including explicit track coverage for every active nucleotide.

960×540, 12fps H.264 preview movies and contact sheets are in
`/tmp/animation-synthesis-previews/`. Each scene was inspected at reaction midpoints
and endpoints. Enzyme labels were moved away from molecule captions after the
first scRRBS contact review. Final 1080p export and decoded-frame audit are owned
by the parent task.
