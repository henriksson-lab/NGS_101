# Assembly, transposition and hybrid cleavage animation review

1 November 2026. Focused source reading for three selected mechanism animations;
this is not a whole-protocol completeness or missing-information audit. Original
files remain in `/husky/henriksson/NGS101_papers/sources/`. No experimental yield,
rate, processivity or population-selection probabilities are inferred.

## Gibson assembly

Source-owned example: `gibson-assembly__10.1038+nmeth.1318/tools/gibson_assembly.py`,
`example_fragments()` and `precursor()`. All base values and the exact six-base
resection boundaries are explicitly invented in that owner. The animation follows
one joint, including both parental duplex bodies and their outer cohesive ends;
it does not claim completed circular vector construction.

Inspected original: published supplementary file
`nmeth.1318_41592_2009_BFnmeth1318_MOESM319_ESM.pdf.txt`, supplementary Figure5 caption
lines498–504 and supplementary Results lines815–925. The caption distinguishes the
T5-exonuclease Figure1 method from the T4-polymerase/ExoIII alternatives. The Results
explain why gap repair and nick sealing are distinct and why incomplete repair
changes assembly recovery. The three-prime alternatives are not used as authority
for five-prime processing or named enzymes in this movie.

The defining main article is a manual manifest row; the direct Nature article route
redirected to a login endpoint in this pass. No new full-main reading is claimed.
The existing source-qualified owner supplies the generic five-prime-exonuclease,
polymerase and ligase mechanism. We deliberately do not assign commercial enzyme
identity or mixture conditions from another supplementary method. Source bundle
primer scan rerun: `/tmp/assembly-gibson-primers.log` (22 candidate hits; pointers,
not evidence). No source sequence was transcribed into the synthetic example.

The exact-end limiting example has no gap: polymerase stops on occupied template
without adding atoms. Showing an arbitrary gap-filling animation would be wrong.
Four terminal resection operations release24 actual nucleotide IDs. All64 retained
IDs are represented by the two32nt joined strands after sealing. Pink highlights
are gated on the actual simulated covalent bonds. This is a selected conceptual
example with openly bounded main-source access, not a newly complete source review.

## ATAC-seq

Inspected published originals:

- `atac-seq__10.1038+nmeth.2688/nmeth.2688_PMC3959825.xml.txt`, Online Methods
  lines73–79:37°C30min transposition, MinElute purification, NEBNext PCR mix and
 72°C5min before denaturation. Primary data processing line101 explicitly describes
  a Tn5 dimer and9bp-separated adapters. The cached XML is the published article's
  author manuscript; the defining DOI and protocol identity are unchanged.
- Its Springer supplementary original
  `nmeth.2688_41592_2013_BFnmeth2688_MOESM200_ESM.pdf.txt`, Supplementary Table1
  lines1–99, especially Ad1_noMX and Ad2.1 lines3/7. These primers support existing
  owner components; they are not a disclosure of each native loaded transposome
  strand. Primer scanner rerun `/tmp/assembly-atac-primers.log` (two collapsed hits),
  followed back to those original table rows.

Generic Tn5 choreography and polarity were already source-reviewed against
Kia2017 and Adey2010 in `animation/tn5-mechanism-research.md`; this pass does not
claim a new reading of those originals. The source supports workup then extension,
not spontaneous disappearance of the protein. The local recipe supplies canonical
Nextera transferred ends and separate ME′ complements with unassigned terminal
chemistry. Target DNA and cut coordinates are synthetic. No nucleosome occupancy,
chromatin-opening transition or native transposome loading stoichiometry is inferred.

Both strand transfers in each homodimer are atomic; a9bp stagger is checked against
the actual target duplex. Each target3′ end copies nine bases, stalls on the bound
ME′ oligo, and continues only after an explicitly conditional contact-release event.
Displacement versus digestion of ME′ and the precise polymerase identity in the
source commercial mixture remain unassigned. No ligase or final PCR library is
invented. The movie ends at a complemented pre-PCR product. Nonselected target
fragments and short oligos retain physical tracks at their preceding offscreen
positions; every active atom has a track in every shot.

## RNase H1

Inspected primary [Lee et al.,2022](https://doi.org/10.1093/nar/gkab1064), Introduction
and Discussion: RNA-selective cleavage in RNA–DNA hybrids; four consecutive RNA
bases around an eligible cut; substrate geometry affects degradation mode. The
movie uses selected endocuts in a blunt synthetic hybrid. It does not transfer the
paper's overhang-dependent processive mechanism to that input.

The source fetcher retrieved the main and supplement via PMC Cloud under
`doi_10.1093+nar+gkab1064/` with manifest records. Main PDF text
`nar_gkab1064_PMC8886854.pdf.txt` lines70–87 and1714–1729 were checked. Primer scan
rerun `/tmp/assembly-rnaseh-primers.log` (19 candidate hits); supplementary original
`nar_gkab1064_supp_gkab1064_supplemental_file.pdf.txt` TableS1 lines150–183 was read
for substrate identities/modifications. Those experimental fluorophore-labelled
substrates are not substituted for our explicitly synthetic sequence. Their other
specific cleavage histories and supplementary figures are outside this reading.

[IUBMB EC3.1.26.4](https://iubmb.qmul.ac.uk/enzyme/EC3/1/26/4.html), reaction/comments,
was checked for RNA–DNA-hybrid cleavage producing5′-phosphomonoester products.
The model retains hydrolytic3′OH/5′P ends, RNA sugars and the complete unchanged DNA.
First-fragment dissociation is a separate illustrative unpairing operation; cutting
alone leaves hybrid contacts untouched. No RNA base vanishes, becomes DNA or is
silently assigned to a new primer. H2 specificity and commercial CEL-Seq vial
identity are outside this example.

## Verification and visual review

`python3 -m animation.selftest_assembly_transposition`:8 tests pass. They perturb
ATAC target sequence, reject a terminal cut, verify actual transfer/copy/cleavage
relations, distinguish cleavage from dissociation, account for every active atom,
and check both molecular equality and nucleotide-coordinate continuity at view
boundaries. Factories roundtrip through the JSON schema.

Native960×540 stills at establish/action-middle/action-end were rendered and
inspected for every scene under `/tmp/assembly-transposition-review/`. No text
layout overflow was reported. Reviewed features include four accessible exonuclease
ends, overlap docking, both ligation seams, paired Tn5 dimers and short oligos,
recessed target3′ repair, final pre-PCR anatomy and retained RNA cleavage fragments.
Final1080p MP4 auditing is performed by the coordinating worker.
