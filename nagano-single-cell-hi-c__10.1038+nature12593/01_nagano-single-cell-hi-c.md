# Nagano single-cell Hi-C

## 1. What it is

The foundational single-cell Hi-C implementation: proximity ligation is performed in a population of fixed nuclei, after which individual nuclei are isolated and their contact libraries are prepared separately.

## 2. Sources and evidence

- 🟢 Nagano et al., *Nature* 2013, DOI [10.1038/nature12593](https://doi.org/10.1038/nature12593), defining paper.
- 🟢 Nagano et al., *Nature Protocols* 2015, DOI [10.1038/nprot.2015.127](https://doi.org/10.1038/nprot.2015.127), detailed implementation used to resolve BglII, biotin-14-dATP, streptavidin capture and PCR ordering.
- 🟢 The [2013 supplementary methods](https://www.nature.com/nature/journal/v502/n7469/extref/nature12593-s1.pdf),
  pp. 1–5, give the AluI digestion, dA tailing and adapter-ligation order, then print twelve customized Illumina adapters, each with a 3-bp identification
  tag, and both PCR primers. The schematic uses the printed CAA adapter member and the
  historical paired-end primer sites; there is no invented i7 read.

The ordered CAA adapter is 🟢 `5′-ACACTCTTTCCCTACACGACGCTCTTCCGATCTCAA*T-3′`
paired with `5′-pTTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG-3′`. Here `p` is a 5′
phosphate and `*` is a terminal phosphorothioate linkage, as defined by the supplement.
These modifications belong to the ligation oligos; the final PCR product is drawn as
ordinary DNA. The two printed library-amplification primers also end with a terminal
phosphorothioate linkage. The paper reports 2 × 40 bp reads beginning with the inline
identification tag, with separate index acquisition unassigned.

## 3. Distinguishing chemistry

🟢 BglII digestion and biotinylated fill-in precede proximity ligation. 🟢 Single nuclei are separated after the bulk ligation reaction. 🟢 After crosslink reversal, streptavidin beads select marked contact junctions; AluI then cuts the bead-bound DNA, Klenow fragment (3′→5′ exo−) adds a 3′ dA, and the printed adapter's 3′ dT overhang is ligated before PCR. The published adapter pairs over a 17-bp core, leaving the two sequencing arms unpaired.

## 4. Bounded source pass — 2026-10-17

Read original published `nature12593-s1.pdf.txt` identity/DOI and
“Preparation of Hi-C library from single cells”15–65; library validation
context67–80; all twelve modified adapter-pair rows88–145; library PCR
primer rows153–158. The manifest gives the original publisher PDF URL and
identity. Ran the source-bundle oligo scraper and checked modifications and
wrapped primer bases at the original-text pointers. The separate cached
article text has no accompanying original in this bundle; it was not used
as a newly inspected defining article. The 2015 protocol was not read or
used to substitute conditions in this pass.

🟢 Source33–35 specifies picking single nuclei under a microscope with
mouth-controlled glass capillaries, followed by overnight 65 C crosslink
reversal. The graph's generic “dilution and sorting” is now corrected to
that operation. Other nuclei have an unresolved fate outside the picked
lineage, rather than an assumed physical discard.

🟢 Fill-in uses Klenow large fragment, proximity and adapter ligation use
T4 DNA ligase, and library amplification uses Platinum Pfx for 25 cycles.
The two printed library primers, each with a terminal phosphorothioate
linkage, now have source-derived graph states. PCR is followed by removal
of streptavidin beads, AMPure XP purification and 300–700 bp recovery from a
2% agarose gel. The graph now includes these recovery boundaries. Both
source-stated 40-base reads are declared, with separate index acquisition unassigned.

🟡 The local CAA adapter and historical paired-end shell are unchanged.
Known adapter geometry does not validate the heterogeneous pre-PCR contact
population; that PCR relation remains unresolved.

Deferred: full defining article original/provenance; remaining source figures,
linked predecessor/vendor originals and detailed 2015 variant; wash states,
precise bead-removal fate and complete AluI fragment population. BglII is
the selected route; the supplement's “DpnII for some libraries” is not
silently substituted. This capped published-supplement pass does not establish
full experimental validation or any accepted absence claim.


## Bounded native migration — 2026-10-18 batch08

Inspected 1 verified originals: Verified publisher nature12593-s1.pdf DOI/title and original Methods15–65, modifiedCAA/other adapter families88–145 and PCR primer rows153–158. Source MANIFEST verifies this original. Cached main-articleTXT53–61 inspected for context only but no corresponding original/manifest identity, so not counted or promoted as fresh defining-source evidence. Mandatory bundle scraper reviewed.

Removed shared ContactWorkflow/historical endpoint constructor; all localreaction choices from independently read2013SI. Explicit5′GATC BglII geometry/fill/biotin contact, postligationpicking, capture/AluI/dA/nativeYligation washes and actual forked source-strand PCR precursor. Derived daughter endpoint includes source dT/insert dA junctions omitted by legacyendpoint. CAA and opposite tag features retained; original biotin/methylation not copied and residual templates remain unresolved. Physical wash/fraction removal distinct from contact QC. Source2×40 positive, separate indexacquisition unassigned.

Deferred: Full defining article original provenance and remaining sourcefigures; unmanifested cachedmainTXT not fresh evidence; Detailed2015revision, linked predecessor/vendor originals and complete native contact/AluIpopulation topology; Physical removal/recovery of original biotin templates and precise first-cycle fork alignment; DpnII-for-some-libraries not assigned to selectedBglIIroute. This source-capped pass is not a complete protocol review or an acceptance of source absence. Historical audit paragraphs above retain their original scope; this entry states the current model.

Independent implementation review retained source terminal phosphorothioate bonds on the physically incorporated forward/reverse PCR-primer footprints. Original adapter backbone modifications and biotin are not copied into daughter DNA. Native Y-ligated originals retain their own source A*T bond; this is distinct from the C*T bond supplied by a PCR primer.
