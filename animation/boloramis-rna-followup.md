# BOLORAMIS RNA retention during RCA

Focused source follow-up, 31 October 2026. Bundle scraper rerun. This is a scoped
positive-evidence check, not an absence audit or complete reaction validation.

## Original BOLORAMIS recipe

The published article, [doi:10.1093/nar/gkab120](https://doi.org/10.1093/nar/gkab120),
was read in cached `nar_gkab120_PMC8191787.pdf.txt`, lines160–183 and214–221.
The non-384-well route transitions from SplintR ligation through hybridization-buffer
rinses into NxGen Phi29 RCA at37°C. The organoid route rinses into Phi29 buffer.
The published SI `BOLORAMIS_NAR_template_Supplementary_Materials_2021-02-08.pdf.txt`,
“High-content384-well BOLORAMIS library preparation”, lines85–96, adds RCA mix
immediately after ligation. None of these inspected transitions introduces an RNase
operation. This does not establish publication-wide absence of every RNA-related
operation; it establishes that an animation must not invent an RNase between these
specified source steps or simply omit the RNA molecule.

## RNA displacement is a distinct capability

The primary single-molecule study
[Kim et al., *Direct visualization of replication and R-loop collision using
single-molecule imaging*, doi:10.1093/nar/gkad1101](https://doi.org/10.1093/nar/gkad1101)
provides positive evidence beyond DNA strand displacement. The original was fetched
by the repository source tool into
`/husky/henriksson/NGS101_papers/sources/doi_10.1093+nar+gkad1101/`; its bundle scraper
was run. Main PDF text lines119–143 specify30°C pre-incubation and25°C curtain experiments. Results
“RNA remains on the R-loop when the replication stalls”, lines407–440, and Figure3G–H
compare labeled RNA after collisions: passing template-strand hybrids lose much more
RNA than stalled molecules. The authors interpret passage as RNA displacement or
hybrid unwinding. Their experiments also show stalling and retained RNA.

This establishes a Phi29 capability on experimentally specified RNA–DNA hybrids,
not guaranteed passage for every RNA target, a quantitative BOLORAMIS yield, or an
independent measurement of the exact Lucigen NxGen preparation used in BOLORAMIS.
The displayed short RNA in our ACTN1 example is already a conditional match, not a
verified full transcript. Its geometry cannot silently stand for complete cellular RNA.

## Simulator consequence

Retain RNA identity and hybrid bonds after circularization. A selected successful
RCA branch may explicitly allow RNA-strand displacement, with the above mechanistic
support and BOLORAMIS applicability marked as an inference. At each collision,
remove RNA–DNA bonds as the new daughter pairs; retain the RNA molecule rather than
inventing cleavage. Do not derive this permission solely from DNA-displacement
activity. If that inference is outside the selected validation scope, keep RNA
unpairing an explicit conditional boundary and report it unresolved.

Phi29's separately reported3′→5′ exoribonuclease activity is not evidence that the
BOLORAMIS recipe digests its RNA splint. No RNase or RNA degradation operation is
added by this source pass.

## Animation source check and staging correction

The new movie was checked directly against these three originals on31October2026:

| Original and external text twin | Inspected location and consequence |
|---|---|
| Published BOLORAMIS paper, `nar_gkab120_PMC8191787.pdf.txt` | Printedp3, non-384-well preparation, text168–180: padlocks and RCA primer enter together,1:1molar ratio; hybridization/wash precedes SplintR; retained complexes enter NxGen Phi29 RCA with dNTPs and aminoallyl-dUTP. |
| Published Supplementary TableS3, `nar_gkab120_supp_gkab120_supplemental_files_Supplementary_Tables_Final.xlsx.txt` | Text365–373: source23nt mRNA RCA primer and source60nt5′-phosphorylated ACTN1A1 probe. Exact workbook cells are also covered by the protocol suite. The bundle primer scraper was rerun without the known-only filter, and the relevant source rows were read directly. |
| Published Phi29 collision study, `nar_gkad1101_PMC10783495.pdf.txt` | Text119–143 identifies the25°C curtain assay;407–440 distinguishes retained RNA on stalled trajectories from displacement on successful template-strand passages. This does not establish identical behavior of BOLORAMIS's NxGen37°C formulation. |

The earlier animation recipe added its RCA primer after sealing the probe. That was
a modelling/staging error, not a publication error. It has been retired. The current
`animation/boloramis.py` co-hybridizes the primer and retains its nucleotide identities
through folding, ligation and RCA. The public adapter and CLI use that same recipe.
The DNA nucleotide pool and its unresolved aminoallyl-dU composition are supplied
only at the post-ligation RCA stage, matching the inspected source order.

The RCA model selects successful RNA-hybrid passage. Disabling that capability stalls
copying at the occupied RNA hybrid and preserves its contacts; there is no desired-product
fallback. Each new DNA T-equivalent carries a structured unresolved dT/aminoallyl-dU
occupancy annotation. Canonical backbone lettering does not assign modification sites
or crosslinks. The RNA remains a conditional25nt matching window, not a newly verified
full ACTN1 transcript. These are positive source checks and explicit modelling limits,
not missing-information absence claims.
