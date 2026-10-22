# MPRA: pooled reporter construction and Tag-Seq

Melnikov A et al. (2012).
[Systematic dissection and optimization of inducible enhancers in human cells using a massively parallel reporter assay](https://doi.org/10.1038/nbt.2137).
Nature Biotechnology 30:271–277.

## Source scope

Bounded readings, 2026-10-21; a scoped construction and Tag-Seq model is implemented.
Read the defining published article's author manuscript, PMC3297981, own DOI
verified, Methods “Oligonucleotide library design and synthesis”, “Plasmid
construction” and “Tag-Seq”: local JATS text lines 152–166 and 176–182.
The manuscript title differs slightly from the publisher title; it identifies the
same published DOI, not a separate preprint chemistry.

Ran `tools/scrape_primers.py` across the retrieved bundle again. Article Methods,
rather than the variant-table hits, establish the constructions below. Read the
publisher's `MOESM10` Supplementary Figure 1 in the original rendered PDF page 2,
including strand ends and the final library orientation. The NIHMS supplement's
corresponding extracted Figure 1 was consulted for source locators. The publisher
original is the authority for the exact drawn boundaries. Workbook headings identify
the available design tables; their individual design/tag assignments remain unselected
and are not treated as checked molecular facts in this pass.
This is not a complete source review or a confirmed-missing-information audit.

## What is new

🟢 Each designed regulatory sequence is linked to an identifying tag before insertion
of the reporter cassette. The model must preserve that association through both cloning
steps and then through distinct plasmid-DNA and RNA-derived sequencing branches.
The tag is a reporter identifier, not a UMI or a sequencing sample index.

## Published molecular path

🟢 The initial synthesis oligo has this ordered structure, from Methods line 154:

| Region | Source value |
|---|---|
| Universal forward site | {{= mpra.UNIVERSAL_F }} |
| Regulatory variant | {{= mpra.VARIANT_LENGTH }} nucleotides, experiment-specific |
| KpnI / XbaI interval | {{= mpra.KPNI + mpra.XBAI }} |
| Identifying tag | {{= mpra.TAG_LENGTH }} nucleotides, design-assigned |
| Universal downstream site | {{= chemdraw.revcomp(mpra.UNIVERSAL_R) }} |

The downstream site is written on the synthesis strand. Its amplification primer is
the reverse complement, `{{= mpra.UNIVERSAL_R }}`; do not assign both strings the same
physical primer orientation.

🟢 Pool amplification and re-amplification add SfiI cloning sites. SfiI cloning
transfers the variable cassette to the source-specific pGL4.10M backbone.
KpnI/XbaI then opens the interval between regulatory variant and tag, and insertion
of the minimal-promoter/luc2 fragment places the tag downstream in the reporter
transcript. A subsequent KpnI-based background-removal step is a separate physical
operation, not a computational read filter.

🟢 Tag-Seq has two inputs: plasmid DNA, or poly(A)-selected, DNase-treated RNA
converted to cDNA. Indexed PCR generates sequencing libraries from either input.
The article reports 36-nucleotide single-end HiSeq2000 sequencing; the first ten
read bases identify the designed reporter tag.

## Source primers used by the model

🟢 Methods line 160 prints the SfiI-adding primers, written 5′→3′:

- `{{= mpra.SFI_F }}`
- `{{= mpra.SFI_R }}`

🟢 Methods line 180 prints the Tag-Seq primers, written 5′→3′:

- `{{= mpra.TAG_F }}`
- `{{= mpra.TAG_R }}`

The paper writes its sample-index placeholder as {{= mpra.SAMPLE_INDEX_LENGTH }} `X` positions;
the model uses the same number of `N` positions. Keep them distinct
from the reporter's ten-base identifying tag. During implementation, import canonical
Illumina components for the matching fixed regions and retain local source justification
for their use here. The forward Tag-Seq oligo is the canonical P5/TruSeq Read1
primer with the shared ACAC region counted once.

## Construction and validation scope

The synthesis architecture has {{= mpra.synthesis_length() }} nucleotides when the
declared variant and tag lengths are counted. The generic view selects no particular
published design/tag pair. It uses the source's design-assigned identifier semantics,
without asserting an observed barcode association or a measured reporter activity.

🟢 The original Supplementary Figure 1 shows SfiI cuts leaving AAC and CGG on
opposite 3′ ends. The explicit input flanks are `{{= mpra.SFI_LEFT }}` and
`{{= mpra.SFI_RIGHT }}`. The selected insert core is sliced from the actual modelled
SfiI-adding PCR product, preserving the variant and tag; it is not retyped as a
separate endpoint sequence. The reverse SfiI-adding primer has an 11-base 3′ match
to the original downstream site. That source-specific short match is explicit.

🟢 Reporter installation uses KpnI GTAC 3′ and XbaI CTAG 5′ ends. Their mixed
polarity is part of the molecular declarations. Matching letters with opposite
polarities are rejected by the neutral cohesive constructor and independently
checked against the actual input states. Each fusion is retained once. These checks
do not prove a complete plasmid digest or establish frequencies of other products.

🟢 The pGL4.10M preparation remains separate: HindIII/XbaI removal of luc2,
Klenow filling at HindIII, Mung Bean removal of the XbaI extension, and self-ligation.
The SfiI regional view is a separate source declaration from Figure 1; self-ligation
keeps the repaired backbone's declared strand order and does not create those sites. Full-map sequence
identity across the preparatory reactions has not been proved from an original vector
record. The later KpnI background digest, gel selection and reclosure are physical
operations, separate from the computational tag/construct read filter. The selected
KpnI-linearized target is duplex DNA with 3′ cohesive ends, not a single-stranded circle.

🟡 The gene-side Tag-Seq site `{{= mpra.TAG_R_SITE }}` is inferred from the written
primer and final-library figure. Its exact coordinate in the full reporter map
has not been independently established. Constructors mark that site inferred by
default and preserve the flag in derived constructs and cDNA. The public
page has the caveat, affected-step captions and marked bases. PCR products incorporate
the complete written primers, so their known terminal bases replace the template-site
inference; the template-binding geometry still retains its source uncertainty. The source supplies
the {{= mpra.REPORTER_INTERVAL }}-nt interval shown in the final-library schematic;
the full sequence is not expanded here. The source describes its reporter fragment
as 1.78 kb; the model does not turn that rounded description into an exact full length.

The RNA branch starts from total RNA collected after the reported assay. It has
poly(A) selection, DNase treatment and first-strand cDNA synthesis. Potential DNA
carryover is a distinct class with unassigned amount and strandedness, and the graph
does not assert complete removal. Exact kit RT primer species remain unassigned.
The cell experiment is not modelled as a proved transformation of the illustrative
plasmid into a verified transcript. Regulatory variant DNA stays out of the cDNA
target and the final Tag-Seq amplicon.

Initial synthesis DNA and first-strand cDNA are drawn single-stranded. Their PCR
steps do not falsely claim that a complementary template strand already exists.
Intended amplicons are derived from the expected duplex-cycle binding sites and
written primer tails; the initial complementary-copy transitions are not expanded.
The SfiI re-amplification and DNA Tag-Seq binding geometry are computed separately.

Read1 reports the reverse complement of the synthesis-strand identifier. A
non-palindromic synthetic fixture checks this orientation in both DNA and RNA
branches, without proposing a study design. Construction requires the identifier
to occupy Read1 cycles 1–{{= mpra.TAG_LENGTH }} in the reported
{{= mpra.READ1_CYCLES }}-cycle acquisition. The sequencing sample index has a
separate structured feature. Index-read cycles and sample-specific index values
remain unassigned; the index-window length alone is not a run recipe. The displayed
plans cover indexed runs; the Methods also allow individual-library sequencing.

Bacterial propagation/retransformation and culture are not expanded into a population
simulation. No transformation efficiency, pool balance, coverage, RNA/DNA activity
ratio or experimental tag fidelity follows from the declared molecular sequence.

## Remaining work

- Read selected source-table design/tag assignments if a concrete published member
  is to be displayed; do not substitute a newly invented design.
- Read original pGL4.10 and pGL4.23 maps to establish full template and cleavage
  coordinates, the inferred gene-side site and exact source-vector identity.
- Resolve the kit RT primer species, ligase identity and any additional sequencing
  acquisition details from their applicable linked originals before assigning them.
- Initial complementary-copy reactions are not expanded beyond the source PCR/RT
  operation boundaries. Detailed population propagation remains outside the graph.
- Keep reporter activity, barcode association, experimental coverage and expression
  separate from the sequence relations that the chemistry validator can check.
- Cell stimulation and phenotype experiments are outside this construct/readout
  model's initial scope. No stimulation recipe is supplied by the graph.
- A second independent source reader and complete linked-source inventory remain
  pending. No missing-sequence conclusion is accepted here.
