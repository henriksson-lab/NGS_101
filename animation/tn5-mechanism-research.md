# Tn5 mechanism evidence for animation

Bounded mechanism review, 31 October 2026. This supports choreography, not a complete
source audit or a claim that any historical kit's oligos are disclosed. No protocol
module or animation code is changed by this note.

## Primary evidence inspected

- [Kia et al., 2017, *Improved genome sequencing using an engineered transposase*](https://doi.org/10.1186/s12896-016-0326-1).
  Publisher HTML: Background paragraphs beginning “Nextera, developed by Epicentre”
  and “The transposase in the Transposome complex”; Methods “Transposome assembly”
  and “Tn5/DNA binding stability assay”; Additional file 4/Figure S3 caption.
  These establish adapter transfer to target 5′ ends, a homodimer, the usual 9-base
  stagger, duplex ME with a 5′ adapter tail, and protein occupancy as an impediment
  to polymerase. The paper's assay distinguishes protein removal from gap fill and
  subsequent denaturation; its release temperatures are assay/variant-specific.
  Do not imply spontaneous Tn5 diffusion away from products or universal release
  at a particular temperature. The existing repository fetcher then retrieved originals through PMC Cloud into
  `/husky/henriksson/NGS101_papers/sources/doi_10.1186+s12896-016-0326-1/`.
  The main PDF text was checked at lines57,76,160–161,251–254,326 and352–354.
  The bundle scraper returned zero hits; that is not a source-absence conclusion.
  Original `s12896-016-0326-1_supp_12896_2016_326_MOESM4_ESM.pdf`, page1, was
  rendered and visually inspected: equal touching Tn5 lobes flank the9bp target
  overlap; the upper adapter exits left, lower adapter right. Its top target is
  5′left→3′right and bottom target3′left→5′right. It does **not** depict repair
  or establish the fate of ME′. Its OCR text is poor; the rendered original,
  not that extraction, supports this graphical observation.
- [Picelli et al., 2014, *Tn5 transposase and tagmentation procedures for massively
  scaled sequencing projects*](https://doi.org/10.1101/gr.177881.114).
  Accessible indexed article text, Results discussing SDS stripping and the opening
  72°C PCR incubation: protein release during that incubation is an interpretation
  in that specific preparation, not a measured universal Tn5 release rule. The live
  PMC article route returned a browser challenge, so this pass does not claim a new
  full-original or figure reading of this paper.

## Visual reference and its limit

[Adey, 2021, *Tagmentation-based single-cell genomics*](https://doi.org/10.1101/gr.275223.121),
Figure 1A–H and its legend, are unusually useful for animation design. This is a
**review**, not the primary chemical authority. The web tool supplied the figure's
visual extraction and indexed text; direct PDF retrieval returned 403. Its drawing
separates the transfer strand from ME′, shows two nearby insertions bracketing the
selected molecule, and identifies ME′ as “displaced or digested”. This supports a
visual comparison but does not establish which ME′ fate occurs in a particular kit.
An animation should label that fate illustrative unless local source evidence assigns it.

## Molecular geometry to preserve

The following is a coordinate derivation from the strand-transfer polarity and
9-base stagger, not an additional experimental finding. Store both strands 5′→3′;
when drawn antiparallel, the bottom runs right to left.

Let the top target be `X`, with top cuts at boundaries `a` and `b` and bottom cuts at
left-to-right coordinates `a+9` and `b+9`. Let transferred adapters `A` and `B` have
ME at their 3′ ends. Follow the central fragment:

- The initial top is `A + X[a:b]`.
- The initial bottom, in its own 5′→3′ order, is `B + revcomp(X[a+9:b+9])`.
- Both original target 3′ ends are recessed. Their first nine copied bases come
  from the retained target overhang on the other strand.
- Each 3′ end subsequently copies the ME and 5′ adapter tail on the opposite
  transferred strand. The completed top is `A + X[a:b+9] + revcomp(B)`; the bottom
  is its reverse complement.

Thus the nine target bases are **not deleted**. Adjacent repaired fragments retain
copies of the overlap. The short nontransfer ME′ oligos are separate molecules;
they must not silently become the completed complementary adapter strands. No new
ligase step follows from this geometry.

## Animation consequences

1. Show a loaded homodimer with two adapter duplexes; each duplex has a long
   transfer strand and a short nontransfer ME′ strand. The single-stranded 5′ tails
   need no arbitrary diagonal bend.
2. Depict cleavage and strand transfer as coupled chemistry. Register each new
   adapter-3′/target-5′ covalent junction, while retaining target 3′ ends.
3. Two proximal events yield the selected central library fragment. Keep outer
   DNA fragments identifiable and move them fully outside the viewport after
   protein removal; this is a presentation selection, not asserted physical waste.
4. Retain the actual recessed DNA ends into the next shot. Polymerase grows from
   those ends across the nine bases and adapter; do not replace the scene with a
   complete duplex. ME′ departure during copying may illustrate displacement,
   but do not present prior thermal dissociation as established.
5. Show indexing PCR as a later operation. Complete flow-cell/index/entry-point/ME
   anatomy must come from the chosen protocol and version, not the mechanism paper.

The Nextera-XT local bundle scraper was run as a retrieval check. No modern XT
manual was used here to fill original SMART-seq/Epicentre kit identities. No missing
sequence or reagent was declared absent.

## Follow-up: RNA displacement is a separate Phi29 capability

The BOLORAMIS audit and primary evidence are recorded in
[boloramis-rna-followup.md](boloramis-rna-followup.md). Inspected BOLORAMIS main
Methods lines160–183/214–221 and SI lines85–96 transition from ligation/rinses into
RCA without an inserted RNase step. This is a scoped recipe observation, not an
absence audit.

The original [Kim et al., doi:10.1093/nar/gkad1101](https://doi.org/10.1093/nar/gkad1101)
PDF was retrieved externally. Its Results lines407–440/Figure3G–H positively support
Phi29 RNA-hybrid displacement on a passing template-strand substrate, while also
showing stalling and retained RNA. Hence a selected successful engine branch may
explicitly allow RNA displacement, preserving RNA backbone and nucleotide IDs.
Do not derive that permission solely from DNA displacement, infer an RNase, or
claim this experiment directly validates native BOLORAMIS occupancy or yield.
