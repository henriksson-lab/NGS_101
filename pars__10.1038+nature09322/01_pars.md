# PARS

**Defining source:** Kertesz *et al.*, *Nature* (2010), doi:10.1038/nature09322.
**Detailed protocol:** Wan *et al.*, *Nature Protocols* (2013),
doi:10.1038/nprot.2013.045. The PMC author manuscript, supplementary workflow and
protocol description were read. The 2010 Supplementary Information was recovered from
Springer Nature's static supplement and the [SOLiD Small RNA Expression Kit guide](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/cms_054973.pdf)
(Applied Biosystems PN 4399434b, 2009) was checked.

## Evidence used

- 🟢 Poly(A)-selected RNA is refolded in vitro and divided between RNase V1 and S1.
- 🟢 V1 cleaves double-stranded and S1 cleaves single-stranded regions, both producing
  informative 5′-phosphate ends.
- 🟢 Random fragmentation precedes selective 5′-adapter ligation: random 5′-OH ends do
  not ligate, whereas nuclease-created 5′-P ends do. Phosphatase treatment and second
  adapter ligation are followed by RT-PCR and SOLiD sequencing.
- 🟢 Supplementary Methods lines 90–127 specify Adaptor Mix A of the SOLiD Small RNA
  Expression Kit, phosphatase treatment, the second adapter ligation and kit PCR
  primers. The Thermo Fisher guide, Appendix A p. 26, prints the common P1/5′ PCR
  primer and all ten indexed 3′ PCR primer sequences. The public oligo list shows
  index set 1 as a **kit example** because the paper does not identify which of the
  ten sets was assigned to each library.
- 🔴 The guide describes Mix A as a set of oligos with a degenerate end and a SOLiD
  arm but does not print its full ordered sequences or terminal modifications.
  Its protocol and the defining paper likewise do not print the instrument
  ligation-sequencing primer. These remain the exact-structure boundary.

## Construct logic

PARS is intrinsically a paired experiment. The molecular endpoint is a cleavage-defined
read start in each branch; the normalized log₂ V1/S1 ratio assigns double- versus
single-stranded character to the nucleotide immediately upstream.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC3847670 JATS68–82, Methods Summary and SOLiD mapping.

Made the original SOLiD branch explicit and recorded truncation to35 bases for mapping without inventing acquired cycle count. Added caveat that original supplementary library methods still need checking; later 2013 protocol cannot substitute for that verification.

Deferred: Original supplementary methods, random fragmentation/order, selective adapter capture geometry, SOLiD vendor original and actual acquisition; no source absence accepted. This pass does not establish a complete source review or publication-wide absence.
