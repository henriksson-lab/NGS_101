# PARS

**Defining source:** Kertesz *et al.*, *Nature* (2010), doi:10.1038/nature09322.
**Detailed protocol:** Wan *et al.*, *Nature Protocols* (2013),
doi:10.1038/nprot.2013.045. The PMC author manuscript, supplementary workflow and
protocol description were read.

## Evidence used

- 🟢 Poly(A)-selected RNA is refolded in vitro and divided between RNase V1 and S1.
- 🟢 V1 cleaves double-stranded and S1 cleaves single-stranded regions, both producing
  informative 5′-phosphate ends.
- 🟢 Random fragmentation precedes selective 5′-adapter ligation: random 5′-OH ends do
  not ligate, whereas nuclease-created 5′-P ends do. Phosphatase treatment and second
  adapter ligation are followed by RT-PCR and SOLiD sequencing.
- 🟡 Exact proprietary SOLiD adapter and sequencing-primer oligos are not printed.

## Construct logic

PARS is intrinsically a paired experiment. The molecular endpoint is a cleavage-defined
read start in each branch; the normalized log₂ V1/S1 ratio assigns double- versus
single-stranded character to the nucleotide immediately upstream.
