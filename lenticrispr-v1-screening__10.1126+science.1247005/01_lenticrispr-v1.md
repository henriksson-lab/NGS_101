# lentiCRISPR v1 screening

Shalem et al. introduced the all-in-one lentiCRISPR v1 vector and used it for pooled loss-of-function screening. The page keeps guide cloning, vector organization, and library readout visually distinct.

## Sources read

- 🟢 Shalem et al. 2014, *Science*, [doi:10.1126/science.1247005](https://doi.org/10.1126/science.1247005) — defining v1 screen. Read as the PubMed Central author manuscript `science.1247005_PMC4089965.xml.txt` (abstract line 42, Results lines 53–69, Supplementary Material pointer lines 71–73, figure captions lines 255–267).
- 🟢 Addgene plasmid #49535 — deposited lentiCRISPR v1 map; retrieval details are in `ref/plasmids/MANIFEST.md`.
- 🟢 Joung et al. 2017, [doi:10.1038/nprot.2017.016](https://doi.org/10.1038/nprot.2017.016) — compatible one-PCR readout shown at the end of the page.
- 🟢 Zhang-lab GeCKO readout-primer workbook, July 2014 (`GeCKO-plasmid-readout-primers-July2014_web20141106225441.xlsx.txt`, sheet `GeCKO plasmid readout`, header notes lines 2–4, forward primers F01–F12 lines 6–18, reverse primers R01–R12 lines 19–32, v2 adaptor pair lines 33–36). The author-hosted original was recovered on 2026-10-17 from the Internet Archive snapshot (2014-11-06) of `genome-engineering.org/gecko/wp-content/uploads/2013/12/`; the live URL still does not respond.

## Boundary

The original Shalem nested readout, the later one-PCR readout, and lentiCRISPR v2 are separate protocols. The older combined investigation remains in `90_combined_screen_research.md` for provenance, not as the public schematic.

## Oligo review boundary

The deposited Addgene map fixes the v1 vector cassette, and Joung's later paper prints
the one-PCR library primers shown on this page. 🔴 The accessible Shalem author
manuscript does not include its gated supplement, so the original experiment's nested
readout cannot be verified from that text. The Zhang-lab workbook
`GeCKO-plasmid-readout-primers-July2014.xlsx` has now been read (next section): it is an
author original for a v1-compatible readout primer set, but it is a *plasmid*-readout
sheet dated July 2014 and is not the 2014 screen's genomic-DNA Methods, so it does not
by itself establish what Shalem et al. did with screen gDNA.

## Screen parameters taken from the article (bounded pass, 2026-10-17)

🟢 From the accessible author manuscript: the GeCKO library targets **18,080 genes with
64,751 unique guide sequences**, 3–4 guides per gene (lines 42 and 59); transduction is at
**MOI 0.3** into A375 melanoma and HUES62 cells, with guide distributions compared **3 and
14 days** post-transduction and a 14-day vemurafenib (PLX) arm (lines 57, 61, 63, 255–261).
The EGFP pilot adds **puromycin selection** after transduction (line 57); the article text
does not restate a selection step for the pooled screen itself, so the page and graph do
not model one.

🔴 The accessible text contains no Methods section. Everything about library amplification
and the nested two-PCR readout is in the Supplementary Information: the PMC supplementary
PDF (`NIHMS609112-supplement-Supplementary_Information.pdf`) sits behind PMC's download
gate and the PMC Cloud bucket `PMC4089965.1/` holds only JSON/TXT/XML, no supplement. A
fresh attempt (2026-10-17) at the live
`genome-engineering.org/gecko/wp-content/uploads/2013/12/GeCKO-plasmid-readout-primers-July2014.xlsx`
again timed out, and the Europe PMC supplementary endpoint rejects PMC4089965 as not open
access. These are retrieval failures and do **not** establish that the original primers
were never published; the Shalem Methods' own primer pair, its barcode scheme and the
acquired read lengths stay unassigned.

## GeCKO readout-primer workbook (recovered 2026-10-17)

The author-hosted workbook was obtained from the Internet Archive copy of the Zhang-lab
page and read directly. 🟢 Its header note states that **for lentiCRISPR v1 and
lentiGuide-Puro a single PCR with the forward and reverse Illumina primers is performed**,
with a mix of forward primers so that the stagger avoids monotemplate issues, while
**lentiCRISPRv2 needs a first adaptor PCR** to add the reverse priming site (lines 2 and
33–36). It is explicit that this note covers *plasmid* sequencing and that readout from
genomic DNA is described in the paper's Supplementary Methods (line 4) — so it does not
settle the screen's own gDNA amplification.

🟢 Architecture of the twelve primer pairs (lines 6–32):

- forward `F01`–`F12`: P5 + TruSeq Read 1 adapter, a 1–9 nt stagger, an **8-nt inline
  barcode**, then the 24-nt priming site `tcttgtggaaaggacgaaacaccg` (the hU6 3′ end plus
  the Pol III +1 G, the same site the drawn Joung forward primer uses), 91–99 nt total,
  33–41 nt of which are read.
- reverse `R01`–`R12`: P7 + an **8-nt inline barcode** + TruSeq Read 2 adapter + a 1–9 nt
  stagger + the 25-nt priming site `TCTACTATTCTTTCCCCTGCACTGT`, which is **not** inside
  the 76-nt guide scaffold (computed against `crispr.SCAFFOLD_V1`) but downstream of it in
  the vector; the historical comparison in `90_combined_screen_research.md` places it in
  cPPT, which is downstream of the guide cassette in v1 only.
- 🟡 Both sample barcodes are therefore **inline, read inside Read 1 / Read 2**, unlike the
  i7-indexed Joung library drawn on the page. The workbook prints no read lengths, no
  thermocycling and no Illumina index, so the acquired read layout of the 2014 screen is
  still unassigned.
