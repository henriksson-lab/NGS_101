# snRandom-seq — source notes

Defining paper: Xu et al., *Nature Communications* 2023, [doi:10.1038/s41467-023-38409-5](https://doi.org/10.1038/s41467-023-38409-5).

## Sources read

- 🟢 PMC10182092 full-text XML, Methods from “In situ DNA blocking” through “Sequencing”, plus preprocessing.
- 🟢 Supplementary Note 1 protocol, Supplementary Table 1 primer sequences and Supplementary Table 2 three-round barcode lists in `41467_2023_38409_MOESM1_ESM.pdf`.

## Evidence trail

- 🟢 FFPE nuclei undergo deparaffinization, crosslink reversal and permeabilization. A handle+N7 block primer is extended first to reduce genomic-DNA capture.
- 🟢 Five handle-bearing random primers and five handle-bearing oligo(dT) primers are used for in-situ RT with twelve annealing ramps from 8 to 42 °C.
- 🟢 TdT plus dATP adds a dA tail to cDNA.
- 🟢 Acrydite hydrogel beads carry a deoxyuridine-containing oligo and are assembled through three split-pool ligations. The published product is `[10-base barcode1]TGGT[10-base barcode2]GAGA[10-base barcode3]N8T21`.
- 🟢 Droplet extension transfers the bead barcode to cDNA; PCR follows emulsion breakage.
- 🟢 A named Vazyme universal Illumina kit performs fragmentation, end repair/adenylation, adapter ligation and library PCR. Paired-end 150-base NovaSeq data use Read 1 for the 30-nt cell barcode plus 8-nt UMI and Read 2 for expression.
- 🟡 The paper does not print the Vazyme adapter sequence. The schematic therefore marks the outer Illumina-compatible arms inferred and does not present them as vendor-transcribed bases.

## Oligos as printed

- Block primer: `GAGAATGTGAGTGAAGATGTATGGTGANNNNNNN`
- Random primers 1–5, respectively: `GGAGTTGGAGTGAGTGGATGAGTGATGGAAGGAATNNNNNNN`, `GGAGTTGGAGTGAGTGGATGAGTGATGGATGAATGNNNNNNN`, `GGAGTTGGAGTGAGTGGATGAGTGATGGATATGTGNNNNNNN`, `GGAGTTGGAGTGAGTGGATGAGTGATGGTGGTGTANNNNNNN`, `GGAGTTGGAGTGAGTGGATGAGTGATGGTAGGTTANNNNNNN`.
- Oligo(dT) primers 1–5, respectively: `GGAGTTGGAGTGAGTGGATGAGTGATGGAAGGAATTTTTTTTTTTTTTTT`, `GGAGTTGGAGTGAGTGGATGAGTGATGGATGAATGTTTTTTTTTTTTTTT`, `GGAGTTGGAGTGAGTGGATGAGTGATGGATATGTGTTTTTTTTTTTTTTT`, `GGAGTTGGAGTGAGTGGATGAGTGATGGTGGTGTATTTTTTTTTTTTTTT`, `GGAGTTGGAGTGAGTGGATGAGTGATGGTAGGTTATTTTTTTTTTTTTTT`.
- Bead anchor as printed: `/5Acryd/ATTATATATAT U GTG AGT GAT GGT TGA GGA TGT GTG GAGATA`; the separated `U` is the deoxyuridine used to release the synthesized bead primer. The three ten-base barcode pools are listed in Supplementary Table 2; their individual members are not interchangeable with one fixed barcode.
- Hydrogel anchor: `/5Acryd/ATTATATATAT U GTG AGT GAT GGT TGA GGA TGT GTG GAGATA`
- Barcoded bead primer: `[10 bases barcode1] TGGT [10 bases barcode2] GAGA [10 bases barcode3] NNNNNNNN TTTTTTTTTTTTTTTTTTTTT`
- PCR primers: `GGAGTTGGAGTGAGTGGATGAGTGATG` and `GTGAGTGATGGTTGAGGATGTGTGGAGATA`.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC10182092 JATS230–256, DNA block, in-situ RT, dA tailing, barcoding and library preparation.

Specified twelve 8–42°C RT annealing ramps followed by 42°C/30min, and TdT tailing at 37°C/30min. Confirmed droplet extension precedes PFO recovery, preamplification and Vazyme library conversion with paired150 reads.

Deferred: Published SupplementaryTable1 modified primer families, exact droplet-release chemistry, physical wash/PFO/AMPure stages, PCR cycle counts and vendor end geometry. This pass does not establish a complete source review or publication-wide absence.
