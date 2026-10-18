# snRandom-seq — source notes

Defining paper: Xu et al., *Nature Communications* 2023, [doi:10.1038/s41467-023-38409-5](https://doi.org/10.1038/s41467-023-38409-5).

## Sources read

- 🟢 PMC10182092 full-text XML, Methods from “In situ DNA blocking” through “Sequencing”, plus preprocessing.
- 🟢 Supplementary Note 1 protocol, Supplementary Table 1 primer sequences and Supplementary Table 2 three-round barcode lists in `41467_2023_38409_MOESM1_ESM.pdf`.

## Evidence trail

- 🟢 FFPE nuclei undergo deparaffinization, crosslink reversal and permeabilization. A handle+N7 block primer is extended first to reduce genomic-DNA capture.
- 🟢 Five handle-bearing random primers and five handle-bearing oligo(dT) primers are used for in-situ RT with twelve annealing ramps from 8 to 42 °C.
- 🟢 TdT plus dATP adds a dA tail to cDNA.
- 🟢 Acrydite hydrogel beads carry a deoxyuridine-containing oligo and are assembled through three split-pool ligations. The published product is `[10-base barcode1]TGGT[10-base barcode2]GAGA[10-base barcode3]N8T20`.
- 🟢 Droplet extension transfers the bead barcode to cDNA; PCR follows emulsion breakage.
- 🟢 A named Vazyme universal Illumina kit receives fragmented DNA and performs end repair/adenylation, adapter ligation and library PCR; the inspected paragraph does not assign the fragmentation mechanism. Paired-end 150-base NovaSeq data use Read 1 for the 30-nt cell barcode plus 8-nt UMI and Read 2 for expression.
- 🟡 Native Vazyme adapter strands and index allocation are unassigned in the current model; the version-matched vendor original remains deferred. The earlier inferred canonical outer arms have been removed.

## Oligos as printed

- Block primer: `GAGAATGTGAGTGAAGATGTATGGTGANNNNNNN`
- Random primers 1–5, respectively: `GGAGTTGGAGTGAGTGGATGAGTGATGGAAGGAATNNNNNNN`, `GGAGTTGGAGTGAGTGGATGAGTGATGGATGAATGNNNNNNN`, `GGAGTTGGAGTGAGTGGATGAGTGATGGATATGTGNNNNNNN`, `GGAGTTGGAGTGAGTGGATGAGTGATGGTGGTGTANNNNNNN`, `GGAGTTGGAGTGAGTGGATGAGTGATGGTAGGTTANNNNNNN`.
- Oligo(dT) primers 1–5, respectively: `GGAGTTGGAGTGAGTGGATGAGTGATGGAAGGAATTTTTTTTTTTTTTTT`, `GGAGTTGGAGTGAGTGGATGAGTGATGGATGAATGTTTTTTTTTTTTTTT`, `GGAGTTGGAGTGAGTGGATGAGTGATGGATATGTGTTTTTTTTTTTTTTT`, `GGAGTTGGAGTGAGTGGATGAGTGATGGTGGTGTATTTTTTTTTTTTTTT`, `GGAGTTGGAGTGAGTGGATGAGTGATGGTAGGTTATTTTTTTTTTTTTTT`.
- Bead anchor as printed: `/5Acryd/ATTATATATAT U GTG AGT GAT GGT TGA GGA TGT GTG GAGATA`; the separated `U` is deoxyuridine. Its precise release chemistry and retained terminal boundary remain unassigned in the inspected source scope. The three ten-base barcode pools are listed in Supplementary Table 2; their individual members are not interchangeable with one fixed barcode.
- Hydrogel anchor: `/5Acryd/ATTATATATAT U GTG AGT GAT GGT TGA GGA TGT GTG GAGATA`
- Barcoded bead primer: `[10 bases barcode1] TGGT [10 bases barcode2] GAGA [10 bases barcode3] NNNNNNNN TTTTTTTTTTTTTTTTTTTT`
- PCR primers: `GGAGTTGGAGTGAGTGGATGAGTGATG` and `GTGAGTGATGGTTGAGGATGTGTGGAGATA`.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC10182092 JATS230–256, DNA block, in-situ RT, dA tailing, barcoding and library preparation.

Specified twelve 8–42°C RT annealing ramps followed by 42°C/30min, and TdT tailing at 37°C/30min. Confirmed droplet extension precedes PFO recovery, preamplification and Vazyme library conversion with paired150 reads.

Deferred: Published SupplementaryTable1 modified primer families, exact droplet-release chemistry, physical wash/PFO/AMPure stages, PCR cycle counts and vendor end geometry. This pass does not establish a complete source review or publication-wide absence.

The current endpoint removes the earlier inferred canonical outer-adapter shell. Native vendor ends and index acquisition are unassigned; inline barcode/UMI architecture remains source-qualified. A representative variant1 RT handle is used to explain strand order; other RT primer variants are separate molecules. No fixed TdT-tail length or deoxyuridine cleavage site is assigned.

## Native graph migration, batch 5 (2026-10-18)

Inspected: Published PMC10182092 JATS Methods230–266: DNA block, RT, dA tailing, bead synthesis, droplet barcoding, library preparation and preprocessing; original MOESM1 SupplementaryInformation Table1 (text190–260).

Removed generic TruSeq dual-index shell, fixed12-A TdT tract and unsupported enzymatic fragmentation inference. Preserved source-specific barcode/UMI features through role-level fragments, native adapter ends and library. Separate genome/RNA/cDNA populations, three washes after each nuclear operation, source deoxyuridine/acrydite bead tether, PFO aqueous recovery, AMPure before/after cDNA PCR, native endprep/ligation/size-selection/PCR/final cleanup. Both PCR geometries use actual precursors. Nonamplified genomic DNA/RNA remains unresolved; cleanup fractions do not assert species-specific removal. Native release site and fragmentation mechanism remain unassigned.

Deferred: Vazyme ND607-01 version-matched adapter and PCR source; exact fragmentation mechanism and native read-primer geometry; Bead predecessor Delley2021, full five-member RT primer population/whitelist scope and source-data workbook; no new absence claim accepted. This capped pass is not a complete source review and does not accept publication-wide absence.

🟢 The bead-primer Table1 line is `NNNNNNNNTTTTTTTTTTTTTTTTTTTT`: N8 followed by T20. This corrects this repository's earlier T21 transcription; it is not a publication discrepancy.
