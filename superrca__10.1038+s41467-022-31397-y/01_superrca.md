# superRCA — source notes

**Defining source:** Chen L, Eriksson A, Weström S, et al. *Ultra-sensitive monitoring of leukemia patients using superRCA mutation detection assays.* Nature Communications 13:4033 (2022). [doi:10.1038/s41467-022-31397-y](https://doi.org/10.1038/s41467-022-31397-y). [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9276831/).

## Scope and evidence

- 🟢 The authors call the method **superRCA**. It is a DNA mutation-detection assay, with fluorescent flow-cytometry or microscopy readout. The targeted NGS and ddPCR procedures in the paper are comparisons, not superRCA outputs (Introduction; Fig. 1; Methods).
- 🟢 The paper provides the reaction order and reagent roles below (Fig. 1 and Methods, “High fidelity PCR pre-amplification” through “Digital recording”).
- 🟢 The original Supplementary Information was retrieved via PMC Cloud and read alongside the article. It contains Supplementary Table 1 (method comparison), figures and supporting notes, but no oligo table or exact assay oligo sequences. The accompanying source-data workbook contains figure data. In the published peer-review file, a reviewer requested all PCR, ligation-template, padlock, RCA-primer and fluorescent-probe sequences, and the authors replied that they added Supplementary Table S2; that table is absent from the available supplementary file. This is a specific unresolved discrepancy between the reply and the published attachment.
- 🔴 No exact PCR primer, ligation-template, padlock, secondary-RCA primer or fluorescent-detection probe sequence is established by the checked defining article and supplement. The public diagram uses named role tokens without claiming nucleotide lengths or bases. A later source may disclose them, but its reagents must not silently replace the 2022 assay.

## Molecular path

1. 🟢 Genomic DNA from bone marrow or whole blood is enriched by 15 cycles of high-fidelity targeted PCR using forward and reverse primers (Methods, “Extraction of genomic DNA” and “High fidelity PCR pre-amplification”). The paper does not print the primer sequences in its main text.
2. 🟢 A ligation-template oligonucleotide complementary to both ends of one denatured amplicon strand brings its 5′ and 3′ ends together. Ampligase seals that strand into a single-stranded DNA circle (Fig. 1B; Methods, “Ligase-mediated circularization”). 🟡 The circularized strand must have ligatable end chemistry; the main text does not specify how the PCR strand acquires its 5′ phosphate, so the diagram does not assert a particular terminal modification.
3. 🟢 The ligation template then serves as primer for Phi29 RCA of the amplicon circle, making a concatemer with repeated copies of the target region (Fig. 1C; Methods, “Target sequence amplification”).
4. 🟢 Allele-specific mutant and wild-type padlock probes anneal to copies in that primary RCA product. Ampligase circularizes successfully matched probes around it (Fig. 1D–E; Methods, “Genotyping of RCA products”). This is direct ligation of padlock ends, not a reported gap-fill MIP reaction.
5. 🟢 Added primers initiate a second Phi29 RCA from the padlock circles. The resulting products stay associated with the primary RCA product, forming one large superRCA cluster per starting circle (Fig. 1E–F; Methods, “Secondary RCA”).
6. 🟢 Fluorophore-labelled oligonucleotides bind allele-specific secondary repeats. A standard flow cytometer counts mutant and wild-type fluorescent objects (Fig. 1F; Methods, “Digital recording”).

## Construct and readout boundary

The first circle contains one strand of the target amplicon. Primary RCA copies that circle repeatedly. Each allele-specific padlock circle supplies a second, allele-labelled repeating template. A superRCA object is the associated primary concatemer plus many secondary RCA products; it is **not** an Illumina library. No sequencing primer, index read, UMI or nucleotide read layout applies to this assay. The main-text schematic is a role-level molecular path until the exact oligos can be checked in the supplement.

## Open source questions

- Locate the Supplementary Table S2 mentioned in the peer-review response, or another author-supplied source explicitly identifying the 2022 oligos, before adding base-by-base circles or probe junctions.
- Establish from the source how the amplicon strand obtains the 5′ phosphate needed for covalent closure; do not infer it from generic padlock chemistry.
