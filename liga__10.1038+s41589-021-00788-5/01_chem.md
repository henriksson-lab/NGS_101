# Genetically encoded multivalent liquid glycan array displayed on M13 bacteriophage

Defining DOI: 10.1038/s41589-021-00788-5. Status: capped source pass; incomplete private model, no finished public schematic.

Sources actually inspected:

- https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41589-021-00788-5/MediaObjects/41589_2021_788_MOESM53_ESM.pdf: `s41589-021-00788-5_41589_2021_788_MOESM53_ESM.pdf.txt` — Published SI cover identity and SupplementaryFig1–3 captions1573–1602; glycosylation/viability experiments only.
- https://pmc-oa-opendata.s3.amazonaws.com/PMC8380037.1/PMC8380037.1.xml: `s41589-021-00788-5_PMC8380037.xml.txt` — Associated accepted manuscript (manifest author manuscript true), Methods1.8–1.15 lines179–243; operational production and selected CD22 ProteinG capture route, not established whole-VOR equivalence.

🟢 Published SI Fig1 explicitly distinguishes a viability experiment: 1mM DBCO-NHS1h →Zeba purification →2mM azidomethyl-thio-mannose8h4°C. Fig2 describes competing DBCO-NHS hydrolysis and nonuniform acylation. Fig3 describes a minority doubly modified pVIII species. Coupling yield, residue occupancy and affinity are not computed from these captions.

🟢 Associated main cache is an **accepted manuscript**, not the final published article. It is retained as a specifically labelled operational interface because the published main method route was not obtained. No complete version equivalence is asserted. Its production1.8 recipe is a separate selected route: DBCO-NHS0.1–1.5mM45min →add glycan-azide directly2–3mM overnight4°C →Zeba purification/MALDI →mix desired glycan-phage clones. It does **not** inherit the intermediate Zeba step from the SI viability experiment.

🟢 Selected MS1.11.3 ProteinG/CD22-Fc magnetic capture route: bind prepared glycan-phage mixture →two HEPES/1%Tween washes →one HEPES wash →transfer bead-associated preparation into30µLwater →save2µL titer aliquot →90°C15min remainder →collect25µL supernatant for PCR. Residual beads/liquid/debris remain explicitly outside that selected PCR fraction. Other source plate acid-elution and cell-pellet routes are not merged.

🟢 MS1.15: 50µL indexed Phusion PCR (full thermal/oligo geometry deferred) →2%agarose quantification →pool10ng of differently indexed products →eGel purification →Qubit →NextSeq500/550 HighOutputv2.5 paired-end2×75. PCR creates new DNA copies; neither phage protein nor original modified glycan automatically becomes copied DNA.

The silently encoded phage genotype maps externally to glycan identity/density. It is neither an acquired current-cell barcode nor a per-molecule UMI. Full source genome and barcode segments are unassembled rather than represented with arbitrary bases. No deterministic glycan/protein affinity, quantitative modification or perfect wash purity is claimed.

Component lineage: source discusses silent chemical encoding and other DNA-coded glycan arrays, but selected sections do not establish importable full recipe identity. No LiGA→PhIP or generic barcode-platform inheritance edge is inferred.

Deferred targets (not claims of absence):

- Published main final Methods or version equivalence of selected accepted-manuscript operations
- Published detailed Supplementary Methods with production, binding and PCR thermal/primer details; printed final SI inspected only in selected figure-caption scope
- Full source silent-genome sequences, barcoded PCR primer bodies/site geometry and barcode→glycan dictionary mapping; additional plate/cell/lectin branches
- Genetic construction/manufacturing route and cited direct predecessors remain unassembled; separate viability condition variant not merged

Verification: mandatory standard oligo scraper was run over each fetched source bundle; source pointers and relevant modifications were checked against the originals stated above. Local tests check transcription and consequential fraction/order distinctions; graph validation is separate from molecular/source completeness. No absence audit or absence acceptance.

Physical QC branch: analytical aliquots remain separate from the production preparation. MALDI matrix-exposed, qPCR-amplified, gel-analysed or Qubit-assayed material is not automatically pooled for sequencing. Selected source quantities are retained; unreviewed analytical allocation amounts stay unspecified.
