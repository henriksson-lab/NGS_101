# PAIso-seq — research notes

## Sources read

- 🟢 Liu et al. 2019, DOI `10.1038/s41467-019-13228-9`, PMC full text and Supplementary Table 1.
- 🟢 Liu et al. 2022, DOI `10.1038/s41596-022-00704-8`,
  [author-uploaded detailed protocol](https://www.researchgate.net/publication/361976122_Transcriptome-wide_measurement_of_polyA_tail_length_and_composition_at_subnanogram_total_RNA_sensitivity_by_PAIso-seq),
  Table 2, Materials and steps 29–30.
- 🟢 [PacBio SMRTbell Template Prep Kit 1.0 insert](https://www.pacb.com/wp-content/uploads/2015/09/Insert-SMRTbell-Template-Prep-Kit-1.0.pdf), p. 1, component list: `Blunt Adapter` PN 100-348-100 and End Repair Mix.
- 🟢 [PacBio SMRTbell Express Template Prep Kit 2.0 overview](https://www.pacb.com/wp-content/uploads/SMRTbell-Express-TPK-2.0-Overview-Large-Insert-gDNA-Library-Prep-Customer-Training.pdf), pp. 6 and 8: Overhang Adapter v3 has a T overhang and the kit uses A/T ligation instead of blunt ligation. These are kit-version comparators; neither is identified as the kit used in the PAIso-seq experiment.

## Construct path

🟢 A barcoded oligo containing strategically placed deoxyuridines supports end extension over the RNA poly(A) tail. USER cleavage removes that oligo before reverse transcription. Template switching adds the opposite PCR handle, and full-length cDNA is converted to a PacBio SMRTbell. This preserves tail length and internal non-A residues in the same CCS read as the transcript isoform.

🟢 The TSO, RT primer, PCR primer, and all fifteen 2dU-BC extension oligos on the schematic are transcribed from Supplementary Table 1. BC3 is the member explicitly named for the bulk library in the defining paper's Methods; BC1–BC15 barcode individual oocytes.

🟢 **The TSO linker is `ACAT`, published by PAIso-seq itself.** Supplementary Table 1 row
`TSO` reads `AAGCAGTGGTATCAACGCAGAGTACATrGrG+G` verbatim, i.e. the 23-nt SMART handle, the
four-base linker `ACAT` and the `rGrG+G` tail. This matters because SMART-seq2 has two
published linkers — `AC` in Picelli et al. 2013 Supplementary Table 2 for the chosen
`rGrG+G` oligo, `ACAT` in the *Nature Protocols* form — and a protocol that only cites
SMART-seq2 cannot tell them apart. PAIso-seq prints its own, so no base is inherited by
assumption. `tools/chem.py` therefore takes the bases from `smart-seq2__10.1038+nmeth.2639/tools/smartseq2` (handle,
`DT_LINKER` for the RT primer's `AC`, `TSO_LINKER`, `TSO`) rather than retyping them, and
`tools/selftest.py` pins all three ordering strings against Supplementary Table 1 so a
change in the shared constant would fail here rather than silently redraw this page.

🟡 The authors' 2022 detailed protocol is a later version: its Table 2 TSO adds
5′ `iCiGiC` and a 3′ `/NH2-C3/` block around the same core sequence. The
2019 defining Supplementary Table 1 does not print those modifications. This
page models the 2019 oligo; the later ordered oligo must not be substituted
for the original experiment.

🟢 Amplification uses the single published PCR primer. Template switching and reverse transcription place its common handle at opposite ends of the cDNA, so the same physical primer species binds both inward-facing sites. The sample barcode is carried by the extension oligo; it is not a second PCR primer.

## SMRTbell junction evidence

🟢 The 2019 PAIso-seq Online Methods, “PAIso-seq library construction,” identify size-selected amplified cDNA and an unversioned `SMRTbell Template Prep Kit`, followed by sequencing-primer annealing and Sequel sequencing. They do not specify dA tailing or the adapter's ligating end. The 2022 detailed protocol explicitly names `SMRTbell template prep kit 1.0-SPv3` (PacBio 100-991-900) in Materials and step 29; it also uses a modified TSO, so its kit choice cannot automatically be assigned to the 2019 experiment. The official PacBio 1.0 insert lists a Blunt Adapter, whereas the contemporaneous Express 2.0 overview describes T-overhang adapters and A/T ligation. 🟡 The model therefore retains only the source-supported closed SMRTbell topology and records the 2019 end-preparation and ligation geometry as unresolved; it does not draw an A/T or blunt junction for that experiment. Exact hairpin and sequencing-primer bases remain unverified.
