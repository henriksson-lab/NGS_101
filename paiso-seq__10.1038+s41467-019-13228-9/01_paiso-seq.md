# PAIso-seq — research notes

## Current native graph scope — 2026-10-18 batch11

Two originals were inspected: published PMC6876564 Methods144–160 (bulk/single
oocyte branches, extension/USER/RT/PCR/size selection/Sequel/read processing) and
original MOESM1 PDF SupplementaryTableS1, all native primer/15guide rows and
modification definitions. Mandatory scraper:64hits in6files. No new reading of
the later2022 protocol, vendor comparators or author repository was performed.
Those existing records remain historical evidence rather than current kit identity.

🟢 The selected bulkBC3 guide is protocol-owned. Klenow with dNTP extends the
originalRNA3′end with DNA, copying the barcode and RT-guide region beyond the
nativeRNA polyA tail. The product is covalent RNA–DNA, not newly transcribedRNA.
USER cleaves guideDNA at its two physicaldU sites; column recovery>200nt does not
by itself establish complete removal of every guidefragment. RT with the printed
primer and original2019 rGrG+G TSO produces actual switchedssDNA, which feeds
single-species KAPA PCR. The originalRNA/guide/RT/TSO lineages and their chemical
modifications are distinct from PCR daughters.

🟡 Two actual portions independently undergo1×PurePB>200bp and0.4×PurePB>2kb
selection, then equimolar contributions are physically pooled. Those overlapping
size classes are not falsely represented as disjoint bins; remaining unallocated
fractions have unresolved fates. The2019 SMRTbell kit version/end-preparation/
hairpin chemistry remains a scoped vendor boundary. A closed circuit retains the
source cDNA/complement and two unknown native hairpin roles; there is no invented
18ntN hairpin or A/T junction. Known sequencing-primer annealing, polymerase binding
and Magbead loading follow the native boundary. CCS≥10pass selection and ambiguous
barcode removal are computational operations, not physical purification.

Existing public Results-vs-TableS1 notices for dU/T7-vs-T8 spacing and nonidentical
TSO/RT/PCR handle descriptions remain. The exact sourceTableS1 ordering strings
win for depicted native reagents, without borrowing Smart-seq2 defaults.

Deferred: cited Smart-seq2 detailed RT conditions/native junctions, unversioned2019
SMRTbell vendor identity and internal operation/ends, hairpin/loaded-primer bases,
class-specific residual recoveries. Single-oocyte15-member pooling and later2022
workflow remain distinct alternatives. This pass accepts no absent information or
historical status upgrade.


## Sources read

- 🟢 Liu et al. 2019, DOI `10.1038/s41467-019-13228-9`, PMC full text and Supplementary Table 1.
- 🟢 Liu et al. 2022, DOI `10.1038/s41596-022-00704-8`,
  [author-uploaded detailed protocol](https://www.researchgate.net/publication/361976122_Transcriptome-wide_measurement_of_polyA_tail_length_and_composition_at_subnanogram_total_RNA_sensitivity_by_PAIso-seq),
  Table 2, Materials and steps 29–30.
- 🟢 [PacBio SMRTbell Template Prep Kit 1.0 insert](https://www.pacb.com/wp-content/uploads/2015/09/Insert-SMRTbell-Template-Prep-Kit-1.0.pdf), p. 1, component list: `Blunt Adapter` PN 100-348-100 and End Repair Mix.
- 🟢 [PacBio SMRTbell Express Template Prep Kit 2.0 overview](https://www.pacb.com/wp-content/uploads/SMRTbell-Express-TPK-2.0-Overview-Large-Insert-gDNA-Library-Prep-Customer-Training.pdf), pp. 6 and 8: Overhang Adapter v3 has a T overhang and the kit uses A/T ligation instead of blunt ligation. These are kit-version comparators; neither is identified as the kit used in the PAIso-seq experiment.

## 2026-10-17 bounded pass

Re-read the published PMC6876564 JATS Results line 54 (Fig. 1a design paragraph) and
Online Methods lines 142–148 (bulk and single-oocyte library construction), Supplementary
Information PDF Supplementary Table 1 (text lines 688–745), and the article PDF text line
158 for the guide-oligo 3′ end. Reaction details now in the graph, all 🟢 from Methods
line 144: end extension with dNTP and Klenow fragment (exo−, NEB) on the BC3-2dU template
at 37 °C for 1 h; USER (NEB) at 37 °C for 30 min; RNA Clean & Concentrator-5 cleanup of
end-extended RNA (>200 nt), which is the represented separation of the short cleaved
guide fragments; KAPA HiFi HotStart ReadyMix PCR; Pure PB size selection at 1× (>200 bp)
and 0.4× (>2 kb) with the two fractions combined equimolar before SMRTbell preparation.

Published discrepancies (shown on the page):

- 🟢 Results (JATS line 54; PDF text line 158) prints the guide 3′ end as
  `5′-dUTTTTTTTdUTTTTTTT-3′` (seven T after each dU); every 2dU-BC row of
  Supplementary Table 1 prints `dUTTTTTTTTdUTTTTTTTT` (eight T). The model uses the table.
- 🟢 Results describes the guide, RT and PCR primers as "TSO minus triple G", but the
  table RT primer ends `…GAGTAC` and PCR primer `…GAGT`, not the `…GAGTACAT` left by
  removing `rGrG+G`. The exact table primers are drawn.

## Construct path

🟢 A barcoded oligo containing strategically placed deoxyuridines supports end extension over the RNA poly(A) tail. USER cleavage removes that oligo before reverse transcription. Template switching adds the opposite PCR handle, and full-length cDNA is converted to a PacBio SMRTbell. This preserves tail length and internal non-A residues in the same CCS read as the transcript isoform.

🟢 The TSO, RT primer, PCR primer, and all fifteen 2dU-BC extension oligos on the schematic are transcribed from Supplementary Table 1. BC3 is the member explicitly named for the bulk library in the defining paper's Methods; BC1–BC15 barcode individual oocytes.

🟢 **The TSO linker is `ACAT`, published by PAIso-seq itself.** Supplementary Table 1 row
`TSO` reads `AAGCAGTGGTATCAACGCAGAGTACATrGrG+G` verbatim, i.e. the 23-nt SMART handle, the
four-base linker `ACAT` and the `rGrG+G` tail. This matters because SMART-seq2 has two
published linkers — `AC` in Picelli et al. 2013 Supplementary Table 2 for the chosen
`rGrG+G` oligo, `ACAT` in the *Nature Protocols* form — and a protocol that only cites
SMART-seq2 cannot tell them apart. PAIso-seq prints its own, so no base is inherited by
assumption. `tools/chem.py` now owns those source-qualified ordering strings locally,
and `tools/selftest.py` compares the native TSO, RT primer, PCR primer and all15
modified guide members directly with the optional original TableS1. The historical
Smart-seq2 constructor imports no longer supply this protocol's chemistry.

🟡 The authors' 2022 detailed protocol is a later version: its Table 2 TSO adds
5′ `iCiGiC` and a 3′ `/NH2-C3/` block around the same core sequence. The
2019 defining Supplementary Table 1 does not print those modifications. This
page models the 2019 oligo; the later ordered oligo must not be substituted
for the original experiment.

🟢 Amplification uses the single published PCR primer. Template switching and reverse transcription place its common handle at opposite ends of the cDNA, so the same physical primer species binds both inward-facing sites. The sample barcode is carried by the extension oligo; it is not a second PCR primer.

## SMRTbell junction evidence

🟢 The 2019 PAIso-seq Online Methods, “PAIso-seq library construction,” identify size-selected amplified cDNA and an unversioned `SMRTbell Template Prep Kit`, followed by sequencing-primer annealing and Sequel sequencing. They do not specify dA tailing or the adapter's ligating end. The 2022 detailed protocol explicitly names `SMRTbell template prep kit 1.0-SPv3` (PacBio 100-991-900) in Materials and step 29; it also uses a modified TSO, so its kit choice cannot automatically be assigned to the 2019 experiment. The official PacBio 1.0 insert lists a Blunt Adapter, whereas the contemporaneous Express 2.0 overview describes T-overhang adapters and A/T ligation. 🟡 The model therefore retains only the source-supported closed SMRTbell topology and records the 2019 end-preparation and ligation geometry as unresolved; it does not draw an A/T or blunt junction for that experiment. Exact hairpin and sequencing-primer bases remain unverified.

## Additional source scope

🟢 Supplementary Figure 2 and Table 3 contain poly(A)-length calibration
spike-ins; Table 2 contains independent PAT validation primers. They are
separate from the library oligos in Table 1. The remaining supplementary
figures and Table 4 concern coverage, correlations and measured tail results.

🟢 The article's linked [author analysis repository](https://github.com/niehu2018/GV_oocyte_PAIsoSeqAnalysis/tree/ed419f406c35eb6812df10531348631892f85361)
prints the TSO/cDNA/poly(A)/barcode/RT-adapter read orientation in its README.
Its trimming code matches the supplied barcode-plus-adapter in either CCS
orientation; these read-side adapters are not SMRTbell hairpins.

The inspected repository revision has two execution inconsistencies: `run.pl`
requires `scripts/basic.pm20190608` while its tree contains `scripts/basic.pm`,
and the README trimming example includes an extra `sample` argument that the
script's argument check does not accept. These are author-code packaging/CLI
issues in the inspected revision, not demonstrated errors in the original
2019 chemical protocol or biological data. The schematic does not depend on
executing that command or guessing the missing module name.
