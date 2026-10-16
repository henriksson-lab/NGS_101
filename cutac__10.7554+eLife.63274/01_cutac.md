# CUTAC — source-owned molecular model

CUTAC uses antibody-tethered pA(G)–Tn5 under low-salt conditions to favor accessible DNA near transcription-associated nucleosomes. The antibody-to-accessible-site connection is spatial, not a covalent molecular bond.

## Original sources and scope

- 🟢 [Henikoff et al., eLife (2020)](https://doi.org/10.7554/eLife.63274), original Methods, identifies in-house pA–Tn5 and commercial Epicypher pAG–Tn5 separately. Low-salt options vary between experiments. Its post-tagmentation route is SDS release, Triton neutralization, gap filling/PCR and post-PCR SPRI cleanup.
- 🟢 Linked [in-house loading protocol 8yrhxv6](https://doi.org/10.17504/protocols.io.8yrhxv6), November 26, 2019, integer29425: original table prints the complete three loaded strands, including the MErev phosphate. The anonymous public API preserves whole table cells that are fragmented by PDF extraction.
- 🟢 Linked [CUT&Tag-direct with CUTAC Version2](https://doi.org/10.17504/protocols.io.bmbfk2jn), September18,2020, integer42055: materials and PCR sections specify custom Buenrostro primers and explicitly say **“Do not use Nextera or NEBNext primers.”** The PCR step permits universal **or** barcoded i5 with uniquely barcoded i7.
- 🟢 Cited [Buenrostro2015](https://doi.org/10.1038/nature14590) original MOESM36 workbook prints the custom indexed PCR families. Their source-specific arms overlap the in-house ME-A/B adapters. Shared canonical sequence names do not make the full custom PCR oligos interchangeable with a vendor kit's primers.

## Loaded oligos: in-house branch only

All strands are written 5′→3′, with source modification retained:

ME-A: `{{= cutac.ME_A }}`

ME-B: `{{= cutac.ME_B }}`

MErev: `/5Phos/{{= cutac.ME_REVERSE }}`

These source-known oligos enter the graph as an annealed pool. The original publication's commercial pAG–Tn5 input is a distinct reagent; this in-house table does not assign its loaded molecular contents or historical kit version.

## Conditional assembly and validation

🟡 The represented branch combines the cited in-house adapter pool with the permitted **barcoded-i5** custom PCR family. The selected i5/i7 members and their sample assignment are unassigned. The universal-i5 alternative is not silently converted into a variable sample barcode.

🟡 A/B, A/A and B/B fragment-end assemblies are inferred possibilities, not measured product frequencies or asserted dimer pairings. Same-end products remain explicit unresolved terminals; no absolute amplification impossibility or physical discard is asserted.

The graph follows the CUTAC release/neutralization route before gap filling and PCR. It does not import the earlier CUT&Tag extraction sequence. Exact inward binding of the declared custom PCR arms is recomputed on the conditional A/B template; this validates that assembly, not its selection in every experiment.

The article reports paired-end HiSeq2500 Turbo sequencing. The model does not assign exact loaded sequencing-primer oligos, acquired cycle counts or index orientation from PCR-site architecture alone.

## Remaining source work

The audit remains in progress. The original dataset workbook supplies metadata and names additional extraction-protocol files; those names need source-identity checks. The transparent-reporting form points to the analysis tutorial. Linked nuclei-preparation and precursor methods also need their bounded original-source checks. A commercial reagent's unassigned contents and selected indexes are model questions, not invented unavailable specification documents or confirmed publication-wide absence.
