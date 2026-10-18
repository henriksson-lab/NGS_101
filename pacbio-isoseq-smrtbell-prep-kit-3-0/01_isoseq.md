# PacBio Iso-Seq express 2.0 with SMRTbell prep kit 3.0

## Source scope

- 🟢 PacBio, [Preparing Kinnex-compatible cDNA using Iso-Seq express 2.0 kit](https://www.pacb.com/wp-content/uploads/Procedure-checklist-%E2%80%93-Preparing-Kinnex-compatible-cDNA-using-Iso-Seq-express-2.0-kit-1.pdf), 102-396-000 Rev10 (September 2026): source for the full-length cDNA front end. Pages 1 and 13 send monomer SMRTbell libraries to the multiplexed amplicon procedure below. Its revision history lists "Remove SPK3 library preparation" at Rev09 (July 2026), so revisions up to Rev08 carried their own SMRTbell 3.0 section.
- 🟢 PacBio, [official Iso-Seq v2 primer FASTA](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta): exact bc01–bc12 cDNA PCR primer family and opposite amplification primer.
- 🟢 PacBio, [Preparing multiplexed amplicon libraries using SMRTbell prep kit 3.0](https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-multiplexed-amplicon-libraries-using-SMRTbell-prep-kit-3.0.pdf), 102-359-000 Rev08 (August 2026), sections 1–5: the downstream monomer route named by Rev10. It gives the repair & A-tailing, adapter ligation, nuclease and annealing reagents and conditions, but prints no adapter or primer bases.
- 🟡 PacBio, [Iso-Seq library preparation using SMRTbell prep kit 3.0 training overview](https://www.pacb.com/wp-content/uploads/Iso-Seq-library-preparation-using-SMRTbell-prep-kit-3.0-training.pdf), PN 102-393-400 v01 (April 2022). It describes the earlier Iso-Seq express workflow (NEBNext Single Cell/Low Input module with oligo-dT30VN RT primer, PacBio express oligo kit 101-737-500, 24 barcoded primer pairs). This is a different cDNA version from express 2.0 and is not used for oligo identities. Its workflow page agrees with the downstream order: repair/A-tailing, ligation, nuclease, each with 1.3X cleanups.
- 🔴 The Iso-Seq-specific 102-396-000 Rev07 (June 2025) PDF URL returns 404 and was not read. No inspected source prints the SMRTbell adapter bases or ligating end, the RT primer mix, TSO 2.0 or the standard sequencing primer.

## Molecular path

From 102-396-000 Rev10, sections 2–3:

1. Mix 300 ng total RNA (RIN ≥7) with the Iso-Seq RT primer mix, then 70 °C for 5 min and hold at 20 °C.
2. Add the RT buffer and RT enzyme mix and incubate 42 °C for 45 min.
3. Add 2 µL Iso-Seq Express TSO 2.0 and incubate 42 °C for 15 min. Template switching is a separate addition after first-strand synthesis.
4. Clean up with 1.3X SMRTbell beads and elute in 21 µL.
5. Run PCR with the Iso-Seq cDNA PCR mix, the cDNA amplification primer and one bc01–bc12 primer: 98 °C 45 s; **12 cycles** of 98 °C 10 s, 60 °C 15 s, 72 °C 3 min for SMRTbell prep kit 3.0 (10 cycles for Kinnex); then 72 °C 5 min.
6. Clean up with 1.1X beads (0.9X for longer cDNA) and elute in 47 µL. Proceed with at least 100 ng. Barcoded cDNAs may be pooled by equal mass up to 500 ng.

From 102-359-000 Rev08, sections 1–4:

1. Clean up the input with 1.3X beads.
2. Repair & A-tailing (Repair buffer, End repair mix, DNA repair mix): 37 °C 30 min, then 65 °C 5 min.
3. Adapter ligation (SMRTbell adapter, Ligation mix, Ligation enhancer): 20 °C 30 min, then a 1.3X cleanup.
4. Nuclease treatment: 37 °C 15 min, then a 1.3X cleanup.
5. Section 5 anneals the Standard sequencing primer before polymerase binding.

The official FASTA establishes complete PCR oligo identities, including the
selected barcode; it does not by itself establish their initial footprints on
the unbarcoded RT/template-switch product. The final dumbbell drawing records closed topology, with inferred hairpin bodies and an explicitly unavailable ligation-end geometry. It does not claim a source-known T overhang: the amplicon procedure does not state one, and the sister WGS kit page is a different procedure.

## Open questions

- SMRTbell adapter bases and ligating end, and the standard sequencing primer bases. These are not printed in either inspected procedure; the Rev07 Iso-Seq SMRTbell section is unread.
- RT primer mix and TSO 2.0 bases. They are not printed in Rev10. The oligo-dT30VN in the 2022 training belongs to the earlier NEB-based workflow.
- The substrate specificity of the nuclease mix. The procedures say only "Nuclease treatment".

## Local chemistry ownership migration — 2026-10-18

Removed the unused shared template-switch illustration. Its fixed14-A/14-T
tracts and invented TSO-terminal bases are not used to fill the current RT
primer or TSO2.0 gaps. Native RT/template-switch states instead declare the
source reagent roles and variable-length RNA/cDNA regions. PCR consumes the
actual unbarcoded single-stranded precursor; its binding geometry remains an
explicit implementation boundary rather than using the already barcoded PCR
product as the template. bc01 is the selected representative primer, while all
twelve alternatives and the opposite primer remain source-listed and checked
against the complete official FASTA.

The local polymerase traversal preserves the sample-barcode feature and unknown
insert/hairpin lengths. Drawing widths are not molecular lengths. Unknown
hairpin end chemistry is now set structurally on the dumbbell, preventing the
renderer from labelling it a direct repaired junction. The inspected Rev08
procedure's1.3X amplicon-input cleanup now precedes repair/A-tailing, distinct
from the preceding Rev10 PCR cleanup. The selected closed topology is a
candidate shown from a ligation mixture; nuclease specificity and other
products are unresolved, and hydrolysis is not equated with physical discard.

Re-read three originals: printed102-396-000 Rev10 SEP2026 RT/template-switch
and PCR sections2–3; official IsoSeq_v2_primers_12.fasta; printed102-359-000
Rev08 AUG2026 input-cleanup and nuclease/cleanup sections1 and4. These are the
already selected, directly linked modern monomer route. The older2022 NEB
workflow and unavailable Rev07 are not substituted for it. The bundle scraper
was run; the original source passages and sequence records were read directly.
This migration does not establish absent publication information.
