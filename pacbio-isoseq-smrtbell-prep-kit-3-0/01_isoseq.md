# PacBio Iso-Seq v2 with SMRTbell prep kit 3.0

## 1. What it is

The current PacBio Iso-Seq v2 workflow converts full-length RNA isoforms to amplified cDNA and then to closed SMRTbell templates for circular-consensus long-read sequencing.

## 2. Sources and evidence

- 🟢 PacBio, *Procedure & checklist — Preparing Kinnex-compatible cDNA using Iso-Seq express 2.0 kit*, document 102-396-000 Rev10 (September 2026): https://www.pacb.com/wp-content/uploads/Procedure-checklist-%E2%80%93-Preparing-Kinnex-compatible-cDNA-using-Iso-Seq-express-2.0-kit-1.pdf
- 🟢 PacBio official `IsoSeq_v2_primers_12.fasta`: https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta
- 🟢 PacBio, *Preparing multiplexed amplicon libraries using SMRTbell prep kit 3.0*, 102-359-000 Rev08 (August 2026): https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-multiplexed-amplicon-libraries-using-SMRTbell-prep-kit-3.0.pdf. This is the separate library-preparation procedure referenced by the Iso-Seq checklist.
- 🟢 The vendor workflow specifies input RNA QC, cDNA synthesis, cDNA amplification, repair/A-tailing, adapter ligation, nuclease treatment, primer annealing and polymerase binding.
- 🟢 Iso-Seq uses full-length cDNA and repeated traversal of a closed SMRTbell to generate a circular-consensus sequence.
- 🟢 The official FASTA publishes the complete bc01–bc12 PCR primers and the opposite cDNA amplification primer. Rev10 confirms that these two primer classes are paired during cDNA amplification.
- 🔴 The RT primer mix, TSO 2.0, SMRTbell hairpin, and separate sequencing-primer bases remain proprietary.
- 🟢 Rev08 confirms the downstream order: repair and dA tailing of the indexed cDNA, SMRTbell adapter ligation, nuclease cleanup, sequencing-primer annealing and polymerase binding. It names the adapter and standard sequencing primer, but gives neither base sequence.

## 3. Molecular path

Template switching creates an amplifiable full-length cDNA. Repair and dA tailing prepare the duplex for two hairpin ligations. The resulting molecule has no free DNA ends and carries the sequencing-primer site in its adapter.
