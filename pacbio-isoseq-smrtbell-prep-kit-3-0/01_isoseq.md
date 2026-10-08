# PacBio Iso-Seq v2 with SMRTbell prep kit 3.0

## 1. What it is

The current PacBio Iso-Seq v2 workflow converts full-length RNA isoforms to amplified cDNA and then to closed SMRTbell templates for circular-consensus long-read sequencing.

## 2. Sources and evidence

- 🟢 PacBio, *Procedure & checklist — Preparing Iso-Seq v2 libraries using SMRTbell prep kit 3.0*, document 102-396-000 Rev07 (June 2025): https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-Iso-Seq-libraries-using-SMRTbell-prep-kit-3.0.pdf
- 🟢 The vendor workflow specifies input RNA QC, cDNA synthesis, cDNA amplification, repair/A-tailing, adapter ligation, nuclease treatment, primer annealing and polymerase binding.
- 🟢 Iso-Seq uses full-length cDNA and repeated traversal of a closed SMRTbell to generate a circular-consensus sequence.
- 🔴 Current cDNA oligo and SMRTbell hairpin bases are proprietary.

## 3. Molecular path

Template switching creates an amplifiable full-length cDNA. Repair and dA tailing prepare the duplex for two hairpin ligations. The resulting molecule has no free DNA ends and carries the sequencing-primer site in its adapter.
