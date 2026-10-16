# PacBio Iso-Seq express 2.0 with SMRTbell prep kit 3.0

## Source scope

- 🟢 PacBio, [Preparing Kinnex-compatible cDNA using Iso-Seq express 2.0 kit](https://www.pacb.com/wp-content/uploads/Procedure-checklist-%E2%80%93-Preparing-Kinnex-compatible-cDNA-using-Iso-Seq-express-2.0-kit-1.pdf), 102-396-000 Rev10 (September 2026): source for the full-length cDNA front end. It explicitly offers SMRTbell prep kit 3.0 as a downstream route; it does not describe its adapter junction.
- 🟢 PacBio, [official Iso-Seq v2 primer FASTA](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta): exact bc01–bc12 cDNA PCR primer family and opposite amplification primer.
- 🟢 PacBio, [Iso-Seq library preparation using SMRTbell prep kit 3.0 training overview](https://www.pacb.com/wp-content/uploads/Iso-Seq-library-preparation-using-SMRTbell-prep-kit-3.0-training.pdf), workflow overview pages 5 and 12: repair and A-tailing followed by hairpin-adapter ligation, then nuclease treatment. This supports the downstream operation order.
- 🔴 The linked kit-specific [Iso-Seq SMRTbell prep kit 3.0 procedure](https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-Iso-Seq-libraries-using-SMRTbell-prep-kit-3.0.pdf) currently returns 404. The linked amplicon guide is a distinct sample class. Neither accessible source verifies the ligating hairpin overhang or complete adapter sequence; RT primer mix, TSO 2.0 and sequencing-primer bases also remain unresolved.

## Molecular path

Template switching creates full-length cDNA. The official FASTA establishes the two cDNA PCR primer footprints, including the selected sample barcode. The PacBio training overview then supports repair, dA tailing and SMRTbell adapter ligation. The final dumbbell drawing records closed topology, with inferred hairpin bodies and an explicitly unavailable ligation-end geometry. It does not claim a source-known T overhang.
