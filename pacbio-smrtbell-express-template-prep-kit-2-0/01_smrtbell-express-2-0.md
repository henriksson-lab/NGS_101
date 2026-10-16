# PacBio SMRTbell Express Template Prep Kit 2.0

## 1. What it is

PacBio's historical large-insert kit repairs mechanically sheared genomic DNA, removes
single-stranded overhangs, dA-tails both duplex ends, and ligates a hairpin adapter to
each end. The product is a covalently closed SMRTbell: one continuous strand folded into
a dumbbell so a polymerase can read both insert strands repeatedly.

Commercial defining source: PacBio “Preparing gDNA Libraries Using the SMRTbell Express
Template Preparation Kit 2.0,” procedure **101-693-800 Version 01 (January 2019)**.

## 2. Sources read

- 🟢 PacBio procedure 101-693-800 Version 01, archived with the publication supplement linked in the manifest.
- 🟢 PacBio customer training PN 101-742-500 Version 01 (February 2019), [large-insert library preparation](https://www.pacb.com/wp-content/uploads/SMRTbell-Express-TPK-2.0-Overview-Large-Insert-gDNA-Library-Prep-Customer-Training.pdf), PDF pp. 7 and 9.
- 🟢 PacBio product note, “SMRTbell express template prep 2.0 for large-insert libraries.”
- 🟢 Travers *et al.*, *Nucleic Acids Research* (2010), doi:[10.1093/nar/gkq543](https://doi.org/10.1093/nar/gkq543), for the published closed SMRTbell topology.

## 3. Insert preparation

🟢 The guide covers 10 kb, greater-than-15 kb and greater-than-30 kb libraries. After
shearing, DNA Prep Enzyme removes single-stranded overhangs for 15 min at 37 °C. Damage
repair follows for 30 min at 37 °C. End Prep Mix acts for 10 min at 20 °C and 30 min at
65 °C, producing repaired, phosphorylated, dA-tailed insert ends.

## 4. Close the SMRTbell

🟢 Overhang Adapter v3, Ligation Mix, Ligation Additive and Ligation Enhancer act for
60 min at 20 °C. The procedure names Overhang Adapter v3 without printing its ends.
The contemporaneous PacBio training document identifies a **T overhang** on this
adapter (PDF p. 7) and explicitly describes **A/T ligation** rather than blunt
ligation (PDF p. 9). The kit procedure's dA-tailed insert therefore pairs with the
reported adapter T end. The terminal strand polarity follows from the paired dA/dT
ends; the training slides do not draw a base-resolved adapter duplex.
After bead purification, nuclease removes linear or incompletely ligated molecules.

🟡 The shared `Dumbbell` constructor guarantees the topology: insert strand 1 → right
hairpin → reverse complement of the insert → left hairpin → origin. The second insert
strand and repeated polymerase circuit are derived rather than transcribed twice.

🔴 The inspected sources do not disclose the complete Overhang Adapter v3 body or
sequencing-primer bases. The schematic shows the source-reported terminal T and known
hairpin topology without inventing the proprietary body sequence or borrowing
current-kit oligos.

## 5. Sequencing

🟢 A PacBio sequencing primer and polymerase bind within a hairpin adapter. The
polymerase traverses the forward insert, opposite hairpin, reverse-complement insert and
first hairpin repeatedly, enabling circular consensus sequencing.
