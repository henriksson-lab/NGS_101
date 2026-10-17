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
shearing and a 0.45X AMPure PB concentration, DNA Prep Enzyme (with NAD and 1:5-diluted
DNA Prep Additive, Reaction Mix 1) removes single-stranded overhangs for 15 min at 37 °C.
DNA Damage Repair Mix v2 is added for 30 min at 37 °C. End Prep Mix ("End-Repair/A-tailing")
acts for 10 min at 20 °C and 30 min at 65 °C, producing dA-tailed insert ends. All four
enzymatic steps are additions to one tube (procedure pp. 13–14).

🟡 5′ phosphorylation of the insert ends is inferred from ligation competence; the
procedure names the step "End-Repair/A-tailing" without stating it.

## 4. Close the SMRTbell

🟢 Overhang Adapter v3, Ligation Mix, Ligation Additive and Ligation Enhancer act for
60 min at 20 °C. The procedure names Overhang Adapter v3 without printing its ends.
The contemporaneous PacBio training document identifies a **T overhang** on this
adapter (PDF p. 7) and explicitly describes **A/T ligation** rather than blunt
ligation (PDF p. 9). The kit procedure's dA-tailed insert therefore pairs with the
reported adapter T end. The terminal strand polarity follows from the paired dA/dT
ends; the training slides do not draw a base-resolved adapter duplex.
🟢 The ligation (97 µL Reaction Mix 4) is purified with 44 µL (0.45X) AMPure PB beads,
eluted at 37 °C (31 µL if size selection follows, 11 µL otherwise; procedure pp. 15–16).
Optional BluePippin size selection (>15 kb or >30 kb) is followed by a 1X AMPure PB
cleanup (pp. 17–18). The library then goes to primer annealing and polymerase binding.

**Model correction (2026-10-17).** Earlier versions of this model and page drew a
post-ligation nuclease digestion of molecules with free ends. Procedure 101-693-800 v01
has no such step: ligation is followed directly by bead purification. The training
document's workflow comparison (PN 101-742-500 v01) places an ExoIII/ExoVII step only in
the older SMRTbell Template Prep Kit 1.0 column. The graph now draws the 0.45X AMPure PB
cleanup; bead purification does not select by topology, so the purified library is not
claimed to be free of incompletely ligated molecules.

🟡 The shared `Dumbbell` constructor guarantees the topology: insert strand 1 → right
hairpin → reverse complement of the insert → left hairpin → origin. The second insert
strand and repeated polymerase circuit are derived rather than transcribed twice.

🔴 The inspected sources do not disclose the complete Overhang Adapter v3 body or
sequencing-primer bases. The schematic shows the source-reported terminal T and known
hairpin topology without inventing the proprietary body sequence or borrowing
current-kit oligos.

## 5. Sequencing

🟢 The procedure (p. 18) routes primer annealing and polymerase binding to SMRT Link
Sample Setup, using Sequencing Primer v4 from the kit (diffusion loading only, not
MagBead loading). A PacBio sequencing primer and polymerase bind within a hairpin adapter. The
polymerase traverses the forward insert, opposite hairpin, reverse-complement insert and
first hairpin repeatedly, enabling circular consensus sequencing.

## 6. Open questions

- Complete Overhang Adapter v3 hairpin body and Sequencing Primer v4 bases (not printed in
  the procedure, training slides pp. 7, 9, 22 or product note).
- Contemporaneous SMRT Link v6 Sample Setup instructions (named by procedure p. 18) remain
  uninspected.
- The optional BluePippin size selection and 1X AMPure PB branch is described here but not
  drawn as a graph branch.
