# PacBio SMRTbell prep kit 3.0

## 1. What it is

PacBio's PCR-free WGS and metagenome preparation caps both ends of a repaired long DNA
duplex with hairpin adapters. The resulting SMRTbell is chemically one covalently closed
single strand folded into a dumbbell, allowing a polymerase to traverse both insert strands
repeatedly.

Commercial defining source: PacBio, “Preparing whole genome and metagenome libraries using
SMRTbell prep kit 3.0,” procedure **102-166-600 REV09 AUG2026**.

## 2. Sources read

- 🟢 PacBio procedure 102-166-600 Rev09, August 2026.
- 🟢 PacBio SMRTbell prep kit 3.0 training, 102-390-900 Rev03, November 2025.
- 🟢 Travers *et al.*, “A flexible and efficient template format for circular consensus
  sequencing and SNP detection,” *Nucleic Acids Research* (2010),
  doi:[10.1093/nar/gkq543](https://doi.org/10.1093/nar/gkq543), for the published SMRTbell
  topology and polymerase path.

## 3. Insert preparation

🟢 The WGS workflow recommends shearing to 15–20 kb (optionally after Short Read
Eliminator on gDNA), followed by a 1X SMRTbell bead cleanup. Repair buffer (8 µL), End
repair mix (4 µL) and DNA repair mix (2 µL) act for 30 min at 37 °C and 5 min at 65 °C
("Repair and A-tailing", procedure section 4). Training p. 6 draws a 3′ A on both insert
strands.

🟡 5′ phosphorylation is inferred from ligation competence; the procedure says only that
the step "repairs sites of DNA damage and prepares the sheared DNA for ligation".

## 4. Dumbbell formation and selection

🟢 The SMRTbell adapter (or, optionally, 4 µL of an adapter from SMRTbell adapter index
plate 96A, replacing it) is ligated with Ligation mix and Ligation enhancer for 30 min at
20 °C (section 5). The ligation is purified with 95 µL (1X) SMRTbell cleanup beads
(steps 5.9–5.22). Nuclease treatment (Nuclease buffer + Nuclease mix, 37 °C for 15 min,
section 6) "removes unligated DNA fragments and leftover SMRTbell adapter". A 3.1X
(155 µL) cleanup with 35% v/v AMPure PB beads then removes the nucleases and progressively
depletes fragments shorter than 5 kb (section 7, step 6.8); libraries <10 kb use a 1X
SMRTbell bead cleanup instead.

🟡 The adapter's 3′ T is inferred: the procedure and training establish a dA-tailed
insert and name the adapter design "Overhang" (training, Revio/Vega application table),
but neither prints the adapter terminal base. Fully closed SMRTbells survive nuclease
treatment because they have no free DNA ends (derived from the topology).

**Model correction (2026-10-17).** The graph previously went directly from ligation to
nuclease and called the nuclease product the final library. It now draws the post-ligation
1X SMRTbell bead cleanup and the post-nuclease 3.1X 35% AMPure PB size selection, and
labels the adapter T as inferred rather than source-stated.

🟡 `Dumbbell` makes the topology an invariant: top insert → right hairpin → reverse
complement insert → left hairpin → origin. Forward and reverse-complement passes and the
covalently closed path are derived from the one insert rather than entered separately.

## 5. Primer, polymerase and repeated passes

🟢 In anneal/bind/cleanup (ABC, section 8, Revio/Vega), the standard sequencing primer is
annealed for 15 min at room temperature, followed by sequencing polymerase for 15 min and
a 100 µL SMRTbell bead cleanup; Sequel II/e libraries go to SMRT Link Sample Setup. Both hairpin adapters carry
primer-binding sites. One productive polymerase follows the closed template through one
insert strand, the opposite hairpin, the other insert strand in reverse-complement
orientation, and back through the first hairpin for another circuit.

🔴 PacBio does not print the full current adapter, standard primer, enzyme-mix or
polymerase sequences/compositions. Historical adapters and primers are not substituted.
The schematic uses explicit proprietary placeholders but shows both binding sites and the
derived molecular trajectory.

## 6. Open questions

- Adapter terminal base and body, adapter index plate 96A–D sequences and standard
  sequencing-primer bases (index plate documents not fetched).
- Optional gel-based size selection, <10 kb 1X cleanup and Short Read Eliminator branches
  are described but not drawn in the graph.
