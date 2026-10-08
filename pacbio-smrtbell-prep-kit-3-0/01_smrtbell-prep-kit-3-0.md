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

🟢 The WGS workflow recommends shearing to 15–20 kb, followed by bead cleanup. Repair
buffer, End Repair Mix and DNA Repair Mix act for 30 min at 37 °C and 5 min at 65 °C,
performing damage repair, end repair and dA-tailing.

## 4. Dumbbell formation and selection

🟢 A T-overhang SMRTbell adapter is ligated for 30 min at 20 °C to both dA-tailed insert
ends. A subsequent nuclease treatment at 37 °C for 15 min removes unligated DNA and excess
adapter. Fully closed SMRTbells survive because they have no free DNA ends.

🟡 `Dumbbell` makes the topology an invariant: top insert → right hairpin → reverse
complement insert → left hairpin → origin. Forward and reverse-complement passes and the
covalently closed path are derived from the one insert rather than entered separately.

## 5. Primer, polymerase and repeated passes

🟢 In anneal/bind/cleanup, the standard sequencing primer is annealed for 15 min at room
temperature, followed by sequencing polymerase for 15 min. Both hairpin adapters carry
primer-binding sites. One productive polymerase follows the closed template through one
insert strand, the opposite hairpin, the other insert strand in reverse-complement
orientation, and back through the first hairpin for another circuit.

🔴 PacBio does not print the full current adapter, standard primer, enzyme-mix or
polymerase sequences/compositions. Historical adapters and primers are not substituted.
The schematic uses explicit proprietary placeholders but shows both binding sites and the
derived molecular trajectory.
