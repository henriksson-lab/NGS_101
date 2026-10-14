# TUC-seq

**Defining paper:** Riml et al., *Angewandte Chemie International Edition* (2017), doi:10.1002/anie.201707465.

## Molecular boundary

🟢 RNA synthesized during a 4sU pulse contains thiouridine. Osmium-tetroxide-mediated chemistry converts 4sU into a cytidine derivative before reverse transcription.

🟢 Newly synthesized molecules are recognized by U-to-C sequence changes distributed through the read. The conversion mark is molecular content, not an indexing sequence.

## Evidence boundary

🟡 TUC-seq defines the conversion and its sequencing readout rather than one immutable adapter kit. Kit-dependent ends are labelled by role.

## Oligo source audit

🟢 The [publisher record](https://onlinelibrary.wiley.com/doi/full/10.1002/anie.201707465) establishes osmium-mediated 4sU-to-C conversion and targeted detection of seven mammalian mRNAs, and lists `anie201707465-sup-0001-misc_information.pdf` as the detailed supplement. 🔴 The corresponding publisher supplement URL returns HTTP 403 through `get_sources.py`; its methods and any gene-specific RT/PCR/sequencing primers remain unread. The model's generic indexed RNA-seq ends are only role placeholders, and a source-verified TUC-seq oligo inventory is still blocked on that supplement.
