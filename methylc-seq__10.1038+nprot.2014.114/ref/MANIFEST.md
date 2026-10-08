# Source manifest

- Urich et al. 2015, DOI `10.1038/nprot.2014.114`, PMC4465251; obtain with
  `python3 tools/get_sources.py --doi 10.1038/nprot.2014.114`.
- Bioo Scientific, NEXTflex Bisulfite-Seq product documentation, cat. 511912. The
  protocol names this exact product; the vendor documentation supplies the single 6-nt
  barcode architecture and NEXTflex PCR-primer sequences.
- Illumina, *Adapter Sequences*, document 1000000002694. Canonical Read 1, Index 1 and
  Read 2 primers are represented through `lib/illumina.py` and located on the modeled
  NEXTflex endpoint by `lib/seqprimers.py`.
