# Source manifest

- Mahat et al. 2024: `python3 tools/get_sources.py --doi 10.1038/s41586-024-07517-7`.
- Zhou et al. 2004 (Mahat ref. 52), original Clinical Chemistry article PDF mirrored at [gene-quantification.de](https://gene-quantification.de/zhou-unlabeled-hrm-2004.pdf): `python3 tools/get_sources.py "scGRO-seq"`. The publisher route currently returns HTTP 403. Its discussion names a 3′ spacer C3 as a nonextendable terminator; the study's own probes used 3′ phosphate.
- Wu et al. 2018 Click-Code-Seq (Mahat ref. 50), original [ACS Figshare Supporting Information](https://acs.figshare.com/articles/journal_contribution/Nucleotide-Resolution_Genome-Wide_Mapping_of_Oxidative_DNA_Damage_by_Click-Code-Seq/6869855): `python3 tools/get_sources.py "scGRO-seq"`. The main ACS article remains purchase-gated; its SI contains the click chemistry methods and oligo table.
