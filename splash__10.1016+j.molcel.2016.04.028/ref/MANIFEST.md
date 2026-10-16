# Source manifest

- Defining article: https://doi.org/10.1016/j.molcel.2016.04.028
- Detailed author protocol: https://doi.org/10.3791/55255. The retrieved JoVE PDF covers RNA end preparation, linker ligation, reverse transcription, ssDNA circularization, and PCR; the web article's materials list prints the modified adapter and RT primer.
- The author materials list names NEB E7335L as the final PCR primer kit. The vendor manual is https://www.neb.com/en/-/media/nebus/files/manuals/manuale7335_e7500_-e7710_e7730.pdf; the exact SPLASH circle-to-amplicon geometry still needs derivation.
- Retrieve with `python3 tools/get_sources.py --slugs splash__10.1016+j.molcel.2016.04.028`.

- Original 2016 supplemental information: https://ars.els-cdn.com/content/image/1-s2.0-S1097276516301046-mmc1.pdf (Figure S2F, Supplemental Experimental Procedures). Additional tables have repeatable candidate routes in `catalogue/extra_sources.tsv`; downloaded sources and reading status are tracked separately.
