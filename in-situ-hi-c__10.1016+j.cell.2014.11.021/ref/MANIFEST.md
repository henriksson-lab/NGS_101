# in situ Hi-C — third-party reference material

Third-party files are deliberately not committed. The protocol model was checked against:

- Rao *et al.*, “A 3D map of the human genome at kilobase resolution reveals principles
  of chromatin looping,” *Cell* (2014),
  doi:[10.1016/j.cell.2014.11.021](https://doi.org/10.1016/j.cell.2014.11.021), PMCID
  **PMC5635824**. The main text defines the in situ Hi-C operations and MboI/DpnII option.
- The paper's Extended Experimental Procedures, section I.a.1. PMC distributes these in
  `NIHMS649130-supplement-Zipped.zip`; its interactive gate may require manual download.
- [ENCODE's faithful transcription of section I.a.1](https://www.encodeproject.org/documents/0dd3624f-d1f9-4505-b97c-da01b90a70f4/),
  which prints the MboI digest, biotin-14-dATP/Klenow fill-in, T4 ligation, 300–500-bp
  shearing, streptavidin selection, end repair, dA tailing and indexed-Illumina-adapter steps.

Run `python3 tools/get_sources.py "in situ Hi-C"`. Files are stored under
`$CHEM_DATA/sources/in-situ-hi-c__10.1016+j.cell.2014.11.021/`, never here.
