# ProteinSeq: High-Performance Proteomic Analyses by Proximity Ligation and Next Generation Sequencing

Defining DOI: 10.1371/journal.pone.0025583. Status: capped source pass; incomplete private model, no finished public schematic.

Sources actually inspected:

- https://pmc-oa-opendata.s3.amazonaws.com/PMC3183061.1/PMC3183061.1.xml: `journal.pone.0025583_PMC3183061.xml.txt` — Published main Methods194–274 and overview112–150; source qPCR versus NGS branches distinguished.
- https://pmc-oa-opendata.s3.amazonaws.com/PMC3183061.1/pone.0025583.s012.docx: `journal.pone.0025583_supp_pone.0025583.s012.docx.txt` — Published TableS8, complete table including Index1–12/reverse/enrichment primers.

🟢 Published main identifies three affinity recognition events: bead capture, followed by two antibody–DNA probes and ligation. DNA attachment is via source 5′ or 3′ thiol-bearing oligos; full conjugation sites and probe polymers are unassembled. Sliding-splint design intends cognate ligation; this graph does not predict affinity specificity, occupancy or ligation yield.

🟢 Source operational route: 1.5h sample capture → two magnetic washes → 1.5h probe binding → two washes → AmpLigase/NAD/splint50°C10min → one bead wash → NGS universal PCR35cycles →24µL aliquot/index program →10µL aliquot/3-cycle enrichment →GFX column →199bp2%E-gel selection →qPCR quantification →groups of12indexed samples →GAIIx single-end95bases with6bp sample index. Source literal thermal programs retained; extra indexing cycles are not inferred.

🟢 The earlier qPCR validation recipe uses dUTP and15universal cycles. The NGS paragraph specifies dNTPmix and35cycles; its dNTP composition is left unassembled rather than copying the qPCR dUTP assumption. Original bead/protein/conjugate species are not automatically discarded by PCR. PCR DNA copies are conceptually distinct from protein-linked input. Only source magnetic/column/gel separations assign removal.

🟢 TableS8 primers, verbatim:

- Index1 forward: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCCGTGATGCGAAACCTGGTCCGGTATC`
- Index reverse: `CAAGCAGAAGACGGCATACGAGATCGATTCGAGAACGTGACTGC`
- Enrich forward: `CAAGCAGAAGACGGCATACGAGAT`
- Enrich reverse: `AATGATACGGCGACCACCGAGATCT`

🟡 Canonical P5/P7 prefixes are imported and checked against the whole printed primers. The Index1 six-base window is a structured sample index; analyte identity derives from source paired PLA arms, not a new UMI. Exact whole reporter/ligation/read-primer geometry is unassembled.

Component lineage: main Methods explicitly cites earlier solid-support preparation and probe conjugation methods; those originals remain deferred, so no recipe inheritance/import is established here. ProteinSeq is PLA ligation, not Olink PEA extension.

Deferred targets (not claims of absence):

- Published TablesS5–S7 (probe bodies and sliding-splint sequence) and FigS1 original geometry; TableS3 source NGS probe concentrations
- Source-cited bead/antibody conjugation predecessors (main refs19/29), full manufacturing and workup
- Full universal primer sites, complete reporter copy/ligation and read geometry; actual sample/index assignment

Verification: mandatory standard oligo scraper was run over each fetched source bundle; source pointers and relevant modifications were checked against the originals stated above. Local tests check transcription and consequential fraction/order distinctions; graph validation is separate from molecular/source completeness. No absence audit or absence acceptance.

Physical QC branch: analytical aliquots remain separate from the production preparation. MALDI matrix-exposed, qPCR-amplified, gel-analysed or Qubit-assayed material is not automatically pooled for sequencing. Selected source quantities are retained; unreviewed analytical allocation amounts stay unspecified.
