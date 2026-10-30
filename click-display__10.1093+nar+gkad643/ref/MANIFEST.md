# Source retrieval and identity

Run `python3 tools/get_sources.py --slugs click-display__10.1093+nar+gkad643`. All originals/text twins remain in the external source store; canonical manifest URLs/hashes, not cache filenames alone, identify them.

- https://pmc-oa-opendata.s3.amazonaws.com/PMC10484664.1/PMC10484664.1.xml — Publishedmain Table1,Methods134–174 including mock-library selection and NGS; source recipe authority except corrected oligo groups/names
- https://pmc-oa-opendata.s3.amazonaws.com/PMC11662653.1/PMC11662653.1.xml — Mandatory published2024correction doi10.1093/nar/gkae1172; whole corrected O1/O2 table and corrected descriptions, originalp1 rendered /tmp/click18-correction.png; source-transcription checks use exact XML table
- https://pmc-oa-opendata.s3.amazonaws.com/PMC10484664.1/gkad643_supplemental_files.zip — PublishedoriginalSI TableS1 and FigS3 template legend; originaloligo1/2 names reversed, correction controls; complete assembly/primer landing unassembled

See [selected source scope and deferrals](../01_chem.md). No source files are checked in.
