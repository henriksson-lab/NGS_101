#!/usr/bin/env python3
from pathlib import Path
import os
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import kinnex as K
from checks import Check, Source, have
c=Check(); c.section("Kinnex source transcription")
c("parallel orientation-specific PCR reactions", K.KINNEX_PCR_REACTIONS, 8)
c("array-formation incubation", K.ARRAY_FORMATION, (45, 60))
c("post-array repair incubation", K.DNA_REPAIR, (45, 30))
source_dir=Path(os.environ.get("CHEM_DATA", H.parents[1]/"_data"))/"sources"/"pacbio-kinnex-full-length-rna"
for name, declared, url in (
    ("mas8_primers.fasta", K.MAS8_REFERENCE_MOTIFS,
     "https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-MAS_adapters/MAS-Seq_Adapter_v3/mas8_primers.fasta"),
    ("Kinnex_barcoded_adapter_indexes.fasta__Kinnex_barcoded_adapter_indexes.fasta",
     K.KINNEX_TERMINAL_BARCODES,
     "https://www.pacb.com/wp-content/uploads/Kinnex_barcoded_adapter_indexes.fasta_.zip"),
):
    source=Source(source_dir/name,url)
    if have(c,source,label=f"source transcription of {name}"):
        lines=source.path.read_text().splitlines()
        c(f"{name} original FASTA agrees with displayed sequences",
          dict(zip((s[1:] for s in lines[::2]), lines[1::2])),declared)
guide=Source(source_dir/"SMRT-Link-User-Guide-v13.1.pdf.txt",
             "https://www.pacb.com/wp-content/uploads/SMRT-Link-User-Guide-v13.1.pdf")
if have(c,guide,label="Iso-Seq v2 trimming reference transcription"):
    guide_text=guide.path.read_text()
    c("Iso-Seq v2 bc01 trimming reference matches original guide",
      f">IsoSeqX_bc01_5p\n{K.ISOSEQ_V2_BC01_5P_REFERENCE}" in guide_text)
    c("Iso-Seq v2 reverse-complemented 3-prime trimming reference matches original guide",
      f">IsoSeqX_3p\n{K.ISOSEQ_V2_3P_REFERENCE}" in guide_text)
c.report()
