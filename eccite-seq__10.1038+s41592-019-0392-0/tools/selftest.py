#!/usr/bin/env python3
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import eccite_seq
from checks import Check, Source, have
from multimodal_page import render
check=Check()
render(eccite_seq)
print("PASS  transcript, antibody, hashtag and direct-guide libraries build")
source_dir=data_root()/"sources"/H.parent.name
table=Source(source_dir/"s41592-019-0392-0_41592_2019_392_MOESM1_ESM.pdf.txt",
             "Mimitou et al. 2019 original Springer Supplementary Tables 3–5")
if have(check,table,label="ECCITE-seq ordered conjugate barcode source"):
    source=table.path.read_text(encoding="utf-8")
    antibody=source.split("Supplementary Table 3: ECCITE-seq antibody-oligo conjugates",1)[1].split("Supplementary Table 4:",1)[0]
    hashtags=source.split("Supplementary Table 4:",1)[1].split("Supplementary Table 5:",1)[0]
    for name,sequence,section in (
        ("15-base antibody",eccite_seq.ANTIBODY_LONG_BARCODE,antibody),
        ("12-base antibody",eccite_seq.ANTIBODY_SHORT_BARCODE,antibody),
        ("mouse Hashtag_A",eccite_seq.HASHTAG_A_BARCODE,hashtags)):
        check(f"{name} barcode is printed as a whole source row",
              sequence in section.splitlines())
    barcode_rows=[line for line in antibody.splitlines() if re.fullmatch("[ACGT]{12,15}",line)]
    check("Table 3 antibody barcode lengths",(sum(map(lambda s:len(s)==15,barcode_rows)),
                                              sum(map(lambda s:len(s)==12,barcode_rows))),(49,3))
protocol=Source(source_dir/"s41592-019-0392-0_41592_2019_392_MOESM3_ESM.pdf.txt",
                "Mimitou et al. 2019 MOESM3 supplementary protocol oligonucleotide list")
if have(check,protocol,label="ECCITE-seq supplementary protocol oligo list"):
    text=re.sub(r"\s+","",protocol.path.read_text(encoding="utf-8")).replace("*","")
    x8="x"*8
    for name,oligo in (
        ("protein-tag example",eccite_seq.PROTEIN_TAG_HANDLE+x8+eccite_seq.TAG_TSO_END),
        ("hashtag example",eccite_seq.HASHTAG_HANDLE+x8+eccite_seq.TAG_TSO_END),
        ("RPI1",eccite_seq.RPI1_P7+x8+eccite_seq.RPI1_HANDLE),
        ("D701_s",eccite_seq.D701S_P7+x8+eccite_seq.D701S_HANDLE),
        ("Next_nst_x",eccite_seq.NEXT_NST_P7+x8+eccite_seq.NEXT_NST_HANDLE+eccite_seq.NEXT_NST_SCAFFOLD),
        ("SI-PCR",eccite_seq.SI_PCR)):
        check(f"{name} bases transcribed from MOESM3",oligo in text)
check.report()
