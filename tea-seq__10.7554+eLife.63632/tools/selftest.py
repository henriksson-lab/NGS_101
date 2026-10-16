#!/usr/bin/env python3
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import tea_seq as M
from checks import Check, Source, have
c=Check()
bundle=data_root()/"sources"/"tea-seq__10.7554+eLife.63632"
guide=Source(bundle/"CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevA.pdf.txt",M.MULTIOME_GUIDE)
if have(c,guide,label="complete source-owned Multiome endpoints"):
    text=re.sub(r"\s+","",guide.path.read_text())
    products=re.findall(r"5[′’ʹ']-([A-Za-z0-9_-]+)-3[′’ʹ']",text)
    expanded=[]
    for product in products:
        product=product.replace("cDNA_Insert","X"*40)
        product=re.sub(r"-+insert-+","X"*42,product)
        product=re.sub(r"N(\d+)",lambda m:"N"*int(m[1]),product)
        expanded.append(product.replace("-",""))
    c("GEX complete molecule matches the original guide",M.RNA.top() in expanded)
    c("ATAC complete molecule matches the original guide",M.ATAC.top() in expanded)
article=Source(bundle/"eLife.63632_PMC8034981.xml.txt",M.SOURCE)
if have(c,article,label="TEA-specific demultiplexing masks"):
    text=article.path.read_text()
    masks=list(dict.fromkeys(re.findall(r"Y\d+n\*,I\d+(?:n\*)?,(?:I\d+n\*|Y\d+|n\*),Y\d+n\*",text)))
    c("all three source masks located",len(masks),3)
    if len(masks)==3:
        for heading,mask in zip(("RNA library","ATAC library","Antibody-tag library"),masks):
            source_cycles=tuple(int(re.search(r"\d+",x)[0]) if re.search(r"\d+",x) else None for x in mask.split(","))
            c(heading+" follows its own published mask",tuple(M.READ_LENGTHS[heading].get(role) for role in M.sp.ROLES),source_cycles)
table=Source(bundle/"eLife.63632_elife-63632-supp5-v2.xlsx.txt","TEA Supplementary file 5")
if have(c,table,label="whole ADT amplification primers from original table"):
    rows={line.split("\t")[0]:line.split("\t")[1] for line in table.path.read_text().splitlines() if "\t" in line}
    def bases(name):
        printed=rows[name].split("-",1)[1].rsplit("-",1)[0]
        return re.sub(r"[^ACGTN]","",printed)
    c("whole SI-PCR-Oligo matches Table 5",M.ADT_P5_PRIMER,bases("SI-PCR-Oligo"))
    c("whole ADT-i7 primer matches Table 5",M.il.P7+"N"*8+M.ADT_I7_ARM,bases("ADT-i7 Primer"))
    c("whole ADT-Rev-AMP matches Table 5",M.ADT_CAPTURE[:-1],bases("ADT-Rev-AMP"))

c.report()
