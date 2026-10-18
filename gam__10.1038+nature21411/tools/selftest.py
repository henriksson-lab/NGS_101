"""Source transcription of WGA end lengths and selected NP index widths."""
from pathlib import Path
import csv,re,sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import data_root
import chem
c=Check();root=data_root()/"sources"/H.parent.name
source=Source(root/"nature21411_PMC5366070.xml.txt","defining article Methods library preparation")
if have(c,source,label="source WGA adapter extent"):
 match=re.search(r"Each library has (\d+)bp WGA adaptors at both ends",source.path.read_text())
 expected=int(match.group(1)) if match else None
 c("source WGA4 end lengths",[s.length_bp for s in chem.FINAL_LIBRARY if "WGA4 adapter" in s.name],[expected,expected])
table=Source(root/"nature21411_supp_NIHMS71147-supplement-Supplementary_Table_1.xlsx.txt","published SupplementaryTable1 NP index metadata")
if have(c,table,label="source selected HT NP index widths"):
 lines=table.path.read_text().splitlines();rows=list(csv.DictReader(lines[1:],delimiter='\t'))
 selected=[r for r in rows if r.get('Type')=='1NP' and r.get('Library prep kit')=='TruSeq HT']
 widths={len(r['Index 1 sequence']) for r in selected}|{len(r['Index 2 sequence']) for r in selected}
 c("selected HT NP index widths",widths,{8})
c.report()
