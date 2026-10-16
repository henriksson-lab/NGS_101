#!/usr/bin/env python3
"""Guard Stereo-seq's local oligos against its actual author protocol."""
from pathlib import Path
from html.parser import HTMLParser
import re
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from storage import data_root
import stereo_seq as model
from checks import Check, Source, have

class Rows(HTMLParser):
    def __init__(self):
        super().__init__(); self.rows=[]; self.row=None; self.cell=None
    def handle_starttag(self, tag, attrs):
        if tag == 'tr': self.row=[]
        elif tag in ('td', 'th') and self.row is not None: self.cell=[]
    def handle_data(self, data):
        if self.cell is not None: self.cell.append(data)
    def handle_endtag(self, tag):
        if tag in ('td', 'th') and self.cell is not None:
            self.row.append(''.join(self.cell).strip()); self.cell=None
        elif tag == 'tr' and self.row is not None:
            self.rows.append(self.row); self.row=None

check=Check()
data=data_root()
source=Source(data/'sources'/HERE.parent.name/'stereo.seq.html',
              'Original MOSTA author-linked protocol https://db.cngb.org/stomics/assets/html/stereo.seq.html')
if have(check, source, label='Original MOSTA protocol and primer table'):
    parser=Rows(); parser.feed(source.path.read_text())
    rows={re.sub(r'\s+', ' ', r[0]): re.sub(r'\s+', '', r[1])
          for r in parser.rows if len(r)==2}
    selected={'DNB library oligo1': model.CID_PRECURSOR,
              'Stereo-TSO': model.TSO,
              'cDNA PCR primer': model.CDNA_PCR_PRIMER,
              'Stereo-seq Library F-primer': '/5phos/'+model.LIBRARY_F[:-2]+'*C*A',
              'Stereo-seq Library R-primer': model.LIBRARY_R,
              'Stereo-seq DNB splint oligo': model.LIBRARY_SPLINT,
              'Stereo-seq sequencing primer1': model.READ1,
              'Stereo-seq sequencing primer2': model.READ2,
              'Stereo-seq MDA primer': model.MDA_PRIMER}
    check('local primer declarations preserve original author table and modifications',
          selected, {k:rows[k] for k in selected})
    flat=re.sub(r'\s+', '', source.path.read_text())
    check('short ligation input is the source phosphorylated capture oligo',
          '/5phos/'+model.CAPTURE_OLIGO in flat)
check.report()
