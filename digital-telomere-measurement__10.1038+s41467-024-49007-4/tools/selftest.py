#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dtm as M
import storage
from checks import Check,Source,have,run_common
c=Check(); c.section("DTM source transcription")
c("seqTether",M.SEQ_TETHER,"AACCTTGGAGATGCACGGAGCAAGCAAT")
c("capture oligo",M.CAPTURE,"TGCTCCGTGCATCTCCAAGGTTCACAAAGACACCGACAACTTTCTTCCTAACC")

# The arm/barcode/telomere-arm partition is not printed; derive it from the
# Supplementary Data 3 rows (shared prefix/suffix of the eight canonical oligos and
# the shared prefix of the TSQ variants with p-t3-nb01) and compare with the model.
SD3=Source(storage.sources_root()/"digital-telomere-measurement__10.1038+s41467-024-49007-4"
           /"s41467-024-49007-4_supp_41467_2024_49007_MOESM6_ESM.xlsx.txt",
           "Nature Communications 15:5148 Supplementary Data 3 (MOESM6); tools/get_sources.py")
if have(c,SD3,label="Supplementary Data 3 partition"):
    rows=dict(l.rstrip("\n").split("\t")[:2] for l in SD3.path.read_text().splitlines()
              if l.startswith("p-"))
    canon=[v for k,v in rows.items() if k.startswith("p-t3-")]
    tsq=[v for k,v in rows.items() if k.startswith("p-tsq")]
    def prefix(xs):
        i=0
        while all(x[i]==xs[0][i] for x in xs): i+=1
        return i
    def suffix(xs):
        i=0
        while all(x[-1-i]==xs[0][-1-i] for x in xs): i+=1
        return i
    c("p-t3-nb01 in Supplementary Data 3",rows.get("p-t3-nb01"),M.CAPTURE)
    c("tether arm = shared canonical prefix",len(M.CAPTURE_TETHER_ARM),prefix(canon))
    c("telomere arm = shared canonical suffix",len(M.CAPTURE_TELOMERE_ARM),suffix(canon))
    c("TSQ variants change only the telomere arm",
      len(M.CAPTURE_TETHER_ARM)+len(M.CAPTURE_BARCODE),prefix(tsq+[M.CAPTURE]))
run_common(c); c.report()
