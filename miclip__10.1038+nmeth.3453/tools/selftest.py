#!/usr/bin/env python3
import re
import build_page
import storage
from checks import Check, Source, have

M = build_page.M
c = Check(); c("published inline identifier pattern", M.INLINE_LAYOUT, (3, 4, 2))

# Source transcription: the inferred 2017 oligos must occur verbatim in Grozhik 2017 §2.3.
grozhik = Source(storage.sources_root() / "miclip__10.1038+nmeth.3453" / "PMC5562447_PMC5562447.xml.txt",
                 "python3 tools/get_sources.py --slugs miclip__10.1038+nmeth.3453 (pmccloud:PMC5562447)")
if have(c, grozhik):
    text = grozhik.path.read_text(encoding="utf-8")
    for label, printed in (("L3 linker", f"/5rApp/{M.L3_LINKER}/3ddC/"),
                           ("cut oligo", f"{M.CUT_OLIGO}/3ddC/"),
                           ("P5 Solexa primer", M.P5_SOLEXA), ("P3 Solexa primer", M.P3_SOLEXA)):
        c(f"Grozhik 2017 prints {label}", printed in text)
    rt = re.findall(r"RT(\d+): /5phos/N[Nn]([acgt]{4})NNN" + M.RT_FIXED[:19] + "gatc" + M.RT_FIXED[23:], text)
    c("Grozhik 2017 RT1–RT16 barcodes", tuple(b for _, b in sorted(rt, key=lambda x: int(x[0]))), M.RT_BARCODES)

build_page.render(M); c.report()
