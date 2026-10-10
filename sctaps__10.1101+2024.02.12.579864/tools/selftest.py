#!/usr/bin/env python3
import build_page
from checks import Check
c=Check(); c("published TAPS conversion paths",
             tuple(build_page.M.TAPS_PATHS[x].states for x in ("5mC","5hmC","C")),
             (("5mC","5caC","DHU","T"),("5hmC","5caC","DHU","T"),("C","C")))
build_page.render(build_page.M); c.report()
