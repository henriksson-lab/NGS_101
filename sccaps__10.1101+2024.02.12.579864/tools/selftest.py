#!/usr/bin/env python3
import build_page
from checks import Check
c=Check(); c("published CAPS+ conversion paths",
             tuple(build_page.M.CAPS_PATHS[x].states for x in ("5hmC","5mC","C")),
             (("5hmC","5fC","5caC","DHU","T"),("5mC","5mC"),("C","C")))
build_page.render(build_page.M); c.report()
