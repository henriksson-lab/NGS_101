#!/usr/bin/env python3
import build_page
from checks import Check
c=Check(); c("published run recipe", build_page.M.READ_LENGTHS, {"Read 1":75,"Index 1 (i7)":6,"Index 2 (i5)":10})
build_page.render(build_page.M); c.report()
