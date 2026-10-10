#!/usr/bin/env python3
import build_page
from checks import Check
c=Check(); c("published inline identifier pattern", build_page.M.INLINE_LAYOUT, (3,4,2))
build_page.render(build_page.M); c.report()
