#!/usr/bin/env python3
import build_page
from checks import Check
c=Check(); c("published UMI and linker length", (build_page.M.UMI_LENGTH,build_page.M.SMART_LINKER_LENGTH), (12,7))
build_page.render(build_page.M); c.report()
