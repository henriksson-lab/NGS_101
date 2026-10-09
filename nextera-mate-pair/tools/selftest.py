#!/usr/bin/env python3
import build_page
from checks import Check
M=build_page.M; c=Check(); c.section("vendor transcription")
c("duplicated mate-pair junction", M.DUPLICATE_JUNCTION,
  "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG")
build_page.render(M); c.report()
