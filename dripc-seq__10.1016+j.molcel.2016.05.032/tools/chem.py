from hybrid_base_protocols import get
SPEC = get("dripc-seq")
SPEC.CAVEAT = ("The authors’ detailed protocol prints the P5/P7 library PCR primers and "
               "example adapter-trimming motifs. It names TruSeq Single index ligation "
               "adapters but does not print the selected index oligos or their end geometry.")
SPEC.oligos = lambda: (
    "<b>Library PCR P5</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA-3′",
    "<b>Library PCR P7</b> — 5′-CAAGCAGAAGACGGCATACGAGAT-3′",
)
