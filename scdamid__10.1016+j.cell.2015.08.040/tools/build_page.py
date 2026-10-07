#!/usr/bin/env python3
"""Build the molecule-focused scDamID schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import scdamid as D
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel
from page import head, info

OUT = HERE.parent / "scdamid.html"


def oligos() -> str:
    preamp = [D.secondary("N4", "N" * D.RANDOM_NT, placeholder=True),
              D.secondary("DamID stub", D.PREAMP_STUB),
              D.secondary("restored junction", "TC")]
    return ("<h2>Key oligos</h2><seq>"
            + oligo("DamID AdRt", [D.secondary("", D.ADRT)])
            + oligo("DamID AdRb", [D.secondary("", D.ADRB)])
            + oligo("single pre-amplification primer", preamp)
            + "</seq>"
            + panel(D.damid_adaptor_scene().rows(), cls="small",
                    caption="INFERRED — adaptor fold reconstructed from the secondary-source "
                            "oligos. AdRb pairs only the 10-nt ligating end of AdRt."))


def construction() -> str:
    final = D.final_library()
    dpni = Scene.duplex(D.dpni_fragment().segments)
    final_sc = Scene.duplex(final.segments)
    final_sc.strands["top"].label = final_sc.strands["bottom"].label = ""
    return f"""<h2>Library construction</h2>
<h3>(1) DpnI cuts Dam-methylated GATC sites</h3>
{panel(dpni.rows(), cls="small",
       caption="The selected blunt fragment lies between two Dam-methylated GATC sites.")}
<h3>(2) Ligate the DamID adaptor to both blunt ends</h3>
{panel(D.damid_adaptor_scene().rows(), cls="small",
       caption="INFERRED — adaptor sequence and strand arrangement come from the secondary "
               "reconstruction; the defining paper establishes adaptor ligation after DpnI.")}
<h3>(3) Amplify with one N4-tailed primer</h3>
{panel(D.preamp_duplex().rows(), cls="small",
       caption="INFERRED — the same restored adaptor/GATC junction at both ends lets one primer "
               "amplify the selected fragment.")}
<h3>(4) End repair, dA-tail, ligate the indexed Illumina adaptor and PCR</h3>
{panel([*final_sc.rows(), *annotation_rows(final)], cls="long",
       caption="INFERRED — final library reconstructed from secondary-source oligos. Dotted "
               "outer regions are not printed in the accessible defining paper.")}"""


def sequencing() -> str:
    return ("<h2>Sequencing</h2>"
            + info("Single-end 51-nt Read 1 plus a 6-nt i7 index. Read 1 begins with N4, "
                   "the DamID stub and the restored GATC before entering genomic DNA.")
            + sp.section(D.final_library(), D.SEQ_PRIMERS,
                         intro="Only the two reads actually used by this single-end, "
                               "single-index run are declared.",
                         required_roles=D.RUN_ROLES))


def render() -> str:
    return "\n".join([
        head("scDamID library chemistry"), '<div class="wrap">',
        '<h1>scDamID &mdash; single-cell nuclear-lamina contact mapping</h1>',
        info('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.08.040">'
             'Kind et al., <i>Cell</i> (2015)</a>.'),
        '<div class="caveat"><b>Sequence source.</b> The accessible paper establishes DpnI '
        'selection, adaptor ligation, single-primer amplification and indexed Illumina '
        'sequencing, but its oligo tables are unavailable. All dotted bases and steps labelled '
        '<b>INFERRED</b> use the upstream scg_lib_structs reconstruction.</div>',
        oligos(), construction(), sequencing(), '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
