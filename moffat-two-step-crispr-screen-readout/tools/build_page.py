#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr as M
import pooled_crispr_page as P
import seqprimers as sp

OUT = HERE.parent / "moffat-two-step-crispr-screen-readout.html"


def render() -> str:
    intermediate, final = M.two_step_intermediate(), M.two_step_library()
    body = f'''<h2>PCR1 — enrich the integrated guide locus</h2>
{P.construct_panel(intermediate, "PCR1 uses ordinary locus primers and produces the ~600-bp enrichment product. No flow-cell adapters are present.")}
<h2>PCR2 — add dual-index TruSeq arms</h2>
{P.construct_panel(final, "PCR2 uses the nested hU6 and scaffold sites and yields the ~200-bp sequencing library.")}
{sp.section(final, M.two_step_sequencing_primers(), heading="Sequencing primers", intro="The protocol uses single-end sequencing with 21 dark cycles, then 26 imaged guide cycles, plus 8-cycle i7 and i5 reads.", required_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)"))}'''
    return P.page("Moffat two-step CRISPR screen readout", "01_moffat-two-step.html",
        'Source: Moffat Lab <a href="https://media.addgene.org/cms/filer_public/44/9f/449f22c6-f818-4af6-8af4-2bab1d20f00a/moffat_lab_crispr_screen_protocol_v2.pdf">Pooled CRISPR Screen Protocol v2</a> (2019), LCV2::TKOv3 branch.', body)


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
