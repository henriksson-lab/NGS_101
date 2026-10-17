"""Version declaration for 10x GEM-X Single Cell 5' V(D)J v3."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "10x-chromium-single-cell-5-vdj" / "tools"),
                str(HERE.parents[1] / "lib"), str(HERE.parents[1])]

from dataclasses import replace

from vdj_common import VDJChemistry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
from protocol_families.tenx_graph import chromium_pcr, partition_material
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Scene

ADAPTER_CAVEAT = (
    "CG000734 Rev A adds the V(D)J Read 2 adaptor as part of Ligation Mix PN-2001109 "
    "(Step 4.2 Adaptor Ligation Mix, p.74); the inspected kit tables, Step 4.2 and "
    "Appendix pp.111-113 print no adaptor oligo, only the Read 2 sequence after the "
    "dA junction in the library product (p.111). The drawn half adaptor is the CG000331 Rev E Adaptor Oligos PN-2000094 "
    "structure (Protocol Step 4.2, Appendix p.74); its use in CG000734, including the 3'-T end implied by "
    "Step 4.1 A-tailing, is inferred.")

RT_PRIMER_CAVEAT = (
    "CG000734 Rev A uses Poly-dT RT Primer B (PN-2001110, kit tables pp.8-9), not the "
    "v2 Poly-dT RT Primer PN-2000007, and its Appendix (pp.111-113) prints no RT-primer "
    "sequence. The drawn handle-AC-dT30-VN sequence is the v2 PN-2000007 sequence printed "
    "in CG000331 Rev E Appendix p.72; its use for Primer B is inferred.")


class GemXV3Chemistry(VDJChemistry):
    """GEM-X v3 provenance: which printed source each reagent sequence comes from.

    CG000734 Rev A prints the gel-bead primer and Sample Index PCR products
    (Appendix pp.111-112) and, on the image-only Additional Sequences page p.113,
    the V(D)J forward primer and every outer/inner Mix 1/2 v2 reverse primer.
    It names cDNA Primers by the same PN-2000089 as v2 (sequence printed in CG000331
    Rev E p.72) and a new Poly-dT RT Primer B whose sequence is not printed there.
    """

    def polydt_segments(self):
        return [replace(s, inferred=True) for s in super().polydt_segments()]

    def adapter_top(self):
        return [replace(s, inferred=True) for s in super().adapter_top()]

    def adapter_bottom(self):
        return [replace(s, inferred=True) for s in super().adapter_bottom()]

    def adapter_scene(self) -> Scene:
        """Predecessor half adaptor, inferred for the premixed Ligation Mix PN-2001109."""
        sc = Scene()
        sc.strand("long oligo", self.adapter_top(), label="long adaptor oligo (inferred)")
        sc.anneal("short oligo", self.adapter_bottom(), to="long oligo",
                  pair=("12-bp stem complement", "Read 2 adapter"),
                  label="short adaptor oligo (inferred)")
        sc.ligation_end("short oligo", "3' T", kind="cohesive", polarity="3-prime",
                        owner="adapter", text="inferred 3'-T ligation overhang")
        return sc


CHEM = GemXV3Chemistry(
    version="v3", umi_nt=12, guide="CG000734 Rev A",
    cdna_pcr_statement=("CG000734 Rev A uses cDNA Primers PN-2000089 (kit tables pp.8-9); "
                        "the forward and reverse sequences of that PN are printed in "
                        "CG000331 Rev E Appendix p.72, not in CG000734."),
    vdj_pcr_statement=("CG000734 Rev A Appendix p.113 (image-only page) prints the "
                       "forward primer and every {stage} reverse primer, identical to "
                       "CG000331 Rev E pp.73-74; neither guide maps each oligo to a "
                       "transcript coordinate."),
    index_pcr_statement=("CG000734 Rev A prints the V(D)J Sample Index PCR product "
                         "(Appendix p.111), not the Dual Index TT Set A primers; the "
                         "drawn primer 3' ends inside Read 1/Read 2 are derived."))


def protocol_graph() -> ProtocolGraph:
    authority=(f'Primary source: <a href="https://www.10xgenomics.com/support">'
               f'10x Genomics User Guide {CHEM.guide}</a>.')
    g=ProtocolGraph(f"10x Chromium Single Cell 5-prime V(D)J {CHEM.version}",
                    authority=authority)
    cells=g.input("cells","cells containing receptor and poly(A) transcripts",
                  MaterialKind.BIOLOGICAL,
                  state=molecule("cells",(role("V(D)J transcripts"),role("poly(A) RNA"))),
                  primary=True)
    reviewed=["cells"]
    occupied=partition_material(g,cells,reviewed,label="cells")
    rt_primer=g.input("rt_primer","Poly-dT RT Primer B (PN-2001110; sequence not printed in CG000734)",MaterialKind.OLIGO)
    rt_enzyme=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
    first=g.product("first","first-strand full-length cDNA",MaterialKind.MOLECULAR,
                    state=CHEM.cdna_template())
    g.react("reverse_transcribe","Reverse-transcribe poly(A) transcripts",
            StepKind.REVERSE_TRANSCRIPTION,
            ((occupied,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
             (rt_enzyme,InputRole.CATALYST)),(first,))
    tso=g.input("tso","barcoded gel-bead template-switch oligo",MaterialKind.OLIGO)
    switched=g.product("switched","cell-barcoded full-length cDNA",MaterialKind.MOLECULAR,
                       state=CHEM.cdna_template())
    g.react("template_switch","Template-switch onto the barcoded gel-bead oligo",
            StepKind.TEMPLATE_SWITCHING,
            ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
             (rt_enzyme,InputRole.CATALYST)),(switched,))
    reviewed.extend(("rt_primer","rt","reverse_transcribe","tso","template_switch"))
    amplified=chromium_pcr(g,"cdna_pcr",switched,CHEM.cdna_template(),reviewed,
                           label="Amplify full-length cDNA",geometry=CHEM.cdna_pcr())
    outer=chromium_pcr(g,"outer_vdj_pcr",amplified,CHEM.vdj_template(inner=True),reviewed,
                       label="Outer V(D)J PCR with Mix 1",geometry=CHEM.vdj_pcr(inner=False))
    inner=chromium_pcr(g,"inner_vdj_pcr",outer,CHEM.vdj_amplicon(),reviewed,
                       label="Inner V(D)J PCR with Mix 2",geometry=CHEM.vdj_pcr(inner=True))
    fragmentase=g.input("fragmentase","V(D)J amplicon fragmentation enzyme",
                        MaterialKind.CATALYST)
    fragments=g.product("fragments","fragmented barcode-bearing V(D)J amplicons",
                        MaterialKind.MOLECULAR)
    g.react("fragment","Fragment nested V(D)J amplicons",StepKind.FRAGMENTATION,
            ((inner,InputRole.TEMPLATE),(fragmentase,InputRole.CATALYST)),(fragments,))
    end_mix=g.input("end_mix","end-repair and dA-tailing enzymes",MaterialKind.CATALYST)
    tailed=g.product("tailed","end-repaired dA-tailed V(D)J fragments",
                     MaterialKind.MOLECULAR)
    g.react("end_prep","End-repair and dA-tail V(D)J fragments",StepKind.END_REPAIR,
            ((fragments,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(tailed,))
    adapter=g.input("read2_adapter","Read 2 adaptor in Ligation Mix PN-2001109 (not printed)",
                    MaterialKind.ADAPTER,state=CHEM.adapter_scene(),inferred=True,
                    source=authority+" (adaptor structure inferred from CG000331 Rev E PN-2000094)")
    ligase=g.input("ligase","adapter ligase",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated barcode-bearing V(D)J fragments",
                      MaterialKind.MOLECULAR,state=CHEM.pre_index_library())
    g.react("ligate","Ligate the Read-2 adapter",StepKind.LIGATION,
            ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
             (ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.known(
                 "cohesive",CHEM.adapter_scene(),statement=ADAPTER_CAVEAT))
    reviewed.extend(("fragmentase","fragment","end_mix","end_prep",
                     "read2_adapter","ligase","ligate"))
    library=chromium_pcr(g,"index_pcr",adapted,CHEM.final_library(),reviewed,
                         label="Select the barcode-bearing fragment by dual-index PCR",
                         library=True,geometry=CHEM.index_pcr())
    lengths={"Read 1":26 if CHEM.umi_nt==10 else 28,
             "Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":90}
    g.sequence(library,SequencingPlan(library=CHEM.final_library(),
        primers=CHEM.seq_primers,acquired_roles=tuple(p.role for p in CHEM.seq_primers),
        read_lengths=lengths))
    reviewed.append(library.id)
    g.declare_reviewed_scope(reviewed,authority)
    return g.require_publishable()
