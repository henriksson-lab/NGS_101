"""Version declaration for 10x Next GEM Single Cell 5' V(D)J v2."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "10x-chromium-single-cell-5-vdj" / "tools"),
                str(HERE.parents[1] / "lib"), str(HERE.parents[1])]

from vdj_common import VDJChemistry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
from protocol_families.tenx_graph import READ2_ADAPTER_END, chromium_pcr, partition_material

CHEM = VDJChemistry(version="v2", umi_nt=10, guide="CG000331 Rev F")


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
    rt_primer=g.input("rt_primer","poly(dT) RT primer",MaterialKind.OLIGO)
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
    adapter=g.input("read2_adapter","3-prime-dT Read-2 adapter",MaterialKind.ADAPTER)
    ligase=g.input("ligase","adapter ligase",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated barcode-bearing V(D)J fragments",
                      MaterialKind.MOLECULAR,state=CHEM.pre_index_library())
    g.react("ligate","Ligate the Read-2 adapter",StepKind.LIGATION,
            ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
             (ligase,InputRole.CATALYST)),(adapted,),geometry=READ2_ADAPTER_END)
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
