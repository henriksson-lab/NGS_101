"""Version declaration for 10x Next GEM Single Cell 5' V(D)J v2."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "10x-chromium-single-cell-5-vdj" / "tools"),
                str(HERE.parents[1] / "lib"), str(HERE.parents[1])]

from vdj_common import CDNA_FORWARD, CDNA_REVERSE, VDJ_FORWARD, VDJChemistry
from chemdraw import Construct, Segment
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
from adapter_ligation import AdapterLigationGeometry

CHEM = VDJChemistry(
    version="v2", umi_nt=10, guide="CG000331 Rev F",
    cdna_pcr_statement=("CG000331 Rev F Appendix, Protocol Step 2.2 prints both cDNA Primers "
                        "(Forward Partial Read 1 and Reverse Non-poly(dT))."),
    vdj_pcr_statement=("The CG000331 Rev F Appendix (Protocol Steps 3.1 and 3.3) prints every {stage} "
                       "reverse-primer sequence, and its V(D)J Amplification overview states "
                       "that the primers are specific to TCR or BCR constant regions; it does "
                       "not map each oligo to a transcript coordinate."),
    index_pcr_statement=("CG000331 Rev F Appendix, Protocol Step 4.4 prints both Dual Index "
                         "TT Set A primer structures."))
# CG000331 Rev F, Appendix p. 74, Human T Cell Mix 1/2 v2 (PN-2000242/2000246).
HUMAN_T_OUTER = ("TGAAGGCGTTTGCACATGCA", "TCAGGCAGTATCTGGAGTCATTGAG")
HUMAN_T_INNER = ("AGTCTCTCAGCTGGTACACG", "TCTGATGGCTCAAACACAGC")


APPENDIX = "CG000331 Rev F Appendix 'Oligonucleotide Sequences', pp. 72–75"


def _oligo(name: str, segments) -> Construct:
    return Construct(list(segments), name=name)


def ligation_ready_adapter():
    """Adaptor Oligos PN-2000094 as printed (Protocol Step 4.2), with its 3'-T end declared."""
    sc = CHEM.adapter_scene()
    sc.ligation_end("short oligo", "3' T", kind="cohesive", polarity="3-prime", owner="adapter")
    return sc


def _pcr(g, prefix, template, product, reviewed, *, label, forward, reverse, polymerase,
         geometry, kind=MaterialKind.MOLECULAR):
    """One local PCR with source-named reagents (replaces the legacy shared recipe)."""
    f = g.input(prefix + "_forward", forward[0], MaterialKind.OLIGO, state=forward[1],
                source=APPENDIX)
    r = g.input(prefix + "_reverse", reverse[0], MaterialKind.OLIGO, state=reverse[1],
                source=APPENDIX)
    pol = g.input(prefix + "_polymerase", polymerase, MaterialKind.CATALYST)
    out = g.product(prefix + "_product", label + " product", kind, state=product)
    g.react(prefix, label, StepKind.PCR,
            ((template, InputRole.TEMPLATE), (f, InputRole.STRUCTURAL),
             (r, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (out,), geometry=geometry)
    reviewed.extend((f.id, r.id, pol.id, prefix))
    return out


def protocol_graph() -> ProtocolGraph:
    authority=(f'Primary source: <a href="https://www.10xgenomics.com/support">'
               f'10x Genomics User Guide {CHEM.guide}</a>.')
    g=ProtocolGraph(f"10x Chromium Single Cell 5-prime V(D)J {CHEM.version}",
                    authority=authority)
    cells=g.input("cells","cells containing receptor and poly(A) transcripts",
                  MaterialKind.BIOLOGICAL,
                  state=molecule("cells",(role("V(D)J transcripts"),role("poly(A) RNA"))),
                  primary=True)
    gems=g.input("gems","Chromium Next GEM partitioning with Single Cell 5' Gel Beads v2",
                 MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied GEMs containing cells",MaterialKind.MOLECULAR,
                       state=cells.state)
    failed=g.product("failed","empty, multiplet or unresolved GEM products",MaterialKind.WASTE)
    g.react("partition","Partition cells with barcoded gel beads",StepKind.PARTITION,
            ((cells,InputRole.TEMPLATE),(gems,InputRole.CONDITION)),(occupied,failed),
            changes_molecule=False)
    g.discard(failed,"does not retain a resolvable single-cell barcode")
    reviewed=["cells","gems","partition"]
    rt_primer=g.input("rt_primer","Poly-dT RT Primer PN-2000007",MaterialKind.OLIGO,
                      state=_oligo("Poly-dT RT Primer",CHEM.polydt_segments()),
                      source=APPENDIX+", Protocol Step 1.5")
    rt_enzyme=g.input("rt","RT Enzyme C (template-switching reverse transcriptase)",
                      MaterialKind.CATALYST)
    first=g.product("first","first-strand full-length cDNA",MaterialKind.MOLECULAR,
                    state=CHEM.cdna_template())
    g.react("reverse_transcribe","Reverse-transcribe poly(A) transcripts",
            StepKind.REVERSE_TRANSCRIPTION,
            ((occupied,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
             (rt_enzyme,InputRole.CATALYST)),(first,))
    tso=g.input("tso","Gel Bead Primer (barcoded template-switch oligo)",MaterialKind.OLIGO,
                state=_oligo("Gel Bead Primer",CHEM.gel_bead_segments()),
                source=APPENDIX+", Protocol Step 1.5")
    switched=g.product("switched","cell-barcoded full-length cDNA (GEM-RT products)",
                       MaterialKind.MOLECULAR,state=CHEM.cdna_template())
    g.react("template_switch","Template-switch onto the barcoded gel-bead oligo",
            StepKind.TEMPLATE_SWITCHING,
            ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
             (rt_enzyme,InputRole.CATALYST)),(switched,))
    recovery=g.input("gem_recovery","Recovery Agent GEM break and Dynabeads MyOne SILANE cleanup (Step 2.1)",
                     MaterialKind.CONDITION)
    recovered=g.product("recovered","recovered first-strand cDNA",MaterialKind.MOLECULAR,
                        state=switched.state)
    g.react("recover","Break GEMs and clean up first-strand cDNA",StepKind.EXTRACTION,
            ((switched,InputRole.TEMPLATE),(recovery,InputRole.CONDITION)),(recovered,),
            changes_molecule=False)
    reviewed.extend(("rt_primer","rt","reverse_transcribe","tso","template_switch",
                     "gem_recovery","recover"))
    seg1=lambda name,seq,tag=None: [Segment(name,seq,tag)]
    amplified=_pcr(g,"cdna_pcr",recovered,CHEM.cdna_template(),reviewed,
                   label="Amplify full-length cDNA (Step 2.2)",
                   forward=("cDNA Primers PN-2000089 forward (Partial Read 1)",_oligo("cDNA forward",seg1("Partial Read 1",CDNA_FORWARD,"r1"))),
                   reverse=("cDNA Primers PN-2000089 reverse (Non-poly(dT))",_oligo("cDNA reverse",seg1("Non-poly(dT)",CDNA_REVERSE,"tso"))),
                   polymerase="Amp Mix",geometry=CHEM.cdna_pcr())
    outer=_pcr(g,"outer_vdj_pcr",amplified,CHEM.vdj_template(inner=True),reviewed,
               label="Outer V(D)J PCR with Mix 1 (Step 3.1)",
               forward=("T/B Cell Mix 1 v2 forward PCR Primer",_oligo("V(D)J forward",seg1("V(D)J forward",VDJ_FORWARD,"r1"))),
               reverse=("T/B Cell Mix 1 v2 reverse outer constant-region pool (human T PN-2000242: 2 oligos)",None),
               polymerase="Amp Mix",geometry=CHEM.vdj_pcr(inner=False))
    inner=_pcr(g,"inner_vdj_pcr",outer,CHEM.vdj_amplicon(),reviewed,
               label="Inner V(D)J PCR with Mix 2 (Step 3.3)",
               forward=("T/B Cell Mix 2 v2 forward PCR Primer",_oligo("V(D)J forward",seg1("V(D)J forward",VDJ_FORWARD,"r1"))),
               reverse=("T/B Cell Mix 2 v2 reverse inner constant-region pool (human T PN-2000246: 2 oligos)",None),
               polymerase="Amp Mix",geometry=CHEM.vdj_pcr(inner=True))
    fragmentase=g.input("fragmentase","Fragmentation Enzyme",MaterialKind.CATALYST)
    fragments=g.product("fragments","fragmented barcode-bearing V(D)J amplicons",
                        MaterialKind.MOLECULAR)
    g.react("fragment","Fragment nested V(D)J amplicons",StepKind.FRAGMENTATION,
            ((inner,InputRole.TEMPLATE),(fragmentase,InputRole.CATALYST)),(fragments,))
    end_mix=g.input("end_mix","Fragmentation Enzyme/Buffer end repair and A-tailing",MaterialKind.CATALYST)
    tailed=g.product("tailed","end-repaired dA-tailed V(D)J fragments",
                     MaterialKind.MOLECULAR)
    g.react("end_prep","End-repair and dA-tail V(D)J fragments",StepKind.END_REPAIR,
            ((fragments,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(tailed,))
    adapter=g.input("read2_adapter","Adaptor Oligos PN-2000094",MaterialKind.ADAPTER,
                    state=CHEM.adapter_scene(),source=APPENDIX+", Protocol Step 4.2")
    ligase=g.input("ligase","DNA Ligase",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated barcode-bearing V(D)J fragments",
                      MaterialKind.MOLECULAR,state=CHEM.pre_index_library())
    g.react("ligate","Ligate the Read-2 adaptor",StepKind.LIGATION,
            ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
             (ligase,InputRole.CATALYST)),(adapted,),
            geometry=AdapterLigationGeometry.known("cohesive",ligation_ready_adapter(),
                statement=APPENDIX+", Protocol Step 4.2 prints the 33-nt top strand paired over 12 nt by 3'-TCTAGCCTTCTCG-5', leaving a single 3'-T; 5' phosphorylation is not stated."))
    reviewed.extend(("fragmentase","fragment","end_mix","end_prep",
                     "read2_adapter","ligase","ligate"))
    library=_pcr(g,"index_pcr",adapted,CHEM.final_library(),reviewed,
                 label="Select the barcode-bearing fragment by dual-index PCR (Step 4.4)",
                 forward=("Dual Index Kit TT Set A forward (P5–i5–Partial Read 1)",_oligo("TT forward",CHEM.index_p5_segments())),
                 reverse=("Dual Index Kit TT Set A reverse (P7–i7–Partial Read 2)",_oligo("TT reverse",CHEM.index_p7_segments())),
                 polymerase="Amp Mix",geometry=CHEM.index_pcr(),kind=MaterialKind.LIBRARY)
    lengths={"Read 1":26 if CHEM.umi_nt==10 else 28,
             "Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":90}
    g.sequence(library,SequencingPlan(library=CHEM.final_library(),
        primers=CHEM.seq_primers,acquired_roles=tuple(p.role for p in CHEM.seq_primers),
        read_lengths=lengths))
    reviewed.append(library.id)
    g.declare_reviewed_scope(reviewed,authority)
    return g.require_publishable()
