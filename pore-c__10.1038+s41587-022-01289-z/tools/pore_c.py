"""Molecular model for restriction-enzyme Pore-C."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      annotation_rows, complement_segments)
from restriction import RestrictionEnzyme
from endprep import LigationJunction, dA_tailed_scene, repair_and_dA_tail
from adapter_ligation import reported_adapter_end
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

def seg(name,top,tag=None,**kw): return Segment(name,top,tag,**kw)

NLAIII=RestrictionEnzyme('NlaIII','CATG',4,0)
def chromatin_locus(letter):
    return Construct([
        seg(f'locus {letter} left','X'*12,placeholder=True),
        seg(f'NlaIII site {letter}',NLAIII.site,'me'),
        seg(f'locus {letter} right','X'*12,placeholder=True),
    ],name=f'crosslinked locus {letter}')

def chromatin_rows():
    out=[]
    for letter in 'ABC':
        out.extend(rows(chromatin_locus(letter),f'proximal locus {letter}'))
    return out
def restriction_fragment_scene(letter):
    """One physical NlaIII product with its enzyme-derived 3-prime overhang."""
    body=seg(f'restriction fragment {letter}','X'*24,placeholder=True)
    over=seg(f'{letter} cohesive {NLAIII.overhang} end',NLAIII.overhang,'me')
    sc=Scene(); sc.strand('top',[body,over],label=f'fragment {letter}')
    sc.anneal('bottom',complement_segments([body]),to='top',
              pair=(body.name+"'",body.name),label=f'fragment {letter}')
    sc.mark('top',over.name,f'NlaIII {NLAIII.end} overhang')
    return sc

def digested_rows():
    out=[]
    for letter in 'ABC':
        out.extend(restriction_fragment_scene(letter).rows())
    return out
CONCATEMER=Construct([
    seg('contact fragment A','X'*24,placeholder=True),seg('A-B NlaIII ligation junction',NLAIII.site,'me'),
    seg('contact fragment B','X'*24,placeholder=True),seg('B-C NlaIII ligation junction',NLAIII.site,'me'),
    seg('contact fragment C','X'*24,placeholder=True),
],name='multiway Pore-C concatemer')
DA_CONCATEMER=repair_and_dA_tail(CONCATEMER)
TA_JUNCTION=LigationJunction('T','A')
FINAL_LIBRARY=Construct([
    seg('left motor-loaded ONT adapter','[ONT ligation adapter]','r1',placeholder=True),
    seg('left T:A junction',TA_JUNCTION.adapter_overhang,
        bottom=TA_JUNCTION.insert_overhang),*CONCATEMER.segments,
    seg('right A:T junction',TA_JUNCTION.insert_overhang,
        bottom=TA_JUNCTION.adapter_overhang),
    seg('right ONT adapter','[ONT ligation adapter]','r2',placeholder=True),
],name='Pore-C nanopore library')

def rows(con,label=''):
    unknown=tuple(s.name+"'" for s in con if s.is_role_token())
    return [*Scene.duplex(list(con),label=label,unpaired=unknown).rows(),*annotation_rows(con)]

def final_rows():
    unknown=tuple(s.name+"'" for s in FINAL_LIBRARY if s.is_role_token())
    sc=Scene.duplex(list(FINAL_LIBRARY),label='Pore-C',unpaired=unknown)
    sc.junction('top','left T:A junction','contact fragment A','adapter ligation')
    sc.junction('top','contact fragment C','right A:T junction','adapter ligation')
    return [*sc.rows(),*annotation_rows(FINAL_LIBRARY)]

INITIAL_NAME='Formaldehyde-crosslinked proximal chromatin'
INITIAL_ROWS=tuple(chromatin_rows())

def workflow():
    wf=Workflow(MolecularState(INITIAL_NAME,INITIAL_ROWS))
    wf.react('Overnight NlaIII restriction digest',digested_rows(),
             note='NlaIII cuts CATG and exposes mutually compatible cohesive ends while spatial contacts remain crosslinked.')
    ligrows=rows(CONCATEMER,'multiway contact')
    # Exact covalent boundaries are owned by Scene junction markers.
    sc=Scene.duplex(list(CONCATEMER),label='multiway contact')
    sc.junction('top','A-B NlaIII ligation junction','contact fragment B','proximity ligation')
    sc.junction('top','B-C NlaIII ligation junction','contact fragment C','proximity ligation')
    ligrows=[*sc.rows(),*annotation_rows(CONCATEMER)]
    wf.react('T4 DNA ligase in intact crosslinked nuclei',ligrows,
             note='Compatible ends from several proximal loci become one chimeric dsDNA polymer.')
    wf.react('Proteinase K reverse-crosslinking and DNA purification',rows(CONCATEMER,'released concatemer'),
             note='The covalent multiway DNA concatemer is released intact; no PCR is performed.')
    wf.react('Optional size selection, then repair and dA-tailing',dA_tailed_scene(DA_CONCATEMER).rows(),
             note='The later ONT info sheet recommends 0.85X SPRI selection when yield permits and a native ligation kit; the original paper\'s kit remains unverified.')
    wf.react('ONT ligation-adapter attachment',final_rows(),
             note='In the later V14 example, Salt-T4 DNA ligase attaches motor-loaded adapters to repaired concatemer ends.')
    return wf

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/s41587-022-01289-z"
    vendor="Later ONT RE-Pore-C info sheet and SQK-LSK114 ligation guidance; the kit used in the defining 2022 paper is unverified"
    graph=ProtocolGraph("Pore-C",authority=source)
    chromatin=graph.input("chromatin","formaldehyde-crosslinked proximal chromatin",
                          MaterialKind.BIOLOGICAL,state=tuple(chromatin_rows()),primary=True)
    nla=graph.input("nlaiii","NlaIII restriction enzyme",MaterialKind.CATALYST)
    fragments=graph.product("restriction_fragments","crosslinked NlaIII restriction fragments",
                            MaterialKind.MOLECULAR,state=tuple(digested_rows()))
    graph.react("digest","Digest crosslinked chromatin with NlaIII",StepKind.DIGESTION,
                ((chromatin,InputRole.TEMPLATE),(nla,InputRole.CATALYST)),(fragments,))
    ligase=graph.input("t4_ligase","T4 DNA ligase",MaterialKind.CATALYST)
    concatemer=graph.product("crosslinked_concatemer","crosslinked multiway proximity concatemer",
                             MaterialKind.MOLECULAR,state=CONCATEMER)
    graph.react("proximity_ligation","Ligate proximal restriction fragments in intact nuclei",
                StepKind.PROXIMITY_LIGATION,
                ((fragments,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(concatemer,))
    proteinase=graph.input("proteinase_k","proteinase K and reverse-crosslink conditions",
                           MaterialKind.CATALYST)
    released=graph.product("released_concatemer","purified intact multiway DNA concatemer",
                           MaterialKind.MOLECULAR,state=CONCATEMER)
    graph.react("deproteinize","Reverse crosslinks and remove proteins",StepKind.DEPROTEINIZATION,
                ((concatemer,InputRole.TEMPLATE),(proteinase,InputRole.CATALYST)),(released,))
    selector=graph.input("long_fragment_selection","optional 0.85X SPRI size-selection beads",
                         MaterialKind.SELECTION)
    selected=graph.product("selected_concatemer","long selected multiway concatemer",
                           MaterialKind.MOLECULAR,state=CONCATEMER)
    rejected=graph.product("short_fragments","rejected short DNA fragments",MaterialKind.WASTE)
    graph.react("size_selection","Enrich long multi-contact concatemers when input permits",StepKind.SELECTION,
                ((released,InputRole.TEMPLATE),(selector,InputRole.SELECTION)),(selected,rejected))
    graph.discard(rejected,"short fragments are removed to preserve multi-contact molecules")
    endprep=graph.input("end_repair_mix","end-repair and dA-tailing enzymes",MaterialKind.CATALYST)
    tailed=graph.product("a_tailed_concatemer","end-repaired dA-tailed concatemer",
                         MaterialKind.MOLECULAR,state=dA_tailed_scene(DA_CONCATEMER))
    graph.react("end_prepare","End-repair and dA-tail the selected concatemer",StepKind.END_REPAIR,
                ((selected,InputRole.TEMPLATE),(endprep,InputRole.CATALYST)),(tailed,))
    adapter=graph.input("ont_adapter","motor-loaded ONT ligation adapter",MaterialKind.ADAPTER,
                        inferred=True,source=vendor)
    adapter_ligase=graph.input("ont_salt_t4_ligase","Salt-T4 DNA ligase for later ONT V14 kit",
                               MaterialKind.CATALYST,inferred=True,source=vendor)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,inferred=True,
                          source=vendor)
    geometry=reported_adapter_end(
        "ONT ligation adapter",kind="cohesive",polarity="3-prime",terminal_sequence="T",
        substrate_scenes=(dA_tailed_scene(DA_CONCATEMER),),
        statement="The later ONT info sheet recommends an SQK-LSK114 native ligation branch; its kit documentation specifies dA-tailed DNA and complementary adapter dT. This does not establish the kit used in the 2022 paper.")
    graph.react("adapter_ligation","Attach motor-loaded ONT ligation adapters",
                StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (adapter_ligase,InputRole.CATALYST)),(library,),geometry=geometry)
    graph.sequence(library,SequencingPlan(
        limitation="Nanopore sequencing is primer-free. The exact kit and adapter used by the defining 2022 paper remain unresolved; the drawn dA/dT attachment is a later vendor-documented example.",
        acquired_roles=("Sequencing output",),read_uses={"Sequencing output":"include"}))
    graph.declare_reviewed_scope((
        "chromatin","nlaiii","digest","t4_ligase","proximity_ligation","proteinase_k",
        "deproteinize","long_fragment_selection","size_selection","end_repair_mix",
        "end_prepare","ont_adapter","ont_salt_t4_ligase","adapter_ligation","library"),vendor)
    return graph.require_publishable()
