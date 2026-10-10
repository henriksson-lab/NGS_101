"""Molecular model for restriction-enzyme Pore-C."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      annotation_rows, complement_segments)
from restriction import RestrictionEnzyme
from endprep import LigationJunction, dA_tailed_scene, repair_and_dA_tail

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
    wf.react('Long-fragment selection, end repair and dA-tailing',dA_tailed_scene(DA_CONCATEMER).rows(),
             note='The published workflow preserves long multi-contact polymers for nanopore sequencing.')
    wf.react('ONT ligation-adapter attachment',final_rows(),
             note='Motor-loaded long-read adapters are ligated to the repaired concatemer ends.')
    return wf
