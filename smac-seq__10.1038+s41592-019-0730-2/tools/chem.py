"""Source-qualified SMAC-seq operations; private incomplete declaration."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan

@dataclass(frozen=True)
class SelectedModifiedStrand:
    """Caller-selected modifications, not an accessibility/enzyme predictor."""
    dna:str
    adenine_m6:tuple=()
    cytosine_m5:tuple=()
    inferred:bool=True
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT'):raise ValueError('Concrete DNA required')
        for positions,base in ((self.adenine_m6,'A'),(self.cytosine_m5,'C')):
            if len(set(positions))!=len(positions) or any(type(i)is not int or not 0<=i<len(self.dna) or self.dna[i]!=base for i in positions):raise ValueError('Explicit unique positions on the correct substrate base required')


def protocol_graph_for(variant='yeast_combined'):
    if variant not in {'yeast_combined','human_eco'}:raise ValueError('Source variant required')
    g=ProtocolGraph('SMAC-seq '+variant,authority='Associated manuscript PMC7968351 Methods134–168; final2020version not merged with later detailed protocol.',representation='published_instance')
    x=g.input('cells','BY4741 yeast cells' if variant=='yeast_combined' else 'GM12878 human cells',K.BIOLOGICAL,primary=True)
    def op(key,label,kind,reagent,product):
        nonlocal x
        is_enzyme=kind==S.CHEMICAL_CONVERSION
        c=g.input(key+'_conditions',reagent,K.CATALYST if is_enzyme else K.CONDITION);y=g.product(key+'_product',product,K.MOLECULAR)
        kw={'reason':'Actual source component; full molecular geometry unassembled.'} if kind==S.MAGIC else {}
        g.react(key,label,kind,((x,R.TEMPLATE),(c,R.CATALYST if is_enzyme else R.CONDITION)),(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_conditions',label,K.CONDITION);y=g.product(key+'_kept','Recovered source preparation; carryover and recovery unquantified',K.MOLECULAR);w=g.product(key+'_removed','Physical source separated fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Source separation; no ideal molecular partition claimed.');x=y
    if variant=='yeast_combined':
        clean('cell_wash','Source pellet/sorbitol wash')
        op('spheroplast','Source yeast spheroplasting30C5min',S.PERMEABILIZATION,'200ul sorbitolbuffer/10mMDTT/0.5mgml100TZymolase','Spheroplast preparation plus cell-wall debris/reagents')
        clean('spheroplast_wash','Source pellet/sorbitol wash')
    op('lysis','Source ice-cold nuclei lysis10min',S.PERMEABILIZATION,'Source NP40 nuclei lysis buffer','Nuclei and source lysis mixture')
    clean('nuclei_wash','Source nuclei pelleting/wash then resuspension in source reactionbuffer')
    if variant=='yeast_combined':
        op('gc_a','Source initial M.CviPI/EcoGII exposure30C7.5min',S.CHEMICAL_CONVERSION,'200U M.CviPI plus200U EcoGII/SAM0.6mM/sucrose300mM','Chromatin with possible source5mC/m6A addition; original DNA/native modifications and reagents retained; no complete marking asserted')
        op('replenish','Source replenishment30C7.5min',S.CHEMICAL_CONVERSION,'128pmolSAM and another100U of enzymes; division between enzymes unassigned','Further exposed chromatin plus original/replenished enzymes and reaction remnants')
        op('cpg','Source subsequent M.SssI30C7.5min',S.CHEMICAL_CONVERSION,'60U M.SssI/128pmolSAM/10mMMgCl2','Combinedm6A/5mC source preparation; existing marks retained, no base-coordinate assignment or replacement')
        op('stop','Source equal-volume stop-buffer addition',S.QUENCHING,'20mMTrisHClpH8.5/600mMNaCl/1percentSDS/10mMEDTA stock','Stopped source preparation; no physical fragment loss or conversion to reported base')
    else:
        op('a_only','Source human EcoGII37C10min',S.CHEMICAL_CONVERSION,'200UEcoGII/0.6mMSAM;1xCutSmart/0.3Msucrose','EcoGII exposed chromatin with possiblem6A additions; endogenous cytosine marks remain')
        op('stop','Source SDS stop',S.QUENCHING,'Source0.2percentSDS','Stopped preparation retaining DNA/native marks and reaction remnants')
    clean('hmw','Source MagAttract HMW DNA isolation; historical kit internals unreviewed')
    op('kit','Source SQK-LSK108 ligation-library component',S.MAGIC,'Historical ONT Ligation Sequencing Kit1D SQK-LSK108; exact original guide unreviewed','Source HMW DNA sequencing library preparation; end repair/adapters/remnants and full modification-preserving geometry unassembled')
    lib=g.product('library','Native DNA source library; full ONT adapter geometry unassembled',K.LIBRARY)
    g.react('library_identity','Recovered source library',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Source R9.4 MinION up to48h; native DNA direct readout. No PCR/bisulfite conversion imposed on the ONT route. Tombo/Albacore processing is computational, not physical loss. Modified base identities are not ordinary base substitutions.'))
    g.add_source_item(SourceItem('pending','Published final article/SI and historical MagAttract/SQK-LSK108 originals','Associated manuscript selected Methods inspected; full source/linked audit and kit geometry remain open.'))
    return g

def protocol_graph():return protocol_graph_for()
def selftest():
    from checks import Check
    c=Check();r=SelectedModifiedStrand('ACGATC',(0,3),(1,5));c('Modifications preserve DNA letters',r.dna,'ACGATC');c('Inference default retained',r.inferred,True)
    for args in [('AC',(1,),()),('AC',(),(0,)),('AC',(0,0),())]:
        try:SelectedModifiedStrand(*args);bad=False
        except ValueError:bad=True
        c('Wrong substrate or duplicate coordinate rejected '+str(args),bad,True)
    a=protocol_graph();b=protocol_graph_for('human_eco');c('Source replenishment feeds subsequent methyltransferase',a.steps['cpg'].inputs[0].material,'replenish_product');c('Human route has no imposed CpG treatment','cpg' in b.steps,False);c('Nuclei not amplified for native readout',any(s.kind==S.PCR for s in a.steps.values()),False);c('Unassembled kit prevents promotion',bool(a.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
