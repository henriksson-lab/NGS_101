"""Associated-source H4S47C mapping; incomplete, no public page marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
@dataclass(frozen=True)
class SelectedCleavage:
    dna:str
    cuts:tuple
    inferred:bool=True
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or tuple(sorted(set(self.cuts)))!=self.cuts or any(type(i)is not int or not 0<i<len(self.dna) for i in self.cuts):raise ValueError('Explicit distinct sorted internal cuts required; not histone-position predictor')
    def fragments(self):
        ends=(0,)+self.cuts+(len(self.dna),);return tuple(self.dna[a:b] for a,b in zip(ends,ends[1:]))
def scene():
    s=Scene();s.strand('unassembled source product',[Segment('geometry','[full genomic and historical adapter geometry unassembled]',placeholder=True)]);return s

def protocol_graph_for(high_mw=False):
    if type(high_mw)is not bool:raise ValueError('Explicit observed gel-band route required')
    g=ProtocolGraph('H4S47C chemical mapping '+('higherMW followup' if high_mw else 'shortest band'),authority='AssociatedmainPMC3786739 Results30–40/Methods75–77; publishedSI MOESM262 sections1–4 lines24–112 controls preparation. Final main and linked protocol scope remain open.',representation='published_instance')
    x=g.input('yeast','Source BY4741 H4S47C strain;HHF1S47C andHHF2deletion established input, no cloning recipe inferred',K.BIOLOGICAL,primary=True)
    def op(key,label,kind,reagent,product):
        nonlocal x
        incorporated=key=='label';binding=kind==S.BINDING;enzyme=key in {'rnase','klenow','terminator'}
        c=g.input(key+'_stock',reagent,K.MOLECULAR if incorporated else K.CATALYST if enzyme else K.SELECTION if binding else K.CONDITION);y=g.product(key+'_product',product,K.MOLECULAR)
        kw={'reason':'Actual source component; internal geometry unassembled.'} if kind==S.MAGIC else {}
        inputs=[(x,R.TEMPLATE),(c,R.STRUCTURAL if incorporated else R.CATALYST if enzyme else R.SELECTION if binding else R.CONDITION)]
        if incorporated:inputs.append((g.input(key+'_buffer','Source labelingbuffer/DMSO and light-sensitive incubation',K.CONDITION),R.CONDITION))
        g.react(key,label,kind,inputs,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Source recovered preparation; full carryover/recovery unquantified',K.MOLECULAR);w=g.product(key+'_removed','Source physically separated fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Physical source separation; no perfect partition assumed.');x=y
    clean('harvest','Source harvest and1Msorbitol wash')
    op('permeabilize','Source immediate permeabilizationRT5min',S.PERMEABILIZATION,'5ml1Msorbitol/5mMBME/0.5mglyticase','Permeabilized source cells, debris and reagents; no formaldehyde crosslinking')
    clean('permeabilized_wash','Source2x5ml1Msorbitol/0.1percentNP40 wash;resuspend2ml labelingbuffer')
    op('label','Source light-sensitive cysteine labelingRT2h then4Covernight',S.CHEMICAL_CONVERSION,'N(1,10phenanthroline-5yl)iodoacetamide1.4mM/20percentDMSO','Covalently labeled histoneCys/protein-containing chromatin plus free label/remnants; no nucleotide label or fixed DNA cut assigned')
    clean('label_wash','Source2x10ml sorbitol/NP40wash;resuspend1.5mlmappingbuffer')
    op('cu','Source copper loading5min',S.BINDING,'CuCl2 final150uM','Copper-loaded labeled chromatin plus free copper/reagents')
    clean('cu_wash','Source2x10mlmappingbuffer wash;resuspend1.5ml')
    op('reduce','Source MPA treatment5min',S.CHEMICAL_CONVERSION,'Beta-mercaptopropionic acid final6mM','Reduced copper source preparation; no DNA loss assumed')
    op('cleave','Source peroxide-initiated chemical mapping20min',S.FRAGMENTATION,'Hydrogenperoxide final6mM','Hydroxyl-radical DNA fragments plus uncut DNA/proteins/reagents; all daughters retained; neither restriction motif nor fixed dyad cut predicted')
    op('quench','Source Neocuprine quench',S.QUENCHING,'Neocuprine final2.8mM','Quenched DNA/chromatin mixture')
    clean('pellet','Source spin/resuspend1ml1mMTrisHClpH7.5/0.1mMEDTA')
    op('disrupt','Source histone-DNA contact disruption1h',S.DEPROTEINIZATION,'SDS2percent/NaCl2M','Source disrupted chromatin mixture; protein and DNA remain until extraction')
    clean('extract','Source phenol-chloroform extraction and ethanol precipitation')
    op('rnase','Source RNase1mgml42Covernight',S.DIGESTION,'Source RNase identity beyond RNase unassigned','DNA plus RNA digest products/reagents; no invented intervening cleanup')
    c=g.input('gel_condition','Source2percentagarose gel and extraction',K.SELECTION);keep=g.product('gel_kept','Source200–500bp higherMW fraction' if high_mw else 'Source shortest DNA band; source nominal~150–200bp mainMethods versus~125–200bp Results',K.MOLECULAR);w=g.product('gel_rest','Other physical gel fractions',K.WASTE);g.react('gel','Source gel-band selection',S.SELECTION,((x,R.TEMPLATE),(c,R.SELECTION)),(keep,w));g.discard(w,'Source physical band choice, not computational nucleosome call.');x=keep
    op('klenow','Source FIRST Klenow large fragment blunting',S.END_REPAIR,'NEB Klenow largefragmentM0210L','Source blunted-end preparation plus residualoverhangs/uncut DNA/reagents; actual end bases unassembled')
    op('terminator','Source THEN DNATerminator end repair/phosphorylation',S.END_REPAIR,'Lucigen DNATerminator40035-2','Further sourceblunted/phosphorylated DNA with reactionremnants; no silent cleanup')
    if high_mw:
        op('capture_shortening','Source Boyle2008 biotin/MmeI shortening component',S.MAGIC,'Biotinylated sequencing adapter/MmeI site;bluntligation/pulldown/MmeI18–20bp distalcut/secondadapter source-described, linked protocol unreviewed','Source shortened end-bearing DNA library intermediate; all daughter/reagent fates and historical geometry unassembled')
    else:
        a=g.input('solid_adapter','Source ABI SOLiD adapter stock; exact oligos/end geometry unreviewed',K.ADAPTER);cat=g.input('ligase','Historical ABI ligation catalytic component unreviewed',K.CATALYST);y=g.product('adapter_product','Source adapter-ligated DNA plus unligated DNA/free adapters/reagents',K.MOLECULAR);g.react('adapter','Source ABI SOLiD adapter ligation',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(scene(),),'Source blunt/phosphorylated insert; historical full adapter geometry unassembled.'));x=y
    op('nick_translation','Source ABI Fragment Preparation nick translation component',S.MAGIC,'Historical ABI Fragment Preparation protocol; internal primer/enzyme geometry unreviewed','Source nick-translated library mixture; exact fates and residual substrates unassembled')
    op('pcr','Source ABI Fragment Preparation PCR component',S.MAGIC,'Historical ABI Fragment Preparation PCR primers/conditions unreviewed','Amplified source preparation plus originals/freeprimers/reactionremnants; no invented cleanup or PCR cyclecount')
    lib=g.product('library','Source SOLiD library preparation; exact whole geometry unassembled',K.LIBRARY);g.react('identity','Source library identity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Actual source ABI SOLiD. Shortest-band experiments fourSE/twoPE; higherMW followupSE. No Illumina adapter/read primer assumed. Cleavage profiles/nucleosome deconvolution are computational; library cut sites not predicted from dyad positions.'))
    g.add_source_item(SourceItem('pending','Final published original/SI, historical ABI Fragment Preparation protocol and source-cited Boyle2008 component','Associated-source selected chemistry only; source versions/linked predecessors and kit geometry remain open.'))
    return g

def protocol_graph():return protocol_graph_for()
def selftest():
    from checks import Check
    c=Check();a=SelectedCleavage('ACGTACGT',(2,5));c('All cleavage daughters conserve template',''.join(a.fragments()),a.dna);c('Coordinate choices inferred by default',a.inferred,True)
    for cuts in ((0,),(8,),(3,2),(2,2)):
        try:SelectedCleavage('ACGTACGT',cuts);bad=False
        except ValueError:bad=True
        c('Invalid cut coordinates rejected '+str(cuts),bad,True)
    g=protocol_graph();h=protocol_graph_for(True);c('Source extra copper washed before reduction',g.steps['reduce'].inputs[0].material,'cu_wash_kept');c('No cleanup invented after RNase',g.steps['gel'].inputs[0].material,'rnase_product');c('Repair order retained',g.steps['terminator'].inputs[0].material,'klenow_product');c('HigherMW shortening distinct from ordinary adapter route',('capture_shortening' in h.steps,'adapter' in h.steps),(True,False));c('Unknown historical geometry prevents promotion',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
