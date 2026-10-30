"""Private neuronal Repair-seq: only associated source Figure1 preparation."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
class B:
    def __init__(self,title,authority,start):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,K.BIOLOGICAL,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture with incomplete/unreacted classes; occupancy and recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; whole source-specific molecular inputs/products unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual retained source fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Actual separation only, not perfect specificity or complete molecular cleanup.');self.x=y

AUTHORITY='Associated AM NIHMS1747554/PMC9179101 main89/Fig1caption270/MatMethdelegation193; no final-equivalence claim.'
def protocol_graph():
    b=B('Repair-seq neuronal synthesis, associated-manuscript coarse preparation',AUTHORITY,'Source prepared ESC-induced postmitotic neurons; differentiation chemistry outside graph');g=b.g
    b.op('incorporate','Source culture neurons with EdU for24h',S.MAGIC,
       (b.stock('edu','Actual5-ethynyl-2prime-deoxyuridine supplied to cellular culture; dose unassembled',K.MOLECULAR,R.STRUCTURAL),b.cond('culture_conditions','Source24h culture; medium/formulation unassembled')),
       reason='Positive cellular nonreplicative incorporation; uptake/phosphorylation/metabolic intermediates, repair pathways and molecular occupancies unassembled. Not a fresh in-vitro EdUTP addition.',product='Source cellular genomic DNA contains EdU-incorporated repair classes; cellular metabolites/unincorporated precursor retained until actual extraction')
    b.op('extract','Source isolate genomes',S.EXTRACTION,(b.cond('extract_conditions','Actual source genome extraction; internal reagents/fractions delegated to gatedMatMeth'),),product='Source extracted genomic DNA mixture; exact extraction molecular classes unassembled')
    b.op('sonicate','Source fragment isolated genomes by sonication',S.FRAGMENTATION,(b.cond('sonication_conditions','Actual source sonication; fragment sizes/end chemistry unassembled'),))
    b.op('click','Source click reaction adds biotin to EdU',S.CHEMICAL_CONVERSION,
       (b.stock('biotin_reagent','Actual source biotin-bearing click reactant; exact linker/formulation unassembled',K.MOLECULAR,R.STRUCTURAL),b.cond('click_conditions','Source click catalyst/conditions delegated to gatedMatMeth')),
       product='Source biotin-labelled EdU DNA, unlabelled DNA, incomplete reaction and free-label classes; no quantitative specificity/yield claim')
    bound=g.product('bound','Source streptavidin-bead retained DNA preparation; recovery/specificity unassigned',K.MOLECULAR)
    unbound=g.product('unbound','Source physically non-retained capture fraction; exact composition/later fate unassigned',K.MOLECULAR)
    g.react('capture','Source enrich biotin-DNA fragments on streptavidin beads',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('streptavidin','Source streptavidin beads',K.SELECTION,R.SELECTION)),(bound,unbound));g.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,'Non-retained fraction is not implicitly sequenced.'));b.x=bound
    b.op('amplify','Source subsequently amplify retained material; whole library preparation delegated to actualMatMeth',S.MAGIC,
       (b.stock('amplification_preparation','Actual prepared source amplification/library stock; whole primers, nucleotide pool, catalysts and quantities unassembled',K.MOLECULAR,R.STRUCTURAL),),
       reason='Associated source Fig1 specifies amplification after capture; exact adapter installation, primer ends and PCR programme remain in gated source, not filled from another library recipe.')
    lib=g.product('library','Source Repair-seq NGS library preparation; full DNA/read geometry unassembled',K.LIBRARY)
    g.react('load','Source NGS of enriched/amplified preparation; loading aliquot/stock remainder unassembled',S.MAGIC,((b.x,R.TEMPLATE),),(lib,),reason='Positive NGS readout; exact platform/run mask/loading allocation and complete molecules delegated to gated MatMeth.')
    g.sequence(lib,SequencingPlan(limitation='Associated manuscript source only. Exact platform, adapters/ends, primers and runmask unassembled; no final-equivalence, completeness or source-absence claim. Other RNA/ATAC/ChIP/HiC comparators are not production steps.'))
    g.add_source_item(SourceItem('matmeth','Actual named associated MatMeth/Figs supplement is gated/TDM.','NIHMS1747554-supplement-Figs__S1-S20___Mat___Meth__Refs_.pdf; source partial only.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Source click follows physical sonication',g.steps['click'].inputs[0].material,'sonicate_product')
    c('Biotin reactant contributes matter',next(x.role for x in g.steps['click'].inputs if x.material=='biotin_reagent'),R.STRUCTURAL)
    c('Amplification uses physically bound fraction',g.steps['amplify'].inputs[0].material,'bound')
    c('Non-retained material not silently sequenced',g.terminals['unbound'].disposition,Disposition.UNRESOLVED)
    c('Cellular precursor not fresh EdUTP',g.materials['edu'].label.find('deoxyuridine')>=0,True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'science.abb9032_PMC9179101.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1126/science.abb9032')):
        t=p.read_text();c('Original positive Fig1 process scope','fragmented with sonication, a click reaction adds biotin to the EdU' in t,True)
        c('Actual whole recipe delegated','Materials and methods are available as supplementary materials.' in t,True)
    c('Partial model cannot claim public readiness',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
