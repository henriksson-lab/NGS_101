"""Private native DAP source route; no ampDAP or newer recipe inheritance."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

AUTHORITY='AssociatedAMNIHMS790647/PMC4907330ExperimentalProcedures148/150/156/178; exactpublishedcorrectionscope/fullSIunread.'
def protocol_graph():
    b=B('DAP-seq nativeDNA associated-manuscript coarse preparation',AUTHORITY,'SourceyoungArabidopsisleaves; nativeDNAmodificationpositionsunassigned');g=b.g
    b.op('extract','SourceextractgDNAfromyoungArabidopsisleaves',S.EXTRACTION,(b.stock('extraction_preparation','ActualsourcegDNAextractionstock, exactformulation/internalphysicalfatesunassembled',K.MOLECULAR,R.STRUCTURAL),b.cond('extraction_conditions','SourcegenomicDNAextractionconditionsunassembled')))
    b.op('fragment','SourcefragmentnativegDNA, exactfragmenter/length/enddistributionsunassembled',S.FRAGMENTATION,(b.cond('fragment_conditions','Sourcefragmentationconditionsunassembled'),))
    adapter=g.input('truncated_truseq','ActualsourcetruncatedIlluminaTruSeqadapterstock; wholeoligos/endsunassembled',K.ADAPTER)
    geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Source truncatedTruSeqadapter name positive; actualwholeoligos/truncation/occupiedDNAends/sealedproducts unassembled, not inferred frommodernTruSeq.')
    library=b.op('adapter_ligate','SourceligatetruncatedTruSeqadapter togenomicfragments',S.LIGATION,((adapter,R.STRUCTURAL),b.stock('ligation_preparation','Actualsourceadapterligationpreparedstock, enzyme/donor/formulationunassembled',K.MOLECULAR,R.STRUCTURAL)),geometry=geometry,product='Native gDNA adapterlibrary plus incomplete/unligatedclasses, native modification positions unassigned; no sourcepre-bindingPCR inselectedroute')
    b.x=g.input('tf_expression_preparation','ActualseparateHaloTF wheatgermexpressionpreparation, template/wholemolecularformulationunassembled',K.MOLECULAR)
    b.op('tf_expression','SourceHalo-taggedTFexpressed invitrowheatgermsystem',S.MAGIC,reason='Positive separateproteinexpression route; completeTF/template/translationproducts/internalmolecularstocks unassembled. Not anew sequencingDNAtemplate.')
    b.op('tf_immobilize','SourceimmobilizeexpressedHaloTFonMagneHaloTagbeads',S.BINDING,(b.stock('halo_beads','SourceMagneHaloTagbeads',K.SELECTION,R.SELECTION),))
    b.clean('tf_wash','SourcewashTFbeads beforeDNA libraryincubation')
    tfbeads=g.product('washed_tf_beads','Same physically washed HaloTF/bead selection preparation, not a fresh addition',K.SELECTION)
    g.react('tf_selection_component','Represent retained washed preparation as the actual selection reagent',S.GENERIC,((b.x,R.TEMPLATE),),(tfbeads,),changes_molecule=False);b.x=library
    bound=g.product('dna_bound','SourceHaloTFbead-retainedDNApreparation; unbound/offtarget/unligated populations andrecoveryunassigned',K.MOLECULAR);unbound=g.product('dna_unbound','ActualnonretainedDNA/workupfraction, laterfateunassigned',K.MOLECULAR)
    g.react('dna_capture','SourceincubatepreparednativeDNA librarywithwashedHaloTFbeads',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(tfbeads,R.SELECTION),b.cond('dna_binding_conditions','Actualsourcebindingconditions/quantitiesunassembled')),(bound,unbound));g.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,'Nonretainedcapturefraction notimplicitlysequenced.'));b.x=bound
    b.clean('dna_wash','SourcewashDNA-boundTFbeads')
    b.op('elute','SourceeluteDNAafterbeadwash',S.MAGIC,(b.cond('elution_conditions','Exactsourceelutionchemistry/fractionsunassembled'),),reason='Sourceelution positive; complete releasedDNA/bead/protein molecularclasses andcarryoverunassembled.')
    b.op('index_pcr','SourceamplifyrecoveredDNAwithindexedTruSeqprimers',S.MAGIC,(b.stock('indexed_pcr_preparation','ActualsourceindexedTruSeqPCRstock; completeprimers/enzymes/NTs/cyclesunassembled',K.MOLECULAR,R.STRUCTURAL),),reason='NativeDAPpost-bindingamplification positive; exactinputstock/endgeometry/wholeprimers/cycles delegated tounavailableSI. Do notborrowampDAP20cycles or2017recipe.')
    b.split('load_aliquot','Actualsequencingloadallocation, amountunassembled','SelectedDAPsequencingaliquot','UntransferredpurifiedDAPpreparation; exactfateunassembled')
    lib=g.product('library','SourceDAP-seqselectednativeDNAcapturelibrary',K.LIBRARY);g.react('load','SourceHiSeq2500100bpSRsequencing',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='AssociatedAMHiSeq2500100bpSR; fulladapter/primerlanding/indexmask unassembled. Nativeinputbeforebinding, notampDAPmethylation-erasedinput orChIProute. Unreadpublishedcorrection/finalbodyequivalence unresolved.'))
    g.add_source_item(SourceItem('unavailable_sources','ActualSupplementExperimentalProcedures PDFsmmc1/mmc5 anddeclaredcorrectionDOI10.1016/j.cell.2016.08.063 unavailable.','Correction scope notassumedchemistry-neutral; no completeaudit or absenceclaim.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('DNAadapterlibrarymadebeforebinding',g.steps['dna_capture'].inputs[0].material,'adapter_ligate_product')
    c('TFbeadfactoryactuallywashed',g.steps['tf_selection_component'].inputs[0].material,'tf_wash_kept')
    c('NoPCRerasurebeforebinding',list(g.steps).index('index_pcr')>list(g.steps).index('dna_capture'),True)
    c('ElutionuseswashedboundDNA',g.steps['elute'].inputs[0].material,'dna_wash_kept')
    c('ProteinpreparationnotprimaryDNAtemplate',g.materials['tf_expression_preparation'].primary,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.cell.2016.04.038_PMC4907330.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.cell.2016.04.038')):
        t=p.read_text();c('OriginalnativeDNAplusseparatewheatgermTF','truncated Illumina TruSeq adaptor' in t and 'in vitro wheat germ system' in t,True)
        c('OriginalsourceHiSeqSRmask','HiSeq 2500 with 100-bp SR reads' in t,True)
    c('Privatepartialremainsunready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
