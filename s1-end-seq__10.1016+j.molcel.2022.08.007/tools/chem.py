"""Private S1-END-seq associated-AM frontend and source outline; not complete END recipe."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
class B:
    def __init__(self,title,authority,start):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,K.MOLECULAR,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source copresent reaction mixture; incomplete/unreacted classes and carryover unassigned'):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']='Positive source operation; historical component formulation/whole molecular products unassembled, not source absence.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def split(self,key,label,branches):
        yy=tuple(self.g.product(key+'_'+name,txt,K.MOLECULAR)for name,txt,disp in branches);self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),yy)
        for y,(_,txt,disp)in zip(yy,branches):
            if disp is not None:self.g.mark_terminal(Terminal(y.id,disp,txt))
        self.x=yy[0];return yy
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Source retained physical fraction; recovery/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated wash/workup fraction',K.WASTE);self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Source separated workup; no perfect specificity/yield.');self.x=y
    def capture(self,key,label,inputs):
        y=self.g.product(key+'_bound','Source bead/membrane bound mixture; purity unassigned',K.MOLECULAR);u=self.g.product(key+'_unbound','Unretained physical capture fraction',K.MOLECULAR);self.g.react(key,label,S.AFFINITY_CAPTURE,((self.x,R.TEMPLATE),)+tuple(inputs),(y,u));self.g.mark_terminal(Terminal(u.id,Disposition.UNRESOLVED,'No later recovery assigned.'));self.x=y
    def finish(self,limitation):
        self.split('loading_split','Source library aliquot for loading/pooling; amounts unassembled',[('aliquot','Selected library loading aliquot',None),('remainder','Untransferred library stock',Disposition.UNRESOLVED)])
        lib=self.g.product('library','Selected source productionlibrary aliquot; other loading contributions unassigned',K.LIBRARY);self.g.react('loading_boundary','Represent source sequencing aliquot without inferred pool recipe',S.GENERIC,((self.x,R.TEMPLATE),),(lib,),changes_molecule=False);self.g.sequence(lib,SequencingPlan(limitation=limitation));self.g.add_source_item(SourceItem('pending','Whole molecular products/history-specific kit internals/read landing unassembled.','Selected bounded pass only; no confirmed absence/full chemistry proof.'));return self.g

from illumina import TRUSEQ_READ1,NEBNEXT_ARM_READ2
AD1_PREFIX=revcomp(TRUSEQ_READ1)[1:-1]
AD1_LOOP='UUTUTUU'
AD1=AD1_PREFIX+AD1_LOOP+TRUSEQ_READ1
AD2='GATCGGAAGAGCACACGTC'+'U'*8+'AGACGTGTGCTCTTCCGATCT'
RAW1='5P-'+AD1_PREFIX+'UU[Biotin-dT]U[Biotin-dT]UU'+TRUSEQ_READ1[:-1]+'*'+TRUSEQ_READ1[-1]
RAW2='5P-'+AD2[:-1]+'*'+AD2[-1]
SOURCE_OLIGO_EXCEPTIONS={'adapter1':'AssociatedAM originalKeyResources821 whole modified ordering read directly inJATS: 5P GATCGGAAGAGCGTCGTGTAGGGAAAGAGTG UU[Biotin-dT]U[Biotin-dT]UU ACACTCTTTCCCTACACGACGCTCTTCCGATC* T. Text matcher fragments Biotin-dT notation. Local normalized-original whole modified pattern plus actual model polymer fixture retain every fixed base/dU/biotinT/PS, no U-to-N/T substitution.'}
def inventory(n):
    s=Scene()
    if n==1:
        segs=[Segment('prefix',AD1_PREFIX),Segment('dU2','UU',note='Source deoxyuridines'),Segment('biotin_T1','T',note='SourcebiotindT'),Segment('dU1','U',note='Sourcedeoxyuridine'),Segment('biotin_T2','T',note='SourcebiotindT'),Segment('dU2b','UU',note='Sourcedeoxyuridines'),Segment('suffix',TRUSEQ_READ1)]
    else:segs=[Segment('prefix',AD2[:19]),Segment('dU8','U'*8,note='Source8deoxyuridines'),Segment('suffix',AD2[27:])]
    s.strand('ordering',segs,mod5='P',mod3='OH; sourcephosphorothioatefinalinternucleotidebond');return s
def unknown_end(n):
    s=Scene();s.strand('end',[Segment('unassembled','[source ligation end]',placeholder=True,inferred=True)],label='Unknown localhistorical END endgeometry, notpublicationabsence');return AdapterLigationGeometry('unavailable','unavailable',(inventory(n),s),'Knownwholeadapterordering5P/dU/biotinT/PS drawnasfold-unassembledinventory. ActualsourceENDligationoutlineknown; source-linkedhistoricalendpreparation/wholeligationgeometryunassembledinthispass.')
def protocol_graph():
    b=B('S1-END-seq source-specific200Ufrontend','AssociatedAMPMC9547894main203/209/Fig1caption636/KeyResources821,823; oneoriginal. Finalequivalenceunverified. LaterENDdetailedrecipeunassembled.','Source10–20millioncells; sampleallocation/applicationunassigned');g=b.g
    b.op('embed','Source1percentagaroseplugs',S.MAGIC,(b.stock('agarose','Sourceagarosephysicalmatrix',K.MOLECULAR,R.STRUCTURAL),),product='Sourceagaroseembeddedcells; exactplugnumber/allocationunassembled')
    b.op('proteinase','SourcePK1h50Cthen7h37C',S.DEPROTEINIZATION,(b.stock('pk','SourceProteinaseK'),b.cond('pk_conditions','Source50C/37C programme, otherhistoricalstockdetailsunassembled')))
    b.clean('pk_wash','Source1xWB/3xTE,15mineachretainplugs')
    b.op('rna_digest','SourceRNaseA1h37C',S.DIGESTION,(b.stock('rnase','SourceRNaseA'),))
    b.clean('rna_wash','Source3xWB/EB/2xS1buffer,15mineach;retainplugs')
    b.split('plug_allocation','Selectedplugcomponent; actualplugcount/cohortunassigned',[('selected','One selectedsourceplugcomponentfor100ulS1reaction',None),('others','Otherpossibleplug/samplematerial; allocation/fateunassembled',Disposition.UNRESOLVED)])
    b.op('s1_digest','Source200US1/100ulS1buffer37C30min perplug; NOT1.8UneuronalSSBrecipe',S.DIGESTION,(b.stock('s1','SourceThermoFisherS1nuclease200U'),b.cond('s1_buffer','40mMsodiumacetatepH4.5/300mMNaCl/2mMZnSO4')),product='SourceS1cleavedDNAstructuralclassesandunreacted/offtargetDNAcopresentinplug; no100percentspecificity')
    b.op('edta','Source10mMEDTAfinal terminatesS1',S.QUENCHING,(b.cond('edta_stock','SourceEDTA10mMfinal'),))
    b.op('end_frontend','Unassembled detailedstandardENDprocessing before firstligation; sourceFig1outlineusedbelow',S.MAGIC,(b.stock('end_frontend_stock','ActualsourcehistoricalENDpreparation stockcomponent; detailsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    # Both source whole-ordering inventories are supplied once; do not infer wholefold/endstate fromthisinventory.
    a1=g.input('adapter1','SourceKeyResourcesENDadapter1'+RAW1,K.ADAPTER,state=inventory(1),chemical_state={'source_ordering':RAW1,'dU_not_T':True,'biotin_dT_count':2,'PS':'finalbond'})
    a2=g.input('adapter2','SourceKeyResourcesENDadapter2'+RAW2,K.ADAPTER,state=inventory(2),chemical_state={'source_ordering':RAW2,'dU_not_T':True,'PS':'finalbond'})
    st=g.input('end_preparation','HistoricalENDligation/capturepreparation componentstock boundary; actualstockallocationsunassembled',K.MOLECULAR);l1=g.product('ligase1','Sourcefirstligationcatalystcomponent, namedidentityunassembled',K.CATALYST);l2=g.product('ligase2','Sourcesecondadapterligationcatalystcomponent, namedidentityunassembled',K.CATALYST);g.react('end_components','Representonesourcehistoricalpreparation; no freshnamedenzymeinferred',S.GENERIC,((st,R.TEMPLATE),),(l1,l2),changes_molecule=False)
    b.op('ligation1','SourceFig1biotinylatedadapterligationtoDSBends; fullendpreparationunassembled',S.LIGATION,((a1,R.STRUCTURAL),(l1,R.CATALYST)),unknown_end(1))
    b.op('sonicate','SourceFig1sonication AFTER biotinadapterinstallation',S.FRAGMENTATION,(b.cond('sonicator','SourceENDsonicationconditionsunassembled'),))
    b.capture('streptavidin','SourceFig1streptavidinmagneticbeadcapture',(b.stock('beads','Sourcephysicalstreptavidinmagneticbeads',K.SELECTION,R.SELECTION),))
    b.op('ligation2','SourceFig1additionalIlluminaadapterinstallation; detailedhistoricalprocessingunassembled',S.LIGATION,((a2,R.STRUCTURAL),(l2,R.CATALYST)),unknown_end(2))
    b.op('end_finish','UnassembledhistoricalENDlibraryfinish/PCR/cleanup betweenknownoutlineandsequencing',S.MAGIC,(b.stock('end_finish_stock','ActualhistoricalENDfinish molecularpreparation; oligos/enzymes/dosesunassembled',K.MOLECULAR,R.STRUCTURAL),))
    return b.finish('AssociatedAMdefaultNextSeq550SE75; also150bppaired-endkits on literalMiSeq200label. CompleteENDPCR/index/readlandingunassembled. IndependentPacBioWGS/EdUcellcyclestainingexcluded.')
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import re
    c=Check();g=protocol_graph();c('PKbeforeRNase',g.steps['rna_digest'].inputs[0].material,'pk_wash_kept');c('S1washbeforecleavage',g.steps['s1_digest'].inputs[0].material,'plug_allocation_selected');c('EDTAquenchesS1',g.steps['edta'].inputs[0].material,'s1_digest_product');c('Adapterbeforefragment',g.steps['sonicate'].inputs[0].material,'ligation1_product');c('Captureafterfragment',g.steps['streptavidin'].inputs[0].material,'sonicate_product');c('Unretainedcaptureopen',g.terminals['streptavidin_unbound'].disposition,Disposition.UNRESOLVED)
    c('Adapter1actualwholepolymer', ''.join(x.top for x in inventory(1).strands['ordering'].segs),'GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGUUTUTUUACACTCTTTCCCTACACGACGCTCTTCCGATCT');c('Adapter2actualwholepolymer',''.join(x.top for x in inventory(2).strands['ordering'].segs),'GATCGGAAGAGCACACGTCUUUUUUUUAGACGTGTGCTCTTCCGATCT');c('UinventorynotN',AD1.count('U'),5)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.molcel.2022.08.007_PMC9547894.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.molcel.2022.08.007')):
        text=re.sub(r'\s+','',p.read_text()).replace('∗','*')
        c('Sourcewholemodifiedadapter1',('5’-phosphate-'+RAW1[3:]+"-3’")in text,True);c('Sourcewholemodifiedadapter2',('5’-phosphate-'+RAW2[3:]+"-3’")in text,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
