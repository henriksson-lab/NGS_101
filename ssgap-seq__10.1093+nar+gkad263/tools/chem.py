"""Private ssGap-seq: retained crosslinked chromatin conversion before reversecrosslinking."""
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

from pcr import PCRGeometry
def unknown_scene():
    s=Scene();s.strand('unknown',[Segment('unassembled','[historical source kit product unassembled]',placeholder=True)]);return s
def protocol_graph():
    b=B('ssGap-seq selectedchickenantiRecAcomponent','PublishedPMC10287960Methods59–85; publishedSI FigS1. Main200–300bp/actualRNase dose selected; SI300–500bp/.2mgml differencesretained. No wholehistoricalkitconstructorinheritance.','Source40mlOD~.5MG1655Ecoli cellcomponent, non-UV selectedbranch; genotype/proteinvariantsnotextrapolated');g=b.g
    b.op('crosslink','Source1percentformaldehyde10minRT',S.CROSS_LINKING,(b.cond('formaldehyde','Source1percentformaldehyde'),))
    b.op('quench','Source125mMglycine10minRT',S.QUENCHING,(b.cond('glycine','Source125mMglycine'),))
    b.clean('cell_wash','Source2xicecoldTBSwash/pelletstored−80C')
    b.op('resuspend','Source2mllysisbuffer; actuallysis/mixedsolubilitynotyetphysicalfractionation',S.MAGIC,(b.cond('lysis_conditions','Source20mMTris8/150mMNaCl/1mMEDTA/1percentTriton/1percentdeoxycholate/2xproteaseinhibitors'),),product='Sourcemixedcell/chromatinlysate; solubleandinsolublecomponentscopresentuntilcentrifugation')
    b.op('sonicate','SourceCovarisS2main200–300bp; SI FigS1states300–500bp, notnormalized',S.FRAGMENTATION,(b.cond('covaris','Source8x60s20percentdutyintensity8/200cyclesburst5–8C'),))
    sup=g.product('clarified_lysate','Sourcesupernatant afterinsolubleremoval',K.MOLECULAR);pel=g.product('clarify_pellet','Sourceinsolublepellet',K.MOLECULAR);g.react('clarify','Source16000g5min4C',S.CLEANUP,((b.x,R.TEMPLATE),b.cond('centrifuge','Sourcecentrifugation')),(sup,pel));g.mark_terminal(Terminal(pel.id,Disposition.UNRESOLVED,'Separatedpelletlaterfateunassigned'));b.x=sup
    b.split('lysate_aliquot','Source250ullysateselectedfrom~2mlclarifiedstock',[('selected','Source250ullysate',None),('rest','Untransferredclarifiedlysatestock; actualvolumelossunassigned',Disposition.UNRESOLVED)])
    b.op('antibody_bind','Source250ullysisbuffer+5ulchickenantiRecA4Covernight',S.BINDING,(b.stock('antibody','SourcephysicalchickenantiRecAprotein5ul',K.MOLECULAR,R.STRUCTURAL),b.cond('antibody_conditions','Source250ullysisbuffer/overnight4C')))
    b.capture('ip','Source100ulPrecipHenprewashedbeads2hRT',(b.stock('preciphen','SourceagarosegoatantichickenIgYresin100ul',K.SELECTION,R.SELECTION),))
    b.clean('ip_wash','Source2x1mleachlysis/highsalt500mMNaCl/LiCl250mM/TE; retainbeads')
    eluates=[]
    for i in [1,2]:
        e=g.product('eluate'+str(i),'Source250uleluate'+str(i),K.MOLECULAR);r=g.product('resin'+str(i),'Sourcepostelutionresin/residualmixture',K.MOLECULAR);g.react('elute'+str(i),'Source250ul.1MNaHCO3/1percentSDS5min',S.CLEANUP,((b.x,R.TEMPLATE),b.cond('elution_conditions'+str(i),'Sourceelutionbuffer')),(e,r));eluates.append(e);b.x=r
    g.mark_terminal(Terminal(b.x.id,Disposition.UNRESOLVED,'Post-secondelutionresinfateunassigned'))
    pool=g.product('eluate_pool','Actual500ulcombinedtwoeluates',K.MOLECULAR);g.react('pool_eluates','Sourcepooltwo250uleluates',S.POOL,tuple((x,R.TEMPLATE)for x in eluates),(pool,));b.x=pool
    b.clean('ip_buffer_exchange','SourceAmiconUltra10k3xTEwashthenBioSpinP6; retain~40–50ulCROSSLINKEDchromatin')
    b.op('bisulfite','Sourceadd460ul5Mbisulfite/20mMhydroquinone37C18hNONdenaturing',S.CHEMICAL_CONVERSION,(b.cond('bisulfite_solution','Fresh5Mstock/20mMhydroquinone addedsolution; finaldilutionnotcalled5M'),),product='Sourcecrosslinkedchromatin withC→UconversionenrichedinssDNA; dsDNAconversion/background/unconvertedclassespersist. ChemicalU, notreportedT.')
    b.clean('bisulfite_wash','SourceAmicon10kfourwaterwashesBEFOREdesulphonation')
    b.op('desulphonate','Source.3MNaOH20minRTdesulphonation',S.CHEMICAL_CONVERSION,(b.cond('naoh','Source.3MNaOH'),),product='SourceconvertedDNAchromatin/unconvertedclasses; no3end/wholeproductinference')
    b.clean('post_naoh','SourceAmiconfourTEwashesthenBioSpinP6/10mMTris8')
    b.op('rnase','Source1ul10ugulRNaseA37C10min; SI FigS1 insteadlabels.2mgml',S.DIGESTION,(b.stock('rnase_enzyme','SourceRNaseA10ugactualaddition; incubationvolumeunassigned'),))
    b.op('reverse_crosslink','Source200ultotal/PK.25mgml/1percentSDS65C4h',S.DEPROTEINIZATION,(b.stock('pk','SourceProteinaseK'),b.cond('reverse_conditions','SourceSDS65C4hthermalreversecrosslink+proteinasecompoundboundary')),product='Sourceproteolyzed/reversecrosslinkedconvertedDNAwithresidualproteins; completecrosslinkreversalnotassumed')
    b.clean('dna_cleanup','SourceMonarchPCR/DNAcleanup/elute20ulEB8.5/store−20C')
    b.split('library_input','Source10ulof20ulIPDNAselectedforlibrary',[('selected','Source10ullibraryinput',None),('rest','SourceuntransferredIPDNAstock',Disposition.UNRESOLVED)])
    b.op('adaptase','SourcehistoricalxGenMethylSeqAdaptaseenzymaticadapterattachment; NOsecondbisulfiteinferredfromkitname',S.MAGIC,(b.stock('adaptase_preparation','SourceactualhistoricalAdaptaseadapter/reactionmolecularpreparation; fullends/oligos/enzymesunassembled',K.MOLECULAR,R.STRUCTURAL),),product='Sourceadapter-bearingDNAplusunadapted/offclasses; fullmolecularendgeometryunassembled')
    st=g.input('pcr_preparation','SourcehistoricalxGen16cyclepreparedPCRreaction',K.MOLECULAR);cat=g.product('pcr_catalyst','PCRcatalystcomponentofsamehistoricalreaction',K.CATALYST);nt=g.product('pcr_nt','PCRnucleotidesofsamepreparedreaction',K.MOLECULAR);f=g.product('pcr_f','SourcePCRforwardprimercomponent, wholeorderingunassembled',K.OLIGO);r=g.product('pcr_r','SourcePCRreverseprimercomponent, wholeorderingunassembled',K.OLIGO);g.react('pcr_components','RepresentonesourcepreparedPCRstock; no freshspecificenzyme/NTdoseinferred',S.GENERIC,((st,R.TEMPLATE),),(cat,nt,f,r),changes_molecule=False)
    b.op('pcr','Source16PCRcycles; exactwholeprimer/runlandingunassembled',S.PCR,((cat,R.CATALYST),(nt,R.STRUCTURAL),(f,R.STRUCTURAL),(r,R.STRUCTURAL)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Historicalwholeprimer/productstatesunassembled, notunreported.'))
    b.clean('spri1','SourcefirstSPRIpurificationround');b.clean('spri2','SourcesecondSPRIpurificationround')
    b.split('library_qc','ActualQubitandqPCRQCaliquotsfromqualifiedlibrary',[('production','Remainingproductionlibrary',None),('qubit','SourceQubitQCaliquot; quantityunassigned',Disposition.QC),('qpcr','SourceqPCRQCaliquot; quantityunassigned',Disposition.QC)])
    b.split('loading','Selectedsequencingaliquot/remainder; amountunassigned',[('selected','Sourceproductionloadingaliquot',None),('rest','Untransferredfinishedlibrarystock',Disposition.UNRESOLVED)])
    lib=g.product('library','SelectedssGap-seqproductionlibraryaliquot',K.LIBRARY);g.react('load','SourceMiniSeqMidOutput300cyclePE150',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceMiniSeq2x150; exacthistoricalindexes/wholeprimerlandingunassembled. C→T/G→AreportlabelsnotchemicalU-to-Trepair. Computationaldedup/windowfiltersnotphysicaldiscard.'))
    g.add_source_item(SourceItem('source_differences','Main200–300bp vsSI300–500bp; main1ul10ugulRNaseA vsSI.2mgml.','Sourcevariantsretained; historicalkitgeometry/proteinbranchdetailsunassembled. Noabsenceaudit.'));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();c('BisulfiteBEFOREproteinase',g.steps['bisulfite'].inputs[0].material,'ip_buffer_exchange_kept');c('RNaseAFTERconversionworkup',g.steps['rnase'].inputs[0].material,'post_naoh_kept');c('PKAFTERrnase',g.steps['reverse_crosslink'].inputs[0].material,'rnase_product');c('NoDNApurificationbeforebisulfite',g.steps['dna_cleanup'].inputs[0].material,'reverse_crosslink_product');c('Secondelutionusesretainedresin',g.steps['elute2'].inputs[0].material,'resin1');c('Pooltwoproducedeluates',len(g.steps['pool_eluates'].inputs),2);c('Kitonlyselected10ul',g.steps['adaptase'].inputs[0].material,'library_input_selected');c('Noautomatedsecondbisulfite',sum(s.kind==S.CHEMICAL_CONVERSION and s.id=='bisulfite'for s in g.steps.values()),1);c('PCRnucleotidesfrompreparedstock',g.producer_of('pcr_nt').id,'pcr_components');c('LibraryQCnotwholeproduction',g.steps['loading'].inputs[0].material,'library_qc_production');c('Antibodyaddedonce',sum(any(p.material=='antibody'for p in s.inputs)for s in g.steps.values()),1);c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
