"""Private source-owned Topo-seq; publishedmain+correction+FigS1, threeoriginals."""
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

def blunt_adapter_end():
    s=Scene();s.strand('adapter',[Segment('unknown_adapter','[source Adapter2 duplex]',placeholder=True)],mod5='Unassembled',mod3='Unassembled');s.ligation_end('adapter','unknown_adapter',kind='blunt',polarity='3-prime',owner='adapter',text='SourceFigS1D Adapter2double-strand blunt end; exact terminalbase/modificationunassembled');return s
def protocol_graph():
    b=B('Topo-seq gyrase cleavage complexes','PublishedPMC6379681main70–78/50,128; FigS1originalp2. Actualgky1312correctionchangesref98only. Historicalkitwholeoligosunassembled.','Source100mlmidlogGyrA-SPAEcoli culture, notgenericnoncovalentChIP');g=b.g
    b.split('culture_split','Sourcebisection into treated/control',[('treated','Source50mldrugtreatedhalf',None),('untreated','Sourceuntreatedhalf; subsequentcontrolpreparationnotassembled',Disposition.UNRESOLVED)])
    b.op('trap','Sourceciprofloxacin/oxolinicacid/microcinB17 sourcechoiceunassigned15min37C',S.CHEMICAL_CONVERSION,(b.cond('poison','Sourcepoisonstock; sourceallowedconcentrationsinnotes, cohortunassigned'),),product='Sourcepoison-stabilizedgyrasecomplexes plusother/untrappedcellularDNA; noquantitativeperfecttrapping')
    b.clean('wash_cells','Source4500g10C/TESwashrepeatedtwice/resuspend1mlTESS')
    b.op('rna_digest','SourceRNaseAinTESS/proteaseinhibitors',S.DIGESTION,(b.stock('rnase','SourceRNaseA'),b.cond('tess','Source10mMTris7.5/1mMEDTA/250mMNaCl/.02percentSDS/.2percentTween20/proteaseinhibitors')))
    b.op('sonicate','Source200–700bpSONOPULS24x10son20soff65percentpower',S.FRAGMENTATION,(b.cond('sonicator','SourceSONOPULSHD3100'),))
    # Source names100ulmocklysate aftersonication, but positionrelativeTESdilution isunassigned.
    b.split('mock_split','Source100ulsonicatedlysate mock; relationtofollowingTESdilutionunassigned',[('production','Remainingcapturelysate; volumeunassigned',None),('mock','Source100ulmockaliquot, independentdeproteinization/extractionbelow',None)])
    production=b.x;mock=g.materials['mock_split_mock'];b.x=mock
    b.op('mock_pk','Source mockdeproteinizedasabove; exactreplicateamountunassembled',S.DEPROTEINIZATION,(b.stock('mock_proteinase','SourceProteinaseKmockreaction'),))
    b.clean('mock_extract','Sourcephenolchloroform/ethanolmockextraction');g.mark_terminal(Terminal(b.x.id,Disposition.UNRESOLVED,'ActualmockDNAcontrol retained; completesequencedmocklibraryunassembled.'));b.x=production
    b.op('dilute','Sourceaddition1mlTES',S.GENERIC,(b.cond('tes_dilution','SourceTES'),))
    b.capture('flag','Source100ulantiFLAGM2affinitygel1.5–2hRTmixing',(b.stock('flag_resin','SourceantiFLAGM2gel',K.SELECTION,R.SELECTION),))
    b.clean('wash_flag','Source2x1mlTESS/1x1mlTES/1x1mlTEretaingel')
    b.op('proteinase','Sourcebringgelto200ulTES/PK.5mgml55C≥3h',S.DEPROTEINIZATION,(b.stock('pk','SourceProteinaseK'),b.cond('pk_buffer','SourceTESproteolysisconditions')))
    sup=g.product('proteinase_supernatant','Sourcecollectedsupernatant',K.MOLECULAR);res=g.product('proteinase_pellet','Source2000g2minresidualgel/pellet',K.MOLECULAR);g.react('clarify','Source2000g2minRT; collectsupernatant',S.CLEANUP,((b.x,R.TEMPLATE),b.cond('centrifuge','Sourcecentrifugation')),(sup,res));g.mark_terminal(Terminal(res.id,Disposition.UNRESOLVED,'Residualgelnotassumedperfectlydepleted'));b.x=sup
    b.clean('extract','Sourcephenol/chloroformextractionandethanolprecipitation')
    b.op('melt','SourceFigS1CdsDNAmelting; conditionsunassembled',S.THERMAL_TREATMENT,(b.cond('melt_conditions','SourcehistoricalAccel1Smeltingconditionsunassembled'),),product='Copresentsourcefreeandpeptide-linkedssDNAchains; neitherphysicallydiscardedbymelting')
    st=g.input('library_preparation','One sourcehistoricalAccelNGS1Spreparedreactionstocks, fullallocation/compositionunassembled',K.MOLECULAR)
    td=g.product('tdt','SourceTdTcatalystcomponent',K.CATALYST);dg=g.product('dgtp','SourceGtailingnucleotidecomponent; amountunassembled',K.MOLECULAR);kl=g.product('klenow','SourceKlenowcomponent',K.CATALYST);nts=g.product('extension_nt','Sourceextensionnucleotidecomponent, historicalsupplyunassembled',K.MOLECULAR);a1=g.product('adapter1','SourceAdapter1oligoChybridizingcomponent, wholeordering/tail-lengthunassembled',K.OLIGO);a2=g.product('adapter2','SourceAdapter2dsbluntcomponent, wholeorderingunassembled',K.ADAPTER,state=blunt_adapter_end());lig=g.product('ligase','SourceAdapter2ligationcatalystcomponent, identityunassembled',K.CATALYST)
    g.react('kit_components','Representonesourcepreparedkitboundary, notnewindependentstocks',S.GENERIC,((st,R.TEMPLATE),),(td,dg,kl,nts,a1,a2,lig),changes_molecule=False)
    b.op('g_tail','SourceFigS1CTdTaddsGtofree3ends; nunknown',S.TAILING,((td,R.CATALYST),(dg,R.STRUCTURAL)),product='SourceG-tailedandunmodifiedssDNAclasses, tailbaseGbutlengthunassigned; noG6default')
    b.op('anneal1','SourceAdapter1oligoCtrackhybridizestooligoG',S.HYBRIDIZATION,((a1,R.STRUCTURAL),),product='CopresentGtailedDNA/annealedAdapter1/unmatchedclasses; noexactpairlengthclaim')
    t=g.product('retained_template','Same-tubeGtailedDNAtemplatecomponent',K.MOLECULAR);p=g.product('retained_adapter1','Same-tubepreviouslysuppliedAdapter1primercomponent',K.OLIGO);g.react('annealed_components','Nonphysicalsame-tubecomponenttracking',S.GENERIC,((b.x,R.TEMPLATE),),(t,p),changes_molecule=False);b.x=t
    b.op('fill','SourceFigS1CKlenowfillsstrandtowardsoriginal5end',S.EXTENSION,((p,R.STRUCTURAL),(kl,R.CATALYST),(nts,R.STRUCTURAL)),product='Sourcefilledduplexclasses withfreeorpeptideblockedoriginal5ends; exactproductsunassembled')
    b.op('ligate2','SourceFigS1DbluntAdapter2ligation; inefficientonpeptideblockedend',S.LIGATION,((a2,R.STRUCTURAL),(lig,R.CATALYST)),AdapterLigationGeometry.known('blunt',blunt_adapter_end(),statement='SourceS1Dexplicitbluntdouble-strandedAdapter2. Baseorder/wholejunctionunassembled; inefficientblockedclassesarenotphysicaldiscard.'))
    b.op('kit_finish','HistoricalAccelNGS1Sremaininglibraryfinish/cleanup/amplificationunassembled',S.MAGIC,(b.stock('finish_stock','Actualhistoricalkitfinish molecularpreparation, identitiesanddosesunassembled',K.MOLECULAR,R.STRUCTURAL),),product='Sourceproductionlibrary plusinefficient/offproductclassesuntilunassembledkitfinishworkup')
    return b.finish('SourceNextSeq150bppairedend; sourceAdapter2readstartFigS1E. Wholeprimers/index/runlandingunassembled, computational4bpgapnotphysicalselection.')
def selftest():
    from checks import Check
    c=Check();g=protocol_graph();c('Controlnotcapture',g.steps['dilute'].inputs[0].material,'mock_split_production');c('MockisphysicalretainedDNA',g.terminals['mock_extract_kept'].disposition,Disposition.UNRESOLVED);c('PKbeforephysicalsupernatant',g.steps['clarify'].inputs[0].material,'proteinase_product');c('NoTDP2default',any('tdp2'in x.lower()for x in g.steps),False);c('Tailbeforeprimer',g.steps['anneal1'].inputs[0].material,'g_tail_product');c('Retainedprimernotexternalreaddition',any(p.material=='retained_adapter1'and p.role==R.STRUCTURAL for p in g.steps['fill'].inputs),True);c('Nucleotidecomponentfromsamekit',g.producer_of('extension_nt').id,'kit_components');c('Adapter2blunt',g.steps['ligate2'].geometry.kind,'blunt');c('Nopeptidephysicaldiscardinligation',len(g.steps['ligate2'].outputs),1);c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
