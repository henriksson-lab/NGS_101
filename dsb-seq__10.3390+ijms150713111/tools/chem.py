"""Private2014DSB-specific TdT/biotin/S1 route, locally source-qualified."""
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

from illumina import TRUSEQ_READ1,P7
FW=TRUSEQ_READ1[:22]
RV=P7[:22]+'C'
LITERAL_TE_ML=200
AUTHORITY='Published2014PMC4139894main3.4DSB/common3.5/3.1/Fig1b/5/3.7; actualPDFp9 verified200mL/Fw/Rv. One original, not later2018chapter.'
def oligo(seq):
    s=Scene();s.strand('primer',[Segment('whole source primer',seq)]);return s

def source_a_end():
    s=Scene();s.strand('top',[Segment('unknown local body','XXXX',placeholder=True,note='Schematic body only; full sequence/length unassigned'),Segment('source protruding dA','A')],mod3='OH')
    s.anneal('bottom',[Segment('same body complement','xxxx',placeholder=True)],to='top',pair=('same body complement','unknown local body'))
    s.overhang('top','source protruding dA',polarity='3-prime',owner='substrate',text='Source3primeA local substrate end only; whole fragment/adapter ends unassembled');return s

def protocol_graph():
    b=B('DSB-Seq2014 source-specific biotin/TdT preparation',AUTHORITY,'Source500ug preparedHCT116HMWgDNA; common3.1/3.2upstream extraction/ETOexposure outside graph',K.MOLECULAR);g=b.g;production=b.x
    # Separate actual500ug control preparation; never a fresh duplicate of production output.
    b.x=g.input('control_dna','Actual separate500ugHMWgDNA source negative-label control preparation',K.MOLECULAR,primary=True)
    b.op('control_tail','Source negative-label tailing: substituteTTP forbiotin16dUTP',S.MAGIC,(b.stock('control_tdt','SourcecontrolTdT24000U',K.CATALYST,R.CATALYST),b.stock('control_dctp','Sourcecontrol0.5mMdCTP',K.MOLECULAR,R.STRUCTURAL),b.stock('control_ttp','SourceTTP substitutesbiotin16dUTP, exactsourceorderingretained',K.MOLECULAR,R.STRUCTURAL),b.cond('control_conditions','Source same3mlTdTbuffer37C30min exceptlabelsubstitution')),
       reason='Positive source negative-label control; whole unlabelled tails/end distributions unassembled, not a zero-background assumption.')
    b.op('control_workup','Source negative-label selection/recovery control; exact downstream processing/quantification allocation unassembled',S.MAGIC,(b.stock('control_preparation','Actual source parallelworkup/selectionpreparation; internalstocks/fractionsunassembled',K.MOLECULAR,R.STRUCTURAL),),reason='SourceResults2.2/Fig2 define a negative control of labeling andselection, but complete processing/terminalS1 allocation ofcontrol material unassembled.')
    g.mark_terminal(Terminal(b.x.id,Disposition.QC,'SourceFig2 recovered-controlDNAquantification; not silently a productionsequencinglibrary.'));b.x=production
    b.op('tdt','Source500ugDNA/3mlRochebuffer/24000UTdT37C30min',S.MAGIC,(b.stock('tdt_enzyme','SourceTdT24000U',K.CATALYST,R.CATALYST),b.stock('dctp','Source0.5mMdCTP tailing substrate',K.MOLECULAR,R.STRUCTURAL),b.stock('biotin_dutp','Source0.02mMBiotin16dUTP covalenttail label',K.MOLECULAR,R.STRUCTURAL),b.cond('tdt_conditions','Source3mlRocheTdTbuffer37C30min')),
       reason='Source3primeendtailing positive; occupiedDNA3ends, mixed C/biotinU tail lengths and residual/unreacted classes unassembled.',product='Source biotin-nucleotide-containing3tails plus unreactedDNA/freeNTPs; no quantitative DSB specificity/yield claim')
    b.op('tdt_stop','SourceEDTA20mM stop',S.QUENCHING,(b.cond('edta','Source20mMfinalEDTA'),))
    b.clean('tail_clean','Sourcephenolchloroform/two2Mammoniumacetateethanolprecipitations; literal200mLTEresuspension(PDFverified), intendedvolumeunresolved')
    b.op('sonicate','SourceBioruptor30x30smedium/30siceintervals to200–400bp',S.FRAGMENTATION,(b.cond('sonication','Source physical sonication conditions'),))
    bound=g.product('streptavidin_bound','Source kilobaseBINDER bead-retained preparation; exactbound/freeDNA classes unassigned',K.MOLECULAR);flow=g.product('streptavidin_unbound','Actual nonretained selection fraction; laterfateunassigned',K.MOLECULAR)
    g.react('capture','SourceDynabeadskilobaseBINDER streptavidin capture',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('streptavidin_beads','Source streptavidin-coatedkilobaseBINDERbeads',K.SELECTION,R.SELECTION)),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Nonretained fraction not implicitly sequenced.'));b.x=bound
    b.clean('capture_wash','Source4washes50degrees+4RT with10mMTrispH7.5/1mMEDTA/2MNaCl')
    b.op('release','Source disruptbiotinstreptavidin75C4h in2MbetaME',S.MAGIC,(b.cond('release_buffer','Source10mMTrispH7.5/1mMEDTA/1MNaCl/2Mbetamercaptoethanol'),),reason='Source positive freeDNA release; complete protein/bead/biotinlinker molecular products and recovery unassembled.')
    b.clean('release_clean','SourceQIAquickPCRpurification offreeDNA')
    b.op('s1','Source30US1nuclease/110ulrecommendedbuffer37C30min',S.DIGESTION,(b.stock('s1_enzyme','SourceFermentasS1nuclease30U',K.CATALYST,R.CATALYST),b.cond('s1_conditions','Source110ulrecommendedbuffer37C30min')),product='Source tail-removal DNA classes plus remaining/incompletelycleaved/damagedDNA and tailfragments; exact molecular ends unassembled')
    b.clean('s1_clean','SourceQIAquickPCRpurification afterS1')
    b.op('end_repair','Local source3.5 EpicentreENDRepair45minRT',S.END_REPAIR,(b.stock('end_enzymes','Source1ulT4DNApolymerase+T4PNK mix',K.CATALYST,R.CATALYST),b.stock('end_dntps','Source0.25mMeachdNTP structuralpool',K.MOLECULAR,R.STRUCTURAL),b.stock('end_atp','Source1mMATP',K.MOLECULAR,R.STRUCTURAL),b.cond('end_buffer','Source33mMTrisacetate/66mMpotassiumacetate/10mMmagnesiumacetate/.5mMDTT')))
    b.clean('end_clean','Source purification afterend-repair; exactinternalworkupunassembled')
    b.op('a_tail','Source15UKlenowexo-minus/.2mMdATP37C30min',S.TAILING,(b.stock('klenow','Source15UKlenowexo-minus',K.CATALYST,R.CATALYST),b.stock('datp','Source.2mMdATP covalentsubstrate',K.MOLECULAR,R.STRUCTURAL),b.cond('atail_conditions','Source37C30min')),product='Source protruding3A candidateDNA classes plus incomplete/unreactedDNA andstocks; fullboth-endDNAunassembled')
    geometry=AdapterLigationGeometry('unavailable','unavailable',(source_a_end(),unknown_scene()),'Known source local substrate3A preserved; wholehistoricalIlluminaadapter/strandends/sealedproducts unassembled, not source-sequenceabsence.')
    b.op('adapter_ligate','Source.1ulIlluminaadaptormix/1000UT4ligaseRT30min',S.LIGATION,(b.stock('adapters','Sourceactualhistorical.1ulIlluminaadaptormix, wholeoligos/endsunassembled',K.ADAPTER,R.STRUCTURAL),b.stock('ligase','Source1000UT4DNAligase',K.CATALYST,R.CATALYST),b.cond('ligation_conditions','SourceRT30min')),geometry=geometry)
    b.clean('postligation','SourceQIAquickPCRpurification')
    b.clean('pre_pcr_gel','Source2percentEGel200–400bpband/MinElutegelcleanup/elute12ul; rejectedgel physicallyseparate')
    stock=g.input('pcr_stock','Sourceprepared18cyclePCRstock, exactquantitiesunassembled',K.MOLECULAR);f=g.product('fw','SourcewholeFw',K.OLIGO,state=oligo(FW));r=g.product('rv','SourcewholeRv23nt',K.OLIGO,state=oligo(RV));nts=g.product('pcr_nts','SamepreparedPCR nucleotidecomponent',K.MOLECULAR);enzyme=g.product('pcr_enzyme','SamePCR enzymecomponent, identityunassembled',K.CATALYST);buf=g.product('pcr_buffer','SamepreparedPCRconditions',K.CONDITION)
    g.react('pcr_components','Same source preparedPCRstockcomponents, not freshduplicateadditions',S.GENERIC,((stock,R.TEMPLATE),),(f,r,nts,enzyme,buf),changes_molecule=False)
    b.op('pcr','Source18cycles98C30s/65C30s/72C30s',S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(nts,R.STRUCTURAL),(enzyme,R.CATALYST),(buf,R.CONDITION)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Wholehistoricaladapter/occupiedprimerlanding/completeDNAproducts unassembled; actualsourcewholeFW/RV retained.'))
    b.clean('post_pcr_gel','Source2.5percentagarose220–500bpband/Qiagengelpurification')
    b.split('load_aliquot','Sourceclusterloadingallocation, exactamountunassigned','ActualGenomeAnalyzerloadingaliquot','Untransferredsourcepurifiedlibrarystock, fateunassigned')
    lib=g.product('library','SourceDSBSeq2014 selectedproductionlibrary',K.LIBRARY);g.react('load','SourceIlluminaGenomeAnalyzerclustergeneration/sequencing',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Source3.7 aligned36bpreads; completephysicalrunmask/primerlanding/indexlayoutunassembled. TdTbiotin/S1classesnotbaseexactDSBcoordinates. Negative-labelcontrolnotmerged withproduction.'))
    g.add_source_item(SourceItem('source_discrepancy','Source3.4actualPDFp9literal200mLTE retained, intendedvolumeunresolved. Methods3.1ETO20min versusFig5caption5min bothretained.','No silentunit/dose repair; no fullsource/reaction/absenceaudit.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Biotintailing precedesphysicalshear',g.steps['sonicate'].inputs[0].material,'tail_clean_kept')
    c('S1usesreleasedpurifiedDNA',g.steps['s1'].inputs[0].material,'release_clean_kept')
    c('LibraryrepairfollowsS1',g.steps['end_repair'].inputs[0].material,'s1_clean_kept')
    c('PCRusesactualselectedband',g.steps['pcr'].inputs[0].material,'pre_pcr_gel_kept')
    c('ControlisactualindependentDNA',g.steps['control_tail'].inputs[0].material,'control_dna')
    c('Controlisnotproductionreads',g.terminals['control_workup_product'].disposition,Disposition.QC)
    c('Sourceunitnotnormalized',LITERAL_TE_ML,200)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'ijms150713111_PMC4139894.xml.txt'
    if have(c,Source(p,'https://doi.org/10.3390/ijms150713111')):
        t=''.join(p.read_text().split());c('ActualwholeFwtranscription',FW in t,True);c('ActualwholeRvtranscription',RV in t,True);c('Sourcebiotin/S1inventory','Biotin-16-dUTP'in t and '30UofS1'in t,True);c('Literal200mLsource','200mLofTEbuffer'in t,True)
    c('Privatepartial remainsunready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
