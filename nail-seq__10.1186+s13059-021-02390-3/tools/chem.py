"""Private NAIL source-owned sequential-label and two-fraction library graph."""
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

from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
AUTHORITY='PublishedPMC8188667main169–183; publishedSI1FigS1c/d actualp1/captionp2; originalAdditionalfile5TableS4 workbook. Three originals only.'
CORE='CCACGCGTGCTCTACA'
BARCODE='NNNTNNNNANNNNTNNN'
BA1_UP=CORE+BARCODE+'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGT'
BA1_DOWN=revcomp(CORE)+'N'*6
BA2_UP=TRUSEQ_READ2[3:]
BA2_DOWN='N'*6+revcomp(BA2_UP)
P5_I5=P5+TRUSEQ_READ1[:22]
P7_I7=P7+'N'*6+TRUSEQ_READ2[:23]
SOURCE_ORDERINGS={'Bridge Adapter-1-up':"5'P-"+BA1_UP+'-NH2','Bridge Adapter-1-down':"5'-"+BA1_DOWN+'-NH2','Bridge Adapter-2-up':"5'-"+BA2_UP,'Bridge Adapter-2-down':"5'-"+BA2_DOWN+'-NH2','P5-I5':"5'-"+P5_I5,'P7-I7-Index':"5'-"+P7_I7}
def adapter(which):
    s=Scene()
    if which==1:
        s.strand('up',[Segment('duplex core',CORE),Segment('source inline barcode',BARCODE,placeholder=True,feature=feature('nail.bridge_barcode','inline_barcode','degenerate',note='Source calls bridge barcoded; identity/UMI use unassigned')),Segment('source tail',BA1_UP[len(CORE+BARCODE):])],mod5='P',mod3='NH2')
        s.anneal('down',[Segment('paired core',BA1_DOWN[:-6]),Segment('random splint','N'*6,placeholder=True,note='Random splint; not automatically an identifier')],to='up',pair=('paired core','duplex core'),unpaired=('random splint',),mod3='NH2')
        s.ligation_end('up','duplex core',kind='splinted',polarity='5-prime',owner='adapter',text='Source Bridge1up5P; adjacent N6 splint supplied, occupied DNA3end/full junction unassembled')
    elif which==2:
        s.strand('up',[Segment('duplex core',BA2_UP)],mod3='OH — unmodified source DNA ordering')
        s.anneal('down',[Segment('random splint','N'*6,placeholder=True,note='Random splint; not automatically an identifier'),Segment('paired core',BA2_DOWN[6:])],to='up',pair=('paired core','duplex core'),unpaired=('random splint',),mod3='NH2')
        s.ligation_end('up','duplex core',kind='splinted',polarity='3-prime',owner='adapter',text='Bridge2 unmodifiedDNA3OH interpretation; N6splint present, occupied substrate5P/junction unassembled')
    else:raise ValueError('Two source bridge types only')
    return s

def primer(seq):
    s=Scene()
    if seq==P7_I7:
        segs=[Segment('P7',P7),Segment('source i7 index','N'*6,placeholder=True,feature=feature('nail.i7','sample_index','unknown',note='Actual study index unassigned')),Segment('source priming tail',TRUSEQ_READ2[:23])]
    elif seq==P5_I5:segs=[Segment('P5',P5),Segment('fixed source tail',TRUSEQ_READ1[:22],note='No variable i5 window is printed despite P5-I5 name')]
    else:raise ValueError('Source primers only')
    s.strand('primer',segs);return s

def protocol_graph():
    b=B('NAIL-seq published sequential EdU/BrdU K562 preparation',AUTHORITY,'Source cultured K562 cells; whole upstream cell-culture preparation outside graph');g=b.g
    b.op('arrest','Source5uMpalbociclib36h G1 arrest',S.GENERIC,(b.cond('arrest_conditions','SourceK5625uMpalbociclib36h'),))
    b.clean('release_wash','Source prewarmedRPMI1640wash and2.5h freshmedium release')
    b.op('edu_pulse','Source10uMEdU15min cellular pulse',S.MAGIC,(b.stock('edu','Source suppliedEdU nucleoside10uM',K.MOLECULAR,R.STRUCTURAL),b.cond('edu_conditions','15min sourcepulse')),
        reason='Positive metabolic incorporation; intracellular nucleotide intermediates and whole nascent-DNA classes unassembled. No fresh in-vitro EdUTP added.')
    b.clean('between_labels','Source prewarmedRPMI1640wash between two pulses')
    b.op('brdu_pulse','Source50uMBrdU15min subsequent pulse',S.MAGIC,(b.stock('brdu','Source suppliedBrdU nucleoside50uM',K.MOLECULAR,R.STRUCTURAL),b.cond('brdu_conditions','15min sourcepulse')),
        reason='Sequential cellular labeling; exact fork incorporation/spatial occupancies and nucleotide intermediates unassembled.')
    b.clean('harvest','Source harvest/PBSwash')
    b.op('extract','Source gDNA extraction with ProteinaseK perref70',S.EXTRACTION,(b.stock('proteinasek','SourceProteinaseK',K.CATALYST,R.CATALYST),b.cond('extraction_conditions','Historicalref70 extraction internal stocks/fractions unassembled')))
    b.op('click_fragment','Source Click25C1h, copper hydrolyzesgDNA to~200–400bp',S.CHEMICAL_CONVERSION,
        (b.stock('biotin_peg_azide','Source100uMbiotinPEGazide molecular label',K.MOLECULAR,R.STRUCTURAL),b.stock('click_copper','Source2mMCuSO4',K.CATALYST,R.CATALYST),b.cond('click_conditions','Source4mMTHPTA/10mMsodiumascorbate;25C1h')),
        product='Source biotinEdU/BrdU/unlabelled fragmentedgDNA classes, free label and residual chemical stocks; length/end distributions unassembled, not purelabelpartition')
    b.clean('precipitate','Source DNA precipitation after Click reaction')
    b.op('denature','Source95C5min DNA denaturation',S.THERMAL_TREATMENT,(b.cond('denature_conditions','95C5min'),))
    kinase_stock=g.input('pnk_stock','Actual prepared source kinase reaction stock; exact source donor/formulation unassembled',K.MOLECULAR)
    kinase=g.product('pnk_enzyme','T4PNK component of same prepared stock',K.CATALYST)
    kinase_rest=g.product('pnk_rest','Other components of same source kinase stock; donor identity/quantity unassembled',K.MOLECULAR)
    g.react('pnk_components','Represent components of same actual source kinase stock, no duplicate fresh additions',S.GENERIC,((kinase_stock,R.TEMPLATE),),(kinase,kinase_rest),changes_molecule=False)
    b.op('pnk','SourceT4PNK37C30min',S.MAGIC,((kinase,R.CATALYST),(kinase_rest,R.STRUCTURAL),b.cond('pnk_conditions','37C30min')),
        reason='Source positive phosphorylation; known intended5P class, complete donor/end distributions and full molecular states unassembled, no source-absence claim.',product='Source phosphorylated5P candidateDNA plus incompletelytreated DNA and remaining stocks; not quantitativeconversion')
    brdu=g.product('brdu_bound','Source antiBrdU/ProteinG bead-retained fraction; mixed/incomplete/off-target classes unassigned',K.MOLECULAR)
    sn=g.product('brdu_supernatant','Actual source supernatant from same antiBrdU capture, not a new externalDNAinput',K.MOLECULAR)
    g.react('brdu_capture','SourceantiBrdU(BU1/75)/ProteinG captureRT4h',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('brdu_antibody','SourceBU1/75',K.SELECTION,R.SELECTION),b.stock('protein_g','SourceProteinGbeads',K.SELECTION,R.SELECTION)),(brdu,sn))
    b.x=sn;b.op('supernatant_denature','Source actualsupernatant denatured again',S.THERMAL_TREATMENT,(b.cond('sn_denaturation','Source re-denaturation; exact repeated programmeunassembled'),))
    edu=g.product('edu_bound','SourceC1 bead-retained EdU-labelled fraction from actualBrdU-supernatant',K.MOLECULAR);flow=g.product('c1_flowthrough','SourceC1 non-retainedfraction; composition/fateunassigned',K.MOLECULAR)
    g.react('edu_capture','SourceStreptavidinC1 captureRT4h',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('c1_beads','SourceStreptavidinC1',K.SELECTION,R.SELECTION)),(edu,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Nonretained fraction not implicitly sequenced.'))
    for name,start in [('brdu',brdu),('edu',edu)]:
        b.x=start
        bridges=tuple(b.stock(name+'_bridge'+str(i),'SourceBridgeAdapter'+str(i)+'1uM wholeduplex stock',K.ADAPTER,R.STRUCTURAL,state=adapter(i),chemical_state={'source_orderings':tuple(SOURCE_ORDERINGS['Bridge Adapter-'+str(i)+'-'+end]for end in ['up','down'])})for i in (1,2))
        geometry=AdapterLigationGeometry.known('splinted',adapter(1),adapter(2),statement='Partial source-owned adapter interfaces only: Bridge1up5P and unmodified Bridge2up3OH with adjacent N6splints. Known3NH2blocks retained. Full occupied targetDNA3ends, exact6ntbinding and sealed products remain unassembled; not a whole ligation or quantitative-product proof.')
        b.op(name+'_ligate','Source simultaneousonbeadBridge1/2ligationRT>=4h',S.LIGATION,bridges+(b.stock(name+'_ligase','Source1600UT4DNAligase',K.CATALYST,R.CATALYST),b.cond(name+'_ligation_conditions','Source15percentPEG8000;1uMeachbridge; actualbufferdonorstockunassembled')),geometry=geometry)
        b.clean(name+'_wash','Source carefullywashligatedDNA')
        stock=g.input(name+'_tagging_stock','ActualsourceP5/P7taggingpreparation perref70; whole chemistry/formulationunassembled',K.MOLECULAR)
        f=g.product(name+'_p5','Source fixedP5-I5 primer component fromsamepreparation',K.OLIGO,state=primer(P5_I5),chemical_state={'source_ordering':SOURCE_ORDERINGS['P5-I5']})
        r=g.product(name+'_p7','SourceP7-I7N6primer component fromsamepreparation',K.OLIGO,state=primer(P7_I7),chemical_state={'source_ordering':SOURCE_ORDERINGS['P7-I7-Index']})
        rest=g.product(name+'_tagging_rest','Other components ofsametaggingstock; internalidentitiesunassembled',K.MOLECULAR)
        g.react(name+'_tagging_components','Samepreparedstockcomponents, not freshduplicateprimers',S.GENERIC,((stock,R.TEMPLATE),),(f,r,rest),changes_molecule=False)
        b.op(name+'_tag','SourceP5-I5/P7-I7tagging delegatedtoref70',S.MAGIC,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(rest,R.STRUCTURAL)),reason='Positive source librarytagging; full operationsequence/PCRproducts/cycles/primersites unassembled, not inherited fromanotherprotocol.')
        lib=g.product(name+'_library','Source'+name+' NAIL library, full molecules/loading allocationunassembled',K.LIBRARY)
        g.react(name+'_load','SourceHiSeq2x150readout; complete loadingallocationunassembled',S.MAGIC,((b.x,R.TEMPLATE),),(lib,),reason='Positive sourceHiSeqreadout; exactloadingaliquot/indexmask/wholeDNAunassembled.')
        g.sequence(lib,SequencingPlan(limitation='SourceHiSeq2x150; actualsourceP7N6 sampleindex, fixedprintedP5-I5, inlinebridgebarcodepurposeunassigned. Exactprimersites/indexmask/fullDNAunassembled. BrdU/EdUlibrariesnotphysicallymerged; E-B iscomputationalsignalsubtraction.'))
    g.add_source_item(SourceItem('unassembled','Mainref70 extraction/tagging historicaloriginal remainsordinaryunread.','Selected3originals positive only; no complete source or source-absence audit.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Click occurs on extractedDNA',g.steps['click_fragment'].inputs[0].material,'extract_product')
    c('No invented sonication stage','sonicate' in g.steps,False)
    c('PNK follows clickedprecipitatedDNA',g.steps['pnk'].inputs[0].material,'denature_product')
    c('Secondcapture usesactualfirstsupernatant',g.steps['supernatant_denature'].inputs[0].material,'brdu_supernatant')
    c('BrdUligation usesitsboundfraction',g.steps['brdu_ligate'].inputs[0].material,'brdu_bound')
    c('EdUligation usesactualC1boundfraction',g.steps['edu_ligate'].inputs[0].material,'edu_bound')
    c('No physical EminusB pool',any(s.kind==S.POOL for s in g.steps.values()),False)
    c('Two distinct library endpoints',set(g.sequencing),{'brdu_library','edu_library'})
    c('PrintedP5 has no variableindex','N' in P5_I5,False)
    c('PrintedP7indexwindow exactly6',P7_I7.count('N'),6)
    c('Bridge1has3primeblockbothstrands',sum(bool(s.mod3=='NH2')for s in adapter(1).strands.values()),2)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s13059-021-02390-3_supp_13059_2021_2390_MOESM5_ESM.xlsx.txt'
    if have(c,Source(p,'https://doi.org/10.1186/s13059-021-02390-3')):
        t=''.join(p.read_text().split())
        for name,order in SOURCE_ORDERINGS.items():c('Original modifiedwholeordering '+name,''.join(order.split()) in t,True)
    c('Privatepartial remainsunready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
