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
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder,disposition=Disposition.UNRESOLVED):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,disposition,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

import re
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
from nextera import ME,ME_RC,ADAPTOR_S5
SOURCE_ORDERS={'21.068C-8N': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCTANNNNNNNN',
 '22tn5.003A': 'A*C*ACTCTTTCCCTACACGACGCTCTTCCGATCTGCTACACATAGATGTGTATAAGAGACAG',
 '22tn5.MOS-3p': 'CTGTCTCTTATACACATCT/3Phos/',
 '22tn5.005': '/5Phos/NNNNNNATGTGTAAATCCATTCGAGACACTGAGAGT*C*A/3Phos/',
 '22tn5.006': '/5Phos/NNNNNNATGTGTAAATCCATAGCAGCGTAAGTGGTA*G*G/3Phos/',
 '21.006G1': '/5Phos/GTTCAGACGTGTGCTCTTCCGATCTTAGCTNWWNNNNWWSWNNNWSWNNNWSWWNNNNWWNTGACTCTCAGTGTCTCGAGGACCATAATACGACTCACTATAGCAAATT',
 '21.007G1': '/5Phos/GCTCTTCCGATCTTAGCTNWWNNNNWWSWNNNWSWNNNWSWWNNNNWWNCCTACCACTTACGCTGCTTCTACTGGACCATAATACGACTCACTATAGCAATACGACCAC',
 'splint6F5': 'GTCTGAACAATTTGCT',
 'splint7F5': 'GGAAGAGCGTGGTCGT',
 '21.004G1-BC': '/5Phos/CCTACCACTTACGCNWNNNWWSWNNNWSWNNNWSWWNNNNWNCTGTCTCTTATACACATCTGACGCTGCCGACGAXXXXXXGTGACTGGAGTTCAGACGTGTGC*T*C',
 '21.004G2-BC': '/5Phos/CCTACCACTTACGCATNNWWWNNWSNNNSWNWNWNWSNWNNNNWCTGTCTCTTATACACATCTGACGCTGCCGACGAXXXXXXGTGACTGGAGTTCAGACGTGTGC*T*C',
 '21.073pt': 'CAATACGACCACGCT*C*T',
 '21.074B': '/5Phos/TGACTCTCAGTGTCTCGAGGTTGCTATAGTGAGTCGTATTAT*A*T/3Phos/',
 '21.075': 'GCTATAATACGACTCACTATAGCAA/3Phos/',
 '21.077': 'CAAGCAGAAGACGGCATACGAGATGCAATACGACCACGCTCT',
 '21.085': 'ACACTCTTTCCCTACACG',
 '21.077-G1': 'CAAGCAGAAGACGGCATACGAGATGCAATACGACCACGCTCTTCCGATCTTAGCT',
 '21.076BB': 'AATGATACGGCGACCACCGAGATCTACACTGCAACCTCGAGACACTGAGAGTCA',
 '4E4.interf1': 'TTTTTTTTTTCCTACCACTTACGC/3Phos/',
 '4E.interf2': 'TTTTTTTTTTGAGCACACGTCTGAACTC/3Phos/',
 '21.081b': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGCTA',
 '21.046G1-BC': 'CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTAGCT',
 'SBS3B': 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'}
AUTHORITY='PublishedPMC12339006.2main40–76/190–240 +publishedSI1S2/S3/S4actualPDFp6–8 +publishedTable3MOESM5actualXLSXZIP/XML;3originals. RNAmodality, no 2019/companionrecipe inheritance.'
CIRCLES={'21.006G1':'splint6F5','21.007G1':'splint7F5'}
UEI_NAMES={'21.004G1-BC','21.004G2-BC'}
def parsed(name):
    order=SOURCE_ORDERS[name];body=re.sub(r'/[^/]+/','',order);ps=tuple(len(body[:i].replace('*',''))for i,c in enumerate(body)if c=='*');return body.replace('*',''),('P'if order.startswith('/5Phos/')else''),('P'if order.endswith('/3Phos/')else''),ps

def oligo_scene(name):
    seq,m5,m3,ps=parsed(name);parts=[]
    for i,block in enumerate(re.findall('[ACGT]+|[NWSRYKMBDHV]+|X+',seq)):
        kw={};placeholder=bool(set(block)-set('ACGT'))
        if placeholder and not block.startswith('X') and name in CIRCLES:kw['feature']=feature('volumetric.umi_'+('i'if name=='21.006G1'else'ii'),'umi','degenerate')
        elif placeholder and not block.startswith('X') and name in UEI_NAMES:kw['feature']=feature('volumetric.event','molecular_complex_barcode','degenerate',note='UEI identifies joining event, not cell or sample')
        elif block.startswith('X'):kw['feature']=feature('volumetric.sample','sample_index','unknown',note='SI1S3 arbitrary6ntsamplebarcode; actualvalueunassigned')
        # Reuse exact canonical bases when a source block equals a canonical sequence.
        block=next((x for x in (P5,P7,TRUSEQ_READ1,TRUSEQ_READ2,ME,ME_RC,ADAPTOR_S5)if x==block),block)
        parts.append(Segment(name+' block'+str(i),block,placeholder=placeholder,**kw))
    s=Scene();s.strand('oligo',parts,mod5=m5,mod3=m3)
    if ps:s.note('oligo','Source phosphorothioate bonds after nucleotide '+','.join(map(str,ps)))
    if name in ('22tn5.005','22tn5.006'):s.ligation_end('oligo',parts[0].name,kind='splinted',polarity='5-prime',owner='adapter',text='Source5P/N6RNA-splintedend; occupiedRNA/junctionunassembled')
    return s

def local_circle_window(name):
    seq,_,_,_=parsed(name);splint=parsed(CIRCLES[name])[0];window=seq[-8:]+seq[:8]
    if revcomp(splint)!=window:raise ValueError('Sourcecircle/splintterminalwindow doesnotpair')
    s=Scene();s.strand('ends',[Segment('source terminal window',window)]);s.anneal('splint',[Segment('source short splint',splint)],to='ends',pair=('source short splint','source terminal window'));s.note('ends','Local last8+first8 closure window only; fullcircle/occupiedligationproduct unassembled');return s

def tn5_duplex():
    seq,_,_,ps=parsed('22tn5.003A');short=parsed('22tn5.MOS-3p')[0]
    if not seq.endswith(ME)or short!=ME_RC:raise ValueError('SourceMEinterface mismatch')
    s=Scene();s.strand('long',[Segment('source prefix',seq[:-len(ME)]),Segment('source ME',ME)]);s.anneal('short',[Segment('source reverse ME',ME_RC)],to='long',pair=('source reverse ME','source ME'),mod3='P');s.note('long','SourcePSafter '+','.join(map(str,ps))+'; shortstrand3P retained, not changedto5P');return s

def oligo_input(b,key,name,label):
    seq,m5,m3,ps=parsed(name);return b.g.input(key,label,K.OLIGO,state=oligo_scene(name),chemical_state={'source_ordering':SOURCE_ORDERS[name],'phosphorothioate_after':ps,'five_prime':m5,'three_prime':m3})
def pcr_inputs(b,key,kind):
    names=('21.046G1-BC','21.081b')if kind=='cdna'else('21.077-G1','21.076BB','4E4.interf1','4E.interf2')
    stock=b.g.input(key+'_master','ONEactualprepared'+kind+'PCRstock:PlatinumTaqHiFi/200uMdNTP/2mMMgSO4/HiFibuffer/300nMprimers'+(';3.3uMblockedinterferenceoligos/5percentDMSO'if kind=='uei'else''),K.MOLECULAR);out=[];ins=[]
    for suffix,k,r in [('enzyme',K.CATALYST,R.CATALYST),('nt',K.MOLECULAR,R.STRUCTURAL),('buffer',K.CONDITION,R.CONDITION)]:
        x=b.g.product(key+'_'+suffix,'Copresent'+suffix+'ofSAMEphysicalPCRstock',k);out.append(x);ins.append((x,r))
    for i,name in enumerate(names):
        x=b.g.product(key+'_oligo'+str(i),('3PblockedinterferenceNOTPCRprimer 'if i>=2 else'ActualPCRprimer ')+name,K.OLIGO,state=oligo_scene(name),chemical_state={'source_ordering':SOURCE_ORDERS[name],'three_prime':parsed(name)[2]});out.append(x);ins.append((x,R.STRUCTURAL))
    b.g.react(key+'_components','CopresentcomponentsofONEPCRstock, notfreshindependentadditions',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False);return tuple(ins)
def unknown_scene():
    s=Scene();s.strand('source',[Segment('WholeoccupiedPCRtemplate/productunassembled','[source molecule unassembled]',placeholder=True)]);return s

def protocol_graph():
    b=B('Volumetric DNA microscopy enhanced-preRCA embryoRNA/proximitypreparation',AUTHORITY,'Actual24hpfzebrafishembryo, upstreamcultureoutsidegraph');g=b.g
    b.op('prepare','Pronase1mgml28C5min/PFA4percent4Cnight/methanolRT15–30min/-80C>=2h/rehydrate75-50-25percent/4xPBST',S.MAGIC,(b.stock('pronase','Actualsource1mgmlPronase'),b.cond('prep_conditions','Actualsourceorderedpreparationprogramme')))
    b.op('permeabilize','ThermolabilePK5e-5UulRT12min/55C15min',S.MAGIC,(b.stock('thermo_pk','ActualsourcePK'),b.cond('pk_program','Actualprogramme')));b.clean('pk_wash','4xPBST5min')
    b.op('pre_rt','20percentformamide/Superase/DTT/albumin4C1h65C10min/4C',S.MAGIC,(b.stock('pre_rt_carriers','ActualSuperase/albumin/DTTpreparation',K.MOLECULAR,R.STRUCTURAL),b.cond('pre_rt_buffer','SourcepreRTbuffer')));b.clean('pre_rt_rinse','Onewaterrinse')
    primer=oligo_input(b,'initial_rt_primer','21.068C-8N','Actual1uMRTprimer,N8primingNOTmoleculeUMI')
    b.op('initial_rt','400uMdNTP/32uMaminoallyldUTP/SSIII10Uul4C1h60C3min37Cnight50C1h',S.MAGIC,((primer,R.STRUCTURAL),b.stock('ssiii_initial','ActualSSIII'),b.stock('initial_nt','Actual400uMdNTP+32uMaminoallyldUTPembryorecipe',K.MOLECULAR,R.STRUCTURAL),b.stock('initial_carriers','Actualalbumin/Superase/DTT',K.MOLECULAR,R.STRUCTURAL),b.cond('initial_rt_buffer','FSbuffer;32nMSummaryinitialculturedcontextnotmerged')),reason='SourceRTmarking; occupiedRNA/cDNA/modifiedpositions/anchoringyieldunassembled.')
    b.op('initial_exo','ExoI1.43Uul4C1h37C1h removesprimer/displacedcDNA',S.DIGESTION,(b.stock('exo_initial','ActualExoI'),b.cond('exo_initial_buffer','SourceExoIbuffer')));b.clean('exo_initial_wash','3xPBST/water');tissue=b.x
    long=oligo_input(b,'tn5_long','22tn5.003A','Actual100uMoligowithfirsttwoPS');short=oligo_input(b,'tn5_short','22tn5.MOS-3p','Actual100uMMEreverse3P, not5P');tn5_duplex()
    b.x=long;b.op('tn5_prepare','Anneal1:195C2min/cool1Cmin70min;25uMduplex/1ugultagmentase23C30min/50percentglycerol-20C',S.MAGIC,((short,R.STRUCTURAL),b.stock('tagmentase','Actualtagmentase'),b.cond('tn5_loading','40mMTris8/50mMNaCl')),reason='SourcePS/3Pduplexinterfacecomputed, fulloccupieddimerunassembled.')
    b.split('tn5_working','Stock3:10then1:50reactionaddition','Actualtransposomeinput','Actualuntransferredstock');tn5=b.x;b.x=tissue
    b.op('tagment','4C1h55C1h afterRT',S.MAGIC,((tn5,R.STRUCTURAL),b.cond('tagment_buffer','MgCl2/Tris/DMF/PEG8000/ATP')),reason='Actualmodifiedtransposome, wholeRNA/cDNA/end/gapdistributionunassembled.');b.clean('tagment_wash','2xPBST/PBS')
    b.op('anchor','5mMBSPEG5PBSRT1h',S.CROSS_LINKING,(b.stock('bspeg','ActualBSPEG5',K.MOLECULAR,R.STRUCTURAL),b.cond('anchor_condition','Exactcrosslinkpositions/yieldunassembled')))
    b.op('anchor_quench','1MTris8rinse/quench30min',S.QUENCHING,(b.cond('anchor_tris','Actual1MTris'),))
    b.op('tn5_denature','4xSSC4C10min70C15min/formamideSSC4C10min50C10min',S.MAGIC,(b.cond('denature_program','Sourceorderedprogramme'),));b.clean('denature_wash','3xPBST5min/water')
    ad1=oligo_input(b,'adapter_i','22tn5.005','Actual500nMtypeIadapter5P/3P/PS/N6bindingwindow');ad2=oligo_input(b,'adapter_ii','22tn5.006','Actual500nMtypeIIadapter5P/3P/PS/N6bindingwindow')
    b.op('adapters','SplintR1.25Uul4C1h23Cnight',S.MAGIC,((ad1,R.STRUCTURAL),(ad2,R.STRUCTURAL),b.stock('splintr','ActualSplintR'),b.cond('splintr_buffer','SplintRbufferwithATP')),geometry=AdapterLigationGeometry.known('splinted',oligo_scene('22tn5.005'),oligo_scene('22tn5.006'),statement='Localsource5P/N6RNAsplintadapterinterfaces only;3P/PSretained, completeoccupiedRNA/targetDNAjunctionunassembled.'),reason='Source3primeadapterligation; N6bindingNOTidentifierUMI; wholeproductsunassembled.');b.clean('adapter_wash','3xPBST5min/water')
    b.op('cip','QuickCIP.5Uul4C1h37C1h removesligatedadapter3P',S.DEPHOSPHORYLATION,(b.stock('quick_cip','ActualQuickCIP'),b.cond('cip_buffer','rCutSmart')),product='Intendedadapter3OHclass plusresidual/unreactedclasses; full5primejunction/enddistributionunassembled');b.clean('cip_wash','3x2xSSCT5min');anchored=b.x
    circle_inputs=[]
    for i,(linear,splint)in enumerate(CIRCLES.items()):
        key='circle'+str(i+1);local_circle_window(linear);b.x=oligo_input(b,key+'_linear',linear,'Actualsource5PlinearUMIcircleprecursor '+linear);sp=oligo_input(b,key+'_splint',splint,'Actualsource16ntshortsplint '+splint)
        b.op(key+'_ligate','SourceT4ligase/shortsplintmakesCirc6G1orCirc7G1',S.MAGIC,((sp,R.STRUCTURAL),b.stock(key+'_ligase','ActualT4DNAligase'),b.cond(key+'_conditions','Detailedcircleligationdose/programmeunassembled')),reason='ActualTable3mappedcircleprecursor/splint; local8+8windowpairs, fulloccupiedproduct/alias24versus21equivalenceunassembled.')
        b.clean(key+'_clean','SourceZymoOligoConcentratorspincolumn');b.split(key+'_qc','ActualTBEureagelQCaliquot','Actualproductioncirclestock','ActualanalyticalQCaliquot',Disposition.QC);b.split(key+'_anneal_input','Actual100nMcircleinputselectedfrompreparedstock','Actualcircleannealinput','Actualuntransferredcirclestock');circle_inputs.append(b.x)
    b.x=anchored;b.op('circle_anneal','100nMeachCirc6G1/7G1hybridize40Cnight',S.HYBRIDIZATION,tuple((x,R.STRUCTURAL)for x in circle_inputs)+(b.cond('circle_hyb','2xSSC/10percentformamide/.1percentTween'),));b.clean('circle_wash','40Chybridizationwash/2xSSCT/1xSSCT/PBST')
    b.op('pre_rca','SelectENHANCEDpreRCA25ngulg32/1mMdNTP/.04UulT4DNApol',S.MAGIC,(b.stock('pre_rca_pol','ActualT4DNApolymerase'),b.stock('pre_rca_nt','Actual1mMdNTP',K.MOLECULAR,R.STRUCTURAL),b.stock('pre_rca_g32','ActualT4gene32protein',K.MOLECULAR,R.STRUCTURAL),b.cond('pre_rca_buffer','NEBuffer2.1; normal100uMdNTPalternativeoutsidegraph')),reason='Sourceinitialonce-copyUMIextensionbeforeRCA, fullheterogeneousoccupiedproductsunassembled.')
    b.op('rca','Phi29.2Uul/d(AUGC)TP250uM/g32/albumin4C1h30Cnight',S.MAGIC,(b.stock('phi29','Actualphi29'),b.stock('rca_nt','ActualsourceDNAprecursorsdA/dU/dG/dC, NOTsubstituteddT',K.MOLECULAR,R.STRUCTURAL),b.stock('rca_carriers','Actualg32/albumin',K.MOLECULAR,R.STRUCTURAL),b.cond('rca_buffer','Sourcephi29buffer, optionalfluoresceindUTPimagingnotadded')),reason='SourceanchoredtandemUMIDNAnanoballs; fullconcatemer/size/diffusion/yieldunassembled, notdeterministicspatialprediction.');b.clean('rca_wash','3x2xSSCT')
    uei=[oligo_input(b,'uei_'+str(i),name,'ActualsourceUEI/bridgingoligo '+name)for i,name in enumerate(('21.004G1-BC','21.004G2-BC','21.073pt','21.074B'))]
    b.op('uei_anneal','SourceG1/2-BCmix100nM/073pt100nM/074B100nM4C1h50C2h',S.HYBRIDIZATION,tuple((x,R.STRUCTURAL)for x in uei)+(b.cond('uei_hyb','2xSSC/5percentformamide/.1percentTween'),));b.clean('uei_wash','50C1hhybwash/2xSSCT/1xSSCT/PBST/water')
    b.op('uei_extend_ligate','T4DNApol.15Uul/T4ligase20Uul/1mMdNTP+ATP4C1hRT40min',S.MAGIC,(b.stock('uei_pol','ActualT4DNApol'),b.stock('uei_ligase','ActualT4ligase'),b.stock('uei_nt','Actual1mMdNTP/ATP sourcecomponents',K.MOLECULAR,R.STRUCTURAL),b.cond('uei_buffer','T4ligasebuffer')),reason='SourceRCA-UEI/intermediate/cDNA-UMIclasses; wholeheterogeneousligation/extensionproductsandactualproximityfidelityunassembled.');b.clean('uei_post_wash','3xPBST/water')
    promoter=oligo_input(b,'promoter','21.075','Actual100nMT7promoteroligo3P, notextendablePCRprimer')
    b.op('ivt','30ulUSER/T7IVT/PEGgelRT2hthen37C20h',S.MAGIC,((promoter,R.STRUCTURAL),b.stock('user','Actual.05UulUSER'),b.stock('t7','Actual10percentT7EnzymeMix'),b.stock('ivt_rntp','Actual7.5mMrNTPmix',K.MOLECULAR,R.STRUCTURAL),b.stock('ivt_g32','Actual100ngulg32',K.MOLECULAR,R.STRUCTURAL),b.stock('hydrogel','Actual6.4ugulETTMP/73.6ugul4armPEGvinylsulfone',K.MOLECULAR,R.STRUCTURAL),b.cond('ivt_buffer','ActualsourceT7buffer/gelprogramme')),reason='USERdigestsdUnanoballs releasingRCA-UEI/intermediate/UMI-cDNAtemplates;T7transcriptionandunanchoredaRNAextensioncreateIVT-UEIs. No fullrelease/diffusion/recombinationyieldor spatialproof.',product='SourceDNA/RNAmixture includingRCA-UEIaRNAs,intermediateUEIaRNAs,UMI-cDNAaRNAsandIVT-UEIaRNAs; intendedT7RNA5PPPclassFigS4, fullproductsunassembled')
    b.op('gel_denature','Add12ulKOH457.5mM/EDTA100mM/DTT42.5mM4C2h',S.MAGIC,(b.cond('gel_denature_buffer','Actualsource12uldenaturationsolution'),))
    b.op('gel_stop','Add12ulTris600mM7.5/.4NHCl→54ultotal',S.QUENCHING,(b.cond('gel_stop_buffer','Actual12ulstopsolution'),))
    b.op('post_gel_pk','Add30ulPKmix→84ul50C1h',S.DEPROTEINIZATION,(b.stock('post_gel_pk_stock','ActualPK.09Uulmix'),b.cond('post_gel_pk_buffer','.28percentTween/8.6mMTris7.5')));b.clean('rna_isolate','1.2xRNAClean/elutewater')
    b.op('dnase','DNaseI.1Uul/Superase.8Uul37C30min',S.DIGESTION,(b.stock('dnase_stock','ActualDNaseI'),b.cond('dnase_buffer','ActualSuperase/DNaseIbuffer')));b.clean('rna_reclean','1.2xRNAClean/elutewater')
    rt1=oligo_input(b,'ivt_rt_uei','21.077','Actual500nMUEI-RTprimer');rt2=oligo_input(b,'ivt_rt_cdna','21.085','Actual500nMcDNA-RTprimer')
    b.op('ivt_rt_anneal','Primers/dNTP500uM65C5min/ice BEFOREenzymeaddition',S.THERMAL_TREATMENT,((rt1,R.STRUCTURAL),(rt2,R.STRUCTURAL),b.stock('ivt_rt_nt','Actual500uMdNTP',K.MOLECULAR,R.STRUCTURAL),b.cond('ivt_rt_anneal_condition','Sourceprogramme')))
    b.op('ivt_rt','AddSSIII10Uul/Superase/5mMDTT/FSbuffer50C1h70C15min',S.MAGIC,(b.stock('ivt_ssiii','ActualSSIII'),b.stock('ivt_rt_inhibitor','ActualSuperase',K.MOLECULAR,R.STRUCTURAL),b.cond('ivt_rt_buffer','FSbuffer/DTT')),reason='SourceRNApostDNase/carryingactualexternalprimers/NTs fromanneal; nofreshprimer/NTaddition orwholeproductsproof.')
    b.op('ivt_exo','ExoI3.3Uul37C30min80C20min',S.DIGESTION,(b.stock('ivt_exo_stock','ActualExoI'),b.cond('ivt_exo_condition','Sourceprogramme')))
    allocations=[g.product('cdna_template_'+str(i),'ActualpostRTstockaliquotforcDNAPCR'+str(i)+' atfinal1:80dilution',K.MOLECULAR)for i in range(4)]+[g.product('uei_template_'+str(i),'ActualpostRTstockaliquotforUEIPCR'+str(i)+' atfinal1:40dilution',K.MOLECULAR)for i in range(2)]
    rem=g.product('template_stock_remaining','ActualuntransferredpostRTstock, amount/fateunassembled',K.MOLECULAR);g.react('allocate_pcr_templates','Source4cDNA+2UEIPCRallocationsfromONEpostRTstock',S.SPLIT,((b.x,R.TEMPLATE),),tuple(allocations)+(rem,));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'UntransferredstocknotclonedintosixPCRs.'))
    products=[]
    for i,template in enumerate(allocations):
        kind='cdna'if i<4 else'uei';key=kind+'_pcr'+str(i if i<4 else i-4);b.x=template
        program='95C2min/5cycles95C30s56C30s68C2min/20cycles95C30s68C2min/68C5min'if kind=='cdna'else'95C2min/1cycle95C30s66C30s68C2min/18cycles95C30s68C2min/68C5min'
        b.op(key,'Source'+kind+'PCR '+program,S.PCR,pcr_inputs(b,key,kind),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Wholeoccupiedprimerlanding/heterogeneousampliconproductsunassembled; sourcewholeprimers/blockersretained.'));products.append(b.x)
    clean_cdna=g.product('clean_cdna_library_stock','ActualsourcecleanedcDNAampliconfamily stock; 4PCRreplicatehandlingunassembled',K.MOLECULAR);clean_uei=g.product('clean_uei_library_stock','ActualsourcecleanedUEIampliconfamily stock;2PCRreplicatehandlingunassembled',K.MOLECULAR);discard=g.product('pcr_cleanup_other','ActualexcludedPCRcleanupfractions',K.WASTE)
    g.react('pcr_cleanup','Source.75xAMPurecleanup ofPCRproducts; actualreplicatepool/transfer scheduleunassembled',S.MAGIC,tuple((x,R.TEMPLATE)for x in products)+(b.stock('ampure_stock','ActualsourceAMPureXPallocations',K.SELECTION,R.SELECTION),),(clean_cdna,clean_uei,discard),reason='Source4+2PCRtubes and2readoutclasses retained; exactreplicatepooling/fraction-transfer amountsunassembled, notassertedphysicalpoolrecipe.');g.discard(discard,'Actualcleanupfractions; no perfectyield/purityclaim.')
    custom=oligo_input(b,'custom_sbs3b_stock','SBS3B','Actualcustomspike-insequencingprimerstock, exactrunallocationunassembled');g.mark_terminal(Terminal(custom.id,Disposition.UNRESOLVED,'SourcecustomSBS3Baddedtosequencerrun; physicalper-library/poolallocationsoutsidepartialgraph, notinventedtwodoses.'))
    for kind,stock in (('cdna',clean_cdna),('uei',clean_uei)):
        b.x=stock;b.split(kind+'_quant','Actualsourceanalyticalquantificationaliquot, method/amountunassembled','Actualremainingproductionlibrary','Actualanalyticalaliquotnotreturned',Disposition.QC)
        b.split(kind+'_load_input','ActualNGSloadingcontribution, pooling/amountunassembled','ActualNGSinput','Actualuntransferredlibrary stock')
        lib=g.product(kind+'_library','Actualsource'+kind+'NGSlibraryclass, exactreplicatepool/primerlandingunassembled',K.LIBRARY);g.react(kind+'_load','SourceNextSeq500150cyclekit;112R1/44R2 andSBS3Bcustomspike',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':112,'R2':44},limitation='Source150cyclekit/112R1/44R2 retained asreported; customSBS3B='+ADAPTOR_S5+'. WholeheterogeneousUMI/UEI/cDNAproducts/readlanding/indexmask/replicatepoolingunassembled. Spatialinferencecomputed, notwetlabvalidation.'))
    g.add_source_item(SourceItem('context_and_names','SelecteddetailedembryoRT32uMaminoallyldUTP, Summaryinitialculturedoptimization32nM, contextsnotmerged. SI1S1/S2labels24.006/007 versusTable3labels21.006/007 retained; localfactoryusesexplicitTable3Circ6/7mapping.','No wholealias/versionequivalence/publicationerrorclaim. Source3PshortTn5/adapters/T7promoter/interferenceoligos/PS preserved; RCA dU/USER releases precedeIVTdiffusion.'));return g
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as ET
    c=Check();g=protocol_graph()
    c('Pronaseexplicitphysicalsourceport',g.steps['prepare'].inputs[1].material,'pronase');c('Adapterlocalendskind',g.steps['adapters'].geometry.kind,'splinted');c('Tn5short3P',g.materials['tn5_short'].chemical_state['three_prime'],'P');c('Tn5PSpositions',g.materials['tn5_long'].chemical_state['phosphorothioate_after'],(1,2));c('CIPafterligation',g.steps['cip'].inputs[0].material,'adapter_wash_kept');c('RCAafteroncecopy',g.steps['rca'].inputs[0].material,'pre_rca_product');c('USERIVTafterUEIassembly',g.steps['ivt'].inputs[0].material,'uei_post_wash_kept');c('PostRTafterDNasecleanup',g.steps['ivt_rt_anneal'].inputs[0].material,'rna_reclean_kept');c('Actual4plus2allocations',len(g.steps['allocate_pcr_templates'].outputs),7);c('SixphysicalPCRs',sum(x.kind==S.PCR for x in g.steps.values()),6);c('OneactualPCRstock',g.steps['cdna_pcr0_components'].inputs[0].material,'cdna_pcr0_master');c('Interference3Pretained',g.materials['uei_pcr0_oligo2'].chemical_state['three_prime'],'P');c('TwoactualNGSclasses',set(g.sequencing),{'cdna_library','uei_library'})
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41587-025-02613-z_supp_41587_2025_2613_MOESM5_ESM.xlsx'
    if have(c,Source(p,'https://doi.org/10.1038/s41587-025-02613-z')):
        z=zipfile.ZipFile(p);ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};ss=[''.join(e.itertext())for e in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)];actual={}
        for row in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:row',ns):
            cells={}
            for cell in row.findall('s:c',ns):
                v=cell.find('s:v',ns)
                if v is not None and cell.attrib.get('t')=='s':cells[''.join(x for x in cell.attrib['r']if x.isalpha())]=ss[int(v.text)]
            if cells.get('B')in SOURCE_ORDERS:actual[cells['B']]=''.join(cells.get('C','').split())
        for name,order in SOURCE_ORDERS.items():c('Originalwholeordering '+name,actual.get(name),''.join(order.split()))
    c('Partialnotpublicready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
