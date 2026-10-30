"""Private partial Paired-Damage-seq: qualified AM recipe, published Table1 orderings.

One explicit R01 tube plus source-qualified parallel remainder; both final modalities.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from nextera import ME,ME_RC,S5
from illumina import TRUSEQ_READ1,P5,P7
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
TITLE='Paired-Damage-seq source-qualified joint preparation'
AUTHORITY='AssociatedAM PMC12182908 main59/125–185/223/227; finalNatureHTMLsubscriptionpreview/maps; publishedMOESM3 originalXLSX SupplementaryTable1 sheets1–6 selectedorderings. Threeoriginals; AMnotclaimedfullyequivalenttofinalMethods.'
PREFIX='AGGCCAGAGCATTCG'
R01_CODES=('CATC','ATGA','AGCT','ACAG','GAAT','TACG','TTAC','GTTG','CCGT','CGAA','TCTA','GGGC')
NOTI='GCGGCCGC';SBFI='CCTGCAGG'
DNA1=PREFIX+'A'+R01_CODES[0]+NOTI+ME
RT_T=PREFIX+'T'+R01_CODES[0]+SBFI+'T'*16+'VN'
RT_N=PREFIX+'T'+R01_CODES[0]+SBFI+'N'*6
TSO='AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G'
LINK2='CGAATGCTCTGGCCTCTCAAGCACGTGGAT'
BLOCK2=revcomp(LINK2)
LINK3='GGTCTGAGTTCGCACCGAAACATCGGCCAC'
QUENCH3=revcomp(LINK3)
PA_F='CAGACGTGTGCTCTTCCGATCT';PA_R='AAGCAGTGGTATCAACGCAGAGT'
ANCHOR=PA_R+'GAAGGATGT'+'G'*9+'*H'
R02='/5Phos/GTGCGAACTCAGACCAAACAAACATCCACGTGCTTGAG'
R03=PA_F+'N'*10+'AAACAAACGTGGCCGATGTTTCG'
P5_TOP=TRUSEQ_READ1;P5_BOTTOM='NNDC'+revcomp(TRUSEQ_READ1)
P5H_TOP=TRUSEQ_READ1+'H';P5H_BOTTOM='NNDCD'+revcomp(TRUSEQ_READ1)
P5_INDEX=P5+'NNNNNNNN'+TRUSEQ_READ1[:-1]+'*'+TRUSEQ_READ1[-1]
N5_INDEX=P5+'NNNNNNNN'+S5
P7_INDEX=P7+'NNNNNN'+'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'
def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('unassembled','[sourcewholegeometryunassembled]',placeholder=True)]);return s
def oligo_scene(seq,*,index=None,mod5='',mod3=''):
    s=Scene()
    if index:
        prefix,n,tail,role=index;parts=[Segment('fixedprefix',prefix),Segment('identifier','N'*n,placeholder=True,feature=feature('pd_'+('i5'if prefix==P5 else'i7'),role,'unknown')),Segment('fixedtail',tail)]
    else:parts=[Segment('source oligo',seq,placeholder=bool(set(seq)-set('ACGTU')))]
    s.strand('oligo',parts,mod5=mod5,mod3=mod3);return s
def dna_scene():
    s=Scene();s.strand('adapter',[Segment('R01fixedprefix',PREFIX+'A'),Segment('R01cellbarcode','N'*4,placeholder=True,feature=feature('pd_r01','cell_barcode','unknown')),Segment('NotI',NOTI),Segment('ME',ME)],mod5='P');return s
def pments_scene():
    s=Scene();s.strand('oligo',[Segment('sourceMEcomplement',ME_RC),Segment('additionalterminalddC','C',note='Source appended ddC is an additional nucleotide, not a replacement/modification of the preceding T')],mod5='P',mod3='ddC');return s
def tn5_duplex():
    s=dna_scene()
    s.anneal('pments',[Segment('sourceMEcomplement',ME_RC),Segment('additionalterminalddC','C',note='Source twentieth nucleotide; kept outside assertedME pairing, no dropped/mismatch-repaired base')],to='adapter',pair=('sourceMEcomplement','ME'),unpaired=('additionalterminalddC',),mod5='P',mod3='ddC');return s

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('nuclei','Source3.6million permeabilizedHeLa nuclei containingDNA/nativeRNA; celltreatmenttimepointunassigned',K.BIOLOGICAL,primary=True)
    def stock(key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(key,label,kind,**kw),role)
    def cond(key,label):return stock(key,label,K.CONDITION,R.CONDITION)
    def op(key,label,kind,ins=(),product='Source copresentDNA/RNA/reagentmixture; residualsubstrates retained until physicalworkup',geometry=None):
        nonlocal x
        kw={}
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.MAGIC:kw['reason']='Source-specific parallel/plate composition or complete reaction geometry remainsunassembled, not assertedabsent.'
        if geometry is not None:kw['geometry']=geometry
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);a=g.product(key+'_kept','Source physically retainedpreparation; recovery/carryoverunassigned',K.MOLECULAR);w=g.product(key+'_waste','Source separatedworkupfractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(a,w));g.discard(w,'Source physicalworkup, no perfect molecular specificity.');x=a
    def pcr(key,label,primers):
        ins=[stock(key+'_enzyme','SourcefreshKAPAHiFiHS'if key!='dna_index'and key!='rna_index'else'SourcefreshNEBQ5'),stock(key+'_nts','Sourcefresh'+('1ul10mMdNTPs'if key in{'dna_index','rna_index'}else'0.5ul10mMdNTPs')+'; per-base/total conventionunassigned',K.MOLECULAR,R.STRUCTURAL)]
        for j,name in enumerate(primers):
            kw={}
            if name.startswith('SourceP5XX '):kw={'state':oligo_scene(P5_INDEX.replace('*',''),index=(P5,8,TRUSEQ_READ1,'sample_index')),'chemical_state':{'source_ordering':P5_INDEX,'phosphorothioate_after':(len(P5_INDEX)-3,)}}
            elif name.startswith('SourceN5XX '):kw={'state':oligo_scene(N5_INDEX,index=(P5,8,S5,'sample_index')),'chemical_state':{'source_ordering':N5_INDEX}}
            elif name.startswith('SourceP7XX '):kw={'state':oligo_scene(P7_INDEX,index=(P7,6,'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC','sample_index')),'chemical_state':{'source_ordering':P7_INDEX}}
            ins.append(stock(key+'_primer'+str(j),name,K.OLIGO,R.STRUCTURAL,**kw))
        op(key,label,S.PCR,ins,geometry=PCRGeometry.unavailable(unknown_scene(),reason='Sourcewholeconvertedmixture/primerproduct geometry unassembled; not confirmedmissing.'))
    # Independent source stock assembly, chosen1ul and explicitremainder.
    da=g.input('dna_adapter','SourcebarcodedDNAadapterfamily, selectedtubeidentityunassigned; Table1 DNA#01 wholeordering '+DNA1,K.OLIGO,state=dna_scene(),chemical_state={'source_representative':'/5Phos/'+DNA1,'source_role':'R01cellbarcode, exactsample/tubeunassigned','source_codes':R01_CODES})
    pm=g.input('pments','Source5Phos/'+ME_RC+'ddC',K.OLIGO,state=pments_scene(),chemical_state={'five_prime':'P','three_prime':'ddC','additional_terminal_modified_base':(len(ME_RC),'ddC'),'source_ordering':'5Phos/'+ME_RC+'ddC'})
    ac=g.input('anneal_condition','Source50uMeach,95C5min/cool20C-0.1Csec',K.CONDITION);duplex=g.product('annealed_adapter','Source1ul selectedannealedDNAadapter/pMENTs aliquot; wholepreparedstockvolume/remainderunassigned',K.MOLECULAR,state=tn5_duplex())
    g.react('anneal','SourceindependentDNAadapterannealing',S.HYBRIDIZATION,((da,R.TEMPLATE),(pm,R.STRUCTURAL),(ac,R.CONDITION)),(duplex,))
    enzyme=g.input('unloaded_patn5','Source6ulunloadedproteinA-Tn50.5mgml',K.MOLECULAR);lc=g.input('load_condition','Source1ulannealedadapter+6ulproteinA-Tn5,RT30min4C10min; labelled0.5mgmlcomplexstock/dilutionbridgeunassigned',K.CONDITION);loaded=g.product('loaded_stock','Source7ulpreparedbarcodedproteinA-Tn5stock',K.MOLECULAR)
    g.react('load_tn5','Source actual independenttransposomeloading',S.BINDING,((enzyme,R.TEMPLATE),(duplex,R.STRUCTURAL),(lc,R.CONDITION)),(loaded,))
    aliquot=g.product('loaded_aliquot','Source1ulpreparedcomplexselectedforthisR01tube',K.MOLECULAR);rem=g.product('loaded_remainder','Untransferredpreparedtransposomestock,nominal6ul;laterfateunassigned',K.MOLECULAR);g.react('loaded_split','Source1ulof7ulstockaliquot',S.SPLIT,((loaded,R.TEMPLATE),),(aliquot,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Untransferredstockisnotaddedtothesample.'))
    ab=g.input('biotin_ab','Source10ulanti-biotin0.2mgml',K.SELECTION);bc=g.input('ab_load_condition','Source20ulMedBuffer1RT60minrotation',K.CONDITION);abcomplex=g.product('ab_complex','Preparedsourceanti-biotin/proteinA-Tn5/adaptercomplex',K.SELECTION)
    g.react('ab_load','Source antibodyprecomplexwiththeONEselectedloadedstockaliquot',S.BINDING,((aliquot,R.TEMPLATE),(ab,R.SELECTION),(bc,R.CONDITION)),(abcomplex,))
    op('fix','Source2.5mlNIB-HEPES+123ul37percentformaldehyde,RT10min',S.CROSS_LINKING,(cond('fix_condition','Sourceformaldehyde'),))
    clean('fix_wash','Source600g4C8min,washcoldNIB-Tris,pellet10min')
    op('nucleosome_deplete','Source800ulr2.1/0.3percentSDS42C30min300rpm',S.PERMEABILIZATION,(cond('deplete_condition','SourceSDSnucleosomedepletion'),),product='SourcefixeddepletednucleiwithDNA/RNAandreleasedmaterialbeforewashing')
    clean('deplete_wash','Source600g10min,2x200ulNIB-Trisand1x200ulr2.1/0.1percentTritonwash')
    op('incise','Source84ulbuffer+8ulFpg+8ulEndoIV37C16h',S.DIGESTION,(stock('fpg','SourceFpg'),stock('endoiv','SourceEndoIV')),product='SourceDNAwithprocessedlesions/AP/SSBendsandresiduallesions,nativeRNAretained; noideal lesion-removalguarantee')
    clean('incise_spin','Source600g10minretainnuclei,resuspend79ulNEBuffer3labelingmixture')
    nts=stock('biotin_nts','Source4ulmixture5uMeachdATP/dCTP/dGTP/biotin-dUTP; dTTPpre-repaircontrolnotselected',K.MOLECULAR,R.STRUCTURAL,chemical_state={'nucleotides':('dATP','dCTP','dGTP','biotin-dUTP')})
    op('repair_label','Sourceadd5ulBstfull-length+10ulTaqDNAligase+2ulNAD+4ulnucleotides,37C6h',S.EXTENSION,(stock('bst','SourceBstfull-length'),stock('taq_ligase','SourceTaqDNAligase'),nts,cond('repair_condition','SourceNAD/labelingbuffer')),product='Sourcebiotin-dUrepair-incorporatedDNAandresidualnicked/unprocessedDNA; freebiotinnucleotides/nativeRNAcopresentuntilwash')
    clean('pretarget_wash','Source600g4C10min,twice200ulMedBuffer1wash')
    chosen=g.product('r01_nuclei','Sourceoneoftwelve300k-nucleialiquots,50ulMedBuffer1',K.MOLECULAR);others=g.product('other_r01_nuclei','Other11sourcealiquotsfromsame3.6millionnucleipreparation',K.MOLECULAR);g.react('r01_split','Source12tubes;oneexplicitrouteplusparallelremainder',S.PARTITION,((x,R.TEMPLATE),),(chosen,others));x=chosen
    op('target_bind','Sourceaddpreparedantibody-transposomecomplex,4Covernightrotation',S.BINDING,((abcomplex,R.SELECTION),))
    clean('target_wash','Source600g4C10min,twice50ulMedBuffer2withoutEDTA')
    template=g.product('tag_template','Same-tubeDNA/RNAnuclearcomponentafterbinding/wash',K.MOLECULAR);cat=g.product('retained_tn5','Same-tube retainedproteinA-Tn5fromloaded_aliquot,notnewstock',K.CATALYST);ad=g.product('retained_adapter','Same-tubeloadedadaptercomponentfromloaded_aliquot',K.ADAPTER,state=tn5_duplex())
    g.react('retained_components','Representcopresentcomponentsofretainedloadedstock;no physicalfractionationornewaddition',S.GENERIC,((x,R.TEMPLATE),),(template,cat,ad),changes_molecule=False);x=template
    op('tag','Source2ul250mMMgCl237C60min550rpm',S.TAGMENTATION,((cat,R.CATALYST),(ad,R.STRUCTURAL),cond('tag_condition','SourceMedBuffer2/Mg')),product='SourcetargetedtagmentedDNA/R01adapterproductsandresidualDNA/nativeRNA; exactloadedadapter/genomicendpaths unassembled')
    op('tag_stop','Source16.5ul40.5mMEDTA',S.QUENCHING,(cond('edta','SourceEDTA'),));clean('tag_spin','Source1000g4C10min;retainnucleiforRT')
    rtprimers=[stock('rt_t','PublishedTable1 RNA#01_RE representative '+RT_T,K.OLIGO,R.STRUCTURAL,chemical_state={'source_ordering':'/5Phos/'+RT_T,'priming':'sourceT16VN;AMlabelT15retainedasversionqualified','R01_role':'cellbarcode'}),stock('rt_n','PublishedTable1 RNA#01_NRE representative '+RT_N,K.OLIGO,R.STRUCTURAL,chemical_state={'source_ordering':'/5Phos/'+RT_N,'priming':'randomN6,notUMI','R01_role':'cellbarcode'}),stock('tso','SourceTSO '+TSO,K.OLIGO,R.STRUCTURAL,chemical_state={'source_ordering':TSO,'three_prime':'rG/rG/LNAG'}),stock('maxima','SourceMaximaReverseHMinus1Uul'),stock('rt_nts','SourceRT0.5mMdNTPs;per-base/total conventionunassigned',K.MOLECULAR,R.STRUCTURAL)]
    op('rt','Source20ulRTmix,T15label/N6each2.5uM,TSO2uM;42C90min thenpublishedAMstagedtemperatureprogramme169; sourcecompletecDNA/template-switchgeometryunassembled',S.REVERSE_TRANSCRIPTION,rtprimers,product='SourcecDNA/R01primers,originalRNA,tagmentedDNAandretainedRTreagents; noRNAdegradationinvented')
    other_ready=g.product('other_r01_ready','Sourceother11tubesafterparallelR01tagging/RT;completeotherstockassembly/geometryunassembled',K.MOLECULAR);g.react('other_r01','Actualparallelother11R01routes,notnewbiologicalmaterial',S.MAGIC,((others,R.TEMPLATE),),(other_ready,),reason='Sourceother11tube-specifictransposome/RTpreparationsknownbutnotfullyassembledinthisselectedpartialgraph.')
    pooled=g.product('r01_pool','Sourceall12R01tubespooled,cooled2min/Triton4.8ul5percent',K.MOLECULAR);g.react('pool_r01','SourcephysicalR01poolbeforecombinatorialligation',S.POOL,((x,R.TEMPLATE),(other_ready,R.TEMPLATE)),(pooled,));x=pooled;clean('r01_pool_spin','Source1000g4C10minretainnuclei')
    for n,raw,link in [(2,R02,LINK2),(3,R03,LINK3)]:
        plate=stock(f'r0{n}_plate',f'SourcepreparedR0{n}barcodedadapterplates,5ulworkingstockperwell;whole50ulstock/remainingnineworkingplatesnotassembled',K.ADAPTER,R.STRUCTURAL,chemical_state={'source_representative':raw,'linker':link,'preparation':'6ul100uMbarcodedoligo+44ul12.5uMlinker,95C5min/cool20C-0.1Csec;ten5ulworkingplates','cell_barcode':'combinatorialR0'+str(n),'R03_random_N10':'UMIassociationunassembled; sourceanalysisusesUMI, noautomaticassignmenttoeveryN'})
        op(f'ligate_r0{n}',f'Sourcesplitintotwo96wellR0{n}plates,20ulreactionperwell/5ulpreparedadapter;freshT4ligase37C30min;exactpercellplateallocation/wholejunctionsunassembled',S.LIGATION,(plate,stock(f'ligase_r0{n}','SourcefreshT4DNAligase')),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'SourceR0'+str(n)+'ligationknown; wholeplate/splint/genomicgeometryunassembled,notconfirmedmissing.'))
        op(f'block_r0{n}','SourceR02blocking30min'if n==2 else'SourceR03terminationwithoutfurtherincubation',S.QUENCHING,(stock(f'blocker_r0{n}',BLOCK2 if n==2 else QUENCH3,K.OLIGO,R.STRUCTURAL),cond(f'block_condition{n}','SourceR02ligasebuffer vsR03EDTAtermination')))
        clean(f'pool_r0{n}',f'SourcephysicalR0{n}pool/pellet;retainnuclei,notcomputationalsorting')
    selected=g.product('selected_sublibrary','Source5000–10000nucleisublibraryin35ulPBS',K.MOLECULAR);other_sub=g.product('other_sublibraries','Otherpooledbarcodednuclei/sublibraries;amount/fateunassigned',K.MOLECULAR);g.react('sublibrary','Sourceactualsublibraryaliquoting',S.SPLIT,((x,R.TEMPLATE),),(selected,other_sub));g.mark_terminal(Terminal(other_sub.id,Disposition.UNRESOLVED,'Unselectednucleiarephysicalremainder,notdiscardedreads.'));x=selected
    op('pk','Source5uleach4MNaCl/10percentSDS/10mgmlProteinaseK55C2h850rpm',S.DEPROTEINIZATION,(stock('pk_enzyme','SourceProteinaseK'),))
    clean('lysate_clean','Source1xSPRI/elute12.5ul;DNA/cDNAretained')
    op('tail','Source12.5ulDNA/cDNA+1.5ul10xTdTbuffer+0.5ul1mMdCTP,95C5min/ice5min;add1ulTdT37C30min75C20min',S.TAILING,(stock('tdt','SourceTdT'),stock('dctp','SourcedCTP',K.MOLECULAR,R.STRUCTURAL)),product='SourceC-tailedDNA/cDNAplusuntail/incompleteclasses,freepolymerase/nucleotidesretained')
    op('anchor','SourcefreshAnchorMix6ul5xKAPAbuffer/0.6uldNTP/0.6ul10uManchor/0.6ulKAPAHS;16cyclesAMprogramme177',S.EXTENSION,(stock('anchor_primer','Source'+ANCHOR,K.OLIGO,R.STRUCTURAL,chemical_state={'source_ordering':ANCHOR,'three_prime':'degenerateH=A/C/T;terminalPS'}),stock('anchor_kapa','SourcefreshKAPAHiFiHS'),stock('anchor_nts','Sourcefresh0.6ul10mMdNTPs;per-base/total conventionunassigned',K.MOLECULAR,R.STRUCTURAL)))
    pcr('preamp','Sourcefresh4ul5xbuffer/0.5uldNTP/2ul10uMPA-F/PA-R/0.5ulKAPA;5cycles98C20s65C20s72C2.5min',('SourcePA-F '+PA_F,'SourcePA-R '+PA_R))
    clean('preamp_clean','SourceSPRI0.2x10ul+0.65x32.5uldouble-size-selection/elute35ul; sourceinternalfractionrecipeunassembled, no defaultkitfatesclaimed')
    common=x;dna=g.product('dna_tube','Source17ulmixedDNA/cDNAintoDNAconstructiontube',K.MOLECULAR);rna=g.product('rna_tube','Source17ulmixedDNA/cDNAintoRNAconstructiontube',K.MOLECULAR);unused=g.product('split_remainder','Sourceuntransferrednominal1ulmixedeluate',K.MOLECULAR);g.react('modal_split','Physical35ul→17ulDNA+17ulRNA+nominal1ulremainder,NOTmolecularspeciespurification',S.SPLIT,((common,R.TEMPLATE),),(dna,rna,unused));g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Untransferredeluate.'))
    for modality,start in [('dna',dna),('rna',rna)]:
        x=start;enzyme='SbfI-HF'if modality=='dna'else'NotI-HF'
        op(modality+'_digest1','Source'+enzyme+'37C60min/rCutSmart',S.DIGESTION,(stock(modality+'_digest1_enzyme','Source'+enzyme),),product='Sourcecleaved/noncleavedDNA/cDNAclassesstillcopresent; digestdoesnotphysicallyremovetheothermodality')
        clean(modality+'_clean1','Source1.8xSPRI/elute11ul')
        pcr(modality+'_second_pcr','SourcefreshKAPA/PA-F+'+('PA-R,5cycles'if modality=='dna'else'TSOprimer,7cycles')+'98C20s65C20s72C2.5min',('SourcePA-F '+PA_F,'Source'+('PA-R '+PA_R if modality=='dna'else'TSOprimer '+PA_R)))
        clean(modality+'_clean2','Source1.8xSPRI/elute'+('20ul'if modality=='dna'else'17ul'))
        op(modality+'_digest2','Source'+('SbfI-HF+FokI'if modality=='dna'else'NotI-HF')+'37C60min/rCutSmart',S.DIGESTION,tuple(stock(modality+'_cut_'+str(j),'Source'+e)for j,e in enumerate(('SbfI-HF','FokI')if modality=='dna'else('NotI-HF',))))
        clean(modality+'_clean3','Source1.25xSPRI/elute10ul')
        if modality=='dna':
            mix=stock('p5_adapter_mix','Source2ulTOTALreadyP5/P5Hcomplexmix1:3;NNDC/D/HsourcebaseclassesnotUMIs',K.ADAPTER,R.STRUCTURAL,chemical_state={'top':P5_TOP,'bottom':'5Phos/'+P5_BOTTOM,'h_top':P5H_TOP,'h_bottom':'5Phos/'+P5H_BOTTOM,'source_preparation':'anneal50uMeach95C5min/cool20C-0.1Csec,mix1:3;selectedstockandremainingstockquantityunassembled'})
            op('dna_adapter2','Source2ulP5mix+2ul10xligasebuffer+4ulwater+2ulT4ligase;4C10min10C15min16C15min25C45min',S.LIGATION,(mix,stock('dna_ligase','SourcefreshT4DNAligase')),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'FokI/sourceNNDCmixed-adapterend/junctiongeometryunassembled,notconfirmedmissing.'))
            clean('dna_adapter_clean','Source1.25xSPRI/elute30ul')
        else:
            op('rna_tagment','Source10ul2xtagmentbuffer+0.5ul0.05mgmlTn5-N5,37C30min550rpm',S.MAGIC,(stock('rna_tn5_loaded','SourcesinglepreparedTn5-N5loadedstock; separateadapter/enzymecomponentsunassembled',K.MOLECULAR,R.STRUCTURAL),),product='SourceN5-taggedcDNAandcopresentresidualDNA/classes;notnativeRNAtransposition')
            clean('rna_tag_clean','SourceQIAquickPCRpurification/elute30ul0.1xelutionbuffer')
        pcr(modality+'_index','Source30ultaggedproduct+10ul5xQ5buffer/1uldNTP/0.5ul50uM'+('P5'if modality=='dna'else'N5')+'/2.5ul10uMP7/5ulwater/1ulQ5;72C5min98C30s,8–10x98C10s63C30s72C1min,72C1min',('Source'+('P5XX '+P5_INDEX if modality=='dna'else'N5XX '+N5_INDEX),'SourceP7XX '+P7_INDEX))
        clean(modality+'_final_clean','Source0.85x42.5ulSPRIcleanup')
        lib=g.product(modality+'_library','Source'+modality+'librarypreparation',K.LIBRARY);g.react(modality+'_load','Sourcesequencerloading',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
        g.sequence(lib,SequencingPlan(limitation='AMsourceNovaSeqX100/8/8/100; exactpublishedP7XXorderingcontainsN6whileindexrun8cycles,landingunassemblednotautomaticallycontradictory. Chemicalbiotin-dUrepairmarkisnotanUMI; R01/R02/R03jointcellidentityknown, exactUMI/wholecyclegeometryunassembled. Computationalfiltersnotphysicalloss.'))
    g.add_source_item(SourceItem('pending','FinalmainMethods subscriptionaccessbarrier;AMrecipeapplicabilitytoentirefinalversionunverified. Fullparallelplate/stocks/end/readgeometry andsourceUMIwindowmappingunassembled.','Source-positive source-qualified private scopeonly,noconfirmedabsence.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    from zipfile import ZipFile
    from xml.etree import ElementTree as E
    c=Check();g=protocol_graph()
    c('pMENTswholemodeledpolymerretainsadditionalddC',''.join(s.top for s in g.materials['pments'].state.strands['oligo'].segs),ME_RC+'C')
    c('LoadedMEinventoryretainsadditionalddC',''.join(s.top for s in tn5_duplex().strands['pments'].segs),ME_RC+'C')
    c('SourcepMENTsterminalmodifiedbaseposition',g.materials['pments'].chemical_state['additional_terminal_modified_base'],(len(ME_RC),'ddC'))
    c('Bothlabelingcatalystscoadded',[p.material for p in g.steps['incise'].inputs if p.role==R.CATALYST],['fpg','endoiv'])
    c('Selectedrepairpoolretainsbiotin-dU',g.materials['biotin_nts'].chemical_state['nucleotides'],('dATP','dCTP','dGTP','biotin-dUTP'))
    c('Loadedstockexplicitselectedaliquot',g.steps['loaded_split'].inputs[0].material,'loaded_stock')
    c('Unusedloadedstocknotconsumed',g.terminals['loaded_remainder'].disposition,Disposition.UNRESOLVED)
    c('Tagmentusesproducedsame-stockcomponents',[p.material for p in g.steps['tag'].inputs if p.role in{R.CATALYST,R.STRUCTURAL}],['retained_tn5','retained_adapter'])
    c('NoexternaladapterreadditionatTAG',all(g.producer_of(p.material) is not None for p in g.steps['tag'].inputs if p.role==R.STRUCTURAL),True)
    c('RTconsumessourcemolecularnucleotides',[(i.material,i.role)for i in g.steps['rt'].inputs if i.material=='rt_nts'],[('rt_nts',R.STRUCTURAL)])
    c('Anchorextensionhasitsfreshsourcenucleotides',[(i.material,i.role)for i in g.steps['anchor'].inputs if i.material=='anchor_nts'],[('anchor_nts',R.STRUCTURAL)])
    c('RTaftertagmentpelleting',g.steps['rt'].inputs[0].material,'tag_spin_kept')
    c('Other11tubesretainprimaryprovenance',g.steps['other_r01'].inputs[0].material,'other_r01_nuclei')
    c('ActualR01poolusesproducedparallelpreparation',[p.material for p in g.steps['pool_r01'].inputs],['rt_product','other_r01_ready'])
    c('Modalbranchisphysicalaliquotnotpurification',g.steps['modal_split'].kind,S.SPLIT)
    c('Untransferred1ulexplicit',g.terminals['split_remainder'].disposition,Disposition.UNRESOLVED)
    c('DNAbranchSbfIbeforeFokI',g.steps['dna_digest1'].inputs[1].material,'dna_digest1_enzyme')
    c('RNAroutetagmentssourcetaggedcDNA',g.steps['rna_tagment'].inputs[0].material,'rna_clean3_kept')
    path=sources_root()/Path(__file__).resolve().parents[1].name/'s41592-025-02632-3_41592_2025_2632_MOESM3_ESM.xlsx'
    if have(c,Source(path,'https://doi.org/10.1038/s41592-025-02632-3')):
        with ZipFile(path)as z:
            shared=[''.join(e.itertext())for e in E.fromstring(z.read('xl/sharedStrings.xml'))]
            for name,seq in [('pMENTs','5Phos/'+ME_RC+'ddC'),('DNA1','/5Phos/'+DNA1),('RTT','/5Phos/'+RT_T),('RTN','/5Phos/'+RT_N),('TSO',TSO),('Anchor',ANCHOR),('R02',R02),('R03',R03),('P5',P5_INDEX),('N5',N5_INDEX),('P7',P7_INDEX)]:c('Wholepublishedsourceordering '+name,seq in shared,True)
    c('Privateincompletewholepreparation',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
