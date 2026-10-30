"""Published LISH Plus-slide branch; source RNA2 reporter and opaque amplification."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from illumina import P5,P7,TRUSEQ_READ2,INDEX1_PRIMER
OF='AGATCGGAAGAGCACAC'
OR='GGAGCTGTCGTTCACTC'
LEFT_TARGET_DNA='TAGAACCAGTTCTCATCG'
LEFT_ORDERING=OR+LEFT_TARGET_DNA+'rUrA'
RIGHT_TARGET='GGGAGCAAGCTCTTTGTGCT'
RIGHT_ORDERING='/5Phos/'+RIGHT_TARGET+OF
CUSTOM_R1='ACACTCTTTCCGATCTGGAGCTGTCGTTCACTC'
CUSTOM_I5='GAGTGAACGACAGCTCCAGATCGGAAAGAGTGT'
I7_SP=INDEX1_PRIMER
P5_ORDERING=P5+'B'*8+CUSTOM_R1
P7_ORDERING=P7+'B'*8+TRUSEQ_READ2
AUTHORITY='PublishedPMC5737328Methods76–112; publishedSI AdditionalFig4p6/TableS1p7 hRPS19_2 andTableS3p11. Two originals only.'

def left_scene():
    s=Scene();s.strand('diribo',[Segment('OR handle',OR),Segment('target DNA',LEFT_TARGET_DNA),Segment('terminal RNA2','UA')],mod3='OH; terminal two residues RNA, not DNA');return s

def right_scene():
    s=Scene();s.strand('phospho',[Segment('target DNA',RIGHT_TARGET),Segment('OF same-sense tail',OF)],mod5='P',mod3='OH');return s

def source_primer(seq):
    s=Scene()
    if 'B'*8 in seq:
        prefix,suffix=seq.split('B'*8);which='i5'if seq==P5_ORDERING else'i7'
        segs=[Segment('source 5prime tail',prefix),Segment('source index window','N'*8,placeholder=True,feature=feature('lish.'+which,'sample_index','unknown',group='lish.sample'),note='Source B8 means index nucleotides; actual study value unassigned, not IUPAC B'),Segment('source annealing tail',suffix)]
    else:segs=[Segment('source whole primer',seq)]
    s.strand('primer',segs);return s

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('whole source products','[unassembled source products]',placeholder=True)]);return s

def reporter_class():
    # Expected sealed subclass, not the complete ligation mixture or a DNA sequence.
    return dict(source_orderings=(LEFT_ORDERING,RIGHT_ORDERING),five_prime_to_three_prime_parts=(OR,LEFT_TARGET_DNA,'rUrA',RIGHT_TARGET,OF),junction='left terminal RNA-A 3OH to right DNA-G 5P',internal_RNA_residues=2,ends=('unassigned left5prime','right3OH'),scope='Source intended hRPS19_2 ligated class; unreacted/off-target probes remain until actual separation')

def protocol_graph():
    g=ProtocolGraph('LISH published Plus-slide Illumina branch',authority=AUTHORITY,representation='published_instance')
    x=g.input('tissue','Source10umFFPE tissue section onPlus slide; >=48h formalin before embedding, actual RNA population unassigned',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return(g.input(k,label,kind,**kw),role)
    def op(k,label,kind,inputs=(),product='Source mixture with incomplete/off-target/unreacted classes; occupancy/yield/carryover unassigned',**kw):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR)
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def wash(k,label,inputs=()):
        nonlocal x
        y=g.product(k+'_kept','Source retained physical fraction, not assumed pure desired molecules',K.MOLECULAR);w=g.product(k+'_waste','Actual removed wash/workup fraction',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),inp(k+'_buffer',label))+tuple(inputs),(y,w));g.discard(w,'Physical separation only; exact free/retained classes and recovery unassigned');x=y
    def split(k,label,selected,remainder):
        nonlocal x
        y=g.product(k+'_selected',selected,K.MOLECULAR);r=g.product(k+'_remainder',remainder,K.MOLECULAR)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));x=y
    def pcr(k,label,forward,reverse,reason):
        bulk=g.input(k+'_stock','Actual prepared sourceHerculase-II PCR stock; internal quantities unassembled',K.MOLECULAR)
        enzyme=g.product(k+'_enzyme','Herculase-II catalyst from same stock',K.CATALYST);nts=g.product(k+'_nt','Nucleotide component of same stock; identities/doses unassembled',K.MOLECULAR)
        f=g.product(k+'_f','Source '+forward,K.OLIGO,state=source_primer(forward),chemical_state={'source_ordering':forward});r=g.product(k+'_r','Source '+reverse,K.OLIGO,state=source_primer(reverse),chemical_state={'source_ordering':reverse})
        buf=g.product(k+'_buffer','Condition component of same PCR stock',K.CONDITION)
        g.react(k+'_components','Same source prepared stock components; no duplicate independent additions',S.GENERIC,((bulk,R.TEMPLATE),),(enzyme,nts,f,r,buf),changes_molecule=False)
        op(k,label,S.PCR,((enzyme,R.CATALYST),(nts,R.STRUCTURAL),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(buf,R.CONDITION)),geometry=PCRGeometry.unavailable(unknown_scene(),reason=reason),product='Source amplified DNA reporters plus residual chimeric reporters, primers and nucleotides; whole population unassembled')
    op('bake','Source60C1h bake',S.THERMAL_TREATMENT,(inp('bake_conditions','60C1h'),))
    wash('deparaffinize','SourceTrilogy95C2x30min,4water rinses,3x10minPBS/.1percentTriton,2x5min2xSSC/20percentformamide/.1percentTween')
    left=g.input('diribo_probe','ActualhRPS19_2 3primeDiribo '+LEFT_ORDERING,K.OLIGO,state=left_scene(),chemical_state={'source_ordering':LEFT_ORDERING,'sugar_map':'DNA until terminal two RNA residues','known_end':'3OH'})
    right=g.input('phospho_probe','ActualhRPS19_2 5primePhospho '+RIGHT_ORDERING,K.OLIGO,state=right_scene(),chemical_state={'source_ordering':RIGHT_ORDERING,'sugar':'DNA','five_prime':'P','three_prime':'OH'})
    op('hybridize','Source250ulhybridizationmix5pMeachprobe60C30min then45C2h',S.HYBRIDIZATION,((left,R.STRUCTURAL),(right,R.STRUCTURAL),inp('other_probe_panel','Actualother members ofequimolarmultiplexpanel; completewhole inventory unassembled',K.OLIGO,R.STRUCTURAL),inp('hyb_carrier','Source1mgmlE.colitRNA carrier inhybridizationmix',K.MOLECULAR,R.STRUCTURAL),inp('hyb_conditions','Source2xSSC/20percentformamide/.2mgmlBSA/2mMvanadyl/.1gml dextransulfate')))
    wash('posthyb','Source45C2xSSC/2mMvanadyl/1mgmltRNA,then37CRnl2buffer', (inp('wash_carrier','Actual1mgmlE.colitRNA suppliedinposthybridizationwash',K.MOLECULAR,R.STRUCTURAL),inp('rnl_buffer_atp','Actual1mMATPcomponentofRnl2bufferwash',K.MOLECULAR,R.STRUCTURAL)))
    op('rnl2','Source200ulRnl2reaction1.5Uul,37C2h; RNA-templated nick sealing',S.MAGIC,(inp('rnl2_enzyme','SourceQiagenRnl2',K.CATALYST,R.CATALYST),inp('ligation_atp','Source1mMATP inactual freshligationreactionmix',K.MOLECULAR,R.STRUCTURAL),inp('ligation_conditions','Source50mMTris/10mMMg/5mMDTT')),
       reason='Source RNA-templated ligation positive; full occupied template/probe/panel geometry unassembled, not unavailable sequence claim.',product='LISH chimeric reporter classes retain two internal RNA residues; free probes/misligation/RNA/tRNA remain. Expected selected sealed class specified separately.')
    wash('postligation','Sourceaspirate200ulreactionmix thenroomtemp2xSSCwash')
    op('rnaseh','Source6URNaseH in150ul50mMTrispH8.3/75mMKCl/3mMMg/fresh10mMDTT37C30min',S.MAGIC,(inp('rnase_h','SourceInvitrogenRNaseH; exact isoform/activity distribution unassigned',K.CATALYST,R.CATALYST),inp('rnase_conditions','SourceRNaseH conditions')),
       reason='Source release by RNA hydrolysis; complete target RNA cleavage, original reporter preservation and all released products unassembled. No assumption of perfect RNaseH specificity.',product='Source release mixture: recovered reporter classes and RNA fragments/residual probes; no forced conversion of internal reporter RNA2 to DNA')
    split('aspirate','SourceaspiratesRNaseHreleasemix','Actualaspirated150ulreleasefraction; recoveryunassigned','Actualretainedslide/tissuefraction; laterfateunassigned')
    split('clarify','Source20,000g15min4C; collectclearedsupernatant','Sourceclearedsupernatant','Actualpellet/insolublefraction; laterfateunassigned')
    wash('concentrate','SourceOligoClean&Concentratorcolumn purification; manufacturerinternalprotocolunassembled')
    pcr('preamp','Source20ulHerculaseII95C2min;20cycles95C20s/59C30s/72C30s;72C3min; routinePlusbranch',OF,OR,
        'TableS3OF is identical to right-probe3prime tail, whereas an opposed primer would require its reverse complement. Source printed OF retained; exact effective production primer/geometry unresolved. DNA copies of RNA2 must be distinguished from original chimeric templates.')
    split('index_input','Source2ulpreampselectedinto20ulindexPCR','Actual2ulpreampaliquot','Remainingnominal18ulpreampstock; exactfateunassigned')
    pcr('index','Source20uldual8merindexPCR10cycles',P5_ORDERING,P7_ORDERING,
        'Whole effective preamplification substrate unresolved. Source indexed primer orderings retained; B8 is unassigned sample index, not IUPAC degeneracy.')
    op('indexed_pool','Sourcepools indexedPCR products; actualotherpreparations/membership/quantitiesunassembled',S.MAGIC,
       (inp('other_indexed_products','ActualotherindexedPCRproductcontributions; upstream separatepreparationsunassembled',K.MOLECULAR,R.STRUCTURAL),),reason='Positive source physical pool; complete panel/preparation aliquots and quantities unassembled, not complete single-input pool.')
    wash('qiaquick1','SourcepooledproductsQIAquickPCRpurification firstcolumn')
    wash('qiaquick2','SourcepooledproductsQIAquickPCRpurification secondcolumn')
    split('load_aliquot','Actualsequencerloadingallocation; volumeunassigned','SourceNextSeq500loadingaliquot','Untransferredpooledlibrarystock; fateunassigned')
    lib=g.product('library','SelectedPlus-slideLISHIlluminapreparation',K.LIBRARY);g.react('load','SourceNextSeq500single-end50cycles/customR1/customi5',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourceNextSeq500SE50/dual8mer; customR1 '+CUSTOM_R1+' andcustomi5 '+CUSTOM_I5+'. Full primerlanding/runmask/preampgeometry unresolved. Source40nttarget-derivedreporters/trim3prime10 readbases; mismatchedligation~5percent filtering is computational, not physical purification. HiSeq2500/LCM/ONTalternatives not merged.'))
    g.add_source_item(SourceItem('discrepancy','PublishedTableS1hRPS19_2 righttail equals printedTableS3OF; opposed preamp orientation unresolved. No silently corrected source primer.',AUTHORITY));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();rc=reporter_class()
    c('ActualreporterRNAratherthanDNA',rc['internal_RNA_residues'],2)
    c('Primerorientationissuepreserved',RIGHT_ORDERING.endswith(OF),True)
    c('Noimplicitreplacementwithreversecomplement',OF!=revcomp(OF),True)
    c('Indexreverseprimerdoesendinopposedtail',P7_ORDERING.endswith(revcomp(OF)),True)
    c('Indexforwardprimerendmatcheslefthead',P5_ORDERING.endswith(OR),True)
    c('IndexPCRusesactual2ulaliquot',g.steps['index'].inputs[0].material,'index_input_selected')
    c('Preampusescolumnfraction',g.steps['preamp'].inputs[0].material,'concentrate_kept')
    c('RNaseHreleasedfractionnotwholefixedtissue',g.steps['clarify'].inputs[0].material,'aspirate_selected')
    c('Twoactualpooledcolumnpurifications',g.steps['qiaquick2'].inputs[0].material,'qiaquick1_kept')
    c('SourceB8retainedliterally',P5_ORDERING.count('B'),8)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkx471_supp_gkx471_supp.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1093/nar/gkx471')):
        t=''.join(p.read_text().split())
        for n,s in [('OF',OF),('OR',OR),('CustomR1',CUSTOM_R1),('Customi5',CUSTOM_I5),('P5index',P5_ORDERING),('P7index',P7_ORDERING)]:c('Sourcewholeordering '+n,s in t,True)
    c('Partialnotpublicready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
