from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        g.react(key,label,kind,inputs,(y,));return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g

from dataclasses import dataclass
import illumina as il
import nextera as nx
from protocol_graph import Terminal,Disposition
TITLE='NEAT-seq — published protocol v1.0 nuclear protein, ATAC and RNA recovery'
SOURCE='PublishedSI MOESM1 NEAT-seqprotocolv1.0 andOligosequencespp4–6/10;publishedMOESM3TableS4;acceptedPMC11192021Methods separatelyidentified'
LIMIT='Private partial model. Published protocolv1.0 uses22uL preamp product,not35uL fromacceptedmanuscript. PublishedSIHTONexterahandle andantibody-type-specificrecipesnotmergedwithacceptedbarnyardTruseqHTOs. Whole historicalMultiomeRevA products/readgeometry andfinalpublishedmain applicabilityremainunassembled.'
ADT_HEAD=il.TRUSEQ_READ2[12:]
HTO_HEAD=nx.ADAPTOR_S7[8:]
PRE_ADT=il.TRUSEQ_READ2[12:-1]
PRE_HTO=HTO_HEAD
@dataclass(frozen=True)
class ModifiedDNA:
    parts:tuple
    source_notation:str
    ps_links:tuple=()
    terminal3:str=''
    terminal_base:str=''
    def scene(self):
        sc=Scene();sc.strand('oligo',list(self.parts),label=self.source_notation,mod3=self.terminal3)
        sc.note('oligo','Physical sourceDNA: '+self.source_notation+';phosphorothioatelink positions='+str(self.ps_links)+'. Endpointmodificationsbelongto suppliedpolymer,not automaticallytoPCRcopies.')
        sc.physical_dna=self;return sc

@dataclass(frozen=True)
class ProteinAssociatedDNA:
    dna:ModifiedDNA
    protein:str='EcoSSB'
    def scene(self):
        sc=self.dna.scene();sc.note('oligo','Source qualitative association with '+self.protein+';unbound protein and buffer remain,no stoichiometry or perfect occupancy asserted.')
        sc.protein_association=self;return sc

def blocking_oligo():
    return ModifiedDNA((Segment('source randomblock extent','[random DNA extent;source calls30mer,ddC count convention unresolved]',placeholder=True),Segment('terminal dideoxycytidine nucleotide','C')),'source30merDNA block withterminalddC',terminal3='ddC sugar;no3primeOH',terminal_base='C')

def tag(which,member=''):
    if which not in ('adt','hto'):raise ValueError(which)
    head=ADT_HEAD if which=='adt' else HTO_HEAD
    role='feature_barcode' if which=='adt' else 'sample_index'
    return ModifiedDNA((Segment('source handle',head),Segment('source reporter UMI','N'*12,placeholder=True,feature=feature('neat.'+which+'_umi','umi','random')),Segment('source antibody orhashidentity','N'*36,placeholder=True,feature=feature('neat.'+which+'_identity'+member,role,'unknown')),Segment('polyA','A'*25)),head+'[12bp UMI][36bp Barcode]'+'A'*23+'*A*A',(-2,-1))

def index_primer(which):
    if which not in ('p5','adt_p7','hto_p7'):raise ValueError(which)
    if which=='p5':prefix=il.P5;arm='ACACTCTTTCCCTACACGACGCTCTTCCGATCT';role='neat.library_i5'
    else:prefix=il.P7;arm=('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT' if which=='adt_p7' else nx.ADAPTOR_S7);role='neat.library_i7'
    return ModifiedDNA((Segment('source flowcellpart',prefix),Segment('8bpindex','N'*8,placeholder=True,feature=feature(role,'sample_index','unknown')),Segment('source PCRarm',arm)),prefix+'[8bp index]'+arm[:-2]+'*'+arm[-2]+'*'+arm[-1],(-2,-1))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    cells=g.input('cells','Source protocol cell preparation;not assigned a particularacceptedmanuscript sample',K.BIOLOGICAL,primary=True)
    x=op('fix',cells,'1.6%PFA2minRT;0.25Mglycine5minonice','Source fixation/quenchconditions')
    x=clean('fix_wash',x,'TwoPBS/3%BSA/0.1%Tweenwashes;retainpellet')
    x=op('perm',x,'0.5%NP40isotonic lysis/permeabilization10minonice','PublishedSITris/NaCl/Mg/NP40/Tween/digitonin/RNaseinhibitor/DTT')
    x=clean('perm_wash',x,'Source800gpelleting andwashbuffer/Tweenwashes')
    block=g.input('blocker','Published source randomDNAblock withactualterminalddC nucleotide',K.OLIGO,state=blocking_oligo().scene())
    blocked=g.product('blocked','Blockednuclei plusfree blocker/buffer components',K.MOLECULAR)
    g.react('block','1mg/mlDNAblock30minRT;boilblock95C2minbeforebufferpreparation',S.GENERIC,((x,R.TEMPLATE),(block,R.STRUCTURAL)),(blocked,))
    x=clean('block_wash',blocked,'Tweenpelleting/twowashes;omitDTTsecondwash;no freshssDNAblockafterthisstep')
    coated={}
    for key in ('adt','hto_1','hto_2','hto_3','hto_4','hto_5'):
        which='adt' if key=='adt' else 'hto'
        member='' if key=='adt' else key[-1]
        dna=g.input(key+'_oligo','Source '+key.upper()+' antibody-linkedDNA;distincthashmember,attachmentposition/chemistryunassembled',K.OLIGO,state=tag(which,member).scene())
        ssb=g.input(key+'_ssb','Source EcoSSB protein;no computed occupancy/resistanceguarantee',K.MOLECULAR)
        c=g.product(key+'_coated','Antibody-linkedDNA qualitativelyEcoSSB-associated;unboundSSB/bufferalso remain',K.MOLECULAR,state=ProteinAssociatedDNA(tag(which,member)).scene())
        g.react(key+'_ssb_bind','Source antibody/EcoSSBpreincubation;addstainingbufferdirectly,NO intermediatepurification',S.GENERIC,((dna,R.TEMPLATE),(ssb,R.STRUCTURAL)),(c,));coated[key]=c
    # Source5tubes carrydistincthasholigos,notfivecopiesofeveryoriginalnucleus.
    tubes=[g.product('tube_'+str(i),'Disjointsource nuclei aliquot '+str(i)+' fordistinctNPChash',K.MOLECULAR) for i in range(1,6)]
    g.react('hash_allocation','Split source nucleipopulationinto5tubes',S.SPLIT,((x,R.TEMPLATE),),tuple(tubes))
    washed=[]
    for i,tube in enumerate(tubes,1):
        y=g.product('hashed_'+str(i),'NPChashed nuclei aliquot '+str(i)+' plusunboundreagents',K.MOLECULAR)
        g.react('hash_'+str(i),'Source distinctantiNPC/HTOstaining30minRT',S.GENERIC,((tube,R.TEMPLATE),(coated['hto_'+str(i)],R.STRUCTURAL)),(y,))
        washed.append(clean('hash_wash_'+str(i),y,'Two stainingbuffer/Tweenwashes;retainhashednuclei'))
    x=g.product('hash_pool','PooledsourceNPChashednuclei;cellallocationsnotduplicated',K.MOLECULAR)
    g.react('pool_hashes','Poolwashedhashednuclei',S.POOL,tuple((y,R.TEMPLATE) for y in washed),(x,))
    y=g.product('tf_stained','TF-stainednuclei plusunboundantibody/reagents',K.MOLECULAR)
    g.react('tf_stain','SourceTFantibodymix30minRT;protocolv1.0genericpanel',S.GENERIC,((x,R.TEMPLATE),(coated['adt'],R.STRUCTURAL)),(y,))
    x=clean('tf_wash',y,'Two stainingbuffer/Tweenwashes;resuspend1xNucleibuffer,40umFlowmifilter')
    x=op('tagment',x,'Source historicalMultiomeTn5transposition;do notimportcurrentkitrecipe','HistoricalMultiomeRevA componentsunreviewed')
    x=op('multiome',x,'HistoricalMultiomepartitioning,RTandrecoverycontinuation;exactintermediatestepsunreviewed','SourceMultiomeRevA reagents')
    preadt=g.input('pre_adt','Published ADT preamp primer',K.OLIGO,state=Construct([Segment('sourceADTpreamp',PRE_ADT)]))
    prehto=g.input('pre_hto','Published HTO preamp primer',K.OLIGO,state=Construct([Segment('sourceHTOpreamp',PRE_HTO)]))
    pre=g.product('preamp_product','Source preamp mixture containinggenomic,cDNA,ADTandHTO-derivedspecies',K.MOLECULAR)
    g.react('preamp','PublishedSIstep4.1add1uL20uMprimer mixtureduringpreamp',S.GENERIC,((x,R.TEMPLATE),(preadt,R.STRUCTURAL),(prehto,R.STRUCTURAL)),(pre,))
    gene=g.product('gene_library_input','Sourcegeneexpressionlibraryallocation;exactvolumeunassigned',K.MOLECULAR)
    atac=g.product('atac_library_input','SourceATACsequencinglibraryallocation;exactvolumeunassigned',K.MOLECULAR)
    adt=g.product('adt_aliquot','22uL preamp aliquotforADTindexPCR',K.MOLECULAR)
    hto=g.product('hto_aliquot','22uL preamp aliquotforHTOindexPCR',K.MOLECULAR)
    remain=g.product('preamp_remaining','Other unallocatedpreampmaterial;notpresumeddiscarded',K.MOLECULAR)
    g.react('library_allocations','Separate manufacturerGEX/ATACcontinuationsand22uLreporterlibraryaliquots;geometryunassembled',S.SPLIT,((pre,R.TEMPLATE),),(gene,atac,adt,hto,remain))
    g.mark_terminal(Terminal(remain.id,Disposition.UNRESOLVED,'Exact totalpreampvolume/remainderfatenotassignedfrominspectedSI.'))
    for key,x in (('gene',gene),('atac',atac)):
        y=op(key+'_library',x,'HistoricalMultiome'+key+'libraryconstruction andsourceQC','Historicalmanufacturerrecipeunreviewed')
        end(g,y,'Sourceprotocolgenerates'+key+'sequencinglibrary;acceptedMSruncontextsconditional,notmergedwithpublishedSIrecipe.Exactprimergeometry/runlengthunassembled.')
    p5=g.input('p5_index','PublishedPS-modifiedP5indexprimer',K.OLIGO,state=index_primer('p5').scene())
    for which,x in (('adt',adt),('hto',hto)):
        p7=g.input(which+'_p7','PublishedPS-modified'+which+'P7indexprimer',K.OLIGO,state=index_primer(which+'_p7').scene())
        y=g.product(which+'_pcr_product','ReporterDNAamplificationmixture;originalmodifiedprimersandunmodifiedcopiesnotconflated',K.MOLECULAR)
        g.react(which+'_index_pcr','Source50uLPCR22uLpreamp/25uLNEBNextHiFi/3uLprimers;98/67/72cycles,source-cyclecountn unassembled',S.GENERIC,((x,R.TEMPLATE),(p5,R.STRUCTURAL),(p7,R.STRUCTURAL)),(y,))
        sel=g.input(which+'_spri06','0.6xSPRI30uL',K.SELECTION)
        sup=g.product(which+'_supernatant','Source0.6xSPRI supernatant includingreporterproducts',K.MOLECULAR)
        beads=g.product(which+'_large_beads','NonretainedlargefragmentSPRIbeadfraction;compositionnotpure',K.WASTE)
        g.react(which+'_size06','FirstSPRIselection;retainsupernatant',S.CLEANUP,((y,R.TEMPLATE),(sel,R.SELECTION)),(sup,beads));g.discard(beads,'Source transferssupernatanttoanother tube.')
        used=g.product(which+'_used75','Source75uLsupernatanttransferred',K.MOLECULAR)
        left=g.product(which+'_leftover','Untransferredsupernatantvolumeunassembled',K.MOLECULAR)
        g.react(which+'_take75','Transfer75uLselectedsupernatant;noassumedquantitativerecovery',S.SPLIT,((sup,R.TEMPLATE),),(used,left));g.mark_terminal(Terminal(left.id,Disposition.UNRESOLVED,'Sourceuntransferredsupernatantnotassignedtolibraryordisposal.'))
        y=clean(which+'_size12',used,'Add28uLSPRI forsource1.2xfinal;80%ethanolwash/elute20uLEB;removebeads')
        y=op(which+'_qc',y,'SourceTapeStation/Bioanalyzerassessment;sequencinglibraryhandoff')
        end(g,y,'PublishedSI'+which.upper()+'reporterlibrary;8bpindexprimerclasses distinctfromreporterUMI12/identity36andkitcellbarcodeunknown. NoMS-only28/48/8/8runassignedtothisgenericrecipe.')
    return incomplete(g)
def render_page():
    g=protocol_graph();return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_neat.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(g,'adt_qc_product'))+'</main>'
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph();block=blocking_oligo()
    c('Block retains actualterminalcytidine',block.terminal_base,'C')
    c('Blockactualstrand hasddCend',g.materials['blocker'].state.strands['oligo'].mod3,'ddC sugar;no3primeOH')
    c('Coatedstate retains sourceprotein association',g.materials['adt_coated'].state.protein_association.protein,'EcoSSB')
    c('No postSSB purification invented',g.steps['tf_stain'].inputs[1].material,'adt_coated')
    c('Hash aliquots receive distinct prepared source antibodies',len({g.steps['hash_'+str(i)].inputs[1].material for i in range(1,6)}),5)
    c('Five disjointsourcehashaliquots',len(g.steps['hash_allocation'].outputs),5)
    c('ReporterPCRsource22uLallocation',g.steps['adt_index_pcr'].inputs[0].material,'adt_aliquot')
    c('Secondselection usestransferredsupernatant',g.steps['adt_size12'].inputs[0].material,'adt_used75')
    c('TransferredSPRIremainder accounted',g.terminals['adt_leftover'].disposition,Disposition.UNRESOLVED)
    c('Privatepartialgraph',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41592-022-01461-y_41592_2022_1461_MOESM1_ESM.pdf.txt'
    if p.exists():
        t=p.read_text();c('PublishedreporterhandleandPSnotation',ADT_HEAD in t and 'AAAAAAAAAAAAAAAAAAAAAAA*A*A' in t,True);c('PublishedterminalddCblockwording','dideoxy' in t and 'cytosine' in t,True)
    else:c.skip('PublishedSItranscription','Externaloriginalabsent')
    c.report()
