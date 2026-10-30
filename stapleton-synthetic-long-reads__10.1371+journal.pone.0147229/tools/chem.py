"""Source-local synthetic-long-read routes; variant identity and molecular gaps explicit."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Scene,Construct,Segment,feature,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry
from page import head,info
PAGE_FILE='chem.html'
TITLE='Stapleton synthetic long reads — source-specific barcode pairing'
SOURCE='Published PMC4720449 DOI10.1371/journal.pone.0147229 Methods183–212,S13Table and S1Note;no BAsE/Loop kit recipe inheritance'
TABLE='journal.pone.0147229_supp_pone.0147229.s030.docx.txt'
LIMIT='Private partial model. Complete dU/USER-cut/barcode/circle/biotin-bound products,end geometry and final sequencing-primer assembly remain unassembled. Operational source continuity is distinct from validated full molecular chemistry. Single-tube and two-tube routes are separate source variants.'

@dataclass(frozen=True)
class PhysicalOligo:
    name:str
    sequence:str
    source_notation:str
    phosphate5:bool=False
    deoxy_u:tuple[int,...]=()
    ribo:tuple[int,...]=()
    lna:tuple[int,...]=()
    ps_after:tuple[int,...]=()
    def __post_init__(self):
        for i in self.deoxy_u:
            if not 0<=i<len(self.sequence) or self.sequence[i]!='U':raise ValueError('deoxyU needs its physical uracil base')
        if any(not 0<=i<len(self.sequence) for i in self.ribo+self.lna):raise ValueError('Sugar annotation outside physical primer')
        if any(not 0<=i<len(self.sequence)-1 for i in self.ps_after):raise ValueError('PS linkage must join two existing bases')
    def as_scene(self):
        if self.name.startswith('Oligo12') or self.name.startswith('Oligo13'):
            a=self.sequence.index('N'*16)
            tail=self.sequence[a+16:];anchored=tail.endswith('VN')
            parts=[Segment(self.name+' 5primebody',self.sequence[:a]),Segment(self.name+' molecule barcode',self.sequence[a:a+16],placeholder=True,feature=feature('stapleton.rt_barcode' if self.name.startswith('Oligo12') else 'stapleton.tso_barcode','linked_read_barcode','random')),Segment(self.name+' 3primebody',tail[:-2] if anchored else tail)]
            if anchored:parts.append(Segment('source anchoredVN','VN',placeholder=True))
        else:parts=[Segment(self.name,self.sequence)]
        sc=Scene();sc.strand('oligo',parts,mod5='phosphate' if self.phosphate5 else '')
        sc.note('oligo','Source physical ordering: '+self.source_notation+';sugar/linkage annotations apply to this reagent,not polymerase copies.')
        sc.source_physical_oligo=self
        return sc
    def rows(self):return self.as_scene().rows()

@dataclass(frozen=True)
class SourceAdapterPreparation:
    oligos:tuple[PhysicalOligo,...]
    def as_scene(self):
        sc=Scene()
        for i,o in enumerate(self.oligos):
            sc.strand('oligo'+str(i),[Segment(o.name,o.sequence)],mod5='phosphate' if o.phosphate5 else '',label='Source ordering component;annealed geometry unassembled')
            sc.note('oligo'+str(i),'Physical source ordering: '+o.source_notation)
        sc.source_adapter_preparation=self
        return sc
    def rows(self):return self.as_scene().rows()

@dataclass(frozen=True)
class IndexedPrimer:
    construct:Construct
    ps_after:tuple[int,...]=()
    def as_scene(self):
        sc=Scene();sc.strand('primer',self.construct.segments,label=self.construct.name)
        if self.ps_after:sc.note('primer','Source physical PS bond immediately before terminal T;not assigned to polymerase-made copies.')
        sc.source_index_primer=self
        return sc
    def rows(self):return self.as_scene().rows()

PRIMER2_SEQ='CCATTCATUAATGCACAUAACTATTCCUG'
PRIMER2=PhysicalOligo('Oligo2 selected unambiguous internal-dU variant',PRIMER2_SEQ,'CCATTCAT/ideoxyU/AATGCACA/ideoxyU/AACTATTCC/ideoxyU/G',deoxy_u=tuple(i for i,b in enumerate(PRIMER2_SEQ) if b=='U'))
OL3=PhysicalOligo('Oligo3','CCGAGAAT TCCAT'.replace(' ',''),'CCGAGAATTCCA*T',ps_after=(11,))
OL4=PhysicalOligo('Oligo4','TGGAATTCTCGGGTGCCAAGG','/5Phos/TGGAATTCTCGGGTGCCAAGG',True)
OL6_SEQ=il.P5+il.TRUSEQ_READ1[4:]
OL6=PhysicalOligo('Oligo6',OL6_SEQ,OL6_SEQ[:-1]+'*'+OL6_SEQ[-1],ps_after=(len(OL6_SEQ)-2,))
OL7=PhysicalOligo('Oligo7','ACACTCTTTCCCTACACGACGCTCTTCC','ACACTCTTTCCCTACACGACGCTCTTCC')
OL8=PhysicalOligo('Oligo8','ATCGGAAGAGCACACGTCT','/5Phos/A*TCGGAAGAGCACACGTCT',True,ps_after=(0,))
OL10=PhysicalOligo('Oligo10',il.P5+'GTTCAGAGTTCTACAGTCCGA',il.P5+'GTTCAGAGTTCTACAGTCCGA')
RT12_SEQ='CCATTCATTAATGCACATAACTATTCCTGG'+'N'*16+'GATCGTCGGACTGTAGAACTCTGAAC'+'T'*30+'VN'
RT12=PhysicalOligo('Oligo12 barcoded oligo-dT',RT12_SEQ,'CCATTCATTAATGCACATAACTATTCCTGG'+'N'*16+'GATCGTCGGACTGTAGAACTCTGAACT30VN')
TSO13_SEQ='GCGGCCATTCATTAATGCACATAACTATTCCTGT'+'N'*16+'AGATCGGAAGAGCGTCGTGTAGG'+'TGGG'
TSO13=PhysicalOligo('Oligo13 mixedDNA/rG/rG/LNAG',TSO13_SEQ,'GCGGCCATTCATTAATGCACATAACTATTCCTGT'+'N'*16+'AGATCGGAAGAGCGTCGTGTAGGTrGrG+G',ribo=(len(TSO13_SEQ)-3,len(TSO13_SEQ)-2),lna=(len(TSO13_SEQ)-1,))

def barcode_adapter(alternate=False):
    return Construct([Segment('unassigned random3 window','N'*3,note='Source N3;not silently assigned a new identifier function'),
        Segment('source sequencing region','GTTCAGAGTTCTACAGTCCGACGATC' if alternate else 'CCTACACGACGCTCTTCCGATCT'),
        Segment('molecule barcode','N'*16,placeholder=True,feature=feature('stapleton.molecule','linked_read_barcode','random')),
        Segment('source adjacent spacer','CC' if alternate else 'AC'),Segment('PCR priming region','AGGAATAGTTATGTGCATTAATGAATGG'),Segment('source terminal sequence','CCGC')],
        'Oligo1 selected source adapter variant;5primephosphate')

def index_primer(pair=False):
    tail=il.TRUSEQ_READ2 if pair else 'GTGACTGGAGTTC'+'CTTGGCACCCGAGAATTCCA'
    return Construct([Segment('P7',il.P7),Segment('source sixbase index','N'*6,'idx',placeholder=True,
        feature=feature('stapleton.pair_index' if pair else 'stapleton.main_index','sample_index','unknown',note='Two-tube branch must choose a different index from main.')),
        Segment('source PCR body',tail)],'Oligo9 pairing index primer' if pair else 'Oligo5 main index primer')

def unknown_scene():
    sc=Scene();sc.strand('unknown',[Segment('source product unassembled','[unassembled molecular product/end geometry]',placeholder=True)])
    sc.note('unknown','This is an assembly boundary,not an accepted source-disclosure absence.');return sc

def variant_graph(variant):
    if variant not in ('dna_two_tube','rna_two_tube','dna_single_tube'):raise ValueError('Choose an explicitly source-scoped variant')
    g=ProtocolGraph(TITLE+' / '+variant,authority=SOURCE)
    def inp(key,label,kind=K.CATALYST,state=None):return g.input(key,label,kind,state=state.as_scene() if hasattr(state,'as_scene') else state)
    def op(key,x,label,kind,reagents=(),geometry=None):
        y=g.product(key+'_product',label+';reaction products persist to workup',K.MOLECULAR);kw={} if geometry is None else {'geometry':geometry}
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(reagents),(y,),**kw);return y
    def clean(key,x,label):
        sel=inp(key+'_selection',label,K.SELECTION);y=g.product(key+'_kept','Retained source DNA fraction',K.MOLECULAR);w=g.product(key+'_waste','Removed physical workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(y,w));g.discard(w,'Physical source separation only;no perfect recovery or specificity assumed.');return y
    def pcr(key,x,primers,label,symmetric=False):
        enzyme=inp(key+'_polymerase','Source LongAmpTaq40cycles for initial amplification' if key=='long_pcr' else 'Source PhusionHF18cyclelcPCR')
        return op(key,x,label,S.PCR,tuple((p,R.STRUCTURAL) for p in primers)+((enzyme,R.CATALYST),),PCRGeometry.unavailable(unknown_scene(),reason='Complete source strands/primer sites unassembled;physical primer inventory is source-specific.',symmetric_species=symmetric))
    def ligate(key,x,adapter,label):
        enzyme=inp(key+'_ligase','Source NEB TA/Blunt ligase master mix')
        return op(key,x,label,S.LIGATION,((adapter,R.STRUCTURAL),(enzyme,R.CATALYST)),AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Source physical oligos known;full donor/substrate end/ligation junction geometry not assembled.'))
    ol2=inp('oligo2','Physical dU-containing Oligo2 selected variant',K.OLIGO,PRIMER2)
    if variant.startswith('rna'):
        rna=g.input('rna','Source mRNA;upstream extraction outside graph',K.MOLECULAR,primary=True)
        rt=inp('rt12','Source barcoded oligo-dT Oligo12',K.OLIGO,RT12);tso=inp('tso13','Source mixed-sugar TSO Oligo13',K.OLIGO,TSO13)
        enzyme=inp('rt_reagents','Source-cited full-length RT procedure[36];exact predecessor version/reagents not inspected')
        x=op('rna_front',rna,'Source modified full-length RT/TSO procedure;products enter library quantification',S.REVERSE_TRANSCRIPTION,((rt,R.STRUCTURAL),(tso,R.STRUCTURAL),(enzyme,R.CATALYST)))
    else:
        dna=g.input('dna','Source genomic DNA',K.MOLECULAR,primary=True)
        shear=inp('gtube','Source g-TUBE about10kb');x=op('input_shear',dna,'Fragment genomic DNA',S.FRAGMENTATION,((shear,R.CATALYST),))
        repair=inp('input_end_repair','NEBNext End Repair Module');x=op('input_repair',x,'End repair',S.END_REPAIR,((repair,R.CATALYST),));x=clean('input_column',x,'Zymo column purification')
        dt=inp('dt_tail_reagents','Klenowexo-/dTTP/PNK37C1h');x=op('dt_tail',x,'Add source dT tails',S.TAILING,((dt,R.CATALYST),))
        adapter_state=Scene();adapter_state.strand('adapter',barcode_adapter().segments,mod5='phosphate')
        if variant.endswith('single_tube'):
            adapter_state.strand('alternate_adapter',barcode_adapter(True).segments,label='Second source sequencing-region variant for single-tube barcode pairing',mod5='phosphate')
        # Reagent formation is not a second substrate genome or copied dU chemistry.
        source1=inp('oligo1','Source5primephosphorylated barcode adapter;selected sequence inventory',K.OLIGO,state=adapter_state)
        assembly=inp('adapter_extension','Klenowexo-/dNTPs37C60min after95C10min anneal')
        made=g.product('assembled_adapter','Source extended barcode-adapter preparation,complete duplex unassembled',K.ADAPTER)
        g.react('make_adapter','Anneal and extend source barcode oligos',S.EXTENSION,((source1,R.TEMPLATE),(ol2,R.STRUCTURAL),(assembly,R.CATALYST)),(made,))
        x=ligate('barcode_ligation',x,made,'Ligate source extended adapters to dT-tailed DNA')
        x=clean('barcoded_recovery',x,'Source optionalgel purification if necessary;elute inEB')
    # Both source front ends explicitly enter quantification/amplification here.
    chosen=g.product('chosen','Approximately100000 doubly barcoded molecule input',K.MOLECULAR);unused=g.product('untransferred','Untransferred quantified preparation',K.MOLECULAR)
    g.react('molecule_bottleneck','SourceqPCR/dilution-series estimate then dilute/select input',S.SPLIT,((x,R.TEMPLATE),),(chosen,unused));g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Untransferred material not presumed disposed;no destructiveQC aliquot is returned to library.'))
    x=pcr('long_pcr',chosen,(ol2,),'Single-primer LongAmpTaq PCR across eight tubes',True);x=clean('long_pcr_gel',x,'Qiagen gel purification')
    user=inp('user','USER enzyme37C30min');x=op('user_digest',x,'Remove dU-containing primer regions with USER;cleavage products remain',S.DIGESTION,((user,R.CATALYST),))
    frag=inp('fragmentase','dsDNAfragmentase/BSA37C15min');x=op('random_break',x,'Source aboutone break on average;not exactlyone break per molecule',S.FRAGMENTATION,((frag,R.CATALYST),))
    stop=inp('edta','0.5MEDTA quench',K.CONDITION);x=op('stop_fragmentase',x,'Quench fragmentase',S.QUENCHING,((stop,R.CONDITION),));x=clean('fragment_cleanup',x,'0.8×AMPure cleanup')
    ecoli=inp('ecoli_ligase','E.coliDNAligase20C20min');x=op('nick_repair',x,'Source E.coliDNA ligase treatment',S.NICK_SEALING,((ecoli,R.CATALYST),))
    fill=inp('biotin_fill_reagents','T4DNApolymerase+Klenow+50uMbiotin-dCTP20C10min');x=op('biotin_fill',x,'Initial biotin-dCTP-only treatment',S.END_REPAIR,((fill,R.CATALYST),))
    dntp=inp('remaining_dntps','Then50uMdGTP,dTTP,dATP20C15min',K.CONDITION);x=op('complete_fill',x,'Add remaining three nucleotides;prior T4/Klenow enzymes retained,no fresh enzyme addition',S.END_REPAIR,((dntp,R.CONDITION),(fill,R.CATALYST)));x=clean('fill_cleanup',x,'1×AMPure cleanup')
    ligase=inp('circle_ligase','T4DNAligase,1ng/uL,roomtemperature16–48h')
    x=op('circularize',x,'Source intramolecular dsDNA circularization;heterogeneous products retained',S.CIRCULARIZATION,((ligase,R.CATALYST),),CircularizationGeometry.known('double-strand intramolecular ligation',statement='Published mainMethods specifies T4DNA ligase on preceding end-treated DNA;exact heterogeneous closure coordinates/unreacted molecules remain unassembled.'))
    t5=inp('t5','T5exonuclease37C60min');x=op('t5_digest',x,'Digest accessible linear DNA;no complete purity assumed',S.DIGESTION,((t5,R.CATALYST),));x=clean('circle_column',x,'Zymo column purification')
    covaris=inp('covaris','CovarisS2 to approximately500bp');x=op('circle_shear',x,'Fragment circular preparation',S.FRAGMENTATION,((covaris,R.CATALYST),))
    beads=inp('streptavidin','DynabeadsM280streptavidin20C15min',K.SELECTION)
    bound=g.product('bound','Retained biotin-DNA/bead complex',K.MOLECULAR);waste=g.product('unbound','UnboundDNA and source wash fractions',K.WASTE)
    g.react('capture','Source streptavidin capture and washes',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,waste));g.discard(waste,'Physical magnet/wash separation only;capture efficiency not predicted.')
    if variant.endswith('two_tube'):
        main=g.product('main_beads','Remaining approximately85% bead-bound fraction',K.MOLECULAR);pair=g.product('pair_beads','Separate approximately15% bead-bound fraction',K.MOLECULAR)
        g.react('pair_split','Move15% beads to separate pairing tube',S.SPLIT,((bound,R.TEMPLATE),),(main,pair))
    else:main=bound;pair=None
    repair=inp('bead_repair_reagents','NEBNextEndRepair20C30min');main=op('bead_repair',main,'Bead-bound end repair',S.END_REPAIR,((repair,R.CATALYST),));main=clean('bead_repair_wash',main,'Magnet-retain beadDNA and wash')
    atail=inp('a_tail_reagents','Klenowexo-/dATP37C30min');main=op('a_tail',main,'Bead-bound dA tailing',S.TAILING,((atail,R.CATALYST),));main=clean('a_tail_wash',main,'Magnet-retain beadDNA and wash')
    adapter34=inp('adapter34','Source annealed Oligo3+4 lcPCR adapter;exact annealed end geometry unassembled',K.ADAPTER,state=SourceAdapterPreparation((OL3,OL4)))
    main=ligate('main_ligation',main,adapter34,'Ligate Oligo3+4,roomtemperature15min');main=clean('main_ligation_wash',main,'Magnet-retain beadDNA and wash')
    five=inp('oligo5','Source mainPCRindex primer;index6 not assigned',K.OLIGO,state=index_primer())
    six=inp('oligo6','Source Oligo6 with terminalPS bond',K.OLIGO,state=OL6)
    primers=(five,six)
    if variant.endswith('single_tube'):
        ten=inp('oligo10','Additional source single-tube pairing primerOligo10',K.OLIGO,state=OL10);primers+=(ten,)
    main=pcr('main_pcr',main,primers,'Source18cyclelcPCR;single-tube addsOligo10 only inits own variant')
    main=clean('main_release',main,'Retain PCR supernatant;discard beads')
    main=clean('main_size',main,'Source0.7×AMPure OR500–900bpgel recovery;selection alternative not merged')
    if pair is not None:
        sap=inp('sapi','SapI37C1h');pair=op('pair_sapi',pair,'Separate pairing tube SapI digest',S.DIGESTION,((sap,R.CATALYST),));pair=clean('pair_sapi_wash',pair,'Magnet-retain beadDNA and wash')
        adapter78=inp('adapter78','Source annealed Oligo7+8 Yadapter;exact annealed end geometry unassembled',K.ADAPTER,state=SourceAdapterPreparation((OL7,OL8)))
        pair=ligate('pair_ligation',pair,adapter78,'Pairingadapterligation4C15min then20C15min');pair=clean('pair_ligation_wash',pair,'Magnet-retain beadDNA and wash')
        nine=inp('oligo9','Pairing index primer uses a DIFFERENT source6base index from Oligo5;terminalPS bond',K.OLIGO,state=IndexedPrimer(index_primer(True),(len(index_primer(True))-2,)))
        pair=pcr('pair_pcr',pair,(six,nine),'Pairing18cyclelcPCR')
        pair=clean('pair_release',pair,'Retain pairingPCR supernatant;discard beads');pair=clean('pair_size',pair,'1.8×AMPure pairing-library cleanup')
        mixed=g.product('run_pool','Physical sequencing pool with1–5% pairing-library molecules',K.MOLECULAR)
        g.react('run_pooling','Mix pairing and main libraries at source1–5% pairing abundance',S.POOL,((main,R.TEMPLATE),(pair,R.TEMPLATE)),(mixed,));main=mixed
    g.sequence(main,SequencingPlan(acquisition='unspecified',limitation=LIMIT+' Source MiSeq/NextSeq/HiSeq standardprimer mixtures;exact final physical priming/read-plan geometry unresolved. E.coliMG1655Results150×150 example belongs to single-tubevariant;not silently assigned to generic two-tuberoute. Computational barcode-pairing does not physically join original genomes.'))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE);g.add_source_item(SourceItem('unassembled_source','Remaining source/variant and exact polymer/read assembly',LIMIT));return g

def protocol_graph():return variant_graph('dna_two_tube')

def render_page():
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',info(LIMIT),'<p><a href="01_stapleton.html">Research notes</a></p>',
        ''.join('<h2>'+v+'</h2>'+workflow_panel(workflow_from_protocol_graph(variant_graph(v),'run_pool' if v.endswith('two_tube') else 'main_size_kept')) for v in ('dna_two_tube','rna_two_tube','dna_single_tube')),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph();single=variant_graph('dna_single_tube');rna=variant_graph('rna_two_tube')
    c('Allvariants private',any(x.publishable for x in (g,single,rna)),False)
    c('No false completion marker','data-protocol-complete=' in render_page(),False)
    c('Initialsingle-primerPCR inputcount',len([x for x in g.steps['long_pcr'].inputs if x.role==R.STRUCTURAL]),1)
    c('USER precedes fragmentation',g.steps['random_break'].inputs[0].material,'user_digest_product')
    c('OtherdNTPs added after initial biotin-dCTP treatment',g.steps['complete_fill'].inputs[0].material,'biotin_fill_product')
    c('Pairing aliquot before main repair',g.steps['bead_repair'].inputs[0].material,'main_beads')
    c('Single-tube has no invented15% split','pair_split' in single.steps,False)
    c('Single-tube retains both source adapter sequence variants',len(single.materials['oligo1'].state.strands),2)
    c('RNA enters explicit common quantified material',rna.steps['molecule_bottleneck'].inputs[0].material,'rna_front_product')
    tso_state=rna.materials['tso13'].state.source_physical_oligo
    c('RNA TSO retains twoRNA and oneLNA residues',(len(tso_state.ribo),len(tso_state.lna)),(2,1))
    p=sources_root()/Path(__file__).resolve().parents[1].name/TABLE
    if p.exists():
        source=''.join(p.read_text().split()).upper()
        for o in (PRIMER2,OL3,OL4,OL6,OL7,OL8,OL10,RT12,TSO13):c('PublishedS13physicaloligo '+o.name,o.source_notation.upper() in source,True)
        import re
        definition=re.search(r'INDEX=([1-9][0-9]*)-BASEILLUMINATRUSEQSMALLRNAMULTIPLEXINGINDEXSEQUENCE',source)
        width=int(definition.group(1)) if definition else None
        for pair in (False,True):
            model_width=len(next(s.top for s in index_primer(pair).segments if s.feature is not None))
            c('Actual model index width matches printed S13 definition '+str(pair),model_width,width)
    else:c.skip('S13transcriptions','External original unavailable')
    c.report()
