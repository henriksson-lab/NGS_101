"""Smart-seq-total v1/v2: artificial RNA tails and typed degradable TSO chemistry."""
from dataclasses import dataclass,replace
from pathlib import Path
import re,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Segment,feature
import nextera as nx
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SourceDiscrepancy,SourceClaim,Highlight
UDG_CATALYSIS='NEB M0280 official current product description,2026-10-29: monofunctional glycosylase hydrolyzes deoxyuridine N-glycosidic bond; current formulation/units not substituted for2021 reagent'
SOURCE='Isakova et al.2021 DOI10.1073/pnas.2113568118; published PMC8713755 Methods111–133, SI FigS2/TableS2 and Supplementary Protocol'
TSO_WRITTEN={1:'dUCGdUCGGCAGCGdUCAGdUdUGdUAdUCAACdUCAGACAdUrGrG+G',2:'AdUGGCdUCGGAGAdUGdUGdUAdUAAGAGACAGdUCdUrGrG+G'}

@dataclass(frozen=True)
class Residue:
    base: str
    sugar: str
    abasic: bool=False
    def __post_init__(self):
        if self.base not in 'ACGTU' or self.sugar not in {'DNA','RNA','LNA'}:raise ValueError('invalid nucleotide')
        if self.abasic and (self.base,self.sugar)!=('U','DNA'):raise ValueError('only selected dU glycosylase product represented')
    def udg(self):return replace(self,abasic=True) if (self.base,self.sugar)==('U','DNA') else self

@dataclass(frozen=True)
class TSO:
    residues: tuple[Residue,...]
    biotin5: bool=True
    backbone_retained: bool=True
    def __post_init__(self):
        if not self.residues or not self.backbone_retained:raise ValueError('TSO backbone cannot disappear through base excision')
    @classmethod
    def source(cls,version):
        raw=TSO_WRITTEN[version];tokens=re.findall(r'dU|rG|\+G|[ACGT]',raw)
        if ''.join(tokens)!=raw:raise ValueError('unparsed source modification')
        return cls(tuple(Residue('U','DNA') if t=='dU' else Residue('G','RNA') if t=='rG' else Residue('G','LNA') if t=='+G' else Residue(t,'DNA') for t in tokens))
    def udg(self):return replace(self,residues=tuple(r.udg() for r in self.residues))
    @property
    def abasic_positions(self):return tuple(i for i,r in enumerate(self.residues) if r.abasic)
    def complementary_copy(self):
        if self.abasic_positions:raise ValueError('abasic-template traversal requires a separate specified path')
        comp={'A':'T','C':'G','G':'C','T':'A','U':'A'}
        return ''.join(comp[r.base] for r in reversed(self.residues))

@dataclass(frozen=True)
class TailPath:
    rna: str
    endogenous_tail: int=0
    added_tail: int=0
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU') or self.endogenous_tail<0 or self.added_tail<0:raise ValueError('explicit RNA and nonnegative tail lengths required')
    def extend_selected(self,number):
        if not isinstance(number,int) or number<=0:raise ValueError('selected positive addition length required; source yield distribution unknown')
        return replace(self,added_tail=self.added_tail+number)
    @property
    def tailed_rna(self):return self.rna+'A'*(self.endogenous_tail+self.added_tail)

def oligo_dt(version):
    if version==1:
        return Construct([Segment('v1 handle','CATA'+nx.ADAPTOR_S7),Segment('dT30','T'*30),Segment('anchor','VN',placeholder=True)])
    if version!=2:raise ValueError('source versions1/2 only')
    return Construct([Segment('v2 handle',nx.S5+'AGTTGTATCAACTCAGAC'),Segment('constrained UMI','NNNHNNNHNNHNNHNN',placeholder=True,feature=feature('smart_total_v2_umi','umi','degenerate')),
       Segment('local linker','GCTGTAAC'),Segment('dT22','T'*22),Segment('anchor','VN',placeholder=True)])

def graph_for_version(version=2,collection=False):
    if version not in {1,2}:raise ValueError('versions1/2 only')
    g=ProtocolGraph('Smart-seq-total v'+str(version)+' — RNA/TSO prefix',authority=SOURCE,representation='published_instance')
    if collection:
        raw=g.input('sorted_cells','Actual sorted single-cell collection across source eight384-well plates; occupied well counts/identities unassigned',K.BIOLOGICAL,primary=True)
        lysis=g.input('lysis_allocations','Per-well0.3ul source lysis buffer allocations:4U RRI/0.12percentTriton/TE; total quantity/occupied well multiplicity unassigned. Biotin oligo(dT) component declared separately below from same lysis inventory.',K.MOLECULAR)
        ercc=g.input('ercc','Actual per-well0.1ul1:300000 ERCC spike allocations, total occupied well count unassigned',K.MOLECULAR)
        prior=g.product('cell_rna','Compartment-resolved RNA collection after per-well72°C3min lysis; no cell/well identities or original RNA end classes invented',K.MOLECULAR)
        g.react('cell_lysis','Per-well72°C3min lysis across actual reaction collection; cool on ice',S.MAGIC,((raw,R.TEMPLATE),(lysis,R.STRUCTURAL),(ercc,R.SPIKE_IN)),(prior,),reason='Array of source single-cell reactions, not a cloned representative. Per-well quantities do not imply all wells occupied.')
        g.declare_reviewed_scope(('sorted_cells','lysis_allocations','ercc','cell_lysis'),SOURCE)
    else:
        prior=g.input('cell_rna','Single-cell total RNA after source72°C3min lysis; original RNA ends/circular species retained as unassigned populations',K.MOLECULAR,primary=True)
    polya=g.input('polya','E.coli poly(A) polymerase, source version-specific activity,37°C15min',K.CATALYST)
    atp=g.input('atp','ATP source version-specific concentration; incorporation length unassigned',K.MOLECULAR)
    tailed=g.product('tailed','Selected eligible linear RNA obtains an artificial3′poly(A) extension; endogenous tail origin no longer identifiable from total tail alone',K.MOLECULAR,chemical_state={'added_tail':'unassigned length','native_RNA_retained':True,'v2_capping':'reported simultaneous treatment; substrate-specific cap chemistry unassembled' if version==2 else 'not part of v1 treatment'})
    cap_inputs=()
    if version==2:
        cap=g.input('capping_enzyme','Vaccinia capping enzyme0.1U in source simultaneous tail/cap mix',K.CATALYST)
        donors=g.input('cap_donors','GTP0.05mM and SAM2.5nM; source cap chemistry/eligible RNA ends unassembled',K.MOLECULAR)
        cap_inputs=((cap,R.CATALYST),(donors,R.STRUCTURAL))
    g.react('tail','Source polyadenylation'+(' with simultaneous capping treatment; cap topology/yields unassigned' if version==2 else '')+'; no guarantee all RNA/circular species acquire an extendable tail',S.TAILING,((prior,R.TEMPLATE),(polya,R.CATALYST),(atp,R.STRUCTURAL))+cap_inputs,(tailed,))
    primer=g.input('oligo_dt','Source biotinylated anchored oligo-dT'+(' with constrained16baseUMI' if version==2 else ''),K.OLIGO,state=oligo_dt(version),chemical_state={'modification':'source5primebiotin; present in initial lysis allocations, not freshly added after tailing'})
    rt=g.input('ssii','Actual source SuperScriptII15U per-well allocations,42°C90min; same reaction across template switching',K.CATALYST)
    extra_rt=()
    if collection:
        rtnt=g.input('rt_dntp','Actual per-well source RT dNTP allocations1.5mM; component of source1ul RT mix, not a fresh additional dose',K.MOLECULAR)
        rtri=g.input('rt_rri','Actual per-well4U RRI allocation within same source1ul RT mixture',K.CATALYST)
        rtb=g.input('rt_buffer','Actual per-well RT allocation components:1.5xFirstStrand,7.5mMDTT,1.5Mbetaine,10mMMgCl2; same1ul source RT mix. Packaging/total volume unassigned.',K.MOLECULAR)
        extra_rt=((rtnt,R.STRUCTURAL),(rtb,R.STRUCTURAL),(rtri,R.CATALYST))
        g.declare_reviewed_scope(('rt_dntp','rt_buffer','rt_rri'),SOURCE)
    cdna=g.product('firststrand','Selected oligo-dT-primed cDNA with RNA template retained; actual priming sites and residual RNA classes unassigned',K.MOLECULAR,chemical_state={'RNA_retained':True,'primer':oligo_dt(version)})
    g.react('rt','Source first-strand synthesis; no undocumented RNaseH removal',S.REVERSE_TRANSCRIPTION,((tailed,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST))+extra_rt,(cdna,))
    tso=g.input('tso','Source5′biotin TSO; TableS2 dU/RNA/LNA identity retained',K.OLIGO,chemical_state=TSO.source(version))
    switched=g.product('switched','Conditional TSO-templated cDNA extension; original RNA and free/extended TSO populations retained',K.MOLECULAR,chemical_state={'RNA_retained':True,'tso':TSO.source(version),'copied_TSO_DNA':TSO.source(version).complementary_copy(),'copy_scope':'selected complete TSO-templated path, efficiency/truncation unassigned'})
    g.react('switch','Source template switching; typed TSO is template, not a DNA molecule silently substituted by its sequence',S.TEMPLATE_SWITCHING,((cdna,R.TEMPLATE),(tso,R.STRUCTURAL),(rt,R.CATALYST)),(switched,))
    udg=g.input('udg','UDG '+('1U37°C30min' if version==1 else '2U37°C60min'),K.CATALYST)
    excised=g.product('abasic_tso','Selected dU base-excision products; TSO backbone retained with AP sites and liberated uracil; cDNA/RNA remain in reaction',K.MOLECULAR,chemical_state={'tso':TSO.source(version).udg(),'liberated_bases':'uracil from selected dU sites; not yet physically removed','cDNA':'copied canonical DNA retained','RNA_retained':True})
    g.react('udg_excision','Source calls this TSO digestion; local model resolves base excision only, backbone breakdown not inferred',S.DIGESTION,((switched,R.TEMPLATE),(udg,R.CATALYST)),(excised,))
    g.mark_terminal(Terminal(excised.id,Disposition.UNRESOLVED,'Remaining TSO/AP processing and source PCR/index/tagmentation branches, physical cleanup, DASH and acquired read plan unassembled. v2 simultaneous capping not yet reconstructed; no universal5′end eligibility assumed.'))
    g.declare_reviewed_scope(['cell_rna','polya','atp','tail','oligo_dt','ssii','rt','tso','switch','udg','udg_excision'],SOURCE)
    if version==2:g.declare_reviewed_scope(['capping_enzyme','cap_donors'],SOURCE)
    g.add_source_item(SourceItem('unreviewed_remaining_preparation','Version-specific capping, cDNA/PCR/index branches, DASH and acquired read geometry','Published Methods123–133; SI Supplementary Protocol and TableS2'))
    if version==1:
        g.annotate(Highlight('tso','Source TSO notation: explicit dU table used',
            'Published main Methods117/PDFp6 prints UCGUCGGCAGCGUCAGUUGUAUCAACUCAGACAUrGrG+G, whereas SI TableS2p37 explicitly marks internal deoxyU as dU. Bare U in the main oligo string does not by itself specify RNA sugar and can abbreviate dU; no RNA-versus-DNA conflict or publication error is established. Local ordered residues follow the explicit TableS2 notation, without claiming an independently verified order sheet.'))
    if version==2:
        g.annotate(SourceDiscrepancy('abasic_tso','Published v2 custom Read1 values differ',
            (SourceClaim('Published main Methods133/PDFp7','TCGTCGGCAGCGTCAGTTGTATCAACTCAGAC'),SourceClaim('Published SI TableS2 originalp37','TCGGCAGCGTCAGTTGTATCAACTCAGACATGGG')),
            'Source read cycles are declared, but occupied primer geometry and actual custom Read1 identity remain unassembled; neither sequence is silently substituted.',
            'Both original PDF locations inspected; v2 row in SI repeats the v1 custom Read1 value.',
            'Actual acquired primer/reagent identity and branch applicability remain unresolved.'))
    return g

def protocol_graph():return full_v2_graph()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();t=TSO.source(2);u=t.udg()
    c('UDG preserves backbone residue positions and terminal RNA/LNA', (len(u.residues),u.residues[-3:]),(len(t.residues),t.residues[-3:]))
    c('only DNA uracils become abasic',all((r.base,r.sugar)==('U','DNA') for r in u.residues if r.abasic),True)
    c('RNA uracil is not silently made UDG-susceptible',Residue('U','RNA').udg().abasic,False)
    try:u.complementary_copy();bad=False
    except ValueError:bad=True
    c('AP traversal cannot masquerade as ordinary complete copying',bad,True)
    a=TailPath('CGU',2).extend_selected(3);b=TailPath('CGU',0).extend_selected(5)
    c('same acquired tail cannot determine endogenous polyA status',(a.tailed_rna==b.tailed_rna,a.endogenous_tail!=b.endogenous_tail),(True,True))
    c('v2 UMI keeps constrainedH synthesis',next(s.top for s in oligo_dt(2) if s.feature),'NNNHNNNHNNHNNHNN')
    original=Source(sources_root()/'smart-seq-total__10.1073+pnas.2113568118'/'pnas.2113568118_supp_pnas.2113568118.sapp.pdf.txt','SI TableS2, rendered originalp37 inspected')
    if have(c,original,label='original oligos'):
        text=''.join(original.path.read_text().split())
        for v in (1,2):c('sourceTSOv'+str(v)+' full modification spelling',TSO_WRITTEN[v] in text,True)
        c('sourcev2 anchored UMI primer sequence',''.join(s.top for s in oligo_dt(2)) in text,True)
    full=protocol_graph()
    # TableS2 v2 FW/RV cDNA Amp and index Amp are produced reagent components,
    # so the external-input oligo guard does not inspect them.
    if have(c,original,label='Produced v2 PCR reagent source transcription'):
        text=''.join(original.path.read_text().split())
        for key,written in (('preamplify_fw','GCTCGGAGATGTGTATAAGAGACAG'),
                            ('preamplify_rv','TCGTCGGCAGCGTCAGTTG'),
                            ('index_pcr_fw','GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTC')):
            actual=''.join(segment.top for segment in full.materials[key].state)
            c('TableS2 whole produced '+key,(actual,written in text),(written,True))
    c('Original acquired geometry obligation remains required',full.publication_problems(),["authoritative source item 'unreviewed_remaining_preparation' is not accounted for"])
    c('Eight distinct normalized source plate aliquots enter final pool',tuple(p.material for p in full.steps['eight_plate_pooling'].inputs),tuple('plate'+str(i)+'_pool_aliquot_selected' for i in range(1,9)))
    c('Final mixed sequencing pool retains untransferred stock',full._consumers['loading_remainder'],[])
    c('ERCC stock contributes through physical primary preparation','ercc' in full._contributing_origins('loading_selected'),True)
    c('UDG treated mixture reaches PCR before physical removal',full._consumers['abasic_tso'],['preamplify'])
    c.report()


# Publication wiring deferred until the original required source/model scope is assembled.

def full_v2_graph():
    from protocol_graph import SequencingPlan,SourceCoverage
    from chemdraw import Scene
    from pcr import PCRGeometry
    g=graph_for_version(2,collection=True)
    # The prefix represents the complete reaction collection, not one cell copied eight times.
    del g.terminals['abasic_tso']
    def inp(k,label,kind=K.MOLECULAR):return g.input(k,label,kind)
    def op(k,label,x,extras=(),kind=S.MAGIC,outkind=K.MOLECULAR,geometry=None):
        out=g.product(k+'_product',label+'; whole occupied library molecules and recovery distributions unassembled',outkind)
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(extras),(out,),geometry=geometry,reason='Source-selected v2 direct-index workflow. Compartment arrays retained; unknown whole source polymers do not supply a canonical kit recipe.' if kind is S.MAGIC else '')
        return out
    def clean(k,x,ratio):
        b=inp(k+'_beads',ratio+' source AMPure allocations per corresponding compartment',K.SELECTION)
        out=g.product(k+'_kept','Recovered DNA collection after '+ratio+' cleanup',K.MOLECULAR);w=g.product(k+'_removed','Corresponding physically removed bead/wash fractions',K.WASTE)
        g.react(k,'Source '+ratio+' cleanup of corresponding compartment collection',S.CLEANUP,((x,R.TEMPLATE),(b,R.SELECTION)),(out,w));g.discard(w,'Physical removed fractions, not catalytic disappearance.');return out
    def select(k,x,label,qc=False):
        take=g.product(k+'_selected',label,K.QC if qc else x.kind);rest=g.product(k+'_remainder','Remaining source stock in matching compartments; volume unassigned',x.kind)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(take,rest),changes_molecule=False)
        if qc:
            g.mark_terminal(Terminal(take.id,Disposition.QC,'Actual source assay fraction; quantity unassigned.'));return rest
        g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferred source stock, no full-stock loading inference.'));return take
    def pcr(k,x,programme,fw,rv,index=False):
        stock=inp(k+'_supplies',('Aggregate actual per-well index-PCR reagent allocations; final4ul including1.5ul cDNA, packaging/stock count unassigned' if index else 'Actual per-well3.2ul source PCR mix allocations, final1x KAPA HiFi plus0.5uM forward/reverse primers; occupied well count/total supply unassigned'))
        polymerase=g.product(k+'_enzyme','KAPA HiFi polymerase component from actual source KAPA supply allocations',K.CATALYST)
        nt=g.product(k+'_nt','Nucleotide component of same actual KAPA supply allocations',K.MOLECULAR);buf=g.product(k+'_buf','Buffer component of same actual KAPA supply allocations',K.CONDITION)
        f=g.product(k+'_fw','Source forward/IndexAmp primer component',K.OLIGO,state=Construct([Segment('source primer',fw)]))
        r=g.product(k+'_rv',('Source reverse primer component' if not index else 'Actual0.4ul each5uM i5/i7 primer allocations; selected index wholeorders unassembled'),K.OLIGO,state=None if index else Construct([Segment('source primer',rv)]))
        g.react(k+'_components','Expose components of actual source PCR allocations',S.MAGIC,((stock,R.TEMPLATE),),(polymerase,nt,buf,f,r),changes_molecule=False,reason='Source PCR component inventory. Same actual KAPA MasterMix component views; primer allocations not fresh additional doses or an invented premix. Index species/count and physical packaging unassigned.')
        sc=Scene();sc.strand('unassembled',[Segment('occupied cDNA sites','[unassembled]',placeholder=True)],label='Source-positive PCR, geometry unassembled')
        geom=PCRGeometry.unavailable(sc,reason='Whole occupied cDNA/selected indexes and primer landing remain unassembled; no canonical finished library is borrowed.')
        return op(k,programme,x,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(polymerase,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION)),S.PCR,geometry=geom)
    x=pcr('preamplify',g.materials['abasic_tso'],'Per-well95°C3min;21cycles98°C20s/65°C15s/72°C5min;72°C5min; AP-TSO/RNA remains until actual cleanup','GCTCGGAGATGTGTATAAGAGACAG','TCGTCGGCAGCGTCAGTTG')
    x=clean('cdna_cleanup',x,'1.8x')
    x=select('picogreen',x,'Per-well actual PicoGreen assay aliquots; volume unassigned',True)
    dilution=inp('dilution_water','Source dilution supplies to0.2ng/ul in corresponding wells',K.CONDITION)
    x=op('dilute','Dilute measured cDNA stocks to0.2ng/ul, keeping distinct occupied wells',x,((dilution,R.CONDITION),))
    x=select('index_template',x,'Actual1.5ul diluted cDNA selected per occupied well,0.3ng calculated template')
    x=pcr('index_pcr',x,'Per-well source direct-index4ul PCR:95°C3min;10cycles98°C20s/62°C15s/72°C5min;72°C5min','GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTC','',True)
    # Explicit actual eight-plate collection boundary, not a physical split of one pooled molecule.
    plates=tuple(g.product('plate_'+str(i),'Distinct per-well indexed reaction products belonging to actual plate '+str(i)+'; occupied well count/index assignments unassigned',K.MOLECULAR) for i in range(1,9))
    g.react('plate_array','Expose eight actual plate groups of the original compartment collection',S.MAGIC,((x,R.TEMPLATE),),plates,changes_molecule=False,reason='Array boundary only. These are distinct reactions from the outset, not a split of one amplified sample or eight copies of one chosen well.')
    finished=[]
    for i,p in enumerate(plates,1):
        key='plate'+str(i)
        p=select(key+'_contributions',p,'Selected per-well library contributions for actual within-plate pool; amounts unassigned')
        p=op(key+'_pool','Physically pool distinct wells within plate'+str(i),p)
        p=clean(key+'_pool_cleanup',p,'0.8x')
        guides=inp(key+'_guides','Actual source guideRNA/tracrRNA allocations for45S+rRNA targets; TableS1 identities, whole orders/selected dose unassembled',K.OLIGO)
        cas=inp(key+'_cas9','Source Cas9 allocation; main guide/tracrRNA/Cas9 ratio wording retained without new stoichiometry',K.CATALYST)
        cb=inp(key+'_cas9_buffer','Source1xCas9 buffer in actual plate-pool reaction',K.CONDITION)
        p=op(key+'_dash','Main selected DASH37°C1h; cut and uncut DNA remain physically together',p,((guides,R.STRUCTURAL),(cas,R.CATALYST),(cb,R.CONDITION)),S.DIGESTION)
        pk=inp(key+'_proteinase_supply','Actual source proteinase K allocation',K.CATALYST)
        p=op(key+'_proteinase','Main selected proteinase K50°C15min; protein breakdown does not physically remove products',p,((pk,R.CATALYST),),S.DEPROTEINIZATION)
        p=clean(key+'_dash_cleanup',p,'1.2x')
        gel=inp(key+'_pippin','Actual source2percent Pippin gel selection210–600bp',K.SELECTION)
        kept=g.product(key+'_size_kept','Size-selected DNA from actual plate pool',K.MOLECULAR);reject=g.product(key+'_size_removed','Physical size-selection rejects',K.WASTE)
        g.react(key+'_size_select','Pippin210–600bp physical selection',S.CLEANUP,((p,R.TEMPLATE),(gel,R.SELECTION)),(kept,reject));g.discard(reject,'Actual size-selection fraction.');p=kept
        p=select(key+'_fragment_qpcr',p,'Actual FragmentAnalyzer and KAPA qPCR aliquot collection; quantities unassigned',True)
        norm=inp(key+'_normalization','Source normalization dilution supplies to2nM',K.CONDITION)
        p=op(key+'_normalize','Normalize measured actual plate pool to2nM',p,((norm,R.CONDITION),))
        p=select(key+'_pool_aliquot',p,'Equal-volume selected contribution from actual plate'+str(i)+'; volume unassigned')
        finished.append(p)
    pool=g.product('eight_plate_pool','Physical pool of equal volumes from all eight actual plate pools',K.LIBRARY)
    g.react('eight_plate_pooling','Mix eight actual normalized plate-pool aliquots',S.MAGIC,tuple((p,R.TEMPLATE) for p in finished),(pool,),changes_molecule=False,reason='Actual eight-plate source pool; no invented pool members or384 occupied-well assumption.')
    phix=inp('phix','Source10percent PhiX sequencing spike-in; stock/input-volume convention unassigned',K.LIBRARY)
    pool=op('phix_mix','Add source10percent PhiX to sequencing pool',pool,((phix,R.TEMPLATE),),outkind=K.LIBRARY)
    pool=select('loading',pool,'Final pooled-library loading aliquot; untransferred mixed pool retained')
    g.sequence(pool,SequencingPlan(acquired_roles=('Read1','Read2','Index1','Index2'),read_lengths={'Read1':30,'Read2':112,'Index1':12,'Index2':12},read_uses={'Read1':'include','Read2':'include','Index1':'include','Index2':'include'},limitation='Source v2 NextSeq500 R1=30/R2=112 and two12cycle index reads. Exact occupied library geometry and Read1 reagent identity unresolved: main and TableS2 custom Read1 differ. UMI16 constrained bases from RT primer; whole read landing not fabricated.'))
    g.declare_reviewed_scope(('preamplify_supplies','preamplify_components','preamplify','cdna_cleanup_beads','cdna_cleanup','picogreen','dilution_water','dilute','index_template','index_pcr_supplies','index_pcr_components','index_pcr','plate_array','eight_plate_pooling','phix','phix_mix','loading','loading_selected'),SOURCE)
    # Fixed source-declared eight plate groups, independently reviewed preparation inventory.
    for i in range(1,9):
        k='plate'+str(i)
        g.declare_reviewed_scope(tuple(k+s for s in ('_contributions','_pool','_pool_cleanup_beads','_pool_cleanup','_guides','_cas9','_cas9_buffer','_dash','_proteinase_supply','_proteinase','_dash_cleanup_beads','_dash_cleanup','_pippin','_size_select','_fragment_qpcr','_normalization','_normalize','_pool_aliquot')),SOURCE)
    g.annotate(Highlight('plate1_dash','Main and general SI DASH recipes differ','Main selected route: proteinase K50°C15min and1.2xAMPure. SI general DASH protocol:56°C10min and0.8x/1xSPRI; its printed88+88+704uL entries sum880uL while total prints1140uL. This model uses main v2 branch conditions; SI general bulk stock recipe is not silently merged.'))
    return g

def render_page():
    from protocol_graph import AnnotationSlot
    from page import head,info
    from protocol_views import assemble_page,WorkflowSection,OligoSection,OligoEntry
    return assemble_page(protocol_graph().require_publishable(),(head('Smart-seq-total v2'),'<div class="wrap" data-protocol-complete="true"><h1>Smart-seq-total v2</h1>',AnnotationSlot.INTRO,
      '<p class="research-notes"><a href="01_smart_total.html">Research notes</a></p>',
      '<p>Single-cell RNA is polyadenylated with simultaneous capping treatment, primed with a biotinylated oligo(dT) bearing a constrained sixteen-base UMI, and copied by template-switching reverse transcription. UDG excises deoxyuracil bases from the TSO; it does not by itself remove its backbone. cDNA amplification, cleanup and direct index PCR precede within-plate pooling, DASH depletion and physical size selection.</p>',
      info('Source v2 direct-index route; no tagmentation branch is imported. The model preserves actual compartment collections and eight plate pools, with occupied well counts unassigned. Whole occupied RNA/cDNA/library structures, eligible capping substrates and amplification efficiency remain unassembled. Unknown structures do not imply absent source information.'),
      OligoSection('Source v2 anchored RT primer', (OligoEntry('5′ biotinylated v2 oligo(dT)',tuple(oligo_dt(2))),),note_html='The ordered primer has a5′biotin. NNNHNNNHNNHNNHNN is the source constrained16-base UMI, followed by source linker andT22VN. Modification is not an extra nucleotide.'),
      '<h2>Modified TSO and read setup</h2><p>Source v2 TSO:5′biotin-AdUGGCdUCGGAGAdUGdUGdUAdUAAGAGACAGdUCdUrGrG+G-3′. Internal dU is DNA uracil; terminal rG/rG/+G retains RNA/RNA/LNA sugar identities. Base excision preserves backbone positions in the local chemical declaration. NextSeq500 acquires R1=30,R2=112 and two12-cycle index reads.</p>',
      '<h2>Published custom Read1 discrepancy</h2><p>Main Methods133 prints TCGTCGGCAGCGTCAGTTGTATCAACTCAGAC; original TableS2p37 prints TCGGCAGCGTCAGTTGTATCAACTCAGACATGGG, repeating its v1 custom-primer value. Actual acquired v2 reagent identity and precise read landing remain unresolved. Neither value is silently substituted.</p>',
      '<h2>Actual pools and source conditions</h2><p>Measured cDNA supplies1.5ul per occupied well for direct indexing. Distinct well products pool within their original plates; plate pools undergo guide/Cas9 digestion, proteinase K, cleanup and210–600bp Pippin selection. Eight normalized2nM plate pools contribute equal volumes;10percent PhiX is added before selecting a final loading aliquot. All QC fractions and untransferred stocks remain explicit.</p>',
      WorkflowSection('Full source-selected v2 preparation and source discrepancies'),AnnotationSlot.AFTER_WORKFLOW,'</div>'))

if __name__=='__main__':selftest()
