"""Source-owned m6A-SAC-seq chemical branches; modified RNA is not a basecall."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, Construct, feature, revcomp, panel
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_views import assemble_page, WorkflowSection
from page import head, info
PAGE_FILE="m6a-sac.html"
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition, SequencingPlan
SOURCE='Hu et al.2022 doi:10.1038/s41587-022-01243-z, accepted manuscript Methods211–227; publisher SupplementaryTable1; VOR procedure applicability unverified'
SOURCE_OLIGO_EXCEPTIONS={'rna_adapter':'Publisher Table 1 whole RNA adapter ordering has /iBiodT/, an actual biotinylated deoxythymidine nucleotide. The model retains its final physical T; the generic source parser drops this internal modified nucleotide. Narrow transcription-format exception only, verified by the whole modified source fixture and physical-backbone comparison; not an absence or chemistry exception.'}
RNA_ADAPTER_WRITTEN='rAppNNNNNATCACGAGATCGGAAGAGCACACGTC/iBiodT//3SpC3/'
CDNA_ADAPTER_WRITTEN='Phos-'+'N'*5+revcomp(il.TRUSEQ_SMALL_RNA_RA5)+'/3SpC3/'
RT_PRIMER=il.TRUSEQ_READ2[-21:]

@dataclass(frozen=True)
class AdenosinePath:
    original: str
    stage: str='native'
    fto_erased: bool=False
    def __post_init__(self):
        if self.original not in {'A','m6A'} or self.stage not in {'native','allyl','cyclized'}:
            raise ValueError('unsupported source adenosine pathway')
        if self.fto_erased and self.original!='m6A':raise ValueError('FTO-erased path must originate at m6A')
    @property
    def substrate(self):return 'A' if self.fto_erased else self.original
    def erase_selected(self):
        if self.original!='m6A' or self.stage!='native' or self.fto_erased:raise ValueError('selected FTO path requires native m6A')
        return replace(self,fto_erased=True)
    def allylate(self):
        if self.stage!='native':raise ValueError('only native A/m6A substrate enters this labeling event')
        return replace(self,stage='allyl')
    def cyclize(self):
        if self.stage!='allyl':raise ValueError('iodine cyclization requires the allylated substrate')
        return replace(self,stage='cyclized')
    @property
    def species(self):
        if self.stage=='native':return self.substrate
        label='a6m6A' if self.substrate=='m6A' else 'a6A'
        return ('cyclized '+label) if self.stage=='cyclized' else label
    def as_canonical_rna_base(self):
        if self.stage!='native' or self.substrate!='A':
            raise ValueError('a covalent RNA adduct/modification is not a canonical sequencing base')
        return 'A'

@dataclass(frozen=True)
class RNAPopulation:
    control: str
    selected_reacted_paths: tuple[AdenosinePath,...]
    residual_native_material: str='unassigned fractions retained; not complete conversion'
    soluble_material: str='present until physical washes; exact composition unassigned'


def published_adapter_inventory_scene():
    sc=Scene()
    sc.strand('stock adapter',[
        Segment('random identifier window','N'*5,placeholder=True,
          feature=feature('sac_random_identifier','inline_barcode','random',note='Source random window; molecule-UMI interpretation remains unassigned')),
        Segment('printed index','ATCACG',feature=feature('sac_printed_index','sample_index','fixed',note='Printed stock example; no actual study sample assignment')),
        Segment('adapter DNA','AGATCGGAAGAGCACACGTC'),
        Segment('internal biotin dT','T',note='Physical nucleotide /iBiodT/')],
        label='RNA 3′ adapter, free stock',mod5='rApp',mod3='SpC3')
    sc.ligation_end('stock adapter','random identifier window',kind='single-strand',polarity='5-prime',owner='adapter',text='5rApp donor; activation AMP is released')
    sc.note('stock adapter','/iBiodT/ modifies the final T; 3SpC3 is a non-nucleotide block.')
    return sc

def cdna_adapter_scene():
    sc=Scene();sc.strand('adapter',[Segment('random window','N'*5,placeholder=True,feature=feature('sac_cdna_random','inline_barcode','random',note='Random adapter window; no automatic UMI or sample assignment')),Segment('source adapter body',revcomp(il.TRUSEQ_SMALL_RNA_RA5))],mod5='P',mod3='SpC3',label='Source Table 1 cDNA adapter; free donor inventory, not mature library')
    sc.ligation_end('adapter','random window',kind='single-strand',polarity='5-prime',owner='adapter',text='Known source donor 5P; complete cDNA acceptor and sealed product unassembled')
    return sc

def unknown_scene():
    sc=Scene();sc.strand('context',[Segment('whole occupied molecule unassembled','[unassembled molecule]',placeholder=True,inferred=True)])
    return sc

class Preparation:
    def __init__(self,control):
        self.g=ProtocolGraph('m6A-SAC-seq — '+control,authority=SOURCE,representation='published_instance')
        self.x=self.g.input('rna','Actual independently prepared '+control+' RNA library input: 30–100 ng poly(A)+ or rRNA-depleted RNA; actual transcript identities unassembled',MaterialKind.MOLECULAR,primary=True)
    def material(self,key,label,kind=MaterialKind.MOLECULAR):return(self.g.input(key,label,kind),InputRole.SELECTION if kind is MaterialKind.SELECTION else InputRole.STRUCTURAL)
    def catalyst(self,key,label):return(self.g.input(key,label,MaterialKind.CATALYST),InputRole.CATALYST)
    def condition(self,key,label):return(self.g.input(key,label,MaterialKind.CONDITION),InputRole.CONDITION)
    def op(self,key,label,kind,extra=(),chemical_state=None,geometry=None,inferred=False):
        out=self.g.product(key+'_product','Actual '+label+' preparation; whole occupied strands, carryover and recovery unassembled',MaterialKind.MOLECULAR,chemical_state=chemical_state,inferred=inferred)
        kw={}
        if geometry is not None:kw['geometry']=geometry
        if kind is StepKind.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,InputRole.TEMPLATE),)+tuple(extra),(out,),**kw);self.x=out;return out
    def clean(self,key,label,extra=(),inferred=False):
        out=self.g.product(key+'_kept','Actual retained '+label+' material; recovery/carryover unassigned',MaterialKind.MOLECULAR,inferred=inferred)
        waste=self.g.product(key+'_waste','Actual excluded '+label+' workup fractions',MaterialKind.WASTE,inferred=inferred)
        self.g.react(key,label,StepKind.CLEANUP,((self.x,InputRole.TEMPLATE),)+tuple(extra),(out,waste));self.g.discard(waste,'Actual source workup fractions; no perfect purification asserted');self.x=out

def source_supplies(b,key,label,enzyme_label,dyed=False):
    g=b.g;stock=g.input(key,label,MaterialKind.MOLECULAR)
    enzyme=g.product(key+'_enzyme',enzyme_label,MaterialKind.CATALYST)
    nts=g.product(key+'_nts','Copresent source nucleotide component; no extra dose',MaterialKind.MOLECULAR)
    buffer=g.product(key+'_buffer','Copresent source buffer component; formulation/stock quantities unassigned',MaterialKind.CONDITION)
    outputs=[enzyme,nts,buffer]
    if dyed:outputs.append(g.product(key+'_dye','Copresent qPCR detection contribution; identity/amount unassigned',MaterialKind.MOLECULAR))
    g.react(key+'_inventory','Same actual supplies: component bookkeeping, not physical separation or a claimed premix factory',StepKind.GENERIC,((stock,InputRole.TEMPLATE),),tuple(outputs),changes_molecule=False)
    return tuple((x,InputRole.CATALYST if x is enzyme else InputRole.CONDITION if x is buffer else InputRole.STRUCTURAL)for x in outputs)

# Source main text219: initial1h; fresh20min recipe followed by 'this labeling step was repeated six times'.
# Local interpretation: six fresh20min allocations TOTAL after initial1h. The prose can also be read as six repeats after the separately described first20min round; public note preserves this ambiguity.
LABEL_ROUNDS=(('label_0','Initial 20 µl labeling mixture: 6 µl bead water +2 µl 10x buffer +2 µl inhibitor +6 µl allylic-SAM +4 µl MjDim1 (1.6 mM); 50°C 1 h'),
 ('label_1','Fresh round 1: 4 µl water +1 µl buffer +1 µl inhibitor +2 µl allylic-SAM +2 µl enzyme; 50°C 20 min'),
 ('label_2','Fresh round 2: same actual source recipe; 50°C 20 min'),('label_3','Fresh round 3: same actual source recipe; 50°C 20 min'),
 ('label_4','Fresh round 4: same actual source recipe; 50°C 20 min'),('label_5','Fresh round 5: same actual source recipe; 50°C 20 min'),('label_6','Fresh round 6: same actual source recipe; 50°C 20 min'))

def graph_for_control(control):
    if control not in {'input','FTO-','FTO+'}:raise ValueError('source controls are input,FTO-,FTO+')
    b=Preparation(control);g=b.g
    b.op('polya_remove','Oligo(dT) annealing; RNase H/DNase I treatment',StepKind.DIGESTION,
         (b.material('oligodt','Actual oligo(dT) reagent; whole ordering unassigned',MaterialKind.OLIGO),b.catalyst('rnaseh_upstream','Actual RNase H'),b.catalyst('dnase','Actual DNase I')))
    b.clean('rna_cleanup','RNA Clean & Concentrator recovery',(b.material('rna_column','Actual source cleanup materials; exact inventory unassembled',MaterialKind.SELECTION),))
    b.op('sonicate','Source Bioruptor 30 cycles 30 s on/off; roughly 150 nt fragments',StepKind.FRAGMENTATION,(b.condition('sonication_condition','Actual source sonication operation; size distribution unassembled'),))
    b.op('pnk','PNK treatment 37°C 30 min; source intended 3OH end class',StepKind.CHEMICAL_CONVERSION,(b.catalyst('pnk_enzyme','Actual PNK; ATP and 5P creation not inferred'),))
    b.op('spike','Add actual 0.6% calibration RNA mix',StepKind.POOL,((g.input('calibration','Actual calibration RNA population; independent of endogenous target methylation/sequence',MaterialKind.MOLECULAR),InputRole.TEMPLATE),))
    adapter=g.input('rna_adapter','Actual source preadenylated RNA 3′ DNA adapter; Table 1 component-qualified, study index member unassigned',MaterialKind.OLIGO,state=published_adapter_inventory_scene(),chemical_state={'raw_source_order':RNA_ADAPTER_WRITTEN,'modifications':'5rApp; internal physical biotin-dT; non-nucleotide 3SpC3','scope':'Known free stock pattern; reference index not assigned to study sample; full occupied RNA/end occupancy unassembled'})
    b.op('rna_ligate','RNA 3′ adapter ligation; activation AMP/coproducts remain in reaction until workup',StepKind.LIGATION,((adapter,InputRole.STRUCTURAL),b.catalyst('rnl2','Actual T4 RNA ligase 2 truncated KQ')),geometry=AdapterLigationGeometry.known('single-strand',published_adapter_inventory_scene(),statement='Source activated 5rApp donor; intended RNA 3OH class. Whole RNA, occupancy and sealed junction unassembled.'))
    b.op('deadenylate','Add 1 µl 5′ deadenylase; 30°C 1 h',StepKind.DIGESTION,(b.catalyst('deadenylase','Actual source 1 µl 5′ deadenylase'),))
    b.op('recjf_digest','Then add 1 µl RecJf; 37°C 1 h',StepKind.DIGESTION,(b.catalyst('recjf','Actual source 1 µl RecJf'),))
    primer=g.input('rt_primer','Actual 1 µl 50 µM RT primer; remains copresent after annealing/capture',MaterialKind.OLIGO,state=Construct([Segment('source RT primer',RT_PRIMER)]))
    b.op('primer_anneal','RT primer annealing:75°C5min/37°C15min/25°C15min',StepKind.BINDING,((primer,InputRole.STRUCTURAL),))
    b.clean('capture','15 µl C1 beads capture/wash of adapter-ligated RNA',(b.material('c1','Actual streptavidin C1 bead/protein contribution',MaterialKind.SELECTION),))
    population=RNAPopulation(control,(AdenosinePath('A'),AdenosinePath('m6A')))
    if control=='FTO+':
        population=replace(population,selected_reacted_paths=(*population.selected_reacted_paths,AdenosinePath('m6A').erase_selected()))
        b.op('demethylate','Actual on-bead FTO control; residual erasure fractions unassigned',StepKind.CHEMICAL_CONVERSION,(b.catalyst('fto','Actual 0.2 nmol FTO'),b.material('fto_supplies','Actual 50 µl demethylation mixture: HEPES/Fe(II)/ascorbate/alpha-ketoglutarate/RNase inhibitor; actual chemical cofactors, no new enzyme dose'),b.condition('fto_condition','Source 70°C30s denaturation/ice before FTO treatment; no complete erasure assumed')),chemical_state=population)
        b.clean('fto_wash','Actual post-FTO bead washing')
    if control!='input':
        for key,label in LABEL_ROUNDS:
            stock=g.input(key+'_supplies','INFERRED — repeat convention: source mixture under the disclosed six-round interpretation: '+label,MaterialKind.MOLECULAR,inferred=True)
            enzyme=g.product(key+'_enzyme','MjDim1 component under the repeat convention; exact source recipe retained',MaterialKind.CATALYST,inferred=True)
            sam=g.product(key+'_sam','Allylic-SAM component under the repeat convention; no extra dose',MaterialKind.MOLECULAR,inferred=True)
            other=g.product(key+'_other','Same repeat-conditioned mixture buffer/RNase inhibitor/water components; no new dose',MaterialKind.MOLECULAR,inferred=True)
            g.react(key+'_inventory','INFERRED — repeat convention: fresh mixture components; no physical isolation or extra additions',StepKind.GENERIC,((stock,InputRole.TEMPLATE),),(enzyme,sam,other),changes_molecule=False)
            if key=='label_0':population=replace(population,selected_reacted_paths=tuple(p.allylate()for p in population.selected_reacted_paths))
            b.op(key,'INFERRED — repeat convention: '+label,StepKind.CHEMICAL_CONVERSION,((enzyme,InputRole.CATALYST),(sam,InputRole.STRUCTURAL),(other,InputRole.STRUCTURAL)),chemical_state=population,inferred=True)
            b.clean(key+'_drain','INFERRED — repeat convention: supernatant removal after '+key+(' and wash/resuspend 25 µl water' if key=='label_6' else '; retained bead population receives next fresh mixture'),inferred=True)
        iodine=g.input('iodine','Actual 1 µl 125 mM molecular I2 reactant',MaterialKind.MOLECULAR)
        population=replace(population,selected_reacted_paths=tuple(p.cyclize()for p in population.selected_reacted_paths))
        b.op('cyclize','Iodine cyclization in dark at room temperature 1 h',StepKind.CHEMICAL_CONVERSION,((iodine,InputRole.STRUCTURAL),b.condition('iodine_environment','Source dark/room-temperature reaction environment')),chemical_state=population)
        b.op('iodine_quench','Source printed 1 µl 40 mM Na2S2SO3 quench; formula retained without normalization',StepKind.QUENCHING,(b.material('iodine_quencher','Actual source quench reagent printed Na2S2SO3; chemical identity beyond printed formula unassigned'),b.condition('quench_condition','Actual source quench operation with the same stated reagent; no separate dose')))
        b.clean('iodine_wash','Actual bead wash/resuspend in 9 µl water')
    # The RT primer was physically added before bead capture; do not add a new primer dose here.
    rt=g.input('hiv_rt','Actual HIV RT: 1 µl for input 1 h; 2 µl for treated samples 3 h',MaterialKind.CATALYST)
    b.op('rt','Source HIV RT at 37°C: '+('1 h input'if control=='input'else'3 h treated'),StepKind.REVERSE_TRANSCRIPTION,
       ((rt,InputRole.CATALYST),b.material('rt_nts','Actual 2 µl 10 mM dNTP'),b.material('rt_other','Actual RT contributions: 2 µl 10× buffer / 2 µl 25 mM Mg / 1.25 µl 0.1 M DTT / 2 µl RNaseOUT; source RT-primer-bearing template, no fresh primer addition')),
       chemical_state={'source':'RNA/cDNA/primer/bead population; whole copied strands and mutation distribution unassembled','template':'Covalent modified RNA is not a canonical basecall; selected-path chemistry retained in prior stages','primer':'Actual preannealed RT primer carried through source preparation, occupancy unassigned'})
    b.clean('rt_wash','Actual post-RT bead wash; resuspend in 8 µl water')
    b.op('rnaseh','Add 1 µl RNase H buffer +1 µl RNase H;37°C30min',StepKind.DIGESTION,(b.catalyst('rnaseh_enzyme','Actual source 1 µl RNase H'),b.condition('rnaseh_buffer','Actual source 1 µl RNase H buffer')))
    b.clean('rnaseh_wash','Actual bead wash before 50 µl water heat elution')
    eluate=g.product('heat_eluate','Actual released cDNA-containing 50 µl heat eluate; residual RNA/protein and recovery unassigned',MaterialKind.MOLECULAR)
    residual=g.product('residual_beads','Actual residual beads after collected heat eluate',MaterialKind.MOLECULAR)
    g.react('heat_release','95°C 10 min and recover actual heat eluate',StepKind.CLEANUP,((b.x,InputRole.TEMPLATE),b.condition('heat_conditions','Actual 50 µl water / 95°C 10 min')),(eluate,residual));g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residual bead fraction separate from collected eluate; no unreported disposal'));b.x=eluate
    b.clean('cdna_cleanup','DNA Clean & Concentrator recovery; elute 10 µl',(b.material('cdna_column','Actual source DNA cleanup materials',MaterialKind.SELECTION),))
    donor=g.input('cdna_adapter','Actual 1 µl 50 µM source cDNA 3′ adapter',MaterialKind.OLIGO,state=cdna_adapter_scene(),chemical_state={'raw_source_order':CDNA_ADAPTER_WRITTEN,'modifications':'5P donor and non-nucleotide 3SpC3 block; no blocker propagated to PCR copies'})
    b.op('cdna_ligate','Actual cDNA adapter ligation 25°C overnight',StepKind.LIGATION,((donor,InputRole.STRUCTURAL),b.catalyst('rnl1','Actual 1 µl T4 RNA ligase 1'),b.material('ligation_atp','Actual 2 µl 10 mM ATP'),b.material('ligation_peg','Actual 10 µl 50% PEG8000'),b.condition('ligation_buffer','Actual 2 µl 10× T4 RNA ligase buffer')),geometry=AdapterLigationGeometry.known('single-strand',cdna_adapter_scene(),statement='Known source 5P/C3 donor; intended cDNA 3OH acceptor is an inferred constructor view. Whole cDNA/sealed junction unassembled.'))
    b.clean('ligation_cleanup','DNA Clean & Concentrator recovery; elute 21 µl',(b.material('ligation_column','Actual source DNA cleanup materials',MaterialKind.SELECTION),))
    qc=g.product('qpcr_template','Actual 1 µl of 21 µl eluate for qPCR testing',MaterialKind.QC)
    production=g.product('pcr_template','Actual 15 µl of 21 µl eluate for production PCR',MaterialKind.MOLECULAR)
    stock=g.product('untransferred','Nominal remaining 5 µl eluate; source fate unassigned',MaterialKind.MOLECULAR)
    g.react('aliquots','Actual 21 µl eluate allocations:1 µl qPCR +15 µl production +5 µl nominal stock',StepKind.SPLIT,((b.x,InputRole.TEMPLATE),),(qc,production,stock));g.mark_terminal(Terminal(stock.id,Disposition.UNRESOLVED,'Nominal untransferred source stock; no normalization or disposal'))
    qc_parts=source_supplies(b,'qpcr_supplies','Actual aggregate qPCR analytical supplies; stock packaging, identities and amounts unassigned','Actual analytical polymerase functional class',dyed=True)
    measured=g.product('qpcr_readout','Actual qPCR analytical readout; no production library input',MaterialKind.QC)
    g.react('qpcr','Source analytical qPCR testing; actual formulation unassembled',StepKind.MAGIC,((qc,InputRole.TEMPLATE),)+qc_parts,(measured,),reason='Source-positive analytical reaction; exact qPCR supply composition/primer orders unassembled, no recipe supplied by default.');g.mark_terminal(Terminal(measured.id,Disposition.QC,'Source actual qPCR analytical fraction'))
    b.x=production;supplies=source_supplies(b,'q5','Actual NEBNext Ultra II Q5 Master Mix allocation; amount unassigned','Actual copresent Q5 polymerase')
    oligos=g.input('pcr_oligo_supplies','Actual source NEBNext amplification oligo inventory; exact bodies/member choices/packaging unassembled',MaterialKind.MOLECULAR)
    pf=g.product('pcr_forward','Forward functional primer class from actual source oligo supplies',MaterialKind.OLIGO)
    pr=g.product('pcr_reverse','Reverse functional primer class from actual source oligo supplies',MaterialKind.OLIGO)
    g.react('pcr_oligo_inventory','Actual copresent functional primer classes; no extra fresh oligo dose',StepKind.GENERIC,((oligos,InputRole.TEMPLATE),),(pf,pr),changes_molecule=False)
    b.op('pcr','Actual source library amplification; cycles unassigned',StepKind.PCR,supplies+((pf,InputRole.STRUCTURAL),(pr,InputRole.STRUCTURAL)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Actual source PCR, known source supply classes; whole historical primer orders, occupied cDNA and landing/product geometry unassembled.'))
    b.clean('ampure','Actual 0.8× AMPure library purification',(b.material('ampure_beads','Actual source AMPure cleanup beads',MaterialKind.SELECTION),))
    loading=g.product('loading','Actual selected HiSeq X Ten library loading contribution; amount unassigned',MaterialKind.LIBRARY)
    library_stock=g.product('library_stock','Untransferred prepared library stock; amount/fate unassigned',MaterialKind.LIBRARY)
    g.react('loading_allocate','Actual library sequencing allocation; no inferred final DNA QC aliquot',StepKind.SPLIT,((b.x,InputRole.TEMPLATE),),(loading,library_stock));g.mark_terminal(Terminal(library_stock.id,Disposition.UNRESOLVED,'Untransferred source library stock'))
    g.sequence(loading,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':150,'Read 2':150},limitation='Source HiSeq X Ten PE150. Whole library, acquired index cycles, study index choices and read-primer landing unassembled. Modified RNA is not directly sequenced.'))
    g.declare_reviewed_scope(REVIEWED_COMMON,SOURCE)
    if control=='FTO+':g.declare_reviewed_scope(REVIEWED_FTO,SOURCE)
    if control!='input':g.declare_reviewed_scope(REVIEWED_LABEL,SOURCE)
    return g.require_publishable()

REVIEWED_COMMON=tuple("""rna oligodt rnaseh_upstream dnase polya_remove polya_remove_product rna_cleanup rna_cleanup_kept rna_cleanup_waste rna_column sonicate sonicate_product sonication_condition pnk pnk_product pnk_enzyme calibration spike spike_product rna_adapter rnl2 rna_ligate rna_ligate_product deadenylase deadenylate deadenylate_product recjf recjf_digest recjf_digest_product rt_primer primer_anneal primer_anneal_product c1 capture capture_kept capture_waste hiv_rt rt_nts rt_other rt rt_product rt_wash rt_wash_kept rt_wash_waste rnaseh rnaseh_enzyme rnaseh_buffer rnaseh_product rnaseh_wash rnaseh_wash_kept rnaseh_wash_waste heat_release heat_eluate residual_beads heat_conditions cdna_column cdna_cleanup cdna_cleanup_kept cdna_cleanup_waste cdna_adapter rnl1 ligation_atp ligation_peg ligation_buffer cdna_ligate cdna_ligate_product ligation_column ligation_cleanup ligation_cleanup_kept ligation_cleanup_waste aliquots qpcr_template pcr_template untransferred qpcr_supplies qpcr_supplies_inventory qpcr_supplies_enzyme qpcr_supplies_nts qpcr_supplies_buffer qpcr_supplies_dye qpcr qpcr_readout q5 q5_inventory q5_enzyme q5_nts q5_buffer pcr_oligo_supplies pcr_oligo_inventory pcr_forward pcr_reverse pcr pcr_product ampure ampure_kept ampure_waste ampure_beads loading_allocate loading library_stock""".split())
REVIEWED_FTO=tuple("""fto fto_supplies fto_condition demethylate demethylate_product fto_wash fto_wash_kept fto_wash_waste""".split())
REVIEWED_LABEL=tuple("""label_0_supplies label_0_inventory label_0_enzyme label_0_sam label_0_other label_0 label_0_product label_0_drain label_0_drain_kept label_0_drain_waste
label_1_supplies label_1_inventory label_1_enzyme label_1_sam label_1_other label_1 label_1_product label_1_drain label_1_drain_kept label_1_drain_waste
label_2_supplies label_2_inventory label_2_enzyme label_2_sam label_2_other label_2 label_2_product label_2_drain label_2_drain_kept label_2_drain_waste
label_3_supplies label_3_inventory label_3_enzyme label_3_sam label_3_other label_3 label_3_product label_3_drain label_3_drain_kept label_3_drain_waste
label_4_supplies label_4_inventory label_4_enzyme label_4_sam label_4_other label_4 label_4_product label_4_drain label_4_drain_kept label_4_drain_waste
label_5_supplies label_5_inventory label_5_enzyme label_5_sam label_5_other label_5 label_5_product label_5_drain label_5_drain_kept label_5_drain_waste
label_6_supplies label_6_inventory label_6_enzyme label_6_sam label_6_other label_6 label_6_product label_6_drain label_6_drain_kept label_6_drain_waste
iodine iodine_environment cyclize cyclize_product iodine_quencher quench_condition iodine_quench iodine_quench_product iodine_wash iodine_wash_kept iodine_wash_waste""".split())

def protocol_graph():return graph_for_control('FTO-')

def render_page():
    graphs=tuple(graph_for_control(c)for c in ('input','FTO-','FTO+'))
    intro=(head('m6A-SAC-seq — RNA chemistry to library preparation'),'<main class="wrap" data-protocol-complete="true"><h1>m6A-SAC-seq</h1><p>Selective allylation and iodine cyclization of RNA adenosine modifications create substrates with different HIV-RT mutational readouts. Input, FTO− and FTO+ libraries are independently prepared controls, not physically pooled here.</p><p><a href="https://doi.org/10.1038/s41587-022-01243-z">Hu et al. (2022)</a> · <a href="01_m6a_sac.html">Research note</a></p>',info('Accepted-manuscript procedure with explicitly cited publisher Table 1 adapter components. This component-level link does not assert whole final-version/manuscript equivalence. INFERRED — intended RNA/cDNA acceptor classes and selected reacted site pathways do not establish actual transcript sequences, occupancy, complete conversion or whole library geometry. No direct RNA sequencing or canonical-base rewriting is implied.'),'<h2>Known free adapter stocks</h2>','<p>RNA adapter, 5′→3′: <seq>'+RNA_ADAPTER_WRITTEN+'</seq></p><p>cDNA adapter, 5′→3′: <seq>'+CDNA_ADAPTER_WRITTEN+'</seq></p>',panel(published_adapter_inventory_scene().rows(),caption='Source Table 1 RNA adapter pattern, not a complete RNA or sequencing product; printed sample-index member is a reference stock choice, not an assigned study sample.'),panel(cdna_adapter_scene().rows(),caption='Source Table 1 free cDNA-adapter donor with 5′ phosphate and non-nucleotide 3′ C3 block; whole occupied acceptor and sealed junction remain unassembled.'),'<h2>Source values and read setup</h2><p><strong>Repeat-count interpretation.</strong> Accepted-manuscript Methods219 first describes a one-hour mixture, then a fresh twenty-minute recipe followed by “this labeling step was repeated six times”. This model interprets that as six total fresh twenty-minute rounds after the initial one-hour round (seven mixture allocations). The wording can also be read as six repeats after the separately described first twenty-minute round; that total remains unresolved here. No publication error is asserted.</p><p>The manuscript prints quench reagent Na2S2SO3; this page retains that formula without identifying a corrected reagent. The 21 µl eluate contributes 1 µl to qPCR and 15 µl to production PCR, leaving 5 µl nominal untransferred stock. HiSeq X Ten PE150 is reported; whole primer landing, acquired index cycles and final library molecules remain unassembled.</p>')
    pieces=list(intro)
    for g in graphs:pieces.extend(('<details><summary>Full '+g.name+' preparation and physical fates</summary>',WorkflowSection(g.name).render(g),'</details>'))
    pieces.append('</main>');return assemble_page(graphs[1],pieces)

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check()
    for origin,species in [('A','cyclized a6A'),('m6A','cyclized a6m6A')]:
        result=AdenosinePath(origin).allylate().cyclize();c(origin+' retains covalent chemical origin',(result.original,result.species),(origin,species))
    try:AdenosinePath('m6A').cyclize()
    except ValueError:blocked=True
    else:blocked=False
    c('Allylation required before cyclization',blocked,True)
    try:AdenosinePath('m6A').allylate().cyclize().as_canonical_rna_base()
    except ValueError:blocked=True
    else:blocked=False
    c('Covalent adduct not flattened into canonical base',blocked,True)
    c('FTO erased origin preserved',AdenosinePath('m6A').erase_selected().allylate().cyclize().species,'cyclized a6A')
    g=protocol_graph()
    c('Seven mixture classes under explicit six-round interpretation',(sum(k in g.steps for k,_ in LABEL_ROUNDS),all(g.materials[k+'_supplies'].inferred and g.materials[k+'_product'].inferred and g.materials[k+'_drain_kept'].inferred for k,_ in LABEL_ROUNDS)),(7,True))
    c('Fresh mixtures have distinct actual supply provenance',len({g.producer_of(k+'_enzyme').id for k,_ in LABEL_ROUNDS}),7)
    c('Each reaction has physical subsequent drain',all(k+'_drain' in g.steps for k,_ in LABEL_ROUNDS),True)
    c('Iodine is molecular reactant',g.materials['iodine'].kind,MaterialKind.MOLECULAR)
    c('Actual aliquots preserved',tuple(g.steps['aliquots'].outputs),('qpcr_template','pcr_template','untransferred'))
    c('Production uses actual15ul class',g.steps['pcr'].inputs[0].material,'pcr_template')
    c('RNA adapter non-nucleotide block not added to backbone',published_adapter_inventory_scene().strands['stock adapter'].mod3,'SpC3')
    c('No experimental labeling on input','label_0' in graph_for_control('input').steps,False)
    c('FTO control retains erasure and residual paths',sorted(p.fto_erased for p in graph_for_control('FTO+').materials['demethylate_product'].chemical_state.selected_reacted_paths if p.original=='m6A'),[False,True])
    for control in ('input','FTO-','FTO+'):c(control+' normal gate',graph_for_control(control).publication_problems(),[])
    source=Source(sources_root()/'m6a-sac-seq__10.1038+s41587-022-01243-z'/'s41587-022-01243-z_41587_2022_1243_MOESM1_ESM.pdf.txt','Source Table1 component inventory only')
    if have(c,source,label='Table1 source fixtures'):
        text=''.join(source.path.read_text().split());c('Whole modified RNA ordering',RNA_ADAPTER_WRITTEN in text,True);c('RT primer exact',RT_PRIMER in text,True)
        physical=RNA_ADAPTER_WRITTEN.removeprefix('rApp').replace('/iBiodT/','T').replace('/3SpC3/','')
        actual=''.join(seg.top for seg in published_adapter_inventory_scene().strands['stock adapter'].segs)
        c('Internal biotin-dT retained as one physical nucleotide',actual,physical)
        c('Whole modified cDNA ordering',CDNA_ADAPTER_WRITTEN.replace('Phos-','Phos') in text,True)
    c.report(exit_on_fail=True)
if __name__=='__main__':selftest()
