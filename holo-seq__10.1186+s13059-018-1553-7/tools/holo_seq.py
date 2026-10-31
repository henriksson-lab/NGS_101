"""Holo-seq: source-local post-library carrier-DNA digestion branches."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SourceDiscrepancy,SourceClaim
SOURCE='Xiao et al.2018 DOI10.1186/s13059-018-1553-7; published PMC6193298 Methods181–197; Additionalfile4 step-by-step library recipes; FigS1 original carrier-template image'

@dataclass(frozen=True)
class SelectedReferenceCut:
    """Conditional reference partition only; NotI stagger/end chemistry not reconstructed."""
    dna: str
    position: int
    origin: str
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or type(self.position) is not int or not 0<self.position<len(self.dna):raise ValueError('explicit DNA with internal selected reference boundary required')
        if self.origin not in {'carrier_derived','endogenous','late_dna_carrier'}:raise ValueError('declared molecular origin required')
    @property
    def reference_pieces(self):return self.dna[:self.position],self.dna[self.position:]
    @property
    def inherited_outer_handles(self):return (frozenset({'left'}),frozenset({'right'}))

def graph_for_branch(branch='mrna'):
    if branch not in {'mrna','total','small'}:raise ValueError('source mRNA/total/small branches only')
    g=ProtocolGraph('Holo-seq — '+branch+' carrier-DNA removal boundary',authority=SOURCE,representation='published_instance')
    prior=g.input('library_dna',('Post-PCR small-RNA library' if branch=='small' else 'Source adapter-ligated '+branch+' cDNA')+' containing native and carrier-derived DNA; full upstream kit architecture unassembled',K.MOLECULAR,primary=True)
    scope=['library_dna']
    def user(key,x):
        enzyme=g.input(key+'_enzyme','Source USER enzyme37°C15min; exact historical substrate/cleavage geometry unassembled',K.CATALYST)
        out=g.product(key+'_product','Source USER-treated DNA mixture; cleavage/degradation products and unreacted material not physically separated yet',K.MOLECULAR)
        g.react(key,'Source USER treatment at branch-specific position',S.DIGESTION,((x,R.TEMPLATE),(enzyme,R.CATALYST)),(out,));scope.extend([key+'_enzyme',key]);return out
    def cleanup(key,x,late=False,ratio='1×'):
        inputs=[(x,R.TEMPLATE)]
        if late:
            carrier=g.input(key+'_dna_carrier','Source100ng~2kb DNA carrier; sequence/ends/recovery unreviewed',K.MOLECULAR)
            inputs.append((carrier,R.STRUCTURAL));scope.append(key+'_dna_carrier')
        beads=g.input(key+'_beads',ratio+'AMPureXP; source physical purification',K.SELECTION)
        keep=g.product(key+'_kept','Physically recovered DNA mixture; native/remaining carrier/late DNA fractions unquantified',K.MOLECULAR)
        waste=g.product(key+'_waste','Physical bead-binding/wash rejects; exact fragment recovery unassigned',K.WASTE)
        g.react(key,'Source bead cleanup; no species disappears before separation',S.CLEANUP,tuple(inputs)+((beads,R.SELECTION),),(keep,waste));g.discard(waste,'Physical purification rejects only');scope.extend([key+'_beads',key]);return keep
    if branch=='mrna':prior=user('user_before_noti',prior)
    if branch!='small':prior=cleanup('ligation_cleanup',prior,late=True)
    enzyme=g.input('noti','NotI2µl in source CutSmart50µl,37°C2h',K.CATALYST)
    cut=g.product('noti_digest','NotI-treated DNA containing cut reference pieces, uncut molecules and unassigned carrier/native recovery; no direct RNA excision',K.MOLECULAR,chemical_state={'substrate':'library DNA','cut_products':'both physical fragments retained until cleanup','specificity_scope':'native NotI-site DNA can also be affected'})
    g.react('digest_noti','Source carrier-library cleavage; do not assert perfect selective depletion',S.DIGESTION,((prior,R.TEMPLATE),(enzyme,R.CATALYST)),(cut,));scope+=['noti','digest_noti']
    prior=cleanup('digest_cleanup',cut,late=branch=='small')
    if branch=='total':prior=user('user_after_noti',prior)
    g.mark_terminal(Terminal(prior.id,Disposition.UNRESOLVED,'Source post-removal PCR'+('/USER mixture' if branch=='total' else '')+' or small-library sequencing handoff; exact historical full molecular/USER/adapter/read chemistry unassembled. Carrier depletion efficiency and residual fractions not predicted.'))
    g.declare_reviewed_scope(scope,SOURCE)
    g.add_source_item(SourceItem('unreviewed_full_library','Historical E7490/E7530/E7420/E7300 polymers, upstream reactions, USER substrate topology and final read plan','Published main Methods177–197 and Additionalfile4 cited kit procedures'))
    if branch=='small':
        g.annotate(SourceDiscrepancy('library_dna','Published small-RNA PCR counts differ',
           (SourceClaim('Main Methods193','15cycles before NotI'),SourceClaim('Additionalfile4 small-RNA PCR block','22cycles before NotI')),
           'Model begins after PCR and assigns no cycle count.',
           'Both original published locations inspected.',
           'Applicability to individual datasets/procedure revisions unresolved; no numerical repair inferred.'))
    return g

def protocol_graph():return full_mrna_graph()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=SelectedReferenceCut('ACGTGCGGCCGCTT',6,'carrier_derived')
    c('selected reference cut preserves all DNA bases',''.join(p.reference_pieces),p.dna)
    c('cut pieces retain opposite original outer handles',p.inherited_outer_handles,(frozenset({'left'}),frozenset({'right'})))
    q=SelectedReferenceCut(p.dna,p.position,'endogenous');c('origin label does not protect a susceptible endogenous fragment',p.reference_pieces,q.reference_pieces)
    c('mRNA USER precedes carrier digest',tuple(graph_for_branch('mrna').steps)[:3],('user_before_noti','ligation_cleanup','digest_noti'))
    c('total-RNA USER follows actual carrier-removal cleanup',tuple(graph_for_branch('total').steps)[-2:],('digest_cleanup','user_after_noti'))
    c('small branch lateDNA carrier is a declared material','digest_cleanup_dna_carrier' in graph_for_branch('small').materials,True)
    original=Source(sources_root()/'holo-seq__10.1186+s13059-018-1553-7'/'s13059-018-1553-7_supp_13059_2018_1553_MOESM4_ESM.docx.txt','Additionalfile4')
    if have(c,original,label='source latecarrier and smallbranch checks'):
        text=''.join(original.path.read_text().split())
        c('source lateDNA carrier retained in recipe','Add100ng~2KbDNAcarrier' in text,True)
        c('source detailed smallRNA count kept distinct','22cyclesof' in text,True)
    full=protocol_graph()
    c('Unreviewed historical polymers and USER/read obligation remains required',full.publication_problems(),["authoritative source item 'unreviewed_full_library' is not accounted for"])
    c('NotI library is cleaned before selected PCR input',full._consumers['digest_noti_product'],['digest_cleanup'])
    c('Source end-template stock remainder is not consumed',full._consumers['end_template_remainder'],[])
    c('Selected loading stock is finished DNA library',full.materials['library_handoff_product'].kind,K.LIBRARY)
    c.report()


# Publication wiring deferred until the original required source/model scope is assembled.

def full_mrna_graph():
    from protocol_graph import SequencingPlan,SourceCoverage
    from adapter_ligation import AdapterLigationGeometry
    from pcr import PCRGeometry
    from chemdraw import Scene,Segment
    g=ProtocolGraph('Holo-seq — single-cell mRNA / NotI preparation',authority=SOURCE,representation='published_instance')
    def inp(k,label,kind=K.MOLECULAR,primary=False):return g.input(k,label,kind,primary=primary)
    def op(k,label,x,extras=(),kind=S.MAGIC,outkind=K.MOLECULAR,geometry=None):
        out=g.product(k+'_product',label+'; occupied molecules, efficiency and residual populations unassembled',outkind)
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(extras),(out,),geometry=geometry,reason='Source Additional4 selected mRNA recipe; explicit local stage, no borrowed historical kit polymer architecture.' if kind is S.MAGIC else '')
        return out
    def aliquot(k,x,label):
        take=g.product(k+'_selected',label,x.kind);rest=g.product(k+'_remainder','Untransferred source stock; amount/fate unassigned unless stated',x.kind)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(take,rest),changes_molecule=False)
        g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferred source stock, not silently consumed.'));return take
    def clean(k,x,ratio,extra=()):
        b=inp(k+'_beads',ratio+' AMPure XP and source wash/elution supply',K.SELECTION)
        out=g.product(k+'_kept','Recovered DNA after '+ratio+' cleanup; composition/size recovery unassigned',K.MOLECULAR)
        waste=g.product(k+'_removed','Physically removed bead/wash fractions',K.WASTE)
        g.react(k,'Source '+ratio+' DNA purification',S.CLEANUP,((x,R.TEMPLATE),(b,R.SELECTION))+tuple(extra),(out,waste));g.discard(waste,'Physical purification rejects.');return out
    cell=inp('single_cell','Source single-cell input; endogenous RNA populations not instantiated',K.BIOLOGICAL,True)
    carrier=inp('carrier_rna','100ng source carrier RNA mixture; IVT-derived polyA/non-polyA carrier stock, full sequence topology not adopted',K.MOLECULAR)
    lysis=inp('lysis_supplies','Source0.1%BSA/PBS or listed Triton/RNase-inhibitor lysis formulation; BSA remains physical component, no free additional dose',K.MOLECULAR)
    lysate=op('lysis','95°C5min then ice, source single-cell/carrier lysis',cell,((carrier,R.STRUCTURAL),(lysis,R.STRUCTURAL)))
    beads=inp('oligodt_beads','Actual10ul source oligo(dT)25 bead stock in25ul2xRNA Binding Buffer',K.SELECTION)
    poly=g.product('polya_selected','Bead-selected cellular and polyA carrier RNA',K.MOLECULAR)
    sup=g.product('polya_supernatant','Nonretained RNA/carrier supernatant in selected mRNA-only route',K.MOLECULAR)
    g.react('polya_capture','Source E7490 polyA selection;65°C5min,RT10min',S.BINDING,((lysate,R.TEMPLATE),(beads,R.SELECTION)),(poly,sup));g.mark_terminal(Terminal(sup.id,Disposition.UNRESOLVED,'Not followed in this selected mRNA route; dual small-RNA preparation is separate.'))
    priming=inp('fragment_mix','Actual8ul First-Strand Buffer/random-primer mix; source recipe4ul5xbuffer+1ul random primers+5ulwater preparation, selected8ul only',K.MOLECULAR)
    random= g.product('random_primers','Random-primer component of same8ul supply; whole species unassembled',K.OLIGO)
    buffer= g.product('fragment_buffer','Buffer/water component of same8ul supply',K.CONDITION)
    g.react('fragment_components','Expose actual selected fragment/primer supply',S.MAGIC,((priming,R.TEMPLATE),),(random,buffer),changes_molecule=False,reason='Two component views of one selected physical stock; no duplicate doses.')
    frag=op('fragment','Elute/fragment and prime at94°C15min,4°C hold',poly,((random,R.STRUCTURAL),(buffer,R.CONDITION)))
    first=aliquot('first_template',frag,'7.5ul fragmented/primed RNA selected from8ul source reaction')
    rt=inp('protoscript','0.5ul ProtoScript II RT',K.CATALYST)
    rri=inp('rri','0.25ul murine RNase inhibitor',K.CATALYST)
    rtwater=inp('rt_water','1.75ul water; source10ul first-strand reaction',K.CONDITION)
    # dNTP availability in historical buffer is not assigned to a new physical stock.
    cdna=op('first_strand','25°C10min,42°C60min,70°C15min first-strand cDNA synthesis',first,((first,R.STRUCTURAL),(rt,R.CATALYST),(rri,R.CATALYST),(rtwater,R.CONDITION)),S.REVERSE_TRANSCRIPTION)
    second=inp('second_mix','Actual2ul Second-Strand Synthesis Enzyme Mix; polymerase/nuclease species and packaging unassembled',K.CATALYST)
    secondbuf=inp('second_buffer','4ul10x Second-Strand Reaction Buffer+24ulwater; source40ul reaction',K.CONDITION)
    ds=op('second_strand','Second-strand synthesis16°C60min',cdna,((second,R.CATALYST),(secondbuf,R.CONDITION)))
    ds=clean('ds_cleanup',ds,'1.8x')
    ds=aliquot('end_template',ds,'27.75ul purified ds-cDNA selected for source end preparation')
    end=inp('end_mix','1.5ul NEBNext End Prep Enzyme Mix',K.CATALYST)
    eb=inp('end_buffer','3.25ul10x End Repair Reaction Buffer',K.CONDITION)
    ds=op('end_prep','Source32.5ul end preparation20°C30min then65°C30min; exact end distribution unassembled',ds,((end,R.CATALYST),(eb,R.CONDITION)),S.END_REPAIR)
    ad=inp('adapter','Actual1ul20-fold diluted NEBNext adaptor; whole polymer/endmods unassembled',K.ADAPTER)
    lig=inp('ligase_mix','Actual7.5ul Blunt/TA Ligase Master Mix; same source ligase supply, composition unresolved',K.CATALYST)
    lw=inp('ligation_water','0.75ulwater in source ligation reaction',K.CONDITION)
    unknown_adapter=Scene();unknown_adapter.strand('unknown',[Segment('unassembled historical source interface','[unassembled]',placeholder=True)],label='Molecular model limit')
    geom=AdapterLigationGeometry('unavailable','unavailable',(unknown_adapter,),'Historical source adaptor and occupied junction remain unassembled; no claim that all cited kit material omits them.')
    ds=op('adapter_ligation','Source adaptor ligation20°C15min',ds,((ad,R.STRUCTURAL),(lig,R.CATALYST),(lw,R.CONDITION)),S.LIGATION,geometry=geom)
    user=inp('user','Actual1.5ul USER enzyme; historical substrate topology unassembled',K.CATALYST)
    ds=op('user_before_noti','USER37°C15min; products remain together until cleanup',ds,((user,R.CATALYST),),S.DIGESTION)
    late=inp('late_dna_carrier','Actual100ng approximately2kb DNA carrier introduced at ligation cleanup; full sequence/end distribution unassembled')
    ds=clean('ligation_cleanup',ds,'1x',((late,R.STRUCTURAL),))
    ds=aliquot('digest_template',ds,'43ul adaptor-ligated DNA selected for50ul NotI reaction')
    noti=inp('noti','Actual2ul NotI enzyme',K.CATALYST);cutbuf=inp('cutsmart','Actual5ul10xCutSmart Buffer',K.CONDITION)
    ds=op('digest_noti','NotI37°C2h: carrier-derived and susceptible native DNA may be cut; both fragments remain until cleanup',ds,((noti,R.CATALYST),(cutbuf,R.CONDITION)),S.DIGESTION)
    ds=clean('digest_cleanup',ds,'1x');ds=aliquot('pcr_template',ds,'10ul recovered DNA selected for25ul source PCR')
    master=inp('q5_master','One actual12.5ul2x NEBNext Q5 Hot Start HiFi PCR Master Mix')
    enz=g.product('q5_enzyme','Polymerase component of same12.5ul stock',K.CATALYST);nt=g.product('q5_nt','Nucleotide component of same12.5ul stock',K.MOLECULAR);buf=g.product('q5_buffer','Buffer component of same12.5ul stock',K.CONDITION)
    g.react('q5_components','Inventory views of one actual PCR Master Mix',S.MAGIC,((master,R.TEMPLATE),),(enz,nt,buf),changes_molecule=False,reason='One physical source stock, not independent enzyme/nucleotide additions.')
    fw=inp('pcr_i7','Actual1.25ul source Index/i7 primer; whole selected order unassembled',K.OLIGO);rv=inp('pcr_i5','Actual1.25ul Universal/i5 primer; whole selected order unassembled',K.OLIGO)
    pg=PCRGeometry.unavailable(unknown_adapter,reason='Source positive indexed PCR; historical complete primer/template geometry unassembled.')
    ds=op('index_pcr','98°C30s;16cycles98°C10s/65°C75s;65°C5min source indexed PCR',ds,((fw,R.STRUCTURAL),(rv,R.STRUCTURAL),(enz,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION)),S.PCR,geometry=pg)
    ds=clean('final_cleanup',ds,'0.9x');ds=aliquot('loading',ds,'Selected finished mRNA library loading contribution; volume unassigned')
    lib=op('library_handoff','Source purified indexed ds-cDNA library for deep sequencing; no early RNA intermediate is sequenced',ds,outkind=K.LIBRARY)
    g.sequence(lib,SequencingPlan(acquisition='unspecified',limitation='Illumina-compatible indexed mRNA DNA library; selected read cycles and occupied primer geometry are unassigned in this source-scoped model. Separate degenerate-primer captured-library experiment is not imported.'))
    g.declare_reviewed_scope(('single_cell','carrier_rna','lysis_supplies','lysis','oligodt_beads','polya_capture','fragment_mix','fragment_components','fragment','first_template','protoscript','rri','rt_water','first_strand','second_mix','second_buffer','second_strand','ds_cleanup_beads','ds_cleanup','end_template','end_mix','end_buffer','end_prep','adapter','ligase_mix','ligation_water','adapter_ligation','user','user_before_noti','late_dna_carrier','ligation_cleanup_beads','ligation_cleanup','digest_template','noti','cutsmart','digest_noti','digest_cleanup_beads','digest_cleanup','pcr_template','q5_master','q5_components','pcr_i7','pcr_i5','index_pcr','final_cleanup_beads','final_cleanup','loading','library_handoff','library_handoff_product'),SOURCE)
    g.add_source_item(SourceItem('unreviewed_full_library','Historical E7490/E7530/E7420/E7300 polymers, upstream reactions, USER substrate topology and final read plan','Published main Methods177–197 and Additionalfile4 cited kit procedures'))
    return g

def render_page():
    from page import head,info
    from protocol_views import assemble_page,WorkflowSection
    return assemble_page(protocol_graph().require_publishable(),(head('Holo-seq'),'<div class="wrap" data-protocol-complete="true"><h1>Holo-seq</h1>',
      '<p class="research-notes"><a href="01_holo_seq.html">Research notes</a></p>',
      '<p>Single-cell RNA is mixed with carrier RNA, selected on oligo(dT) beads, fragmented and converted to double-stranded cDNA. Adaptor ligation and USER treatment precede NotI digestion. Cleanup, sixteen-cycle indexed PCR and final purification produce the mRNA sequencing library.</p>',
      info('Selected published mRNA / NotI route. Carrier-derived DNA is cut, not RNA directly. Native DNA containing susceptible sites may also be affected; neither perfect depletion nor fragment recovery is predicted. Historical whole adaptors, USER substrate topology and occupied read geometry remain unassembled. No canonical kit recipe is borrowed.'),
      '<h2>Oligos and readout</h2><p>The selected recipe names bead-bound oligo(dT)25, random primers, NEBNext adaptor and i7/index plus universal/i5 primers. Their whole historical orders are unassembled here. The separate captured-library experiment’s degenerate primers are not substituted. Deep sequencing is declared for the finished DNA library; read cycles remain unassigned.</p>',
      '<h2>Published branch boundaries</h2><p>The mRNA route uses sixteen PCR cycles. The separate small-RNA route prints fifteen cycles in main Methods193 and twenty-two in Additional file4. Neither small-RNA count is applied to this mRNA route. Source carrier mixture polyA-minus:polyA-plus1:20 wording is retained without inversion.</p>',
      WorkflowSection('Complete selected mRNA preparation and explicit molecular limits'),'</div>'))

if __name__=='__main__':selftest()
