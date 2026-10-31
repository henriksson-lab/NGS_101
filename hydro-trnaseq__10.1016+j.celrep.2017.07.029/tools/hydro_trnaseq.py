"""hydro-tRNAseq size enrichment, hydrolysis and source-specific end treatment."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Segment, feature
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
SOURCE='Gogakos et al.2017 doi:10.1016/j.celrep.2017.07.029, accepted-manuscript Methods122–124'
PREDECESSOR='Hafner et al.2012 doi:10.1016/j.ymeth.2012.07.030, accepted-manuscript adapter table92–148'
BOUNDARY=('Barcoded small-RNA library is explicitly inherited from Hafner2012; actual adapter '
          'member, exact fragment-selection stage and acquired read geometry remain deferred. '
          'Gated defining supplements and VOR verification remain open. CIP treatment does '
          'not establish complete resolution of every possible cyclic/2′/3′ phosphate species.')

def barcoded_adapter_inventory():
    """Unassigned member of the cited twenty-adapter family, not fabricated Ad01 use."""
    return Construct([Segment('unassigned source pentamer barcode','N'*5,'cbc',placeholder=True,
                       feature=feature('hydro_sample','sample_index','whitelist',
                                       whitelist='Hafner2012 twenty adapters; actual member unassigned')),
                      Segment('common3′adapterDNA','TCGTATGCCGTCTTCTGCTTG')],
                     name='inherited-family inventory only;5′rApp/3′aminohexyl; actual member unassigned')

def graph_for_replicate(replicate):
    if replicate not in (1,2,3,4): raise ValueError("source replicates are1–4")
    duration=10 if replicate==1 else 60
    g=ProtocolGraph('hydro-tRNAseq — size enrichment and hydrolysis',authority=SOURCE,
                    representation='published_instance')
    total=g.input('total_rna','HEK293 total RNA after TRIzol extraction',MaterialKind.MOLECULAR,primary=True)
    gel=g.input('urea_page','12%urea-PAGE;60–100nt window',MaterialKind.SELECTION)
    enriched=g.product('enriched_rna','60–100nt RNA including mature and precursor tRNAs',MaterialKind.MOLECULAR)
    rejected=g.product('outside_window','RNA outside recovered60–100nt window',MaterialKind.WASTE)
    g.react('size_enrich','Resolve20µg totalRNA; recover60–100nt RNA',StepKind.SELECTION,
            ((total,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(enriched,rejected))
    g.discard(rejected,'Physical PAGE fraction outside the collected window')
    carbonate=g.input('carbonate','10mMNa2CO3+10mMNaHCO3,pH9.7,65°C',MaterialKind.CONDITION)
    mixture=g.product('hydrolysis_mixture','all hydrolysis fragments and any deacylation-derived molecular products',MaterialKind.MOLECULAR)
    g.react('hydrolyze',f'Limited alkaline hydrolysis:{duration}min, source replicate{replicate}',
            StepKind.FRAGMENTATION,((enriched,InputRole.TEMPLATE),(carbonate,InputRole.CONDITION)),(mixture,))
    cip=g.input('cip','calf intestinal phosphatase37°C1h',MaterialKind.CATALYST)
    treated=g.product('cip_treated','CIP-treated fragment mixture; exact residual cyclic/2′/3′phosphate distribution unassigned',MaterialKind.MOLECULAR)
    g.react('dephosphorylate','Source CIP dephosphorylation of hydrolysis products',
            StepKind.DEPHOSPHORYLATION,((mixture,InputRole.TEMPLATE),(cip,InputRole.CATALYST)),(treated,))
    prepared=g.product('pnk_input','RNA contribution entering source20µlPNK reaction; preceding50µlCIP transfer/preparation interface unassembled',MaterialKind.MOLECULAR)
    retained=g.product('cip_interface_residual','Other CIP-stage material, if any; allocation/partition/fate unassigned',MaterialKind.MOLECULAR)
    g.react('cip_pnk_interface','Source50µlCIP→20µlPNK preparation interface, not an inferred whole-stock transfer',StepKind.MAGIC,((treated,InputRole.TEMPLATE),),(prepared,retained),reason='Source names50µlCIP and20µlPNK reactions without assigning the intervening transfer/reconstitution here. No new cleanup/aliquot or negative volume is invented.')
    g.mark_terminal(Terminal(retained.id,Disposition.UNRESOLVED,'Unassembled source preparation interface, not known discard.'))
    pnk=g.input('pnk','T4PNK37°C1h',MaterialKind.CATALYST)
    atp=g.input('atp','Source1mMATP molecular phosphate-donor supply; stockvolume unassigned',MaterialKind.MOLECULAR)
    repaired=g.product('rephosphorylated','PNK-treated RNA-fragment mixture; no undeclared physical debris removal',MaterialKind.MOLECULAR)
    g.react('rephosphorylate','Rephosphorylate with T4PNK and ATP',StepKind.PHOSPHORYLATION,
            ((prepared,InputRole.TEMPLATE),(pnk,InputRole.CATALYST),(atp,InputRole.STRUCTURAL)),(repaired,))
    inherited_library_trace(g,repaired)
    g.declare_reviewed_scope(('cip_pnk_interface',),SOURCE)
    g.declare_reviewed_scope(('total_rna','urea_page','size_enrich','carbonate','hydrolyze',
                             'cip','dephosphorylate','pnk','atp','rephosphorylate'),SOURCE)
    return g.require_publishable()

def protocol_graph():
    """Default source instance: replicate1,10min alkaline hydrolysis."""
    return graph_for_replicate(1)

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('molecular hydrolysis products remain before physical separation',
      g.materials['hydrolysis_mixture'].kind,MaterialKind.MOLECULAR)
    c('actual sample identifier remains unassigned',
      list(barcoded_adapter_inventory())[0].feature.role,'sample_index')
    c('source-supported DNA finish passes normal gate',g.publication_problems(),[])
    ports=lambda i:tuple(p.material for p in g.steps[i].inputs)
    c('sequencing material is actual produced DNA',g.materials['library_loaded'].chemical_state['polymer'],'DNA')
    c('DNA produced by typedPCR',g.steps['library_pcr'].kind,StepKind.PCR)
    c('DNA produced after actual reverse transcription',g.steps['library_rt'].kind,StepKind.REVERSE_TRANSCRIPTION)
    c('RT carries source primer from elution','retained_carrier' in ports('library_rt'),True)
    c('RT does not re-add carrier stock','reverse_carrier' in ports('library_rt'),False)
    c('PCR has actual molecular nucleotide supplies','library_pcr_nucleotides' in ports('library_pcr'),True)
    c('untransferred cDNA remains accounted',g.terminals['pcr_template_rest'].disposition,Disposition.UNRESOLVED)
    c('ancestral polymer applicability marked inferred',all(g.materials[i].inferred for i in ('three_adapter','five_adapter','reverse_carrier','pcr_forward','pcr_reverse')),True)
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/LIBRARY_FOLDER/LIBRARY_FILE,LIBRARY_SOURCE)
    if have(c,source,label='Inherited library original'):
        text=''.join(source.path.read_text().split())
        c('whole ancestral forward ordering',FORWARD in text,True)
        c('whole ancestral reverse ordering',REVERSE in text,True)
        c('whole ancestral RNA5 ordering',(RA5 if BARCODED else ''.join('r'+b for b in RA5)) in text,True)
        c('whole ancestral activated blockedDNA3 ordering',('rAppTCACT'+COMMON3+'-L' if BARCODED else 'App'+COMMON3+'-L') in text,True)
    c('PNK ATP molecular donor',g.materials['atp'].kind,MaterialKind.MOLECULAR)
    c.report()



LIBRARY_FOLDER='doi_10.1016+j.ymeth.2012.07.030'
LIBRARY_FILE='j.ymeth.2012.07.030_PMC3508525.xml.txt'
BARCODED=True
PAGE_FILE='hydro_trnaseq.html'
NOTE_FILE='01_hydro_trnaseq.html'
PAGE_TITLE='hydro-tRNAseq'
LIBRARY_SOURCE=PREDECESSOR+'; common-core barcoded smallRNA ligation/RT/PCR inheritance explicitly cited2017Methods124; individual-study choices unresolved'
READ_LIMIT='Source2017 reports IlluminaHiSeq2500; acquired cycles, read-end count and exact sequencing primer landing remain unassigned. Do not import the2012 guide run setup as a2017 measurement.'
PAGE_SOURCE='<a href="https://doi.org/10.1016/j.celrep.2017.07.029">Gogakos et al.2017</a>, associatedAMPMC5564215, and explicitly cited Hafner2012 library method, associatedAMPMC3508525.'
PAGE_SUMMARY='Size-enrich60–100nt HEK293 RNA, partially hydrolyze it, treat ends with CIP and ATP/PNK, then produce the source-reported barcoded small-RNA cDNA library. The default graph is source replicate1(10min); replicates2–4 use60min.'
PAGE_LIMITS='The50µlCIP→20µlPNK preparation interface is unassembled, not a new invented cleanup. CIP is not proof that every cyclic/2′/3′phosphate has become ligation-compatible. The19–35nt fragment class does not locate a distinct preligationPAGE step: cited-guide adapter-bearing-product size selection is retained at a qualified recovery boundary. Actual barcode/member, calibrator/PmeI and pool choices remain unassigned. Gated defining supplements and finalVOR comparison remain open; no source-absence claim.'

# Source-owned common-core inheritance only. No study barcode/variant/run defaults.
# References are inferred by constructor: target-study exact applicability remains open.
import illumina as il
from chemdraw import Scene,revcomp,panel
from protocol_model import role
from protocol_graph import MaterialKind as K,InputRole as R,StepKind as S,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_views import assemble_page,WorkflowSection
from page import head,info
RA5=il.TRUSEQ_SMALL_RNA_RA5.replace('T','U')
REVERSE=il.P7[:21]  # Whole historical source primer verified, not modern full P7.
COMMON3=revcomp(REVERSE)
FORWARD='AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA'

def reagent_reference(name,body,*,inferred=True):
    return Construct([Segment(name,body,inferred=inferred)],name='INFERRED inherited source reference: '+name)

def library_end_scene(three=True):
    sc=Scene()
    if three:
        body=three_reference()
        sc.strand('adapter',list(body),label='INFERRED inherited 3DNAadapter reference; actual member/variant unassigned',mod5='App',mod3='L (blocked)')
        sc.strand('rna',[role('eligible RNA3OH substrate class; residual ends unassigned')],label='INFERRED eligible class only',mod3='OH')
        sc.ligation_end('adapter',list(body)[0].name,kind='single-strand',polarity='5-prime',owner='adapter')
        sc.ligation_end('rna','eligible RNA3OH substrate class; residual ends unassigned',kind='single-strand',polarity='3-prime',owner='substrate')
    else:
        sc.strand('adapter',list(reagent_reference('RNA5 adapter',RA5)),label='INFERRED ancestral RNA5adapter',mod3='OH')
        sc.strand('rna',[role('adapter-bearing eligible RNA5P class')],label='INFERRED eligible class only',mod5='P')
        sc.ligation_end('adapter','RNA5 adapter',kind='single-strand',polarity='3-prime',owner='adapter')
        sc.ligation_end('rna','adapter-bearing eligible RNA5P class',kind='single-strand',polarity='5-prime',owner='substrate')
    return sc

def three_reference(*,inferred=True):
    parts=[]
    if BARCODED:
        parts.append(Segment('unassigned source pentamer','N'*5,placeholder=True,inferred=inferred,feature=feature('hydro_sample','sample_index','whitelist',whitelist='Hafner2012 twenty adapters; actual member unassigned')))
    parts.append(Segment('ancestral adapterDNA body',COMMON3,inferred=inferred))
    return Construct(parts,name='INFERRED ancestral adapter family, not assigned actual study member')

def inherited_library_trace(g,source_rna):
    """Actual source-positive DNA producer; inherited core, not whole-recipe identity."""
    def mixed(label):
        sc=Scene();sc.strand('dna',[role('actual RNA-derived DNA component; occupied sequence unassembled')],label=label)
        sc.strand('other',[role('remaining RNA/reagent or digestion material; physical partition unassigned')],label='same tube until separation')
        return sc
    def recover(i,x,label,carrier=None):
        c=g.input(i+'_conditions',label,K.CONDITION)
        y=g.product(i+'_rna','Recovered source adapter-bearing RNA class; yield/carryover unassigned',K.MOLECULAR)
        waste=g.product(i+'_removed','Source physical recovery rejects; composition unassigned',K.WASTE)
        ports=[(x,R.TEMPLATE),(c,R.CONDITION)]
        if carrier:ports.append((carrier,R.STRUCTURAL))
        outputs=[y,waste]
        if carrier:
            retained=g.product('retained_carrier','Same-tube retained reverse-primer carrier for RT, not new dose',K.OLIGO,state=carrier.state,inferred=True)
            outputs.append(retained)
        g.react(i,label,S.CLEANUP,ports,outputs);g.discard(waste,'Physical source-guide recovery; no full contaminant removal or study-specific schedule inferred.')
        return y
    adapter=g.input('three_adapter','INFERRED inherited preadenylated blockedDNA3adapter reference; study member/variant unassigned',K.ADAPTER,state=three_reference(),chemical_state={'5prime':'source App activating group','3prime':'source L blocker; 2012 aminohexyl, 2008 aminolinker/blocker'},inferred=True)
    ligase3=g.input('ligase3',('Source2012 Rnl2(1–249)K227Q reference; target-study exact enzyme/dose applicability unassigned' if BARCODED else 'Source2008 truncatedRnl2 orK227Q alternatives; target-study choice unassigned'),K.CATALYST)
    buffer3=g.input('ligase3_buffer','Source-guide ATP-free ligation buffer/DMSO; source-study stock volumes unassigned',K.CONDITION)
    a3=g.product('rna_three_adapter','RNA with source3adapter-bearing eligible class; all occupied products unassembled',K.MOLECULAR)
    g.react('ligate_three','INFERRED reagent reference — source-inherited RNA3adapter ligation',S.LIGATION,((source_rna,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase3,R.CATALYST),(buffer3,R.CONDITION)),(a3,),geometry=AdapterLigationGeometry.known('single-strand',library_end_scene(),statement='Ancestral AppDNA5donor/blocked3end and eligible RNA3OH class; full target-study polymer/member/end distribution unresolved.'))
    a3=recover('recover_three',a3,'Source inherited adapter-product gel/ethanol recovery; actual barcode-pool composition and19–35nt selection placement unassigned')
    five=g.input('five_adapter','INFERRED ancestral RNA5adapter reference; no modern tail extension',K.ADAPTER,state=reagent_reference('RNA5 adapter',RA5),inferred=True)
    ligase5=g.input('ligase5','Source-guide T4RNA ligase1; target-study dose/variant applicability unassigned',K.CATALYST)
    atp5=g.input('five_atp','Actual ATP-containing molecular donor supply for inherited RNA5ligation; study dose unassigned',K.MOLECULAR)
    buffer5=g.input('five_buffer','Source-guide RNA5ligation buffer/DMSO; stock formulation unassigned here',K.CONDITION)
    a5=g.product('rna_two_adapters','Two-adapter-bearing sourceRNA class, not yet DNA',K.MOLECULAR)
    g.react('ligate_five','INFERRED reagent reference — source-inherited RNA5adapter ligation',S.LIGATION,((a3,R.TEMPLATE),(five,R.STRUCTURAL),(ligase5,R.CATALYST),(atp5,R.STRUCTURAL),(buffer5,R.CONDITION)),(a5,),geometry=AdapterLigationGeometry.known('single-strand',library_end_scene(False),statement='Ancestral RNA5adapter3OH accepts eligible RNA5P; exact target-study occupied junction unassembled.'))
    carrier=g.input('reverse_carrier','INFERRED inherited reverse-PCR-primer carrier reference added during RNA-product gel elution; actual study amount unassigned',K.OLIGO,state=reagent_reference('reverse primer',REVERSE),inferred=True)
    recovered=recover('recover_five',a5,'Source-guide RNA5ligation-product gel elution/ethanol recovery with reverse primer carrier; retained components same tube',carrier)
    rt=g.input('library_rt_enzyme','Source-inherited reverse-transcriptase supplies; cited guides use SuperScriptIII, exact target-study dose/variant unassigned',K.CATALYST)
    nt=g.input('rt_nucleotides','Actual nucleotide-containing molecular supplies for inherited RNA→DNA synthesis; study formulation/dose unassigned',K.MOLECULAR)
    cdna=g.product('library_cdna','Actual source RT-produced DNA component plus RNA template/reagent mixture; occupied DNA sequence unassembled',K.MOLECULAR,state=mixed('RNA-derived cDNA, not RNA sequenced directly'),chemical_state={'polymer':'DNA component in RNA/cDNA mixture','source_production':'RNA→DNA reverse transcription'})
    g.react('library_rt','Actual source-inherited reverse transcription produces cDNA; INFERRED carrier-primer reference',S.REVERSE_TRANSCRIPTION,((recovered,R.TEMPLATE),(g.materials['retained_carrier'],R.STRUCTURAL),(rt,R.CATALYST),(nt,R.STRUCTURAL)),(cdna,))
    koh=g.input('library_koh','Source-guide KOH/Tris RNA-template hydrolysis afterRT; study volumes unassigned',K.CONDITION)
    hydro=g.product('library_hydrolysate','cDNA plus RNA hydrolysis products; fragments not automatically discarded',K.MOLECULAR,state=mixed('DNA survives in template-hydrolysis mixture'),chemical_state={'polymer':'DNA component plus RNA hydrolysate'})
    g.react('library_rna_hydrolysis','Source-inherited postRT RNA hydrolysis; keep digestion products in mixture',S.CHEMICAL_CONVERSION,((cdna,R.TEMPLATE),(koh,R.CONDITION)),(hydro,))
    acid=g.input('library_neutralization','Source-guide HCl neutralization/pH requirement before PCR; exact study implementation unassigned',K.CONDITION)
    neutral=g.product('neutral_cdna','Neutralized source cDNA-containing mixture',K.MOLECULAR,state=hydro.state,chemical_state=hydro.chemical_state)
    g.react('library_neutralize','Source-inherited neutralization before DNA amplification',S.QUENCHING,((hydro,R.TEMPLATE),(acid,R.CONDITION)),(neutral,))
    aliquot=g.product('pcr_template','Selected actual source cDNA PCR contribution; guide10µl, study transfer volume unassigned',K.MOLECULAR,state=neutral.state,chemical_state=neutral.chemical_state)
    rest=g.product('pcr_template_rest','Untransferred neutralized cDNA stock, if any; amount/future use unassigned',K.MOLECULAR,state=neutral.state)
    g.react('library_pcr_transfer','Source cDNA template contribution and remainder',S.SPLIT,((neutral,R.TEMPLATE),),(aliquot,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No whole-stock consumption or invented remaining amount.'))
    f=g.input('pcr_forward','INFERRED historical forwardPCR primer reference, no modernP5 extension',K.OLIGO,state=reagent_reference('historical forward',FORWARD),inferred=True)
    r=g.input('pcr_reverse','INFERRED fresh source PCR reverse-primer reference; separate from earlier carrier dose',K.OLIGO,state=reagent_reference('historical reverse',REVERSE),inferred=True)
    pol=g.input('library_pcr_enzyme','Source-inherited PCR catalyst supplies; guides nameTaq, target-study exact dose/cycle program unassigned',K.CATALYST)
    pcrnt=g.input('library_pcr_nucleotides','Actual molecular nucleotide supplies for source DNA amplification; study doses unassigned',K.MOLECULAR)
    amplified=g.product('library_amplified','Actual PCR-produced source DNA library class; full amplicons/side products unassembled',K.MOLECULAR,state=Construct([role('source PCR DNA, not an RNA endpoint; occupied library unassembled')]),chemical_state={'polymer':'DNA','source_production':'RNA→RTcDNA→PCRDNA'})
    g.react('library_pcr','Actual inherited DNA PCR; INFERRED historical primer references, exact cycle/whole products unassigned',S.PCR,((aliquot,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST),(pcrnt,R.STRUCTURAL)),(amplified,),geometry=PCRGeometry.unavailable(aliquot.state,reason='Exact target-study adapter/primer variants, occupied template and binding/product geometry remain unassembled; real source PCR DNA production is declared.'))
    finish=g.input('library_finish_conditions','Source-inherited library recovery/finishing interface; calibrator/PmeI/cycleQC/NTC and actual pool plan not assigned to target study',K.CONDITION)
    library=g.product('library','Actual source sequenced DNA-library class after inherited finishing; full physical components unresolved',K.LIBRARY,state=amplified.state,chemical_state=amplified.chemical_state)
    residue=g.product('library_finish_residual','Residual/untransferred finishing material, if any; partition/fate unassigned',K.MOLECULAR)
    g.react('library_finish','Source-positive DNA-library finishing; study-specific marker/control/pooling/size schedule unresolved',S.MAGIC,((amplified,R.TEMPLATE),(finish,R.CONDITION)),(library,residue),reason='Defining article explicitly converts RNA into cDNA libraries using this cited method and sequences them. Common-core RT/PCR DNA production is explicit; source-study barcode/adapter alternatives, calibrator and exact finishing/pooling schedule are not invented.')
    g.mark_terminal(Terminal(residue.id,Disposition.UNRESOLVED,'No complete removal, exact yield or fictional side-product amount asserted.'))
    loaded=g.product('library_loaded','Selected actual sequencing-library DNA contribution; loading amount unassigned',K.LIBRARY,state=library.state,chemical_state=library.chemical_state)
    remaining=g.product('library_remaining','Untransferred source DNA library stock, if any; amount/future use unassigned',K.LIBRARY,state=library.state)
    g.react('library_loading','Source DNA sequencing contribution and stock remainder',S.SPLIT,((library,R.TEMPLATE),),(loaded,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Actual loading quantities unassigned.'))
    g.sequence(loaded,SequencingPlan(limitation=READ_LIMIT,acquisition='unspecified'))
    g.declare_reviewed_scope(LIBRARY_SCOPE,LIBRARY_SOURCE)

LIBRARY_SCOPE=('three_adapter','ligase3','ligase3_buffer','ligate_three','recover_three_conditions','recover_three','five_adapter','ligase5','five_atp','five_buffer','ligate_five','reverse_carrier','recover_five_conditions','recover_five','library_rt','library_rt_enzyme','rt_nucleotides','library_cdna','library_koh','library_rna_hydrolysis','library_neutralization','library_neutralize','library_pcr_transfer','pcr_forward','pcr_reverse','library_pcr_enzyme','library_pcr_nucleotides','library_pcr','library_finish_conditions','library_finish','library','library_loading','library_loaded','library_remaining')

def render_page():
    g=protocol_graph().require_publishable()
    return assemble_page(g,(head(PAGE_TITLE),'<div class="wrap" data-protocol-complete="true"><h1>'+PAGE_TITLE+'</h1>',info(PAGE_SOURCE),'<p>'+PAGE_SUMMARY+'</p>',
      info('INFERRED — exact target-study adapter members/variants and historical primer applicability remain unassigned. Displayed ancestral reagent polymers are source-verified references, inferred by constructor, not established individual-study orders. Whole occupied RNA/DNA strands, yield and read geometry remain unassembled. All inspected defining/predecessor sources are author manuscripts; full final-version equivalence is not asserted.'),
      '<h2>Inherited component references</h2><p>INFERRED source-family 3DNA adapter: '+('unassigned pentamer + ' if BARCODED else '')+COMMON3+'; 5′App and terminal L block retained. '+('The family has twenty pentamer sample identifiers; no Ad01 is assigned to this study. ' if BARCODED else 'Specific 2010 adapter/ligase alternatives remain unassigned. ')+'</p>',
      '<p>RNA5 adapter: '+RA5+'<br>Historical PCR forward: '+FORWARD+'<br>Historical PCR reverse / source carrier: '+REVERSE+'</p>',
      '<h2>Actual DNA production and read setup</h2><p>The source explicitly reports cDNA-library generation using the cited small-RNA method. The graph includes source-inherited adapter ligations, reverse transcription, RNA-template hydrolysis, DNA PCR and finishing before its DNA sequencing terminal. It never sequences the RNA front end directly. '+READ_LIMIT+'</p>',
      info(PAGE_LIMITS),'<p>The predecessor guides prescribe reverse-primer carrier retention, cycle titration and no-template PCR controls. Target-study calibrator/marker use, exact control schedule and pooling design remain unassigned; no new wet-control samples or fixed all-sample pool are invented.</p>',
      '<details><summary>Full source-qualified preparation and unresolved molecular checks</summary>',WorkflowSection('RNA front end → source RT/PCR DNA library'),'</details>',
      '<p><a href="'+NOTE_FILE+'">Source notes and remaining questions</a></p></div>'))

if __name__=='__main__':selftest()
