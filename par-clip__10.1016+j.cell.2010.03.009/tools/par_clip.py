"""Original PAR-CLIP front end; no guessed library adapters or mutation chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Segment,feature
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
SOURCE = 'Hafner et al. 2010 doi:10.1016/j.cell.2010.03.009, Methods PAR-CLIP (JATS205)'
BOUNDARY = ('Original chemistry supplement is gated. Library preparation explicitly cites '
            'Hafner2008 (10.1016/j.ymeth.2007.09.009), whose Solexa route is read, but '
            'source-specific adapter/variant applicability and downstream molecular model '
            'remain deferred. No sequencing-ready terminal is declared.')


def rna(name):
    return Construct([Segment(name, '[variable 4SU-labelled RNA]', placeholder=True)], name=name)


def protocol_graph():
    g = ProtocolGraph('PAR-CLIP — 4SU crosslink and RNA recovery', authority=SOURCE,
                      representation='published_instance')
    sample = g.input('labelled_cells', 'HEK293 cells grown overnight with 100 µM 4SU',
                     MaterialKind.BIOLOGICAL, primary=True)
    uv = g.input('uv365', '365 nm UV irradiation', MaterialKind.CONDITION)
    crosslinked = g.product('crosslinked_rnp', '4SU-containing RNA crosslinked to RBP',
                            MaterialKind.MOLECULAR, state=rna('crosslinked RNA'))
    g.react('crosslink', 'Irradiate live 4SU-labelled cells at 365 nm', StepKind.CROSS_LINKING,
            ((sample, InputRole.TEMPLATE), (uv, InputRole.CONDITION)), (crosslinked,))
    nt1 = g.input('rnase_t1', 'RNase T1', MaterialKind.CATALYST)
    trimmed = g.product('trimmed_rnp', 'RNase-T1-trimmed crosslinked RNA–protein complexes',
                        MaterialKind.MOLECULAR, state=rna('trimmed RNA'))
    g.react('digest_rna', 'NP40 lysis and RNase T1 digestion', StepKind.DIGESTION,
            ((crosslinked, InputRole.TEMPLATE), (nt1, InputRole.CATALYST)), (trimmed,))
    antibody = g.input('anti_flag', 'anti-FLAG Protein G Dynabeads', MaterialKind.SELECTION)
    nt2=g.input('rnase_t1_onbead','Fresh source RNaseT1 added ON immunoprecipitate; dose unassigned',MaterialKind.CATALYST)
    ip = g.product('ip_rnp', 'immunoprecipitated FLAG/HA-tagged RBP–RNA complexes',
                   MaterialKind.MOLECULAR, state=rna('immunopurified RNA'))
    reject = g.product('ip_waste', 'unretained lysate and wash material', MaterialKind.WASTE)
    g.react('immunoprecipitate', 'Anti-FLAG immunoprecipitation and on-bead RNase T1 digestion',
            StepKind.AFFINITY_CAPTURE, ((trimmed, InputRole.TEMPLATE),
                                      (antibody, InputRole.SELECTION), (nt2, InputRole.CATALYST)), (ip, reject))
    g.discard(reject, 'Immunoprecipitation flow-through and wash fractions')
    cip = g.input('cip', 'calf intestinal alkaline phosphatase', MaterialKind.CATALYST)
    dephos = g.product('dephosphorylated', 'dephosphorylated complex-associated RNA',
                      MaterialKind.MOLECULAR, state=rna('dephosphorylated RNA'))
    g.react('dephosphorylate', 'CIP dephosphorylation on beads', StepKind.DEPHOSPHORYLATION,
            ((ip, InputRole.TEMPLATE), (cip, InputRole.CATALYST)), (dephos,))
    pnk = g.input('pnk', 'polynucleotide kinase', MaterialKind.CATALYST)
    atp = g.input('radioactive_atp','Actual radioactive ATP molecular phosphate donor',MaterialKind.MOLECULAR)
    labelled = g.product('labelled_rnp', 'radiolabelled crosslinked RNA–protein complex',
                         MaterialKind.MOLECULAR, state=rna('radiolabelled RNA'))
    g.react('radiolabel', 'PNK and radioactive ATP label complex-associated RNA', StepKind.PHOSPHORYLATION,
            ((dephos, InputRole.TEMPLATE), (pnk, InputRole.CATALYST), (atp, InputRole.STRUCTURAL)), (labelled,))
    pk = g.input('proteinase_k', 'proteinase K', MaterialKind.CATALYST)
    recovered = g.product('recovered_rna', 'RNA recovered after gel/electroelution and proteinase K extraction',
                           MaterialKind.MOLECULAR, state=rna('recovered RNA'))
    g.react('recover', 'SDS-PAGE, electroelution, proteinase K, acidic phenol/chloroform and precipitation',
            StepKind.DEPROTEINIZATION, ((labelled, InputRole.TEMPLATE), (pk, InputRole.CATALYST)), (recovered,))
    inherited_library_trace(g,recovered)
    g.declare_reviewed_scope(('rnase_t1_onbead',),SOURCE)
    g.declare_reviewed_scope(('labelled_cells','uv365','crosslink','rnase_t1','digest_rna','anti_flag',
                             'immunoprecipitate','cip','dephosphorylate','pnk','radioactive_atp',
                             'radiolabel','proteinase_k','recover'), SOURCE)
    return g.require_publishable()


def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('RNA recovery is consumed by actual downstream library generation','recovered_rna' in tuple(p.material for p in g.steps['ligate_three'].inputs),True)
    c('source-qualified DNA finish passes normal gate',g.publication_problems(),[])
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
    c('fresh source on-bead RNaseT1','rnase_t1_onbead' in ports('immunoprecipitate'),True)
    c('old lysate RNaseT1 not re-added','rnase_t1' in ports('immunoprecipitate'),False)
    c('radiolabel ATP molecular donor',g.materials['radioactive_atp'].kind,MaterialKind.MOLECULAR)
    c.report()



LIBRARY_FOLDER='doi_10.1016+j.ymeth.2007.09.009'
LIBRARY_FILE='j.ymeth.2007.09.009_PMC2847350.html.txt'
BARCODED=False
PAGE_FILE='par_clip.html'
NOTE_FILE='01_par_clip.html'
PAGE_TITLE='PAR-CLIP'
LIBRARY_SOURCE='Hafner2008 doi:10.1016/j.ymeth.2007.09.009, associatedAMPMC2847350 Solexa table/ligation/RT/PCR; explicitly cited2010Methods205, exact2010 variants unresolved'
READ_LIMIT='Defining source reports Solexa sequencing. The32nt representative sequence-context and40nt mutation-control context are not a universal acquired run length; exact cycles/read-end count, run primer and index/UMI architecture remain unassigned.'
PAGE_SOURCE='<a href="https://doi.org/10.1016/j.cell.2010.03.009">Hafner et al.2010</a>, associatedAMPMC2861495, and explicitly cited Hafner2008 Solexa library method, associatedAMPMC2847350.'
PAGE_SUMMARY='365nm crosslinking of4SU-labelled HEK293 RNA to tagged RBP, RNase/IP and radiolabelling recover RNA for the explicitly cited small-RNA cDNA library route. Source-inherited RT/PCR now produces DNA before sequencing.'
PAGE_LIMITS='4SUoccupancy/crosslinking and backgroundRNA are not uniform. The T→C sequencing signature is not a literal conversion of every RNAU toC. The6SG/G→A route is separate. Source QKI oligo mutation controls are distinct from the selected cellular-RBP preparation, not invented paired-control lysates. Exact source2010 adapter/enzyme variants, defining chemistry supplement and tables remain unreconciled; cited2008 WT/K227Qalternatives are not silently resolved. No finalVOR or full-component equivalence/absence claim.'

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
