"""Protocol-owned source-limited model; no finished-publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, Construct, feature, revcomp, oligo, panel, workflow_from_protocol_graph, workflow_panel
from protocol_model import role, molecule
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, Terminal, Disposition, SourceItem, SequencingPlan
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry
from page import head, info

def mixture(*species):
    sc = Scene()
    for n, (label, parts) in enumerate(species):
        sc.strand('s' + str(n), list(parts), label=label)
    return sc

def cleanup(g, ident, source, template, state):
    reagent = g.input(ident + '_condition', source, K.CONDITION)
    out = g.product(ident, 'recovered source preparation', K.LIBRARY if template.kind is K.LIBRARY else K.MOLECULAR, state=state)
    waste = g.product(ident + '_waste', 'physically separated workup fractions; exact molecular partition unquantified', K.WASTE)
    g.react(ident + '_step', source, S.CLEANUP, ((template, R.TEMPLATE), (reagent, R.CONDITION)), (out, waste))
    g.discard(waste, 'Actual source physical separation; no perfect recovery/depletion asserted.')
    return out

def end_boundary(label):
    sc = molecule('unassembled ligation-ready geometry', (role(label),))
    return AdapterLigationGeometry('unavailable', 'unavailable', (sc,), label + '; no source-absence claim.')

def render_page():
    import html
    g = protocol_graph()
    parts = [head(TITLE), '<main class="wrap"><h1>' + html.escape(TITLE) + '</h1>',
             '<p><a href="' + NOTES + '">Research notes</a></p>', info(SOURCE),
             info('Private incomplete model. ' + CAVEAT)]
    if callable(globals().get('oligos')):
        parts.extend(['<h2>Source ordering oligos</h2>', *oligos()])
    for label, terminal in PREVIEWS:
        parts.extend(['<h2>' + html.escape(label) + '</h2>', workflow_panel(workflow_from_protocol_graph(g, terminal))])
    parts.extend([info(READOUT), '</main>'])
    return '\n'.join(parts)

def pop(label):
    return mixture(('source DNA',(role(label),)),('other physical material',(role('other original DNA/RNA/protein or reaction remnants; exact partition and recovery unquantified'),)))
def op(g,key,label,kind,x,out_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION,additional=()):
    reagent=g.input(key+'_stock',reagent_label,reagent_kind)
    product=g.product(key+'_product',label,K.MOLECULAR,state=out_state)
    kw={'reason':'Source operation known; exact molecular substrate/product relation and mixture fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(key,label,kind,((x,R.TEMPLATE),(reagent,reagent_role))+tuple(additional),(product,),**kw)
    return product


def native(g,x,limitation):
    # Source acquisition may be stated without fabricated sequencing primers/cycle counts.
    g.sequence(x,SequencingPlan(limitation=limitation))

def pending(g,label):
    g.add_source_item(SourceItem('unreviewed_assembly',label,'Named source operations are represented; historical end/kit/primer assembly and complete original-version applicability unreviewed.'))



PAGE_FILE='pu-seq.html'
TITLE='Pu-seq'
NOTES='01_pu-seq.html'
SOURCE='Associated manuscript EMS61634 Methods115 Libraryproductionandsequencing/Fig1d–e; publishedSI MOESM16 identity/captions scope; DOI10.1038/nsmb.2962.'
CAVEAT='Original Pu-seq alkali-fragmentation plus dUTP-complement/delayedUSER route. Full final-VOR applicability, source BioPrime/NEBNextUltra polymer/end geometry remain unreviewed. Source native rNMP boundary differs from canonical readout; no perfect selection/strand erasure/copying or biological polymerase assignment guaranteed.'
READOUT='Source100bp paired-end HiSeq2500 or150bp paired-end NextSeq500. Source mapsR1 start to downstream5′end generatedbyalkali; the embedded upstream rNMP and polymerase usage are computationally assigned, not first sequenced nucleotide identity. No read length/priming geometry fabricated beyond these source acquisition branches.'
PREVIEWS=(('Original alkaline cleavage/dUTPcopy route','library'),)
@dataclass(frozen=True)
class RiboseCleavage:
    reference:str
    site:int
    def __post_init__(self):
        if not self.reference or set(self.reference)-set('ACGT') or type(self.site)is not int or not 0<=self.site<len(self.reference)-1:raise ValueError('explicit embeddedrNMP and downstreamfragment required')
    @property
    def upstream(self):return self.reference[:self.site+1]
    @property
    def downstream(self):return self.reference[self.site+1:]
    @property
    def upstream_end(self):return 'rNMP2′,3′cyclicphosphate class'
    @property
    def downstream_end(self):return '5′OH'
    def selected_complement(self):return revcomp(self.downstream).replace('T','U')

def select_source_input(g,key,current,label):
    selected=g.product(key+'_selected',label,K.MOLECULAR,state=current.state)
    rest=g.product(key+'_stock','Untransferred source material, if any; amount and later use unassigned',K.MOLECULAR,state=current.state)
    g.react(key,label,S.SPLIT,((current,R.TEMPLATE),),(selected,rest))
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source input allocation; no positive remaining amount or later use inferred.'))
    return selected

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Individual sourceearly-logS.pombe strains, rnh201 background/polymerasemutants not pooledbiologicaldefaults',K.MOLECULAR,state=pop('source genomicDNA with embeddedrNMP sites and othercellularmaterial'),primary=True)
    x=op(g,'extract','Source QIAGEN100/GGenomic-tip DNA preparation',S.MAGIC,cells,pop('recovered genomicDNA and unquantifiedcarryover; nativeends/rNMP sites unassigned'),'source historicalGenomic-tip kit, exact internals unreviewed')
    x=select_source_input(g,'hydrolysis_input',x,'Actual20ug genomicDNA input to alkaline treatment')
    naoh=g.input('naoh','Actual0.23M NaOH molecular reagent',K.MOLECULAR)
    x=op(g,'alkali','Source20µgDNA/0.23MNaOH55C2h; cleavage3′ofembeddedrNMP',S.CHEMICAL_CONVERSION,x,pop('upstream rNMP2′,3′cyclicphosphate classes plusdownstreamDNA5′OH fragments and all other cleavage/background/remnant classes'),'Source55C2h reaction environment; molecular NaOH supplied separately',additional=((naoh,R.STRUCTURAL),))
    x=select_source_input(g,'gel_input',x,'Actual10ug source ssDNA loaded for gel separation')
    x=cleanup(g,'gel','Source10µgssDNA2%TBEgel100V2h;acridineorange5µg/ml2h,destainovernight;300–500bp gelregion/NucleoSpin extraction',x,pop('recovered selectedssDNAfragmentpopulation and possiblecarryover; exactendcombinations unassigned'))
    x=select_source_input(g,'copy_input',x,'Actual100ng purified source ssDNA for BioPrime duplex conversion')
    supplies=g.input('copy_supplies','Actual BioPrime/dNTP inventory; packaging, number of stocks and kit internals unassembled',K.MOLECULAR)
    copy_nt=g.product('copy_nt','Actual dATP,dCTP,dGTP,dUTP contribution; dTTP substituted by dUTP in source copy reaction',K.MOLECULAR)
    copy_enzyme=g.product('copy_enzyme','Actual BioPrime DNA-copying catalyst class; identity unreviewed',K.CATALYST)
    copy_priming=g.product('copy_priming','Actual DNA-priming functional class; source whole sequence/priming positions unreviewed',K.OLIGO)
    copy_buffer=g.product('copy_buffer','Actual DNA-copy reaction-buffer class; recipe unreviewed',K.CONDITION)
    g.react('copy_components','Expose source DNA-copy supply classes; inventory grouping, not a documented premix',S.GENERIC,((supplies,R.STRUCTURAL),),(copy_nt,copy_enzyme,copy_priming,copy_buffer),changes_molecule=False)
    copied=g.product('copy_product','Original ssDNA with dU-bearing complementary DNA and incompletely copied/reagent classes; exact geometry unassembled',K.MOLECULAR,state=pop('original DNA/dU-copy, primers and reaction remnants'))
    g.react('copy','Actual source100ng BioPrime duplex conversion with dUTP substitution',S.MAGIC,((x,R.TEMPLATE),(copy_nt,R.STRUCTURAL),(copy_enzyme,R.CATALYST),(copy_priming,R.STRUCTURAL),(copy_buffer,R.CONDITION)),(copied,),reason='Positive source DNA-copy operation; actual enzyme/primer identities, copied fraction and ends unreviewed.')
    x=copied
    x=cleanup(g,'copy_cleanup','Source AMPureXP physical cleanup before separate analytical allocations',x,x.state)
    pico=g.product('picogreen','Source PicoGreen analytical allocation; amount/assay internals unassembled',K.QC)
    bio=g.product('bioanalyzer','Source Bioanalyzer analytical allocation; amount/assay internals unassembled',K.QC)
    all_dna=g.product('all_library_dna','ALL source20–60ng converted DNA after analytical allocations for library preparation',K.MOLECULAR,state=x.state)
    g.react('copy_qa_allocations','Separate source analytical material; ALL remaining20–60ngDNA enters NEBNext',S.SPLIT,((x,R.TEMPLATE),),(pico,bio,all_dna))
    g.mark_terminal(Terminal(pico.id,Disposition.QC,'Analytical material is not production library input.'))
    g.mark_terminal(Terminal(bio.id,Disposition.QC,'Analytical material is not production library input.'))
    x=all_dna
    x=op(g,'end_clean','Source20–60ng NEBNextUltra DNA endcleaning, historicalexactreaction unreviewed',S.MAGIC,x,x.state,'source historicalNEBNextUltra; no modernkitmodule default')
    ds=mixture(('sourceadapter',(role('originalNEBNextUltra adapterpolymers/end modifications unreviewed'),)))
    donor=g.input('adapter','sourceNEBNextUltra adapter',K.ADAPTER,state=ds);ligase=g.input('ligase','sourcekitligationreagents',K.CATALYST)
    ligated=g.product('ligated','source adapter-ligatedduplexDNA preparation beforeUSER',K.MOLECULAR,state=pop('source adapter-linked originaltemplate/dUcopy DNA and unligatedspecies/reagentremnants'))
    g.react('ligate','Source endclean/adaptorligation expresslyWITHOUT USERtreatment',S.LIGATION,((x,R.TEMPLATE),(donor,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),geometry=AdapterLigationGeometry('unavailable','unavailable',(ds,),'Original historicaladapter/end assembly unreviewed; no sourceabsenceclaim.'))
    x=cleanup(g,'insert_selection','Source AMPureXP insertselection250–600bp BEFORE USER',ligated,ligated.state)
    x=op(g,'user','Source purifiedDNA thenUSERenzyme; dU-copy cleavage doesnot physically discard fragments',S.MAGIC,x,pop('retained originalstrands, USERcleaveddU-copy/remnants, possible incompletelydigested/carryover molecules; exactsourceadaptor opening unassembled'),'Actual source USER enzyme; exact kitdUadaptorversuscopy susceptibility unresolved',K.CATALYST,R.CATALYST)
    f=g.input('primer_forward','Source multiplexPCRforward primerbody/selectedmemberunreviewed',K.OLIGO)
    r=g.input('primer_reverse','Source multiplexPCRindexprimer, study member/wholepolymerunreviewed',K.OLIGO)
    pcr_supply=g.input('pcr_supply','Actual source PCR supply inventory; packaging and recipe unreviewed',K.MOLECULAR)
    cat=g.product('pcr_enzyme','Actual source PCR catalyst class; identity unreviewed',K.CATALYST)
    nts=g.product('pcr_nt','Actual PCR nucleotide class; composition unreviewed, no dUTP-copy recipe inherited',K.MOLECULAR)
    buffer=g.product('pcr_buffer','Actual PCR buffer class; composition unreviewed',K.CONDITION)
    g.react('pcr_components','Expose actual PCR functional supplies without inventing a premix',S.GENERIC,((pcr_supply,R.STRUCTURAL),),(cat,nts,buffer),changes_molecule=False)
    amp=g.product('amp','source13cyclePCRpreparation with templates/remnants retained',K.LIBRARY,state=pop('source amplifiedsequencingDNA plusoriginaltemplates/USERfragments andPCRreagentremnants'))
    g.react('pcr','Source13cycles after delayedUSER',S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='Historicalsourceprimer/adaptor/end geometryunreviewed; no sourceabsenceclaim.'))
    x=cleanup(g,'library','Source finalAMPureXPpurification beforesequencing',amp,amp.state)
    native(g,x,'Source100bp paired-endHiSeq2500 or150bp paired-endNextSeq500;exactsourceindexmember/priminggeometryunreviewed.')
    g.add_source_item(SourceItem('historical_assembly','Original BioPrime/NEBNextUltra and Zhang2012 methoddocuments; finalVORapplicability','Associatedmain115 names actualdirectoriginals; preciseversions/assemblynotreviewedincap.'))
    return g

def selftest():
    from checks import Check
    c=Check();x=RiboseCleavage('ACGATCA',2);g=protocol_graph()
    c('HydrolyzedrNMP belongsupstream,notfirstdownstreambase',(x.upstream,x.downstream),('ACG','ATCA'))
    c('Alkalineendclassesaredistinct',(x.upstream_end,x.downstream_end),('rNMP2′,3′cyclicphosphate class','5′OH'))
    c('SelecteddUTPcopyiscomplementaryandphysicalU-bearing',x.selected_complement(),'UGAU')
    c('USERconsumesphysicallysizeselectedligatedDNA',g.steps['user'].inputs[0].material,'insert_selection')
    c('USERremnantsretaineduntilPCR/physicalcleanup',g.steps['user'].outputs,('user_product',))
    c('PCRusespostUSERpreparation',g.steps['pcr'].inputs[0].material,'user_product')
    c('Hydrolysis feeds actual10ug gel allocation',g.steps['gel_input'].inputs[0].material,'alkali_product')
    c('DNA-copy reaction consumes selected100ng',g.steps['copy'].inputs[0].material,'copy_input_selected')
    c('ALL post-QA20–60ng DNA feeds NEBNext without another loading remainder',g.steps['end_clean'].inputs[0].material,'all_library_dna')
    c('USER is supplied as an actual enzyme catalyst',(g.materials['user_stock'].kind,g.steps['user'].inputs[1].role),(K.CATALYST,R.CATALYST))
    c('Analytical fractions have QC fates',(g.terminals['picogreen'].disposition,g.terminals['bioanalyzer'].disposition),(Disposition.QC,Disposition.QC))
    c('PCR supplies do not inherit copy dUTP composition',g.producer_of('pcr_nt').id,'pcr_components')
    c('Libraryassemblyboundarypreventspromotion',bool(g.publication_problems()))
    c('Privatepreviewhasnofinishedmarker','data-protocol-complete=' not in render_page())
    c.report()
