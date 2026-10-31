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
    out = g.product(ident, 'recovered source preparation', K.MOLECULAR, state=state)
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
def op(g,key,label,kind,x,out_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION,extra_inputs=()):
    reagent=g.input(key+'_stock',reagent_label,reagent_kind)
    product=g.product(key+'_product',label,K.MOLECULAR,state=out_state)
    kw={'reason':'Source operation known; exact molecular substrate/product relation and mixture fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(key,label,kind,((x,R.TEMPLATE),(reagent,reagent_role))+tuple(extra_inputs),(product,),**kw)
    return product


def native(g,x,limitation):
    # Source acquisition may be stated without fabricated sequencing primers/cycle counts.
    g.sequence(x,SequencingPlan(limitation=limitation))

def pending(g,label):
    g.add_source_item(SourceItem('unreviewed_assembly',label,'Named source operations are represented; historical end/kit/primer assembly and complete original-version applicability unreviewed.'))



PAGE_FILE='fcab-seq.html'
TITLE='fCAB-seq'
NOTES='01_fcab-seq.html'
SOURCE='AssociatedmanuscriptNIHMS463944 mainResultsbase-resolutionsection/Methods146–148/Fig6; DOI10.1016/j.cell.2013.04.001. ChemistrySI01/02gated.'
CAVEAT='Main148 says bisulfite thermal programme run twice or high-throughput bisulfite amplicon sequencing (see Extended Procedures); exact applicability of twice-program to this selectedamplicon route is unreviewed, not silently assumed or declared a publication error. SelectedsourceampliconfCABroute; ChIP-fCAB and separatefC-Sealaffinityassay aredistinctsourcebranches. O-ethylhydroxylamine protection isnot sodiumborohydridereduction. No wholelibrary/sourceprimerbody or idealconversiondistribution inferred.'
READOUT='SourceFig6deepampliconNGS ofprotectedandordinarybisulfitepreparations; exactplatform/cycles/primergeometry unreviewed.5fCis computationallyinferred bycomparison ofprotectedC-likecalls withordinarybisulfite; protectiondoesnotchemicallyreplace5fCwithC.'
PREVIEWS=(('ProtectedfCABampliconpreparation','library_protected'),)
@dataclass(frozen=True)
class FormylDNA:
    reference:str
    sites:frozenset[int]
    protected:bool=False
    def __post_init__(self):
        if not self.reference or set(self.reference)-set('ACGT') or any(type(i)is not int or not 0<=i<len(self.reference) or self.reference[i]!='C' for i in self.sites):raise ValueError('selected5fC requiresactualC-originreferencepositions')
    def protect_selected(self):return FormylDNA(self.reference,self.sites,True)
    @property
    def modifications(self):return {i:'EtONH2-protected5fCoxime' if self.protected else '5fC' for i in self.sites}
    def source_reported_class(self):return {i:'C-like' if self.protected else 'T-like' for i in self.sites}

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    dna=g.input('dna','SourcegenomicDNA, selecteddeepampliconroute; native5fCsitesunassigned',K.MOLECULAR,state=pop('source5fC/othercytosine-containinggenomicDNA'),primary=True)
    x=op(g,'sonicate','Source sonicatedgenomicDNA average400bp; synthetic76mer/ChIPDNA areseparatesubstrates',S.FRAGMENTATION,dna,dna.state,'sourcesonication')
    prot=g.product('protected_input','SourcealiquotselectedforEtONH2',K.MOLECULAR,state=x.state);bs=g.product('bs_input','SourceordinaryBScomparisonaliquot',K.MOLECULAR,state=x.state)
    g.react('compare_split','Sourceprotectedversusordinarybisulfitecomparisonpreparations',S.SPLIT,((x,R.TEMPLATE),),(prot,bs))
    protection=g.input('protection_programme','Source100mMMESpH5 and37C2h condition; physical EtONH2 supplied separately',K.CONDITION)
    prot=op(g,'protect','Source10mMO-ethylhydroxylamine/100mMMESpH5,37C2h',S.CHEMICAL_CONVERSION,prot,pop('5fCoxime-protectedDNA pluspossibleunprotected5fC/othercytosines/reagentremnants'),'Actual source O-ethylhydroxylamine contribution10mM; no NaBH4reduction/5hmCsubstitution',K.MOLECULAR,R.STRUCTURAL,extra_inputs=((protection,R.CONDITION),))
    prot=cleanup(g,'protect_cleanup','SourceQiagennucleotideremovalkit',prot,prot.state)
    for branch,current in (('protected',prot),('ordinary',bs)):
        current=op(g,'bisulfite_'+branch,'SourceEpiTect/delegated highthroughputbisulfite route; general twice-program wording, selectedamplicon applicability unreviewed',S.MAGIC,current,pop('sourcebisulfiteprocessedDNA,'+('protected5fC-derivativeclasses/' if branch=='protected' else 'ordinaryunprotected5fC-derivedclasses/')+'othercytosineproducts andpossibledegradation/carryover; exactdistributionsunassembled'),'HistoricalEpiTect/amplicon supplies; main148 says twice-program or highthroughputamplicon sequencing, exact selectedbranch programme/fullphysicalbaseproducts unreviewed')
        kit=g.input('amplicon_'+branch,'SourceExtendedProcedureshighthroughputbisulfiteampliconandlibraryrecipe; originalgated',K.CONDITION)
        lib=g.product('library_'+branch,'Source'+branch+'deepampliconsequencinglibrary',K.LIBRARY,state=current.state)
        g.react('library_'+branch+'_prep','Source '+branch+' amplicon/librarycontinuation',S.MAGIC,((current,R.TEMPLATE),(kit,R.CONDITION)),(lib,),reason='OriginalFig6establishesNGS butpreciseprimer/PCR/adapter/recoverygeometrybehindunreadSI; no inheritedrecipe.')
        native(g,lib,'SourceFig6deepampliconNGS; exactplatform/run/primer/endgeometryunreviewed.')
    g.add_source_item(SourceItem('gated_si','ActualExtendedExperimentalProceduresSupplement01/02','https://pmc.ncbi.nlm.nih.gov/articles/instance/3657391/bin/NIHMS463944-supplement-01.pdf and supplement-02.pdf; originalgatesrecorded.'))
    return g

def selftest():
    from checks import Check
    c=Check();x=FormylDNA('ACCG',frozenset({1}));p=x.protect_selected();g=protocol_graph()
    c('Protection consumes actual chemical contributor rather than only a programme',g.materials['protect_stock'].kind,K.MOLECULAR)
    c('Original ExtendedProcedures obligation remains uncovered','gated_si' in g.source_coverage,False)
    c('ProtectionretainsC-originDNAsequence',p.reference,x.reference)
    c('PhysicaloximederivativediffersfromreportedC-likebase',p.modifications[1],'EtONH2-protected5fCoxime')
    c('Comparativecallsdonotreplacephysicaltemplate',p.source_reported_class()[1],'C-like')
    c('Protectedandordinarypreparationsremainseparate',set(g.sequencing),{'library_protected','library_ordinary'})
    c('SourceEtONH2cleanupoccursbeforebisulfite',g.steps['bisulfite_protected'].inputs[0].material,'protect_cleanup')
    c('Gatedlibrarypreventsfinishedpublication',bool(g.publication_problems()))
    c('Privatepreviewhasnofinishedmarker','data-protocol-complete=' not in render_page())
    c.report()
