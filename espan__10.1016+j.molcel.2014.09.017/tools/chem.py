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
def op(g,key,label,kind,x,out_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
    reagent=g.input(key+'_stock',reagent_label,reagent_kind)
    product=g.product(key+'_product',label,K.MOLECULAR,state=out_state)
    kw={'reason':'Source operation known; exact molecular substrate/product relation and mixture fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(key,label,kind,((x,R.TEMPLATE),(reagent,reagent_role)),(product,),**kw)
    return product


def native(g,x,limitation):
    # Source acquisition may be stated without fabricated sequencing primers/cycle counts.
    g.sequence(x,SequencingPlan(limitation=limitation))

def pending(g,label):
    g.add_source_item(SourceItem('unreviewed_assembly',label,'Named source operations are represented; historical end/kit/primer assembly and complete original-version applicability unreviewed.'))


PAGE_FILE='espan.html'
TITLE='eSPAN'
NOTES='01_espan.html'
SOURCE='Associated manuscript NIHMS667987 mainMethods152–156/Fig1A; DOI10.1016/j.molcel.2014.09.017. DetailedSI gated; actual erratum identity/disposition unreviewed.'
CAVEAT='Source-known serial proteinChIP and BrdU-IP retain strand-linked nascentDNA before the citedMeyer2012 ssDNA library. No generic ChIP/ancientDNA recipe, oligo or end chemistry is substituted. Original correction and detailed published preparation remain unresolved; no complete final-version equivalence assumed.'
READOUT='Source HiSeq2000/2500 paired-end acquisition; exact cycles/primers/library polymers unassembled. Physical BrdU is retained through source-known selection; Watson/Crick counts and origin-relative leading/lagging bias are computational, not direct chemical strand labels or perfect target occupancy.'
PREVIEWS=(('Original yeast serial selection','library'),)
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Original individually prepared yeastW303 strains, not pooledgenotypes',K.MOLECULAR,state=pop('yeastchromatin/cellular material; exact DNA/protein populations unassigned'),primary=True)
    x=op(g,'label','Source Fig1A G1-release intoS with BrdU, with or withoutHU as distinctexperimentalconditions',S.MAGIC,cells,pop('physical BrdU-containingnascentDNA, parental/unlabelledDNA and cellularcomponents'),'source cell-cycle/analogue labeling; exact parameters gatedSI, no full labeling guarantee')
    x=op(g,'chip','Source target-protein chromatin immunoprecipitation, original internal fixation/fragmentation/wash/recovery unreviewed',S.MAGIC,x,pop('source protein-associated chromatin preparation with possible otherDNA/protein/reagentcarryover'),'source standardChIP and targetantibody; no invented purification recipe')
    x=op(g,'proteolyze','Source proteins digested AFTER proteinChIP; resulting peptides remain physically present until actual sourceworkup',S.DEPROTEINIZATION,x,pop('source protein-associatedDNA with physical proteolyticremnants; no nativeDNA removal or assumed crosslink reversal'),'source protein-digestionreagent/conditions unreviewed',K.CATALYST,R.CATALYST)
    x=op(g,'denature','Source protein-associatedDNA denatured before BrdU-IP',S.MAGIC,x,pop('source single-strandDNA classes includingBrdU-bearingnascentDNA and unlabelledDNA, plus remnants'),'source DNAdenaturation conditions; exact heat/crosslink chemistry gatedSI')
    ab=g.input('brdu_antibody','Source antiBrdU enrichment ofnascentssDNA; antibody/beadimplementation unreviewed',K.SELECTION)
    selected=g.product('nascent','source enrichedprotein-associatednascentssDNA preparation',K.MOLECULAR,state=pop('selected BrdU-bearingnascentDNA and possible nonspecificDNA/protein/reagentcarryover; exact recovery/strand ends unassembled'))
    other=g.product('nonselected','companion nonselected BrdU-IP fraction',K.MOLECULAR,state=pop('unboundDNA/proteolyticremnants and possible targetloss; exact source partition unassembled'))
    g.react('brdu_ip','Source BrdU immunoprecipitation afterDNAdenaturation',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(ab,R.SELECTION)),(selected,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Original workup/disposal unreviewed; no perfect depletion assumed.'))
    kit=g.input('meyer_method','Source cited Meyer2012 ssDNA library method, precise components/version unreviewed',K.CONDITION)
    lib=g.product('library','source strand-specific ssDNA-derived sequencinglibrary',K.LIBRARY,state=pop('source ssDNA-derivedlibrary and possible remnants; complete primer/adapter/strand geometry unassembled'))
    g.react('library_prep','Source cited ssDNA library continuation',S.MAGIC,((selected,R.TEMPLATE),(kit,R.CONDITION)),(lib,),reason='Original main names Meyer2012; SI/predecessor not read, no family default or inferred internal operations.')
    native(g,lib,'Source paired-end IlluminaHiSeq2000/2500; cycles and whole-library/priming geometry unreviewed. No nativeONT or fixed readlength substitution.')
    g.add_source_item(SourceItem('gated_si','Actual SupplementalExperimentalProcedures,sevenfigures/twotables','https://pmc.ncbi.nlm.nih.gov/articles/instance/4362665/bin/NIHMS667987-supplement-supplement_1.pdf; preparing-download gate.'))
    g.add_source_item(SourceItem('erratum','PubMed-linked MolCell2014Nov20;56(4):609 erratum, precise DOI/contentunresolved','OriginalPubMed25449133 identifieserratum; candidate2014.11.007 not proven correctionidentity.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Protein digestion consumes selected proteinChIP preparation',g.steps['proteolyze'].inputs[0].material,'chip_product')
    c('BrdU selection follows DNAdenaturation',g.steps['brdu_ip'].inputs[0].material,'denature_product')
    c('Library continuation consumes selectednascentDNA',g.steps['library_prep'].inputs[0].material,'nascent')
    c('Companion capturematerial remains molecular',g.materials['nonselected'].kind,K.MOLECULAR)
    c('Actualgatedoriginal/correction prevents finishedpublication',bool(g.publication_problems()))
    c('Private preview has no finishedmarker','data-protocol-complete=' not in render_page())
    c.report()
