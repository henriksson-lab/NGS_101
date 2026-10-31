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


PAGE_FILE='sns-seq.html'
TITLE='SNS-seq'
NOTES='01_sns-seq.html'
SOURCE='Published original PMC3205567 Resultsinitialoverview/Methods Nascentstrandabundanceassay source1368–1380; DOI10.1101/gr.124644.111.'
CAVEAT='Original2011 lambda-exonuclease SNS-NGS route with400–600bp sucrosefractions; no later2012saturation/2020optimized recipe merged. Wang2004 isolation and exact historicalGAII library assembly remain unreviewed. RNA-primer protection is source enrichment principle, not perfect resistance/specificity. Source also sequenced a separately sheared genomicDNA control; that independent library is not expanded in this selected nascent-strand graph or physically pooled here.'
READOUT='Source IlluminaGenomeAnalyzerII/Solexa acquisition of double-stranded random-primednascentDNA; source gives no assembled adapter/run geometry here. Read-density ratio to separately preparedshearedgenomicDNA and origin calling are computational; no sequencingbarcode/UMI assigned to nativeRNAleader.'
PREVIEWS=(('Original2011 nascent-DNA route','library'),)
from dataclasses import dataclass
@dataclass(frozen=True)
class NascentChain:
    rna:str
    dna:str
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU') or not self.dna or set(self.dna)-set('ACGT'):raise ValueError('explicit RNAprimer and DNApayload required')
    def rows(self):
        return molecule('selected originalRNA-DNA chain',(Segment('nativeRNAprimer',self.rna),Segment('nativeDNApayload',self.dna))).rows()
    def complementary_dna(self):
        # Conditional complete-copy class only; reaction efficiency/priming sites unassigned.
        return revcomp(self.dna)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    cells=g.input('cells','source K562 orMCF7 individualcellpreparation, not pooledcelltypes',K.BIOLOGICAL,state=pop('genomicDNA and RNA-primednascentDNA pluscellularmaterial'),primary=True)
    x=op(g,'isolation','Source genomic/nascentDNA isolated asWang2004 withmodifications; exact innersteps unreviewed',S.MAGIC,cells,pop('source RNA-primednascentDNA/brokengenomicDNA preparation; nativeends/othermaterial unassigned'),'source linkedWang2004 isolation')
    x=cleanup(g,'gradient','Source neutral5–30%sucrose/TNE,20h21000rpmSW40at25C;400–600bp fractions pooled,dialyzedagainstTE,precipitated',x,x.state)
    a=g.product('carrier_aliquot','source companion aliquot BEFORE plasmid addition',K.MOLECULAR,state=x.state)
    x2=g.product('selected_aliquot','source separated aliquot without addedplasmidcontrol',K.MOLECULAR,state=x.state)
    carrier=g.input('plasmid','source plasmidDNAcarrier/positivecontrol, exactidentity unreviewed',K.MOLECULAR,state=pop('sourceplasmidcontrol'))
    g.react('aliquot','Source sample divided into twoaliquots before adding plasmid toone',S.SPLIT,((x,R.TEMPLATE),),(x2,a))
    control=g.product('carrier_branch','source aliquot with addedplasmidcarrier/positivecontrol',K.MOLECULAR,state=pop('source nascentDNA plusaddedplasmidcontrol; exact latercontrol useunreviewed'))
    g.react('carrier_addition','Add source plasmid only to the companioncontrolaliquot',S.MAGIC,((a,R.TEMPLATE),(carrier,R.STRUCTURAL)),(control,),reason='Source addition without covalent joining; exact controlstudy member unreviewed.')
    g.mark_terminal(Terminal(control.id,Disposition.UNRESOLVED,'Source control aliquot retained; exact later control fate outside selected branch.'))
    x=cleanup(g,'organic','Source phenol-chloroform extraction thenethanolprecipitation/resuspend22µlDEPCwater',x2,x2.state)
    x=op(g,'lambda_digest','Source overnight37C lambdaexonuclease digestion in32µlbuffer plus2µlenzyme; fragmentsremainuntilworkup',S.DIGESTION,x,pop('lambda-treatedRNA-primedDNA, possible survivingbrokenDNA and digestionproducts; no perfect enrichment'),'source lambdaexonuclease',K.CATALYST,R.CATALYST)
    x=op(g,'heat','Source95C10min afterlambda digestion',S.MAGIC,x,x.state,'source heat condition; no perfectRNA/primer preservation or ideal inactivation assumed')
    x=cleanup(g,'post_digest','Source ethanolprecipitation/resuspend50µlsterilewater afterheat',x,x.state)
    qc=g.product('qc','source aliquot for knownoriginprimer validation',K.QC,state=x.state);rest=g.product('copy_input','remaining source nascentDNA forrandompriming',K.MOLECULAR,state=x.state)
    g.react('qc_split','Source validationPCR aliquot separated frompreparativeDNA',S.SPLIT,((x,R.TEMPLATE),),(qc,rest));g.mark_terminal(Terminal(qc.id,Disposition.QC,'KnownoriginPCR validation only; not nativeNGSlibrary primeraddition.'))
    supplies=g.input('copy_stock','Aggregate source Klenow/DNAPrime random-priming supplies; physical packaging, number of stocks and doses unassigned',K.MOLECULAR)
    enzyme=g.product('copy_enzyme','Actual named DNApolymeraseI Klenow contribution',K.CATALYST)
    primers=g.product('copy_primers','Source random-priming oligo contribution; whole bodies and distribution unassembled',K.OLIGO)
    nts=g.product('copy_nucleotides','DNA synthesis nucleotide supplies, inferred-required functional class; formulation/doses unassigned',K.MOLECULAR)
    g.react('copy_inventory','Expose functional inventory of aggregate source supplies; not a physical separation or asserted premix',S.GENERIC,((supplies,R.STRUCTURAL),),(enzyme,primers,nts),changes_molecule=False)
    x=g.product('copy_product','Source duplexcopiednascentDNA and reaction remnants',K.MOLECULAR,state=pop('source duplexcopiednascentDNA, originalRNA-DNA chains and primer/nucleotide/reagentremnants; exact productpositions unassembled'))
    g.react('copy','Source randompriming withDNApolymeraseIKlenow/DNAPrimeLabelingSystem1h37C',S.MAGIC,((rest,R.TEMPLATE),(enzyme,R.CATALYST),(primers,R.STRUCTURAL),(nts,R.STRUCTURAL)),(x,),reason='Actual named source catalyst and randompriming contributions exposed without stock-packaging inference; nucleotide class is inferred-required, source formulation and whole copied strands unassembled. No unsourced RNA hydrolysis,dUTP orUSER.')
    x=cleanup(g,'copy_cleanup','Source ethanolprecipitation/resuspend30µlsterilewater',x,x.state)
    kit=g.input('ga_library','Source historicalGAII library continuation; exact reagents/oligos/workup unreviewed',K.CONDITION)
    lib=g.product('library','source double-strandednascentDNA sequencinghandoff,1µg reported',K.LIBRARY,state=x.state)
    g.react('library_prep','Source GAII library handoff afterduplexconversion',S.MAGIC,((x,R.TEMPLATE),(kit,R.CONDITION)),(lib,),reason='Acquisition stated but historical libraryassembly remains unreviewed; not a modernkitrecipe.')
    native(g,lib,'Original IlluminaGAII/Solexa acquisition; cycles/priminglayout unreviewed, no paired-end orstrand-specific run guessed.')
    g.add_source_item(SourceItem('supplement','Actual citedSupplementalMethods/Table1/DataAnalysis, publisherlinkedoriginalidentityunresolved','MainMethods/footnote1411links publisherarticleDOI; actualfile not recovered bystandardfetcher inthispass.'))
    g.add_source_item(SourceItem('linked_isolation','Source-citedWang2004 isolation and originalDNAPrime/GAII documents','Actualversions and complete physical/strand assembly not reviewed.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();chain=NascentChain('ACGU','GATTACA')
    c('Analytical validation fraction is typed QC',g.materials['qc'].kind,K.QC)
    c('Copy catalyst derives from actual aggregate supplies',g.steps['copy_inventory'].inputs[0].material,'copy_stock')
    c('Copy consumes the produced catalyst rather than a fresh duplicated stock',g.steps['copy'].inputs[1].material,'copy_enzyme')
    c('Original supplement and precursor obligations remain uncovered',any(k in g.source_coverage for k in ('supplement','linked_isolation')),False)
    c('No plasmid appears in the split before carrier addition',g.materials['carrier_aliquot'].state is g.materials['selected_aliquot'].state,True)
    c('NativeRNAprimer remains explicit ratherthan a sequencingidentifier',chain.rna,'ACGU')
    c('Selected fullDNAcopy derives from actualDNApayload',chain.complementary_dna(),'TGTAATC')
    c('Nucleasedebris not discarded atdigestion',g.steps['lambda_digest'].outputs,('lambda_digest_product',))
    c('Physicalpostdigestworkup follows sourceheat',g.steps['post_digest_step'].inputs[0].material,'heat_product')
    c('Preparativecopy consumes remainingDNA, notvalidationPCR product',g.steps['copy'].inputs[0].material,'copy_input')
    c('Libraryremains sourcegeometrylimited',bool(g.publication_problems()))
    c('Privatepreview has nofinishedmarker','data-protocol-complete=' not in render_page())
    c.report()
