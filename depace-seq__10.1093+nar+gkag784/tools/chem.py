"""Protocol-owned source-limited model; no finished-publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, Construct, feature, revcomp, oligo, panel, workflow_from_protocol_graph, workflow_panel
from protocol_model import role, molecule
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, Terminal, Disposition, SourceItem
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry
from page import head, info
import illumina as il
RHP5_SEQUENCE=il.TRUSEQ_READ2[2:]

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

def selftest():
    from checks import Check
    c = Check(); g = protocol_graph()
    c('Pending source/assembly boundary prevents finished publication', bool(g.publication_problems()))
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())
    c('Every physical DNA output has a declared downstream operation or terminal', not any(m.kind in {K.MOLECULAR,K.LIBRARY} and m.id not in g.terminals and not any(port.material==m.id for step in g.steps.values() for port in step.inputs) for m in g.materials.values()))
    c('R release consumes residual pellet after L rather than an independent aliquot',g.steps['r_pellet_wash_step'].inputs[0].material=='residual_pellet')
    c('RNaseH2/Q5 remnants remain molecular until column cleanup', all(g.materials[mid].kind is not K.WASTE for mid in g.steps['r_activation_extension'].outputs))
    c('L/R amplified libraries remain separate physical outputs', 'l_library' in g.terminals and 'r_library' in g.terminals)
    rh=g.materials['rhp5'];base=''.join(seg.top for seg in rh.state.strands['rhp5'].segs)
    c('rhP5 physicalnucleotides preserve actual riboguanosine',base,rhp5_chemistry()['sequence'])
    from checks import Source,have
    from storage import sources_root
    import re
    source=Source(sources_root()/HERE.parent.name/'nar_gkag784_PMC13443447.xml.txt','Published main106 wholemodifiedrhP5')
    if have(c,source):
        match=re.search(r'5′NH2C6-([ACGT]+)/rG/([ACGT]+)-3′C3 Spacer',source.path.read_text())
        c('Original modifiedrhP5 ordering located',match is not None,True)
        if match:
            c('Whole modeled rhP5 matches source including riboG',base,match.group(1)+'G'+match.group(2))
            c('RNA sugar coordinate follows original rG position',rh.chemical_state['sugar_blocks'][1]['start'],len(match.group(1)))
    c.report()

PAGE_FILE='depace-seq.html'
TITLE='DEPACE-seq'
NOTES='01_depace-seq.html'
SOURCE='<a href="https://doi.org/10.1093/nar/gkag784">Defining article</a>; Published main Methods94–120; SI Adapter preparationp7, genomic control validationpp13–14, rendered TableS1pp17–18'
CAVEAT='Selected published main genomic route. PAB conjugates AP sites; cleavage does not remove fragments from the tube. SI control experiments use NEB/NextSeq conditions and do not replace main VAHTS/HiSeq-or-DNBSEQ recipe. Left/right release identity and exact adduct-bearing end geometry require further assembly; no ideal molecular partition assumed.'
READOUT='Main HiSeq2000 or DNBSEQ-T7 PE150. Left/right libraries are physically separate; their mapped coordinate intersection is computational. One-library putative AP calls are a separate analysis branch.'
PREVIEWS=(("L-library route","l_library"),("Sequential residual-pellet R-library route","r_library"))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    current=g.input('dna','source genomic DNA for selected main route; genomic bases unassigned',K.MOLECULAR,state=molecule('source DNA',(role('genomic DNA; coordinates and sequence unassigned'),)),primary=True)
    reagent=g.input('shear_reagent','Sonicate genomicDNA200–500bp',K.CONDITION)
    product=g.product('shear_product','Sonicate genomicDNA200–500bp',K.MOLECULAR,state=molecule('source product classes',(role('fragmented genomicDNA with AP sites at unassigned positions'),)))
    g.react('shear','Sonicate genomicDNA200–500bp',S.FRAGMENTATION,((current,R.TEMPLATE),(reagent,R.CONDITION)),(product,))
    current=product
    reagent=g.input('label_reagent','Source2mMPAB molecular AP-site conjugation reagent; occupancyunassigned',K.MOLECULAR)
    label_conditions=g.input('label_conditions','Source4xPBS overnight20-25C; no reactionyieldassumption',K.CONDITION)
    product=g.product('label_product','PAB2mM4xPBS overnight room temperature',K.MOLECULAR,state=molecule('source product classes',(role('PAB-AP covalent conjugates on DNA plus unlabelled/background DNA'),)))
    g.react('label','PAB2mM4xPBS overnight room temperature',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(reagent,R.STRUCTURAL),(label_conditions,R.CONDITION)),(product,))
    current=product
    current=cleanup(g,'label_clean_product','Source Zymo-Spin cleanup after PAB',current,molecule('recovered DNA',(role('recovered PAB-conjugated and background DNA'),)))
    kit=g.input('end_repair_kit','Actual VAHTSND607v3 molecular end-preparation supplies; main genomic branch, internalrecipeunreviewed',K.MOLECULAR)
    prepared=g.product('prepared','source end-repaired/dA-tailed DNA',K.MOLECULAR,state=downstream_state('PAB-associated DNA with source dA-tailed ends; exact kit transition unassembled'))
    g.react('end_repair','Main source end repair and dA-tail',S.MAGIC,((current,R.TEMPLATE),(kit,R.STRUCTURAL)),(prepared,),reason='Source names end repair/dA-tail kit; exact proprietary operation and termini not reconstructed.')
    current=source_ligation(g,'p5','CustomP5 ligation20C3h',prepared,p5_scene())
    current=cleanup(g,'p5_cleanup','Main1xVAHTS DNA Clean Beads after P5 ligation',current,current.state)
    reagent=g.input('cleave_reagent','EndoIV37C30min at source PAB-labelled sites',K.CATALYST)
    product=g.product('cleave_product','EndoIV37C30min at source PAB-labelled sites',K.MOLECULAR,state=molecule('source product classes',(role('source free3′OH cleavage class and corresponding PAB-associated end fragments, complementary DNA and possible uncleaved DNA; all retained'),)))
    g.react('cleave','EndoIV37C30min at source PAB-labelled sites',S.DIGESTION,((current,R.TEMPLATE),(reagent,R.CATALYST)),(product,))
    current=product
    current=cleanup(g,'cleave_clean_product','Source Zymo-Spin cleanup after EndoIV',current,molecule('recovered DNA',(role('recovered cleaved/uncleaved DNA population; fragment recovery unquantified'),)))
    reagent=g.input('capture_reagent','Prewashed streptavidin beads, binding buffer and six source washes',K.SELECTION)
    blocker=g.input('capture_blocker','Source2ugpoly dI:dC DNAblocker contribution; boundcarryoverunquantified',K.MOLECULAR)
    product=g.product('captured','Prewashed streptavidin beads, binding buffer and six source washes',K.MOLECULAR,state=molecule('source product classes',(role('captured biotin-associated DNA complexes and possible blocker/background carryover'),)))
    other=g.product('capture_other','physically separated source flowthrough/washes; partition unquantified',K.WASTE)
    g.react('capture','Prewashed streptavidin beads, binding buffer and six source washes',S.AFFINITY_CAPTURE,((current,R.TEMPLATE),(reagent,R.SELECTION),(blocker,R.STRUCTURAL)),(product,other))
    g.discard(other,'Actual source isolation/washes; no perfect specificity asserted.')
    current=product
    alkali=g.input('alkali','100mMNaOH roomtemperature5min; source magnetic separation',K.CONDITION)
    left=g.product('left','L-source supernatant after alkaline release',K.MOLECULAR,state=downstream_state('source L-eluate DNA; exact cleavage-side strand/end identity unassembled'))
    pellet=g.product('residual_pellet','residual bead-associated DNA after L release',K.MOLECULAR,state=downstream_state('source residual biotin-associated DNA; possible residual complementary/background DNA'))
    g.react('l_release','Sequential L-release: alkaline denaturation then collect supernatant; retain residual pellet for R',S.STRAND_SEPARATION,((current,R.TEMPLATE),(alkali,R.CONDITION)),(left,pellet))
    left=cleanup(g,'l_neutralize','L-source HCl neutralization and Zymo column workup',left,left.state)
    heat=g.input('l_heat_condition','95C5min then immediate ice',K.CONDITION)
    denatured=g.product('l_denatured','denatured L-eluate DNA',K.MOLECULAR,state=left.state)
    g.react('l_denature','Denature L DNA before customP7 ligation',S.THERMAL_TREATMENT,((left,R.TEMPLATE),(heat,R.CONDITION)),(denatured,))
    left=source_ligation(g,'l_p7','L-P7 Instant Sticky-End ligase25Covernight',denatured,p7_scene())
    left=cleanup(g,'l_ligation_cleanup','Main1xVAHTS DNA Clean Beads after L ligation',left,left.state)
    left=source_pcr(g,'l_pcr','Main VAHTSHiFi/i5/i7PCR15cycles',left)
    left=cleanup(g,'l_library','Main1xVAHTS size selection after L PCR',left,left.state)
    washed=cleanup(g,'r_pellet_wash','Residual pellet1xBufferI wash after L release',pellet,pellet.state)
    formamide=g.input('formamide','95%formamide10mMEDTA90C3min; magnetic separation',K.CONDITION)
    right=g.product('r_eluted','R-source formamide-eluted DNA',K.MOLECULAR,state=downstream_state('source R-eluate DNA; tagged-end/adduct partition unassembled'))
    residual=g.product('final_pellet','residual bead-associated material after R elution',K.MOLECULAR,state=downstream_state('possible residual DNA/blocker/adduct material; final disposal not source-resolved'))
    g.react('r_release','Formamide release from residual pellet after L',S.STRAND_SEPARATION,((washed,R.TEMPLATE),(formamide,R.CONDITION)),(right,residual))
    g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Do not assume complete elution or disposal of residual material.'))
    right=cleanup(g,'r_elution_cleanup','R-source Zymo column purification',right,right.state)
    primer=g.input('rhp5','source5amine/rG/3C3 blocked rhP5',K.OLIGO,state=rhp5_scene(),chemical_state=rhp5_chemistry())
    enzyme=g.input('r_extension_enzymes','ThermostableRNaseH2 and Q5; main temperature program',K.CATALYST)
    extension_nt=g.input('r_extension_dntps','Source Rextension1ul dNTPs,300uMfinal; per-base/total conventionunassigned',K.MOLECULAR)
    extension=g.product('r_extended','rhP5/RNaseH2/Q5 products and cleaved primer remnants',K.MOLECULAR,state=mixture(('extension products',(role('source R-primer extension DNA; activated primer-end identity unassembled'),)),('reaction remnants',(role('RNaseH2-cleaved blocked-primer remnants and residual templates; retained until column workup'),))))
    g.react('r_activation_extension','rhP5 RNaseH2 activation and Q5 primer extension; cleavage-product geometry unreviewed',S.MAGIC,((right,R.TEMPLATE),(primer,R.STRUCTURAL),(enzyme,R.CATALYST),(extension_nt,R.STRUCTURAL)),(extension,),reason='Source reaction/modified primer established; exact RNaseH2 cut and extension end geometry not yet assembled.')
    right=cleanup(g,'r_extension_cleanup','R-source Zymo column after primer extension',extension,downstream_state('recovered R-extension DNA; remnant carryover unquantified'))
    right=source_ligation(g,'r_p7','R-P7 Instant Sticky-End ligase25Covernight',right,p7_scene())
    right=cleanup(g,'r_ligation_cleanup','Main1xVAHTS DNA Clean Beads after R ligation',right,right.state)
    right=source_pcr(g,'r_pcr','Main R-PCR15cycles followsLcyclingprogram; exactRformulationunassigned',right)
    right=cleanup(g,'r_library','Main1xVAHTS size selection after R PCR',right,right.state)
    g.mark_terminal(Terminal(left.id,Disposition.UNRESOLVED,'L library source operations complete, molecular/run geometry unassembled.'))
    current=right
    g.mark_terminal(Terminal(current.id,Disposition.UNRESOLVED,'Exact modified P5/P7/rhP5 assembly, exact molecular partition through sequential L/R release and full read geometry remain model work; vendor ND607 manual directly named but uninspected.'))
    g.declare_reviewed_scope(('dna','shear_reagent','shear','label_reagent','label_conditions','label','end_repair_kit','end_repair','capture_reagent','capture_blocker','capture','rhp5','r_extension_enzymes','r_extension_dntps','r_activation_extension'),SOURCE)
    g.add_source_item(SourceItem('pending_assembly','Unreviewed source/assembly boundary','Exact modified P5/P7/rhP5 assembly, exact molecular partition through sequential L/R release and full read geometry remain model work; vendor ND607 manual directly named but uninspected.'))
    return g

def downstream_state(label):
    return mixture(('source-derived DNA',(role(label),)),('possible carryover',(role('background DNA, residual templates and reaction components; exact recovery unquantified'),)))

def source_ligation(g,ident,label,template,adapter):
    donor=g.input(ident+'_adapter','Source adapter strands, fullassemblyunreviewed',K.ADAPTER,state=adapter,chemical_state=p7_chemistry() if ident.endswith('p7') else {'scope':'SourceP5 modifiedstrandinventory; wholejunctionunassembled'})
    ligase=g.input(ident+'_ligase',label,K.CATALYST)
    out=g.product(ident+'_product',label,K.MOLECULAR,state=downstream_state('source adapter-ligated DNA roles; exact junction assembly unreviewed'))
    geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter,), 'Source-known adapter identity; substrate end and sealed junction assembly unreviewed, no source-absence claim.')
    g.react(ident,label,S.LIGATION,((template,R.TEMPLATE),(donor,R.STRUCTURAL),(ligase,R.CATALYST)),(out,),geometry=geometry)
    return out

def source_pcr(g,ident,label,template):
    forward=g.input(ident+'_forward','Source first PCR primerclass; precise indexed ordering/studyassignmentunreviewed',K.OLIGO,state=molecule('first PCR primer',(role('source primer; exact study-specific identity unassigned'),)))
    reverse=g.input(ident+'_reverse','Source second PCR primerclass; precise indexed ordering/studyassignmentunreviewed',K.OLIGO,state=molecule('second PCR primer',(role('source primer; exact study-specific identity unassigned'),)))
    supplies=g.input(ident+'_supplies','Source25ulVAHTSHiFiAmplificationMix forL' if ident=='l_pcr' else 'Actual source R-PCR supplies; sourceonlystatesLcyclingprogram, formulation/doseunassigned',K.MOLECULAR)
    polymerase=g.product(ident+'_enzyme','Polymerase functional component ofsameactualPCRsupplies; no extra freshaddition',K.CATALYST)
    nts=g.product(ident+'_nts','Nucleotide component ofsameactualPCRsupplies; dose/packagingunassigned',K.MOLECULAR)
    g.react(ident+'_inventory','Components ofactualsourcePCRsupplies; RformulationnotassumedidenticaltoL',S.GENERIC,((supplies,R.TEMPLATE),),(polymerase,nts),changes_molecule=False)
    out=g.product(ident+'_product',label,K.LIBRARY,state=downstream_state('source amplicon library roles; primer landing and exact end geometry unassembled'))
    g.react(ident,label,S.PCR,((template,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(polymerase,R.CATALYST),(nts,R.STRUCTURAL)),(out,),geometry=PCRGeometry.unavailable(template.state,reason='Source PCR known; exact selected primer identities and full library geometry unreviewed.'))
    return out

def p7_scene():
    sc=Scene()
    sc.strand('top',[Segment('P7 top','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCG')],label='P7 adapter top',mod5='phosphate',mod3='C3 spacer')
    sc.strand('bottom',[Segment('P7 bottom fixed','CTACACGACGCTCTTCCGATCT'),Segment('six-base degenerate ligation window','NNNNNN',placeholder=True)],label='P7 adapter bottom; random window not assigned UMI',mod5='C3 spacer',mod3='C6 amino; source phosphorothioate between final two Ns')
    return sc

def rhp5_scene():
    sc=Scene()
    sc.strand('rhp5',[Segment('rhP5DNAbody',RHP5_SEQUENCE[:-5]),Segment('source riboguanosine','G',note='Source/rG/ sugar isRNA; structuredsugarblocks retained oninputchemistry'),Segment('rhP5DNAtail',RHP5_SEQUENCE[-4:])],label='source rhP5',mod5='C6 amino',mod3='C3 spacer')
    return sc

def p5_scene():
    sc=Scene()
    sc.strand('top',[Segment('P5 top',RHP5_SEQUENCE)],label='P5 adapter top',mod5='C6 amino')
    sc.strand('bottom',[Segment('P5 bottom','GATCGGAAGAGCACACGTCTGAACTCCAGTC')],label='P5 adapter bottom',mod5='phosphate',mod3='C3 spacer')
    return sc

def oligos():
    return [panel(p5_scene().rows(),caption='Source P5 inventory; duplex/junction assembly unreviewed'),panel(p7_scene().rows(),caption='Source P7 inventory; degenerate window not assigned UMI'),panel(rhp5_scene().rows(),caption='Source blocked rhP5 inventory')]


def rhp5_chemistry():
    body=RHP5_SEQUENCE[:-5]
    return {'polymer':'mixed DNA/RNA','sequence':body+'G'+'ATCT','five_prime':'C6 amino','three_prime':'C3 spacer','sugar_blocks':({'start':0,'end':len(body),'sugar':'DNA'},{'start':len(body),'end':len(body)+1,'sugar':'RNA'},{'start':len(body)+1,'end':len(body)+5,'sugar':'DNA'}),'source_ordering':'5NH2C6-'+body+'/rG/ATCT-3C3Spacer'}

def p7_chemistry():
    bottom='CTACACGACGCTCTTCCGATCT'+'N'*6
    return {'source':'PublishedSI TableS1p17','top':{'five_prime':'P','three_prime':'C3 spacer'},'bottom':{'sequence':bottom,'five_prime':'C3 spacer','three_prime':'C6 amino','phosphorothioate_after':(len(bottom)-2,)},'scope':'Knownsourcephysicalstrands/modifications only; no wholeligationjunction'}

if __name__=='__main__':selftest()
