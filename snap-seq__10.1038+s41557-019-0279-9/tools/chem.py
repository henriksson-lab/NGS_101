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
    c.report()

import illumina as il
PAGE_FILE='snap-seq.html'
TITLE='snAP-seq — hot AP-site cleavage followed by negative biotin selection'
NOTES='01_snap-seq.html'
SOURCE='<a href="https://doi.org/10.1038/s41557-019-0279-9">Liu et al. (2019)</a>; associated accepted manuscriptPMC6589398 Methods105–115 and publisherSI Table1/adapter and extraction Methods.'
CAVEAT='Native HeLa AP route; glycosylase pretreatments, syntheticN10 controls and TruSeq input route remain separate. HIPS/CuAAC reacts with source aldehyde classes before selection; no perfect reaction yield/specificity asserted. Room-temperature alkali washes away complementary DNA, while HOT alkali releases AP-cleavage products. The second capture keeps the UNBOUND supernatant. Me notation and source P5top leadingG are preserved; full chemical structures, historical kit/end chemistry, final strand/run geometry remain unassembled.'
READOUT='Source MiSeq/NextSeq; exact run lengths and final strand/index-cycle layout unassembled. Source reads accumulate one nucleotide downstream of captured AP sites, not at an unchanged tagged base. Source controlN10 deduplication is not assigned to genomic libraries.'
PREVIEWS=(('Native HeLa AP enrichment and library route','library'),('Prepare first capture beads','first_beads'),('Prepare negative-selection beads','second_ready'))
P7_TOP=il.TRUSEQ_READ2
P7_CORE=revcomp(P7_TOP)[1:]
INDEX=feature('snap_sample_index','sample_index','unknown')
P7_BOTTOM=P7_CORE+'N'*6+revcomp(il.P7)
P5_TOP='G'+il.NEBNEXT_UNIVERSAL_PRIMER
P5_BOTTOM='GATCGGAAGAGCG'
PCR_FORWARD=il.NEBNEXT_UNIVERSAL_PRIMER[:44]
PCR_REVERSE=il.P7

def p7_adapter():
    core=revcomp(P7_CORE)
    if not P7_TOP.startswith(core):raise ValueError('source P7 oligos lack their supplied terminal overlap')
    sc=Scene();sc.strand('top',[Segment('P7 paired body',core),Segment('P7 source terminal base',P7_TOP[len(core):])],label='source P7top',mod5='Me(source notation; chemistry unassigned)')
    sc.anneal('bottom',[Segment('P7 body complement',P7_CORE),Segment('source TruSeqLT index','N'*6,placeholder=True,length_bp=6,feature=INDEX),Segment('P7 reverse flowcell sequence',revcomp(il.P7))],to='top',pair=('P7 body complement','P7 paired body'),label='source P7bottom',mod5='phosphate',mod3='C3 spacer')
    sc.note('top','Input oligo overlap only; historical kit/genomic ligation ends and full strand assembly unreviewed.')
    return sc

def reagent(name,seq,*,mod5='',mod3=''):
    sc=Scene();sc.strand('oligo',[Segment(name,seq)],label=name,mod5=mod5,mod3=mod3);return sc

def oligos():
    return [oligo('P7top/extension primer',[Segment('P7top',P7_TOP)],five='5′-Me-'),
            oligo('P7bottom',[Segment('P7body',P7_CORE),Segment('source index','N'*6,placeholder=True,length_bp=6,feature=INDEX),Segment('P7tail',revcomp(il.P7))],five='5′-p-',three='-spacerC3-3′'),
            oligo('P5top',[Segment('source P5top includes initialG',P5_TOP)]),oligo('P5bottom',[Segment('P5bottom',P5_BOTTOM)],five='5′-p-'),
            oligo('PCR forward',[Segment('PCR forward',PCR_FORWARD)]),oligo('PCR reverse',[Segment('PCR reverse',PCR_REVERSE)])]
@dataclass(frozen=True)
class APHistory:
    stage:str
    def __post_init__(self):
        if self.stage not in {'input','hips','biotin','p7','captured','cold','hot_eluent','bead_residue','second_bound','supernatant','copied','p5_ligated','library'}:
            raise ValueError('unsupported snAP-seq stage')
    def rows(self):
        sc=Scene()
        if self.stage in {'input','hips','biotin','p7','captured','cold'}:
            for key,label in [('ap','AP-site DNA'),('offtarget','possible competing aldehyde classes such as5-fU;5-fC reaction condition-dependent')]:
                sc.strand(key,[role(label+(';HIPS-conjugated' if self.stage=='hips' else ';HIPS/CuAAC-biotin conjugate' if self.stage not in {'input','hips'} else ';reaction state unassigned'))],label=label)
            if self.stage not in {'cold'}:sc.strand('complement',[role('independent complementary DNA strands; association geometry unassigned')],label='other DNA strands, not covalent appendages')
        elif self.stage=='bead_residue':
            sc.strand('residue',[role('bead-retained biotin/tag-bearing AP-cleavage residue and other retained material; exact chemistry unassembled')],label='physical bead-associated products')
        elif self.stage=='second_bound':
            sc.strand('retained',[role('still-biotinylated off-target DNA and other second-bead-associated material')],label='negative-selection bound fraction, not library substrate')
        elif self.stage in {'hot_eluent','supernatant'}:
            sc.strand('target',[role('AP-cleavage-derived genomic fragment with source P7-derived end; AP/tag region not retained here')],label='source selected DNA',mod5='phosphate from AP-cleavage route')
            sc.strand('background',[role('other eluted/unbound DNA; purity and tag state not quantitative')],label='possible background class')
        elif self.stage=='copied':
            sc.strand('template',[role('AP-cleavage-derived genomic template; source P7 end context unassembled'),role('source P7 sample index; full strand position unassembled',molecular_feature=INDEX)],label='original template retained during copying',mod5='source AP-cleavage5-prime phosphate')
            sc.strand('new_copy',[role('P7-primer-derived complementary DNA; full duplex/end assembly unassigned')],label='newly synthesized DNA strand; no P5 adapter installed yet',mod5='Me(source primer notation)')
        else:
            sc.strand('library',[role('source P5-derived end; full geometry unassembled'),role('AP-boundary-derived genomic DNA'),role('P7-derived end with source6nt index; exact placement unassigned',length_bp=None,molecular_feature=INDEX)],label='source copied/PCR library roles')
            sc.strand('carryover',[role('residual templates/off-target DNA; complete removal not asserted')],label='possible physical carryover')
        if self.stage in {'captured','cold','hot_eluent','supernatant','second_bound','copied','p5_ligated','library'}:
            sc.strand('blocker',[role('possible poly(dI:dC) blocking DNA carryover; partition/recovery unquantified')],label='independent blocking DNA, not joined to genomic sequence')
        return sc.rows()

def wash_beads(g,prefix,label):
    raw=g.input(prefix+'_raw',label,K.MOLECULAR)
    buffer=g.input(prefix+'_wash_buffer','Three source1xBindingBuffer washes; resuspend in source2xbuffer',K.CONDITION)
    ready=g.product(prefix+'_beads','physically prewashed source streptavidin beads',K.SELECTION)
    waste=g.product(prefix+'_wash_waste','removed prewash fractions',K.WASTE)
    g.react(prefix+'_wash','Prewash source streptavidin beads before DNA addition',S.CLEANUP,((raw,R.TEMPLATE),(buffer,R.CONDITION)),(ready,waste));g.discard(waste,'Actual reagent prewash, before genomic DNA addition.')
    return ready

def protocol_graph():
    g=ProtocolGraph(TITLE,authority='Associated PMC6589398 Methods105–115 plus publisherSI custom adapters/native-HeLa extraction; no complete VOR comparison')
    dna=g.input('dna','source purified/sonicated HeLa genomicDNA,near450bp; no glycosylase pretreatment in selected route',K.MOLECULAR,state=APHistory('input'),primary=True)
    hips=g.input('hips','probe1(10mM),40mMNaPi pH7.4,room temperature2h',K.CONDITION)
    tagged=g.product('hips_tagged','source HIPS-treated aldehyde classes',K.MOLECULAR,state=APHistory('hips'))
    g.react('hips_tag','HIPS reaction with source aldehyde DNA classes',S.CHEMICAL_CONVERSION,((dna,R.TEMPLATE),(hips,R.CONDITION)),(tagged,))
    clean=cleanup(g,'hips_clean','Mini-quick spin oligo column purification',tagged,tagged.state)
    click=g.input('click','CuBr250µM/THPTA1.25mM/biotinPEG3azide500µM,37C2h',K.CONDITION)
    biotin=g.product('biotin','CuAAC-biotinylated DNA classes',K.MOLECULAR,state=APHistory('biotin'))
    g.react('biotinylate','CuAAC installs biotin on reacted probe-bearing DNA',S.CHEMICAL_CONVERSION,((clean,R.TEMPLATE),(click,R.CONDITION)),(biotin,))
    clean=cleanup(g,'click_clean','Prewashed10kDaAmicon filter; water/Tris washes and elution',biotin,biotin.state)
    kit=g.input('kit','source NEBNextUltraII preparation; historical end-processing recipe unreviewed',K.CONDITION)
    prepared=g.product('prepared','source kit-prepared labeled DNA for customP7 ligation',K.MOLECULAR,state=clean.state)
    g.react('kit_prepare','Historical kit preparation before customP7 adapter attachment',S.MAGIC,((clean,R.TEMPLATE),(kit,R.CONDITION)),(prepared,),reason='Source names NEBNextUltraII with customP7 substitution; exact historical end-processing and strand products unassembled.')
    p7=g.input('p7_adapter','source custom P7 duplex with six-base TruSeqLT index',K.ADAPTER,state=p7_adapter())
    lig=g.input('p7_ligation','source kit P7 ligation conditions; original kit revision unreviewed',K.CONDITION)
    adapted=g.product('p7_adapted','customP7-bearing labeled DNA',K.MOLECULAR,state=APHistory('p7'))
    g.react('ligate_p7','Attach source customP7 duplex',S.LIGATION,((prepared,R.TEMPLATE),(p7,R.STRUCTURAL),(lig,R.CONDITION)),(adapted,),geometry=AdapterLigationGeometry('unavailable','unavailable',(p7_adapter(),),'Source input oligo overlap known; kit/genomic ligation-ready ends and complete sealed product not reconstructed.'))
    sap=g.input('sap','3U SAP/CutSmart,37C30min',K.CATALYST)
    dephos=g.product('dephosphorylated','SAP-treated preparation before AP-specific hot cleavage',K.MOLECULAR,state=adapted.state)
    g.react('dephosphorylate','Dephosphorylate source preparation',S.DEPHOSPHORYLATION,((adapted,R.TEMPLATE),(sap,R.CATALYST)),(dephos,))
    clean=cleanup(g,'p7_clean','1.4xAMPure purification',dephos,dephos.state)
    beads=wash_beads(g,'first','50µL source Magnesphere streptavidin beads')
    blocker=g.input('first_blocker','2µg poly(dI:dC) added to source ligatedDNA',K.MOLECULAR,state=molecule('blocking DNA',(role('poly(dI:dC) reagent DNA'),)))
    captured=g.product('captured','source first-bead-retained DNA and possible associated blocker',K.MOLECULAR,state=APHistory('captured'))
    wash=g.product('capture_waste','physically removed binding supernatant and6buffer washes',K.WASTE)
    g.react('first_capture','Capture biotin-bearing aldehyde-tagged DNA and wash6times',S.AFFINITY_CAPTURE,((clean,R.TEMPLATE),(blocker,R.SPIKE_IN),(beads,R.SELECTION)),(captured,wash));g.discard(wash,'Source capture and6physical washes; no perfect molecular specificity or blocker depletion assumed.')
    cold=g.input('cold_naoh','100mMNaOH,room temperature10min',K.CONDITION)
    retained=g.product('cold_retained','tagged strands remain bead-associated after cold alkali',K.MOLECULAR,state=APHistory('cold'))
    cold_eluate=g.product('cold_eluate','complementary DNA and other cold-eluted/wash material',K.WASTE)
    g.react('cold_release','Remove complementary strands with room-temperature alkali',S.STRAND_SEPARATION,((captured,R.TEMPLATE),(cold,R.CONDITION)),(retained,cold_eluate));g.discard(cold_eluate,'Source cold-alkali fraction is removed before the hotAP-elution fraction; not a chemical destruction of complementary DNA.')
    ready=cleanup(g,'cold_washes','Additional100mMNaOH wash then3BindingBuffer washes',retained,retained.state)
    hot=g.input('hot_naoh','100mMNaOH70C15min',K.CONDITION)
    eluted=g.product('hot_eluted','source hotAP-cleavage eluate; new5-prime phosphate',K.MOLECULAR,state=APHistory('hot_eluent'))
    residue=g.product('hot_beads','tag-bearing residue and other bead-retained material',K.MOLECULAR,state=APHistory('bead_residue'))
    g.react('hot_release','Hot alkali cleaves tagged AP sites; recover source eluate separately from bead residue',S.CHEMICAL_CONVERSION,((ready,R.TEMPLATE),(hot,R.CONDITION)),(eluted,residue))
    g.mark_terminal(Terminal(residue.id,Disposition.UNRESOLVED,'Original capture beads/tag residues are not the recovered hot eluate; exact final disposal/chemical constitution unassembled.'))
    tris=g.input('tris','Immediately25µL500mMTrisHCl pH7.0',K.CONDITION)
    neutral=g.product('neutral','neutralized hot-cleavage eluate',K.MOLECULAR,state=eluted.state)
    g.react('quench','Neutralize hot eluate immediately',S.QUENCHING,((eluted,R.TEMPLATE),(tris,R.CONDITION)),(neutral,))
    beads=wash_beads(g,'second','75µL fresh source streptavidin beads')
    blocker=g.input('second_blocker','2µg poly(dI:dC) used with fresh capture beads',K.MOLECULAR,state=molecule('fresh blocking DNA',(role('poly(dI:dC) reagent DNA'),)))
    blocked_beads=g.product('second_ready','fresh prewashed beads prepared with blocker',K.SELECTION)
    g.react('prepare_negative','Prepare source negative-selection bead/blocker reagent',S.BINDING,((blocker,R.TEMPLATE),(beads,R.SELECTION)),(blocked_beads,),changes_molecule=False)
    bound=g.product('second_bound','still-biotinylated off-target/other bound material',K.MOLECULAR,state=APHistory('second_bound'))
    supernatant=g.product('supernatant','KEEP unbound AP-cleavage-derived DNA population',K.MOLECULAR,state=APHistory('supernatant'))
    g.react('negative_selection','Second streptavidin capture15min; retain UNBOUND supernatant',S.AFFINITY_CAPTURE,((neutral,R.TEMPLATE),(blocked_beads,R.SELECTION)),(bound,supernatant))
    g.mark_terminal(Terminal(bound.id,Disposition.UNRESOLVED,'Second-bead bound fraction is not the source library substrate; complete background partition and disposal unassigned.'))
    clean=cleanup(g,'ss_clean','Zymo ssDNAclean&concentrator,omitIICcolumn step',supernatant,supernatant.state)
    primer=g.input('p7_primer','source P7top also used as extension primer; Me notation retained',K.OLIGO,state=reagent('P7 primer',P7_TOP,mod5='Me(source notation)'))
    klenow=g.input('klenow','Source P7/dNTP/NEBuffer2:95C1min,65C30s,37C30min THEN add2U Klenow(exo−)',K.CATALYST)
    copied=g.product('copied','source dsDNA synthesis products and possible template/background carryover',K.MOLECULAR,state=APHistory('copied'))
    g.react('copy','P7-primed copying with source Klenow(exo−) ordering',S.EXTENSION,((clean,R.TEMPLATE),(primer,R.STRUCTURAL),(klenow,R.CATALYST)),(copied,))
    clean=cleanup(g,'ds_clean','Zymo DNAclean&concentrator5 purification',copied,copied.state)
    top=g.input('p5_top','source customP5top retains leadingG',K.OLIGO,state=reagent('P5top',P5_TOP))
    bottom=g.input('p5_bottom','source customP5bottom5-phosphate',K.OLIGO,state=reagent('P5bottom',P5_BOTTOM,mod5='phosphate'))
    lig=g.input('p5_ligase','source Blunt/TA ligase master mix,20C30min',K.CATALYST)
    ligated=g.product('p5_ligated','source customP5-bearing DNA; full sealed geometry unassembled',K.MOLECULAR,state=APHistory('p5_ligated'))
    g.react('ligate_p5','Attach customP5 without importing another library geometry',S.LIGATION,((clean,R.TEMPLATE),(top,R.STRUCTURAL),(bottom,R.STRUCTURAL),(lig,R.CATALYST)),(ligated,),geometry=end_boundary('Source P5 strands and copied DNA supplied; full end geometry and sealed strand assembly not reconstructed'))
    clean=cleanup(g,'p5_clean','1.5xAMPure cleanup',ligated,ligated.state)
    forward=g.input('pcr_forward','source library forwardPCR primer',K.OLIGO,state=reagent('ForwardPCR',PCR_FORWARD))
    reverse=g.input('pcr_reverse','source library reversePCR primer',K.OLIGO,state=reagent('ReversePCR',PCR_REVERSE))
    q5=g.input('q5','source Q5HotStart high-fidelity master mix; exact cycles/run unassembled',K.CATALYST)
    library=g.product('library','source amplified library roles',K.LIBRARY,state=APHistory('library'))
    g.react('pcr','Source libraryPCR with verified oligos; full product unassembled',S.PCR,((clean,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(q5,R.CATALYST)),(library,),geometry=PCRGeometry.unavailable(molecule('source library template',(role('complete source strand assembly unreviewed'),)),reason='Exact oligos known locally, but full adapter/template/product geometry and runmask remain unresolved.'))
    g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'MiSeq/NextSeq endpoint source-positive; full product/read/index geometry not assembled.'))
    g.declare_reviewed_scope(('dna','hips','hips_tag','hips_clean_condition','hips_clean_step','click','biotinylate','click_clean_condition','click_clean_step','kit','kit_prepare','p7_adapter','p7_ligation','ligate_p7','sap','dephosphorylate','p7_clean_condition','p7_clean_step','first_raw','first_wash_buffer','first_wash','first_blocker','first_capture','cold_naoh','cold_release','cold_washes_condition','cold_washes_step','hot_naoh','hot_release','tris','quench','second_raw','second_wash_buffer','second_wash','second_blocker','prepare_negative','negative_selection','ss_clean_condition','ss_clean_step','p7_primer','klenow','copy','ds_clean_condition','ds_clean_step','p5_top','p5_bottom','p5_ligase','ligate_p5','p5_clean_condition','p5_clean_step','pcr_forward','pcr_reverse','q5','pcr'),g.authority)
    g.add_source_item(SourceItem('vor_and_kit_assembly','Defining VOR and full historical kit/strand/run assembly','Associated manuscript and publishedSI used; DOI10.1038/s41557-019-0279-9 complete VOR comparison and kit originals remain unreviewed.'))
    return g

def selftest():
    import re
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();source=Source(sources_root()/HERE.parent.name/'s41557-019-0279-9_41557_2019_279_MOESM1_ESM.pdf.txt','PublisherSI Table1 oligos')
    if have(c,source):
        text=re.sub(r'\s+','',source.path.read_text())
        for name,seq in [('P7top',P7_TOP),('P5top-leadingG',P5_TOP),('P5bottom',P5_BOTTOM),('ForwardPCR',PCR_FORWARD),('ReversePCR',PCR_REVERSE)]:c(name+' source transcription',seq in text)
        c('P7bottom6-base index and terminalC3 preserved',P7_CORE+'-NNNNNN-'+revcomp(il.P7)+'spacerC3' in text)
    g=protocol_graph();c('Negative selection keeps source supernatant',g.steps['ss_clean_step'].inputs[0].material,'supernatant')
    c('Cold and hot alkali fractions are separate',g.steps['hot_release'].inputs[0].material,'cold_washes')
    c('Source incomplete geometry blocks finished publication',bool(g.publication_problems()))
    c('Private page lacks finished marker','data-protocol-complete=' not in render_page())
    c.report()
