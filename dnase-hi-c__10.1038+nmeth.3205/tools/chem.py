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

PAGE_FILE='dnase-hi-c.html'
TITLE='DNase Hi-C — biotin bridges and whole-genome in-gel ligation'
NOTES='01_dnase-hi-c.html'
SOURCE='<a href="https://doi.org/10.1038/nmeth.3205">Ma et al. (online2014)</a>; associated accepted manuscriptPMC4281301 Online Methods114–158; publisher Supplementary Note2; source Table23.'
CAVEAT='Selected K562 whole-genome route. H1 fixation and targeted SeqCap/pre-capture PCR are separate source branches. Genomic loci, ligation yields and full bridge/genomic duplex chemistry are unassigned. Blunt-blocker DNA, dangling/unligated DNA and reaction products do not vanish at enzymatic treatment. Historical library end chemistry, exact source index member and full PCR/run assembly remain unreviewed. Known input adapter pairing is not validation of a chromatin-contact product.'
READOUT='Source HiSeq2000 paired-end150bp,9–11-cycle KAPA libraryPCR and eight-base source sample indexes. Exact strand/primer/index-cycle mapping remains unassembled; no other Hi-C adapter or restriction junction is imported.'
PREVIEWS=(('Whole-genome K562 route','library'),('Prepare the streptavidin capture reagent','washed_beads'))
BIO_TOP='GCTGAGGGATC'
BIO_BOTTOM='CCTCAGCT'
BLOCK_TOP='GCTGAGGGAC'
BLOCK_BOTTOM='CCTCAGC'
INDEX=feature('dnase_hic_sample_index','sample_index','unknown')
def bridge(blocker=False):
    core=revcomp(BLOCK_BOTTOM)
    top=BLOCK_TOP if blocker else BIO_TOP
    bottom=BLOCK_BOTTOM if blocker else BIO_BOTTOM
    if not top.startswith(core) or not bottom.startswith(BLOCK_BOTTOM):
        raise ValueError('source bridge/blocker lacks supplied antiparallel core')
    sc=Scene();parts=[Segment('paired core',core)]
    if blocker:
        parts.append(Segment('source blocker exposed end',top[len(core):]))
    else:
        tail=top[len(core):]
        if len(tail)<2 or tail[-2]!='T':
            raise ValueError('source internal biotin dT position must remain a thymidine')
        parts.extend([Segment('bridge exposed prefix',tail[:-2]),Segment('internal biotin dT',tail[-2]),Segment('bridge final base',tail[-1])])
    sc.strand('top',parts,label='nonbiotin blunt blocker' if blocker else 'biotin T-bridge top',mod5='' if blocker else 'phosphate')
    paired=[Segment('paired core complement',BLOCK_BOTTOM)]
    if not blocker:paired.append(Segment('source3-prime T overhang',bottom[len(BLOCK_BOTTOM):]))
    sc.anneal('bottom',paired,to='top',pair=('paired core complement','paired core'),label='source bottom oligo')
    sc.note('top','Computed input overlap only; genomic/proximity junction remains unassembled.' + ('' if blocker else ' Biotin is on the internal dT, not the5-prime terminus.'))
    return sc

def oligos():
    return [oligo('Bridge top',[Segment('core plus GA',BIO_TOP[:9]),Segment('biotinylated T','/iBiodT/',placeholder=True),Segment('last C','C')],five='5′-/5Phos/'),
            oligo('Bridge bottom',[Segment('Bridge bottom',BIO_BOTTOM)]),
            oligo('Blunt blocker top',[Segment('Blocker top',BLOCK_TOP)]),oligo('Blunt blocker bottom',[Segment('Blocker bottom',BLOCK_BOTTOM)])]
@dataclass(frozen=True)
class ChromatinState:
    stage:str
    def __post_init__(self):
        if self.stage not in {'input','fixed','digested','a_tail','bridged','joined','deproteinized','recovered','dangling','sheared','library'}:
            raise ValueError('unsupported DNase Hi-C stage')
    def rows(self):
        sc=Scene()
        if self.stage in {'input','fixed','digested','a_tail'}:
            if self.stage=='a_tail':
                sc.strand('dna',[role('genomic body; full opposite strand unassembled'),Segment('source terminal A','A')],label='source dA-tailed end class')
                sc.strand('residual_blunt',[role('residual blunt genomic fragments; source explicitly supplies blocking adapter')],label='separate non-A-tailed class; no population fraction assigned')
            else:
                sc.strand('dna',[role({'input':'chromatin genomic DNA','fixed':'formaldehyde-crosslinked genomic DNA','digested':'DNase genomic fragment; cut coordinates unassigned'}[self.stage])],label='source genomic population',mod5='phosphate' if self.stage=='digested' else '',mod3='hydroxyl' if self.stage=='digested' else '')
        elif self.stage=='bridged':
            sc.strand('end_tagged',[role('one genomic fragment with source T-bridge at its end; no proximity join yet')],label='biotin bridge-bearing genomic end class')
            sc.strand('blocked',[role('blunt blocker-associated DNA and unligated DNA; exact end classes unassigned')],label='other physical DNA end classes')
        else:
            sc.strand('contact',[role('genomic arm A'),role('source bridge-associated junction; exact strand assembly unreviewed'),role('genomic arm B')],label='source contact-product class; no assigned pair/yield')
            sc.strand('other_dna',[role('dangling/unligated/nonbiotin-blocker-associated DNA remains; exact enrichment and terminal labels unassigned')],label='other physical DNA classes, not automatically discarded')
        if self.stage in {'fixed','digested','a_tail','bridged','joined'}:
            sc.strand('associated',[role('source protein/chromatin association and preparation components; not joined as DNA')],label='independent associated components')
        if self.stage=='deproteinized':
            sc.strand('debris',[role('protein/agarose digestion products remain until actual workup')],label='separate reaction debris')
        if self.stage=='library':
            sc.strand('index',[role('source8-base sample index; PCR member and strand position unassembled',length_bp=8,molecular_feature=INDEX)],label='identifier role, not an assigned covalent library position')
        sc.note('dna' if self.stage in {'input','fixed','digested','a_tail'} else 'end_tagged' if self.stage=='bridged' else 'contact','Source population roles only. No complete fragmentation, ligation, background depletion or quantitative recovery asserted.')
        return sc.rows()

def treatment(g,ident,label,kind,template,reagent_label,reagent_kind,reagent_role,state):
    reagent=g.input(ident+'_reagent',reagent_label,reagent_kind)
    out=g.product(ident+'_product',label,K.MOLECULAR,state=state)
    g.react(ident,label,kind,((template,R.TEMPLATE),(reagent,reagent_role)),(out,))
    return out

def protocol_graph():
    g=ProtocolGraph(TITLE,authority='Associated PMC4281301 Online Methods114–158 and publisher Supplementary Note2/Table23; whole-genome K562 route')
    cells=g.input('cells','K562 chromatin',K.MOLECULAR,state=ChromatinState('input'),primary=True)
    fixed=treatment(g,'fix','Fix K562 with1% formaldehyde,10min room temperature',S.CROSS_LINKING,cells,'1% formaldehyde',K.CONDITION,R.CONDITION,ChromatinState('fixed'))
    quenched=treatment(g,'quench','Quench10min with0.125M glycine',S.QUENCHING,fixed,'0.125M glycine',K.CONDITION,R.CONDITION,fixed.state)
    washed=cleanup(g,'cell_wash','PBS wash; cold lysis bufferI10min; source freeze/thaw',quenched,quenched.state)
    sds=treatment(g,'sds','SDS/TE lysis,37C10min',S.THERMAL_TREATMENT,washed,'TE lysis buffer with1%SDS',K.CONDITION,R.CONDITION,washed.state)
    nuclei=cleanup(g,'nuclei_wash','Collect nuclei;0.15%Igepal wash; resuspend in DNase/RNase buffer',sds,sds.state)
    rnase=treatment(g,'rnase_buffer','Source RNase-containing DNase buffer,37C10min; RNA products remain until bead workup',S.DIGESTION,nuclei,'0.5xDNase buffer/MnCl2/CaCl2 with RNaseA',K.CATALYST,R.CATALYST,nuclei.state)
    digested=treatment(g,'dnase','DNaseI5min25C; genomic cut coordinates unassigned',S.DIGESTION,rnase,'source0.015U/µL DNaseI',K.CATALYST,R.CATALYST,ChromatinState('digested'))
    stopped=treatment(g,'stop_dnase','Stop with source SDS/EDTA solution',S.QUENCHING,digested,'10xstop:10%SDS/250mMEDTA,one-tenth volume',K.CONDITION,R.CONDITION,digested.state)
    ready=cleanup(g,'digest_clean','2xAMPure recovery; two80%ethanol washes; distribute3tubes',stopped,stopped.state)
    repaired=treatment(g,'repair','End-repair1h room temperature',S.END_REPAIR,ready,'T4DNA polymerase/Klenow/dNTPs in source ligation buffer',K.CATALYST,R.CATALYST,ready.state)
    ready=cleanup(g,'repair_clean','SDS stop; equal volume20%PEG/2.5MNaCl; bead recovery and80%ethanol washes',repaired,repaired.state)
    tailed=treatment(g,'a_tail','Klenow(exo−)/dATP37C1h',S.TAILING,ready,'source Klenow(exo−)/dATP in NEBbuffer2',K.CATALYST,R.CATALYST,ChromatinState('a_tail'))
    ready=cleanup(g,'tail_clean','SDS stop and source bead workup',tailed,tailed.state)
    bio=g.input('biotin_bridge','10µM source internally biotinylated T-bridge',K.ADAPTER,state=bridge())
    block=g.input('blunt_blocker','20µM source nonbiotin blunt blocker; not an index',K.ADAPTER,state=bridge(True))
    lig=g.input('bridge_ligase','25U T4DNA ligase in Rapid ligation buffer;1h room temperature',K.CATALYST)
    bridged=g.product('bridged','bridge-bearing and blocker/other genomic populations',K.MOLECULAR,state=ChromatinState('bridged'))
    g.react('bridge_ligation','Ligate T-bridge; source blunt blocker handles residual non-A-tailed fragments',S.LIGATION,((ready,R.TEMPLATE),(bio,R.STRUCTURAL),(block,R.STRUCTURAL),(lig,R.CATALYST)),(bridged,),geometry=AdapterLigationGeometry('unavailable','unavailable',(bridge(),bridge(True)),'Exact source input oligos paired; genomic ligation-ready duplexes, sealing pattern and full strand-resolved product not assembled.'))
    ready=cleanup(g,'bridge_clean','SDS stop and source PEG/NaCl bead recovery',bridged,bridged.state)
    phosphorylated=treatment(g,'phosphorylate','T4PNK37C1h',S.PHOSPHORYLATION,ready,'source T4PNK in ligation buffer',K.CATALYST,R.CATALYST,ready.state)
    gel=g.input('gel','0.4%low-melting agarose in source T4ligation buffer; solidify on ice',K.CONDITION)
    ligase=g.input('proximity_ligase','200U T4DNA ligase;16C4h then25C1h',K.CATALYST)
    joined=g.product('joined','source in-gel contact-product population and other DNA',K.MOLECULAR,state=ChromatinState('joined'))
    g.react('proximity','In-gel proximity ligation; no exact genomic pairing/yield assigned',S.PROXIMITY_LIGATION,((phosphorylated,R.TEMPLATE),(gel,R.CONDITION),(ligase,R.CATALYST)),(joined,))
    melted=treatment(g,'melt','Melt agarose70C10min',S.THERMAL_TREATMENT,joined,'70C10min',K.CONDITION,R.CONDITION,joined.state)
    agarase=treatment(g,'agarase','Preincubate42C10min; agarase overnight',S.DIGESTION,melted,'source agarase; polymer digestion does not physically remove DNA',K.CATALYST,R.CATALYST,melted.state)
    concentrated=cleanup(g,'concentrate','AmiconUltra15,30kDa concentration to near1mL',agarase,agarase.state)
    deproteinized=treatment(g,'proteinase','ProteinaseK/SDS65C overnight; morePK55C2h',S.DEPROTEINIZATION,concentrated,'source ProteinaseK/SDS; de-crosslinking conditions',K.CATALYST,R.CATALYST,ChromatinState('deproteinized'))
    recovered=cleanup(g,'dna_recovery','GlycoBlue/Na-acetate/isopropanol precipitation;QIAquickGelExtraction recovery',deproteinized,ChromatinState('recovered'))
    dangling=treatment(g,'dangling_treatment','Source T4DNA polymerase/dATP+dGTP,25C30min',S.END_REPAIR,recovered,'50U T4DNApolymerase,0.2mM each dATP/dGTP',K.CATALYST,R.CATALYST,ChromatinState('dangling'))
    stopped=treatment(g,'stop_dangling','EDTA stop',S.QUENCHING,dangling,'25µL0.5MEDTA',K.CONDITION,R.CONDITION,dangling.state)
    concentrated=cleanup(g,'dangling_concentrate','AmiconUltra0.5,30kDa concentrate; source TE/SDS buffer',stopped,stopped.state)
    sheared=treatment(g,'covaris','CovarisS2 shear to100–300bp',S.FRAGMENTATION,concentrated,'source Covaris parameters; resulting length distribution, not an exact insert',K.CONDITION,R.CONDITION,ChromatinState('sheared'))
    recovered=cleanup(g,'shear_clean','AMPure recovery and two80%ethanol washes',sheared,sheared.state)
    repaired=treatment(g,'library_repair','Library end repair16C10min',S.END_REPAIR,recovered,'Fermentas EndRepair mix',K.CATALYST,R.CATALYST,recovered.state)
    recovered=cleanup(g,'library_repair_clean','0.5%SDS stop and source bead recovery',repaired,repaired.state)
    tailed=treatment(g,'library_a_tail','Library Klenow(exo−)/dATP37C1h',S.TAILING,recovered,'source Klenow(exo−),2mMdATP',K.CATALYST,R.CATALYST,recovered.state)
    recovered=cleanup(g,'library_tail_clean','0.5%SDS stop and source bead recovery',tailed,tailed.state)
    adapter=g.input('library_adapter','3µM Table23 T-tailed sequencing adapter; full duplex unassembled',K.ADAPTER,state=molecule('sequencing adapter roles',(role('Table23 adapter duplex; full geometry unassembled'),)))
    condition=g.input('library_ligation_condition','source Rapid ligation buffer,room temperature30min; historical formulation unreviewed',K.CONDITION)
    adapted=g.product('library_adapted','sequencing-adapter-bearing source DNA classes',K.MOLECULAR,state=recovered.state)
    g.react('library_ligation','Ligate source sequencing adapters',S.LIGATION,((recovered,R.TEMPLATE),(adapter,R.STRUCTURAL),(condition,R.CONDITION)),(adapted,),geometry=end_boundary('Source T-tailed library adapter and A-tailed DNA; full Table23 duplex/junction unassembled'))
    recovered=cleanup(g,'library_adapter_clean','0.5%SDS stop/bead recovery; elute200µLwater',adapted,adapted.state)
    beads=g.input('raw_beads','25µL MyOneC1 streptavidin beads before source washing',K.MOLECULAR)
    bw=g.input('bead_wash_condition','B&W wash,1xBSA wash,1xB&W wash; resuspend200µL2xB&W',K.CONDITION)
    washed_beads=g.product('washed_beads','physically washed streptavidin beads',K.SELECTION)
    bead_waste=g.product('bead_waste','removed bead-preparation wash fractions',K.WASTE)
    g.react('wash_beads','Prepare capture beads by source physical washes',S.CLEANUP,((beads,R.TEMPLATE),(bw,R.CONDITION)),(washed_beads,bead_waste));g.discard(bead_waste,'Source bead-preparation wash occurs before adding DNA.')
    captured=g.product('captured','biotin-enriched source DNA classes; background not assumed zero',K.MOLECULAR,state=recovered.state)
    flow=g.product('capture_flow','physically removed binding supernatant',K.WASTE)
    g.react('capture','Bind biotin-bearing DNA25C20min and remove supernatant',S.AFFINITY_CAPTURE,((recovered,R.TEMPLATE),(washed_beads,R.SELECTION)),(captured,flow));g.discard(flow,'Source physically removes capture supernatant; no exact class purity assigned.')
    selected=cleanup(g,'capture_wash','Source B&W/TE transfer;2xB&W,NEBuffer2,EB washes;resuspend80µLEB',captured,captured.state)
    forward=g.input('forward','Table23 universal PCR primer; full library landing geometry unassembled',K.OLIGO,state=molecule('forward PCR role',(role('source universal PCR primer; assembly deferred'),)))
    reverse=g.input('reverse','Table23 indexed reverse PCR primer; actual sample member unassigned',K.OLIGO,state=molecule('indexed reverse PCR role',(role('source indexed reverse primer; sample-code boundaries/strand assembly unassigned',molecular_feature=INDEX),)))
    kapa=g.input('kapa','KAPA HiFiHotStart;whole-genome9–11cycles',K.CATALYST)
    amplified=g.product('amplified','whole-genome PCR products and possible background classes',K.MOLECULAR,state=ChromatinState('library'))
    g.react('pcr','Whole-genome9–11-cycle PCR; no SeqCap branch merged',S.PCR,((selected,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(kapa,R.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(molecule('source enriched library',(role('exact post-library adapter/template assembly unreviewed'),)),reason='Table23 primers known as source reagents; complete model assembly and source member assignment remain unresolved.'))
    library=cleanup(g,'library','Pool source PCR products;AMPure size fractionation',amplified,amplified.state)
    g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'Source HiSeq2000 PE150 endpoint; exact molecule, sample member and run-primer/index cycle geometry unassembled.'))
    g.declare_reviewed_scope(('cells','fix_reagent','fix','quench_reagent','quench','cell_wash_condition','cell_wash_step','sds_reagent','sds','nuclei_wash_condition','nuclei_wash_step','rnase_buffer_reagent','rnase_buffer','dnase_reagent','dnase','stop_dnase_reagent','stop_dnase','digest_clean_condition','digest_clean_step','repair_reagent','repair','repair_clean_condition','repair_clean_step','a_tail_reagent','a_tail','tail_clean_condition','tail_clean_step','biotin_bridge','blunt_blocker','bridge_ligase','bridge_ligation','bridge_clean_condition','bridge_clean_step','phosphorylate_reagent','phosphorylate','gel','proximity_ligase','proximity','melt_reagent','melt','agarase_reagent','agarase','concentrate_condition','concentrate_step','proteinase_reagent','proteinase','dna_recovery_condition','dna_recovery_step','dangling_treatment_reagent','dangling_treatment','stop_dangling_reagent','stop_dangling','dangling_concentrate_condition','dangling_concentrate_step','covaris_reagent','covaris','shear_clean_condition','shear_clean_step','library_repair_reagent','library_repair','library_repair_clean_condition','library_repair_clean_step','library_a_tail_reagent','library_a_tail','library_tail_clean_condition','library_tail_clean_step','library_adapter','library_ligation_condition','library_ligation','library_adapter_clean_condition','library_adapter_clean_step','raw_beads','bead_wash_condition','wash_beads','capture','capture_wash_condition','capture_wash_step','forward','reverse','kapa','pcr','library_condition','library_step'),g.authority)
    g.add_source_item(SourceItem('unreviewed_vor_and_assembly','Complete defining VOR comparison and post-ligation strand assembly','Defining DOI10.1038/nmeth.3205; associated manuscript used, full VOR comparison/end product validation not completed.'))
    return g

def selftest():
    import re
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();source=Source(sources_root()/HERE.parent.name/'nmeth.3205_supp_NIHMS641840-supplement-S23.xlsx.txt','Associated Table23 bridge/adaptor originals')
    if have(c,source):
        text=re.sub(r'\s+','',source.path.read_text())
        for label,written in [('biotin bridge','/5Phos/GCTGAGGGA/iBiodT/C'),('bridge bottom',BIO_BOTTOM),('blunt blocker top',BLOCK_TOP),('blunt blocker bottom',BLOCK_BOTTOM)]:
            c(label+' complete source transcription',written in text)
    g=protocol_graph();c('Private source boundary blocks finished publication',bool(g.publication_problems()))
    c('No finished marker','data-protocol-complete=' not in render_page())
    c('Capture uses physically washed beads',[port.material for port in g.steps['capture'].inputs if port.role is R.SELECTION],['washed_beads'])
    c.report()
