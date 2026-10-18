"""Nano3P source-owned starter duplex, native3′ capture and AMX boundary."""
from chemdraw import Construct,Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import Disposition,Terminal,InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,Highlight,AnnotationSlot
from protocol_model import molecule,role
from nano3p_chem import R_RNA,D_DNA,COMPA_DNA,final_library,oligos
TITLE='Nano3P-seq — Nanopore3′-end capture'
NOTES='01_nano3p_seq.html'
SOURCE='Begik et al., Nature Methods(2023), <a href="https://doi.org/10.1038/s41592-022-01714-w">doi:10.1038/s41592-022-01714-w</a>; published Methods and SupplementaryProtocol.'
SUMMARY='TGIRT starts on an RNA/DNA duplex and template-switches onto native RNA3′ends; CompA provides a local duplex/A-overhang end for PCR-free AMX ligation.'
CAVEAT='Selected Methods route uses50ngRNA/1µLTGIRT and100µL AMX ligation. Source100ng/2µLTGIRT and smaller supplementary ligation are distinct alternatives. Native RNA and blocked starterRNA persist through RT until RNase; cDNA remains single stranded apart from CompA-bound handle. The SQK-DCS109 AMX complementary end/motor geometry remains unassigned; no RNA-adapter ligation orPCR is inserted.'
FINAL_LIBRARY=final_library()
SEQ_PRIMERS=()
FINAL_CAPTION='PCR-free cDNA contains the actual D_DNA handle and RNA-derived copy. AMX joining details remain unassigned.'
SEQUENCING_MECHANISM='nanopore_motor'
MOTOR_ADAPTER_NAME='SQK-DCS109 AMX motor adapter'
READ_LENGTHS=None
SEQUENCING_INTRO='Source uses R9.4.1 Flongle orMinION flowcells, MinKNOW3.5.5; no fixed read-cycle plan.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE='The AMX physical adapter/motor end is source-named but not reconstructed.'

def dna_parts(*,copied=False):
    p=[Segment('D_DNA core',D_DNA[:-1]),Segment('starter3′N','N',placeholder=True,length_bp=1,note='overhang/capture base; not UMI')]
    if copied:p.append(role('RNA-derived cDNA after matched terminalN; nativeRNA3′→5′'))
    return p

def starter_scene():
    sc=Scene();sc.strand('RNA',[Segment('R_RNA starter',R_RNA)],label='source RNA starter',mod3='3SpC3')
    sc.anneal('DNA',dna_parts(),to='RNA',pair=('D_DNA core','R_RNA starter'),unpaired=('starter3′N',),label='5′phosphorylated DNA starter',mod5='P')
    return sc

def cdna_scene():
    sc=Scene();sc.strand('cDNA',dna_parts(copied=True),label='tracked first-strand cDNA target',mod5='P',mod3='OH');return sc

def ready_scene():
    sc=cdna_scene()
    sc.anneal('CompA',[Segment('CompA paired core',COMPA_DNA[:-1]),Segment('CompA3′A overhang',COMPA_DNA[-1])],to='cDNA',pair=('CompA paired core','D_DNA core'),label='source CompA local complement')
    sc.overhang('CompA','CompA3′A overhang',polarity='3-prime',owner='substrate',text='reported CompA3′A end')
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input('rna','native50ngRNA; selected source baseline',MaterialKind.BIOLOGICAL,state=molecule('nativeRNA',(role('RNA with native3′end; tail status variable'),)),primary=True)
    starter=g.input('rt_starter','preannealed3′C3-blocked R_RNA and5′P D_DNA;94°C1min→25°C',MaterialKind.OLIGO,state=starter_scene(),source='Published Methods243/SuppProtocol oligo table')
    rt=g.input('reverse_transcriptase','TGIRT1µL;RT30min beforedNTP,60°C60min,75°C15min',MaterialKind.CATALYST)
    cdna=g.product('full_cdna','single-strand template-switched cDNA target',MaterialKind.MOLECULAR,state=cdna_scene())
    templates=Scene();templates.strand('nativeRNA',[role('nativeRNA templates')],label='independent nativeRNA');templates.strand('starterRNA',[Segment('R_RNA starter',R_RNA)],label='independent starterRNA',mod3='3SpC3')
    residual=g.product('rt_template_rna','native/starter RNA after RT',MaterialKind.MOLECULAR,state=templates)
    g.react('template_switch_rt','TGIRT switches from starter duplex to nativeRNA3′end',StepKind.TEMPLATE_SWITCHING,((rna,InputRole.TEMPLATE),(starter,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,residual))
    rnase=g.input('rnase_cocktail','RNaseCocktail AM2286;37°C10min',MaterialKind.CATALYST)
    digested=g.product('digested_cdna','cDNA target after RNase',MaterialKind.MOLECULAR,state=cdna.state)
    fragments=g.product('rna_fragments','hydrolysed native/starter RNA products',MaterialKind.MOLECULAR)
    g.react('digest_rna','Digest native and starterRNA',StepKind.DIGESTION,((cdna,InputRole.TEMPLATE),(residual,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),(digested,fragments))
    g.mark_terminal(Terminal(fragments.id,Disposition.UNRESOLVED,'RNase produces fragments, not physical removal. Individual recovery through the following0.8×bead cleanup is unassigned; graph follows cDNA targets.'))
    beads=g.input('first_cleanup','0.8×AMPureXP;physical supernatant/70%ethanol removal;16µLwater elution',MaterialKind.SELECTION)
    purified=g.product('purified_cdna','purified cDNA target component',MaterialKind.MOLECULAR,state=digested.state)
    wash=g.product('first_cleanup_wash','excluded small molecules/adapters and ethanol washes',MaterialKind.WASTE)
    g.react('purify_cdna','Purify cDNA targets on0.8×AMPureXP',StepKind.CLEANUP,((digested,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(purified,wash),changes_molecule=False);g.discard(wash,'physical magnet/supernatant and ethanol removals in source')
    comp=g.input('compa','CompA_DNA1µL10µM;Methods94°C1min→25°C',MaterialKind.OLIGO,state=COMPA_DNA,source='Published Methods243/SuppProtocol')
    ready=g.product('ligation_ready','partially duplex CompA-bound cDNA;3′A on CompA',MaterialKind.MOLECULAR,state=ready_scene())
    g.react('anneal_compa','Anneal CompA to actual cDNA handle',StepKind.HYBRIDIZATION,((purified,InputRole.TEMPLATE),(comp,InputRole.STRUCTURAL)),(ready,))
    amx=g.input('motor_adapter','SQK-DCS109 AMX motor adapter; physical end unassigned',MaterialKind.ADAPTER)
    ligase=g.input('adapter_ligase','Blunt/TA LigaseMix50µL;22.5µLcDNA,5µLAMX,22.5µLwater;RT10min',MaterialKind.CATALYST)
    ligated=g.product('amx_ligated','AMX-ligated cDNA target',MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
    g.react('ligate_motor','Ligate AMX withoutPCR',StepKind.LIGATION,((ready,InputRole.TEMPLATE),(amx,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.unavailable('Actual CompA partial duplex/A-overhang is explicit, but SQK-DCS109 AMX complementary end and complete joining geometry are unassigned; no compatible end is borrowed from newer ONT kits.'))
    post=g.input('motor_cleanup','0.8×AMPureXP;WSB washes,elutionbuffer',MaterialKind.SELECTION)
    lib=g.product('library','purified PCR-free Nano3P library',MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True)
    waste=g.product('motor_cleanup_waste','WSB/supernatant and unbound adapter material',MaterialKind.WASTE)
    g.react('motor_cleanup_step','Purify AMX-ligated targets before loading',StepKind.CLEANUP,((ligated,InputRole.TEMPLATE),(post,InputRole.SELECTION)),(lib,waste),changes_molecule=False);g.discard(waste,'source physical postligation bead/WSB wash removal')
    g.sequence(lib,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_INTRO+' '+SEQUENCING_UNAVAILABLE))
    g.annotate(Highlight('ligate_motor','Published condition alternatives','This graph selects article Methods100µLligation. SupplementaryProtocol lists22.5µLcDNA+2.5µLAMX+25µLligasemix and90°CCompA anneal; article lists100µLligation and94°Canneal. These are retained as separate stated recipes, not silently merged.',slot=AnnotationSlot.AFTER_WORKFLOW))
    g.declare_reviewed_scope(('rna','rt_starter','reverse_transcriptase','template_switch_rt','rnase_cocktail','digest_rna','first_cleanup','purify_cdna','compa','anneal_compa','motor_adapter','adapter_ligase','ligate_motor','motor_cleanup','motor_cleanup_step','library'),SOURCE)
    return g.require_publishable()
