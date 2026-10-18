"""HiRES joint Hi-C and RNA library, as described by the author protocol."""
from pcr import PCRGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp, complement_segments
from protocol_model import role
from protocol_graph import Terminal, Disposition
import seqprimers as sp
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
TITLE = 'HiRES — single-cell Hi-C and RNA-seq'
NOTES = '01_hires.html'
SOURCE = 'Defining <a href="https://doi.org/10.1126/science.adg3797">Science paper</a>; <a href="https://doi.org/10.21769/BioProtoc.4886">author protocol</a>; authors\' <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223917">GEO GSE223917</a> sample records.'
SUMMARY = 'The same fixed nucleus is reverse-transcribed, MboI-digested and proximity-ligated, then sorted into its own well. A joint amplification feeds tagmentation, an RNA-P7/i5 enrichment and i7 indexing PCR.'
CAVEAT = "The defining Science article and supplement were not inspected in this capped pass; this page follows the authors' Bio-protocol (doi:10.21769/BioProtoc.4886). The authors' GEO GSE223917 sample records repeat the tagmentation, two-stage index PCR, RNA-P7 sequence and paired 150-base NovaSeq 6000 reads, and state that RNA reads are recognized by the GAT5-RT arm at the 5′ end of Read 2 and deduplicated by UMI. The author protocol prints four HiRES primers and specifies per-well i5/i7 indexing. The TruePrep manual prints a standard finished-library shell, but its loaded adapter oligos, the HiRES RNA-P7 junction and selected index-member bases remain unresolved. The two views below are read classes from one library, not separate index-ligation rounds."
UMI_FEATURE = feature('hires_umi', 'umi', 'random', note='GAT5-RT N8; sourceGEO deduplication association, not GAT5-7N priming window')
I5_FEATURE = feature('hires_i5', 'sample_index', 'unknown', note='well-specificUD primer; member/physicalwidth unassigned')
I7_FEATURE = feature('hires_i7', 'sample_index', 'unknown', note='well-specificUD primer; member/physicalwidth unassigned')
CONTACT = Construct([role('native TruePrep/PCR endA'), role('MboI-derived DNAcontact insert'), role('native i5 endcontext', molecular_feature=I5_FEATURE), role('native i7 endcontext', molecular_feature=I7_FEATURE)], name='HiRES contact read-class reference')
RNA = Construct([role('native RNA-P7/PCR endA'), role('GAT5-RT N8UMI lineage', length_bp=8, molecular_feature=UMI_FEATURE), role('RNA-derived insert'), role('native i5 endcontext', molecular_feature=I5_FEATURE), role('native i7 endcontext', molecular_feature=I7_FEATURE)], name='HiRES transcript read-class reference')
FINAL_LIBRARY = None
FINAL_LIBRARIES = (('Contact read-class reference', CONTACT, (), 'Computational class from one physicaljointlibrary; native strandassembly unassigned.', 'PE150; native index acquisition unassigned.'), ('Transcript read-class reference', RNA, (), 'Computational class from samejointlibrary; native RNA-P7junction unassigned.', 'RNA reads recognized by sourceGAT5RT arm nearRead2; RNARead1aligned.'))
SEQ_PRIMERS = ()
SHOW_READ_CONTENT = False
READ_SPECS = (sp.read_spec('Read 1', 150, 'contact or RNA-derived sequence; class-specific mapping', 'include'), sp.read_spec('Read 2', 150, 'contact mapping or GAT5RTarm/UMI-based RNAclassrecognition; nativecyclepositions unassigned', 'extract_then_include'), sp.read_spec('Index 1 (i7)', 'unassigned', 'WellUD member/physicalwidth/acquisition unassigned', 'undetermined'), sp.read_spec('Index 2 (i5)', 'unassigned', 'WellUD member/physicalwidth/acquisition unassigned', 'undetermined'))
SEQUENCING_ENDING = 'BioProtocol/GSM6998595: NovaSeq6000PE150; DNA/RNA are computationalclasses of onejointphysicallibrary. NativeUDindex acquisition remains unassigned.'
GAT5_RT = 'GTAGGTGTGAGTGATGGTTGAGGTAGTATTGCGCAATGNNNNNNNNTTTTTTTTTTTTTTTVN'
GAT5_7N = 'GTAGGTGTGAGTGATGGTTGAGGTAGTNNNNNNN'
GAT5 = 'GTAGGTGTGAGTGATGGTTGAGGTAGT'
RNA_P7 = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGTTGAGGTAGTATTGCGCAATG'

def oligos():
    return [f'<b>GAT5-RT</b> — 5′-{GAT5_RT}-3′', f'<b>GAT5-7N</b> — 5′-{GAT5_7N}-3′', f'<b>GAT5</b> — 5′-{GAT5}-3′', f'<b>RNA-P7</b> — 5′-{RNA_P7}-3′']

def _oligo(name, seq):
    if name == 'GAT5-RT':
        return Construct(rt_parts(), name=name)
    if name == 'GAT5-7N':
        return Construct([Segment('GAT5common', GAT5), Segment('random7priming tract', 'N' * 7, placeholder=True, note='random priming, not an assigned UMI')], name=name)
    return Construct([Segment(name, seq)], name=name)

def rt_parts():
    n = GAT5_RT.index('N')
    tail = GAT5_RT[n + 8:]
    return [Segment('GAT5RT fixedhandle', GAT5_RT[:n]), Segment('GAT5RT UMI', 'N' * 8, placeholder=True, feature=UMI_FEATURE), Segment('RT polyT tract', tail[:-2]), Segment('VN anchor', 'VN', placeholder=True)]

def mixture(stage):
    sc = Scene()
    if stage in ('raw', 'fixed'):
        sc.strand('DNA', [role('nuclear genomicDNA; native modifications retained')], label='independent DNA')
        sc.strand('RNA', [role('nuclear RNA; physicalU')], label='independent RNA')
        return sc
    dna = [role('native genomicDNA/contact extent unassigned')]
    if stage in ('digested',):
        dna = [Segment('MboI5′GATC', 'GATC'), role('restrictionfragmentbody; independentopposite fragments unassigned')]
    elif stage == 'copied':
        dna = [role('original nuclear genomicDNA; not yet digested or randomly primed')]
    elif stage in ('ligated', 'lysed'):
        dna = [role('locusA'), Segment('unfilled MboI religation junction', 'GATC'), role('locusB')]
    else:
        dna = [role('native random-primer-derived contactDNA; full strand/endassembly unassigned')]

    if stage in ('enriched','indexed','pooled'): dna.append(role('native i5-bearing DNA/carryover classes',molecular_feature=I5_FEATURE))
    if stage in ('indexed','pooled'): dna.append(role('native i7-bearing DNA/carryover classes',molecular_feature=I7_FEATURE))
    sc.strand('DNA', dna, label='actual DNA/contact lineage; ' + stage)
    if stage == 'digested':
        opposite=complement_segments([dna[1]])
        sc.anneal('DNA opposite',opposite,to='DNA',pair=(opposite[0].name,dna[1].name),unpaired=(opposite[0].name,),mod5='p',label='opposite restriction strand; native extents unassigned')
        sc.overhang('DNA',dna[0].name,polarity='5-prime',owner='substrate',text='sourceMboI5′GATC cohesive end')
    cdna = [*rt_parts(), role('antisense transcript cDNA')]
    if stage in ('tagged', 'extended', 'enriched', 'indexed', 'pooled'):
        cdna.append(role('native TruePrep/RNA-P7extension footprint; gap/strand/retention unassigned'))
    if stage in ('enriched', 'indexed', 'pooled'):
        cdna.append(role('installed i5 primer lineage; member/extent unassigned', molecular_feature=I5_FEATURE))
    if stage in ('indexed', 'pooled'):
        cdna.append(role('installed i7 primer lineage; member/extent unassigned', molecular_feature=I7_FEATURE))
    sc.strand('cDNA', cdna, label='independent RNAderivedDNA/copy population; ' + stage)
    sc.strand('originals', [role('possible original RNA, genomic/cDNAtemplates, randomprimers and nonproductive tag products; recovery unresolved')], label='independent residualorigins; copying/nonamplification not disposal')
    return sc

def protocol_graph():
    source = 'Chen et al. 2023: https://doi.org/10.21769/BioProtoc.4886'
    g = ProtocolGraph(TITLE, authority=source)
    cells = g.input('cells', 'embryo single-cell suspension (~100,000 in 200 µl PBS) or cortex nuclei (~250,000 in 500 µl storage buffer)', MaterialKind.BIOLOGICAL, primary=True, state=mixture('raw'))
    formaldehyde = g.input('formaldehyde', '1.75% formaldehyde 10 min RT, BSA quench; Igepal lysis; 0.5% SDS 65 °C 10 min, Triton quench', MaterialKind.CONDITION)
    nuclei = g.product('nuclei', 'fixed, permeabilized nuclei containing chromatin and RNA', MaterialKind.MOLECULAR, state=mixture('fixed'))
    g.react('fix_permeabilize', 'Fix, lyse and SDS-permeabilize nuclei (procedure C)', StepKind.CROSS_LINKING, ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (nuclei,))
    rt_primer = g.input('gat5_rt', 'published GAT5-RT anchored primer', MaterialKind.OLIGO, state=_oligo('GAT5-RT', GAT5_RT))
    rt = g.input('rt', 'SuperScript IV reverse transcriptase (10 U/µl, 55 °C 15 min)', MaterialKind.CATALYST)
    copied = g.product('copied', 'nuclei containing in-situ cDNA and chromatin', MaterialKind.MOLECULAR, state=mixture('copied'))
    g.react('reverse_transcribe', 'Reverse-transcribe RNA in situ', StepKind.REVERSE_TRANSCRIPTION, ((nuclei, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (copied,))
    mboi = g.input('mboi', 'MboI restriction nuclease (37 °C overnight; 65 °C 20 min inactivation)', MaterialKind.CATALYST)
    digested = g.product('digested', 'MboI-cut chromatin with retained cDNA', MaterialKind.MOLECULAR, state=mixture('digested'))
    g.react('digest', 'Digest chromatin with MboI', StepKind.DIGESTION, ((copied, InputRole.TEMPLATE), (mboi, InputRole.CATALYST)), (digested,))
    ligase = g.input('ligase', 'T4 DNA ligase (Invitrogen, 16 °C 4 h, diluted to ~1 ml)', MaterialKind.CATALYST)
    joined = g.product('joined', 'nuclei with proximity-ligated DNA and cDNA', MaterialKind.MOLECULAR, state=mixture('ligated'))
    g.react('proximity_ligate', 'Ligate nearby genomic ends', StepKind.PROXIMITY_LIGATION, ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (joined,))
    sorter = g.input('facs', 'FACSAria single-nucleus sort into dry 96-well plates (FSC/SSC gating)', MaterialKind.CONDITION)
    sorted_nuclei = g.product('sorted', 'one ligated nucleus per well', MaterialKind.MOLECULAR, state=mixture('ligated'))
    unsorted = g.product('unsorted', 'debris, clumps and unallocated nuclei; disposition unassigned', MaterialKind.MOLECULAR)
    g.react('sort', 'Sort single nuclei into wells (procedure F)', StepKind.PARTITION, ((joined, InputRole.TEMPLATE), (sorter, InputRole.CONDITION)), (sorted_nuclei, unsorted), changes_molecule=False)
    g.mark_terminal(Terminal(unsorted.id,Disposition.UNRESOLVED,'Outside selectedsingle-nucleus lineage; unallocated suspension is not assumed destroyed.'))
    protease = g.input('qiagen_protease', 'Qiagen protease in SC lysis buffer with GAT5 (50 °C 3 h, 70 °C 1 h)', MaterialKind.CATALYST)
    lysed = g.product('lysed', 'deproteinized single-nucleus contact DNA and cDNA', MaterialKind.MOLECULAR, state=mixture('lysed'))
    g.react('lyse', 'Lyse each nucleus with protease (procedure G1-3)', StepKind.DEPROTEINIZATION, ((sorted_nuclei, InputRole.TEMPLATE), (protease, InputRole.CATALYST)), (lysed,))
    gat5 = g.input('gat5', 'published GAT5 primer', MaterialKind.OLIGO, state=_oligo('GAT5', GAT5))
    gat5_7n = g.input('gat5_7n', 'published GAT5-7N primer', MaterialKind.OLIGO, state=_oligo('GAT5-7N', GAT5_7N))
    polymerase = g.input('polymerase', 'Deep Vent (exo-) polymerase, ThermoPol buffer', MaterialKind.CATALYST)
    amplified = g.product('amplified', 'joint single-nucleus amplicons', MaterialKind.MOLECULAR, state=mixture('amplified'))
    g.react('preamp', 'Amplify contact DNA and transcript cDNA in one well (9 stepped random-priming cycles, then 14 cycles)', StepKind.PCR, ((lysed, InputRole.TEMPLATE), (gat5, InputRole.STRUCTURAL), (gat5_7n, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (amplified,), geometry=PCRGeometry.unavailable(lysed.state, reason='Actualmixed contactDNA/ss-cDNA/possibleRNA precursor; native randompriming/firstcycle geometry bounded.'))
    transposome = g.input('transposome', 'TruePrep V2 TTE Mix V50 in TTBL', MaterialKind.ADAPTER, inferred=True, source=source + ' (native adapter oligos unassigned)')
    tn5 = g.input('tn5', 'TruePrep transposase (55 °C 10 min; 0.2% SDS stop)', MaterialKind.CATALYST)
    tagged = g.product('tagged', 'tagmented joint amplicons', MaterialKind.MOLECULAR, inferred=True, source=source + ' (native adapter oligos unassigned)', state=mixture('tagged'))
    tagment_waste = g.product('tagment_waste', 'nonproductive tagmentation classes; recovery unassigned', MaterialKind.MOLECULAR)
    g.react('tagment', 'Tagment 2 µl diluted preamplification product', StepKind.TAGMENTATION, ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (tagged, tagment_waste))
    g.mark_terminal(Terminal(tagment_waste.id, Disposition.UNRESOLVED, 'Nonproductive/end-ineligible molecules are not physically disposed by tagmentation; actual cleanup recovery unassigned.'))
    stop = g.input('sds_stop', '1.25µl0.2%SDS;RT10min', MaterialKind.CONDITION)
    stopped = g.product('stopped', 'SDS-stopped native taggedmixture', MaterialKind.MOLECULAR, state=tagged.state)
    g.react('stop_tagment', 'Stop source tagmentation withSDS', StepKind.QUENCHING, ((tagged, InputRole.TEMPLATE), (stop, InputRole.CONDITION)), (stopped,))
    tagged = stopped
    rna_p7 = g.input('rna_p7', 'published RNA-P7 enrichment primer', MaterialKind.OLIGO, state=_oligo('RNA-P7', RNA_P7))
    i5 = g.input('i5', 'well-specific i5 UD index primer', MaterialKind.OLIGO, inferred=True, source=source + ' (index member bases unavailable)')
    pcr_pol = g.input('pcr_polymerase', 'KAPA HiFi polymerase, GC buffer', MaterialKind.CATALYST)
    extended = g.product('extended', 'actual72°Cinitial-extension mixture', MaterialKind.MOLECULAR, state=mixture('extended'))
    g.react('initial_extension', 'Native72°C3min precycleextension', StepKind.EXTENSION, ((tagged, InputRole.TEMPLATE), (rna_p7, InputRole.STRUCTURAL), (i5, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)), (extended,))
    tagged = extended
    enriched = g.product('enriched', 'i5-extended products with RNA-P7-enriched cDNA reads', MaterialKind.MOLECULAR, inferred=True, source=source + ' (TruePrep adapter bases unavailable)', state=mixture('enriched'))
    g.react('enrich_pcr', '72 °C 3 min, then 5 cycles with RNA-P7 and i5 only (slight cDNA enrichment)', StepKind.PCR, ((tagged, InputRole.TEMPLATE), (rna_p7, InputRole.STRUCTURAL), (i5, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)), (enriched,), geometry=PCRGeometry.unavailable(tagged.state, reason='Actual native extended mixedDNA precursor; selectedTruePrep/i5 landing geometry unassigned.'))
    i7 = g.input('i7', 'well-specific i7 UD index primer', MaterialKind.OLIGO, inferred=True, source=source + ' (index member bases unavailable)')
    indexed = g.product('indexed', 'per-well dual-indexed HiRES library', MaterialKind.LIBRARY, inferred=True, source=source + ' (TruePrep adapter bases unavailable)', state=mixture('indexed'))
    g.react('index_pcr', 'Add i7 to the same well (i5, RNA-P7 and KAPA HiFi carried over) and run 5 more cycles', StepKind.PCR, ((enriched, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL), (i7, InputRole.STRUCTURAL), (rna_p7, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)), (indexed,), geometry=PCRGeometry.unavailable(enriched.state, reason='Actual priorenrichment mixture, not finalduplex; nativei7landing/completeendgeometry unassigned.'))
    others = g.input('other_well_libraries', 'other independently indexed wells; same nativejointclass', MaterialKind.LIBRARY, state=mixture('indexed'))
    pooled = g.product('pooled', 'physicaljoint library pool', MaterialKind.LIBRARY, state=mixture('pooled'))
    g.react('pool_wells', 'Pool indexed wells', StepKind.POOL, ((indexed, InputRole.TEMPLATE), (others, InputRole.TEMPLATE)), (pooled,), changes_molecule=False)
    indexed = pooled
    ampure = g.input('ampure', '0.6x then 0.15x AMPure XP double-sided selection', MaterialKind.SELECTION)
    joint = g.product('joint_library', 'pooled, size-selected joint HiRES library', MaterialKind.LIBRARY, inferred=True, source=source + ' (TruePrep adapter bases unavailable)', state=mixture('pooled'))
    size_waste = g.product('size_waste', 'fragments outside the bead size window and leftover primers', MaterialKind.WASTE)
    g.react('pool_select', 'Pool wells and size-select with AMPure XP', StepKind.CLEANUP, ((indexed, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (joint, size_waste), changes_molecule=False)
    g.discard(size_waste, 'removed by the double-sided AMPure selection')
    g.sequence(joint, SequencingPlan(acquired_roles=('Read 1', 'Read 2'), read_lengths={'Read 1': 150, 'Read 2': 150}, limitation=SEQUENCING_ENDING))
    g.declare_reviewed_scope(('cells', 'formaldehyde', 'fix_permeabilize', 'gat5_rt', 'rt', 'reverse_transcribe', 'mboi', 'digest', 'ligase', 'proximity_ligate', 'facs', 'sort', 'qiagen_protease', 'lyse', 'gat5', 'gat5_7n', 'polymerase', 'preamp', 'transposome', 'tn5', 'tagment', 'sds_stop', 'stop_tagment', 'rna_p7', 'i5', 'pcr_polymerase', 'initial_extension', 'enrich_pcr', 'i7', 'index_pcr', 'other_well_libraries', 'pool_wells', 'ampure', 'pool_select', 'joint_library'), source)
    return g.require_publishable()

CAVEAT += " OriginalRNA/DNA/cDNA classes remain independent; exactUDindexwidth and acquiredindexcycles are unassigned. RNA/readclass separation is computational only. Sourceauthorprotocol defines the selected native chemistry; definingScienceoriginal remains uninspected."
