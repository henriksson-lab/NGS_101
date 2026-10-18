"""Source-owned fragmented, random-nonamer SHAPE-MaP preparation."""
from chemdraw import Scene, Segment
from protocol_model import role
from protocol_graph import (ProtocolGraph, MaterialKind as M, InputRole as I,
                            StepKind as K, SequencingPlan)
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

TITLE = 'SHAPE-MaP — RNA structure by mutational profiling'
NOTES = '01_shape-map.html'
SOURCE = '<a href="https://doi.org/10.1038/nmeth.3029">Siegfried et al. (2014)</a>, published author manuscript PMC4259394, Online Methods pp. 8–10.'
SUMMARY = 'SHAPE RNA 2′-O-adducts cause mutations during MaP reverse transcription. This graph selects the fragmented ribosome/group II intron/HCV IRES branch with random nonamer DNA primers.'
CAVEAT = 'Native component states show independent RNA, DNA and residual reagent populations. Ordered TruSeq forked-adapter strands, library PCR primer sites, index identities and acquisition cycles remain unassigned in this bounded review. No canonical Illumina endpoint supplies those values. Untreated/denatured controls, HIV LNA priming and targeted PCR branches are outside this selected preparation.'
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SHOW_READ_CONTENT = False
SEQUENCING_ENDING = 'The source names MiSeq or HiSeq. Exact run lengths and primer sites remain unassigned. Demultiplexing, quality trimming (retain reads ≥25 nt after Q10 trimming), alignment and mutation counting are computational operations; they do not physically purify the library.'


def native(stage):
    sc = Scene()
    if stage in ('rna', 'modified', 'fragmented'):
        sc.strand('RNA', [role('original RNA; physical ribonucleotides; '+stage)],
                  label='RNA; 2′-O-adducts are substituents, not added bases' if stage != 'rna' else 'folded input RNA')
    else:
        sc.strand('original_RNA', [role('possible carried original SHAPE-modified RNA/fragments; recovery unresolved')], label='original RNA lineage; adducts not copied into DNA')
        if stage == 'rt':
            sc.strand('cDNA', [Segment('incorporated random DNA primer', 'N'*9), role('new DNA copy; mutations at some adduct sites; exact positions unassigned')], label='independent primer-extended DNA; not modified RNA')
        else:
            body = [role('MaP DNA insert; sequence and termini unassigned')]
            if stage in ('adapted', 'library'):
                body = [role('forked-adapter-derived DNA left end; assembly unassigned'), *body,
                        role('forked-adapter-derived DNA right end; assembly unassigned')]
            sc.strand('DNA_one', body, label='DNA '+stage+' strand population; exact ends unassigned')
            sc.strand('DNA_two', [role('independent second DNA strand; pairing/overhangs unassigned')], label='second-strand DNA; no fabricated duplex columns')
        sc.strand('residual_oligos', [role('possible free nonamers, off-target DNA and native reagents; class recovery unresolved')])
        if stage in ('adapted', 'library'):
            sc.strand('residual_adapters', [role('possible free forked adapters, adapter dimers and unligated DNA; physical fates unresolved')])
        if stage == 'library':
            sc.strand('PCR_daughters', [role('new DNA amplification products; selected primer/end assembly unassigned')], label='PCR daughters; original RNA adducts are not copied chemical groups')
            sc.strand('residual_primers', [role('possible original PCR primers and original ligated templates; no post-PCR purification inferred')])
    return sc


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    def inp(i, label, kind=M.CONDITION, state=None):
        return g.input(i, label, kind, state=state)
    def op(i, label, kind, x, extra, stage):
        y = g.product(i+'_out', label+' mixture', M.MOLECULAR, state=native(stage))
        g.react(i, label, kind, ((x,I.TEMPLATE), *extra), (y,))
        return y
    def cleanup(i, label, x):
        selection = inp(i+'_selection', label, M.SELECTION)
        y = g.product(i+'_out', label+' retained mixture', x.kind, state=x.state)
        waste = g.product(i+'_waste', label+' nonretained fraction', M.WASTE)
        g.react(i, label, K.CLEANUP, ((x,I.TEMPLATE),(selection,I.SELECTION)), (y,waste))
        g.discard(waste, 'Physical cleanup fraction; exact RNA/DNA/oligo recoveries are unresolved, not assumed class-pure.')
        return y
    raw = g.input('rna', 'folded RNA aliquot; selected model-RNA branch', M.MOLECULAR, state=native('rna'), primary=True)
    reagent = inp('shape', 'SHAPE reagent condition; untreated/denatured controls independently prepared')
    x = op('modify', 'Acylate flexible RNA 2′ hydroxyls', K.CHEMICAL_CONVERSION, raw, ((reagent,I.CONDITION),), 'modified')
    x = cleanup('modification_cleanup', 'RNeasy MinElute or G-50 purification named by source; selection unresolved', x)
    frag = inp('fragment_condition', 'Mg-containing fragmentation condition; source Online Methods p. 9')
    x = op('fragment', 'Fragment original modified RNA', K.FRAGMENTATION, x, ((frag,I.CONDITION),), 'fragmented')
    x = cleanup('fragment_desalt', 'G-50 spin-column desalting', x)
    primer_scene = Scene(); primer_scene.strand('primer', [Segment('random nonamer DNA','N'*9)], label='source random nonamer DNA pool')
    primer = inp('rt_primer', 'random nonamer DNA primers; selected non-HIV fragmented branch', M.OLIGO, primer_scene)
    rt = inp('map_rt', 'SuperScript II under source Mn-containing MaP conditions', M.CATALYST)
    x = op('reverse_transcribe', 'Copy RNA into mutationally encoded DNA', K.REVERSE_TRANSCRIPTION, x, ((primer,I.STRUCTURAL),(rt,I.CATALYST)), 'rt')
    x = cleanup('rt_desalt', 'G-50 spin-column desalting after RT', x)
    second = inp('second_mix', 'NEBNext E6111 second-strand synthesis module', M.CATALYST)
    x = op('second_strand', 'Synthesize independent second DNA strand', K.EXTENSION, x, ((second,I.CATALYST),), 'duplex')
    x = cleanup('second_cleanup', 'PureLink Micro PCR cleanup', x)
    repair = inp('repair_mix', 'NEBNext E6050 end repair', M.CATALYST)
    x = op('end_repair', 'Repair double-stranded cDNA ends', K.END_REPAIR, x, ((repair,I.CATALYST),), 'repaired')
    x = cleanup('repair_cleanup', 'AMPure XP 1.6:1 bead cleanup', x)
    tail = inp('tail_mix', 'NEBNext E6053 dA tailing', M.CATALYST)
    x = op('da_tail', 'Add DNA terminal dA; strand-specific geometry unassigned', K.TAILING, x, ((tail,I.CATALYST),), 'tailed')
    ads = Scene(); ads.strand('adapter_a',[role('ordered TruSeq forked-adapter strand A; bases/ends unassigned')]); ads.strand('adapter_b',[role('independent ordered forked-adapter strand B; bases/ends unassigned')])
    adapter = inp('truseq_adapters', 'Illumina-compatible forked TruSeq adapters named by source', M.ADAPTER, ads)
    ligase = inp('ligation_mix', 'NEB M2200 quick ligation module', M.CATALYST)
    ligated = g.product('adapted', 'native adapter-ligated DNA mixture', M.MOLECULAR, state=native('adapted'))
    g.react('ligate_adapters', 'Ligate forked adapters onto dA-tailed DNA', K.LIGATION, ((x,I.TEMPLATE),(adapter,I.STRUCTURAL),(ligase,I.CATALYST)), (ligated,), geometry=AdapterLigationGeometry.unavailable('Ordered forked-adapter end bases, polarity and exact junction assembly remain unassigned in this source scope.'))
    f = inp('pcr_forward', 'fragmented-library forward PCR primer; identity unassigned', M.OLIGO)
    r = inp('pcr_reverse', 'fragmented-library reverse PCR primer; identity unassigned', M.OLIGO)
    pol = inp('pcr_polymerase', 'Q5 hot-start high-fidelity polymerase; source 30-cycle emulsion PCR', M.CATALYST)
    library = g.product('library', 'PCR DNA daughters plus possible native residuals', M.LIBRARY, state=native('library'))
    g.react('library_pcr', 'Emulsion PCR of actual ligated DNA mixture', K.PCR, ((ligated,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)), (library,), geometry=PCRGeometry.unavailable(ligated.state, reason='Actual incoming ligated material is the template; selected primer identities/sites and full amplified assembly remain unassigned.'))
    g.sequence(library, SequencingPlan(limitation=SEQUENCING_ENDING+' Source reports Qubit quantification, Bioanalyzer verification and pooling before sequencing; sampling volumes, physical QC aliquots and companion libraries are not reconstructed.'))
    g.declare_reviewed_scope(REVIEWED, SOURCE)
    return g.require_publishable()

REVIEWED = ('rna','shape','modify','modification_cleanup_selection','modification_cleanup','fragment_condition','fragment','fragment_desalt_selection','fragment_desalt','rt_primer','map_rt','reverse_transcribe','rt_desalt_selection','rt_desalt','second_mix','second_strand','second_cleanup_selection','second_cleanup','repair_mix','end_repair','repair_cleanup_selection','repair_cleanup','tail_mix','da_tail','truseq_adapters','ligation_mix','ligate_adapters','pcr_forward','pcr_reverse','pcr_polymerase','library_pcr','library')
