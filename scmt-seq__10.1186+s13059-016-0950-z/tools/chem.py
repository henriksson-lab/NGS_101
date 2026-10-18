"""Hu et al. 2016: native scMT-seq branches; unassigned inherited strands stay roles."""
from chemdraw import Construct, Scene, Segment
from protocol_model import role
from endprep import EndPreparedDuplex, dA_tailed_scene
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind, InputRole, StepKind,
                            SequencingPlan, Terminal, Disposition)
import seqprimers as sp

TITLE = 'scMT-seq — transcriptome and DNA methylome'
NOTES = '01_scmt-seq.html'
SOURCE = ('Hu et al., Genome Biology (2016), '
          '<a href="https://doi.org/10.1186/s13059-016-0950-z">doi:10.1186/s13059-016-0950-z</a>; '
          'published PMC4858893 Methods205,215,219.')
SUMMARY = ('Physically pipetted nucleus and cytosol enter separate scRRBS and Smart2-seq branches. '
           'SuperScript III is the study-specific RT substitution; methylome PCR has two rounds, '
           '20 then 12 cycles, followed by 180–500 bp PAGE.')
CAVEAT = ('Native Smart2-seq primer/TSO strands, Nextera index primers and scRRBS adapters remain '
          'unassigned in this bounded source scope. The named predecessors are Picelli2014 '
          'nprot.2014.006 and Guo2015 nprot.2015.039; an older generic RRBS recipe and a chosen '
          'Smart-seq2 linker are not substituted. DNA 5′ phosphorylation and native adapter '
          'end geometry remain unassigned. PCR round counts do not establish nestedness or an '
          'inter-round purification.')

RNA = Construct([role('native RNA-library end A'), role('transcript/ERCC-derived cDNA insert'),
                 role('native RNA-library end B')], name='scMT transcriptome endpoint')
DNA = Construct([role('native DNA-library end A'),
                 role('nuclear/lambda DNA insert; protected C versus converted-site T after PCR'),
                 role('native DNA-library end B')], name='scMT methylome endpoint')
RP = DP = ()
FINAL_LIBRARIES = (
    ('Transcriptome library', RNA, RP, 'Native Nextera endpoint roles; exact primer strands unassigned.',
     'The study names Nextera library construction; read/index acquisition is unassigned.'),
    ('DNA methylome library', DNA, DP, 'Native scRRBS endpoint roles after two PCR rounds and PAGE.',
     'HiSeq2500 is reported; read/index cycles and exact primer landing remain unassigned.'),
)
READ_SPECS = {name: tuple(sp.read_spec(read, 'unassigned',
               'Study read/index acquisition unassigned in the inspected source scope', 'undetermined')
               for read in ('Read 1', 'Read 2', 'Index 1 (i7)', 'Index 2 (i5)'))
              for name in ('Transcriptome library', 'DNA methylome library')}


def rna_population(stage):
    sc = Scene()
    for identity in ('cytoplasmic RNA', 'ERCC RNA spike-in'):
        sc.strand(identity, [role(identity)], label=identity + '; independent RNA population')
        if stage != 'raw':
            parts = [role('native RT-primer-derived end: ' + identity),
                     role('first-strand copy: ' + identity)]
            if stage == 'switched':
                parts.append(role('template-switch-derived end: ' + identity))
            sc.strand(identity + ' cDNA', parts, label='single-strand cDNA; RNA remains a separate molecule')
    return sc


def duplex(parts, label):
    return Scene.duplex(parts, label=label,
                        unpaired=tuple(part.name + "'" for part in parts if part.is_role_token()))


def dna_population(stage):
    sc = Scene()
    for identity in ('nuclear genomic DNA', 'lambda DNA spike-in'):
        sc.strand(identity, [role(identity + ': ' + stage)], label=identity + '; independent population')
    return sc


def prepared_dna():
    rows = []
    for identity in ('nuclear genomic DNA', 'lambda DNA spike-in'):
        insert = Construct([role(identity + ' end-repaired body')], name=identity)
        # False suppresses a phosphate drawing; the source has not assigned 5′ P status.
        declared = EndPreparedDuplex(insert, five_prime_phosphate=False,
                                     three_prime_overhang='A', overhang_inferred=False)
        rows.extend(dA_tailed_scene(declared, label=identity + '; 5′ chemistry unassigned').rows())
    return tuple(rows)


def converted_dna():
    sc = Scene()
    for identity in ('nuclear genomic DNA', 'lambda DNA spike-in'):
        for strand in ('strand 1', 'strand 2'):
            sc.strand(identity + strand, [role('native adapter end A: ' + identity),
                      role(identity + ' ' + strand + ': unprotected C→U; protected methyl-C retained'),
                      role('native adapter end B: ' + identity)],
                      label=identity + '; converted physical strand, U before PCR')
    return sc


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cell = g.input('cell', 'one cell', MaterialKind.BIOLOGICAL, primary=True)
    lysis = g.input('lysis', '2% Triton, 20mM NaCl/Tris and RNase inhibitor', MaterialKind.CONDITION)
    ercc = g.input('ercc', 'ERCC RNA spike-in', MaterialKind.MOLECULAR)
    cytosol = g.product('cytosol', 'cytosol and ERCC RNA', MaterialKind.MOLECULAR,
                        state=rna_population('raw'))
    nucleus = g.product('nucleus', 'nucleus pipetted in 0.2µl into 4µl RRBS lysis buffer',
                        MaterialKind.BIOLOGICAL)
    g.react('separate', 'Lyse and physically pipette the nucleus away from cytosol', StepKind.SPLIT,
            ((cell, InputRole.TEMPLATE), (lysis, InputRole.CONDITION), (ercc, InputRole.SPIKE_IN)),
            (cytosol, nucleus))
    rtp = g.input('rt_primer', 'Smart2-seq RT primer; native strand unassigned', MaterialKind.OLIGO)
    rt = g.input('rt', 'SuperScript III replaces SuperScript II', MaterialKind.CATALYST)
    copied = g.product('copied', 'first-strand RNA/cDNA mixture', MaterialKind.MOLECULAR,
                       state=rna_population('copied'))
    g.react('reverse_transcribe', 'Reverse-transcribe cytosolic poly(A) RNA and RNA spike-ins',
            StepKind.REVERSE_TRANSCRIPTION, ((cytosol, InputRole.TEMPLATE),
            (rtp, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (copied,))
    tso = g.input('tso', 'Smart2-seq TSO; native bases/modifications unassigned', MaterialKind.OLIGO)
    switched = g.product('switched', 'template-switched single-strand cDNA reaction',
                         MaterialKind.MOLECULAR, state=rna_population('switched'))
    g.react('template_switch', 'Template-switch nascent cDNA; inherited component sequence unassigned',
            StepKind.TEMPLATE_SWITCHING, ((copied, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
            (rt, InputRole.CATALYST)), (switched,))
    pre_f = g.input('pre_f', 'cDNA PCR forward-binding primer role; identity unassigned', MaterialKind.OLIGO)
    pre_r = g.input('pre_r', 'cDNA PCR reverse-binding primer role; identity relation unassigned', MaterialKind.OLIGO)
    pre_pol = g.input('pre_pol', 'Smart2-seq cDNA PCR polymerase', MaterialKind.CATALYST)
    amplicon = duplex([role('native preamplification end A'), role('transcript/ERCC cDNA copy'),
                      role('native preamplification end B')], 'cDNA PCR product')
    amplified = g.product('amplified', 'amplified cDNA', MaterialKind.MOLECULAR, state=amplicon)
    residual = g.product('rna_residual', 'RNA templates/nonamplified products; physical fate unassigned',
                         MaterialKind.MOLECULAR)
    g.react('preamp', 'Amplify template-switched cDNA', StepKind.PCR,
            ((switched, InputRole.TEMPLATE), (pre_f, InputRole.STRUCTURAL),
             (pre_r, InputRole.STRUCTURAL), (pre_pol, InputRole.CATALYST)), (amplified, residual),
            geometry=PCRGeometry.unavailable(switched.state,
                     reason='Actual single-strand switched-cDNA mixture; native primer strands unassigned.'))
    g.mark_terminal(Terminal(residual.id, Disposition.UNRESOLVED,
                    reason='Nonamplification does not establish RNA destruction or physical removal.'))
    recovery = g.input('cdna_recovery', 'cDNA purification conditions; implementation unassigned', MaterialKind.CONDITION)
    purified = g.product('purified_cdna', 'purified amplified cDNA', MaterialKind.MOLECULAR, state=amplicon)
    waste = g.product('cdna_waste', 'cDNA purification non-retained fraction; composition unassigned', MaterialKind.WASTE)
    g.react('purify_cdna', 'Purify amplified cDNA', StepKind.CLEANUP,
            ((amplified, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)), (purified, waste))
    g.mark_terminal(Terminal(waste.id, Disposition.DISCARD, reason='Physical purification non-retained fraction.'))
    amount = g.input('cdna_amount', '0.1ng purified cDNA', MaterialKind.CONDITION)
    aliquot = g.product('cdna_aliquot', '0.1ng cDNA tagmentation input', MaterialKind.MOLECULAR, state=amplicon)
    rest = g.product('remaining_cdna', 'amplified cDNA outside selected aliquot', MaterialKind.MOLECULAR)
    g.react('take_cdna', 'Take 0.1ng amplified/purified cDNA', StepKind.SAMPLING,
            ((purified, InputRole.TEMPLATE), (amount, InputRole.CONDITION)), (aliquot, rest))
    g.mark_terminal(Terminal(rest.id, Disposition.UNRESOLVED, reason='Remaining aliquot is outside selected library route.'))
    transposome = g.input('transposome', 'Nextera loaded adapters; native strands unassigned', MaterialKind.ADAPTER)
    tn5 = g.input('tn5', 'Nextera Tn5', MaterialKind.CATALYST)
    tagged_state = duplex([role('transferred native end A'), role('tagmented cDNA body'),
                          role('transferred native end B; gap topology unassigned')], 'actual tagmented cDNA precursor')
    tagged = g.product('tagged', 'Nextera-tagged cDNA before library PCR', MaterialKind.MOLECULAR,
                       state=tagged_state)
    g.react('tagment', 'Tagment selected cDNA with Nextera', StepKind.TAGMENTATION,
            ((aliquot, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST)), (tagged,))
    rna_f = g.input('rna_f', 'native RNA-library PCR primer A', MaterialKind.OLIGO)
    rna_r = g.input('rna_r', 'native RNA-library PCR primer B', MaterialKind.OLIGO)
    rna_pol = g.input('rna_pol', 'Nextera library PCR polymerase', MaterialKind.CATALYST)
    rna_library = g.product('rna_library', RNA.name, MaterialKind.LIBRARY, state=RNA)
    g.react('rna_pcr', 'Complete native Nextera library; cycle count unassigned', StepKind.PCR,
            ((tagged, InputRole.TEMPLATE), (rna_f, InputRole.STRUCTURAL),
             (rna_r, InputRole.STRUCTURAL), (rna_pol, InputRole.CATALYST)), (rna_library,),
            geometry=PCRGeometry.unavailable(tagged.state, reason='Actual tagged cDNA input; native gap/PCR geometry unassigned.'))

    rrbs_lysis = g.input('rrbs_lysis', 'RRBS nuclear DNA-release conditions', MaterialKind.CONDITION)
    genome = g.product('genome', 'released nuclear genomic DNA', MaterialKind.MOLECULAR,
                       state=duplex([role('released nuclear genomic DNA')], 'genomic DNA before lambda addition'))
    g.react('release_dna', 'Release DNA from separated nucleus', StepKind.EXTRACTION,
            ((nucleus, InputRole.TEMPLATE), (rrbs_lysis, InputRole.CONDITION)), (genome,))
    lambda_dna = g.input('lambda', 'lambda DNA spike-in', MaterialKind.MOLECULAR)
    mspi = g.input('mspi', 'MspI', MaterialKind.CATALYST)
    fragments = g.product('fragments', 'MspI-digested genomic and lambda DNA', MaterialKind.MOLECULAR,
                          state=dna_population('MspI-digested DNA'))
    g.react('digest', 'MspI digest with independent lambda DNA spike-in', StepKind.DIGESTION,
            ((genome, InputRole.TEMPLATE), (lambda_dna, InputRole.SPIKE_IN),
             (mspi, InputRole.CATALYST)), (fragments,))
    prep = g.input('endprep', 'source-named end-repair/dA-tailing reagents; enzyme identities unassigned', MaterialKind.CATALYST)
    prepared = g.product('prepared', 'end-repaired and 3′-dA-tailed DNA', MaterialKind.MOLECULAR,
                         state=prepared_dna())
    g.react('prepare_ends', 'End repair and dA tailing', StepKind.END_REPAIR,
            ((fragments, InputRole.TEMPLATE), (prep, InputRole.CATALYST)), (prepared,))
    adapters = g.input('adapters', 'native scRRBS adaptors; strands/index/terminal geometry unassigned', MaterialKind.ADAPTER)
    ligase = g.input('ligase', 'adapter ligation enzyme', MaterialKind.CATALYST)
    ligated = g.product('ligated', 'adapter-ligated genomic/lambda fragments', MaterialKind.MOLECULAR,
                        state=dna_population('native adapter end A–DNA insert–native adapter end B'))
    g.react('ligate', 'Ligate native scRRBS adaptors', StepKind.LIGATION,
            ((prepared, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (ligated,),
            geometry=AdapterLigationGeometry.unavailable('Native adapter strands and ligation-facing end geometry unassigned; no canonical Y/T-overhang substituted.'))
    bs = g.input('bisulfite', 'bisulfite conversion conditions', MaterialKind.CONDITION)
    converted = g.product('converted', 'physical bisulfite-converted single strands', MaterialKind.MOLECULAR,
                          state=converted_dna())
    g.react('convert', 'Bisulfite convert adapted DNA; unprotected C→U before PCR', StepKind.CHEMICAL_CONVERSION,
            ((ligated, InputRole.TEMPLATE), (bs, InputRole.CONDITION)), (converted,))
    bs_recovery = g.input('bs_recovery', 'converted-DNA purification conditions', MaterialKind.CONDITION)
    recovered = g.product('recovered_bs', 'purified converted DNA', MaterialKind.MOLECULAR, state=converted.state)
    waste = g.product('bs_waste', 'conversion purification non-retained fraction', MaterialKind.WASTE)
    g.react('purify_bs', 'Purify bisulfite-converted DNA', StepKind.CLEANUP,
            ((converted, InputRole.TEMPLATE), (bs_recovery, InputRole.CONDITION)), (recovered, waste))
    g.mark_terminal(Terminal(waste.id, Disposition.DISCARD, reason='Source-described purification non-retained fraction.'))
    dna_f = g.input('dna_f', 'first-round native scRRBS PCR primer A', MaterialKind.OLIGO)
    dna_r = g.input('dna_r', 'first-round native scRRBS PCR primer B', MaterialKind.OLIGO)
    dna_pol = g.input('dna_pol', 'scRRBS PCR polymerase; identity unassigned', MaterialKind.CATALYST)
    round1_state = duplex([role('native DNA PCR end A'), role('PCR-copied DNA: converted U positions yield T/A pairs'),
                           role('native DNA PCR end B')], 'first-round PCR DNA')
    round1 = g.product('round1', 'first-round DNA amplicons; second-round handoff unassigned',
                       MaterialKind.MOLECULAR, state=round1_state)
    g.react('dna_pcr1', 'First methylome PCR:20 cycles', StepKind.PCR,
            ((recovered, InputRole.TEMPLATE), (dna_f, InputRole.STRUCTURAL),
             (dna_r, InputRole.STRUCTURAL), (dna_pol, InputRole.CATALYST)), (round1,),
            geometry=PCRGeometry.unavailable(recovered.state, reason='Actual converted single-strand input; native primers unassigned.'))
    second_f = g.input('second_f', 'second-round PCR primer A; relation to round1 primer unassigned', MaterialKind.OLIGO)
    second_r = g.input('second_r', 'second-round PCR primer B; relation to round1 primer unassigned', MaterialKind.OLIGO)
    round2 = g.product('round2', 'second-round methylome amplicons', MaterialKind.MOLECULAR, state=DNA)
    g.react('dna_pcr2', 'Second methylome PCR:12 cycles; inter-round handling unassigned', StepKind.PCR,
            ((round1, InputRole.TEMPLATE), (second_f, InputRole.STRUCTURAL),
             (second_r, InputRole.STRUCTURAL), (dna_pol, InputRole.CATALYST)), (round2,),
            geometry=PCRGeometry.unavailable(round1.state, reason='Actual first-round product; no nestedness, new primer identity or inter-round cleanup inferred.'))
    page = g.input('page', 'PAGE selection of180–500bp', MaterialKind.SELECTION)
    dna_library = g.product('dna_library', DNA.name, MaterialKind.LIBRARY, state=DNA)
    rejected = g.product('page_rejected', 'DNA outside selected PAGE fraction', MaterialKind.WASTE)
    g.react('page_select', 'Select180–500bp after both PCR rounds', StepKind.SELECTION,
            ((round2, InputRole.TEMPLATE), (page, InputRole.SELECTION)), (dna_library, rejected))
    g.mark_terminal(Terminal(rejected.id, Disposition.DISCARD, reason='PAGE fraction outside180–500bp.'))
    g.sequence(rna_library, SequencingPlan(library=RNA, primers=(), limitation='Native RNA sequencing/index primers and acquired cycle allocation unassigned.'))
    g.sequence(dna_library, SequencingPlan(library=DNA, primers=(), limitation='HiSeq2500 known; native DNA sequencing/index primers and cycle allocation unassigned.'))
    g.declare_reviewed_scope((
        'cell','lysis','ercc','separate','cytosol','nucleus','rt_primer','rt','reverse_transcribe',
        'copied','tso','template_switch','switched','pre_f','pre_r','pre_pol','preamp','amplified',
        'cdna_recovery','purify_cdna','purified_cdna','cdna_amount','take_cdna','cdna_aliquot',
        'transposome','tn5','tagment','tagged','rna_f','rna_r','rna_pol','rna_pcr','rna_library',
        'rrbs_lysis','release_dna','genome','lambda','mspi','digest','fragments','endprep',
        'prepare_ends','prepared','adapters','ligase','ligate','ligated','bisulfite','convert',
        'converted','bs_recovery','purify_bs','recovered_bs','dna_f','dna_r','dna_pol',
        'dna_pcr1','round1','second_f','second_r','dna_pcr2','round2','page','page_select','dna_library'), SOURCE)
    return g.require_publishable()
