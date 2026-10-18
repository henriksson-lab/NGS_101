"""Clark et al. 2018: selected EL16 scNMT-seq branch with native molecular roles."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'lib')]
from chemdraw import Construct, Scene, feature
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind, InputRole, StepKind,
                            SequencingPlan, Terminal, Disposition)
import seqprimers as sp

TITLE = 'scNMT-seq — RNA, DNA methylation and accessibility'
NOTES = '01_scnmt-seq.html'
SOURCE = ('Clark et al., Nature Communications (2018), '
          '<a href="https://doi.org/10.1038/s41467-018-03149-4">doi:10.1038/s41467-018-03149-4</a>; '
          'published PMC5823944 Methods136–144 and interpretation Results.')
SUMMARY = ('M.CviPI marks accessible GpC before RNA capture; separate Smart-seq2/Nextera XT '
           'and native 6NF-based bisulfite-DNA branches retain matching plate-well identity. '
           'This graph selects the EL16 paired125-base cohort.')
CAVEAT = ('The source names 6NF for both DNA synthesis stages. Its native sequence, handle and '
          'tail width remain unassigned; the earlier borrowed N9 scBS-seq scaffold is removed. '
          'Smart-seq2 TSO identity also remains unassigned: Picelli2013 SupplementaryTable2 '
          '(rGrG+G row) prints AC, while Picelli2014 ReagentSetup prints ACAT. These are '
          'different published options; citing both does not choose '
          'the study oligo. No chosen linker is inherited here. Native index members, cycle '
          'allocation and PCR geometry remain unassigned. GCG ambiguity and CCG off-target '
          'marking are computational interpretation exclusions, not DNA purification fractions.')

RNA_INDEX = feature('scnmt_rna_index', 'sample_index', 'unknown', group='plate_well',
                    note='Sample/plate-well identity; native index member, count and width unassigned.')
DNA_INDEX = feature('scnmt_dna_itag', 'sample_index', 'unknown', group='plate_well',
                    note='iTAG-associated sample/plate-well identity; member and width unassigned.')
RNA_LIBRARY = Construct([role('native RNA end A'), role('transcript-derived DNA insert'),
                         role('native indexed RNA end region; geometry/count unassigned',
                              molecular_feature=RNA_INDEX)], name='scNMT EL16 RNA endpoint')
DNA_LIBRARY = Construct([role('native DNA end A; PE1.0-related strand unassigned'),
                         role('DNA insert; PCR-copied methylation/accessibility signal'),
                         role('native iTAG-associated end; geometry unassigned',
                              molecular_feature=DNA_INDEX)], name='scNMT EL16 DNA endpoint')
RNA_PRIMERS = DNA_PRIMERS = ()
FINAL_LIBRARIES = (
    ('RNA branch', RNA_LIBRARY, RNA_PRIMERS,
     'Native Nextera XT endpoint roles; selected EL16 RNA paired125.',
     'Native sequencing-primer landing and index acquisition unassigned.'),
    ('DNA branch', DNA_LIBRARY, DNA_PRIMERS,
     'Native 6NF/PE1.0/iTAG endpoint roles; selected EL16 DNA paired125.',
     'No historical scBS-seq adapter shell or index width is substituted.'),
)
READ_SPECS = {
    name: (sp.read_spec('Read 1', 125, content, 'include'),
           sp.read_spec('Read 2', 125, content, 'include'),
           sp.read_spec('Index 1 (i7)', 'unassigned', 'Native index acquisition unassigned', 'undetermined'),
           sp.read_spec('Index 2 (i5)', 'unassigned', 'Native index acquisition unassigned', 'undetermined'))
    for name, content in (('RNA branch', 'EL16 transcript-derived sequence'),
                           ('DNA branch', 'EL16 bisulfite-DNA sequence; context-filtered computationally'))
}
READ_SPECS_NOTE = {
    'RNA branch': 'Selected EL16: HiSeq2500 PE125. Separate E14 EB RNA cohort uses75-base single-end reads; it has no acquired Read2. Index acquisition unassigned.',
    'DNA branch': 'Selected EL16: HiSeq2500 PE125. Separate initial20-library QC cohort including3controls uses MiSeq50-base single-end; E14 EB DNA uses NextSeq PE75. These are different runs; index acquisition unassigned.',
}


def duplex(parts, label):
    return Scene.duplex(parts, label=label,
                        unpaired=tuple(p.name + "'" for p in parts if p.is_role_token()))


def cellular_state(marked=False):
    sc = Scene()
    sc.strand('RNA', [role('poly(A) RNA retained independently')], label='RNA population')
    description = ('DNA: accessible GpC can carry M.CviPI 5mC; endogenous CpG methylation retained'
                   if marked else 'chromatin DNA with endogenous CpG methylation')
    sc.strand('DNA', [role(description)], label='DNA population; no covalent RNA linkage')
    return sc


def rna_state(stage):
    sc = Scene()
    sc.strand('RNA', [role('captured poly(A) RNA')], label='RNA template')
    if stage == 'captured':
        sc.strand('DNA carryover', [role('marked DNA carried with bead fraction before recovery washes')], label='independent DNA carryover; three washes transfer it to DNA fraction')
    if stage in ('copied', 'switched'):
        parts = [role('native oligo-dT-derived end'), role('first-strand transcript copy')]
        if stage == 'switched':
            parts.append(role('native TSO-derived end; exact linker/modifications unassigned'))
        sc.strand('cDNA', parts, label='single-strand cDNA; RNA is a separate molecule')
    return sc


def converted_state():
    sc = Scene()
    for n in (1, 2):
        sc.strand('DNA' + str(n), [role('converted original strand' + str(n) +
                   ': unprotected C→U; methyl-C retained')], label='original converted DNA; U before copying')
    return sc


def post_exoi_state():
    sc = Scene()
    sc.strand('DNA population', [role('post-ExoI DNA; surviving original/copy assignments and strand topology unassigned')], label='retained DNA population; not a completed duplex library')
    return sc


def first_copy_state():
    sc = converted_state()
    sc.strand('copies', [role('native 6NF-derived first-copy end; handle/tail width unassigned'),
                        role('new DNA copies of converted originals; copied U can yield A')],
              label='independent first-stage DNA copies; not a completed library')
    return sc


def second_copy_state():
    sc = Scene()
    sc.strand('templates', [role('recovered first-stage DNA templates after ExoI/AMPure; original/copy identity unassigned')],
              label='DNA template population')
    sc.strand('second copies', [role('native second-stage 6NF-derived end; same named reagent'),
                               role('new complementary DNA copies'),
                               role('first-stage 6NF-related terminal boundary; topology unassigned')],
              label='actual copied PCR precursor; no completed PE1.0/iTAG shell')
    return sc


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cell = g.input('cell', 'single cell assigned to a plate well', MaterialKind.BIOLOGICAL,
                   state=cellular_state(), primary=True)
    methylase = g.input('mcvi_pi', 'M.CviPI', MaterialKind.CATALYST)
    sam = g.input('sam', 'SAM-containing2.5µl marking reaction', MaterialKind.CONDITION)
    marked = g.product('marked', 'GpC-marked cell with RNA retained', MaterialKind.BIOLOGICAL,
                       state=cellular_state(True))
    g.react('mark_accessibility', 'M.CviPI marking:37°C15min', StepKind.CHEMICAL_CONVERSION,
            ((cell, InputRole.TEMPLATE), (methylase, InputRole.CATALYST),
             (sam, InputRole.CONDITION)), (marked,))
    stop = g.input('rlt_plus', '5µl RLTplus', MaterialKind.CONDITION)
    lysate = g.product('lysate', 'stopped lysate containing independent RNA and DNA',
                       MaterialKind.MOLECULAR, state=marked.state)
    g.react('stop_marking', 'Stop marking and lyse with RLTplus', StepKind.QUENCHING,
            ((marked, InputRole.TEMPLATE), (stop, InputRole.CONDITION)), (lysate,))
    beads = g.input('oligodt_beads', 'Smart-seq2 oligo-dT pre-annealed to magnetic beads', MaterialKind.SELECTION)
    rna = g.product('rna_bound', 'bead-bound poly(A) RNA', MaterialKind.MOLECULAR, state=rna_state('captured'))
    dna_initial = g.product('dna_lysate', 'transferred DNA-containing lysate', MaterialKind.MOLECULAR,
                            state=duplex([role('marked genomic DNA')], 'DNA fraction'))
    g.react('capture_rna', 'Capture RNA on beads; transfer DNA-containing lysate', StepKind.AFFINITY_CAPTURE,
            ((lysate, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (rna, dna_initial))
    wash_buffer = g.input('wash_buffer', 'three15µl bead washes', MaterialKind.CONDITION)
    washed_rna = g.product('washed_rna', 'washed bead-captured RNA', MaterialKind.MOLECULAR, state=rna_state('washed'))
    washes = g.product('dna_washes', 'three transferred DNA-recovery washes', MaterialKind.MOLECULAR,
                        state=dna_initial.state)
    g.react('wash_capture', 'Wash beads three times; transfer all washes to DNA fraction', StepKind.PARTITION,
            ((rna, InputRole.TEMPLATE), (wash_buffer, InputRole.CONDITION)), (washed_rna, washes))
    dna = g.product('dna_fraction', 'combined DNA lysate plus three washes', MaterialKind.MOLECULAR,
                    state=dna_initial.state)
    g.react('combine_dna', 'Combine transferred DNA lysate and washes', StepKind.POOL,
            ((dna_initial, InputRole.TEMPLATE), (washes, InputRole.TEMPLATE)), (dna,))

    rt = g.input('rt', '100U SuperScript II', MaterialKind.CATALYST)
    rt_primer = g.input('rt_primer', 'bead-associated Smart-seq2 oligo-dT; native strand unassigned', MaterialKind.OLIGO)
    first = g.product('rna_first', 'first-strand RNA/cDNA reaction', MaterialKind.MOLECULAR, state=rna_state('copied'))
    g.react('reverse_transcribe', 'RT incubation:42°C60min,50°C30min,60°C10min', StepKind.REVERSE_TRANSCRIPTION,
            ((washed_rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,))
    tso = g.input('tso', 'named Smart-seq2 TSO; ordered variant unassigned', MaterialKind.OLIGO)
    switched = g.product('rna_switched', 'switched single-strand cDNA reaction', MaterialKind.MOLECULAR,
                         state=rna_state('switched'))
    g.react('template_switch', 'Template switch during RT; native linker unassigned', StepKind.TEMPLATE_SWITCHING,
            ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (switched,))
    ispcr = g.input('ispcr', 'ISPCR primer; native strand unassigned', MaterialKind.OLIGO)
    kapa = g.input('kapa', 'KAPA PCR reagents', MaterialKind.CATALYST)
    cdna_state = duplex([role('native ISPCR end A'), role('transcript-derived cDNA'),
                        role('native ISPCR end B')], 'cDNA PCR product')
    amplified = g.product('rna_amplified', '18-cycle cDNA amplicons', MaterialKind.MOLECULAR, state=cdna_state)
    rna_residual = g.product('rna_residual', 'RNA and nonamplified reaction products; partition unassigned', MaterialKind.MOLECULAR)
    g.react('preamp_rna', 'Amplify switched cDNA:18cycles with ISPCR/KAPA', StepKind.PCR,
            ((switched, InputRole.TEMPLATE), (ispcr, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)),
            (amplified, rna_residual), geometry=PCRGeometry.unavailable(switched.state,
               reason='Actual switched single-strand cDNA mixture; native ISPCR binding geometry unassigned.', symmetric_species=True))
    g.mark_terminal(Terminal(rna_residual.id, Disposition.UNRESOLVED, reason='Nonamplification does not establish physical RNA removal.'))
    ampure = g.input('ampure', 'AMPure beads', MaterialKind.SELECTION)
    cdna_clean = g.product('cdna_clean', 'cDNA recovered in20µl water', MaterialKind.MOLECULAR, state=cdna_state)
    cdna_waste = g.product('cdna_waste', '1:1AMPure non-retained fraction', MaterialKind.WASTE)
    g.react('clean_cdna', 'Purify cDNA1:1AMPure; elute20µl water', StepKind.CLEANUP,
            ((amplified, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (cdna_clean, cdna_waste))
    g.mark_terminal(Terminal(cdna_waste.id, Disposition.DISCARD, reason='Physical cDNA cleanup fraction; species composition unassigned.'))
    amount = g.input('cdna_amount', '100–400pg recovered cDNA', MaterialKind.CONDITION)
    aliquot = g.product('cdna_aliquot', 'selected100–400pg cDNA', MaterialKind.MOLECULAR, state=cdna_state)
    rest = g.product('cdna_rest', 'cDNA outside selected input', MaterialKind.MOLECULAR)
    g.react('sample_cdna', 'Select100–400pg cDNA for one-fifth-volume Nextera XT', StepKind.SAMPLING,
            ((cdna_clean, InputRole.TEMPLATE), (amount, InputRole.CONDITION)), (aliquot, rest))
    g.mark_terminal(Terminal(rest.id, Disposition.UNRESOLVED, reason='Remaining aliquot handling is outside selected route.'))
    loaded = g.input('transposome', 'Nextera XT loaded adapters; native strands unassigned', MaterialKind.ADAPTER)
    tn5 = g.input('tn5', 'Nextera XT Tn5', MaterialKind.CATALYST)
    tagged_state = duplex([role('transferred native end A'), role('tagmented transcript-derived DNA'),
                           role('transferred native end B; gap topology unassigned')], 'tagmented precursor')
    tagged = g.product('rna_tagged', 'native tagged cDNA before indexing PCR', MaterialKind.MOLECULAR, state=tagged_state)
    g.react('tagment_rna', 'Tagment selected cDNA in one-fifth-volume Nextera XT prep', StepKind.TAGMENTATION,
            ((aliquot, InputRole.TEMPLATE), (loaded, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (tagged,))
    index_f = g.input('rna_index_f', 'native RNA indexed PCR primer A; member unassigned', MaterialKind.OLIGO)
    index_r = g.input('rna_index_r', 'native RNA indexed PCR primer B; member unassigned', MaterialKind.OLIGO)
    nextera_pcr = g.input('nextera_pcr', 'native Nextera XT library PCR reagents; enzyme identity unassigned', MaterialKind.CATALYST)
    rna_library = g.product('rna_library', RNA_LIBRARY.name, MaterialKind.LIBRARY, state=RNA_LIBRARY)
    g.react('rna_index_pcr', 'Nextera XT library indexing; cycle count unassigned', StepKind.PCR,
            ((tagged, InputRole.TEMPLATE), (index_f, InputRole.STRUCTURAL),
             (index_r, InputRole.STRUCTURAL), (nextera_pcr, InputRole.CATALYST)), (rna_library,),
            geometry=PCRGeometry.unavailable(tagged.state, reason='Actual tagged cDNA; native adapter/gap/index primer geometry unassigned.'))

    dna_clean = g.product('dna_clean', 'DNA recovered in10µl water', MaterialKind.MOLECULAR, state=dna.state)
    dna_waste = g.product('dna_waste', '0.8×AMPure non-retained DNA fraction', MaterialKind.WASTE)
    g.react('clean_dna', 'Purify combined DNA0.8×AMPure; elute10µl water', StepKind.CLEANUP,
            ((dna, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (dna_clean, dna_waste))
    g.mark_terminal(Terminal(dna_waste.id, Disposition.DISCARD, reason='Physical DNA cleanup non-retained fraction.'))
    bisulfite = g.input('bisulfite', 'EZ Methylation Direct MagBead kit, half-volume', MaterialKind.CONDITION)
    converted = g.product('converted', 'converted original DNA strands', MaterialKind.MOLECULAR, state=converted_state())
    g.react('bisulfite_convert', 'Bisulfite-convert original marked DNA in half-volume kit', StepKind.CHEMICAL_CONVERSION,
            ((dna_clean, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)), (converted,))
    six_nf_first = g.input('six_nf_first', '6NF first-stage primer; sequence/handle/tail width unassigned', MaterialKind.OLIGO)
    klenow_first = g.input('klenow_first', '50U Klenow first-stage reaction', MaterialKind.CATALYST)
    copied = g.product('dna_first_copy', 'first-stage copied DNA reaction', MaterialKind.MOLECULAR, state=first_copy_state())
    g.react('first_synthesis', 'First synthesis:37°C30min, repeated four more times; five total', StepKind.EXTENSION,
            ((converted, InputRole.TEMPLATE), (six_nf_first, InputRole.STRUCTURAL),
             (klenow_first, InputRole.CATALYST)), (copied,))
    exoi = g.input('exoi', '20U Exonuclease I', MaterialKind.CATALYST)
    digested = g.product('exoi_treated', 'ExoI-treated first-stage reaction; selected substrate spectrum unassigned',
                        MaterialKind.MOLECULAR, state=post_exoi_state())
    hydrolysate = g.product('exoi_products', 'ExoI hydrolysis products; physical partition unassigned', MaterialKind.MOLECULAR)
    g.react('digest_exoi', 'Exonuclease I treatment at37°C', StepKind.DIGESTION,
            ((copied, InputRole.TEMPLATE), (exoi, InputRole.CATALYST)), (digested, hydrolysate))
    g.mark_terminal(Terminal(hydrolysate.id, Disposition.UNRESOLVED, reason='Chemical digestion is distinct from subsequent physical cleanup.'))
    clean_first = g.product('clean_first', '0.75×AMPure-recovered first-stage copies', MaterialKind.MOLECULAR, state=digested.state)
    waste_first = g.product('waste_first', 'first0.75×AMPure non-retained fraction', MaterialKind.WASTE)
    g.react('purify_first', 'Purify ExoI-treated reaction0.75×AMPure', StepKind.CLEANUP,
            ((digested, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (clean_first, waste_first))
    g.mark_terminal(Terminal(waste_first.id, Disposition.DISCARD, reason='Physical first-stage cleanup fraction; species composition unassigned.'))
    six_nf_second = g.input('six_nf_second', '6NF again for second synthesis; full native sequence unassigned', MaterialKind.OLIGO)
    heat = g.input('second_heat', '98°C2min then ice', MaterialKind.CONDITION)
    denatured = g.product('denatured', 'second-stage template/6NF annealing mixture', MaterialKind.MOLECULAR,
                          state=clean_first.state)
    g.react('second_anneal', 'Second-stage6NF setup:98°C2min then ice', StepKind.THERMAL_TREATMENT,
            ((clean_first, InputRole.TEMPLATE), (six_nf_second, InputRole.STRUCTURAL),
             (heat, InputRole.CONDITION)), (denatured,))
    klenow_second = g.input('klenow_second', 'second-stage Klenow', MaterialKind.CATALYST)
    second = g.product('second_copy', 'native second-stage copied DNA', MaterialKind.MOLECULAR, state=second_copy_state())
    g.react('second_synthesis', 'Second-stage Klenow synthesis:37°C90min', StepKind.EXTENSION,
            ((denatured, InputRole.TEMPLATE), (klenow_second, InputRole.CATALYST)), (second,))
    pre_pcr = g.product('dna_pre_pcr', '0.75×AMPure-recovered copied PCR precursor', MaterialKind.MOLECULAR, state=second.state)
    waste_second = g.product('waste_second', 'second0.75×AMPure non-retained fraction', MaterialKind.WASTE)
    g.react('clean_second', 'Purify second-stage copies0.75×AMPure', StepKind.CLEANUP,
            ((second, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (pre_pcr, waste_second))
    g.mark_terminal(Terminal(waste_second.id, Disposition.DISCARD, reason='Physical second-stage cleanup fraction.'))
    pe1 = g.input('pe1', 'PE1.0 PCR primer; native strand unassigned', MaterialKind.OLIGO)
    itag = g.input('itag', 'iTAG PCR primer; family/member/width unassigned', MaterialKind.OLIGO)
    amplified_dna = g.product('dna_amplified', '14-cycle native DNA library amplicons', MaterialKind.LIBRARY, state=DNA_LIBRARY)
    g.react('dna_pcr', 'Amplify native copied DNA:14cycles with PE1.0+iTAG/KAPA', StepKind.PCR,
            ((pre_pcr, InputRole.TEMPLATE), (pe1, InputRole.STRUCTURAL),
             (itag, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (amplified_dna,),
            geometry=PCRGeometry.unavailable(pre_pcr.state, reason='Actual copied DNA input; native6NF and PE1.0/iTAG strand geometry unassigned.'))
    dna_library = g.product('dna_library', DNA_LIBRARY.name, MaterialKind.LIBRARY, state=DNA_LIBRARY)
    final_waste = g.product('final_waste', '0.7×AMPure non-retained library fraction', MaterialKind.WASTE)
    g.react('final_dna_cleanup', 'Purify DNA library0.7×AMPure', StepKind.CLEANUP,
            ((amplified_dna, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (dna_library, final_waste))
    g.mark_terminal(Terminal(final_waste.id, Disposition.DISCARD, reason='Physical final library cleanup fraction.'))
    g.sequence(rna_library, SequencingPlan(library=RNA_LIBRARY, primers=(),
               limitation='Selected EL16 RNA PE125; native primer landing and index acquisition unassigned.'))
    g.sequence(dna_library, SequencingPlan(library=DNA_LIBRARY, primers=(),
               limitation='Selected EL16 DNA PE125; native primer landing and index acquisition unassigned. GCG and CCG interpretation exclusions are computational.'))
    g.declare_reviewed_scope((
        'cell','mcvi_pi','sam','mark_accessibility','marked','rlt_plus','stop_marking','lysate',
        'oligodt_beads','capture_rna','rna_bound','dna_lysate','wash_buffer','wash_capture',
        'washed_rna','dna_washes','combine_dna','dna_fraction','rt','rt_primer','reverse_transcribe',
        'rna_first','tso','template_switch','rna_switched','ispcr','kapa','preamp_rna','rna_amplified',
        'ampure','clean_cdna','cdna_clean','cdna_amount','sample_cdna','cdna_aliquot','transposome',
        'tn5','tagment_rna','nextera_pcr','rna_tagged','rna_index_f','rna_index_r','rna_index_pcr','rna_library',
        'clean_dna','dna_clean','bisulfite','bisulfite_convert','converted','six_nf_first',
        'klenow_first','first_synthesis','dna_first_copy','exoi','digest_exoi','exoi_treated',
        'purify_first','clean_first','six_nf_second','second_heat','second_anneal','denatured','klenow_second',
        'second_synthesis','second_copy','clean_second','dna_pre_pcr','pe1','itag','dna_pcr',
        'dna_amplified','final_dna_cleanup','dna_library'), SOURCE)
    return g.require_publishable()
