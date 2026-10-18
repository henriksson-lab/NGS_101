"""Protocol-owned oxBS declarations; historical and later implementations stay distinct."""
from types import SimpleNamespace
from chemdraw import Construct, Scene
from pcr import PCRGeometry
from base_conversion import BasePath
from protocol_model import role, molecule
from protocol_graph import (Disposition, Terminal, Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
import seqprimers as sp


def oxbs_path(base: str) -> BasePath:
    """Chemical intermediates followed by PCR readout, not direct base detection."""
    if base == "5hmC":
        return BasePath(("5hmC", "5fC", "U", "T"))
    if base == "5mC":
        return BasePath(("5mC", "5mC", "C"))
    if base == "C":
        return BasePath(("C", "U", "T"))
    raise ValueError("oxBS input must be C, 5mC or 5hmC")


SOURCE = ('<a href="https://doi.org/10.1126/science.1220671">Booth et al., Science (2012)</a>; '
          'separately labelled reference implementation from '
          '<a href="https://doi.org/10.1038/nprot.2013.115">Booth et al., Nature Protocols (2013)</a>.')
DETAIL = 'Published 2013 protocol PMC3919000, Steps 6–22 and 32–36'


def _state(name, *, controls=False, digestion=False):
    sc = molecule(name, (role(name),))
    if controls:
        sc.strand('sequencing_controls', [role('separate sequencing-control DNA molecules')],
                  label='spiked sequencing controls; strand-specific chemistry unassigned')
    if digestion:
        sc.strand("digestion_control", [role("separate digestion-control DNA")])
    return sc


def _library(name, signal):
    return Construct([role('native TruSeq end: selected historical strands unassigned', inferred=True),
                      role(signal),
                      role('opposite native TruSeq end: index architecture unassigned', inferred=True)], name=name)


OX_LIBRARY = _library('2013 oxBS reference library', 'original C→T; 5mC→C; 5hmC→T after PCR')
BS_LIBRARY = _library('2013 matched BS reference library', 'original C→T; 5mC/5hmC→C after PCR')
FINAL_LIBRARY = OX_LIBRARY
SEQ_PRIMERS = ()


def protocol_graph():
    g = ProtocolGraph('oxBS: separately labelled detailed 2013 reference', authority=DETAIL)
    dna = g.input('genomic_dna', '2013 reference: 100 ng–1 µg genomic DNA', MaterialKind.BIOLOGICAL,
                  state=_state('genomic DNA with coexisting C, 5mC and 5hmC sites'), primary=True)
    shear = g.input('sonication', 'sonication, recommended whole-genome fragment range 100–400 bp', MaterialKind.CONDITION)
    fragmented = g.product('fragmented_dna', 'sonicated genomic DNA', MaterialKind.MOLECULAR,
                           state=_state('sonicated genomic DNA'))
    g.react('fragment', 'Sonicate genomic DNA (Step 6)', StepKind.FRAGMENTATION,
            ((dna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragmented,))
    controls = g.input('sequencing_controls', 'sequencing-control DNA spike-in, 0.5% by weight (Step 7)',
                       MaterialKind.BIOLOGICAL, state=_state('separate sequencing-control DNA molecules'))
    mixed = g.product('control_mix', 'genomic DNA plus sequencing controls', MaterialKind.MOLECULAR,
                      state=_state('genomic DNA', controls=True))
    g.react('spike_controls', 'Add sequencing controls before library preparation', StepKind.GENERIC,
            ((fragmented, InputRole.TEMPLATE), (controls, InputRole.SPIKE_IN)), (mixed,), changes_molecule=False)
    kit = g.input('truseq_kit', 'Illumina TruSeq library preparation (Step 8)', MaterialKind.CONDITION)
    adapted = g.product('adapted_dna', 'adapted DNA before replacement bead cleanup', MaterialKind.MOLECULAR,
                        state=_state('native TruSeq-adapted genomic DNA: strand geometry unassigned', controls=True))
    g.react('native_library_prep', 'Prepare TruSeq libraries; omit final kit bead cleanup (Step 8)', StepKind.MAGIC,
            ((mixed, InputRole.TEMPLATE), (kit, InputRole.CONDITION)), (adapted,),
            reason='Step 8 delegates construction to a historical kit manual. Selected adapter strands, '
                   'end geometry and index member are unassigned; no canonical Y adapter is substituted.')
    beads = g.input('step9_beads', 'AMPure XP: 1:1 for sonicated DNA; two 80% ethanol washes, dry 30 min, 22.5 µl water', MaterialKind.SELECTION)
    ox = g.product('ox_aliquot', 'Step 9 oxBS library aliquot', MaterialKind.MOLECULAR, state=adapted.state)
    bs = g.product('bs_aliquot', 'Step 9 nonoxidized BS aliquot set aside', MaterialKind.MOLECULAR, state=adapted.state)
    waste = g.product('step9_waste', 'Step 9 supernatant and washes', MaterialKind.WASTE)
    g.react('clean_and_split', 'Purify adapted DNA and reserve matched BS sample (Step 9)', StepKind.CLEANUP,
            ((adapted, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (ox, bs, waste))
    g.discard(waste, 'AMPure supernatant and washes physically removed')
    digestion = g.input('digestion_control', 'digestion-control DNA, 0.5% by weight (Step 10)', MaterialKind.BIOLOGICAL)
    oxmix_scene = _state('native adapted genomic DNA', controls=True)
    oxmix_scene.strand('digestion_control', [role('separate digestion-control DNA')])
    oxmix = g.product('ox_control_mix', 'oxBS aliquot plus digestion control', MaterialKind.MOLECULAR, state=oxmix_scene)
    g.react('spike_digestion_control', 'Add digestion control to oxBS aliquot (Step 10)', StepKind.GENERIC,
            ((ox, InputRole.TEMPLATE), (digestion, InputRole.SPIKE_IN)), (oxmix,), changes_molecule=False)
    p6 = g.input('p6_column', 'water-equilibrated Bio-Spin P-6 column', MaterialKind.SELECTION)
    exchanged = g.product('water_exchanged', 'DNA in water, P-6 eluate', MaterialKind.MOLECULAR, state=oxmix.state)
    p6w = g.product('p6_retained', 'column-retained buffer components', MaterialKind.WASTE)
    g.react('exchange_buffer', 'Remove buffer/ethanol contamination; retain eluate (Steps 11–13)', StepKind.CLEANUP,
            ((oxmix, InputRole.TEMPLATE), (p6, InputRole.SELECTION)), (exchanged, p6w))
    g.discard(p6w, 'buffer components retained by column; DNA eluate recovered')
    naoh = g.input('naoh', 'NaOH denaturation: 37°C 30 min, then ice 5 min', MaterialKind.CONDITION)
    ss = g.product('denatured_dna', 'denatured genomic library strands', MaterialKind.MOLECULAR,
                   state=_state('population of separate native library strands: C, 5mC, 5hmC', controls=True, digestion=True))
    opposite = g.product('denatured_opposite', 'opposite denatured genomic library strands', MaterialKind.MOLECULAR,
                         state=_state('opposite native library strand population'))
    g.react('denature', 'Alkaline-denature adapted DNA (Steps 14–16)', StepKind.STRAND_SEPARATION,
            ((exchanged, InputRole.TEMPLATE), (naoh, InputRole.CONDITION)), (ss, opposite))
    oxidant = g.input('perruthenate', 'potassium perruthenate; on ice 60 min', MaterialKind.CONDITION)
    oxidized = g.product('oxidized_dna', 'oxidized library strands', MaterialKind.MOLECULAR,
                         state=_state('native library strands: C and 5mC retained; 5hmC→5fC', controls=True, digestion=True))
    g.react('oxidize', 'Oxidize 5hmC to 5fC (Steps 17–18)', StepKind.CHEMICAL_CONVERSION,
            ((ss, InputRole.TEMPLATE), (opposite, InputRole.TEMPLATE), (oxidant, InputRole.CONDITION)), (oxidized,))
    recovered = g.product('oxidized_eluate', 'postoxidation P-6 eluate', MaterialKind.MOLECULAR, state=oxidized.state)
    ow = g.product('oxidation_column_waste', 'postoxidation column-retained components', MaterialKind.WASTE)
    g.react('recover_oxidized', 'Recover DNA through water-equilibrated P-6 (Step 19)', StepKind.CLEANUP,
            ((oxidized, InputRole.TEMPLATE), (p6, InputRole.SELECTION)), (recovered, ow))
    g.discard(ow, 'column-retained components separated from recovered DNA eluate')
    bisulfite = g.input('bisulfite', 'EpiTect FFPE kit with modified conversion program (Steps 20–22)', MaterialKind.CONDITION)
    forward = g.input('truseq_pcr_forward', 'TruSeq PCR primer cocktail: forward species unassigned', MaterialKind.OLIGO, inferred=True)
    reverse = g.input('truseq_pcr_reverse', 'TruSeq PCR primer cocktail: reverse species unassigned', MaterialKind.OLIGO, inferred=True)
    polymerase = g.input('pfuturbo_cx', 'Agilent PfuTurbo Cx DNA polymerase', MaterialKind.CATALYST)
    cycles = g.input('pcr_cycles', '12–15 cycles: 95°C 30 s / 60°C 20 s / 72°C 45 s', MaterialKind.CONDITION)
    cleanup = g.input('post_pcr_beads', '50 µl AMPure XP; 22.5 µl water eluate (Step 34)', MaterialKind.SELECTION)
    for branch, precursor, signal, final in (
            ('oxbs', recovered, 'original C→U; 5mC retained; original 5hmC→5fC→U', OX_LIBRARY),
            ('bs', bs, 'original C→U; 5mC retained; 5hmC→stable methylsulfonate adduct (reads C)', BS_LIBRARY)):
        converted = g.product(branch + '_converted', branch + ' physical product before PCR', MaterialKind.MOLECULAR,
                              state=_state('native single-strand library: ' + signal, controls=True, digestion=branch == 'oxbs'))
        g.react(branch + '_bisulfite', 'Bisulfite-convert ' + branch + ' aliquot (Steps 20–22)', StepKind.CHEMICAL_CONVERSION,
                ((precursor, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)), (converted,))
        amplified = g.product(branch + '_amplified', branch + ' amplified library', MaterialKind.MOLECULAR, state=final)
        control_products = g.product(branch + '_control_products', branch + ' control-derived products: native identities unassigned',
                                     MaterialKind.MOLECULAR,
                                     state=_state('separate control-derived product population; exact amplification unassigned'))
        g.react(branch + '_pcr', 'Amplify 2 µl converted DNA in 50 µl (Steps 32–33)', StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST), (cycles, InputRole.CONDITION)), (amplified, control_products),
                geometry=PCRGeometry.unavailable(converted.state,
                    reason='Actual converted strands are the PCR substrate. Native TruSeq primer/end/index identities remain unassigned.'))
        g.mark_terminal(Terminal(control_products.id, Disposition.UNRESOLVED,
            'Controls are not silently removed: their strand-specific amplification and Step 23–31 QC partition remain bounded. The final construct depicts the genomic library signal.'))
        library = g.product(branch + '_library', branch + ' purified reference library', MaterialKind.LIBRARY, state=final)
        w = g.product(branch + '_post_pcr_waste', branch + ' bead supernatant/washes', MaterialKind.WASTE)
        g.react(branch + '_cleanup', 'AMPure-purify PCR products (Step 34)', StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (cleanup, InputRole.SELECTION)), (library, w))
        g.discard(w, 'supernatant and washes physically removed during bead recovery')
        g.sequence(library, SequencingPlan(final, (), limitation='2013 reference MiSeq/HiSeq endpoint; acquired roles, cycles and native primer sites unassigned.'))
    g.annotate(Highlight('native_library_prep', 'Implementation boundary',
        'The detailed graph is the separately labelled 2013 sonicated TruSeq reference, not the 2012 oxRRBS recipe. '
        'Control DNA species remain independent molecules. Their exact modification patterns and Step 23–31 digestion QC '
        'are outside this bounded molecular reconstruction; QC primers are not genomic-library PCR primers.'))
    g.annotate(Highlight('bs_bisulfite', 'Physical product versus sequencing signal',
        'Ordinary BS converts 5hmC to a stable methylsulfonate adduct, read as C after PCR. '
        '5hmC abundance requires matched BS-minus-oxBS comparison; one oxBS molecule does not distinguish original C from 5hmC.'))
    g.annotate(SourceDiscrepancy('oxidize', 'Published Supplementary Fig. 1b substrate wording',
        (SourceClaim('2013 Supplementary Fig. 1b, opening/final clauses', 'DNA made with d5fCTP; 5fC conversion to U'),
         SourceClaim('2013 Supplementary Fig. 1b, middle clause', 'names 5hmC as the conversion substrate')),
        'Model uses 5fC→U after oxidation, following the surrounding caption and main protocol.',
        'The surrounding figure describes the 5fC experiment.',
        'The inconsistent 5hmC clause appears erroneous; intended wording remains unresolved.'))
    g.declare_reviewed_scope((
        'genomic_dna', 'sonication', 'fragment',
        'sequencing_controls', 'spike_controls', 'truseq_kit', 'native_library_prep', 'step9_beads',
        'clean_and_split', 'digestion_control', 'spike_digestion_control', 'p6_column', 'exchange_buffer',
        'naoh', 'denature', 'perruthenate', 'oxidize', 'recover_oxidized', 'bisulfite', 'truseq_pcr_forward', 'truseq_pcr_reverse',
        'pfuturbo_cx', 'pcr_cycles', 'post_pcr_beads', 'oxbs_bisulfite', 'bs_bisulfite',
        'oxbs_pcr', 'bs_pcr', 'oxbs_cleanup', 'bs_cleanup', 'oxbs_library', 'bs_library'), DETAIL + '; defining 2012 main article; 2013 Supplementary Fig. 1b')
    return g.require_publishable()


SPEC = SimpleNamespace(
    TITLE='oxBS-seq — oxidative bisulfite sequencing', NOTES='01_oxbs_seq.html', SOURCE=SOURCE,
    SUMMARY='Oxidation converts 5hmC to 5fC before bisulfite. PCR reports original C and 5hmC as T, and 5mC as C. Matched ordinary BS reports 5mC plus 5hmC as C.',
    CAVEAT='INFERRED native adapter end roles — historical TruSeq strands, index architecture and primer sites remain unassigned. The graph is the 2013 reference implementation. The original 2012 oxRRBS library construction remains outside this graph pending its original SOM. No default index length, ligation geometry or acquired read is inferred from kit compatibility.',
    oligos=lambda: (), FINAL_LIBRARY=OX_LIBRARY, SEQ_PRIMERS=(),
    FINAL_CAPTION='INFERRED native end roles — 2013 oxBS PCR endpoint; cytosine outcomes are population interpretations, not separated base-class molecules.',
    EXTRA_MOLECULAR_STATES=(('Matched ordinary BS reference endpoint', BS_LIBRARY,
        'INFERRED native end roles — original 5mC and 5hmC both report C after PCR.'),),
    SEQUENCING_ENDING='2013 Step 35 names MiSeq or HiSeq. Native library primer sites, read/index acquisition and cycles remain unassigned within the inspected scope; the original 2012 read layout is separately unresolved.',
    READ_SPECS=tuple(sp.read_spec(read, 'unassigned', 'Acquisition and cycles unassigned for this reference scope', 'undetermined')
                     for read in ('Read 1', 'Read 2', 'Index 1 (i7)', 'Index 2 (i5)')),
    protocol_graph=protocol_graph)
