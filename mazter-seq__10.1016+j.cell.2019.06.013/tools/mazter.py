"""MAZTER-seq: source-qualified cleavage events and library preparation.

The selected ACA event is not a population cleavage or methylation predictor.
Historical Table S7 adapter polymers and end chemistry remain unassembled.
"""
from dataclasses import dataclass
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Scene, Segment, workflow_from_protocol_graph, workflow_panel
from page import head, info
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R,
                            StepKind as S, Terminal, Disposition, SourceItem)

PAGE_FILE = 'mazter-seq.html'
SOURCE = ('Garcia-Campos et al. 2019, DOI10.1016/j.cell.2019.06.013; published '
          'Results p732 and STAR Methods pp e3–e4, institutional VOR PDF')
ORIGINAL = ('https://weizmann.elsevierpure.com/ws/portalfiles/portal/111306132/'
            'ss_Cell_DecipheringThe-m6ACode_PV2019.pdf')
PCR_CYCLES = 13


@dataclass(frozen=True)
class SelectedCut:
    """One caller-selected unmodified ACA scission, preserving every RNA base.

    Ends are deliberately unassigned. Unmodified ACA is not guaranteed to cleave,
    and non-ACA scission is not ruled out; the source discusses both limitations.
    """
    rna: str
    position: int
    first_a_methylated: bool = False

    def __post_init__(self):
        if not self.rna or set(self.rna) - set('ACGU'):
            raise ValueError('An explicit RNA substrate must use A/C/G/U')
        if not isinstance(self.position, int) or not 0 < self.position <= len(self.rna) - 3:
            raise ValueError('A selected internal ACA cut must leave both RNA fragments')
        if self.rna[self.position:self.position + 3] != 'ACA':
            raise ValueError('This selected-event model represents an ACA site only')
        if self.first_a_methylated:
            raise ValueError('Do not model a protected first-A m6A site as this susceptible event')

    @property
    def fragments(self):
        return self.rna[:self.position], self.rna[self.position:]


def selected_event_scene():
    """Symbolic source event; no study transcript or quantitatively complete digest."""
    scene = Scene()
    scene.strand('upstream', [Segment('RNA before selected ACA', '[upstream RNA]', placeholder=True)],
                 label='upstream RNA; new 3′ end chemistry unassigned')
    scene.strand('downstream', [Segment('recognition motif', 'ACA'),
                               Segment('remaining RNA', '[downstream RNA]', placeholder=True)],
                 label='downstream RNA; new 5′ end chemistry unassigned')
    scene.strand('uncut', [Segment('uncut molecules', '[uncut RNA]', placeholder=True)],
                 label='uncut/protected RNA remains in the reaction mixture')
    scene.note('downstream', 'Selected cleavage lies immediately before ACA; all RNA bases are retained across the two fragments.')
    return tuple(scene.rows())


def unreviewed_ligation(reason):
    scene = Scene()
    scene.strand('pending_adapter', [Segment('unreviewed Table S7 adapter end',
                  '[adapter end unreviewed; Table S7 pending]', placeholder=True)])
    scene.note('pending_adapter', reason)
    return AdapterLigationGeometry('unavailable', 'unavailable', (scene,), reason)


def protocol_graph():
    g = ProtocolGraph('MAZTER-seq — published no-IP library branch', authority=SOURCE,
                      representation='published_instance')
    total = g.input('total_rna', 'Extracted total RNA; source-specific extraction is upstream of this route', K.MOLECULAR, primary=True)
    beads = g.input('polya_selection', 'Source oligo(dT) mRNA enrichment; small-input mRNA-DIRECT or large-input GenElute', K.SELECTION)
    rna = g.product('mrna', 'Poly(A)-selected RNA containing modification/sequence/structure-dependent populations', K.MOLECULAR)
    reject = g.product('polya_reject', 'Physical mRNA-selection reject and washes', K.WASTE)
    g.react('select_mrna', 'Select polyadenylated RNA using the source input-dependent kit', S.SELECTION,
            ((total, R.TEMPLATE), (beads, R.SELECTION)), (rna, reject))
    g.discard(reject, 'Physical selection fractions; no claim of perfect depletion')
    heat = g.input('heat', '70°C for 2 min, then ice; source 100 ng poly(A) RNA', K.CONDITION)
    denatured = g.product('denatured', 'Heat-treated mRNA; subsequent structure and cleavage probabilities unassigned', K.MOLECULAR)
    g.react('denature', 'Heat-denature before MazF', S.THERMAL_TREATMENT,
            ((rna, R.TEMPLATE), (heat, R.CONDITION)), (denatured,))
    mazf = g.input('mazf', 'Takara 2415A MazF, 20 U; source buffer, DMSO and RNase inhibitor', K.CATALYST)
    digest = g.product('digest', 'Cleaved plus uncut/protected RNA population', K.MOLECULAR,
                       state=selected_event_scene())
    g.react('cleave', '37°C for 30 min; one ACA event shown, not a complete digest prediction', S.DIGESTION,
            ((denatured, R.TEMPLATE), (mazf, R.CATALYST)), (digest,))

    def cleanup(identifier, template, label):
        retained = g.product(identifier, label, K.MOLECULAR)
        waste = g.product(identifier + '_waste', 'Physical purification supernatant and washes', K.WASTE)
        g.react(identifier + '_cleanup', 'Source MyOne SILANE/RLT/ethanol purification', S.CLEANUP,
                ((template, R.TEMPLATE),), (retained, waste))
        g.discard(waste, 'Physical cleanup; exact population recovery not predicted')
        return retained

    clean_rna = cleanup('clean_rna', digest, 'Purified RNA including both fragments and uncut molecules')
    repair = g.input('repair_enzymes', 'FastAP, T4 PNK and Turbo DNase in source FastAP/PNK buffer mixture', K.CATALYST)
    repaired = g.product('repaired', 'Source-described end-repaired/dephosphorylated RNA; exact end chemistry unassembled', K.MOLECULAR)
    g.react('repair', 'End-repair RNA before adapter ligation', S.END_REPAIR,
            ((clean_rna, R.TEMPLATE), (repair, R.CATALYST)), (repaired,))
    repaired = cleanup('repair_clean', repaired, 'Purified end-repaired RNA')
    adapter = g.input('rna_adapter', '3iLL RNA adapter, 20 pmol; Table S7 polymer/modifications unreviewed', K.ADAPTER)
    preligation_heat = g.input('preligation_heat', 'RNA plus 3iLL: 70°C for 2 min, then snap-cool on ice', K.CONDITION)
    heated_rna = g.product('preligation_rna', 'Heat-treated RNA in the adapter-containing mixture; no isolation or binding geometry asserted', K.MOLECULAR)
    heated_adapter = g.product('preligation_adapter', 'Heat-treated 3iLL in the same mixture; no isolation or binding geometry asserted', K.ADAPTER)
    g.react('preligation_denature', 'Heat and snap-cool RNA together with its adapter before adding ligase', S.THERMAL_TREATMENT,
            ((repaired, R.TEMPLATE), (adapter, R.STRUCTURAL), (preligation_heat, R.CONDITION)),
            (heated_rna, heated_adapter))
    ligase = g.input('rna_ligase', 'T4 RNA ligase 1, NEB M0437M', K.CATALYST)
    tagged = g.product('rna_tagged', 'RNA ligated at its 3′ end to source 3iLL adapter; possible unligated populations unquantified', K.MOLECULAR)
    g.react('rna_ligation', 'Source RNA adapter ligation, 23°C for 75 min', S.LIGATION,
            ((heated_rna, R.TEMPLATE), (heated_adapter, R.STRUCTURAL), (ligase, R.CATALYST)), (tagged,),
            geometry=unreviewed_ligation('Table S7 adapter polymer/modifications and exact ligation-ready ends have not been inspected.'))
    tagged = cleanup('rna_ligation_clean', tagged, 'Purified adapter-ligated RNA population')
    others = g.input('other_barcoded_rna', 'Corresponding source adapter-ligated samples for pooling; exact barcode membership unassigned', K.MOLECULAR)
    pooled = g.product('pooled_rna', 'Pooled source barcoded RNA, no-IP branch', K.MOLECULAR)
    g.react('pool', 'Pool ligated RNA; proceed directly to library preparation for no-IP route', S.POOL,
            ((tagged, R.TEMPLATE), (others, R.TEMPLATE)), (pooled,))
    primer = g.input('rt_primer', 'rTd DNA cDNA primer complementary to 3ILL; Table S7 polymer unreviewed', K.OLIGO)
    rt = g.input('rt', 'SuperScript III, 50°C for 1 h; omit heat inactivation', K.CATALYST)
    hybrid = g.product('hybrid', 'RNA/cDNA hybrids plus residual primers; full geometry unassembled', K.MOLECULAR)
    g.react('reverse_transcribe', 'Adapter-primed reverse transcription; preserve RNA/cDNA hybrid', S.REVERSE_TRANSCRIPTION,
            ((pooled, R.TEMPLATE), (primer, R.STRUCTURAL), (rt, R.CATALYST)), (hybrid,))
    exosap = g.input('exosap', 'ExoSAP-IT, 37°C for 12 min; stop with EDTA', K.CATALYST)
    treated = g.product('primer_digested', 'RNA/cDNA hybrid plus primer-digestion products; no physical debris removal yet', K.MOLECULAR)
    g.react('digest_primers', 'Digest remaining primers without heat denaturation', S.DIGESTION,
            ((hybrid, R.TEMPLATE), (exosap, R.CATALYST)), (treated,))
    alkali = g.input('alkali', 'Source NaOH, 70°C for 12 min, then HCl neutralization', K.CONDITION)
    hydrolysed = g.product('hydrolysed', 'Single-stranded cDNA plus hydrolysed RNA and primer-derived products', K.MOLECULAR)
    g.react('remove_rna', 'Hydrolyse RNA; neutralize without calling hydrolysis products discarded', S.CHEMICAL_CONVERSION,
            ((treated, R.TEMPLATE), (alkali, R.CONDITION)), (hydrolysed,))
    cdna = cleanup('clean_cdna', hydrolysed, 'Purified single-stranded cDNA')
    dna_adapter = g.input('dna_adapter', '5iLL-22 DNA adapter, 50 pmol; Table S7 polymer/modifications unreviewed', K.ADAPTER)
    cdna_ligase = g.input('cdna_ligase', 'T4 RNA ligase 1; 4% DMSO, 1 mM ATP, 21% PEG8000', K.CATALYST)
    adapted = g.product('adapted_cdna', 'cDNA with second adapter at its 3′ end; exact junction unassembled', K.MOLECULAR)
    g.react('cdna_ligation', 'Ligate source DNA adapter to cDNA, 23°C for 3 h', S.LIGATION,
            ((cdna, R.TEMPLATE), (dna_adapter, R.STRUCTURAL), (cdna_ligase, R.CATALYST)), (adapted,),
            geometry=unreviewed_ligation('Table S7 DNA adapter polymer/modifications and cDNA terminal geometry remain unreviewed.'))
    adapted = cleanup('adapted_clean', adapted, 'Purified dual-adapter cDNA population')
    forward = g.input('pcr_f', 'Library PCR primer F, Table S7', K.OLIGO)
    reverse = g.input('pcr_r', 'Library PCR primer R, Table S7', K.OLIGO)
    kapa = g.input('kapa', 'KAPA HiFi HotStart ReadyMix KK2601', K.CATALYST)
    amplified = g.product('amplified', 'Amplified library population; exact template/primer/amplicon geometry unassembled', K.MOLECULAR)
    g.react('pcr', f'Library PCR, {PCR_CYCLES} cycles', S.PCR,
            ((adapted, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL), (kapa, R.CATALYST)), (amplified,),
            geometry=PCRGeometry.unavailable(Scene(), reason='Source Table S7 primers and complete adapter-derived template are unassembled.'))
    library = g.product('library', 'AMPure XP-purified source library for NextSeq 500 paired-end sequencing', K.LIBRARY)
    pcr_waste = g.product('pcr_waste', 'Physical AMPure supernatant/washes', K.WASTE)
    g.react('final_cleanup', 'AMPure cleanup, Qubit and TapeStation assessment', S.CLEANUP,
            ((amplified, R.TEMPLATE),), (library, pcr_waste))
    g.discard(pcr_waste, 'Physical final cleanup; no perfect removal of particular molecular classes asserted')
    g.mark_terminal(Terminal(library.id, Disposition.UNRESOLVED,
                            'Published paired-end NextSeq500 named; Table S7 polymers, exact run-primer sites and cycle/index plan remain unresolved.'))
    g.declare_reviewed_scope(('total_rna', 'polya_selection', 'select_mrna', 'heat', 'denature',
                             'mazf', 'cleave', 'clean_rna_cleanup', 'repair_enzymes', 'repair',
                             'repair_clean_cleanup', 'rna_adapter', 'preligation_heat', 'preligation_denature', 'rna_ligase', 'rna_ligation',
                             'rna_ligation_clean_cleanup', 'other_barcoded_rna', 'pool',
                             'rt_primer', 'rt', 'reverse_transcribe', 'exosap', 'digest_primers',
                             'alkali', 'remove_rna', 'clean_cdna_cleanup', 'dna_adapter',
                             'cdna_ligase', 'cdna_ligation', 'adapted_clean_cleanup',
                             'pcr_f', 'pcr_r', 'kapa', 'pcr', 'final_cleanup'), SOURCE)
    g.add_source_item(SourceItem('table_s7', 'Table S7 exact adapter, primer and barcode polymers',
                                'Published STAR Methods explicitly refer to Table S7; publisher supplement index is inaccessible.'))
    return g


def render_page():
    return ('<!doctype html><html lang="en">' + head('MAZTER-seq — source-limited model') +
            '<main class="wrap"><h1>MAZTER-seq — source-limited model</h1>' +
            info('Private incomplete model. One ACA event and the no-IP route are represented. '
                 'Cleavage is sequence/structure dependent; protected or uncut molecules remain. '
                 'Source Table S7 polymers, end chemistry and run geometry are unassembled.') +
            workflow_panel(workflow_from_protocol_graph(protocol_graph(), 'library')) +
            '<p>Optional two-round m6A IP and calibration branches remain outside this route. '
            '<a href="' + ORIGINAL + '">Published original, Results p732 and STAR Methods e3–e4</a>.</p></main></html>')


def selftest():
    from checks import Check
    c = Check()
    cut = SelectedCut('GGGACAUU', 3)
    c('selected internal cut preserves RNA order and begins downstream with ACA',
      (''.join(cut.fragments), cut.fragments[1]), (cut.rna, 'ACAUU'))
    for label, values in [('protected m6A', ('GGGACAUU', 3, True)),
                          ('non-ACA event', ('GGGACAUU', 2)),
                          ('empty upstream fragment', ('ACAUU', 0)),
                          ('DNA substituted for RNA', ('GGGACATT', 3))]:
        try:
            SelectedCut(*values)
        except ValueError:
            rejected = True
        else:
            rejected = False
        c('selected-event boundary rejects ' + label, rejected, True)
    c('unreviewed polymers and readout prevent public finished status', protocol_graph().publishable, False)
    c.report()


if __name__ == '__main__':
    selftest()
