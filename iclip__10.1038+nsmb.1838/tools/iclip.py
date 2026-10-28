"""Source-scoped partial iCLIP model; intentionally not a publishable library."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Scene, Segment, revcomp, feature
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (ProtocolGraph, MaterialKind, InputRole, StepKind,
                            Terminal, Disposition, MolecularRelation)
from circular import circularize_ssdna, CircularizationGeometry
from illumina import TRUSEQ_P5_FULL, PE_PCR_PRIMER_2

SOURCE = 'König et al. 2010 doi:10.1038/nsmb.1838, publisher Supplementary Methods pp1–4'
L3_WRITTEN = '5′-phosphate-UGAGAUCGGAAGAGCGGTTCAG-3′-Puromycin'
L3 = 'UGAGAUCGGAAGAGCGGTTCAG'
RT_PRIMERS = ('NNNCAAGATCGGAAGAGCGTCGTGGATCCTGAACCGCTC',
              'NNNGAAGATCGGAAGAGCGTCGTGGATCCTGAACCGC',
              'NNNTGAGATCGGAAGAGCGTCGTGGATCCTGAACCGCTC')
LINEARIZATION_PRIMER = 'GTTCAGGATCCACGACGCTCTTCAAAA'
PCR_PRIMERS = (PE_PCR_PRIMER_2, TRUSEQ_P5_FULL)
BOUNDARY = ('This capped model ends at the original replicate-1 cDNA circle. '
            'BamHI reopening, source-specific PCR reconstruction, variant discrepancies '
            'and exact acquired sequencing-primer geometry remain deferred.')


def linear_cdna():
    # The terminal nine bases pair L3 positions 12–20 (1-based); L3's terminal AG remains beyond the paired interval.
    complement = revcomp(L3.replace('U', 'T'))
    assert RT_PRIMERS[0][-9:] == complement[2:11]
    return Construct([
        Segment('random barcode', 'NNN', 'umi', placeholder=True,
                feature=feature('iclip_umi', 'umi', 'random')),
        Segment('experiment CA', 'CA', 'cbc',
                feature=feature('iclip_replicate', 'sample_index', 'fixed')),
        Segment('inverted PCR handles and BamHI site', RT_PRIMERS[0][5:-9]),
        Segment('RT landing suffix', RT_PRIMERS[0][-9:]),
        Segment('copied L3 upstream', revcomp(L3[:11].replace('U', 'T'))),
        Segment('RNA copy to crosslink stop', '[variable cDNA]', placeholder=True),
    ], name='replicate-1 single-stranded cDNA; 5′ phosphate retained')



def l3_scene():
    sc = Scene()
    sc.strand('RNA adapter', [Segment('L3', L3)], label='L3',
              mod5='phosphate', mod3='Puromycin')
    sc.ligation_end('RNA adapter', 'L3', kind='single-strand',
                    polarity='5-prime', owner='adapter')
    return sc


def protocol_graph():
    """Starts at the explicitly scoped post-IP RNA; stops before unsupported RT pairing."""
    g = ProtocolGraph('iCLIP — source-limited RNA linker ligation', authority=SOURCE,
                      representation='published_instance')
    rna = g.input('bound_rna', 'post-IP RNase-I-trimmed hnRNP C-bound RNA',
                  MaterialKind.MOLECULAR, primary=True,
                  state=Construct([Segment('RNA fragment', '[variable RNA]', placeholder=True)]))
    sap = g.input('sap', 'shrimp alkaline phosphatase', MaterialKind.CATALYST)
    oh = g.product('rna_3oh', 'RNA with dephosphorylated 3-prime end', MaterialKind.MOLECULAR,
                   state=Construct([Segment('RNA fragment', '[variable RNA]', placeholder=True)]))
    g.react('dephosphorylate', 'Remove RNA 3-prime phosphate (37 °C, 10 min)',
            StepKind.DEPHOSPHORYLATION, ((rna, InputRole.TEMPLATE), (sap, InputRole.CATALYST)), (oh,))
    adapter = g.input('l3', L3_WRITTEN, MaterialKind.ADAPTER, state=tuple(l3_scene().rows()))
    ligase = g.input('rna_ligase', 'RNA ligase (NEB)', MaterialKind.CATALYST)
    ligated = g.product('ligated_rna', 'RNA–L3 covalent junction; L3 3-prime Puromycin retained',
                        MaterialKind.MOLECULAR,
                        state=Construct([Segment('RNA fragment', '[variable RNA]', placeholder=True),
                                         Segment('L3', L3)]))
    g.react('ligate_l3', 'Ligate L3 5-prime phosphate to RNA 3-prime OH (16 °C overnight)',
            StepKind.LIGATION, ((oh, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                                (ligase, InputRole.CATALYST)), (ligated,),
            geometry=AdapterLigationGeometry.known('single-strand', l3_scene(),
                statement='Supplementary Methods p2 prints 5-prime phosphate L3 and 3-prime Puromycin.'))
    pk = g.input('proteinase_k', 'proteinase K / urea extraction', MaterialKind.CATALYST)
    recovered = g.product('recovered_rna', 'peptide-bearing RNA–L3 recovered from membrane',
                          MaterialKind.MOLECULAR, state=ligated.state)
    g.react('recover_rna', 'Gel/membrane recovery followed by proteinase K digestion',
            StepKind.DEPROTEINIZATION, ((ligated, InputRole.TEMPLATE), (pk, InputRole.CATALYST)), (recovered,))
    primer = g.input('rt_primer', 'replicate-1 5-prime phosphorylated RT primer',
                     MaterialKind.OLIGO, state=Construct(list(linear_cdna())[:4]))
    rt = g.input('superscript_iii', 'SuperScript III', MaterialKind.CATALYST)
    cdna = g.product('cdna', 'crosslink-truncated cDNA after RNase A and size recovery',
                     MaterialKind.MOLECULAR, state=linear_cdna())
    g.react('reverse_transcribe', 'SuperScript III; heat/RNase A; recover 100–175 or 175–350 nt cDNA',
            StepKind.REVERSE_TRANSCRIPTION,
            ((recovered, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (cdna,))
    ligase2 = g.input('circligase', 'CircLigase II with MnCl2', MaterialKind.CATALYST)
    closed = circularize_ssdna(linear_cdna(), five_prime_phosphate=True)
    circle = g.product('cdna_circle', 'replicate-1 cDNA circle', MaterialKind.MOLECULAR,
                       chemical_state=closed, state=tuple(closed.rows()))
    g.react('circularize', 'CircLigase II intramolecular closure (60 °C, 1 h)',
            StepKind.CIRCULARIZATION,
            ((cdna, InputRole.TEMPLATE), (ligase2, InputRole.CATALYST)), (circle,),
            geometry=CircularizationGeometry.known('single-strand intramolecular ligation',
                        statement='Supplementary Methods p3: phosphorylated RT primer and CircLigase II closure.'),
            relations=(MolecularRelation.CIRCLE_CLOSURE,))
    g.mark_terminal(Terminal(circle.id, Disposition.UNRESOLVED, BOUNDARY))
    g.declare_reviewed_scope(('bound_rna', 'sap', 'dephosphorylate', 'l3', 'rna_ligase',
                             'ligate_l3', 'proteinase_k', 'recover_rna', 'rt_primer', 'superscript_iii',
                             'reverse_transcribe', 'circligase', 'circularize'), SOURCE)
    return g


def selftest():
    from checks import Check
    c = Check()
    # Regression for the exact, non-terminal RT landing interval in the printed L3.
    expected = revcomp(L3.replace('U', 'T'))
    c('printed RT1 lands exactly within L3',
      RT_PRIMERS[0][-9:] == expected[2:11], True)
    c('partial graph cannot claim sequencing completeness', protocol_graph().publishable, False)
    c('unresolved molecular endpoint retained',
      protocol_graph().terminals['cdna_circle'].disposition, Disposition.UNRESOLVED)
    c.report()

if __name__ == '__main__':
    selftest()
