"""Check lariat bond ports and geometry without a splicing/RT prediction."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from dataclasses import replace
import xml.etree.ElementTree as ET
from rna_lariat import RNALariatTopology, RNAOxygenBond
from protocol_graph import ProtocolGraph, MaterialKind
from chemistry_validation import check_graph, Status


def run(check):
    check.section('RNA lariat backbone and topology')
    rna = RNALariatTopology('GUGCGAACUUAGU', 6)
    check('loop and tail retain supplied nucleotide order',
          (rna.loop_sequence, rna.tail_sequence), ('GUGCGAA', 'CUUAGU'))
    check('branch attaches first phosphate at a distinct 2prime port',
          rna.bonds[-1], RNAOxygenBond(0, 6, 2))
    check('ordinary backbone continues from branch into tail',
          RNAOxygenBond(7, 6, 3) in rna.bonds, True)
    check('RNA five-prime splice site is covalently occupied', rna.five_prime_is_free, False)
    check('tail retains its terminal three-prime hydroxyl', rna.free_three_prime_index, 12)
    check('loop plus tail does not duplicate a nucleotide',
          rna.loop_sequence + rna.tail_sequence, rna.sequence)
    check('each supplied phosphate has one backbone attachment',
          sorted(b.phosphate_nt for b in rna.bonds), list(range(13)))
    check('unknown bases do not invent a branch base preference',
          RNALariatTopology('NGUC', 2).sequence[2], 'U')
    for label, sequence, branch in [('empty', '', 1), ('DNA thymine', 'GTAC', 2),
                                   ('invalid base', 'GUZC', 2), ('first nucleotide', 'GUAC', 0),
                                   ('out of range', 'GUAC', 4), ('negative', 'GUAC', -1),
                                   ('boolean', 'GUAC', True)]:
        try:
            RNALariatTopology(sequence, branch)
        except ValueError:
            rejected = True
        else:
            rejected = False
        check('reject invalid lariat ' + label, rejected, True)
    final_branch = RNALariatTopology('GUGA', 3)
    check('terminal branch has empty downstream tail, not missing three-prime OH',
          (final_branch.tail_sequence, final_branch.free_three_prime_index), ('', 3))
    for name, obj in [('tailed', rna), ('no tail', final_branch),
                      ('long loop', RNALariatTopology('G' + 'U'*45 + 'ACU', 46))]:
        layout = obj.layout()  # shared planar validation runs in construction
        check(name + ' layout has one coordinate per nucleotide',
              len(layout.strands[obj.name].points), len(obj.sequence))
    svg = rna.svg(caption='Invented RNA reference')
    root = ET.fromstring(svg)
    ns = {'s': 'http://www.w3.org/2000/svg'}
    check('SVG retains nucleotide identities, including uracil',
          ''.join(n.text for n in root.findall('.//s:text[@class="nt"]', ns)), rna.sequence)
    check('SVG marks the real free end without a false free-five-prime label',
          [n.attrib['data-end'] for n in root.findall('.//*[@data-end]')], ['3-prime-free'])
    check('SVG bond metadata distinguishes 2prime oxygen',
          [n.attrib['data-oxygen-position'] for n in root.findall('.//*[@data-oxygen-position]')], ['2'])
    graph = ProtocolGraph('declared lariat reference')
    material = graph.input('rna', 'supplied RNA', MaterialKind.MOLECULAR, chemical_state=rna)
    check('checker establishes only supplied connectivity', check_graph(graph)[0].status, Status.ESTABLISHED)
    class MissingBranch(RNALariatTopology):
        @property
        def bonds(self):
            return super().bonds[:-1]
    damaged = object.__new__(MissingBranch)
    object.__setattr__(damaged, 'sequence', rna.sequence)
    object.__setattr__(damaged, 'branch_index', rna.branch_index)
    object.__setattr__(damaged, 'name', rna.name)
    graph.materials['rna'] = replace(material, chemical_state=damaged)
    check('checker rejects a lost branch even without redrawing',
          check_graph(graph)[0].status, Status.CONTRADICTORY)
    class WrongPort(RNALariatTopology):
        @property
        def bonds(self):
            return super().bonds[:-1] + (RNAOxygenBond(0, self.branch_index, 3),)
    damaged = object.__new__(WrongPort)
    object.__setattr__(damaged, 'sequence', rna.sequence)
    object.__setattr__(damaged, 'branch_index', rna.branch_index)
    object.__setattr__(damaged, 'name', rna.name)
    graph.materials['rna'] = replace(material, chemical_state=damaged)
    check('checker rejects a 3prime bond that consumes the tail attachment port',
          check_graph(graph)[0].status, Status.CONTRADICTORY)


if __name__ == '__main__':
    from checks import Check
    check = Check()
    run(check)
    check.report()
