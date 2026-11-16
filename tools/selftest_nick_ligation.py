#!/usr/bin/env python3
"""Neutral end-chemistry checks: no protocol-dependent primer or adapter recipe."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'lib')]
from animation.simulator import Simulator,Polymerase,Ligase,complement


def fixture(primary_five='OH'):
    s=Simulator();src='synthetic asymmetric duplex fixture'
    s.supply_nt('DNA','ACGT',source=src)
    a=s.add('primary','ACGTGCTAGTCA',five=primary_five,source=src)
    primer=s.add('random',complement('ACGTGCTAGTCA'[-6:][::-1]),five='P',source=src)
    s.anneal('random','primary',binding=6,source=src)
    s.extend('random',Polymerase('copy',('DNA',),'DNA'),source=src)
    b=s.chain('random');before=s.snapshot()
    tail=Polymerase('dA only',('DNA',),'DNA',terminal_bases=('A',))
    s.terminal_add('primary','A',1,tail,source=src)
    s.terminal_add('random','A',1,tail,source=src)
    a,b=s.chain('primary'),s.chain('random')
    for side in ('left','right'):
        s.add(side+'_top','CCGGT',source=src)
        s.add(side+'_bottom','CCGG',five='P',source=src)
        s.anneal(side+'_bottom',side+'_top',binding=4,template_site=3,source=src)
    s.bind_spans(s.chain('left_top')[-1:],b[-1:],source=src)
    s.bind_spans(a[-1:],s.chain('right_top')[-1:],source=src)
    ligase=Ligase('test ligase',('DNA',),('DNA',),('P',),('ATP',),('DNA',))
    return s,a,b,primer,before,ligase


class NickLigationTests(unittest.TestCase):
    def test_copy_and_da_tail_preserve_inherited_five_prime_groups(self):
        s,a,b,primer,before,_=fixture()
        self.assertEqual(before.five[a[0]],'OH')
        self.assertEqual(before.five[b[0]],'P')
        self.assertEqual(s.five[a[0]],'OH');self.assertEqual(s.five[b[0]],'P')
        self.assertEqual(b[:len(primer)],primer)

    def test_three_of_four_nicks_seal_without_repair(self):
        s,a,b,_,_,ligase=fixture()
        blocked=(s.tail('left_top'),a[0])
        sites=s.ligation_sites(ligase,['primary','random'],cofactors=('ATP',))
        self.assertEqual(len(sites),3);self.assertNotIn(blocked,sites)
        operations=s.ligate_available(ligase,['primary','random'],cofactors=('ATP',),source='test')
        self.assertEqual(len(operations),3)
        self.assertFalse(s.ligation_sites(ligase,['primary','random'],cofactors=('ATP',)))
        self.assertEqual(s.five[a[0]],'OH');self.assertNotIn(a[0],s.previous)
        self.assertEqual(len(s.chain('left_top')),5)
        self.assertGreater(len(s.chain('random')),len(b))

    def test_oh_nick_cannot_be_silently_repaired(self):
        s,a,_,_,_,ligase=fixture();before=s.snapshot()
        with self.assertRaises(ValueError):
            s.ligate(s.tail('left_top'),a[0],ligase,mode='nick',cofactors=('ATP',),source='test')
        self.assertEqual(s.snapshot(),before)

    def test_both_phosphorylated_has_four_sealable_nicks(self):
        s,_,_,_,_,ligase=fixture('P')
        self.assertEqual(len(s.ligate_available(ligase,['primary','random'],cofactors=('ATP',),source='test')),4)

    def test_missing_cofactor_does_not_partially_react(self):
        s,_,_,_,_,ligase=fixture();before=s.snapshot()
        with self.assertRaises(ValueError):s.ligate_available(ligase,['primary'],cofactors=(),source='test')
        self.assertEqual(s.snapshot(),before)

    def test_da_only_constructor_retains_asymmetric_phosphates(self):
        from endprep import dA_tail_only,dA_tailed_scene,ThreePrimeOverhangAdapter,ligate_both_ends,cohesive_ligation_geometry
        from chemdraw import Construct,Segment
        product=dA_tail_only(Construct([Segment('insert','ACGT')]),
                             top_five_prime_phosphate=False,bottom_five_prime_phosphate=True)
        scene=dA_tailed_scene(product)
        self.assertEqual(scene.strands['top'].mod5,'OH')
        self.assertEqual(scene.strands['bottom'].mod5,'p')
        adapter=ThreePrimeOverhangAdapter('toy','CCGGT','CCGG')
        self.assertEqual(cohesive_ligation_geometry(product,adapter,statement='test').availability,'known')
        with self.assertRaises(ValueError):ligate_both_ends(product,adapter)

    def test_da_only_does_not_accept_undeclared_end_states(self):
        from endprep import dA_tail_only
        from chemdraw import Construct,Segment
        with self.assertRaises(ValueError):
            dA_tail_only(Construct([Segment('insert','ACGT')]),
                         top_five_prime_phosphate=None,bottom_five_prime_phosphate=True)


if __name__=='__main__':unittest.main()
