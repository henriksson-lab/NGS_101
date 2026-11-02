"""Causal checks and identity/continuity audit for three selected mechanisms."""
import math
import unittest
from dataclasses import replace
from .capture_ligation import crispr_mip,smallseq,hic,smallseq_simulation,MIP,SMALL,HIC
from .other_simulations import mip as simulate_mip,hic as simulate_hic
from .protocols import load
from .evaluate import evaluate,presentation
from .simulator import Ligase


class CaptureLigationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.movies=[crispr_mip(),smallseq(),hic()]

    def test_generated_frames_and_conserved_transition_boundaries(self):
        for movie in self.movies:
            movie.validate()
            for shot in movie.shots:
                for t in [0,shot.timing()[2]]:
                    states,_=evaluate(shot,t);frame=next(iter(states.values())).molecular
                    drawn={id for st in states.values()for j,id in enumerate(st.strand.nucleotide_ids)if st.fraction(j)}
                    self.assertEqual(drawn,set(frame.atoms)-frame.released,(movie.title,shot.id,t))
            for previous,shot in zip(movie.shots,movie.shots[1:]):
                self.assertTrue(shot.transition)
                a=evaluate(previous,previous.timing()[2])[0]
                b=evaluate(shot,0)[0]
                self.assertEqual(next(iter(a.values())).molecular,next(iter(b.values())).molecular,(movie.title,shot.id))
                actual,_,camera=presentation(movie,shot,0)
                self.assertEqual(camera,previous.camera)
                locations={id:(st,j)for st in a.values()for j,id in enumerate(st.strand.nucleotide_ids)if st.fraction(j)}
                for st in actual.values():
                    for j,id in enumerate(st.strand.nucleotide_ids):
                        if st.fraction(j):
                            prior,index=locations[id]
                            self.assertLess(math.dist(st.point(j+.5),prior.segment_point(index+.5,index)),1e-7)

    def test_duplex_contacts_land_on_their_partners_after_each_shot(self):
        for movie in self.movies:
            for shot in movie.shots:
                states,_=evaluate(shot,shot.timing()[2]);frame=next(iter(states.values())).molecular
                loc={id:st.segment_point(i+.5,i)for st in states.values()for i,id in enumerate(st.strand.nucleotide_ids)if st.fraction(i)}
                for a,b in frame.duplex:
                    if a in loc and b in loc:self.assertAlmostEqual(math.dist(loc[a],loc[b]),24,places=6,msg=(movie.title,shot.id,a,b))

    def test_mip_copy_is_substrate_derived_and_only_sealing_closes(self):
        r=simulate_mip(load(MIP,'crisprmip.py'));self.assertEqual(len(r['extend'].added),112)
        self.assertEqual(r['extend'].details['stop'],'occupied_template')
        p=r['probe']
        self.assertNotIn(p[0],r['extend'].after.previous)
        self.assertEqual(r['seal'].after.next[r['filled'][-1]],p[0])
        self.assertEqual(r['seal'].before.atoms,r['seal'].after.atoms)
        self.assertEqual(r['seal'].before.paired,r['seal'].after.paired)
        self.assertFalse(r['circle_exoi'].details['removed']);self.assertFalse(r['circle_exoiii'].details['removed'])
        self.assertEqual(len(r['unused_digest'].details['removed']),len(r['unused'])-2)
        self.assertTrue(set(r['unused_digest'].details['removed'])<=r['closed'].released)

    def test_smallseq_uses_real_adenylated_and_blocked_source_orders(self):
        r=smallseq_simulation(load(SMALL,'smallseq.py'));rna,ra3,ra5=r['rna'],r['ra3'],r['ra5']
        self.assertEqual(r['initial'].five[ra3[0]],'App');self.assertEqual(r['initial'].three[ra3[-1]],'ddC')
        self.assertEqual(r['initial'].five[ra5[0]],'NH2')
        self.assertEqual(r['lig3'].after.next[rna[-1]],ra3[0]);self.assertEqual(r['lig5'].after.next[ra5[-1]],rna[0])
        self.assertEqual(r['lig5'].details['cofactors'],['ATP']);self.assertEqual(r['lig3'].details['cofactors'],[])
        self.assertEqual(len(r['digest'].details['removed']),len(r['free'])-1)
        self.assertTrue(all(id not in r['digest'].after.released for id in r['rtp']))
        self.assertEqual(sum('umi.H8'in r['initial'].atoms[id].features for id in ra5),8)

    def test_smallseq_second_ligation_requires_atp(self):
        r=smallseq_simulation(load(SMALL,'smallseq.py'));s=r['s']
        # A fresh pair of eligible RNA substrates cannot be joined without ATP.
        a=s.add('controla','ACGU',sugar='RNA',source='synthetic causal control')
        b=s.add('controlb','UGCA',sugar='RNA',five='P',source='synthetic causal control')
        with self.assertRaisesRegex(ValueError,'cofactor'):
            s.ligate(a[-1],b[0],Ligase('T4 RNA ligase1',('RNA',),('RNA',),('P',),('ATP',),()),mode='single-strand',cofactors=(),source='synthetic causal control')

    def test_hic_two_strand_sealing_and_biotin_inheritance(self):
        r=simulate_hic(load(HIC,'insitu_hic.py'))
        self.assertEqual(''.join(r['seals'][-1].after.atoms[id].base for id in r['joined'])[14:22],'GATCGATC')
        before=r['seal_before'];after=r['seals'][-1].after
        self.assertEqual(before.atoms,after.atoms);self.assertEqual(len(set(after.next.items())-set(before.next.items())),2)
        added=list(r['fill_a'].added)+list(r['fill_b'].added)
        self.assertEqual(len(added),8)
        self.assertEqual(sum(after.atoms[id].modification=='biotin'for id in added),2)
        self.assertEqual(len(after.atoms),128+8)
        for original in r['originals'].values():
            for ids in original:self.assertTrue(set(ids)<=set(after.atoms))

    def test_native_preview_has_no_text_overflows(self):
        try:from .render import Renderer,Style;r=Renderer(Style(640,360))
        except RuntimeError:self.skipTest('optional pycairo unavailable')
        for movie in self.movies:
            offset=0
            for shot in movie.shots:
                for t in [0,shot.transition/2,shot.establish+shot.action*.75,shot.timing()[2]-.01]:
                    r.image(movie,(offset+t)/1.5)
                    self.assertFalse(r.overflows,(movie.title,shot.id,r.overflows))
                offset+=shot.timing()[2]


if __name__=='__main__':unittest.main(verbosity=2)
