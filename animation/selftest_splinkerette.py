"""Molecular state, smooth transitions and native layout of the splinkerette film."""
import math
import unittest
from .splinkerette import movie
from .evaluate import evaluate,presentation
from .molecular_frame import molecular_frame
from .protocols import load


class SplinkeretteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.movie=movie()

    def test_every_transition_preserves_the_full_molecular_state(self):
        for previous,shot in zip(self.movie.shots,self.movie.shots[1:]):
            self.assertEqual(molecular_frame(previous,previous.timing()[2]),molecular_frame(shot,0),shot.id)
            states,_,_=presentation(self.movie,shot,0)
            old,_=evaluate(previous,previous.timing()[2])
            places={id:(st,j)for st in old.values()for j,id in enumerate(st.strand.nucleotide_ids)if st.fraction(j)}
            for st in states.values():
                for j,id in enumerate(st.strand.nucleotide_ids):
                    if st.fraction(j):
                        if id not in places:
                            self.assertFalse(0<=st.point(j+.5)[0]<=1600 and 0<=st.point(j+.5)[1]<=900,shot.id)
                            continue
                        oldst,k=places[id]
                        self.assertLess(math.dist(st.point(j+.5),oldst.segment_point(k+.5,k)),1e-7,shot.id)

    def test_all_final_pairs_have_close_aligned_coordinates(self):
        for shot in self.movie.shots:
            states,_=evaluate(shot,shot.timing()[2]);frame=molecular_frame(shot,shot.timing()[2])
            loc={id:st.segment_point(j+.5,j)for st in states.values()for j,id in enumerate(st.strand.nucleotide_ids)if st.fraction(j)}
            for a,b in frame.duplex:
                if a in loc and b in loc:self.assertAlmostEqual(math.dist(loc[a],loc[b]),24,places=6,msg=shot.id)

    def test_intact_stock_primers_arrive_before_their_extension(self):
        for shot in self.movie.shots:
            frame=molecular_frame(shot,0)
            for event in shot.events:
                if event.kind!='introduce':continue
                new={n['id']for n in event.molecular_operation['added']}
                self.assertFalse(new & frame.atoms.keys(),shot.id)
                self.assertTrue(new <= molecular_frame(shot,shot.establish+.1).atoms.keys(),shot.id)

    def test_both_nested_complexes_receive_late_polymerases(self):
        for key in ['nested1','nested2']:
            shot=next(q for q in self.movie.shots if q.id==key)
            self.assertEqual(len(evaluate(shot,shot.establish+shot.action*.4)[1]),0)
            self.assertEqual(len(evaluate(shot,shot.establish+shot.action*.7)[1]),2)

    def test_one_restriction_complex_cuts_both_backbones(self):
        shot=self.movie.shots[0]
        self.assertEqual(len([e for e in shot.events if e.kind=='cleave']),2)
        self.assertEqual(len(evaluate(shot,shot.establish+shot.action*.4)[1]),1)

    def test_unphosphorylated_opposite_junction_never_becomes_a_bond(self):
        m=load(self.movie.protocol,'splinkerette.py');r=m.simulate()
        forbidden=(r['genome_top'][-1],r['adapter_top'][0])
        for shot in self.movie.shots:
            frame=molecular_frame(shot,shot.timing()[2])
            self.assertNotIn(forbidden,frame.backbone)

    def test_shared_engine_remains_the_compatibility_engine(self):
        from .simulator import Simulator as old
        from lib.nucleotide_simulator import Simulator as shared
        self.assertIs(old,shared)

    def test_native_frames_have_no_text_overflow(self):
        try:
            from .render import Renderer,Style
            renderer=Renderer(Style(960,540))
        except RuntimeError:self.skipTest('optional Cairo not installed')
        offset=0
        for shot in self.movie.shots:
            for t in [0,.9,shot.establish+shot.action*.75,shot.timing()[2]-.01]:
                renderer.image(self.movie,offset+t,speed=1)
                self.assertFalse(renderer.overflows,(shot.id,renderer.overflows))
            offset+=shot.timing()[2]
