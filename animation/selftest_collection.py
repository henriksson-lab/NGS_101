"""Shared acceptance checks across ten additional mechanism movies."""
import unittest,math
from .examples import EXAMPLES
from .evaluate import evaluate,presentation
from .batch import FOCUS
NAMES=('crispr-mip','smallseq','hic','celseq','dropseq','scrrbs','gibson','atac','rnaseh','golden-gate')

class CollectionTests(unittest.TestCase):
    def test_recipe_serialization_is_repeatable_across_processes(self):
        import json,subprocess,sys
        code="import hashlib,json;from animation.examples import EXAMPLES;from animation.selftest_collection import NAMES;print(json.dumps({n:hashlib.sha256(json.dumps(EXAMPLES[n]().to_dict(),sort_keys=True).encode()).hexdigest()for n in NAMES},sort_keys=True))"
        a=subprocess.check_output([sys.executable,'-c',code])
        b=subprocess.check_output([sys.executable,'-c',code])
        self.assertEqual(json.loads(a),json.loads(b))

    def test_each_example_has_causal_chemistry_and_a_real_focus(self):
        for name in NAMES:
            m=EXAMPLES[name]();m.validate()
            self.assertIn(FOCUS[name],{s.id for s in m.shots},name)
            self.assertTrue(m.provenance.get('molecular_simulation'),name)
            for shot in m.shots:
                self.assertTrue(shot.chemistry_generated,(name,shot.id))
                self.assertTrue(shot.source)
                self.assertTrue(shot.scope)
    def test_every_active_atom_has_one_track_during_reactions(self):
        for name in NAMES:
            m=EXAMPLES[name]()
            for shot in m.shots:
                ids=[id for s in shot.strands for id in s.nucleotide_ids]
                self.assertEqual(len(ids),len(set(ids)),(name,shot.id,'duplicate physical atom'))
                for f in [0,.2,.5,.8,1]:
                    states,_=evaluate(shot,shot.establish+shot.action*f)
                    frame=next(iter(states.values())).molecular
                    self.assertFalse((set(frame.atoms)-frame.released)-set(ids),(name,shot.id,f))
    def test_all_incoming_views_reframe_identical_chemistry_continuously(self):
        for name in NAMES:
            m=EXAMPLES[name]()
            for previous,shot in zip(m.shots,m.shots[1:]):
                self.assertGreater(shot.transition,0,(name,shot.id))
                old,_=evaluate(previous,previous.timing()[2]);new,_,camera=presentation(m,shot,0)
                self.assertEqual(camera,previous.camera)
                locations={id:(st,j)for st in old.values()for j,id in enumerate(st.strand.nucleotide_ids)if st.fraction(j)}
                for st in new.values():
                    for j,id in enumerate(st.strand.nucleotide_ids):
                        if st.fraction(j):
                            before,k=locations[id]
                            self.assertLess(math.dist(st.point(j+.5),before.segment_point(k+.5,k)),1e-7,(name,shot.id,id))
                # Every intermediate and endpoint keeps the exact incoming frame.
                for fraction in [.25,.5,.75,1]:
                    view,_,_=presentation(m,shot,shot.transition*fraction)
                    self.assertEqual(next(iter(view.values())).molecular,next(iter(new.values())).molecular,(name,shot.id))

if __name__=='__main__':unittest.main()
