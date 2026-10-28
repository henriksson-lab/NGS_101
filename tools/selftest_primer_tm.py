"""Offline Tm parameter export and page generation."""
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
sys.path[:0]=[str(Path(__file__).resolve().parents[1]),str(Path(__file__).resolve().parents[1]/'lib')]
import build_tm
from chemdraw import tm_model_parameters


def run(check):
    check.section('primer Tm page and shared thermodynamic data')
    model=build_tm.model_data()
    check('browser receives the existing canonical thermodynamic data',model['thermodynamics'],tm_model_parameters())
    check('preset identifies the ReadyMix product and the original revision',
          'ReadyMix PCR Kit' in model['preset']['name'] and 'KR0370 v16.25' in model['preset']['revision'] and model['preset']['source'].startswith('https://elabdoc-prod.roche.com/'),True)
    with TemporaryDirectory() as directory:
        output=Path(directory);build_tm.build(output);page=(output/'tm.html').read_text()
        payload=page.split('<script type="application/json" id="tm-model">')[1].split('</script>')[0]
        check('page contains its offline model and source-qualified preset',json.loads(payload),model)
        check('monovalent assumption is visible beside the inputs',
              'Monovalent ions (mM, assumed)' in page and 'not Roche’s formulation' in page,True)
        check('page uses its own shared navigation and reusable model script',
              'href="tm.html"' in page and 'src="assets/primer_tm.js"' in page,True)
        check('page separates experimental annealing guidance from calculated Tm',
              'no fixed Tm-to-Ta offset' in page and '60–72°C' in page,True)


if __name__=='__main__':
    from checks import Check
    c=Check();run(c);c.report()
