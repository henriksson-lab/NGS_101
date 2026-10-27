"""Offline restriction catalogue, provenance and site-build checks."""
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
sys.path[:0]=[str(Path(__file__).resolve().parents[1]),str(Path(__file__).resolve().parents[1]/'lib')]
import build_restriction


def run(check):
    check.section('licensed restriction snapshot and browser build')
    data=build_restriction.catalogue();enzymes={e['name']:e for e in data['enzymes']}
    check('reference snapshot has provenance, a pinned revision and the original permission notice',
          len(enzymes)>1000 and data['revision'] in data['source'] and len(data['source_sha256'])==64
          and 'Permission to use, copy, modify, and distribute' in data['license_text'],True)
    check('EcoRI coordinates preserve the declared opposite-strand offset',enzymes['EcoRI']['cuts'],[[1,5]])
    check('PstI is represented with a three-prime overhang',enzymes['PstI']['cuts'],[[5,1]])
    check('BsaI type IIS cuts are outside its recognition sequence',enzymes['BsaI']['cuts'],[[7,11]])
    check('two-cut enzyme retains both cut pairs',len(enzymes['BcgI']['cuts']),2)
    check('HpaII/MspI have separate source-qualified CpG annotations',
          (enzymes['HpaII']['methylation']['cpg'],enzymes['MspI']['methylation']['cpg']),('blocked','unaffected'))
    check('DpnI requirement is not confused with methylation resistance',enzymes['DpnI']['methylation']['dam'],'requires')
    check('unannotated variants do not inherit a different reagent annotation',enzymes['BsaI']['methylation']['cpg'],'unknown')
    check('all assigned sensitivity annotations identify primary sources',
          all(e.get('methylation_sources') for e in enzymes.values() if any(v!='unknown' for v in e['methylation'].values())),True)
    with TemporaryDirectory() as tmp:
        output=Path(tmp);build_restriction.build(output)
        page=(output/'restriction.html').read_text();bundle=json.loads((output/'assets/restriction_catalogue.json').read_text())
        check('normal offline builder publishes the complete licensed snapshot',bundle,data)
        check('page provides topology, counts, contexts and independent sensitivity filters',
              all(f'id="restriction-{key}"' in page for key in ('topology','count-mode','count','context','cpg','dam','dcm','custom','export','export-sites')),True)
        check('restriction tab has shared navigation', 'href="restriction.html"' in page,True)
        check('builder does not install or ship the third-party Python package',list(output.rglob('*.py')),[])


if __name__=='__main__':
    from checks import Check
    c=Check();run(c);c.report()
