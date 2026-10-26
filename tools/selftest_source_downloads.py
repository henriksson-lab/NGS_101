"""Source-request links retain identity without entering public navigation/search."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_source_downloads import parse_requests, render_page, ROOT, URL

FIXTURE = '''Outstanding documents

1. Example manual,
   Revision E1
   Official PDF:
   https://example.org/manual.pdf?file=a&version=E1
   DOI: https://doi.org/10.1016/0378-1119(89)90359-4
   Status: <unavailable> & needs review.
   Directory: example__10.1+example
'''


def run(check):
    from html.parser import HTMLParser
    class Links(HTMLParser):
        def __init__(self):
            super().__init__(); self.hrefs = []
        def handle_starttag(self, tag, attrs):
            if tag == 'a': self.hrefs.append(dict(attrs).get('href'))
    check.section('unlisted source downloads preserve requested document identities')
    request, = parse_requests(FIXTURE)
    check('wrapped manual revision stays visible in its title', request.title, 'Example manual, Revision E1')
    page = render_page(FIXTURE); links = Links(); links.feed(page)
    check('parenthesized DOI and query parameters survive as exact clickable URLs',
          links.hrefs, [url for url, _ in request.links])
    check('request text is escaped rather than interpreted as HTML',
          '&lt;unavailable&gt; &amp; needs review.' in page, True)
    check('page excludes source-root preamble and requests search-engine exclusion',
          ('/husky/' in page, 'content="noindex,nofollow"' in page), (False, True))
    check('fulfilled request removal removes its links', 'example.org' in render_page(''), False)
    check.raises('an unidentified source request cannot disappear silently',
                 lambda: parse_requests('1. Broken request\n   https://example.org/file.pdf\n'))
    check.raises('unsafe destination path is rejected',
                 lambda: parse_requests(FIXTURE.replace('example__10.1+example', '../outside')))
    actual = (ROOT / 'geturl.txt').read_text()
    parser = Links(); parser.feed(render_page(actual))
    check('all curated request URLs remain available on the helper page',
          set(parser.hrefs), set(URL.findall(actual)))
    import build_index, build_site
    index, _ = build_index.build(run_checks=False)
    check('helper page stays out of public index navigation and search', 'source-downloads.html' in index, False)
    check('raw request ledger is not a publishable asset', build_site.published(ROOT / 'geturl.txt'), False)


if __name__ == '__main__':
    from checks import Check
    check = Check(); run(check); check.report()
