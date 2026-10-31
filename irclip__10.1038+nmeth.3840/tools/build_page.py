"""Source-qualified irCLIP page; legacy build entry point retained."""
from pathlib import Path
from irclip import render_page
OUT=Path(__file__).resolve().parents[1]/'irclip_partial.html'
def build():
    OUT.write_text(render_page());print(OUT)
if __name__=='__main__':build()
