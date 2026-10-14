from dsb_capture import SOURCE_OLIGO_EXCEPTIONS, get  # noqa: F401

_SPEC = get("dsbcapture")
globals().update(vars(_SPEC))


def protocol_graph():
    """Exposed as a module-level def so tools/selftest_source_oligos.py finds this graph.

    The discovery scan reads the AST, so the namespace lambda copied in above is invisible
    to it and the published oligo states would never be compared against source text.
    """
    return _SPEC.protocol_graph()
