import pytest
from ltrk2p.cigar import aligned_pair_from_cigartuples, extended_cigar, cs_string

def test_extended_cigar_uses_eq_and_X_not_M():
    a = "ACGTACGT"
    b = "ACGTACCT"
    assert extended_cigar(a, b) == "6=1X1="
    assert "M" not in extended_cigar(a, b)

def test_extended_cigar_insertion_and_deletion_orientation():
    #        query has extra A (I) ; query missing G (D)
    a = "ACGTA-CGT"
    b = "ACGT-GCGT"
    assert extended_cigar(a, b) == "4=1I1D3="

def test_cs_string_matches_minimap2_short_form():
    a = "ACGTACGT"
    b = "ACGTACCT"
    # substitution encoded as *<ref><query>, lowercase
    assert cs_string(a, b) == ":6*cg:1"

def test_cs_string_indels():
    a = "ACGTA-CGT"
    b = "ACGT-GCGT"
    assert cs_string(a, b) == ":4+a-g:3"

def test_pywfa_cigartuples_expand_correctly():
    # pins the op-code mapping: 0 = equal (pywfa emits M for true matches), 8 = X
    tuples = [(0, 4), (8, 1), (0, 3)]
    a, b = aligned_pair_from_cigartuples(tuples, "ACGTACGT", "ACGTCCGT")
    assert (a, b) == ("ACGTACGT", "ACGTCCGT")
    assert extended_cigar(a, b) == "4=1X3="

def test_roundtrip_reconstructs_both_sequences():
    a = "ACGTA-CGTTGCA"
    b = "ACGT-GCGTTGCA"
    cig = extended_cigar(a, b)
    q = a.replace("-", "")
    r = b.replace("-", "")
    # walk the cigar and rebuild
    import re
    qi = ri = 0
    qa, rb_ = [], []
    for n, op in re.findall(r"(\d+)([=XID])", cig):
        n = int(n)
        if op in "=X":
            qa.append(q[qi:qi+n]); rb_.append(r[ri:ri+n]); qi += n; ri += n
        elif op == "I":
            qa.append(q[qi:qi+n]); rb_.append("-"*n); qi += n
        else:
            qa.append("-"*n); rb_.append(r[ri:ri+n]); ri += n
    assert "".join(qa) == a and "".join(rb_) == b

def test_rejects_unequal_lengths():
    with pytest.raises(ValueError):
        extended_cigar("ACGT", "ACG")

def test_pywfa_real_orientation_is_query_first():
    """Guard: pywfa's pattern is the query, text is the ref. If a pywfa
    upgrade flips this, I/D orientation silently inverts -- catch it here."""
    from pywfa import WavefrontAligner
    from ltrk2p.cigar import aligned_pair_from_cigartuples, extended_cigar
    query = "ACGTAACGT"   # one extra A relative to ref
    ref = "ACGTACGT"
    al = WavefrontAligner(query, scope="full", span="end-to-end")
    al(ref)
    a, b = aligned_pair_from_cigartuples(al.cigartuples, query, ref)
    assert a.replace("-", "") == query
    assert b.replace("-", "") == ref
    assert "I" in extended_cigar(a, b)   # extra base is in the QUERY -> I
