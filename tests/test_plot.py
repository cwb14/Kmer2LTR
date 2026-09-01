import random

import numpy as np
import pytest

from kmer2ltr.k2p import insertion_time
from kmer2ltr.plot import _kde, density_plot, myr_per_unit_k2p, read_k2p
from kmer2ltr.runner import COLUMNS

matplotlib = pytest.importorskip("matplotlib")


def _tsv(tmp_path, rows, name="out.tsv"):
    """A results TSV with the real header and only the fields under test set."""
    p = tmp_path / name
    lines = ["\t".join(COLUMNS)]
    for status, k2p in rows:
        row = dict.fromkeys(COLUMNS, "NA")
        row["status"] = status
        row["k2p"] = k2p
        lines.append("\t".join(row[c] for c in COLUMNS))
    p.write_text("\n".join(lines) + "\n")
    return p


def test_read_k2p_takes_only_passing_rows(tmp_path):
    p = _tsv(tmp_path, [("pass", "0.1"), ("weak_pair", "0.2"), ("pass", "0.3"),
                        ("no_pair", "NA"), ("k2p_undefined", "NA")])
    assert list(read_k2p(p)) == [0.1, 0.3]


def test_read_k2p_skips_unparseable_values(tmp_path):
    """A `k2p_undefined` row carries NA, and a hand-edited file can carry
    anything. Neither should abort a figure over the other 10,000 rows."""
    p = _tsv(tmp_path, [("pass", "0.1"), ("pass", "NA"), ("pass", "nan"),
                        ("pass", "inf"), ("pass", "0.2")])
    assert list(read_k2p(p)) == [0.1, 0.2]


def test_the_top_axis_is_the_same_clock_as_the_k2p_time_column():
    """The figure's age axis and the TSV's `k2p_time` must not disagree -- a
    factor of two here is exactly the classic LTR-dating error."""
    mu, d = 7e-9, 0.0625
    assert d / myr_per_unit_k2p(mu) * 1e6 == pytest.approx(insertion_time(d, mu))


def test_kde_integrates_to_one():
    rng = random.Random(3)
    x = np.array([rng.gauss(0.1, 0.02) for _ in range(2000)])
    grid = np.linspace(-0.1, 0.3, 4001)
    dens = _kde(x, grid)
    # Riemann sum rather than np.trapezoid: np.trapz is gone in numpy 2 and
    # np.trapezoid absent before it, and the grid here is fine enough that the
    # difference is far below the tolerance.
    assert float(dens.sum() * (grid[1] - grid[0])) == pytest.approx(1.0, abs=0.01)


def test_kde_declines_on_input_with_no_spread():
    """A bandwidth of zero is a division by zero, and a one-element set has no
    density to estimate. Both must return None rather than a NaN curve."""
    grid = np.linspace(0, 1, 10)
    assert _kde(np.array([0.1]), grid) is None
    assert _kde(np.array([0.1] * 50), grid) is None


def test_density_plot_writes_a_pdf_and_reports_the_element_count(tmp_path):
    rng = random.Random(5)
    rows = [("pass", f"{rng.gauss(0.08, 0.02):.6f}") for _ in range(300)]
    rows += [("weak_pair", "0.4"), ("no_pair", "NA")]
    tsv = _tsv(tmp_path, rows)
    out = tmp_path / "d.pdf"
    assert density_plot(tsv, out, mutation_rate=7e-9) == 300
    assert out.read_bytes()[:5] == b"%PDF-"


def test_density_plot_works_without_a_mutation_rate(tmp_path):
    tsv = _tsv(tmp_path, [("pass", f"{0.05 + i * 0.001:.4f}") for i in range(60)])
    out = tmp_path / "d.pdf"
    assert density_plot(tsv, out) == 60
    assert out.exists()


def test_density_plot_writes_nothing_when_no_element_passed(tmp_path, capsys):
    tsv = _tsv(tmp_path, [("no_pair", "NA"), ("weak_pair", "0.2")])
    out = tmp_path / "d.pdf"
    assert density_plot(tsv, out) == 0
    assert not out.exists()
    assert "no passing element" in capsys.readouterr().err


def test_a_saturated_tail_does_not_squash_the_axis(tmp_path):
    """A handful of near-saturated elements would otherwise compress every real
    burst into the first few percent of the plot. They are counted in the
    annotation, never dropped from the data."""
    rows = [("pass", "0.05")] * 400 + [("pass", "0.9")] * 3
    tsv = _tsv(tmp_path, rows)
    out = tmp_path / "d.pdf"
    assert density_plot(tsv, out) == 403
    x = read_k2p(tsv)
    assert float(np.percentile(x, 99.5)) < 0.9


def test_reflected_kde_still_integrates_to_one_over_the_non_negative_half():
    """K2P cannot be negative. A plain kernel density puts mass below zero and
    understates the density just above it -- exactly where the recent-insertion
    peak sits. Reflection moves that mass back without changing the total."""
    rng = random.Random(11)
    x = np.abs(np.array([rng.gauss(0.0, 0.03) for _ in range(3000)]))
    grid = np.linspace(0.0, 0.4, 4001)
    dens = _kde(x, grid, reflect=True)
    dx = grid[1] - grid[0]
    assert float(dens.sum() * dx) == pytest.approx(1.0, abs=0.02)
    plain = _kde(x, grid)
    assert dens[0] > 1.8 * plain[0], "reflection should roughly double density at the boundary"
    assert dens[-1] == pytest.approx(plain[-1], abs=1e-6), "far from the boundary, unchanged"
