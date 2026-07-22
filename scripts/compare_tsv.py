#!/usr/bin/env python3
"""Golden-master TSV comparison for the AsterX test suite.

Compares the TSV output produced by a test run against the committed
reference data under ``AsterX/test/<testname>/``, at a tight
(near-machine-precision) tolerance.  This is stricter than the default
Cactus testsuite comparison (``ABSTOL/RELTOL = 1e-8`` in
``AsterX/test/test.ccl``) and is used to confirm that a refactor -- in
particular the fused single-sweep flux calculation -- reproduces the
current results to floating-point precision rather than merely to 1e-8.

We deliberately reuse the *existing* tests and their *committed* ``.tsv``
data as the golden master: that data was produced by the current code, so
a faithful refactor must reproduce it to TSV precision.  No new tests or
data are introduced.

CarpetX writes TSV with ``setprecision(digits10 + 1)`` = 16 significant
digits for ``double`` (hardcoded, not tunable), so a bitwise-identical
result reproduces the committed strings exactly and a real regression
shows up far above the ~1e-16 noise floor.  We compare at a configurable
relative tolerance (default 1e-13) rather than exact 0 to absorb the
last-digit rounding of the 16-digit text representation and any benign
cross-platform fp differences.

Rows are keyed by their grid locator ``(iteration, patch, level, i, j,
k)`` -- NOT by position -- because the suite runs at different MPI-rank
counts (1 proc and 2 procs) which reorder the rows of a file.

Usage::

    compare_tsv.py REFERENCE_ROOT PRODUCED_ROOT [PRODUCED_ROOT ...]
        [--rtol 1e-13] [--atol 1e-13] [--rel-floor 1e-9]
        [--allow-missing]

``REFERENCE_ROOT`` is the committed test tree (e.g. ``AsterX/test``).
Each ``PRODUCED_ROOT`` is a test-run output tree (e.g.
``<output-dir>/TEST/sim``); every reference ``.tsv`` is matched to the
produced file with the same trailing ``<testname>/<filename>`` path.

Exit status is 0 iff every compared value in every produced tree agrees
with the reference within tolerance.  Standard library only -- no numpy.
"""

import argparse
import math
import os
import sys

# Grid-locator columns: matched exactly, used as the per-row key so that
# differently-ordered outputs (e.g. 1-proc vs 2-proc) still line up.
KEY_COLUMNS = ("iteration", "patch", "level", "i", "j", "k")


def parse_tsv(path):
    """Parse one TSV file.

    Returns (colnames, rows) where ``rows`` maps the grid-locator key
    tuple to a dict {colname: float} of the non-key numeric columns.
    """
    colnames = None
    rows = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                if colnames is None:
                    # Header like: "# 1:iteration\t2:time\t...\t11:Bvecx"
                    header = line.lstrip()[1:]  # drop leading '#'
                    cols = []
                    for tok in header.split("\t"):
                        tok = tok.strip()
                        if not tok:
                            continue
                        # strip the "N:" ordinal prefix
                        cols.append(tok.split(":", 1)[-1])
                    colnames = cols
                continue
            if colnames is None:
                continue
            fields = line.split("\t")
            if len(fields) != len(colnames):
                # tolerate whitespace-separated variants
                fields = line.split()
            if len(fields) != len(colnames):
                raise ValueError(
                    "%s: row has %d fields, header has %d"
                    % (path, len(fields), len(colnames)))
            record = dict(zip(colnames, fields))
            try:
                key = tuple(int(float(record[c])) for c in KEY_COLUMNS
                            if c in record)
            except (KeyError, ValueError):
                raise ValueError("%s: cannot build grid-locator key" % path)
            values = {}
            for c in colnames:
                if c in KEY_COLUMNS:
                    continue
                try:
                    values[c] = float(record[c])
                except ValueError:
                    # non-numeric column (unexpected) -- skip
                    continue
            rows[key] = values
    if colnames is None:
        raise ValueError("%s: no header found" % path)
    return colnames, rows


def index_produced(root):
    """Map every *.tsv under ``root`` to its full path.

    Keyed both by ``<parent>/<basename>`` and by bare ``<basename>`` so a
    reference file can be matched regardless of the produced tree's exact
    nesting.
    """
    by_pair = {}
    by_base = {}
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if not fn.endswith(".tsv"):
                continue
            full = os.path.join(dirpath, fn)
            base = fn
            pair = os.path.join(os.path.basename(dirpath), fn)
            by_pair.setdefault(pair, []).append(full)
            by_base.setdefault(base, []).append(full)
    return by_pair, by_base


def match_produced(ref_rel, by_pair, by_base):
    """Find the produced file corresponding to reference rel-path
    ``<testname>/<filename>``.  Returns a path or None.

    The produced tree's nesting is a flesh-testsuite internal detail
    (e.g. ``<testname>/output-0000/<file>``), so we match progressively:
    (1) exact ``<testname>/<file>`` immediate-parent match; (2) among
    files with the right basename, the unique one whose path contains
    ``<testname>`` as a directory component; (3) a globally-unique
    basename.
    """
    testname = os.path.basename(os.path.dirname(ref_rel))
    filename = os.path.basename(ref_rel)

    pair = os.path.join(testname, filename)
    if pair in by_pair:
        # multiple identical candidates are fine -- pick the first.
        return by_pair[pair][0]

    candidates = by_base.get(filename, [])
    tagged = [testname + os.sep, os.sep + testname + os.sep]
    scoped = [p for p in candidates
              if any(t in p for t in tagged) or ("/" + testname + "/") in p]
    if len(scoped) == 1:
        return scoped[0]
    if len(scoped) > 1:
        return scoped[0]

    if len(candidates) == 1:
        return candidates[0]
    return None


def compare_file(ref_path, prod_path, rtol, atol, rel_floor):
    """Compare one file pair.  Returns (ok, worst_abs, worst_rel, where,
    n_missing_rows).

    The reported relative deviation is normalised by each column's *peak*
    reference magnitude, not the point's own value.  A point whose value is
    small but non-trivial (e.g. velx ~ 1e-4 in a shock) would otherwise turn an
    ordinary ~1e-9 absolute fp difference into a scary ~1e-5 "relative" number;
    dividing by the column peak measures the error against the field's actual
    scale.  Columns whose peak is <= rel_floor (genuinely-zero fields) are
    excluded from the relative statistic and judged by the absolute deviation.
    """
    _, ref_rows = parse_tsv(ref_path)
    _, prod_rows = parse_tsv(prod_path)

    # Per-column peak |reference|, used as the relative-deviation denominator.
    col_peak = {}
    for vals in ref_rows.values():
        for col, a in vals.items():
            if math.isfinite(a):
                av = abs(a)
                if av > col_peak.get(col, 0.0):
                    col_peak[col] = av

    ok = True
    worst_abs = 0.0
    worst_rel = 0.0
    where = ""
    n_missing_rows = 0
    for key, ref_vals in ref_rows.items():
        prod_vals = prod_rows.get(key)
        if prod_vals is None:
            n_missing_rows += 1
            ok = False
            continue
        for col, a in ref_vals.items():
            b = prod_vals.get(col)
            if b is None:
                ok = False
                continue
            ad = abs(a - b)
            if not math.isfinite(a) or not math.isfinite(b):
                # NaN/Inf must match exactly (both NaN or both same Inf)
                if not (math.isnan(a) and math.isnan(b)) and a != b:
                    ok = False
                    if ad != ad or ad > worst_abs:  # NaN or larger
                        worst_abs = float("inf")
                        where = "%s[%s]" % (col, key)
                continue
            if ad > worst_abs:
                worst_abs = ad
            # relative deviation measured against the column's characteristic
            # scale (its peak), so small-magnitude points don't dominate
            peak = col_peak.get(col, 0.0)
            if peak > rel_floor:
                rd = ad / peak
                if rd > worst_rel:
                    worst_rel = rd
                    where = "%s[%s]" % (col, key)
            # numpy.allclose-style pass/fail tolerance
            if ad > atol + rtol * abs(b):
                ok = False
    return ok, worst_abs, worst_rel, where, n_missing_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference", help="committed test tree, e.g. AsterX/test")
    ap.add_argument("produced", nargs="+",
                    help="test-run output tree(s), e.g. <out>/TEST/sim")
    ap.add_argument("--rtol", type=float, default=1e-13,
                    help="relative tolerance (default 1e-13)")
    ap.add_argument("--atol", type=float, default=1e-13,
                    help="absolute tolerance (default 1e-13)")
    ap.add_argument("--rel-floor", type=float, default=1e-9,
                    help="exclude a column from the relative-deviation report "
                         "when its peak |value| <= this (reporting only; "
                         "default 1e-9)")
    ap.add_argument("--allow-missing", action="store_true",
                    help="do not fail when a reference file has no produced "
                         "counterpart (default: fail)")
    args = ap.parse_args()

    ref_root = os.path.abspath(args.reference)
    ref_files = []
    for dirpath, _dirnames, filenames in os.walk(ref_root):
        for fn in filenames:
            if fn.endswith(".tsv"):
                ref_files.append(os.path.join(dirpath, fn))
    ref_files.sort()
    if not ref_files:
        print("ERROR: no reference .tsv files under %s" % ref_root,
              file=sys.stderr)
        return 2

    overall_status = 0
    for produced_root in args.produced:
        produced_root = os.path.abspath(produced_root)
        print("=== comparing produced tree: %s" % produced_root)
        by_pair, by_base = index_produced(produced_root)
        if not by_pair:
            print("ERROR: no produced .tsv files under %s" % produced_root,
                  file=sys.stderr)
            overall_status = 1
            continue

        n_ok = n_fail = n_missing = 0
        tree_worst_abs = tree_worst_rel = 0.0
        tree_where = ""
        for ref_path in ref_files:
            ref_rel = os.path.relpath(ref_path, ref_root)
            prod_path = match_produced(ref_rel, by_pair, by_base)
            if prod_path is None:
                n_missing += 1
                msg = "MISSING produced file for %s" % ref_rel
                if args.allow_missing:
                    print("  [skip] " + msg)
                else:
                    print("  [FAIL] " + msg)
                    overall_status = 1
                continue
            try:
                ok, wa, wr, where, miss_rows = compare_file(
                    ref_path, prod_path, args.rtol, args.atol, args.rel_floor)
            except ValueError as exc:
                print("  [FAIL] %s: %s" % (ref_rel, exc))
                overall_status = 1
                n_fail += 1
                continue
            if wa > tree_worst_abs:
                tree_worst_abs = wa
            if wr > tree_worst_rel:
                tree_worst_rel = wr
                tree_where = "%s :: %s" % (ref_rel, where)
            if ok:
                n_ok += 1
            else:
                n_fail += 1
                overall_status = 1
                extra = (" (%d rows missing in produced)" % miss_rows
                         if miss_rows else "")
                print("  [FAIL] %s: worst |abs|=%.3e worst |rel|=%.3e (%s)%s"
                      % (ref_rel, wa, wr, where, extra))
        print("  files: %d ok, %d fail, %d missing" % (n_ok, n_fail, n_missing))
        print("  tree worst |abs| = %.3e, worst |rel| = %.3e (%s)"
              % (tree_worst_abs, tree_worst_rel, tree_where))

    if overall_status == 0:
        print("GOLDEN-MASTER: PASS (all values within rtol=%g atol=%g)"
              % (args.rtol, args.atol))
    else:
        print("GOLDEN-MASTER: FAIL")
    return overall_status


if __name__ == "__main__":
    sys.exit(main())
