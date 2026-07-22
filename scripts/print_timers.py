#!/usr/bin/env python3
"""Pretty-print AsterX per-routine timings from TimerReport output.

Scans the given test-run directories for TimerReport ``AllTimersReadable.txt``
files and, for each test, prints the AsterX per-function timers from the final
(largest-iteration) report block, sorted by wall time -- so ``AsterX_Fluxes``
and the other hot routines are easy to read in the CI log, instead of the raw
multi-column timer dump.

Usage:
    print_timers.py DIR [DIR ...]      # e.g. the ONEPROC/TWOPROC output roots

Standard library only.
"""

import os
import re
import sys

# Per-function timers are named like "[0044] AsterX: AsterX_Fluxes in ODESolvers_RHS"
FUNC_RE = re.compile(r"\[\d+\]\s*AsterX:\s*(\S+)\s+in\s+(\S+)")
TOP_N = 8


def parse_rows(path):
    """Yield (iteration, avg_seconds, timer_name) from an AllTimersReadable.txt.

    Row format: '<it> <time>\\t<num>\\t<avg> <min> <max>\\t<name>'.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            try:
                it = int(float(parts[0].split()[0]))
                avg = float(parts[2].split()[0])
            except (ValueError, IndexError):
                continue
            rows.append((it, avg, parts[3]))
    return rows


def proc_label(path):
    if "temp_1" in path:
        return "np1"
    if "temp_2" in path:
        return "np2"
    return "?"


def main():
    files = []
    for root in sys.argv[1:]:
        for dirpath, _dirnames, filenames in os.walk(root):
            if "AllTimersReadable.txt" in filenames:
                files.append(os.path.join(dirpath, "AllTimersReadable.txt"))

    printed = False
    for path in sorted(files):
        rows = parse_rows(path)
        if not rows:
            continue
        maxit = max(r[0] for r in rows)
        # out_every may write a block at the final iteration AND the terminate
        # report fires there too, so the same (routine, bin) can appear twice at
        # maxit -- keep the largest (cumulative) time per (routine, bin).
        best = {}
        for it, avg, name in rows:
            if it != maxit:
                continue
            m = FUNC_RE.match(name)
            if m:
                key = (m.group(1), m.group(2))
                if avg > best.get(key, -1.0):
                    best[key] = avg
        if not best:
            continue
        test = os.path.basename(os.path.dirname(path))
        print("--- %-34s [%s, it %d] ---" % (test, proc_label(path), maxit))
        ranked = sorted(((avg, r, b) for (r, b), avg in best.items()),
                        reverse=True)
        for avg, routine, sched_bin in ranked[:TOP_N]:
            print("    %10.4f s  %-30s (%s)" % (avg, routine, sched_bin))
        printed = True

    if not printed:
        print("(no TimerReport AllTimersReadable.txt found)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
