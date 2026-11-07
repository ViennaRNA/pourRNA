#!/usr/bin/env python3
"""Reproduce the sv11 refolding results of the pourRNA paper and compare them.

Entzian and Raden, pourRNA - a time- and memory-efficient approach for the guided
exploration of RNA energy landscapes, Bioinformatics 36(2), 2020, 462-469:
  * Table 1 (turner-99 parameters, shift moves) and Section 4.5,
  * Supplementary Table 2 (turner-04 parameters, shift moves).

All runs start in the co-transcriptionally trapped meta-stable structure of sv11
(Supplementary Section C) and use its minimum free energy structure as final structure,
maxE = -56 kcal/mol and shift moves, as stated in the paper.

The paper's Delta_bp (Section 3.3.4) is the number of base pairs by which the local minimum
of an explored basin may deviate from a direct path between start and final minimum.
Runs that use Delta_bp can be made in three modes (--modes, default add,detour):

  add     --max-bp-dist-add=Delta_bp   (the paper's values passed to the old option; this
                                        bounds the extra *moves*, i.e. 2 * deviating pairs,
                                        and allows only floor(Delta_bp/2) deviating pairs)
  detour  --max-bp-detour=Delta_bp     (Delta_bp deviating base pairs, as the paper defines it;
                                        the same as --max-bp-dist-add=2*Delta_bp)
  legacy  --max-bp-dist-add=Delta_bp --max-bp-dist-legacy
                                       (as pourRNA 1.0 computed the published results: the
                                        distance of the flooded basin decides, which adds one
                                        layer of basins beyond the bound; reproduces the
                                        published barriers and numbers of flooded basins)

For each run the script reports the energy barrier (lowest possible highest saddle on any
path from the start minimum to the final minimum through the explored basins, computed from
--saddle-file), the number of flooded basins |B_p| and the runtime, next to the published
values. Which mode reproduces the published barrier and |B_p| shows how the paper's runs were
parameterized.

Notes:
  * pourRNA allows only one thread when a final structure is given, the paper used up to
    8 threads; runtimes are therefore not directly comparable.
  * Bug fixes since the published version (e.g. inclusive and rounded energy limits, kBest
    tie-breaking, start basin, see ChangeLog) can change |B_p| slightly.
  * Several published runs took hours to days. By default only the fast runs are made
    (--runs quick); use --runs all, table1, table2 or a comma-separated list of run ids.

Example:
  scripts/reproduce_sv11_paper.py --pourrna src/pourRNA --runs quick --out sv11-repro
"""

import argparse
import csv
import heapq
import os
import re
import subprocess
import sys
import time

SV11 = ("GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGC"
        "GAGGUGACCCCCCGAAAAGGGGGGUUUCCCA")
# meta-stable conformation (Supplementary Section C and G.1)
META = ("(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........))))"
        ")))).((((((((.....)))))))).....")
# minimum free energy structures (last lines of Supplementary Listings 1 and 2)
MFE = {
    1: "(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))"
       "))))))))))))))))..))))))))..))).",  # turner-99, -97.60
    0: "(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))"
       "))))))))))))))))..))))))))..))).",  # turner-04, -97.70
}
MAX_E = -56

# id, source, energy model (0 = turner-04, 1 = turner-99), deltaE, Delta_bp, kBest,
# published barrier (None = not connected), published time, published |B_p|
RUNS = [
    ("t1-dE5",          "Table 1",   1, 5,  None, None, -62.3, "6.3 h",     143032),
    ("t1-dE6",          "Table 1",   1, 6,  None, None, -62.3, "7.8 days",  714359),
    ("t1-dE9",          "Table 1",   1, 9,  None, None, -62.3, "10.6 days", 771300),
    ("t1-dE5-k8",       "Table 1",   1, 5,  None, 8,    -62.3, "4.9 h",     111777),
    ("t1-dE6-k5",       "Table 1",   1, 6,  None, 5,    -62.3, "2.2 h",     40880),
    ("t1-dE6-bp0",      "Sec. 4.5.2", 1, 6, 0,    None, None,  "",          None),
    ("t1-dE6-bp6",      "Table 1",   1, 6,  6,    None, -59.4, "1.2 h",     20607),
    ("t1-dE6-bp8",      "Table 1",   1, 6,  8,    None, -62.3, "2 h",       44130),
    ("t1-dE6-bp10",     "Table 1",   1, 6,  10,   None, -62.3, "4.3 h",     89598),
    ("t1-dE6-bp10-k6",  "Table 1",   1, 6,  10,   6,    -62.3, "2.8 min",   2361),
    ("t2-dE5-bp6",      "Supp. Tab. 2", 0, 5,  6,  None, -56.6, "3.5 min",  7135),
    ("t2-dE5-bp10",     "Supp. Tab. 2", 0, 5,  10, None, -57.5, "14 min",   22496),
    ("t2-dE5",          "Supp. Tab. 2", 0, 5,  None, None, -57.5, "1.9 h",  77480),
    ("t2-dE6-bp5",      "Supp. Tab. 2", 0, 6,  5,  None, -58.5, "15 min",   9165),
    ("t2-dE8-bp5",      "Supp. Tab. 2", 0, 8,  5,  None, -58.5, "1.74 h",   16226),
    ("t2-dE10-bp0",     "Supp. Tab. 2", 0, 10, 0,  None, -58.5, "1.5 h",    3316),
]
PRESETS = {
    "quick": ["t1-dE6-bp0", "t1-dE6-bp10-k6", "t2-dE5-bp6", "t2-dE6-bp5", "t2-dE5-bp10"],
    "table1": [r[0] for r in RUNS if r[0].startswith("t1-")],
    "table2": [r[0] for r in RUNS if r[0].startswith("t2-")],
    "all": [r[0] for r in RUNS],
}


def read_saddles(path):
    """Undirected graph {minimum: {neighbor: saddle energy}} from a --saddle-file CSV."""
    graph = {}
    with open(path) as f:
        next(f, None)  # header
        for line in f:
            fields = [x.strip() for x in line.split(",")]
            if len(fields) != 8:
                continue
            a, b, saddle = fields[1], fields[4], float(fields[7])
            for x, y in ((a, b), (b, a)):
                old = graph.setdefault(x, {}).get(y)
                graph[x][y] = saddle if old is None else min(old, saddle)
    return graph


def barrier(graph, start, target):
    """Lowest possible highest saddle on a path from start to target (minimax Dijkstra)."""
    if start == target:
        return None
    best = {start: float("-inf")}
    queue = [(float("-inf"), start)]
    while queue:
        height, node = heapq.heappop(queue)
        if node == target:
            return height
        if height > best.get(node, float("inf")):
            continue
        for neighbor, saddle in graph.get(node, {}).items():
            h = max(height, saddle)
            if h < best.get(neighbor, float("inf")):
                best[neighbor] = h
                heapq.heappush(queue, (h, neighbor))
    return None  # not connected


def parse_output(text):
    def field(label):
        m = re.search(r"^" + re.escape(label) + r"\s*(\S*)\s*$", text, re.M)
        return m.group(1) if m else ""
    flooded = re.search(r"^Number of flooded minima:\s*(\d+)", text, re.M)
    return (field("The start state ends in basin:"), field("The final minimum is:"),
            int(flooded.group(1)) if flooded else None)


def run_one(args, run, mode):
    rid, _, model, delta_e, dbp, kbest = run[:6]
    name = rid if mode is None else f"{rid}-{mode}"
    wdir = os.path.join(args.out, name)
    os.makedirs(wdir, exist_ok=True)
    saddle_file = os.path.join(wdir, "saddles.csv")
    cmd = [args.pourrna, f"--sequence={SV11}", f"--start-structure={META}",
           f"--final-structure={MFE[model]}", f"--energy-model={model}", "--move-set=1",
           f"--max-energy={MAX_E}", f"--delta-e={delta_e}", "--max-threads=1",
           f"--saddle-file={saddle_file}"]
    if kbest is not None:
        cmd.append(f"--filter-best-k={kbest}")
    if mode == "add":
        cmd.append(f"--max-bp-dist-add={dbp}")
    elif mode == "detour":
        cmd.append(f"--max-bp-detour={dbp}")
    elif mode == "legacy":
        cmd += [f"--max-bp-dist-add={dbp}", "--max-bp-dist-legacy"]
    with open(os.path.join(wdir, "command.txt"), "w") as f:
        f.write(" ".join(f"'{c}'" if "(" in c else c for c in cmd) + "\n")

    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
        rc, out, err = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        rc, out, err = "timeout", "", ""
    seconds = time.time() - t0
    with open(os.path.join(wdir, "pourRNA.out"), "w") as f:
        f.write(out + err)

    result = {"run": name, "rc": rc, "seconds": round(seconds, 1), "barrier": "",
              "flooded": "", "start_min": "", "final_min": ""}
    if rc != 0:
        return result
    start_min, final_min, flooded = parse_output(out)
    result.update(start_min=start_min, final_min=final_min, flooded=flooded)
    if os.path.exists(saddle_file):
        b = barrier(read_saddles(saddle_file), start_min, final_min)
        result["barrier"] = "not connected" if b is None else f"{b:.2f}"
    return result


def fmt_time(seconds):
    if seconds == "" or seconds is None:
        return ""
    if seconds < 120:
        return f"{seconds:.0f} s"
    if seconds < 7200:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.2f} h"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__.split("\n", 2)[2])
    p.add_argument("--pourrna", default="pourRNA",
                   help="pourRNA binary (installed, so that --energy-model=1 finds its parameter file)")
    p.add_argument("--runs", default="quick",
                   help="quick, table1, table2, all, or comma-separated run ids (see --list)")
    p.add_argument("--modes", default="add,detour",
                   help="for runs with Delta_bp: add, detour and/or legacy (default: add,detour)")
    p.add_argument("--timeout", type=float, default=None, help="seconds per pourRNA run")
    p.add_argument("--out", default="sv11-repro", help="output directory")
    p.add_argument("--list", action="store_true", help="list the runs and exit")
    args = p.parse_args()

    by_id = {r[0]: r for r in RUNS}
    if args.list:
        for r in RUNS:
            dbp = "" if r[4] is None else f" Delta_bp={r[4]}"
            kb = "" if r[5] is None else f" kBest={r[5]}"
            pub = "not connected" if r[6] is None else f"{r[6]} kcal/mol, {r[7]}, |B_p|={r[8]}"
            print(f"{r[0]:16s} {r[1]:13s} turner-{'99' if r[2] else '04'} deltaE={r[3]}{dbp}{kb}"
                  f"  -> published: {pub}")
        return 0
    ids = PRESETS.get(args.runs, args.runs.split(","))
    unknown = [i for i in ids if i not in by_id]
    if unknown:
        p.error("unknown run id(s): " + ", ".join(unknown) + " (see --list)")
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    if any(m not in ("add", "detour", "legacy") for m in modes):
        p.error("--modes takes add, detour and/or legacy")

    os.makedirs(args.out, exist_ok=True)
    rows = []
    for rid in ids:
        run = by_id[rid]
        for mode in (modes if run[4] is not None else [None]):
            label = rid if mode is None else f"{rid} ({mode})"
            print(f"running {label} ...", file=sys.stderr, flush=True)
            res = run_one(args, run, mode)
            dbp = run[4]
            if mode == "add":
                option, pairs = f"--max-bp-dist-add={dbp}", dbp // 2
            elif mode == "detour":
                option, pairs = f"--max-bp-detour={dbp}", dbp
            elif mode == "legacy":
                option, pairs = f"--max-bp-dist-add={dbp} --max-bp-dist-legacy", f"{dbp // 2} (+1 layer)"
            else:
                option, pairs = "", ""
            pub_barrier = "not connected" if run[6] is None else f"{run[6]:.2f}"
            match = ""
            if res["barrier"]:
                match = "yes" if res["barrier"] == pub_barrier else "no"
            rows.append({
                "run": rid, "source": run[1], "params": f"turner-{'99' if run[2] else '04'}",
                "deltaE": run[3], "Delta_bp": "" if dbp is None else dbp,
                "kBest": "" if run[5] is None else run[5], "mode": mode or "",
                "option": option, "deviating_pairs_allowed": pairs,
                "published_barrier": pub_barrier, "barrier": res["barrier"] or f"rc={res['rc']}",
                "barrier_matches": match,
                "published_Bp": "" if run[8] is None else run[8], "Bp": res["flooded"],
                "published_time": run[7], "time": fmt_time(res["seconds"]),
            })
            print(f"  barrier {rows[-1]['barrier']} (published {pub_barrier}), |B_p| "
                  f"{res['flooded']} (published {rows[-1]['published_Bp']}), "
                  f"{rows[-1]['time']}", file=sys.stderr, flush=True)

    with open(os.path.join(args.out, "results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    head = ["run", "params", "deltaE", "Delta_bp", "kBest", "option", "pairs",
            "barrier (pub.)", "barrier", "|B_p| (pub.)", "|B_p|", "time (pub.)", "time"]
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for r in rows:
        lines.append("| " + " | ".join(str(x) for x in (
            r["run"], r["params"], r["deltaE"], r["Delta_bp"], r["kBest"], r["option"],
            r["deviating_pairs_allowed"], r["published_barrier"], r["barrier"],
            r["published_Bp"], r["Bp"], r["published_time"], r["time"])) + " |")
    table = "\n".join(lines)
    with open(os.path.join(args.out, "results.md"), "w") as f:
        f.write(table + "\n\npairs = deviating base pairs allowed per local minimum; "
                "time = single thread (the paper used up to 8 threads).\n")
    print(table)
    return 0


if __name__ == "__main__":
    sys.exit(main())
