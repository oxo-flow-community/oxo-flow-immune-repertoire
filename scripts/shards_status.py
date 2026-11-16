#!/usr/bin/env python3
"""shards_status.py — one-screen production dashboard for the 4 shards.

Answers what `oxo-flow cluster status` does not (issue #1085): how many rule
instances are planned/submitted/completed/failed per shard and per rule, which
samples are in flight on which node for how long, and the node temperatures
that enforce the <=20 CPU / thermal policy.

Data sources (all local, no scheduler accounting needed):
  shards/shardN/manifests/groups.tsv            sample count per shard
  shards/shardN/.oxo-flow/runs/latest/events.jsonl   SUBMITTED/COMPLETED/FAILED
  shards/shardN/.oxo-flow/runs/latest/jobs/*/status.json   per-instance state
  squeue -u wsx                                  node / elapsed / state
  sensors (via cluster_health_probe.sh, run as root)

Usage: python3 scripts/shards_status.py [--shards 1,2,3,4] [--no-temps]
Stdlib only.
"""
import argparse
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RULES = ("mixcr", "trust4", "standardize", "cleanup")


def sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True,
                              timeout=30).stdout
    except subprocess.TimeoutExpired:
        return ""


def planned_samples(shard):
    path = os.path.join(BASE, "shards", f"shard{shard}", "manifests", "groups.tsv")
    if not os.path.exists(path):
        return 0
    total = 0
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 1:
                total += len([s for s in parts[1].split(",") if s])
    return total


def latest_run(shard):
    root = os.path.join(BASE, "shards", f"shard{shard}", ".oxo-flow", "runs")
    if not os.path.isdir(root):
        return None
    runs = [os.path.join(root, d) for d in os.listdir(root)
            if os.path.isdir(os.path.join(root, d))]
    return max(runs, key=os.path.getmtime) if runs else None


def events_summary(run_dir):
    """Per-rule counters from events.jsonl (append-only, small)."""
    counters = defaultdict(Counter)
    jobs = {}
    path = os.path.join(run_dir, "events.jsonl")
    if not os.path.exists(path):
        return counters, jobs
    with open(path) as f:
        for line in f:
            try:
                e = json.loads(line)
            except json.JSONDecodeError:
                continue
            rule = e.get("rule", "")
            kind = rule.split("_")[0]
            if kind not in RULES:
                continue
            event = e.get("t")
            counters[kind][event] += 1
            if e.get("job"):
                jobs[e["job"]] = rule
    return counters, jobs


def squeue_map():
    """job id -> (state, elapsed, node, cpus)"""
    out = sh('squeue -u wsx -h -o "%i|%T|%M|%R|%C"')
    mapping = {}
    for line in out.splitlines():
        parts = line.split("|")
        if len(parts) == 5:
            mapping[parts[0]] = parts[1:]
    return mapping


def temps():
    out = sh("bash " + os.path.join(BASE, "scripts", "cluster_health_probe.sh"))
    rows = []
    for line in out.splitlines():
        rows.append(line.strip())
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", default="1,2,3,4")
    ap.add_argument("--no-temps", action="store_true")
    args = ap.parse_args()

    queue = squeue_map()
    print(f"{'shard':<6}{'planned':>8}{'submitted':>10}{'completed':>10}"
          f"{'failed':>7}{'inflight':>9}   per-rule (done/submitted)")
    totals = Counter()
    for s in args.shards.split(","):
        n = planned_samples(s)
        run = latest_run(s)
        if not run:
            print(f"shard{s:<2}  no run directory")
            continue
        counters, job_rule = events_summary(run)
        sub = sum(c["SUBMITTED"] for c in counters.values())
        comp = sum(c["COMPLETED"] for c in counters.values())
        fail = sum(c["FAILED"] + c["CANCELLED"] for c in counters.values())
        inflight = sum(1 for j in job_rule if j in queue)
        detail = "  ".join(
            f"{r}:{counters[r]['COMPLETED']}/{counters[r]['SUBMITTED']}"
            for r in RULES if counters[r])
        print(f"shard{s:<2}{n*4:>8}{sub:>10}{comp:>10}{fail:>7}{inflight:>9}   {detail}")
        totals["planned"] += n * 4
        totals["sub"] += sub
        totals["comp"] += comp
        totals["fail"] += fail

    print(f"{'TOTAL':<6}{totals['planned']:>8}{totals['sub']:>10}"
          f"{totals['comp']:>10}{totals['fail']:>7}")

    if not args.no_temps:
        print("\nnode health (temp °C / 1-min load / project CPUs):")
        for row in temps():
            print("  " + row)

    print("\nin flight:")
    for jid, (state, elapsed, node, cpus) in sorted(queue.items()):
        print(f"  {jid:<10}{state:<10}{elapsed:>10}  {node:<8}{cpus:>3} cpu")


if __name__ == "__main__":
    sys.exit(main())
