#!/usr/bin/python3
import os
import sys
import math
import re
import datetime
from datetime import datetime
import json
import glob
import subprocess

from joblib import Parallel, delayed


class Graph:
    def __init__(self, n, edges):
        self.num_vertices = n
        self.edges = edges

class DisjointSet:
    def __init__(self, n):
        self.parents = list(range(n))
        self.ranks = [0] * n
    def find(self, x):
        if self.parents[x] == x: return x
        self.parents[x] = self.find(self.parents[x])
        return self.parents[x]
    def unite(self, x, y):
        x = self.find(x)
        y = self.find(y)
        if x == y: return x
        if self.ranks[x] < self.ranks[y]:
            self.parents[x] = y
            return y
        if self.ranks[x] > self.ranks[y]:
            self.parents[y] = x
            return x
        self.parents[y] = x
        self.ranks[x] += 1
        return x
    def same(self, x, y):
        return self.find(x) == self.find(y)

def parse_graph(lines):
    n, m = map(int, lines[0].split())
    edges = set()
    for line in lines[1:]:
        if line.strip() == "":
            continue
        u, v = map(int, line.split())
        u -= 1
        v -= 1
        if u > v:
            u, v = v, u
        edges.add((u, v))
    return Graph(n, edges)

def parse_problem(in_data):
    lines = in_data.split("\n")
    gn, gm = map(int, lines[0].split())
    g = parse_graph(lines[0:gm+1])
    g_emb = parse_graph(lines[gm+1:])
    return g, g_emb

def parse_solution(n, stdout):
    result = [None] * n
    for u, line in enumerate(stdout.split("\n")):
        if line.strip() == "":
            continue
        tokens = list(map(lambda x: int(x) - 1, line.split()))
        result[u] = set(tokens[1:])
    return result

def evaluate(g, g_emb, solution):
    # validate solution
    disjoint_set = DisjointSet(g_emb.num_vertices)
    mapping = [-1] * g_emb.num_vertices
    for i, phi in enumerate(solution):
        if phi is None or len(phi) == 0:
            raise ValueError("empty group")
        for j in phi:
            if mapping[j] != -1:
                raise ValueError("duplicated mapping")
            mapping[j] = i
    for u, v in g_emb.edges:
        if mapping[u] == mapping[v]:
            disjoint_set.unite(u, v)
    for i, phi in enumerate(solution):
        phi_list = list(phi)
        for j in phi_list:
            if not disjoint_set.same(phi_list[0], j):
                raise ValueError("disjointed group")
    # compute score
    remaining = set(g.edges)
    total = 5000 + len(g.edges) * 100 + 100000
    score = 5000
    for u_emb, v_emb in g_emb.edges:
        u = mapping[u_emb]
        v = mapping[v_emb]
        if u < 0 or v < 0:
            continue
        if u > v:
            u, v = v, u
        if (u, v) in remaining:
            score += 100
            remaining.remove((u, v))
    if len(remaining) == 0:
        score += 100000
    for phi in solution:
        score -= len(phi) - 1
    return score, total

def run(executable, fname):
    with open(fname, "r") as f:
        in_data = f.read()
    start_dt = datetime.now()
    proc = subprocess.Popen(
            [executable, fname],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate(in_data.encode("utf-8"))
    for line in stderr.decode("utf-8").strip().split("\n"):
        sys.stderr.write("{}: {}\n".format(fname, line))
    end_dt = datetime.now()
    exec_time = (end_dt - start_dt).total_seconds()
    g, g_emb = parse_problem(in_data)
    emb_wh = int(math.sqrt(g_emb.num_vertices))
    solution = parse_solution(g.num_vertices, stdout.decode("utf-8"))
    score, total = evaluate(g, g_emb, solution)
    return {
        "name": fname,
        "score": score,
        "total": total,
        "g": g,
        "g_emb": {"width": emb_wh, "height": emb_wh},
        "solution": solution,
        "time": exec_time
    }

def to_json(obj):
    if isinstance(obj, Graph):
        return obj.__dict__
    if isinstance(obj, set):
        return list(obj)
    else:
        raise TypeError(repr(obj) + " is not JSON serializable")

def main():
    outdir = "stats"
    os.makedirs(outdir, exist_ok=True)
    results = Parallel(n_jobs=-1)([
        delayed(run)("./a.out", fname)
        for fname in list(glob.glob("../dataset/*.in"))])
    summary = []
    score_sum = 0
    total_sum = 0
    max_time = 0
    for result in results:
        fname = os.path.basename(result["name"])
        score, total, time = result["score"], result["total"], result["time"]
        print("{}: {} / {} ({}s)".format(
            result["name"], score, total, result["time"]))
        score_sum += score
        total_sum += total
        max_time = max(max_time, time)
        with open(os.path.join(outdir, fname + ".json"), "w") as f:
            json.dump(result, f, separators=(",", ":"), default=to_json)
        summary.append({
            "name": fname,
            "score": score,
            "total": total,
            "time": time,
            "num_vertices": result["g"].num_vertices,
            "num_edges": len(result["g"].edges),
            "embed_width": result["g_emb"]["width"],
            "embed_height": result["g_emb"]["height"],
        })
    print("{} / {} (max: {}s)".format(score_sum, total_sum, max_time))

if __name__ == "__main__":
    main()
