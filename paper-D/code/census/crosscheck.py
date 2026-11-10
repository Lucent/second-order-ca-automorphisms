"""Cross-check of the censuses of commuting block maps with other SAT solvers.

For a rule f and a radius R the census is the set of all solutions of the CNF of centralizer.py
(commutation of a block map of radius R with F_f). Its exhaustiveness rests on the final answer of
the SAT solver that the CNF together with the blocking clauses of the solutions found is
unsatisfiable. This program repeats the complete enumeration with several solvers of PySAT, which
share the encoding and nothing else, and compares

  (a) the solution sets of the solvers with one another, as sets of tables;
  (b) the solution set with the shipped log of the census: the numbers of clauses and of solutions,
      the counts named / const / other, and the set itself whenever the log determines it (it does
      when the log lists the tables, and when it reports no solution other than the maps
      sigma^a F^b with |a|+|b| <= R and the constants, whose tables are built here directly).

Every table returned by a solver is rechecked on all 4^(2R+3) windows with NumPy.

With --rup, the unsatisfiable CNF (the encoding and one blocking clause for each solution) is in
addition given to Glucose 4 with proof logging, and the clausal proof is checked by the forward
reverse-unit-propagation checker below, which is independent of the solvers (instances with very
many solutions, whose blocking clauses are too large for it, are skipped and marked in the output).

    python3 crosscheck.py --R 1,2,3 --jobs 12                 # all rules of logs/census_R123_summary.txt
    python3 crosscheck.py --R 4 --rules 110,30 --solvers glucose4,cadical153 --jobs 4
    python3 crosscheck.py --R 1,2 --rules 110 --rup
"""
import sys, os, re, ast, time, argparse, hashlib
from multiprocessing import Pool
import numpy as np
from pysat.solvers import Solver
from centralizer import make_rule, encode, named_maps
from classify_r4 import check_commutes_np_R

HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = os.path.join(HERE, "logs")
RUP_LIMIT = 10 ** 7      # --rup is skipped when (number of solutions) x (table entries) exceeds this
HEAD = re.compile(r"=== rule (\d+)R,? radius (?:R=)?(\d+)(?: \[\w+\])?: (\d+) variables, (\d+) clauses,(?:.*,)? (\d+) solutions( \(CAPPED\))?")


def shipped(R):
    """(rule, R) -> dict(nclauses, nsol, capped, counts or None, tables or None), from the shipped logs."""
    files = ["census_R123_summary.txt", "census_110_30_R123.txt"] if R <= 3 else \
            ["census_R4_summary.txt", "census_110_R4.txt", "census_30_R4.txt"]
    out = {}
    for fn in files:
        cur = None
        for line in open(os.path.join(LOGS, fn)):
            m = HEAD.search(line)
            if m:
                key = (int(m.group(1)), int(m.group(2)))
                cur = out.setdefault(key, {"counts": None, "tables": None})
                cur.update(nclauses=int(m.group(4)), nsol=int(m.group(5)), capped=bool(m.group(6)))
                continue
            if cur is None: continue
            d = re.search(r"\{'named'.*?\}", line)
            if d: cur["counts"] = ast.literal_eval(d.group(0))
            t = re.search(r"table=([0-3]+)\s*$", line)
            if t:
                if cur["tables"] is None: cur["tables"] = set()
                cur["tables"].add(bytes(int(ch) for ch in t.group(1)))
    return out


class Collect:
    def __init__(self): self.clauses = []
    def add_clause(self, c): self.clauses.append(c)


def enumerate_with(name, rule, R, cap):
    """all solutions (as bytes, one byte per table entry) with the solver `name`; (set, capped, nclauses)."""
    f, F_loc, _ = make_rule(rule)
    n = 2 * R + 1; nvar = 2 * 4 ** n
    sols = set(); capped = False
    with Solver(name=name) as s:
        nv, nc = encode(R, f, F_loc, s)
        while s.solve():
            model = s.get_model()
            bits = np.array(model[:nvar]) > 0
            sols.add((bits[0::2] * 2 + bits[1::2]).astype(np.uint8).tobytes())
            s.add_clause([-lit for lit in model[:nvar]])
            if len(sols) >= cap: capped = True; break
    return sols, capped, nc


def rup_check(nvar, clauses, proof):
    """Forward check of a clausal proof by reverse unit propagation, with two watched literals.
    Deletion lines are ignored (keeping a clause cannot make a later lemma fail to be implied).
    Returns (ok, number of lemmas checked); ok means that every lemma, and at the end the empty
    clause, follows from the earlier clauses by unit propagation."""
    val = [0] * (nvar + 1)                 # 0 unassigned, 1 true, -1 false
    watch = {}; db = []; trail = []
    def value(l): return val[l] if l > 0 else -val[-l]
    def assign(l):
        val[abs(l)] = 1 if l > 0 else -1; trail.append(l)
    def propagate(start):
        """returns True on conflict"""
        i = start
        while i < len(trail):
            fl = -trail[i]; i += 1
            ws = watch.get(fl, []); keep = []; conflict = False
            for k, ci in enumerate(ws):
                c = db[ci]
                if c[0] == fl: c[0], c[1] = c[1], c[0]
                if value(c[0]) == 1: keep.append(ci); continue
                for j in range(2, len(c)):
                    if value(c[j]) != -1:
                        c[1], c[j] = c[j], c[1]; watch.setdefault(c[1], []).append(ci); break
                else:
                    keep.append(ci)
                    if value(c[0]) == -1:
                        keep.extend(ws[k + 1:]); conflict = True; break
                    assign(c[0])
            watch[fl] = keep
            if conflict: return True
        return False
    def add(c):
        """add a clause at the top level; returns True if the empty clause is derived"""
        c = sorted(set(c), key=lambda l: value(l) == -1)       # non-false literals first
        s = set(c)
        if any(-l in s for l in c): return False
        if not c or value(c[0]) == -1: return True
        if len(c) == 1 or value(c[1]) == -1:
            if value(c[0]) == 1: return False
            n0 = len(trail); assign(c[0])
            if len(c) > 1:
                db.append(c); watch.setdefault(c[0], []).append(len(db) - 1); watch.setdefault(c[1], []).append(len(db) - 1)
            return propagate(n0)
        db.append(c); watch.setdefault(c[0], []).append(len(db) - 1); watch.setdefault(c[1], []).append(len(db) - 1)
        return False
    done = False
    for c in clauses:
        if add(list(c)): done = True; break
    nlem = 0
    if not done:
        for line in proof:
            line = line.strip()
            if not line or line.startswith("d"): continue
            lemma = [int(x) for x in line.split()][:-1]
            nlem += 1
            n0 = len(trail); conflict = False
            for l in lemma:
                if value(l) == 1: conflict = True; break
                if value(l) == 0: assign(-l)
            if not conflict: conflict = propagate(n0)
            for l in trail[n0:]: val[abs(l)] = 0
            del trail[n0:]
            if not conflict: return False, nlem
            if add(lemma): done = True; break
    return done, nlem


def task(args):
    rule, R, solvers, cap, ship, rup = args
    t0 = time.time()
    f, F_loc, Finv_loc = make_rule(rule)
    res = {}
    for name in solvers:
        res[name] = enumerate_with(name, rule, R, cap)
    first = res[solvers[0]]; S = first[0]; nc = first[2]
    msg = [f"rule {rule}R R={R}: {2 * 4 ** (2 * R + 1)} variables, {nc} clauses |"]
    msg.append(" ".join(f"{name}: {len(res[name][0])}{' (CAPPED)' if res[name][1] else ''}" for name in solvers) + " |")
    capped = any(r[1] for r in res.values())
    ok = True
    if capped:
        msg.append("cap reached, sets not compared |")
    else:
        same = all(res[name][0] == S for name in solvers)
        ok &= same
        msg.append(f"sets identical: {same} |")
        commute = all(check_commutes_np_R(np.frombuffer(t, dtype=np.uint8), F_loc, R) for t in S)
        ok &= commute
        names = {bytes(k) for k in named_maps(R, F_loc, Finv_loc)}
        consts = {bytes([s]) * 4 ** (2 * R + 1) for s in range(4)}
        fixed = {bytes([3 * b]) * 4 ** (2 * R + 1) for b in (0, 1) if f(b, b, b) == 0}
        n_named = len(S & names); n_const = len(S & consts - names); n_other = len(S) - n_named - n_const
        ok &= names <= S and (S & consts) == fixed
        msg.append(f"all commute: {commute}; named {n_named}/{len(names)} const {n_const} other {n_other} |")
        h = hashlib.sha256(b"".join(sorted(S))).hexdigest()[:16]
        sh = ship.get((rule, R))
        if sh is None:
            msg.append("no shipped log |")
        elif sh["capped"]:
            msg.append("shipped log capped |")
        else:
            agree = sh["nclauses"] == nc and sh["nsol"] == len(S)
            if sh["counts"] is not None:
                c = sh["counts"]
                agree &= (c["named"], c["const"], c["other"]) == (n_named, n_const, n_other)
            if sh["tables"] is not None:
                agree &= sh["tables"] == S; how = "set equal (tables in the log)"
            elif sh["counts"] is not None and sh["counts"]["other"] == 0:
                agree &= S == names | fixed; how = "set equal (determined by the log)"
            elif sh["counts"] is not None:
                how = "counts equal"
            else:
                how = "number of solutions equal"
            ok &= agree
            msg.append(f"shipped log: {how if agree else 'DISAGREE'} |")
        msg.append(f"sha256 {h} |")
        if rup and len(S) * len(next(iter(S))) > RUP_LIMIT:
            msg.append(f"UNSAT proof not checked ({len(S)} solutions) |")
        elif rup:
            col = Collect(); encode(R, f, F_loc, col)
            nvar = 2 * 4 ** (2 * R + 1)
            block = [[-(2 * e + 1) if t[e] & 2 else (2 * e + 1) for e in range(len(t))] +
                     [-(2 * e + 2) if t[e] & 1 else (2 * e + 2) for e in range(len(t))] for t in sorted(S)]
            cnf = col.clauses + block
            with Solver(name="glucose4", bootstrap_with=cnf, with_proof=True) as s:
                sat = s.solve(); proof = s.get_proof() or []
            good, nlem = (False, 0) if sat else rup_check(nvar, cnf, proof)
            ok &= good
            msg.append(f"glucose4 UNSAT proof: {len(proof)} lines, {nlem} lemmas, RUP check {'passed' if good else 'FAILED'} |")
    msg.append(f"{'OK' if ok else 'MISMATCH'} [{time.time() - t0:.0f}s]")
    return rule, R, ok, " ".join(msg)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", default=None, help="default: all rules of the shipped summary log")
    ap.add_argument("--R", default="1,2,3")
    ap.add_argument("--solvers", default="glucose4,minisat22,cadical153")
    ap.add_argument("--max", type=int, default=5000)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--rup", action="store_true")
    a = ap.parse_args()
    Rs = [int(x) for x in a.R.split(",")]
    solvers = a.solvers.split(",")
    tasks = []
    for R in Rs:
        ship = shipped(R)
        rules = [int(x) for x in a.rules.split(",")] if a.rules else sorted({r for r, RR in ship if RR == R})
        for rule in rules:
            sh = ship.get((rule, R))
            if sh is not None and sh["capped"] and not a.rules:
                print(f"rule {rule}R R={R}: cap reached in the shipped log, skipped", flush=True); continue
            tasks.append((rule, R, solvers, a.max, ship, a.rup))
    print(f"# solvers: {', '.join(solvers)}; {len(tasks)} instances", flush=True)
    bad = 0
    with Pool(a.jobs) as pool:
        for rule, R, ok, line in pool.imap(task, tasks):
            print(line, flush=True); bad += not ok
    print(f"# {len(tasks)} instances, {bad} mismatches", flush=True)
