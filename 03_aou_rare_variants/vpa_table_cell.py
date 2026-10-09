# %% CELL 10 - thesis results tables for the variant analysis
# Reads the saved results files, so it does not need CELLS 3-6 rerun and is light enough to run
# while the PRS job is going. Prints a shareable version (counts 1-19 suppressed) that can be copied out.
import os, subprocess, numpy as np, pandas as pd

OUT = "gs://vpa-full-YOUR-WORKSPACE-ID/vpa_rerun"
os.makedirs("/tmp/vpa", exist_ok=True)
for f in ["results_all.csv", "results_mechanisms.csv"]:
    subprocess.run(["gsutil", "-q", "cp", f"{OUT}/{f}", f"/tmp/vpa/{f}"], check=True)   # /tmp, never "."
r_all = pd.read_csv("/tmp/vpa/results_all.csv")
r_mech = pd.read_csv("/tmp/vpa/results_mechanisms.csv")
print("subsets:", r_all["subset"].unique())
print("tiers:  ", r_all["tier"].unique())
print("groups: ", r_all["comparison"].unique(), "\n")

PRIMARY = "Primary (coding, >=1*)"
BLOCKS = [   # (label in the table, tier, panel name in results_all.csv)
    ("A. Revised panel (167 genes), primary tier",            PRIMARY,                  "Revised wide"),
    ("B. Revised panel, extended tier (+LoF, REVEL >= 0.75)", "Extended (+LoF, REVEL)", "Revised wide"),
    ("C. Negative control genes (10), primary tier",          PRIMARY,                  "Negative control"),
    ("D. Original panel (146 genes), primary tier",           PRIMARY,                  "Wide-147 (original)"),
    ("E. Original panel + GWAS-locus genes (170), primary",   PRIMARY,                  "Wide-171 (original)"),
]

def pct(n, d):
    return "<20" if 0 < n < 20 else f"{n:,} ({100 * n / d:.1f}%)"

def pfmt(p):
    return "<0.001" if p < 0.001 else (f"{p:.3f}" if p < 0.01 else f"{p:.2f}")

def row(r, extra=None):
    small = 0 < r["case_carriers"] < 20 or 0 < r["comp_carriers"] < 20
    d = {"comparison": r["comparison"], "n_comp": f"{r['n_comp']:,}",
         "cases": pct(r["case_carriers"], r["n_cases"]),
         "comp": pct(r["comp_carriers"], r["n_comp"]),
         "OR (95% CI)": "suppressed" if small else f"{r['OR']:.2f} ({r['lo']:.2f}-{r['hi']:.2f})",
         "p": "" if small else pfmt(r["p"])}
    if extra:
        d[extra] = "" if small else pfmt(r[extra])
    return d

sub = r_all[r_all["subset"] == r_all["subset"].iloc[0]]
table = []
for label, tier, panel in BLOCKS:
    blk = sub[(sub["tier"] == tier) & (sub["panel"] == panel)]
    assert len(blk), f"no rows for {tier} / {panel}"
    table += [{"block": label, **row(r)} for r in blk.to_dict("records")]

# Optional block: primary tier minus any variant with a benign / likely benign submission.
# Needs the CELL 3-5 objects in memory; skipped (with a message) if they are not there.
if all(k in globals() for k in ["v", "IDS", "carriers", "case_key", "PANELS", "COL", "or_ci"]):
    from scipy.stats import fisher_exact
    nb = v[v["plp"] & v["coding"] & (v["stars"] >= 1) & ~v["calls"].map(lambda s: any("benign" in c for c in s))]
    vids = set(nb.loc[nb[COL["gene"]].isin(PANELS["Revised wide"]), COL["vid"]])
    carrying = set(carriers.loc[carriers["vid"].isin(vids), "person_id"])
    a, n1 = len(IDS[case_key] & carrying), len(IDS[case_key])
    for grp in [k for k in IDS if k != case_key]:
        c, n0 = len(IDS[grp] & carrying), len(IDS[grp])
        o, lo, hi = or_ci(a, n1 - a, c, n0 - c)
        p = fisher_exact([[a, n1 - a], [c, n0 - c]])[1]
        table.append({"block": "A2. Revised panel, primary tier, no benign/likely benign submissions",
                      **row({"comparison": grp, "case_carriers": a, "n_cases": n1, "comp_carriers": c,
                             "n_comp": n0, "OR": o, "lo": lo, "hi": hi, "p": p})})
else:
    print("No-benign block skipped: run CELL 1 -> RESTORE -> 3 -> 4 -> 5a -> 5b first if you want it.\n")

main = pd.DataFrame(table)
pd.set_option("display.max_rows", 500, "display.width", 250, "display.max_colwidth", 80)
print(main.to_string(index=False))

# Mechanisms (primary tier, revised panel genes), BH q already in results_mechanisms.csv
mech = pd.DataFrame([{"mechanism": r["panel"].replace("Mechanism: ", ""), **row(r, extra="q")}
                     for r in r_mech.to_dict("records")])
print("\n", mech.to_string(index=False))

main.to_csv("/tmp/vpa/thesis_table_main_shareable.csv", index=False)
mech.to_csv("/tmp/vpa/thesis_table_mech_shareable.csv", index=False)
for f in ["thesis_table_main_shareable.csv", "thesis_table_mech_shareable.csv"]:
    subprocess.run(["gsutil", "-q", "cp", f"/tmp/vpa/{f}", f"{OUT}/{f}"], check=True)
