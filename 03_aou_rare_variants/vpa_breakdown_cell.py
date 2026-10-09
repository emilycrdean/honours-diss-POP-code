# %% CELL 9 - full breakdown at the primary tier: panels, mechanisms and every gene, for each comparison group
# Run in the PVA notebook after cells 1-6 (if the kernel restarted: cells 1, 3, 4, 5b, 6 - 5b reads the saved
# carriers, no query). Prints FULL counts for you to look at inside the Workbench, and saves two files:
#   breakdown_full.csv        - every count (Workbench only, never download or paste)
#   breakdown_shareable.csv   - counts 1-19 shown as "<20" and their ORs hidden (safe to take out)
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

for name in ["QUAL", "IDS", "carriers", "PANELS", "MECHANISMS", "case_key", "or_ci", "save"]:
    assert name in globals(), f"{name} is missing - run cells 1 to 6 first"

TIER = "Primary (coding, >=1*)"
groups = [k for k in IDS if k != case_key]
q = QUAL[TIER][["vid", "gene"]].drop_duplicates()
by_gene = carriers.merge(q, on="vid").groupby("gene")["person_id"].apply(set)   # who carries a qualifying variant, per gene

def stats(level, name, genes, mechanism=""):
    carrying = set().union(*[by_gene.get(g, set()) for g in genes]) if genes else set()
    a, n1 = len(IDS[case_key] & carrying), len(IDS[case_key])
    rows = []
    for grp in groups:
        c, n0 = len(IDS[grp] & carrying), len(IDS[grp])
        o, lo, hi = or_ci(a, n1 - a, c, n0 - c)
        rows.append({"level": level, "name": name, "mechanism": mechanism, "comparison": grp,
                     "case_carriers": a, "n_cases": n1, "case_pct": 100 * a / n1,
                     "comp_carriers": c, "n_comp": n0, "comp_pct": 100 * c / n0,
                     "OR": o, "lo": lo, "hi": hi, "p": fisher_exact([[a, n1 - a], [c, n0 - c]])[1]})
    return rows

gene_mech = {g: m for m, gs in MECHANISMS.items() for g in gs}
rows = []
for p_name in ["Wide-147 (original)", "Revised wide", "Core-10", "Negative control"]:
    rows += stats("Panel", p_name, PANELS[p_name])
for m_name, genes in MECHANISMS.items():
    rows += stats("Mechanism", m_name, genes, m_name)
all_genes = sorted(set(PANELS["Wide-147 (original)"]) | set(PANELS["Revised wide"]) | set(PANELS["Negative control"]))
for g in all_genes:
    mech = gene_mech.get(g, "Negative control" if g in PANELS["Negative control"] else "Removed in revised panel")
    rows += stats("Gene", g, [g], mech)
bd = pd.DataFrame(rows)

# Benjamini-Hochberg correction for each comparison: across all 9 mechanisms, and across all genes
bd["q"] = np.nan
for level in ["Mechanism", "Gene"]:
    for grp in groups:
        m = (bd["level"] == level) & (bd["comparison"] == grp)
        bd.loc[m, "q"] = multipletests(bd.loc[m, "p"], method="fdr_bh")[1]

save(bd, "breakdown_full.csv", index=False)
share = bd.copy()
small = share["case_carriers"].between(1, 19) | share["comp_carriers"].between(1, 19)
for col in ["case_carriers", "comp_carriers"]:
    share[col] = share[col].map(lambda n: "<20" if 0 < n < 20 else n)
share.loc[small, ["case_pct", "comp_pct", "OR", "lo", "hi", "p", "q"]] = np.nan
save(share, "breakdown_shareable.csv", index=False)

# full view for you, inside the Workbench only
def fmt(r):
    q_txt = "" if pd.isna(r.q) else f"  q={r.q:.3f}"
    return (f"{r.case_carriers}/{r.n_cases} ({r.case_pct:.1f}%) vs {r.comp_carriers}/{r.n_comp} ({r.comp_pct:.2f}%)  "
            f"OR {r.OR:.2f} ({r.lo:.2f}-{r.hi:.2f})  p={r.p:.3g}{q_txt}")
pd.set_option("display.width", 250, "display.max_colwidth", 200, "display.max_rows", 2000)
for level in ["Panel", "Mechanism", "Gene"]:
    print(f"\n==================== {level.upper()} (WORKBENCH ONLY - full counts) ====================")
    sub = bd[bd["level"] == level]
    if level == "Gene":
        sub = sub[sub["case_carriers"] + sub["comp_carriers"] > 0].sort_values(["mechanism", "name"])
    for (mech, name), part in sub.groupby(["mechanism", "name"], sort=False):
        head = f"{name}" + (f"   [{mech}]" if level == "Gene" else "")
        print(f"\n{head}")
        for r in part.itertuples():
            print(f"   vs {r.comparison:28s} {fmt(r)}")
print("\nsaved breakdown_full.csv (Workbench only) and breakdown_shareable.csv (safe to take out) to", OUT)
