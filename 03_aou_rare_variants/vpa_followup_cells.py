# =====================================================================================
# Follow-up analyses for Aim 2. Run after vpa_pipeline.py (needs v, QUAL, IDS, carriers,
# PANELS, MECHANISMS, COL, or_ci and save in memory). Each '# %%' block is one notebook cell.
# =====================================================================================

# %% CELL 11 - variants driving the comparison-group carrier rate
# DIAGNOSTIC - which variants drive the comparison-group carrier rate? (output safe to export: <20 suppressed)
q = QUAL["Primary (coding, >=1*)"]
q147 = q[q["gene"].isin(PANELS["Wide-147 (original)"])]
c = carriers[carriers["vid"].isin(set(q147["vid"]))]
comp_ids, case_ids = IDS["Complete comparison group"], IDS[case_key]
per = c.groupby("vid")["person_id"].agg(
    comp=lambda s: len(set(s) & comp_ids), cases=lambda s: len(set(s) & case_ids)).reset_index()
info = v[[COL["vid"], COL["gene"], COL["af"], COL["clinvar"], COL["stars"], COL["consequence"]]] \
         .drop_duplicates(COL["vid"]).rename(columns={COL["vid"]: "vid"})
per = per.merge(info, on="vid", how="left").sort_values("comp", ascending=False)
per["comp % of group"] = (100 * per["comp"] / len(comp_ids)).round(2)
for col in ["comp", "cases"]:
    per[col] = per[col].map(lambda n: "<20" if 0 < n < 20 else n)
pd.set_option("display.width", 250, "display.max_colwidth", 45)
print(per.head(25).to_string(index=False))

print("\nTop 10 variants' share of all comparison carriers:",
      round(100 * c[c["vid"].isin(per.head(10)["vid"])]["person_id"].pipe(lambda s: len(set(s) & comp_ids))
            / len(set(c["person_id"]) & comp_ids), 1), "%")

# %% CELL 12 - mechanisms x all four comparison groups, with and without benign-submission variants
# MECHANISMS x ALL FOUR comparison groups, two variant definitions (output safe to export)
from statsmodels.stats.multitest import multipletests
benign_vids = set(v.loc[v["calls"].map(lambda s: len(s & {"benign", "likely benign"}) > 0), COL["vid"]])
q = QUAL["Primary (coding, >=1*)"]
rows = []
for defn, drop in [("All P/LP", set()), ("No benign calls", benign_vids)]:
    for mech, genes in MECHANISMS.items():
        vids = set(q.loc[q["gene"].isin(genes), "vid"]) - drop
        carrying = set(carriers.loc[carriers["vid"].isin(vids), "person_id"])
        a, n1 = len(IDS[case_key] & carrying), len(IDS[case_key])
        for grp in [k for k in IDS if k != case_key]:
            c, n0 = len(IDS[grp] & carrying), len(IDS[grp])
            o, lo, hi = or_ci(a, n1 - a, c, n0 - c)
            rows.append(dict(definition=defn, mechanism=mech[:40], comparison=grp, a=a, c=c, n1=n1, n0=n0,
                             OR=o, lo=lo, hi=hi, p=fisher_exact([[a, n1 - a], [c, n0 - c]])[1]))
mt = pd.DataFrame(rows)
mt["q"] = mt.groupby(["definition", "comparison"])["p"].transform(lambda p: multipletests(p, method="fdr_bh")[1])
small = mt["a"].between(1, 19) | mt["c"].between(1, 19)
mt["cases %"] = np.where(small, "", (100 * mt["a"] / mt["n1"]).round(1).astype(str))
mt["comp %"]  = np.where(small, "", (100 * mt["c"] / mt["n0"]).round(1).astype(str))
mt["OR (95% CI)"] = np.where(small, "suppressed", mt.apply(lambda r: f"{r.OR:.2f} ({r.lo:.2f}-{r.hi:.2f})", axis=1))
fmt = lambda x: "<0.001" if x < 0.001 else f"{x:.3f}"
mt["p"] = np.where(small, "", mt["p"].map(fmt)); mt["q (BH)"] = np.where(small, "", mt["q"].map(fmt))
pd.set_option("display.width", 250, "display.max_rows", 200)
print(mt[["definition", "mechanism", "comparison", "cases %", "comp %", "OR (95% CI)", "p", "q (BH)"]].to_string(index=False))

# %% CELL 13 - gene-level results for each comparison group
# GENE-LEVEL x EACH comparison group - FULL COUNTS (view inside the Workbench only)
from statsmodels.stats.multitest import multipletests
pd.set_option("display.width", 250, "display.max_rows", 500)
q = QUAL["Primary (coding, >=1*)"]
cases_s = IDS[case_key]; n1 = len(cases_s)
gene_carriers = {g: set(carriers.loc[carriers["vid"].isin(set(q.loc[q["gene"] == g, "vid"])), "person_id"])
                 for g in PANELS["Revised wide"]}
allg = []
for grp in [k for k in IDS if k != case_key]:
    comp_s = IDS[grp]; n0 = len(comp_s); rows = []
    for g, carrying in gene_carriers.items():
        a, c = len(cases_s & carrying), len(comp_s & carrying)
        if a + c == 0:
            continue
        o, lo, hi = or_ci(a, n1 - a, c, n0 - c)
        rows.append(dict(comparison=grp, gene=g, cases=a, case_pct=round(100*a/n1, 2), comp=c,
                         comp_pct=round(100*c/n0, 3), OR=round(o, 2), lo=round(lo, 2), hi=round(hi, 2),
                         p=fisher_exact([[a, n1 - a], [c, n0 - c]])[1]))
    d = pd.DataFrame(rows); d["q"] = multipletests(d["p"], method="fdr_bh")[1]
    allg.append(d.sort_values("p"))
    print(f"\n=== {grp} (n = {n0:,}) - top 20 by p ===")
    print(d.sort_values("p").head(20).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("q < 0.05:", list(d.loc[d["q"] < 0.05, "gene"]) or "none")
save(pd.concat(allg), "results_gene_level_by_group.csv", index=False)

