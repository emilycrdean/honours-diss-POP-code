# %% PRS-G3 - genome-wide score (PGS002288), one chromosome at a time - has a real tendency to crash or run into billing errors
# Each chromosome is saved to the bucket as soon as it finishes. Re-running this cell skips
# chromosomes already saved, so if the kernel, cluster or connection dies, just run it again and it'll pick back up.
# It will run happily in the background but if going too slow add more workers. 
import os, subprocess, time
import pandas as pd
import hail as hl

P = "YOUR-WORKSPACE-ID"
PRS_OUT = "gs://YOUR-WORKSPACE-ID/prs"
GW_OUT = f"{PRS_OUT}/pgs002288_by_chr"
MT_PATH = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt"
LOCAL_PGS = "/tmp/PGS002288_hmPOS_GRCh38.txt.gz"                      # /tmp, never the notebook folder

def gs_exists(path):
    return subprocess.run(["gsutil", "-q", "stat", path]).returncode == 0

hl.init(idempotent=True, gcs_requester_pays_configuration=P)

# weights (same filtering as before)
if not os.path.exists(LOCAL_PGS):
    subprocess.run(["gsutil", "-q", "cp", f"{PRS_OUT}/PGS002288_hmPOS_GRCh38.txt.gz", LOCAL_PGS], check=True)
pgs = pd.read_csv(LOCAL_PGS, sep="\t", comment="#", dtype=str)
other = "other_allele" if "other_allele" in pgs.columns else "hm_inferOtherAllele"
pgs = pgs.dropna(subset=["hm_chr", "hm_pos", "effect_allele", other, "effect_weight"])
pgs = pgs[pgs[other].str.fullmatch(r"[ACGT]+") & pgs["effect_allele"].str.fullmatch(r"[ACGT]+")]
pgs["effect_weight"] = pgs["effect_weight"].astype(float)
order = [str(i) for i in range(1, 23)] + ["X"]
chroms = [c for c in order if c in set(pgs["hm_chr"])]
print(f"{len(pgs):,} weighted variants on chromosomes: {', '.join(chroms)}")
pd.Series(chroms).to_csv("/tmp/pgs_chroms.txt", index=False, header=False)
subprocess.run(["gsutil", "-q", "cp", "/tmp/pgs_chroms.txt", f"{GW_OUT}/expected_chromosomes.txt"], check=True)

# everyone in the three PRS cohorts
ids = set()
for name in ["early_pop", "later_pop", "controls"]:
    if not os.path.exists(f"/tmp/prs_{name}.csv"):
        subprocess.run(["gsutil", "-q", "cp", f"{PRS_OUT}/{name}.csv", f"/tmp/prs_{name}.csv"], check=True)
    ids |= set(pd.read_csv(f"/tmp/prs_{name}.csv", dtype=str, usecols=["person_id"])["person_id"])
print(f"women to score: {len(ids):,}")

# chr22 was already scored in the pilot - telling the thing to reuse it instead of running it again
if gs_exists(f"{PRS_OUT}/pgs002288_chr22.tsv.bgz") and not gs_exists(f"{GW_OUT}/chr22.tsv.bgz"):
    subprocess.run(["gsutil", "-q", "cp", f"{PRS_OUT}/pgs002288_chr22.tsv.bgz", f"{GW_OUT}/chr22.tsv.bgz"], check=True)

total_w = len(pgs)
done_w = sum(len(pgs[pgs["hm_chr"] == c]) for c in chroms if gs_exists(f"{GW_OUT}/chr{c}.tsv.bgz"))
t_start = time.time(); ran_w = 0
for c in chroms:
    out = f"{GW_OUT}/chr{c}.tsv.bgz"
    if gs_exists(out):
        print(f"chr{c}: already saved - skipped")
        continue
    t0 = time.time()
    w = pgs[pgs["hm_chr"] == c]
    both = pd.concat([                                                 
        pd.DataFrame({"pos": w["hm_pos"].astype(int).values, "a0": w[other].values,
                      "a1": w["effect_allele"].values, "w": w["effect_weight"].values, "flip": False}),
        pd.DataFrame({"pos": w["hm_pos"].astype(int).values, "a0": w["effect_allele"].values,
                      "a1": w[other].values, "w": w["effect_weight"].values, "flip": True})])
    wt = hl.Table.from_pandas(both)
    wt = wt.annotate(locus=hl.locus(f"chr{c}", wt.pos, reference_genome="GRCh38"),
                     alleles=hl.array([wt.a0, wt.a1]))
    wt = wt.key_by("locus", "alleles").select("w", "flip")

    mt = hl.read_matrix_table(MT_PATH)
    sample = list(mt.col_key)[0]
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(f"chr{c}", reference_genome="GRCh38")])
    mt = mt.filter_cols(hl.literal(ids).contains(mt[sample]))
    mt = mt.annotate_rows(pgs=wt[mt.row_key])
    mt = mt.filter_rows(hl.is_defined(mt.pgs))
    mt = mt.select_entries(d=hl.if_else(mt.pgs.flip, 2 - mt.GT.n_alt_alleles(), mt.GT.n_alt_alleles()))
    mt = mt.annotate_rows(mean_d=hl.agg.mean(mt.d))                    # missing calls get the variant's mean
    mt = mt.annotate_cols(score=hl.agg.sum(hl.or_else(mt.d, mt.mean_d) * mt.pgs.w),
                          n_called=hl.agg.count_where(hl.is_defined(mt.d)))
    cols = mt.cols().key_by()
    cols = cols.select(person_id=hl.str(cols[sample]), score=cols.score, n_called=cols.n_called)
    cols.export(out)

    ran_w += len(w); done_w += len(w)
    hrs = (time.time() - t_start) / 3600
    left = hrs / ran_w * (total_w - done_w) if ran_w else float("nan")
    print(f"chr{c}: {len(w):,} weights, {(time.time() - t0) / 60:.0f} min | "
          f"{100 * done_w / total_w:.0f}% of the score done, about {left:.1f} h to go")
print("all chromosomes saved" if done_w == total_w else "stopped before the end - run the cell again to continue")


# %% PRS-G4 - add up the chromosomes and rerun the PRS analysis with the genome-wide score
# Self-contained: works in a fresh kernel -- means you don't have to reinit hail yay
import os, subprocess
import numpy as np, pandas as pd
import statsmodels.api as sm

P = "YOUR-WORKSPACE-ID"
PRS_OUT = "gs://vpa-full-YOUR-WORKSPACE-ID/prs"
GW_OUT = f"{PRS_OUT}/pgs002288_by_chr"
os.makedirs("/tmp/gw", exist_ok=True)

subprocess.run(["gsutil", "-q", "cp", f"{GW_OUT}/expected_chromosomes.txt", "/tmp/gw/"], check=True)
chroms = [l.strip() for l in open("/tmp/gw/expected_chromosomes.txt") if l.strip()]
subprocess.run(["gsutil", "-q", "-m", "cp", f"{GW_OUT}/chr*.tsv.bgz", "/tmp/gw/"], check=True)
missing = [c for c in chroms if not os.path.exists(f"/tmp/gw/chr{c}.tsv.bgz")]
assert not missing, f"chromosomes not finished yet: {missing} - run PRS-G3 again first"

total = None
for c in chroms:
    s = pd.read_csv(f"/tmp/gw/chr{c}.tsv.bgz", sep="\t", compression="gzip", dtype={"person_id": str})
    s = s.set_index("person_id").iloc[:, 0].astype(float)          # first data column = that chromosome's score
    total = s if total is None else total.add(s, fill_value=0)
gw = -total.rename("prs_gw")      s
print(f"genome-wide score built for {len(gw):,} women from {len(chroms)} chromosomes (sign corrected)")

# cohorts, the 29-SNP score for comparison, age and ancestry
for f in ["early_pop.csv", "later_pop.csv", "controls.csv", "ctd_list.csv", "prs_scores.csv"]:
    subprocess.run(["gsutil", "-q", "cp", f"{PRS_OUT}/{f}", f"/tmp/gw/{f}"], check=True)
subprocess.run(["gsutil", "-u", P, "-q", "cp",
                "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv",
                "/tmp/gw/ancestry_preds.tsv"], check=True)
coh = {n: pd.read_csv(f"/tmp/gw/{n}.csv", dtype=str) for n in ["early_pop", "later_pop", "controls", "ctd_list"]}
early, later = set(coh["early_pop"]["person_id"]), set(coh["later_pop"]["person_id"])
controls = set(coh["controls"]["person_id"]) - early - later
later = later - early
ctd = set(coh["ctd_list"]["person_id"])
people = pd.concat([coh[n][["person_id", "birth_datetime"]] for n in ["early_pop", "later_pop", "controls"]]) \
           .drop_duplicates("person_id").set_index("person_id")
age = (pd.Timestamp("2026-10-04", tz="UTC") - pd.to_datetime(people["birth_datetime"], utc=True)).dt.days / 365.25
anc = pd.read_csv("/tmp/gw/ancestry_preds.tsv", sep="\t", dtype=str)
pcs = anc["pca_features"].str.strip("[]").str.split(",", expand=True).iloc[:, :10].astype(float)
pcs.columns = [f"PC{i}" for i in range(1, 11)]
pcs.index = anc["research_id"].astype(str)
label = anc.set_index("research_id")["ancestry_pred"]
prs29 = pd.read_csv("/tmp/gw/prs_scores.csv", dtype={"person_id": str}).set_index("person_id")["prs29"]

d = pd.concat([gw, prs29], axis=1).join(age.rename("age")).join(pcs).join(label.rename("ancestry"))
d["ctd"] = d.index.isin(ctd)
corr = d[["prs_gw", "prs29"]].corr().iloc[0, 1]
print(f"correlation between the 29-SNP and genome-wide scores: {corr:.2f}  (should now be about +0.59)")
assert corr > 0, "the genome-wide score still correlates negatively with the 29-SNP score - stop and check the sign"

def sup(n):
    return "<20" if 0 < n < 20 else f"{n:,}"

def fit(case_ids, ctrl_ids, score, subset, extra=""):
    x = d[d.index.isin(case_ids | ctrl_ids)].dropna(subset=["age", "PC1", "prs_gw", "prs29"])
    if subset == "eur":
        x = x[x["ancestry"] == "eur"]
    if extra == "no_ctd":
        x = x[~x["ctd"]]
    y = x.index.isin(case_ids).astype(int)
    ctrl = x[y == 0]
    z = (x[score] - ctrl[score].mean()) / ctrl[score].std()
    top = (x[score] >= ctrl[score].quantile(0.8)).astype(int)
    out = {"cases": int(y.sum()), "controls": int((y == 0).sum())}
    for name, v in [("per SD", z), ("top 20% vs rest", top)]:
        X = sm.add_constant(pd.concat([v.rename("score"), x[["age"] + [f"PC{i}" for i in range(1, 11)]]], axis=1))
        m = sm.Logit(y, X).fit(disp=0)
        b, se = m.params["score"], m.bse["score"]
        out[name] = (np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), m.pvalues["score"])
    return out

comparisons = {"Early POP vs no-POP controls": (early, controls),
               "Later POP vs no-POP controls": (later, controls),
               "Early POP vs later POP":       (early, later)}
analyses = [("Primary: European ancestry, genome-wide score", "prs_gw", "eur", ""),
            ("Same women, 29-SNP score (for comparison)",     "prs29",  "eur", ""),
            ("All ancestries, genome-wide score",             "prs_gw", "all", ""),
            ("European, CTD diagnoses removed, genome-wide",  "prs_gw", "eur", "no_ctd")]
rows = []
for a_name, score, subset, extra in analyses:
    for c_name, (ca, co) in comparisons.items():
        r = fit(ca, co, score, subset, extra)
        for est in ["per SD", "top 20% vs rest"]:
            o, lo, hi, p = r[est]
            rows.append({"analysis": a_name, "comparison": c_name, "estimate": est,
                         "cases": r["cases"], "controls": r["controls"], "OR": o, "lo": lo, "hi": hi, "p": p})
res = pd.DataFrame(rows)
res.to_csv("/tmp/gw/prs_gw_results.csv", index=False)
subprocess.run(["gsutil", "-q", "cp", "/tmp/gw/prs_gw_results.csv", f"{PRS_OUT}/prs_gw_results.csv"], check=True)

show = res.copy()
small = show["cases"].between(1, 19) | show["controls"].between(1, 19)
show["n"] = show["cases"].map(sup) + " / " + show["controls"].map(sup)
show["OR (95% CI)"] = np.where(small, "suppressed",
                               show.apply(lambda r: f"{r.OR:.2f} ({r.lo:.2f}-{r.hi:.2f})", axis=1))
show["p"] = np.where(small, "", show["p"].map(lambda v: "<0.001" if v < 0.001 else f"{v:.3f}"))
pd.set_option("display.width", 250, "display.max_colwidth", 60)
for a_name in show["analysis"].unique():
    print(f"\n=== {a_name}")
    print(show[show["analysis"] == a_name][["comparison", "estimate", "n", "OR (95% CI)", "p"]].to_string(index=False))
