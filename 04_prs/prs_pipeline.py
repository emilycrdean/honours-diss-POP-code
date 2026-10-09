# %% PRS-1 - pull the 30 Pujol-Gualdo lead SNPs for every woman in the PRS cohorts (Hail)
# Needs the same Hail (Dataproc) setup as the cell that built the PVA. Saves one genotype file to the bucket 
# that way the big thing only has to run once and you can refer to it again. 
import io, subprocess
import pandas as pd
import hail as hl

P = "YOUR-WORKSPACE-ID"
PRS_OUT = "gs://vpa-full-YOUR-WORKSPACE-ID/prs"
MT_PATH = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt"

# 30 lead variants, Pujol-Gualdo et al. 2022 (GWAS Catalog GCST90102470), GRCh38 positions
SNPS = pd.read_csv(io.StringIO("""rsid,chrom,pos,effect_allele,eaf,or
rs3820282,1,22141722,T,0.16,0.86
rs9306894,2,20678345,G,0.38,1.10
rs11899888,2,55875609,G,0.11,1.11
rs3791675,2,55884174,T,0.25,0.92
rs77648136,2,176056916,T,0.15,0.93
rs58170120,3,128002490,A,0.18,1.08
rs201194999,4,65801177,T,0.30,0.89
rs28403275,4,126003529,C,0.18,1.12
rs10013769,4,126155033,G,0.65,1.07
rs42400,5,4978582,G,0.36,0.94
rs251217,5,128176372,G,0.61,1.06
rs72624976,7,128392779,T,0.01,0.84
rs1493202,8,70993352,G,0.52,1.05
rs10810888,9,18023979,G,0.65,1.05
rs10762631,10,74609065,A,0.09,0.92
rs6484161,11,10287444,T,0.31,1.06
rs35166569,11,32324851,C,0.09,0.93
rs11031796,11,32458261,A,0.31,0.93
rs4944936,11,74681469,C,0.72,0.93
rs12314243,12,12515476,T,0.12,1.09
rs73197353,12,114232191,C,0.08,1.12
rs1247943,12,114235616,A,0.54,0.91
rs4779517,15,31344221,G,0.49,1.07
rs12915554,15,32731285,A,0.32,0.95
rs4886778,15,73933047,A,0.46,1.05
rs12325192,16,51454218,T,0.18,0.89
rs1874008,16,84906873,C,0.77,0.94
rs72839768,17,7226579,A,0.02,1.19
rs235929,21,27264365,C,0.39,0.93
rs2267372,22,38202227,G,0.61,0.93"""))
assert len(SNPS) == 30

# everyone in the three main PRS cohorts (saved earlier with save_cohort)
ids = set()
for name in ["early_pop", "later_pop", "controls"]:
    subprocess.run(["gsutil", "-q", "cp", f"{PRS_OUT}/{name}.csv", f"/tmp/prs_{name}.csv"], check=True)
    d = pd.read_csv(f"/tmp/prs_{name}.csv", dtype=str, usecols=["person_id"])
    print(f"{name}: {len(d):,}")
    ids |= set(d["person_id"])
print(f"women to genotype: {len(ids):,}")

hl.init(default_reference="GRCh38", idempotent=True, gcs_requester_pays_configuration=P)
mt = hl.read_matrix_table(MT_PATH)
print("row key:", list(mt.row_key), "| column key:", list(mt.col_key), "| entry fields:", list(mt.entry))
assert "GT" in mt.entry, "no GT field"
sample = list(mt.col_key)[0]

# a 3-bp window round each position because this is genetics so there's errors and indels
intervals = [hl.parse_locus_interval(f"chr{c}:{p - 1}-{p + 2}", reference_genome="GRCh38")
             for c, p in zip(SNPS["chrom"], SNPS["pos"])]
mt = hl.filter_intervals(mt, intervals)
mt = mt.filter_cols(hl.literal(ids).contains(mt[sample]))
mt = mt.select_rows()
mt = mt.select_entries(d=mt.GT.n_alt_alleles())

found = mt.rows().collect()
print(f"\n{len(found)} variant rows found at the 30 positions:")
for r in found:
    print(f"  {r.locus}  {r.alleles}")
print(f"women in the genotype file: {mt.count_cols():,}")

ent = mt.entries().key_by()
ent = ent.select(person_id=hl.str(ent[sample]),
                 variant=hl.str(ent.locus.contig) + ":" + hl.str(ent.locus.position) + ":" + hl.delimit(ent.alleles, ":"),
                 dosage=ent.d)
ent.export(f"{PRS_OUT}/genotypes_30snps.tsv.bgz")
print("saved", f"{PRS_OUT}/genotypes_30snps.tsv.bgz")


# %% match each SNP to its effect allele, check the genotypes, and build the score
import numpy as np
LOCAL_G = "/tmp/genotypes_30snps.tsv.bgz"                       # /tmp, never the notebook folder
subprocess.run(["gsutil", "-q", "cp", f"{PRS_OUT}/genotypes_30snps.tsv.bgz", LOCAL_G], check=True)

# 1. pick, for each lead SNP, the row at its exact position whose alleles contain the effect allele
variants = pd.DataFrame([{"variant": f"{r.locus.contig}:{r.locus.position}:{':'.join(r.alleles)}",
                          "chrom": r.locus.contig.replace("chr", ""), "pos": r.locus.position,
                          "ref": r.alleles[0], "alt": r.alleles[1]} for r in found])
match = []
for s in SNPS.to_dict("records"):
    rows = variants[(variants["chrom"] == str(s["chrom"])) & (variants["pos"] == s["pos"]) &
                    (variants["ref"].str.len() == 1) & (variants["alt"].str.len() == 1) &
                    ((variants["alt"] == s["effect_allele"]) | (variants["ref"] == s["effect_allele"]))]
    alt_hits = rows[rows["alt"] == s["effect_allele"]]
    pick = alt_hits.iloc[0] if len(alt_hits) == 1 else (rows.iloc[0] if len(rows) == 1 else None)
    match.append({"rsid": s["rsid"], "effect_allele": s["effect_allele"], "reported_eaf": s["eaf"], "or": s["or"],
                  "variant": None if pick is None else pick["variant"],
                  "effect_is": None if pick is None else ("alt" if pick["alt"] == s["effect_allele"] else "ref"),
                  "palindromic": None if pick is None else {pick["ref"], pick["alt"]} in ({"A", "T"}, {"C", "G"}),
                  "rows_at_position": len(rows)})
match = pd.DataFrame(match)
assert match["variant"].notna().all(), f"unmatched: {match.loc[match['variant'].isna(), 'rsid'].tolist()}"

# 2. read only the 30 chosen variants (missing genotypes come out of Hail as NA, ignore it) 
keep = set(match["variant"])
parts = []
for ch in pd.read_csv(LOCAL_G, sep="\t", compression="gzip", dtype={"person_id": str, "variant": str},
                      chunksize=2_000_000):
    parts.append(ch[ch["variant"].isin(keep)])
g = pd.concat(parts).pivot(index="person_id", columns="variant", values="dosage").astype(float)
print(f"genotype matrix: {g.shape[0]:,} women x {g.shape[1]} variants")

# 3. dosage of the EFFECT allele, and checks per SNP (frequencies only - no participant counts - if this is all of us, do not print participant counts)
eff = pd.DataFrame(index=g.index)
for m in match.itertuples():
    eff[m.rsid] = g[m.variant] if m.effect_is == "alt" else 2 - g[m.variant]

subprocess.run(["gsutil", "-u", P, "cp",
                "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv",
                "/tmp/ancestry_preds.tsv"], check=True)
anc = pd.read_csv("/tmp/ancestry_preds.tsv", sep="\t", dtype=str)
print("ancestry file columns:", list(anc.columns))
assert "ancestry_pred" in anc.columns, f"no ancestry_pred column: {list(anc.columns)}"
label = anc.set_index("research_id")["ancestry_pred"]
eur = eff.index[eff.index.map(lambda i: label.get(i, "") == "eur")]

match["missing_pct"] = [round(100 * eff[r].isna().mean(), 2) for r in match["rsid"]]
match["aou_eaf_all"] = [round(eff[r].mean() / 2, 3) for r in match["rsid"]]
match["aou_eaf_eur"] = [round(eff.loc[eur, r].mean() / 2, 3) for r in match["rsid"]]
pd.set_option("display.width", 250)
print(match.drop(columns=["variant"]).to_string(index=False))
assert match["missing_pct"].max() < 5, "a SNP is missing in over 5% of women"

# 4. weighted score: sum of log(OR) x effect-allele dosage --- missing genotypes set to that SNP's mean
w = np.log(match.set_index("rsid")["or"])
filled = eff.fillna(eff.mean())
prs = pd.DataFrame({"person_id": eff.index, "prs_raw": filled[w.index].values @ w.values,
                    "n_missing_snps": eff.isna().sum(axis=1).values})
prs.to_csv("/tmp/prs_scores.csv", index=False)
match.to_csv("/tmp/prs_snp_match.csv", index=False)
for f in ["prs_scores.csv", "prs_snp_match.csv"]:
    subprocess.run(["gsutil", "-q", "cp", f"/tmp/{f}", f"{PRS_OUT}/{f}"], check=True)
print(f"\nscore built for {len(prs):,} women; saved prs_scores.csv and prs_snp_match.csv")


# %% PRS-2d - build the score from 29 SNPs (rs201194999 dropped: not called in most genomes)
DROP = ["rs201194999"]
use = match[~match["rsid"].isin(DROP)]
w = np.log(use.set_index("rsid")["or"])
filled = eff[w.index].fillna(eff[w.index].mean())                 # the rare missing calls get the SNP's mean
no_amb = [r for r in w.index if r != "rs4779517"]                  # sensitivity: drop the ambiguous palindromic SNP
prs = pd.DataFrame({"person_id": eff.index,
                    "prs29": filled.values @ w.values,
                    "prs28_no_rs4779517": filled[no_amb].values @ w[no_amb].values})
prs.to_csv("/tmp/prs_scores.csv", index=False)
match.to_csv("/tmp/prs_snp_match.csv", index=False)
for f in ["prs_scores.csv", "prs_snp_match.csv"]:
    subprocess.run(["gsutil", "-q", "cp", f"/tmp/{f}", f"{PRS_OUT}/{f}"], check=True)
print(f"{len(w)} SNPs in the score; built for {len(prs):,} women; saved prs_scores.csv")


# %% PRS-3 - does the score predict early-onset POP? (logistic regression, OR per SD of the score)
import statsmodels.api as sm

def sup(n):                       # counts 1-19 are never shown
    return "<20" if 0 < n < 20 else f"{n:,}"

subprocess.run(["gsutil", "-q", "cp", f"{PRS_OUT}/ctd_list.csv", "/tmp/prs_ctd_list.csv"], check=True)
coh = {n: pd.read_csv(f"/tmp/prs_{n}.csv", dtype=str) for n in ["early_pop", "later_pop", "controls", "ctd_list"]}
early, later = set(coh["early_pop"]["person_id"]), set(coh["later_pop"]["person_id"])
controls = set(coh["controls"]["person_id"]) - early - later      # no-POP controls, any overlap removed
later = later - early
ctd = set(coh["ctd_list"]["person_id"])

people = pd.concat([coh[n][["person_id", "birth_datetime"]] for n in ["early_pop", "later_pop", "controls"]]) \
           .drop_duplicates("person_id").set_index("person_id")
age = (pd.Timestamp("2026-10-04", tz="UTC") - pd.to_datetime(people["birth_datetime"], utc=True)).dt.days / 365.25

pcs = anc["pca_features"].str.strip("[]").str.split(",", expand=True).iloc[:, :10].astype(float)
pcs.columns = [f"PC{i}" for i in range(1, 11)]
pcs.index = anc["research_id"].astype(str)

d = prs.set_index("person_id").join(age.rename("age")).join(pcs).join(label.rename("ancestry"))
d["ctd"] = d.index.isin(ctd)

def fit(case_ids, ctrl_ids, score, subset, extra=""):
    x = d[d.index.isin(case_ids | ctrl_ids)].dropna(subset=["age", "PC1", score])
    if subset == "eur":
        x = x[x["ancestry"] == "eur"]
    if extra == "no_ctd":
        x = x[~x["ctd"]]
    y = x.index.isin(case_ids).astype(int)
    ctrl = x[y == 0]
    z = (x[score] - ctrl[score].mean()) / ctrl[score].std()          # SD units, scaled on the controls
    top = (x[score] >= ctrl[score].quantile(0.8)).astype(int)        # top 20% of the control distribution
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
analyses = [("Primary: European ancestry, 29 SNPs",     "prs29",              "eur", ""),
            ("All ancestries, 29 SNPs",                 "prs29",              "all", ""),
            ("European, without rs4779517",             "prs28_no_rs4779517", "eur", ""),
            ("European, CTD diagnoses removed",         "prs29",              "eur", "no_ctd")]
rows = []
for a_name, score, subset, extra in analyses:
    for c_name, (ca, co) in comparisons.items():
        r = fit(ca, co, score, subset, extra)
        for est in ["per SD", "top 20% vs rest"]:
            o, lo, hi, p = r[est]
            rows.append({"analysis": a_name, "comparison": c_name, "estimate": est,
                         "cases": r["cases"], "controls": r["controls"], "OR": o, "lo": lo, "hi": hi, "p": p})
res = pd.DataFrame(rows)
res.to_csv("/tmp/prs_results.csv", index=False)
subprocess.run(["gsutil", "-q", "cp", "/tmp/prs_results.csv", f"{PRS_OUT}/prs_results.csv"], check=True)

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


# %% PRS-G1 - get the genome-wide Pujol-Gualdo score (PGS002288, GRCh38 version) into the bucket
import re, time
PGS_DIR = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002288/ScoringFiles/Harmonized/"
idx = subprocess.run(["curl", "-sSL", "-f", "--max-time", "60", PGS_DIR], capture_output=True, text=True)
assert idx.returncode == 0, f"The Workbench could not reach the PGS Catalog ({idx.stderr[:200]})"
names = sorted({h.split("/")[-1] for h in re.findall(r'href="([^"]+)"', idx.stdout) if "PGS002288" in h})
print("files in the folder:", names)
target = [n for n in names if "GRCh38" in n and n.endswith(".txt.gz")]
assert len(target) == 1, f"expected one GRCh38 file, found {target}"
LOCAL_PGS = f"/tmp/{target[0]}"                                    # /tmp, never the notebook folder
subprocess.run(["curl", "-sSL", "-f", "-o", LOCAL_PGS, PGS_DIR + target[0]], check=True)
subprocess.run(["gsutil", "-q", "cp", LOCAL_PGS, f"{PRS_OUT}/{target[0]}"], check=True)
print("saved", f"{PRS_OUT}/{target[0]}")
print(subprocess.run(f"zcat {LOCAL_PGS} | grep '^#' | head -30", shell=True, capture_output=True, text=True).stdout)


# %% PRS-G2 - pilot: score chromosome 22 only, time it, and estimate the full genome-wide run
pgs = pd.read_csv(LOCAL_PGS, sep="\t", comment="#", dtype=str)
print("columns:", list(pgs.columns))
need = {"hm_chr", "hm_pos", "effect_allele", "effect_weight"}
assert need <= set(pgs.columns), f"missing {need - set(pgs.columns)}"
other = "other_allele" if "other_allele" in pgs.columns else "hm_inferOtherAllele"
assert other in pgs.columns, "no other-allele column"
pgs = pgs.dropna(subset=["hm_chr", "hm_pos", "effect_allele", other, "effect_weight"])
pgs = pgs[pgs[other].str.fullmatch(r"[ACGT]+") & pgs["effect_allele"].str.fullmatch(r"[ACGT]+")]
pgs["effect_weight"] = pgs["effect_weight"].astype(float)
print(f"{len(pgs):,} weighted variants with a GRCh38 position")

w22 = pgs[pgs["hm_chr"] == "22"]
both = pd.concat([                                                  # both orientations, so ref/alt order never matters
    pd.DataFrame({"pos": w22["hm_pos"].astype(int).values, "a0": w22[other].values,
                  "a1": w22["effect_allele"].values, "w": w22["effect_weight"].values, "flip": False}),
    pd.DataFrame({"pos": w22["hm_pos"].astype(int).values, "a0": w22["effect_allele"].values,
                  "a1": w22[other].values, "w": w22["effect_weight"].values, "flip": True})])
wt = hl.Table.from_pandas(both)
wt = wt.annotate(locus=hl.locus("chr22", wt.pos, reference_genome="GRCh38"), alleles=hl.array([wt.a0, wt.a1]))
wt = wt.key_by("locus", "alleles").select("w", "flip")

mt = hl.read_matrix_table(MT_PATH)
sample = list(mt.col_key)[0]
mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr22", reference_genome="GRCh38")])
mt = mt.filter_cols(hl.literal(ids).contains(mt[sample]))
mt = mt.annotate_rows(pgs=wt[mt.row_key])
mt = mt.filter_rows(hl.is_defined(mt.pgs))
matched = mt.count_rows()
print(f"chr22: {len(w22):,} weighted variants, {matched:,} found in the genotype file ({100 * matched / len(w22):.1f}%)")

mt = mt.select_entries(d=hl.if_else(mt.pgs.flip, 2 - mt.GT.n_alt_alleles(), mt.GT.n_alt_alleles()))
mt = mt.annotate_rows(mean_d=hl.agg.mean(mt.d))                    # missing calls get the variant's mean
mt = mt.annotate_cols(score=hl.agg.sum(hl.or_else(mt.d, mt.mean_d) * mt.pgs.w))

t0 = time.time()
cols = mt.cols().key_by()
cols = cols.select(person_id=hl.str(cols[sample]), score_chr22=cols.score)
cols.export(f"{PRS_OUT}/pgs002288_chr22.tsv.bgz")
minutes = (time.time() - t0) / 60

workers = hl.spark_context()._jsc.sc().getExecutorMemoryStatus().size() - 1
print(f"\nchr22 took {minutes:.1f} min on this cluster ({workers} executors, {hl.spark_context().defaultParallelism} cores)")
print(f"chr22 holds {100 * len(w22) / len(pgs):.1f}% of the score's variants")
print(f"rough estimate for the whole genome: {minutes * len(pgs) / len(w22) / 60:.1f} hours on this cluster")
