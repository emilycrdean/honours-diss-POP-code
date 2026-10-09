# %% PRS-G5 - why is the genome-wide score pointing the wrong way? (no Hail needed, runs in a minute or two)
# essentially - this entire thing is looking at the PRS - it was not needed for the other analyses, exploratory info only
import os, subprocess
import numpy as np, pandas as pd

P = "YOUR-WORKSPACE-ID"
PRS_OUT = "gs://YOUR-WORKSPACE-ID/prs"
GW_OUT = f"{PRS_OUT}/pgs002288_by_chr"
LOCAL_PGS = "/tmp/PGS002288_hmPOS_GRCh38.txt.gz"
os.makedirs("/tmp/gw", exist_ok=True)

# ---------------------------------------------------------------------------------------------
# PART 1 - asking if the score's weights agree in sign with the published odds ratios at the 30 lead SNPs?

LEAD = pd.DataFrame([
    ("rs3820282", 1, 22141722, "T", 0.86), ("rs9306894", 2, 20678345, "G", 1.10),
    ("rs11899888", 2, 55875609, "G", 1.11), ("rs3791675", 2, 55884174, "T", 0.92),
    ("rs77648136", 2, 176056916, "T", 0.93), ("rs58170120", 3, 128002490, "A", 1.08),
    ("rs201194999", 4, 65801177, "T", 0.89), ("rs28403275", 4, 126003529, "C", 1.12),
    ("rs10013769", 4, 126155033, "G", 1.07), ("rs42400", 5, 4978582, "G", 0.94),
    ("rs251217", 5, 128176372, "G", 1.06), ("rs72624976", 7, 128392779, "T", 0.84),
    ("rs1493202", 8, 70993352, "G", 1.05), ("rs10810888", 9, 18023979, "G", 1.05),
    ("rs10762631", 10, 74609065, "A", 0.92), ("rs6484161", 11, 10287444, "T", 1.06),
    ("rs35166569", 11, 32324851, "C", 0.93), ("rs11031796", 11, 32458261, "A", 0.93),
    ("rs4944936", 11, 74681469, "C", 0.93), ("rs12314243", 12, 12515476, "T", 1.09),
    ("rs73197353", 12, 114232191, "C", 1.12), ("rs1247943", 12, 114235616, "A", 0.91),
    ("rs4779517", 15, 31344221, "G", 1.07), ("rs12915554", 15, 32731285, "A", 0.95),
    ("rs4886778", 15, 73933047, "A", 1.05), ("rs12325192", 16, 51454218, "T", 0.89),
    ("rs1874008", 16, 84906873, "C", 0.94), ("rs72839768", 17, 7226579, "A", 1.19),
    ("rs235929", 21, 27264365, "C", 0.93), ("rs2267372", 22, 38202227, "G", 0.93),
], columns=["rsid", "chrom", "pos", "pub_effect", "pub_or"])
LEAD["pub_logor"] = np.log(LEAD["pub_or"])

if not os.path.exists(LOCAL_PGS):
    subprocess.run(["gsutil", "-q", "cp", f"{PRS_OUT}/PGS002288_hmPOS_GRCh38.txt.gz", LOCAL_PGS], check=True)
pgs = pd.read_csv(LOCAL_PGS, sep="\t", comment="#", dtype=str)
print("score file columns:", list(pgs.columns))
other = "other_allele" if "other_allele" in pgs.columns else "hm_inferOtherAllele"
pgs["w"] = pgs["effect_weight"].astype(float)
print(f"{len(pgs):,} rows; weight sign: {(pgs['w'] > 0).mean():.1%} positive; "
      f"mean weight {pgs['w'].mean():.2e}")

rows = []
for r in LEAD.itertuples():
    hit = pgs[(pgs["hm_chr"] == str(r.chrom)) & (pgs["hm_pos"] == str(r.pos))]
    if hit.empty:
        rows.append({"rsid": r.rsid, "pub_effect": r.pub_effect, "pub_logor": round(r.pub_logor, 3),
                     "in_score": False, "pgs_effect": "", "pgs_other": "", "pgs_w": np.nan, "verdict": "not in score"})
        continue
    h = hit.iloc[0]
    if h["effect_allele"] == r.pub_effect:
        verdict = "agree" if np.sign(h["w"]) == np.sign(r.pub_logor) else "OPPOSITE SIGN"
    elif h[other] == r.pub_effect:
        verdict = "agree (alleles swapped)" if np.sign(h["w"]) == -np.sign(r.pub_logor) else "OPPOSITE SIGN (alleles swapped)"
    else:
        verdict = "alleles do not match (strand?)"
    rows.append({"rsid": r.rsid, "pub_effect": r.pub_effect, "pub_logor": round(r.pub_logor, 3), "in_score": True,
                 "pgs_effect": h["effect_allele"], "pgs_other": h[other], "pgs_w": float(f"{h['w']:.3g}"),
                 "verdict": verdict})
chk = pd.DataFrame(rows)
pd.set_option("display.width", 250, "display.max_rows", 200)
print("\n=== lead SNPs: published direction vs the score's weights")
print(chk.to_string(index=False))
print("\nverdict counts:\n", chk["verdict"].value_counts().to_string())

# ---------------------------------------------------------------------------------------------
# PART 2 - which chromosomes carry the inverted signal?  Standardised difference in each chromosome's
# partial score, later POP cases vs no-POP controls, European ancestry. One number per chromosome.
for f in ["later_pop.csv", "controls.csv", "early_pop.csv"]:
    if not os.path.exists(f"/tmp/gw/{f}"):
        subprocess.run(["gsutil", "-q", "cp", f"{PRS_OUT}/{f}", f"/tmp/gw/{f}"], check=True)
if not os.path.exists("/tmp/gw/ancestry_preds.tsv"):
    subprocess.run(["gsutil", "-u", P, "-q", "cp",
                    "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv",
                    "/tmp/gw/ancestry_preds.tsv"], check=True)
subprocess.run(["gsutil", "-q", "cp", f"{GW_OUT}/expected_chromosomes.txt", "/tmp/gw/"], check=True)
subprocess.run(["gsutil", "-q", "-m", "cp", f"{GW_OUT}/chr*.tsv.bgz", "/tmp/gw/"], check=True)
chroms = [l.strip() for l in open("/tmp/gw/expected_chromosomes.txt") if l.strip()]

coh = {n: pd.read_csv(f"/tmp/gw/{n}.csv", dtype=str, usecols=["person_id"]) for n in ["early_pop", "later_pop", "controls"]}
early, later = set(coh["early_pop"]["person_id"]), set(coh["later_pop"]["person_id"])
controls = set(coh["controls"]["person_id"]) - early - later
later = later - early
anc = pd.read_csv("/tmp/gw/ancestry_preds.tsv", sep="\t", dtype=str)
eur = set(anc.loc[anc["ancestry_pred"] == "eur", "research_id"].astype(str))

out = []
for c in chroms:
    s = pd.read_csv(f"/tmp/gw/chr{c}.tsv.bgz", sep="\t", compression="gzip", dtype={"person_id": str})
    s = s.set_index("person_id").iloc[:, 0].astype(float)          # first data column = that chromosome's score (chr22 pilot names it differently)
    s = s[s.index.isin(eur)]
    ctrl = s[s.index.isin(controls)]
    case = s[s.index.isin(later)]
    out.append({"chr": c, "std_diff_later_vs_controls": round((case.mean() - ctrl.mean()) / ctrl.std(), 3),
                "sd_in_controls": round(ctrl.std(), 4)})
res = pd.DataFrame(out)
print("\n=== standardised difference in each chromosome's partial score (later POP minus controls)")
print(res.to_string(index=False))
print(f"\nchromosomes with a negative difference: {(res['std_diff_later_vs_controls'] < 0).sum()} of {len(res)}")
