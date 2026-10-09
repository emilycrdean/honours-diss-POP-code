# =====================================================================================
# Rare pathogenic variant analysis (Aim 2) - full re-run pipeline for a new workspace
# Python notebook. Each "# %%" block is one notebook cell. Run them in order.
#
# Lessons from the September run are built in:
#   - read every extracted file with dtype=str (the ClinVar star column silently became NaN once)
#   - star rating = maximum across submitters (the field is comma-separated)
#   - P/LP is judged per submitter call, so "clean" and "conflicting" are separated explicitly
#   - every filter step is reported per gene, so genes that drop out (ESR1, BPNT2, RMP64) are visible
#   - BigQuery person_id comes back as INT64; everything is cast to str before set operations
#   - case and comparison cohorts are checked for overlap before any test
#   - every printed participant count of 1-19 is replaced with "<20", so output is safe to copy out
# =====================================================================================

# %% CELL 0a - recover workspace settings the kernel cannot see (run before anything else)
import os, subprocess
P = os.environ["GOOGLE_PROJECT"]

# 1. The Jupyter server that started this kernel usually still holds them
found = {}
for pid in {os.getppid(), 1}:
    try:
        for kv in open(f"/proc/{pid}/environ", "rb").read().split(b"\0"):
            k, _, val = kv.decode(errors="ignore").partition("=")
            if any(s in k for s in ["CDR", "WORKSPACE", "BUCKET", "STORAGE", "WGS"]) or val.startswith("gs://"):
                found[k] = val
    except Exception as e:
        print(f"could not read process {pid}: {e}")
for k, val in found.items():
    os.environ.setdefault(k, val)
os.environ.setdefault("WORKSPACE_CDR", "YOUR-CDR-PROJECT.C2025Q4R6")   # so Dataset Builder code also runs
print("Recovered from the Jupyter server:", found or "nothing")

# 2. If still missing, list where they could be
def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print((r.stdout or r.stderr).strip()[:4000])
if "WORKSPACE_BUCKET" not in os.environ:
    print("\n--- buckets in this workspace's project:"); sh(f"gsutil ls -p {P}")
if "CDR_STORAGE_PATH" not in os.environ:
    print("\n--- All of Us controlled-tier releases:"); sh(f"gsutil -u {P} ls gs://fc-aou-datasets-controlled/")
for k in ["WORKSPACE_CDR", "WORKSPACE_BUCKET", "CDR_STORAGE_PATH"]:
    print(f"{k}: {os.environ.get(k, '(still missing)')}")

# %% CELL 0 - discovery (run first, then fill in CELL 1 from its output)
import os, re, subprocess
from google.cloud import bigquery

client = bigquery.Client()
print("Workspace variables this kernel can see:")
for k in sorted(os.environ):
    if any(s in k for s in ["CDR", "WORKSPACE", "BUCKET", "GOOGLE_PROJECT"]):
        print(f"  {k} = {os.environ[k]}")

CDR = os.environ.get("WORKSPACE_CDR") or "YOUR-CDR-PROJECT.C2025Q4R6"   # C2025Q4R6 CDR
STORAGE = os.environ.get("CDR_STORAGE_PATH", "")
BUCKET = os.environ.get("WORKSPACE_BUCKET", "")
print("\nCDR:", CDR, "\nSTORAGE:", STORAGE or "(not set)", "\nBUCKET:", BUCKET or "(not set)")
print("CDR reachable:", client.get_dataset(CDR).dataset_id)

tables = [t.table_id for t in client.list_tables(CDR)]
for t in [t for t in tables if re.search(r"(^|_)vat(_|$)|variant|clinvar", t.lower())]:   # not "obserVATion"
    print("\nBQ table:", t, [f.name for f in client.get_table(f"{CDR}.{t}").schema])

if STORAGE:
    ls = subprocess.run(["gsutil", "ls", "-r", f"{STORAGE}/wgs/short_read/snpindel/aux/**"],
                        capture_output=True, text=True).stdout.splitlines()
    print("\nVAT / ancestry files:")
    print("\n".join(l for l in ls if "vat" in l.lower() or "ancestry" in l.lower()))

# %% CELL 0d - find the VAT / ancestry paths used by earlier notebooks, and any non-bucket resources
import json
def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return (r.stdout or r.stderr).strip()

print("--- workspace resources that are not your own buckets")
print(sh("wb resource list | grep -v GCS_BUCKET") or "(none)")

print("\n--- data paths mentioned in your notebooks")
P = os.environ["GOOGLE_PROJECT"]
nbs = [l for l in sh(f"gsutil ls 'gs://*-{P}/**.ipynb'").splitlines() if l.startswith("gs://") and "checkpoint" not in l]
print(len(nbs), "notebooks searched")
paths = {}
for nb in nbs:
    text = sh(f"gsutil cat '{nb}'")
    for m in re.findall(r"""(gs://[^\s"'\\,)]+|[\w-]+\.C2025Q4R6[\w.]*)""", text):
        if any(k in m.lower() for k in ["vat", "ancestry", "wgs", "aux", "index"]):
            paths.setdefault(m, set()).add(nb.split("/")[2])
for m, where in sorted(paths.items()):
    print(f"  {m}\n      in: {', '.join(sorted(where))}")

# %% CELL 1 - configuration
import os, subprocess
VAT_SOURCE = "hail"     # no VAT table in BigQuery (CELL 0), so the VAT file is read with Hail
VAT_PATH   = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"
VAT_TABLE  = ""                                    # bigquery route only - not used
ANCESTRY_PATH = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
P = os.environ["GOOGLE_PROJECT"]                  # billing project for the All of Us data bucket
BUCKET = "gs://vpa-full-YOUR-WORKSPACE-ID"    # results go in a vpa_rerun/ folder, apart from the export

COL = {                       # VAT column names (from the September extraction; confirm in CELL 0)
    "vid": "vid", "gene": "gene_symbol", "af": "gvs_all_af", "canonical": "is_canonical_transcript",
    "consequence": "consequence", "clinvar": "clinvar_classification", "stars": "clinvar_rcv_num_stars",
    "revel": "revel",         # optional - used only for the extended tier
}
OUT = f"{BUCKET}/vpa_rerun"

def save(df, name, **kw):      # write in /tmp, then copy to the bucket (never into the notebook folder)
    os.makedirs("/tmp/vpa", exist_ok=True)
    df.to_csv(f"/tmp/vpa/{name}", **kw)
    subprocess.run(["gsutil", "-q", "cp", f"/tmp/vpa/{name}", f"{OUT}/{name}"], check=True)

CORE10 = "FBN1 FBLN5 COL1A1 COL1A2 COL3A1 COL5A1 COL5A2 TGFBR1 TGFBR2 TNXB".split()
WIDE147 = """ABCC6 ABL1 ACTA2 ACVR1 ADAMTS10 ADAMTS17 ADAMTS2 ADAMTSL2 ADAMTSL4 AEBP1 ALDH18A1 ALPL ARIH1 ASPH ATP6AP1
ATP6V0A2 ATP6V1A ATP6V1E1 ATP7A B3GALT6 B3GAT3 B4GALT7 BGN BPNT2 C1R C1S CBS CHD4 CHST14 CHST3 COG7 COL11A1 COL11A2
COL12A1 COL1A1 COL1A2 COL2A1 COL3A1 COL4A1 COL5A1 COL5A2 COL6A2 COL6A3 COL9A1 COL9A2 COL9A3 CREB3L1 CRTAP CSGALNACT1
DCHS1 DLG4 DSE EFEMP2 ELN ERBIN FBLN5 FBN1 FBN2 FKBP14 FLCN FLNA FLNB FN1 FOXE3 GATA5 GGCX GORAB HCN4 IPO8 KIF22 KLK1
KLK10 KLK11 KLK12 KLK13 KLK14 KLK15 KLK2 KLK3 KLK4 KLK5 KLK6 KLK7 KLK8 LEMD3 LOX LOXL3 LRP2 LTBP2 LTBP3 LTBP4 LZTS1
MAT2A MBTPS2 MED12 MFAP5 MOCS1 MYH11 MYLK NEK9 NEPRO NKAP NOG NOTCH1 NT5E P3H1 P4HA1 PKD2 PLOD1 PLOD3 PLP1 PRDM5 PRKG1
PYCR1 RIN2 ROBO3 ROBO4 SERPING1 SGMS2 SHOC2 SKI SLC26A2 SLC2A10 SLC39A13 SMAD2 SMAD3 SMAD4 SMAD6 SMS SP7 SPARC TAB2
TALDO1 TGFB1 TGFB2 TGFB3 TGFBR1 TGFBR2 TLN1 TNFRSF1A TNXB UPF3B VCAN WNT1 ZDHHC9 ZNF469""".split()
ADDON24 = """ACADVL CHRDL2 CRISPLD2 DUSP16 DVL2 EFEMP1 ESR1 FAT4 GDF7 HOXA11 HOXD13 IMPDH1 KLF13 LAMC1 LDAH LOXL1
LOXL1-AS1 MAFF SALL1 TBX5 TIMP2 VCL WNT4 WT1""".split()
NEW = """COL6A1 SERPINH1 FKBP10 PPIB PLOD2 P4HA2 SERPINF1 BMP1 TMEM38B PLS3 IFITM5 LOXL2 LOXL4 MFAP4 MMP1 MMP2 MMP3 MMP9
TIMP1 TIMP3 IL1B IL6 NLRP3 ESR2 PGR HOXA13 DES""".split()
REMOVE = ("KLK1 KLK2 KLK3 KLK4 KLK5 KLK6 KLK7 KLK8 KLK10 KLK11 KLK12 KLK13 KLK14 DLG4 PLP1 ROBO3 "
          "ACADVL LOXL1-AS1 CHRDL2 CRISPLD2 DUSP16 DVL2 FAT4 GDF7 IMPDH1 KLF13 LAMC1 LDAH MAFF TBX5").split()
# The last 14 were included only because they lie at POP GWAS loci. GWAS signals come from common, mostly
# regulatory variants; rare pathogenic coding variants in the nearest gene test a different hypothesis.
# GWAS-locus genes with an independent mechanism (LOXL1, EFEMP1, TIMP2, VCL, ESR1, HOXA11, WNT4...) are kept.
NEG_CONTROL = "CFTR GJB2 PAH HEXA ATP7B SLC26A4 ACADM BTD GALT MMUT".split()

MECHANISMS = {   # every revised-panel gene in exactly one mechanism, assigned before the rerun results
    "Structural collagens and collagen-associated matrix":
        "AEBP1 COL11A1 COL11A2 COL12A1 COL1A1 COL1A2 COL2A1 COL3A1 COL4A1 COL5A1 COL5A2 COL6A1 COL6A2 COL6A3 COL9A1 COL9A2 COL9A3 FN1 TNXB".split(),
    "Collagen biosynthesis, modification and cross-linking":
        "ADAMTS2 ALDH18A1 ATP7A BMP1 CREB3L1 CRTAP FKBP10 FKBP14 IFITM5 LOX LOXL1 LOXL2 LOXL3 LOXL4 MBTPS2 P3H1 P4HA1 P4HA2 PLOD1 PLOD2 PLOD3 PLS3 PPIB PRDM5 PYCR1 SERPINF1 SERPINH1 SP7 SPARC TMEM38B WNT1 ZNF469".split(),
    "Proteoglycans and glycosaminoglycans":
        "B3GALT6 B3GAT3 B4GALT7 BGN BPNT2 CHST14 CHST3 CSGALNACT1 DSE SLC26A2 SLC39A13 VCAN".split(),
    "Elastic fibres and microfibrils":
        "ABCC6 ADAMTS10 ADAMTS17 ADAMTSL2 ADAMTSL4 ATP6AP1 ATP6V0A2 ATP6V1A ATP6V1E1 EFEMP1 EFEMP2 ELN FBLN5 FBN1 FBN2 GGCX GORAB LTBP2 LTBP3 LTBP4 MFAP4 MFAP5 NT5E RIN2 SLC2A10".split(),
    "TGF-beta and BMP signalling":
        "ACVR1 ERBIN IPO8 LEMD3 NOG NOTCH1 SKI SMAD2 SMAD3 SMAD4 SMAD6 TAB2 TGFB1 TGFB2 TGFB3 TGFBR1 TGFBR2".split(),
    "Smooth muscle, cytoskeleton and vessel wall":
        "ACTA2 ARIH1 DES FLNA FLNB FOXE3 GATA5 MAT2A MYH11 MYLK PKD2 PRKG1 ROBO4 TLN1 VCL".split(),
    "Matrix remodelling, complement and inflammation":
        "C1R C1S IL1B IL6 MMP1 MMP2 MMP3 MMP9 NLRP3 SERPING1 TIMP1 TIMP2 TIMP3 TNFRSF1A".split(),
    "Hormonal signalling and pelvic development":
        "ESR1 ESR2 HOXA11 HOXA13 HOXD13 PGR SALL1 WNT4 WT1".split(),
    "Syndromic or other, no direct connective tissue mechanism":
        "ABL1 ALPL ASPH CBS CHD4 COG7 DCHS1 FLCN HCN4 KIF22 KLK15 LRP2 LZTS1 MED12 MOCS1 NEK9 NEPRO NKAP SGMS2 SHOC2 SMS TALDO1 UPF3B ZDHHC9".split(),
}
PANELS = {
    "Core-10": CORE10,
    "Wide-147 (original)": WIDE147,
    "Wide-171 (original)": WIDE147 + ADDON24,
    "Revised wide": sorted((set(WIDE147) | set(ADDON24) | set(NEW)) - set(REMOVE)),
    **{f"Mechanism: {k}": v for k, v in MECHANISMS.items()},
    "Negative control": NEG_CONTROL,
}
ALL_GENES = sorted(set().union(*PANELS.values()))
print(len(ALL_GENES), "genes to extract;", {k: len(v) for k, v in PANELS.items()})

# %% CELL 2 - extract VAT rows for all genes (one-off; saves to the workspace bucket)
EXTRACT = f"{OUT}/vat_extract.tsv"
if VAT_SOURCE == "bigquery":
    job = client.query(f"SELECT * FROM `{VAT_TABLE}` WHERE {COL['gene']} IN UNNEST(@genes)",
                       job_config=bigquery.QueryJobConfig(query_parameters=[
                           bigquery.ArrayQueryParameter("genes", "STRING", ALL_GENES)]))
    save(job.to_dataframe().astype(str), "vat_extract.tsv", sep="\t", index=False)
else:
    # Needs a Dataproc (Hail) cluster. Restart any other kernel using the cluster first -
    # last time an abandoned job blocked it.
    import hail as hl
    hl.init(default_reference="GRCh38", idempotent=True, gcs_requester_pays_configuration=P)
    ht = hl.import_table(VAT_PATH, force_bgz=True, impute=False)       # every column as string
    ht = ht.filter(hl.literal(set(ALL_GENES)).contains(ht[COL["gene"]]))
    ht.export(EXTRACT)
print("saved", EXTRACT)

# %% CELL 3 - qualifying variants per tier, with a per-gene audit of every filter step
# Reads the extract in chunks and keeps only rows that could qualify for any tier (canonical, AF<1%,
# and either a ClinVar entry or a protein-altering consequence). Every tier is defined within that
# subset, so the results are identical to reading the whole file - just far faster and lighter.
import re
import numpy as np, pandas as pd

os.makedirs("/tmp/vpa", exist_ok=True)
LOCAL = "/tmp/vpa/vat_extract.tsv"
subprocess.run(["gsutil", "-q", "cp", EXTRACT, LOCAL], check=True)
print(f"extract: {os.path.getsize(LOCAL) / 1e9:.1f} GB")

header = pd.read_csv(LOCAL, sep="\t", nrows=0).columns
missing = [c for k, c in COL.items() if c not in header and k != "revel"]
assert not missing, f"Columns not found: {missing} - fix COL in CELL 1. Available: {list(header)}"
use = [c for c in COL.values() if c in header]

PLP = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
CODING = {"missense_variant", "stop_gained", "frameshift_variant", "splice_donor_variant",
          "splice_acceptor_variant", "start_lost", "stop_lost", "inframe_insertion", "inframe_deletion"}
LOF = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost"}
MAYBE = "missense_variant|stop_gained|frameshift_variant|splice_donor_variant|splice_acceptor_variant|start_lost"

g = COL["gene"]
rows_all = pd.Series(0, index=ALL_GENES, dtype=float)    # VAT rows per gene (all transcripts)
rows_can = pd.Series(0, index=ALL_GENES, dtype=float)    # canonical-transcript rows = variants
kept, n = [], 0
for ch in pd.read_csv(LOCAL, sep="\t", dtype=str, usecols=use, keep_default_na=False,
                      na_values=[""], chunksize=1_000_000):
    n += len(ch)
    rows_all = rows_all.add(ch[g].value_counts(), fill_value=0)
    can = ch[COL["canonical"]].str.lower().isin(["true", "1", "yes"])
    rows_can = rows_can.add(ch.loc[can, g].value_counts(), fill_value=0)
    rare = pd.to_numeric(ch[COL["af"]], errors="coerce") < 0.01
    maybe = ch[COL["clinvar"]].notna() | ch[COL["consequence"]].fillna("").str.contains(MAYBE)
    kept.append(ch[can & rare & maybe])
    print(f"  {n:,} rows read", end="\r")
v = pd.concat(kept, ignore_index=True)
print(f"\n{n:,} rows read; {len(v):,} kept as possible qualifiers")

def calls(x):         # individual submitter classifications, normalised
    return {c.strip().lower().replace("_", " ") for c in re.split(r"[,|;/]", str(x)) if c.strip()} - {"nan"}

def max_stars(x):
    s = [float(t) for t in re.split(r"[,|;]", str(x)) if t.strip() not in ("", "nan", "None")]
    return max(s) if s else np.nan

def csq(x):
    return {c.strip() for c in re.split(r"[,&|]", str(x)) if c.strip()}

v["calls"]    = v[COL["clinvar"]].map(calls)
v["plp"]      = v["calls"].map(lambda s: len(s & PLP) > 0)
v["clean"]    = v["calls"].map(lambda s: len(s) > 0 and s <= PLP)
v["stars"]    = v[COL["stars"]].map(max_stars)
v["csq"]      = v[COL["consequence"]].map(csq)
v["coding"]   = v["csq"].map(lambda s: len(s & CODING) > 0)
v["lof"]      = v["csq"].map(lambda s: len(s & LOF) > 0)
v["revel_num"] = pd.to_numeric(v[COL["revel"]], errors="coerce") if COL["revel"] in v.columns else np.nan
v["damaging_missense"] = v["csq"].map(lambda s: "missense_variant" in s) & (v["revel_num"] >= 0.75)

TIERS = {
    "Unrestricted":           v["plp"],
    "Primary (coding, >=1*)": v["plp"] & v["coding"] & (v["stars"] >= 1),
    "Strict (>=2*, clean)":   v["plp"] & v["coding"] & (v["stars"] >= 2) & v["clean"],
    "Extended (+LoF, REVEL)": (v["plp"] & v["coding"] & (v["stars"] >= 1)) | v["lof"] | v["damaging_missense"],
}

# per-gene audit (variant counts, not people - safe to share)
audit = pd.DataFrame({"VAT rows": rows_all.astype(int), "canonical variants": rows_can.astype(int)})
steps = {"rare + ClinVar/coding": pd.Series(True, index=v.index), "+ P/LP": v["plp"],
         "+ coding": v["plp"] & v["coding"], **TIERS}
for name, mask in steps.items():
    audit[name] = v[mask].groupby(g)[COL["vid"]].nunique().reindex(audit.index).fillna(0).astype(int)
audit.index.name = "gene"
save(audit, "gene_filter_audit.csv")
print(audit.to_string())
print("\nGenes with no VAT rows:", list(audit.index[audit["VAT rows"] == 0]))
print("Genes with 0 primary-tier variants:", list(audit.index[audit["Primary (coding, >=1*)"] == 0]))

QUAL = {t: v[m][[COL["vid"], g]].drop_duplicates().rename(columns={COL["vid"]: "vid", g: "gene"})
        for t, m in TIERS.items()}

# %% CELL 4 - cohorts, read from the Dataset Builder exports in each bucket
P = "YOUR-WORKSPACE-ID"
COHORT_BUCKETS = {
    "Cases (POP 18-35)":         f"gs://rerun-pva-case-{P}",
    "Complete comparison group": f"gs://vpa-full-{P}",
    "Aged 18-35, no menopause":  f"gs://1835-nomnop-{P}",
    "Aged >40, no menopause":    f"gs://40nomnop-{P}",
    "Aged >40, menopause":       f"gs://40mnop-{P}",
}

def read_export(key, bucket):
    files = subprocess.run(["gsutil", "ls", f"{bucket}/*_data_person-*.csv.gz"],
                           capture_output=True, text=True).stdout.split()
    assert files, f"{key}: no person export found in {bucket}"
    stamps = {f.split("/")[-1].split("_data_person")[0] for f in files}
    assert len(stamps) == 1, f"{key}: {bucket} holds more than one export {sorted(stamps)} - choose which export to use"
    local = f"/tmp/cohorts/{bucket.split('/')[-1]}"          # always /tmp, never the notebook folder
    os.makedirs(local, exist_ok=True)
    subprocess.run(["gsutil", "-q", "-m", "cp", *files, local + "/"], check=True)
    d = pd.concat([pd.read_csv(f"{local}/{f.split('/')[-1]}", dtype=str) for f in files], ignore_index=True)
    print(f"{key}: export {stamps.pop()}, {len(files)} file(s)")
    return d

COHORTS = {k: read_export(k, b) for k, b in COHORT_BUCKETS.items()}
for k, d in COHORTS.items():
    need = {"person_id", "T_DISP_race", "T_DISP_sex_at_birth"}
    assert need <= set(d.columns), f"{k}: missing {need - set(d.columns)}. Columns are: {list(d.columns)}"
    assert d["person_id"].is_unique, f"{k}: duplicate person_ids"

# Short-read WGS filter. Some cohorts allowed short- OR long-read WGS, but carriers come from
# cb_variant_to_person (short-read), so long-read-only people would be wrongly counted as non-carriers.
# The flag is taken from whichever table has it filled in (in C2025Q4R6 it was empty in cb_search_person).
def flag_table(t):
    cols = {f.name for f in client.get_table(t).schema}
    lr = "has_lr_whole_genome_variant" if "has_lr_whole_genome_variant" in cols else "NULL"
    q = (f"SELECT CAST(has_whole_genome_variant AS STRING) AS short_read, CAST({lr} AS STRING) AS long_read, "
         f"COUNT(DISTINCT person_id) AS n FROM `{t}` GROUP BY 1, 2 ORDER BY 1, 2")
    return client.query(q).to_dataframe()

chosen = None
for t in [f"{CDR}_index_061026.T_ENT_person", f"{CDR}.cb_search_person"]:
    try:
        ft = flag_table(t)
    except Exception as e:
        print(t, "- not readable:", str(e)[:150], "\n"); continue
    shown = ft.copy(); shown["n"] = shown["n"].map(lambda n: "<20" if 0 < n < 20 else f"{n:,}")
    print(t); print(shown.to_string(index=False)); print()
    n_sr = ft.loc[ft["short_read"].str.lower().isin(["1", "true"]), "n"].sum()
    if chosen is None and n_sr > 100000:
        chosen = t
assert chosen, "Neither table has the short-read flag filled in"
print("Using:", chosen)

ftype = {f.name: f.field_type for f in client.get_table(chosen).schema}["has_whole_genome_variant"]
flag_true = "true" if ftype in ("BOOL", "BOOLEAN") else "1"
sr_ids = set(client.query(f"SELECT DISTINCT CAST(person_id AS STRING) AS person_id FROM `{chosen}` "
                          f"WHERE has_whole_genome_variant = {flag_true}").to_dataframe()["person_id"])
print(f"Participants with short-read WGS: {len(sr_ids):,}")

for k in list(COHORTS):
    d = COHORTS[k]
    before = len(d)
    COHORTS[k] = d[d["person_id"].astype(str).isin(sr_ids)].copy()
    print(f"  {k}: {before:,} -> {len(COHORTS[k]):,} (removed {before - len(COHORTS[k]):,} without short-read WGS)")   # Workbench only

IDS  ={k: set(d["person_id"].astype(str)) for k, d in COHORTS.items()}
RACE = pd.concat([d[["person_id", "T_DISP_race"]].astype(str) for d in COHORTS.values()]) \
         .drop_duplicates("person_id").set_index("person_id")["T_DISP_race"]

# Checks (inside the Workbench only - these print small counts, do not copy them out)
print("Sex at birth values per cohort (should be female only):")
for k, d in COHORTS.items():
    print(f"  {k}: {sorted(d['T_DISP_sex_at_birth'].astype(str).unique())}")
print("Race labels:", sorted(RACE.unique()))   # CELL 7 filters on the exact label "White"

case_key = "Cases (POP 18-35)"
for k in IDS:
    if k != case_key:
        ov = len(IDS[case_key] & IDS[k])
        assert ov == 0, f"{ov} participants are in both the cases and {k}"
for a, b in [("Aged 18-35, no menopause", "Aged >40, no menopause"),
             ("Aged 18-35, no menopause", "Aged >40, menopause"),
             ("Aged >40, no menopause", "Aged >40, menopause")]:
    assert not (IDS[a] & IDS[b]), f"{a} and {b} overlap"
print({k: len(s) for k, s in IDS.items()})

# %% CELL 5a - safety test: does a literal IN list let BigQuery skip most of the table?
# cb_variant_to_person is ~20 TB. Filtering with IN UNNEST(@param) did NOT use the vid clustering and billed the
# whole table. A literal IN list normally does. Every query here has a hard cap (maximum_bytes_billed): if it
# would bill more than the cap it fails and costs nothing.
all_vids = sorted(set().union(*[set(q["vid"]) for q in QUAL.values()]))

def unpack(p):   # person_ids arrives as a nested STRUCT/list; handle every shape it has come back in
    if isinstance(p, dict):
        p = p.get("list", p.get("element", []))
    out = []
    for x in (p if p is not None else []):
        out.append(x.get("element") if isinstance(x, dict) else x)
    return out

def lookup(vids, cap_gb):
    assert all("'" not in v and "\\" not in v for v in vids)
    sql = (f"SELECT vid, person_ids FROM `{CDR}.cb_variant_to_person` "
           f"WHERE vid IN ({', '.join(chr(39) + v + chr(39) for v in vids)})")
    job = client.query(sql, job_config=bigquery.QueryJobConfig(maximum_bytes_billed=int(cap_gb * 1e9)))
    rows = [(r["vid"], str(pid)) for r in job.result() for pid in unpack(r["person_ids"])]
    return rows, (job.total_bytes_billed or 0) / 1e9

test_rows, gb = lookup(all_vids[:20], cap_gb=100)
print(f"test: 20 variants -> {len({v for v, _ in test_rows})} found, {gb:,.1f} GB billed")

# %% CELL 5b - carriers in sorted batches, each capped; saved so it never has to be queried again
CARRIERS_FILE = f"{OUT}/carriers_flat.csv"
BATCH, CAP_GB = 2000, 300            # stop if any batch would bill more than CAP_GB
if subprocess.run(["gsutil", "-q", "stat", CARRIERS_FILE]).returncode == 0:   # already done: no query
    subprocess.run(["gsutil", "-q", "cp", CARRIERS_FILE, "/tmp/vpa/carriers_flat.csv"], check=True)
    carriers = pd.read_csv("/tmp/vpa/carriers_flat.csv", dtype=str)
    print("read saved carriers - no query run")
else:
    rows, total = [], 0.0
    for i in range(0, len(all_vids), BATCH):
        r, gb = lookup(all_vids[i:i + BATCH], cap_gb=CAP_GB)
        rows += r; total += gb
        print(f"  {min(i + BATCH, len(all_vids)):,}/{len(all_vids):,} variants, {total:,.1f} GB billed so far")
    carriers = pd.DataFrame(rows, columns=["vid", "person_id"]).drop_duplicates()
    save(carriers, "carriers_flat.csv", index=False)
assert carriers["person_id"].map(type).eq(str).all()
print(len(all_vids), "target variants;", carriers["vid"].nunique(), "with at least one carrier")

# %% CELL 6 - carrier comparisons: every panel x tier x comparison group
from scipy.stats import fisher_exact

def sup(n):
    return "<20" if 0 < n < 20 else f"{n:,}"

def or_ci(a, b, c, d):          # Woolf CI, Haldane-Anscombe correction if any cell is 0
    if 0 in (a, b, c, d):
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    o = a * d / (b * c); se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    return o, np.exp(np.log(o) - 1.96 * se), np.exp(np.log(o) + 1.96 * se)

def compare(ids, label=""):
    out = []
    for tier, q in QUAL.items():
        for panel, genes in PANELS.items():
            vids = set(q.loc[q["gene"].isin(genes), "vid"])
            carrying = set(carriers.loc[carriers["vid"].isin(vids), "person_id"])
            a = len(ids[case_key] & carrying); n1 = len(ids[case_key])
            for grp in [k for k in ids if k != case_key]:
                c = len(ids[grp] & carrying); n0 = len(ids[grp])
                o, lo, hi = or_ci(a, n1 - a, c, n0 - c)
                p = fisher_exact([[a, n1 - a], [c, n0 - c]])[1]
                out.append({"subset": label, "tier": tier, "panel": panel, "comparison": grp,
                            "case_carriers": a, "n_cases": n1, "comp_carriers": c, "n_comp": n0,
                            "OR": o, "lo": lo, "hi": hi, "p": p})
    return pd.DataFrame(out)

res = compare(IDS, "All participants")
save(res, "results_all.csv", index=False)       # full counts stay in the workspace

def show(df):   # printable version: counts 1-19 suppressed, and ORs/percentages for those rows hidden
    d = df.copy()
    small = d["case_carriers"].between(1, 19) | d["comp_carriers"].between(1, 19)
    d["cases"] = d["case_carriers"].map(sup) + " / " + d["n_cases"].map(sup)
    d["comp"]  = d["comp_carriers"].map(sup) + " / " + d["n_comp"].map(sup)
    d["OR (95% CI)"] = np.where(small, "suppressed",
                                d.apply(lambda r: f"{r.OR:.2f} ({r.lo:.2f}-{r.hi:.2f})", axis=1))
    d["p"] = np.where(small, "", d["p"].map(lambda x: "<0.001" if x < 0.001 else f"{x:.3f}"))
    return d[["subset", "tier", "panel", "comparison", "cases", "comp", "OR (95% CI)", "p"]]

pd.set_option("display.max_rows", 500, "display.width", 250)
print(show(res).to_string(index=False))

# %% CELL 6b - mechanism-level results (primary tier), BH-corrected across the 9 mechanisms
from statsmodels.stats.multitest import multipletests
m = res[res["panel"].str.startswith("Mechanism:") & (res["tier"] == "Primary (coding, >=1*)")].copy()
m["q"] = m.groupby("comparison")["p"].transform(lambda p: multipletests(p, method="fdr_bh")[1])
save(m, "results_mechanisms.csv", index=False)               # full counts stay in the workspace
out = show(m)
out["q (BH)"] = np.where(out["OR (95% CI)"] == "suppressed", "",
                         m["q"].map(lambda x: "<0.001" if x < 0.001 else f"{x:.3f}").values)
print(out.drop(columns=["subset", "tier"]).to_string(index=False))

# %% CELL 7 - sensitivity analyses
# (a) White participants only - the direct check on ancestry
white = {k: {i for i in s if RACE.get(i, "") == "White"} for k, s in IDS.items()}
res_w = compare(white, "White only"); save(res_w, "results_white_only.csv", index=False)
print(show(res_w[res_w["tier"] == "Primary (coding, >=1*)"]).to_string(index=False))

# (b) leave one gene out: ABCC6, and the four genes with indirect mechanisms
for drop in [["ABCC6"], ["HCN4", "MOCS1", "ROBO4", "TNFRSF1A"]]:
    keep = PANELS.copy(); PANELS = {"Wide-147 minus " + "+".join(drop): [x for x in WIDE147 if x not in drop]}
    r = compare(IDS, "leave-out"); PANELS = keep
    print(show(r[r["tier"] == "Primary (coding, >=1*)"]).to_string(index=False))

# (c) dose-response: 0 vs 1 vs >=2 qualifying variants (primary tier, revised wide panel)
q = QUAL["Primary (coding, >=1*)"]; vids = set(q.loc[q["gene"].isin(PANELS["Revised wide"]), "vid"])
per_person = carriers[carriers["vid"].isin(vids)].groupby("person_id")["vid"].nunique()
for grp in [k for k in IDS if k != case_key]:
    for level, test in [("1 variant", lambda n: n == 1), (">=2 variants", lambda n: n >= 2)]:
        ca = sum(test(per_person.get(i, 0)) for i in IDS[case_key]); c0 = sum(per_person.get(i, 0) == 0 for i in IDS[case_key])
        ga = sum(test(per_person.get(i, 0)) for i in IDS[grp]);     g0 = sum(per_person.get(i, 0) == 0 for i in IDS[grp])
        o, lo, hi = or_ci(ca, c0, ga, g0)
        shown = "suppressed (<20)" if (0 < ca < 20 or 0 < ga < 20) else f"OR {o:.2f} ({lo:.2f}-{hi:.2f})"
        print(f"{grp:28s} {level:13s} vs none: {shown}")

# %% CELL 8 - ancestry-adjusted logistic regression (if the ancestry file exists)
import statsmodels.api as sm
subprocess.run(["gsutil", "-u", P, "cp", ANCESTRY_PATH, "/tmp/ancestry_preds.tsv"], check=True)   # /tmp, never "."
anc = pd.read_csv("/tmp/ancestry_preds.tsv", sep="\t", dtype=str)
pcs = anc["pca_features"].str.strip("[]").str.split(",", expand=True).iloc[:, :10].astype(float)
pcs.columns = [f"PC{i}" for i in range(1, 11)]; pcs.index = anc["research_id"].astype(str)
q = QUAL["Primary (coding, >=1*)"]; carrying = set(carriers.loc[carriers["vid"].isin(
        set(q.loc[q["gene"].isin(PANELS["Wide-147 (original)"]), "vid"])), "person_id"])
for grp in [k for k in IDS if k != case_key]:
    people = [i for i in IDS[case_key] | IDS[grp] if i in pcs.index]
    X = pcs.loc[people].copy(); X["carrier"] = [int(i in carrying) for i in people]
    y = np.array([int(i in IDS[case_key]) for i in people])
    fit = sm.Logit(y, sm.add_constant(X)).fit(disp=0)
    b, se = fit.params["carrier"], fit.bse["carrier"]
    print(f"{grp:28s} ancestry-adjusted OR {np.exp(b):.2f} ({np.exp(b-1.96*se):.2f}-{np.exp(b+1.96*se):.2f}), "
          f"p = {fit.pvalues['carrier']:.4f}  (n = {len(people):,})")
