# Builds vpa_gene_list.csv: every gene, which panel(s) it belongs to, and whether it is kept, added or proposed for removal.
import csv

WIDE147 = """ABCC6 ABL1 ACTA2 ACVR1 ADAMTS10 ADAMTS17 ADAMTS2 ADAMTSL2 ADAMTSL4 AEBP1 ALDH18A1 ALPL ARIH1 ASPH ATP6AP1
ATP6V0A2 ATP6V1A ATP6V1E1 ATP7A B3GALT6 B3GAT3 B4GALT7 BGN BPNT2 C1R C1S CBS CHD4 CHST14 CHST3 COG7 COL11A1 COL11A2
COL12A1 COL1A1 COL1A2 COL2A1 COL3A1 COL4A1 COL5A1 COL5A2 COL6A2 COL6A3 COL9A1 COL9A2 COL9A3 CREB3L1 CRTAP CSGALNACT1
DCHS1 DLG4 DSE EFEMP2 ELN ERBIN FBLN5 FBN1 FBN2 FKBP14 FLCN FLNA FLNB FN1 FOXE3 GATA5 GGCX GORAB HCN4 IPO8 KIF22 KLK1
KLK10 KLK11 KLK12 KLK13 KLK14 KLK15 KLK2 KLK3 KLK4 KLK5 KLK6 KLK7 KLK8 LEMD3 LOX LOXL3 LRP2 LTBP2 LTBP3 LTBP4 LZTS1
MAT2A MBTPS2 MED12 MFAP5 MOCS1 MYH11 MYLK NEK9 NEPRO NKAP NOG NOTCH1 NT5E P3H1 P4HA1 PKD2 PLOD1 PLOD3 PLP1 PRDM5 PRKG1
PYCR1 RIN2 ROBO3 ROBO4 SERPING1 SGMS2 SHOC2 SKI SLC26A2 SLC2A10 SLC39A13 SMAD2 SMAD3 SMAD4 SMAD6 SMS SP7 SPARC TAB2
TALDO1 TGFB1 TGFB2 TGFB3 TGFBR1 TGFBR2 TLN1 TNFRSF1A TNXB UPF3B VCAN WNT1 ZDHHC9 ZNF469""".split()

CORE10 = "FBN1 FBLN5 COL1A1 COL1A2 COL3A1 COL5A1 COL5A2 TGFBR1 TGFBR2 TNXB".split()

ADDON24 = """ACADVL CHRDL2 CRISPLD2 DUSP16 DVL2 EFEMP1 ESR1 FAT4 GDF7 HOXA11 HOXD13 IMPDH1 KLF13 LAMC1 LDAH LOXL1
LOXL1-AS1 MAFF SALL1 TBX5 TIMP2 VCL WNT4 WT1""".split()

ADD = {  # gene: (pathway, reason)
    "COL6A1":   ("Structural collagens", "Third chain of collagen VI (Bethlem/Ullrich); COL6A2/COL6A3 already included"),
    "SERPINH1": ("Collagen biosynthesis", "HSP47 collagen chaperone; recessive OI"),
    "FKBP10":   ("Collagen biosynthesis", "Collagen folding; recessive OI / Bruck syndrome"),
    "PPIB":     ("Collagen biosynthesis", "Collagen prolyl isomerase; recessive OI"),
    "PLOD2":    ("Collagen biosynthesis", "Telopeptide lysyl hydroxylase; Bruck syndrome"),
    "P4HA2":    ("Collagen biosynthesis", "Collagen prolyl 4-hydroxylase; P4HA1 already included"),
    "SERPINF1": ("Collagen biosynthesis", "Recessive OI type VI"),
    "BMP1":     ("Collagen biosynthesis", "Procollagen C-proteinase; recessive OI"),
    "TMEM38B":  ("Collagen biosynthesis", "Recessive OI"),
    "PLS3":     ("Collagen biosynthesis", "X-linked osteoporosis/OI-like bone fragility"),
    "IFITM5":   ("Collagen biosynthesis", "OI type V; main pathogenic variant is non-coding (5'UTR) - expect it to drop at the coding filter"),
    "LOXL2":    ("Collagen cross-linking", "Lysyl oxidase family; completes LOX/LOXL1/LOXL3/LOXL4"),
    "LOXL4":    ("Collagen cross-linking", "Re-added on mechanism (lysyl oxidase family); was dropped for a null POP citation"),
    "MFAP4":    ("Elastic fibre assembly", "Microfibril-associated protein; MFAP5 already included"),
    "MMP1":     ("Matrix remodelling", "Interstitial collagenase; MMP/TIMP imbalance in prolapsed tissue"),
    "MMP2":     ("Matrix remodelling", "Gelatinase A; multicentric osteolysis"),
    "MMP3":     ("Matrix remodelling", "Stromelysin; POP candidate-gene studies"),
    "MMP9":     ("Matrix remodelling", "Gelatinase B; POP candidate-gene studies"),
    "TIMP1":    ("Matrix remodelling", "MMP inhibitor; TIMP2 already included"),
    "TIMP3":    ("Matrix remodelling", "MMP inhibitor; Sorsby fundus dystrophy"),
    "IL1B":     ("Inflammation", "Drives MMP activation and matrix breakdown"),
    "IL6":      ("Inflammation", "Inflammatory matrix degradation"),
    "NLRP3":    ("Inflammation", "Cryopyrin-associated periodic syndrome; autoinflammatory like TNFRSF1A/TRAPS"),
    "ESR2":     ("Hormonal regulation", "Oestrogen receptor beta; ESR1 already included"),
    "PGR":      ("Hormonal regulation", "Progesterone receptor; POP candidate-gene studies"),
    "HOXA13":   ("Pelvic floor development", "Hand-foot-genital syndrome (uterine malformations); HOXA11 already included"),
    "DES":      ("Smooth muscle", "Re-added on mechanism (pelvic floor smooth/skeletal muscle); was dropped for lack of a POP citation"),
}

REMOVE = {  # proposed - decide before running
    **{f"KLK{i}": "Kallikrein cluster bystander: no connective tissue role and no ClinVar P/LP entries (KLK15 kept as hEDS candidate)"
       for i in [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14]},
    "DLG4":  "Neuronal (PSD-95); entered via the POP GWAS locus it shares with ACADVL/DVL2, all prior signal was upstream non-coding calls",
    "PLP1":  "CNS myelin protein (Pelizaeus-Merzbacher); no connective tissue or pelvic mechanism",
    "ROBO3": "Axon guidance (gaze palsy with scoliosis); scoliosis is secondary to neurological disease",
    "CHRDL2": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "CRISPLD2": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "DUSP16": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "DVL2": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "FAT4": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "GDF7": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "IMPDH1": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "KLF13": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "LAMC1": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "LDAH": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "MAFF": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "TBX5": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "LOXL1-AS1": "Included only because it lies at a POP GWAS locus; rare coding variants test a different hypothesis from the common GWAS signal",
    "ACADVL": "Entered via the POP GWAS locus; its ClinVar P/LP variants cause VLCAD deficiency (recessive, metabolic), not a connective tissue mechanism",
}

NEG_CONTROL = {  # recessive carrier-screening genes with no connective tissue role
    "CFTR": "Cystic fibrosis", "GJB2": "Hearing loss", "PAH": "Phenylketonuria", "HEXA": "Tay-Sachs",
    "ATP7B": "Wilson disease", "SLC26A4": "Pendred syndrome", "ACADM": "MCAD deficiency",
    "BTD": "Biotinidase deficiency", "GALT": "Galactosaemia", "MMUT": "Methylmalonic acidaemia",
}

rows = []
for g in sorted(set(WIDE147) | set(ADDON24) | set(ADD) | set(NEG_CONTROL)):
    panels = [p for p, s in [("Core-10", CORE10), ("Wide-147", WIDE147), ("Add-on-24", ADDON24)] if g in s]
    if g in NEG_CONTROL:
        status, pathway, reason = "Negative control", "Negative control", NEG_CONTROL[g]
    elif g in REMOVE:
        status, pathway, reason = "Proposed removal", "", REMOVE[g]
    elif g in ADD:
        status, (pathway, reason) = "New", ADD[g]
    else:
        status, pathway, reason = "Keep", "", ""
    rows.append({"gene": g, "existing_panels": ";".join(panels), "status": status, "pathway": pathway, "reason": reason})

with open("vpa_gene_list.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

from collections import Counter
print(Counter(r["status"] for r in rows))
print("Revised wide panel size:",
      len((set(WIDE147) | set(ADDON24) | set(ADD)) - set(REMOVE)))
