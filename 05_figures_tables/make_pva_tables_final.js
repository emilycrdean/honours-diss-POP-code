// Aim 2 tables from the final (167-gene) rerun. Numbers typed from the SHAREABLE Workbench output only.
const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType, BorderStyle, AlignmentType,
        PageOrientation } = require("docx");
const FONT = "Cambria", SIZE = 17, NOTE = 15;
const rule = { style: BorderStyle.SINGLE, size: 6, color: "000000" }, none = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const run = (t, o = {}) => new TextRun({ text: t, font: FONT, size: o.size || SIZE, bold: o.bold, italics: o.italics });
const para = (t, o = {}) => new Paragraph({ alignment: o.align || AlignmentType.LEFT, spacing: { before: o.before || 0, after: o.after || 0 },
  children: Array.isArray(t) ? t : [run(t, o)] });
function table(widths, header, rows, leftCols = 1, W, pre = null) {
  const n = rows.length;
  const preRow = pre ? [new TableRow({ tableHeader: true, children: pre.map(([t, span]) => new TableCell({ columnSpan: span,
      margins: { top: 40, bottom: 20, left: 70, right: 70 }, borders: { top: rule, bottom: t ? rule : none, left: none, right: none },
      children: [para([run(t, { bold: true })], { align: AlignmentType.CENTER })] })) })] : [];
  const mk = (cells, i, head) => { let col = 0; return new TableRow({ tableHeader: head, cantSplit: true, children: cells.map((c, j) => {
    const txt = typeof c === "object" ? c.t : c, bold = head || (typeof c === "object" && c.bold), ital = typeof c === "object" && c.it;
    const span = typeof c === "object" && c.span ? c.span : 1;
    const w = widths.slice(col, col + span).reduce((a, b) => a + b, 0); const j0 = col; col += span;
    const left = j0 < leftCols || (typeof c === "object" && c.left);
    return new TableCell({ width: { size: w, type: WidthType.DXA }, columnSpan: span, margins: { top: 40, bottom: 40, left: 70, right: 70 },
      borders: { top: head ? rule : none, bottom: head || i === n - 1 ? rule : none, left: none, right: none },
      children: [para([run(txt, { bold, italics: ital })], { align: left ? AlignmentType.LEFT : AlignmentType.CENTER })] }); }) }); };
  return new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: widths, rows: [...preRow, mk(header, -1, true), ...rows.map((r, i) => mk(r, i, false))] });
}
const caption = (num, text) => para([run(num + " ", { bold: true }), run(text)], { after: 120, before: 360 });
const note = t => para(t, { size: NOTE, before: 80 });
const pct = (a, n) => `${a.toLocaleString("en-US")} (${(100 * a / n).toFixed(1)})`;

// ---------------- Table X: main result ----------------
const G = [["Early-onset POP (cases)", 1330], ["Complete comparison group", 320028], ["Aged 18–35, no menopause", 50191],
           ["Aged >40, no menopause", 185861], ["Aged >40, menopause", 50301]];
const panel  = [[135, "—", "—"], [30532, "1.07 (0.90–1.28)", "0.45"], [4803, "1.07 (0.89–1.28)", "0.48"],
                [17651, "1.08 (0.90–1.29)", "0.43"], [4858, "1.06 (0.88–1.27)", "0.54"]];
const negc   = [[319, "—", "—"], [68472, "1.16 (1.02–1.32)", "0.023"], [10627, "1.17 (1.03–1.33)", "0.014"],
                [39190, "1.18 (1.04–1.34)", "0.011"], [11305, "1.09 (0.96–1.24)", "0.20"]];
const W1 = 9071;
const t1 = table([2000, 780, 1130, 1400, 560, 1130, 1400, 560],
  ["Group", "n", "Carriers, n (%)", "OR (95% CI)", "p", "Carriers, n (%)", "OR (95% CI)", "p"],
  G.map(([g, n], i) => [i === 0 ? { t: g, bold: true } : g, n.toLocaleString("en-US"),
    pct(panel[i][0], n), panel[i][1], panel[i][2], pct(negc[i][0], n), negc[i][1], negc[i][2]]), 1, W1,
  [["", 2], ["Connective tissue panel (167 genes)", 3], ["Negative-control genes (10 genes)", 3]]);

// ---------------- Supplementary Table: sensitivity analyses ----------------
const S = [
  ["167-gene panel", null],
  ["Primary tier (main analysis)", ["1.07 (0.90–1.28)", "1.07 (0.89–1.28)", "1.08 (0.90–1.29)", "1.06 (0.88–1.27)"]],
  ["Strict tier", ["1.01 (0.72–1.42)", "1.09 (0.77–1.54)", "1.01 (0.72–1.42)", "0.94 (0.67–1.33)"]],
  ["Extended tier", ["0.93 (0.83–1.03)", "0.91 (0.82–1.02)", "0.92 (0.82–1.02)", "0.95 (0.85–1.06)"]],
  ["Primary tier, White participants only*", ["0.99 (0.79–1.25)", "0.95 (0.75–1.20)", "1.00 (0.80–1.25)", "1.00 (0.80–1.26)"]],
  ["Primary tier, one qualifying variant vs none†", ["1.08 (0.90–1.30)", "1.08 (0.90–1.29)", "1.09 (0.91–1.30)", "1.06 (0.88–1.28)"]],
  ["Original 146-gene panel", null],
  ["Primary tier", ["1.07 (0.89–1.28)", "1.07 (0.89–1.29)", "1.08 (0.90–1.29)", "1.05 (0.88–1.26)"]],
  ["Primary tier, adjusted for 10 ancestry PCs", ["1.07 (0.89–1.28)", "1.03 (0.86–1.24)", "1.07 (0.90–1.29)", "1.07 (0.89–1.29)"]],
  ["Primary tier, excluding ABCC6", ["1.06 (0.86–1.31)", "1.05 (0.85–1.29)", "1.06 (0.86–1.31)", "1.07 (0.87–1.32)"]],
  ["Primary tier, excluding HCN4, MOCS1, ROBO4, TNFRSF1A‡", ["1.03 (0.85–1.25)", "1.03 (0.85–1.25)", "1.03 (0.85–1.25)", "1.01 (0.83–1.23)"]],
  ["Extended tier", ["0.86 (0.77–0.97)", "0.85 (0.76–0.95)", "0.86 (0.77–0.96)", "0.90 (0.81–1.01)"]],
  ["Original panel plus 24 GWAS-locus genes (170 genes)", null],
  ["Primary tier", ["1.06 (0.89–1.26)", "1.06 (0.89–1.27)", "1.07 (0.89–1.27)", "1.04 (0.87–1.24)"]],
  ["Negative-control genes (10 genes)", null],
  ["Primary tier", ["1.16 (1.02–1.32)", "1.17 (1.03–1.33)", "1.18 (1.04–1.34)", "1.09 (0.96–1.24)"]],
  ["Primary tier, White participants only*", ["1.02 (0.87–1.19)", "1.02 (0.87–1.19)", "1.02 (0.87–1.19)", "1.01 (0.86–1.18)"]],
  ["Strict tier", ["1.04 (0.85–1.27)", "1.06 (0.86–1.30)", "1.06 (0.87–1.30)", "0.95 (0.77–1.16)"]],
  ["Extended tier", ["1.10 (0.97–1.23)", "1.09 (0.97–1.23)", "1.11 (0.98–1.25)", "1.07 (0.94–1.20)"]],
];
const W2 = 9071;
const t2 = table([3071, 1500, 1500, 1500, 1500],
  ["Analysis", "Complete comparison group", "Aged 18–35, no menopause", "Aged >40, no menopause", "Aged >40, menopause"],
  S.map(([lab, v]) => v ? ["   " + lab, ...v] : [{ t: lab, bold: true, span: 5, left: true }]), 1, W2);

// ---------------- Supplementary Table: mechanisms ----------------
const N0 = { "Complete comparison group": 320028, "Aged 18–35, no menopause": 50191, "Aged >40, no menopause": 185861, "Aged >40, menopause": 50301 };
const M = [
  ["Structural collagens and collagen-associated matrix", 19, "32 (2.4)", [
    ["Complete comparison group", 7366, "1.05 (0.74–1.49)", "0.78", "1.00"], ["Aged 18–35, no menopause", 1192, "1.01 (0.71–1.45)", "0.93", "1.00"],
    ["Aged >40, no menopause", 4301, "1.04 (0.73–1.48)", "0.78", "1.00"], ["Aged >40, menopause", 1104, "1.10 (0.77–1.57)", "0.57", "1.00"]]],
  ["Elastic fibres and microfibrils", 25, "49 (3.7)", [
    ["Complete comparison group", 11098, "1.06 (0.80–1.42)", "0.65", "1.00"], ["Aged 18–35, no menopause", 1658, "1.12 (0.84–1.50)", "0.44", "0.79"],
    ["Aged >40, no menopause", 6375, "1.08 (0.81–1.43)", "0.60", "1.00"], ["Aged >40, menopause", 1852, "1.00 (0.75–1.34)", "1.00", "1.00"]]],
  ["Syndromic or other, no direct connective tissue mechanism", 24, "24 (1.8)", [
    ["Complete comparison group", 4875, "1.19 (0.79–1.78)", "0.37", "0.88"], ["Aged 18–35, no menopause", 728, "1.25 (0.83–1.88)", "0.30", "0.79"],
    ["Aged >40, no menopause", 2853, "1.18 (0.79–1.77)", "0.43", "0.97"], ["Aged >40, menopause", 808, "1.13 (0.75–1.70)", "0.58", "1.00"]]],
  ["Collagen biosynthesis, modification and cross-linking", 32, "<20", null],
  ["Proteoglycans and glycosaminoglycans", 12, "<20", null],
  ["TGF-β and BMP signalling", 17, "<20", null],
  ["Smooth muscle, cytoskeleton and vessel wall", 15, "<20", null],
  ["Matrix remodelling, complement and inflammation", 14, "<20", null],
  ["Hormonal signalling and pelvic development", 9, "0", null],
];
const rows3 = [];
for (const [mech, ng, cc, groups] of M) {
  if (!groups) { rows3.push([mech, String(ng), cc, { t: cc === "0" ? "No case carriers" : "Tested; estimates not reported (fewer than 20 case carriers)", span: 5, left: true }]); continue; }
  groups.forEach(([g, c, or, p, q], k) => rows3.push([k === 0 ? { t: mech, bold: true } : "", k === 0 ? String(ng) : "", k === 0 ? cc : "",
    { t: g, left: true }, pct(c, N0[g]), or, p, q]));
}
const W3 = 13958;   // landscape
const t3 = table([3300, 700, 1300, 2500, 1700, 1900, 1200, 1358],
  ["Mechanism", "Genes", "Case carriers, n (%)", "Comparison group", "Carriers, n (%)", "OR (95% CI)", "p", "q (BH)"], rows3, 1, W3);

const portrait = { page: { size: { width: 11906, height: 16838 }, margin: { top: 1417, bottom: 1417, left: 1417, right: 1417 } } };
const landscape = { page: { size: { width: 11906, height: 16838, orientation: PageOrientation.LANDSCAPE }, margin: { top: 1417, bottom: 1417, left: 1417, right: 1417 } } };
const doc = new Document({ styles: { default: { document: { run: { font: FONT, size: SIZE } } } },
  sections: [
    { properties: portrait, children: [
      caption("Table X.", "Carriers of rare pathogenic or likely pathogenic variants in early-onset POP cases and four comparison groups, connective tissue panel and negative-control genes (primary tier)."),
      t1,
      note("Primary tier: coding variants (missense, stop-gained, frameshift, splice donor or acceptor, start-lost, stop-lost and in-frame insertions or deletions) with allele frequency below 1% and at least one ClinVar review star, classified pathogenic or likely pathogenic. All groups are people assigned female at birth with short-read whole-genome sequencing. OR (95% CI) for cases vs each comparison group (Woolf method); p from two-sided Fisher's exact test. Bonferroni-corrected threshold for the connective tissue panel: p < 0.0125 across four comparisons."),
      caption("Supplementary Table X.", "Odds ratios (95% CI) for carrying a qualifying variant, early-onset POP cases vs each comparison group, across variant tiers, panels and sensitivity analyses."),
      t2,
      note("Strict tier: primary-tier variants with at least two ClinVar review stars and only pathogenic or likely pathogenic submissions. Extended tier: primary-tier variants plus loss-of-function variants and missense variants with REVEL ≥ 0.75, regardless of ClinVar classification. *857 cases. †Carriers of exactly one qualifying variant vs non-carriers; estimates for carriers of two or more variants rested on fewer than 20 cases and are not reported. ‡Genes carrying relatively common variants with conflicting ClinVar classifications. Ancestry-adjusted estimates from logistic regression with 10 genetic principal components. Analyses in the original panel were run before panel curation and are shown for completeness."),
    ] },
    { properties: landscape, children: [
      caption("Supplementary Table X.", "Carriers of qualifying variants by biological mechanism, connective tissue panel (primary tier)."),
      t3,
      note("Each gene was assigned to one mechanism before the results were examined. All nine mechanisms were tested; estimates for mechanisms with 1–19 case carriers are not reported because they would allow those counts to be derived. q: Benjamini–Hochberg false discovery rate across the nine mechanisms within each comparison group. Counts of 1–19 are shown as <20, in line with the All of Us Data and Statistics Dissemination Policy."),
    ] },
  ] });
Packer.toBuffer(doc).then(b => { fs.writeFileSync("thesis_figures/pva_tables_final.docx", b); console.log("ok"); });
