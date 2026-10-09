const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType,
        BorderStyle, AlignmentType, PageOrientation } = require("docx");

const FONT = "Cambria", SIZE = 20, NOTE = 17;          // 10 pt body, 8.5 pt notes
const W = 9071;                                       // 16 cm content width (A4, 2.5 cm margins)
const rule = { style: BorderStyle.SINGLE, size: 6, color: "000000" };
const none = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };

const run = (t, o = {}) => new TextRun({ text: t, font: FONT, size: o.size || SIZE, bold: o.bold, italics: o.italics });
const para = (t, o = {}) => new Paragraph({
  alignment: o.align || AlignmentType.LEFT, spacing: { before: o.before || 0, after: o.after || 0 },
  indent: o.indent ? { left: o.indent } : undefined,
  children: Array.isArray(t) ? t : [run(t, o)] });

function table(widths, header, rows, opts = {}) {
  const n = rows.length;
  const mk = (cells, i, isHead) => new TableRow({
    tableHeader: isHead,
    children: cells.map((c, j) => {
      const txt = typeof c === "object" ? c.t : c;
      const sub = typeof c === "object" && c.sub;
      const bold = isHead || (typeof c === "object" && c.bold);
      return new TableCell({
        width: { size: widths[j], type: WidthType.DXA },
        margins: { top: 50, bottom: 50, left: 80, right: 80 },
        borders: {
          top: isHead ? rule : none,
          bottom: isHead || i === n - 1 ? rule : none,
          left: none, right: none,
        },
        children: [para(txt, { bold, align: j === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
                               indent: sub ? 240 : 0 })],
      });
    }),
  });
  return new Table({
    width: { size: W, type: WidthType.DXA }, columnWidths: widths,
    rows: [mk(header, -1, true), ...rows.map((r, i) => mk(r, i, false))],
  });
}

const caption = (num, text) => para([run(num + " ", { bold: true }), run(text)], { after: 120, before: 360 });
const note = (t) => para(t, { size: NOTE, before: 80 });

// ---------------- Table: UK Biobank matched sample ----------------
const ukb = table([4071, 1900, 1900, 1200],
  ["", "Marfan syndrome (n = 37)", "Matched controls (n = 148)", "SMD"],
  [
    ["Year of birth, mean (SD)", "1952.1 (8.1)", "1951.7 (8.2)", "0.05"],
    ["BMI, mean (SD)", "25.5 (5.7)", "24.8 (4.3)", "0.15"],
    ["Townsend deprivation index, mean (SD)", "−0.34 (2.90)", "−0.57 (3.35)", "0.07"],
    ["Number of live births, mean (SD)", "1.27 (1.15)", "1.26 (1.11)", "0.01"],
    ["Ever smoked, n (%)", "", "", "0.12"],
    [{ t: "No", sub: true }, "15 (40.5)", "56 (37.8)", ""],
    [{ t: "Yes", sub: true }, "21 (56.8)", "90 (60.8)", ""],
    [{ t: "Unknown", sub: true }, "1 (2.7)", "2 (1.4)", ""],
    ["Menopausal status, n (%)", "", "", "0.04"],
    [{ t: "No", sub: true }, "8 (21.6)", "33 (22.3)", ""],
    [{ t: "Yes", sub: true }, "21 (56.8)", "83 (56.1)", ""],
    [{ t: "Hysterectomy", sub: true }, "6 (16.2)", "23 (15.5)", ""],
    [{ t: "Unsure/unknown", sub: true }, "2 (5.4)", "9 (6.1)", ""],
    [{ t: "POP, n (%)", bold: true }, { t: "6 (16.2)", bold: true }, { t: "14 (9.5)", bold: true }, ""],
  ]);

// ---------------- Table: All of Us matched samples and POP ----------------
const aou = table([2250, 1350, 1100, 1150, 1421, 1800],
  ["Condition", "Cases / controls", "POP, cases n (%)", "POP, controls n (%)", "Difference, % (95% CI)", "OR (95% CI)"],
  [
    ["Ehlers–Danlos syndrome", "1,816 / 7,264", "210 (11.6)", "250 (3.4)", "8.1 (6.6–9.7)", { t: "3.78 (3.11–4.61)", bold: true }],
    ["Marfan syndrome", "124 / 496", "<20", "<20", "—", "1.96 (0.90–4.23)"],
    ["Hypermobility syndrome*", "1,001 / 4,004", "88 (8.8)", "158 (3.9)", "4.8 (2.9–6.8)", { t: "2.38 (1.81–3.13)", bold: true }],
    ["Osteogenesis imperfecta*", "73 / 292", "<20", "<20", "—", "1.07 (0.34–3.39)"],
  ]);

// ---------------- Table: All of Us matching balance ----------------
const bal = table([3571, 1500, 1500, 1250, 1250],
  ["Condition", "Cases", "Controls", "SMD age", "SMD BMI"],
  [
    ["Ehlers–Danlos syndrome", "1,816", "7,264", "0.005", "−0.028"],
    ["Marfan syndrome", "124", "496", "−0.058", "0.006"],
    ["Hypermobility syndrome", "1,001", "4,004", "−0.007", "0.022"],
    ["Osteogenesis imperfecta", "73", "292", "0.011", "−0.013"],
  ]);


// ---------------- Table: E-values ----------------
const ev = table([3071, 1700, 2000, 2300],
  ["Condition (cohort)", "OR (95% CI)", "E-value, estimate", "E-value, CI limit"],
  [
    ["Ehlers–Danlos syndrome (All of Us)", "3.78 (3.11–4.61)", "7.02", "5.67"],
    ["Hypermobility syndrome (All of Us)", "2.38 (1.81–3.13)", "4.19", "3.02"],
    ["Marfan syndrome (UK Biobank)", "1.89 (0.65–5.47)", "3.19", "1.00"],
    ["Marfan syndrome (All of Us)", "1.96 (0.90–4.23)", "3.33", "1.00"],
  ]);

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: SIZE } } } },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 },
                          margin: { top: 1417, bottom: 1417, left: 1417, right: 1417 } } },
    children: [
      caption("Table X.", "Characteristics of women with Marfan syndrome and matched controls, UK Biobank."),
      ukb,
      note("Controls were matched 1:4 by propensity score on Townsend deprivation, BMI, smoking status, menopausal status and number of live births, with exact matching on 5-year birth cohort. SMD, standardised mean difference (pooled standard deviation)."),

      caption("Table X.", "POP in matched cases and controls, All of Us."),
      aou,
      note("Cases were matched 1:4 on age and BMI. Difference is the percentage-point difference in POP prevalence. OR from conditional logistic regression stratified by matched set; bold indicates a 95% CI excluding 1. *Exploratory condition. <20: count suppressed in line with the All of Us Data and Statistics Dissemination Policy."),

      caption("Supplementary Table X.", "Covariate balance after matching, All of Us."),
      bal,
      note("SMD, standardised mean difference; values below 0.1 indicate adequate balance."),

      caption("Supplementary Table X.", "E-values for the primary matched estimates."),
      ev,
      note("The E-value is the minimum strength of association, on the risk ratio scale, that an unmeasured confounder would need to have with both the condition and POP to fully explain away the observed association. The CI-limit E-value is the strength needed to shift the confidence interval to include 1; it is 1.00 where the interval already includes 1. ORs were treated as risk ratios because POP prevalence among matched controls was below 10%."),
    ],
  }],
});

Packer.toBuffer(doc).then(b => { fs.writeFileSync("thesis_figures/aim1_tables.docx", b); console.log("ok"); });
