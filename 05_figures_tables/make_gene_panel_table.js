const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType,
        BorderStyle, AlignmentType } = require("docx");

const FONT = "Cambria", SIZE = 20, NOTE = 17;
const W = 9071;
const rule = { style: BorderStyle.SINGLE, size: 6, color: "000000" };
const none = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const run = (t, o = {}) => new TextRun({ text: t, font: FONT, size: o.size || SIZE, bold: o.bold,
                                          italics: o.italics, superScript: o.sup });
const para = (children, o = {}) => new Paragraph({ alignment: o.align || AlignmentType.LEFT,
  spacing: { before: o.before || 0, after: o.after || 0 }, children });

// ---- data ----
const csv = (f) => fs.readFileSync(f, "utf8").trim().split(/\r?\n/).map(l => {
  const out = []; let cur = "", q = false;
  for (const ch of l) { if (ch === '"') q = !q; else if (ch === "," && !q) { out.push(cur); cur = ""; } else cur += ch; }
  out.push(cur); return out; });
const gl = csv("vpa_gene_list.csv"); const head = gl.shift();
const G = Object.fromEntries(gl.map(r => [r[0], Object.fromEntries(head.map((h, i) => [h, r[i]]))]));
const mech = csv("vpa_mechanisms.csv"); mech.shift();
const order = [];
const byMech = {};
for (const [g, m] of mech) { if (!byMech[m]) { byMech[m] = []; order.push(m); } byMech[m].push(g); }

const isCore = g => (G[g]?.existing_panels || "").includes("Core-10");
const marker = g => G[g]?.status === "New" ? "‡" : (G[g]?.existing_panels === "Add-on-24" ? "†" : "");

function geneRuns(genes, marks = true) {
  const out = [];
  genes.sort().forEach((g, i) => {
    out.push(run(g, { italics: true, bold: isCore(g) }));
    const m = marks ? marker(g) : ""; if (m) out.push(run(m, { sup: true }));
    if (i < genes.length - 1) out.push(run(", "));
  });
  return out;
}

function table(widths, header, rows) {
  const n = rows.length;
  const cell = (content, j, i, isHead) => new TableCell({
    width: { size: widths[j], type: WidthType.DXA },
    margins: { top: 60, bottom: 60, left: 80, right: 80 },
    borders: { top: isHead ? rule : none, bottom: isHead || i === n - 1 ? rule : none, left: none, right: none },
    children: [para(Array.isArray(content) ? content : [run(content, { bold: isHead })],
                    { align: j === 1 && widths.length === 3 ? AlignmentType.CENTER : AlignmentType.LEFT })],
  });
  return new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: widths,
    rows: [new TableRow({ tableHeader: true, children: header.map((h, j) => cell(h, j, -1, true)) }),
           ...rows.map((r, i) => new TableRow({ cantSplit: true, children: r.map((c, j) => cell(c, j, i, false)) }))] });
}

// ---- Part A: final panel by mechanism ----
const total = order.reduce((s, m) => s + byMech[m].length, 0);
const partA = table([2900, 700, 5471], ["Mechanism", "Genes", "Gene symbols"],
  order.map(m => [m, String(byMech[m].length), geneRuns(byMech[m].slice())]));

// ---- Part B: genes removed during curation ----
const rm = Object.values(G).filter(r => r.status === "Proposed removal").map(r => r.gene);
const klk = rm.filter(g => g.startsWith("KLK"));
const neuro = ["DLG4", "PLP1", "ROBO3"];
const gwas = rm.filter(g => G[g].existing_panels === "Add-on-24" && g !== "ACADVL");
const partB = table([3400, 5671], ["Genes removed", "Reason"], [
  [geneRuns(klk, false), "Kallikrein gene cluster: no connective tissue role and no pathogenic or likely pathogenic ClinVar entries (KLK15 retained as a candidate gene for hypermobile EDS)"],
  [geneRuns(neuro.slice(), false), "Nervous system genes (postsynaptic density protein, central myelin, axon guidance) with no connective tissue or pelvic mechanism"],
  [geneRuns(gwas, false), "Included only because they lie near a POP genome-wide association locus; a common non-coding association signal does not imply that rare coding variants in the nearest gene affect POP risk"],
  [geneRuns(["ACADVL"], false), "Near a POP association locus; its pathogenic variants cause very-long-chain acyl-CoA dehydrogenase deficiency, a recessive metabolic disorder"],
]);

// ---- Part C: negative controls ----
const neg = Object.values(G).filter(r => r.status === "Negative control").sort((a, b) => a.gene.localeCompare(b.gene));
const partC = table([3400, 5671], ["Gene", "Associated condition"],
  neg.map(r => [[run(r.gene, { italics: true })], r.reason]));

const nRm = rm.length, nNew = Object.values(G).filter(r => r.status === "New").length;
const caption = (num, text) => para([run(num + " ", { bold: true }), run(text)], { before: 360, after: 120 });
const note = (t) => para([run(t, { size: NOTE })], { before: 80 });
const sub = (t) => para([run(t, { bold: true })], { before: 240, after: 80 });

const doc = new Document({ sections: [{
  properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1418, bottom: 1418, left: 1418, right: 1418 } } },
  children: [
    caption("Supplementary Table X.", `Connective tissue gene panel used in the pathogenic variant analysis (${total} genes).`),
    sub("A. Final panel, grouped by biological mechanism"),
    partA,
    note(`Genes were compiled from the union of eight clinical diagnostic panels for heritable connective tissue disorders (Invitae, Fulgent, GeneDx, PreventionGenetics and Greenwood Genetic Centre connective tissue panels, a Marfan syndrome panel and two Ehlers–Danlos syndrome panels) and genes reported in POP genome-wide association and candidate gene studies, then curated as described in Section 2.3. Each gene was assigned to one mechanism before analysis. Bold, core panel genes with established roles in Marfan syndrome, Ehlers–Danlos syndrome and Loeys–Dietz syndrome. † Gene from the POP association literature retained because it has an independent connective tissue or pelvic mechanism. ‡ Gene added during curation (${nNew} genes), including recessive forms of osteogenesis imperfecta and genes for collagen processing and matrix remodelling. Unmarked genes are from the clinical diagnostic panels.`),
    sub(`B. Genes removed during curation (${nRm} genes)`),
    partB,
    sub(`C. Negative-control genes (${neg.length} genes)`),
    partC,
    note("Negative-control genes cause recessive disorders with no known connective tissue or pelvic involvement. They were analysed with the same variant filters to estimate the background difference in carrier frequency between groups."),
  ] }] });

Packer.toBuffer(doc).then(b => { fs.writeFileSync("thesis_figures/supp_gene_panel.docx", b);
  console.log("saved", total, "genes;", nRm, "removed;", nNew, "added;", neg.length, "negative controls"); });
