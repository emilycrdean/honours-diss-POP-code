const fs = require('fs');
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType, BorderStyle,
        ShadingType, PageOrientation } = require('docx');

const FONT = 'Cambria', SZ = 18, LINE = 360, CHECK = 'FFF2B3';
const C = { confirm: true };

function cell(content, w, opts = {}) {
  const parts = Array.isArray(content) ? content : [content];
  const paras = []; let runs = [];
  for (const p of parts) {
    const [t, o] = Array.isArray(p) ? p : [p, {}];
    t.split('\n').forEach((ln, i) => {
      if (i > 0) { paras.push(runs); runs = []; }
      if (ln) runs.push(new TextRun({ text: ln, font: FONT, size: SZ, bold: opts.bold || o.bold, italics: o.italic,
        shading: o.confirm ? { type: ShadingType.CLEAR, color: 'auto', fill: CHECK } : undefined }));
    });
  }
  paras.push(runs);
  const none = { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' };
  const line = s => ({ style: BorderStyle.SINGLE, size: s, color: '000000' });
  return new TableCell({
    width: { size: w, type: WidthType.DXA }, columnSpan: opts.span,
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    borders: { top: opts.top ? line(8) : none, bottom: opts.bottom ? line(opts.bottom) : none, left: none, right: none },
    children: paras.map(r => new Paragraph({ spacing: { line: LINE, before: 0, after: 0 }, children: r })),
  });
}

const W = [2300, 4200, 4700, 2758];
const TW = W.reduce((a, b) => a + b);
const rows = [
  // [concept set, parent concepts selected, what they capture, exclusions / use]
  ['Study population', 'Sex at birth: female\nAt least one EHR record\nShort-read whole-genome sequence data passing quality control (Aim 2 only)', '—', 'All analyses'],
  ['Pelvic organ prolapse (diagnosis)', 'ICD-10-CM: N81 Female genital prolapse; N99.3 Prolapse of vaginal vault after hysterectomy\nICD-9-CM: 618 Genital prolapse\nSNOMED CT: Prolapse of female genital organs',
    'Uterine prolapse (all degrees); incomplete and complete uterovaginal prolapse; anterior and posterior vaginal wall prolapse; cystocele (midline, lateral); rectocele; vaginal enterocele; urethrocele; vaginal vault prolapse, including after hysterectomy; cervical stump prolapse; incompetence or weakening of rectovaginal tissue; incontinence or overactive bladder due to genital prolapse',
    'Removed: rectal and anal prolapse (K62.2, K62.3; 569.1); prolapse in pregnancy; prolapse of ovary or fallopian tube\nUse: Aim 1 outcome; excluded from Aim 2 comparison groups'],
  ['Pelvic organ prolapse (surgical repair)', 'CPT-4 and SNOMED CT: Repair of vaginal wall prolapse',
    'Anterior, posterior and combined colporrhaphy; repair of cystocele, rectocele, urethrocele and enterocele; colpopexy and sacrocolpopexy (vaginal, abdominal, laparoscopic, robot-assisted); colpocleisis; vaginal hysterectomy with enterocele repair; mesh insertion for pelvic floor defect',
    'Use: Aim 2 case definition (POP diagnosis or repair at age 18–35, no CTD diagnosis)'],
  ['Menopause', 'ICD-10-CM: E28.31 Premature menopause; Z78.0 Asymptomatic menopausal state\nSNOMED CT: Post-hysterectomy menopause; postartificial menopausal syndrome; postsurgical menopause\nObservation: Menopause; menopause finding; menopause present; premature menopause',
    'Natural, premature (symptomatic and asymptomatic) and surgical menopause; postmenopausal status recorded without symptoms',
    'Use: Aim 1 sensitivity covariate; defines Aim 2 comparison groups'],
  ['Ehlers–Danlos syndrome', 'SNOMED CT: Ehlers–Danlos syndrome\nICD-10-CM: Q79.6\nICD-9-CM: 756.83', 'All EDS subtypes', 'Marfan syndrome excluded\nUse: Aim 1 exposure; Aim 2 exclusion'],
  ['Marfan syndrome', 'SNOMED CT: Marfan syndrome\nICD-10-CM: Q87.4\nICD-9-CM: 759.82', 'Marfan syndrome, including with aortic or ocular involvement', 'Use: Aim 1 exposure; Aim 2 exclusion'],
  ['Hypermobility syndrome', 'SNOMED CT: Hypermobility syndrome\nICD-10-CM: M35.7\nICD-9-CM: 728.5', 'Hypermobility syndrome; familial articular hypermobility', 'Use: Aim 1 exposure (exploratory); Aim 2 exclusion'],
  ['Osteogenesis imperfecta', 'SNOMED CT: Osteogenesis imperfecta\nICD-10-CM: Q78.0\nICD-9-CM: 756.51', 'All OI types', 'Use: Aim 1 exposure (exploratory); Aim 2 exclusion'],
  ['Body mass index', 'Measurement: Body mass index', 'All recorded values; median per participant', 'Missing BMI excluded from the primary match'],
  ['Smoking', 'Lifestyle survey: "Have you smoked at least 100 cigarettes in your entire life?"\nObservation: Tobacco smoking status', 'Ever / never smoked', 'Missing retained as unknown'],
  ['Deprivation', 'Zip code socioeconomic status data', 'Composite area-level index (American Community Survey; three-digit zip code)', 'Higher values, greater deprivation'],
  ['Race and ethnicity', 'Basics survey: race; ethnicity', 'Combined into a single variable', 'Not reported retained as a category'],
  ['Vaginal delivery', 'CPT-4: Vaginal delivery, antepartum and postpartum care procedures', 'Vaginal deliveries, including assisted delivery', 'Explored to derive parity; not used (recorded for <1% of participants)'],
];

const header = new TableRow({ tableHeader: true, children: ['Concept set', 'Parent concepts selected', 'Captures (descendant concepts)', 'Exclusions and use']
  .map((h, i) => cell(h, W[i], { bold: true, top: true, bottom: 6 })) });
const body = rows.map((r, k) => new TableRow({ cantSplit: true, children: r.map((c, i) =>
  cell(i === 0 ? [[c, { bold: true }]] : c, W[i], { bottom: k === rows.length - 1 ? 8 : 2 })) }));

const p = (runs, after = 120) => new Paragraph({ spacing: { line: LINE, after }, children: runs });
const doc = new Document({ sections: [{
  properties: { page: { size: { width: 11906, height: 16838, orientation: PageOrientation.LANDSCAPE },
                        margin: { top: 1134, bottom: 1134, left: 1134, right: 1134 } } },
  children: [
    p([new TextRun({ text: 'Supplementary Table S2. ', bold: true, font: FONT, size: 22 }),
       new TextRun({ text: 'Concept sets used to define study populations, outcomes, exposures and covariates in All of Us.', font: FONT, size: 22 })]),
    new Table({ width: { size: TW, type: WidthType.DXA }, columnWidths: W, rows: [header, ...body] }),
    p([new TextRun({ text: 'Concept sets were built in the All of Us Cohort Builder (Controlled Tier, CDR C2025Q4R6) and include all descendant concepts of each concept listed. EHR data are mapped to the OMOP common data model; source ICD-10-CM, ICD-9-CM and CPT-4 codes are shown alongside SNOMED CT standard concepts. CPT-4, Current Procedural Terminology, 4th edition; CTD, connective tissue disorder; EHR, electronic health record; POP, pelvic organ prolapse.', font: FONT, size: 18 })], 0),
  ] }] });
Packer.toBuffer(doc).then(b => { fs.writeFileSync('thesis_figures/supp_aou_concepts.docx', b); console.log('ok'); });
