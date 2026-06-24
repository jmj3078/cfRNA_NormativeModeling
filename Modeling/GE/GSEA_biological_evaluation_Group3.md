# GSEA Biological Evaluation — Group 3: Pancreas/Lung Diseases

**Analysis date:** 2026-06-24  
**Model:** GAMLSS normative model trained on 693 healthy controls (HC), plasma cfRNA  
**GSEA method:** Prerank on mean Z-score per phenotype (gene-level deviation from HC normative distribution)  
**Gene sets:** KEGG 2021 Human, Reactome 2022, GO Biological Process 2023  
**Significance threshold:** FDR q-val < 0.05  
**NES interpretation:** NES > 0 = enriched among high-Z (upregulated) genes; NES < 0 = enriched among low-Z (downregulated) genes

---

## Source Papers

| Disease | Study | DOI | n |
|---------|-------|-----|---|
| Pancreatic Cancer (Moore) | Moore et al., *Nature Communications* 2025 | https://doi.org/10.1038/s41467-025-62685-y | 72 |
| Pancreatic Cancer (Reggiardo) | Reggiardo et al., *Nature Biomedical Engineering* 2023 | https://doi.org/10.1038/s41551-023-01081-7 | 6 |
| Pancreatitis (Moore) | Moore et al., *Nature Communications* 2025 | https://doi.org/10.1038/s41467-025-62685-y | 79 |
| Lung Cancer | Chen et al., *eLife* 2022 | https://doi.org/10.7554/elife.75181 | 30 |

---

## Study Design Notes

**Moore et al. 2025 (Nature Communications):** Collected plasma cfRNA from patients referred for endoscopic ultrasound-guided fine-needle aspiration (EUS-FNA). Discovery cohort: 153 samples (CEDAR); validation cohort: 95 samples (BCC). Standard short-read total RNA-seq of plasma. Derived 29 cfRNA biomarkers for PDAC, many of which relate to liver function (elevated even in early-stage PDAC without liver metastasis). Both Pancreatic Cancer (n=72) and Pancreatitis (n=79) come from this same cohort, making direct comparison particularly meaningful.

**Reggiardo et al. 2023 (Nature Biomedical Engineering):** COMPLETE-seq — a repeat-element-aware liquid biopsy approach combining custom annotation of >5 million repeat element insertions (~15,000 repeat subfamilies) with Illumina short-read and Oxford Nanopore long-read sequencing. The key innovation is profiling transposable element-derived cfRNA (SINEs, LINEs, Alu elements) alongside canonical gene-level expression. Pancreatic cancer samples were enriched for specific Alu subfamilies (AluY, AluSc, AluSg7, AluSc8, AluSx3, AluSg). Only 6 pancreatic cancer samples in the normative modeling cohort.

**Chen et al. 2022 (eLife):** SMART-seq total RNA-seq of plasma cfRNA (>50 nt) profiling both human and microbial cfRNA from ~300 samples across five cancer types including lung cancer (n=35 in original study; n=30 in this cohort). Identified lung cancer-specific biomarkers including IL1R2 and CLEC4E related to immune regulation.

---

## 1. Pancreatic Cancer (Moore, n=72)

### Summary Statistics
- Total significant terms (FDR < 0.05): **554**
- NES > 0 (upregulated): **240**
- NES < 0 (downregulated): **314**

### Top Upregulated Pathways (NES > 0, high-Z genes)

| NES | FDR | Pathway |
|-----|-----|---------|
| +2.495 | 0.000 | KEGG: Complement and coagulation cascades |
| +2.484 | 0.000 | Reactome: RHOC GTPase Cycle |
| +2.456 | 0.000 | Reactome: Post-translational Protein Phosphorylation |
| +2.445 | 0.000 | GO: Hippo Signaling |
| +2.400 | 0.000 | Reactome: Regulation of IGF Transport and Uptake by IGFBPs |
| +2.390 | 0.000 | KEGG: ECM-receptor interaction |
| +2.370 | 0.000 | GO: Negative Regulation of Blood Coagulation |
| +2.315 | 0.001 | GO: Regulation of Endothelial Cell Migration |
| +2.299 | 0.001 | GO: Collagen Fibril Organization |
| +2.275 | 0.001 | Reactome: Platelet Degranulation |
| +2.272 | 0.001 | Reactome: Non-integrin membrane-ECM Interactions |
| +2.236 | 0.001 | Reactome: Response to Elevated Platelet Cytosolic Ca2+ |
| +2.206 | 0.002 | Reactome: Extracellular Matrix Organization |
| +2.194 | 0.002 | GO: Negative Regulation of Angiogenesis |
| +2.181 | 0.003 | Reactome: RHOB GTPase Cycle |
| +2.021 | 0.005 | GO: ERBB Signaling Pathway |
| +1.992 | — | GO: Epidermal Growth Factor Receptor Signaling Pathway |
| +1.941 | 0.008 | Reactome: Signaling By Hippo |
| +1.878 | — | GO: TGF-beta Receptor Signaling Pathway |

### Top Downregulated Pathways (NES < 0, low-Z genes)

| NES | FDR | Pathway |
|-----|-----|---------|
| -2.705 | 0.000 | KEGG: Oxidative phosphorylation |
| -2.646 | 0.000 | GO: Aerobic Electron Transport Chain |
| -2.638 | 0.000 | GO: Cellular Respiration |
| -2.564 | 0.000 | GO: Mitochondrial ATP Synthesis Coupled Electron Transport |
| -2.552 | 0.000 | Reactome: Respiratory Electron Transport / ATP Synthesis / Uncoupling Proteins |
| -2.450 | 0.000 | Reactome: NIK to Noncanonical NF-kB Signaling |
| -2.393 | 0.000 | Reactome: Dectin-1 Mediated Noncanonical NF-kB Signaling |
| -2.386 | 0.000 | Reactome: Degradation of DVL (Wnt) |
| -2.377 | 0.000 | Reactome: Downstream TCR Signaling |
| -2.353 | 0.000 | Reactome: Activation of NF-kappaB in B Cells |
| -2.331 | 0.000 | Reactome: Hh Mutants Abrogate Ligand Secretion |
| -2.330 | 0.000 | Reactome: Hh Mutants Are Degraded by ERAD |
| -2.324 | 0.000 | KEGG: Primary immunodeficiency |
| -2.304 | 0.000 | Reactome: Translation |
| -2.301 | 0.000 | Reactome: ZAP-70 Translocation to Immunological Synapse |

### Biological Interpretation

**Cancer-related pathways — strong signal:**  
Complement and coagulation cascade activation (NES +2.495) is the single strongest upregulated signal and is well-established in PDAC biology. PDAC tumors secrete tissue factor and promote a hypercoagulable state; elevated plasma complement components (C3, C4, factor VII, fibrinogen) are reported in multiple PDAC liquid biopsy studies, including via plasma proteomics and cfRNA. The co-enrichment of Platelet Degranulation (NES +2.275) and Response to Elevated Platelet Cytosolic Ca2+ (NES +2.236) is consistent: PDAC is known to cause paraneoplastic thrombocytosis and platelet activation. These terms in the cfRNA context likely reflect platelet-derived mRNA contamination — platelets package tumor-educated mRNAs during degranulation — but also genuinely reflect systemic coagulation dysregulation in PDAC.

**ECM remodeling — PDAC-characteristic desmoplasia:**  
ECM-receptor interaction (NES +2.390), Extracellular Matrix Organization (NES +2.206), Collagen Fibril Organization (NES +2.299), and Non-integrin membrane-ECM Interactions (NES +2.272) are all significantly upregulated. PDAC is defined by its dense desmoplastic stroma, with type I collagen, fibronectin, and laminin constituting the bulk of the tumor microenvironment. The upregulation of these signatures in plasma cfRNA likely reflects both tumor cell shedding of ECM-associated transcripts and stromal/cancer-associated fibroblast contributions. This is one of the most biologically coherent signals observed across all disease groups.

**Hippo/YAP pathway activation:**  
Hippo Signaling (NES +2.445), Signaling By Hippo (NES +1.941), and YAP1/WWTR1 (TAZ)-stimulated Gene Expression are all upregulated. YAP/TAZ activation is a canonical driver of PDAC progression, stroma activation, and drug resistance. Elevated cfRNA carrying YAP/TAZ target gene signatures is biologically plausible given that tumor cells and cancer-associated fibroblasts both exhibit high YAP activity in PDAC. This signal is notable because Hippo activation is also present in Pancreatitis (NES +2.423), but the PDAC signal is broader (3 vs 2 terms) and includes more downstream effectors.

**RhoGTPase signaling — invasiveness marker:**  
RHOC GTPase Cycle (NES +2.484) is the second strongest upregulated term. RHOC is a master regulator of cytoskeletal reorganization, tumor cell invasion, and metastasis in PDAC. RHOB (NES +2.181) and RHOA (NES +2.102) are also enriched, collectively indicating active cytoskeletal remodeling programs. The Rho GTPase family drives KRAS-mediated oncogenic signaling, which is mutated in >90% of PDAC. This signal is largely unique to Pancreatic Cancer (Moore) and absent from Pancreatitis at the same NES magnitude, making it a candidate discriminatory feature.

**EGFR/ERBB and IGF axis:**  
ERBB Signaling (NES +2.021), EGFR Signaling (NES +1.992), and IGF Transport/IGFBPs (NES +2.400) are upregulated. These growth factor axes are amplified or activated in a substantial fraction of PDACs and are recognized therapeutic targets. The IGF axis in particular has been associated with PDAC desmoplasia and cancer-associated fibroblast crosstalk. TGF-beta Receptor Signaling is present but modest (NES +1.878), consistent with the known dual role of TGF-beta in PDAC (early tumor suppressor, later pro-invasive/immunosuppressive).

**Angiogenesis — a nuanced picture:**  
The enrichment is in "Negative Regulation of Angiogenesis" (NES +2.194) and "Negative Regulation of Blood Vessel Morphogenesis" (NES +2.165), which is counterintuitive at first glance. However, PDAC is paradoxically hypovascular despite having angiogenic drivers. The stroma physically constrains vessel ingrowth; tumor-secreted anti-angiogenic factors (thrombospondin-1, endostatin) mediate this. The simultaneous enrichment of "Regulation of Angiogenesis" (NES +2.157) and "Positive Regulation of VEGF" (NES +2.025) alongside the anti-angiogenic terms captures this tension. Overall the VEGF/endothelial axis is upregulated in cfRNA from PDAC patients, but the tissue biology suppresses productive vascularization.

**Downregulated — mitochondrial/OxPhos suppression:**  
Oxidative phosphorylation (NES -2.705) is the most strongly downregulated pathway, followed by the full hierarchy of mitochondrial electron transport chain complexes (I, II, III, IV), TCA cycle, and ATP synthesis. This is the dominant low-Z signal. In cfRNA biology, this reflects: (1) Warburg effect — PDAC cells reprogram toward aerobic glycolysis, downregulating OxPhos gene expression; (2) A shift in the cellular composition of cfRNA donors — fewer normal epithelial/acinar cells (which are metabolically active) and more tumor/stromal cells. The magnitude (NES -2.705) is the strongest single-term signal in the entire PDAC dataset.

**Downregulated — NF-kB and immune suppression:**  
NIK to Noncanonical NF-kB Signaling (NES -2.450), Dectin-1 Mediated Noncanonical NF-kB Signaling (NES -2.393), Activation of NF-kappaB in B Cells (NES -2.353), and TNFR2 Non-Canonical NF-kB Pathway (NES -2.307) are all strongly downregulated. This is biologically important: PDAC creates an immunosuppressive microenvironment that excludes T cells and B cells and suppresses innate immune signaling. The suppression of canonical and non-canonical NF-kB in cfRNA reflects systemic immune exhaustion and the low lymphocyte cfRNA contribution in PDAC patients. TCR signaling (ZAP-70 Translocation, NES -2.301; Downstream TCR Signaling, NES -2.377; Phosphorylation of CD3/TCR Zeta Chains, NES -2.175) is also uniquely downregulated in PDAC compared to Pancreatitis, consistent with immune exhaustion rather than acute inflammation.

**Downregulated — Hedgehog/Wnt:**  
Multiple Hedgehog terms (Hh Mutants Degraded by ERAD, NES -2.330; Hh Mutants Abrogate Ligand Secretion, NES -2.301; Hedgehog Ligand Biogenesis, NES -2.123) and Wnt terms (Degradation of DVL, NES -2.386; Degradation of Beta-Catenin by Destruction Complex, NES -2.122) are downregulated. The Hedgehog pathway terms refer to the degradation/clearance of mutant Hh proteins — downregulation of these suggests impaired proteasomal clearing of dysfunctional Hh ligands, potentially indicating Hh pathway dysregulation. This is relevant since pancreatic stellate cells respond to Hh ligands from PDAC cells, driving desmoplasia. The Wnt degradation complex terms being low-Z is consistent with stabilized (active) beta-catenin.

**Signal quality: HIGH.** The PDAC (Moore) Z-score GSEA shows a highly coherent biology: (1) ECM/desmoplasia, (2) Hippo/YAP activation, (3) Rho GTPase/invasion, (4) Complement/coagulation, (5) EGFR/ERBB, with a clean mirror of OxPhos suppression and T/B cell immune exhaustion. 554 significant terms with 240 upregulated indicates a broad transcriptome-wide signal. The n=72 provides adequate statistical power.

---

## 2. Pancreatic Cancer (Reggiardo, n=6)

### Summary Statistics
- Total significant terms (FDR < 0.05): **558**
- NES > 0 (upregulated): **117**
- NES < 0 (downregulated): **441**

### Top Upregulated Pathways (NES > 0)

| NES | FDR | Pathway |
|-----|-----|---------|
| +2.218 | 0.000 | GO: Monocyte Chemotaxis |
| +2.115 | 0.003 | KEGG: Oxidative phosphorylation |
| +2.109 | 0.002 | Reactome: Respiratory Electron Transport / ATP Synthesis |
| +2.095 | 0.003 | Reactome: Respiratory Electron Transport |
| +2.057 | 0.005 | Reactome: Vif-mediated Degradation of APOBEC3G |
| +2.057 | 0.005 | GO: Mitochondrial ATP Synthesis Coupled Electron Transport |
| +2.044 | 0.005 | Reactome: APC/C:Cdc20 Mediated Degradation of Mitotic Proteins |
| +2.030 | 0.005 | GO: Cardiac Muscle Cell Contraction |
| +2.024 | 0.005 | Reactome: APC:Cdc20 Mediated Degradation of Cell Cycle Proteins |
| +2.013 | 0.005 | Reactome: Activation of APC/C and APC/C:Cdc20 |
| +2.015 | 0.005 | GO: Eosinophil Migration / Eosinophil Chemotaxis |

### Top Downregulated Pathways (NES < 0)

| NES | FDR | Pathway |
|-----|-----|---------|
| -2.984 | 0.000 | GO: Nucleosome Organization |
| -2.788 | 0.000 | Reactome: SUMOylation of Transcription Cofactors |
| -2.777 | 0.000 | Reactome: RUNX1 Interacts with Co-Factors |
| -2.724 | 0.000 | GO: Regulation of Erythrocyte Differentiation |
| -2.678 | 0.000 | GO: Nucleosome Disassembly |
| -2.614 | 0.000 | Reactome: CDC42 GTPase Cycle |
| -2.612 | 0.000 | Reactome: SUMO E3 Ligases SUMOylate Target Proteins |
| -2.609 | 0.000 | Reactome: Transcriptional Regulation of Granulopoiesis |
| -2.589 | 0.000 | Reactome: Signaling By Cytosolic FGFR1 Fusion Mutants |
| -2.583 | 0.000 | GO: Heterochromatin Formation |
| -2.582 | 0.000 | Reactome: Deadenylation of mRNA |
| -2.578 | 0.000 | Reactome: IL-3, IL-5 and GM-CSF Signaling |
| -2.572 | 0.000 | Reactome: Signaling by Erythropoietin |
| -2.569 | 0.000 | GO: mRNA Metabolic Process |
| -2.566 | 0.000 | GO: Chromatin Remodeling |
| -2.563 | 0.000 | Reactome: SUMOylation |
| -2.556 | 0.000 | GO: Microvillus Assembly |
| -2.550 | 0.000 | Reactome: RAC1 GTPase Cycle |
| -2.532 | 0.000 | Reactome: Signaling to ERKs |
| -2.522 | 0.000 | Reactome: Regulation of TP53 Activity through Acetylation |
| -2.514 | 0.000 | KEGG: Adherens junction |

### Biological Interpretation

**Critical caveat — sample size n=6:**  
With only 6 patients, the normative model Z-scores reflect the mean of an extremely small sample. The GSEA statistics are dominated by the pathway-level aggregation of gene Z-scores rather than patient-to-patient variance, but the biological interpretation must be made with extreme caution. Any pathway signal could be driven by a single outlier patient, batch effects specific to the Reggiardo cohort, or technical artifacts of the COMPLETE-seq platform.

**Reggiardo platform — repeat-element aware cfRNA:**  
Reggiardo et al. used COMPLETE-seq, a hybrid Illumina + Nanopore approach targeting both annotated genes and the full landscape of repetitive elements (Alu, LINE, SINE, endogenous retroviruses). Pancreatic cancer was specifically enriched for Alu subfamilies (AluY, AluSc, AluSg7, AluSc8, AluSx3, AluSg). Critically, the GSEA-prerank used in the normative model is applied to canonical gene-level expression, not the repeat elements themselves. Therefore, the GSEA results here reflect the canonical transcriptomic Z-score profile of Reggiardo's 6 PDAC samples, which may differ from Moore's 72 samples due to: (1) different sample preparation (COMPLETE-seq vs standard cfRNA-seq), (2) Nanopore-specific biases affecting gene quantification, (3) extremely small n, (4) potentially different patient population characteristics.

**Inverted OxPhos signal — a major divergence from Moore:**  
Strikingly, Oxidative phosphorylation (NES +2.115) is *upregulated* in Reggiardo's PDAC samples, the exact opposite of Moore (NES -2.705). This divergence is the most important cross-cohort discrepancy and warrants careful consideration. Possible explanations: (a) The COMPLETE-seq platform may have differential capture efficiency for mitochondrial transcripts (nanopore sequencing captures full-length mitochondrial mRNAs more efficiently than short-read approaches, potentially elevating their apparent abundance); (b) Reggiardo's 6 patients may have been enriched for a metabolically active PDAC subtype; (c) With n=6, a single patient with high mitochondrial gene expression drives the group mean; (d) Technical differences in cfRNA extraction/normalization between the two studies. This OxPhos divergence is a significant red flag for cross-platform generalizability.

**Dominant downregulated signal — chromatin remodeling and SUMOylation:**  
The strongest downregulated terms are Nucleosome Organization (NES -2.984), Nucleosome Disassembly (NES -2.678), Chromatin Remodeling (NES -2.566), Heterochromatin Formation (NES -2.583), and a cascade of SUMOylation terms (NES -2.788 to -2.563). These collectively point to profound epigenetic dysregulation. PDAC is characterized by extensive epigenetic reprogramming, including altered nucleosome positioning, loss of heterochromatin, and SUMO pathway dysregulation. RUNX1 co-factor interactions (NES -2.777) and Transcriptional Regulation of Granulopoiesis (NES -2.609) being downregulated may reflect the low hematopoietic cell cfRNA fraction in these patients (or COMPLETE-seq-specific biases in capturing granulocyte-associated cfRNA). RUNX1 is also expressed in pancreatic progenitors and may reflect tumor-specific downregulation.

**Hematopoietic suppression — erythropoiesis/granulopoiesis:**  
Regulation of Erythrocyte Differentiation (NES -2.724), Signaling by Erythropoietin (NES -2.572), IL-3/IL-5/GM-CSF Signaling (NES -2.578), and Transcriptional Regulation of Granulopoiesis (NES -2.609) are all strongly downregulated. This pattern is consistent with cancer-associated anemia and bone marrow suppression in PDAC. Whether this reflects genuine hematopoietic suppression or reduced hematopoietic cell cfRNA contamination cannot be determined from this data alone.

**GTPase suppression — contrast with Moore:**  
Whereas Moore showed RHOC, RHOB, RHOA upregulation (invasion program), Reggiardo shows CDC42 (NES -2.614), RAC1 (NES -2.550), and many other GTPase cycles as downregulated. This is a second major divergence — Moore's PDAC patients show active invasion signaling while Reggiardo's do not. Again, n=6 and platform differences are the most parsimonious explanation.

**TGF-beta and Wnt — downregulated in Reggiardo:**  
TGF-beta signaling (NES -2.213, -1.847, -1.781) and Wnt pathway terms (Repression of WNT Target Genes, NES -2.476; Formation of beta-catenin:TCF Transactivating Complex, NES -2.322) are downregulated, which is the opposite of Moore for TGF-beta (NES +1.878 upregulated). This inversion likely reflects the platform differences: COMPLETE-seq may capture different cfRNA species than standard total RNA-seq, and with n=6, group-level pathway inference is unreliable.

**Upregulated monocyte chemotaxis:**  
Monocyte Chemotaxis (NES +2.218) is the top upregulated term. PDAC extensively recruits monocytes and tumor-associated macrophages (TAMs). Elevated monocyte chemotactic cfRNA (likely including CCL2, CCL7, CXCL1 derivatives) in 6 patients may reflect genuine immune remodeling in these PDAC patients.

**Signal quality: LOW — use with extreme caution.**  
The combination of n=6, COMPLETE-seq platform specificity, and multiple inversions of established PDAC biology (OxPhos, GTPase, TGF-beta all reversed relative to Moore) makes biological interpretation unreliable for this group. The 558 significant terms with high NES magnitudes are almost certainly artifacts of small-n instability — when 6 patients share an unusual expression profile, the mean Z-score can produce extreme GSEA enrichment even when the signal is unrepresentative of PDAC biology. Results should not be used for biological conclusions without replication in a larger independent COMPLETE-seq cohort.

---

## 3. Pancreatitis (Moore, n=79)

### Summary Statistics
- Total significant terms (FDR < 0.05): **264**
- NES > 0 (upregulated): **30**
- NES < 0 (downregulated): **234**

### Top Upregulated Pathways (NES > 0)

| NES | FDR | Pathway |
|-----|-----|---------|
| +2.423 | 0.000 | GO: Hippo Signaling |
| +2.249 | 0.006 | GO: Positive Regulation of Transcription by RNA Pol III |
| +2.242 | 0.004 | GO: Positive Regulation of Peptidyl-Lysine Acetylation |
| +2.198 | 0.003 | Reactome: RHOC GTPase Cycle |
| +2.196 | 0.002 | GO: Positive Regulation of Focal Adhesion Assembly |
| +2.195 | 0.004 | GO: Positive Regulation of Cell-Substrate Junction Organization |
| +2.179 | 0.003 | GO: Camera-Type Eye Development |
| +2.154 | 0.005 | Reactome: RHOB GTPase Cycle |
| +2.129 | 0.011 | Reactome: Signaling By Hippo |
| +2.116 | 0.010 | Reactome: RUNX3 Regulates NOTCH Signaling |
| +2.082 | 0.017 | KEGG: ECM-receptor interaction |
| +2.064 | 0.024 | GO: Blood Vessel Endothelial Cell Migration |
| +2.027 | 0.033 | GO: Regulation of Epithelial Cell Apoptotic Process |
| +2.024 | 0.035 | Reactome: NOTCH4 Intracellular Domain Regulates Transcription |
| +1.981 | 0.041 | GO: Cellular Response to Glucocorticoid Stimulus |
| +1.974 | 0.044 | GO: Hematopoietic Stem Cell Proliferation |
| +1.965 | 0.049 | Reactome: VEGFR2 Mediated Vascular Permeability |

### Top Downregulated Pathways (NES < 0)

| NES | FDR | Pathway |
|-----|-----|---------|
| -2.658 | 0.000 | Reactome: Respiratory Electron Transport / ATP Synthesis |
| -2.540 | 0.000 | GO: Proton Motive Force-Driven ATP Synthesis |
| -2.537 | 0.000 | KEGG: Oxidative phosphorylation |
| -2.526 | 0.000 | Reactome: Respiratory Electron Transport |
| -2.523 | 0.000 | GO: Oxidative Phosphorylation |
| -2.410 | 0.000 | Reactome: Citric Acid (TCA) Cycle and Respiratory Electron Transport |
| -2.327 | 0.000 | GO: Mitochondrial Electron Transport, NADH to Ubiquinone |
| -2.302 | 0.000 | GO: Mitochondrial Translation |
| -2.295 | 0.000 | GO: Mitochondrial Gene Expression |
| -2.219 | 0.000 | Reactome: Defective CFTR Causes Cystic Fibrosis |
| -2.231 | 0.000 | Reactome: NIK to Noncanonical NF-kB Signaling |
| -2.208 | 0.000 | KEGG: Parkinson disease |
| -2.193 | 0.000 | GO: Immunoglobulin Mediated Immune Response |

### Biological Interpretation

**Expected inflammatory pathways — notably absent or weak:**  
Acute pancreatitis is typically associated with strong innate immune activation: NFkB, IL-6/JAK/STAT3, acute phase response, neutrophil and macrophage infiltration. Strikingly, in the cfRNA Z-score landscape, these canonical inflammatory signals are NOT upregulated. IL-1 Signaling (NES -1.851), IL-4 regulation (NES -1.780), and TNFR2 NF-kB (NES -2.005) are all downregulated — the same direction as in PDAC. Only 30 upregulated terms total. This initially appears paradoxical for an inflammatory disease.

**Biological explanation for low inflammatory signal:**  
The normative model's Z-scores represent deviation from healthy controls. Pancreatitis generates a systemic cfRNA signature where canonical inflammatory gene transcripts (cytokines, NF-kB targets, IL-6 pathway genes) are paradoxically lower relative to HC norms or are not significantly elevated. Two explanations: (1) Chronic pancreatitis (most of Moore's cohort were likely symptomatic referrals for EUS-FNA, including chronic pancreatitis cases) — the inflammatory state may have subsided by blood draw time; (2) cfRNA from plasma in pancreatitis may be dominated by acinar cell death products (enzymatic release) rather than immune cell-derived transcripts. The immune transcripts that would indicate NF-kB/IL-6 activity may be diluted by structural pancreatic cell cfRNA.

**Shared signal with PDAC (Moore) — the overlap problem:**  
21 upregulated terms are shared between Pancreatitis and PDAC (Moore). These include: Hippo Signaling (NES +2.423 in pancreatitis vs +2.445 in PDAC), RHOC GTPase Cycle (NES +2.198 vs +2.484), ECM-receptor interaction (NES +2.082 vs +2.390), Focal Adhesion Assembly, RHOB GTPase Cycle, and Endothelial Cell Migration. This extensive overlap represents a significant diagnostic challenge: the pan-pancreatic disease tissue remodeling response (fibrosis, Hippo activation, RhoGTPase-mediated acinar/ductal restructuring) is indistinguishable between benign and malignant disease at the pathway level in cfRNA.

The 171 shared downregulated terms — dominated by the full OxPhos hierarchy — further underscore that pancreatitis and PDAC produce nearly identical metabolic reprogramming in plasma cfRNA. This convergence is likely driven by the loss of normal acinar cell-derived OxPhos transcripts in both diseases (acinar cell destruction in pancreatitis via autodigestion; acinar cell replacement by tumor/stroma in PDAC).

**Unique upregulated in PDAC (Moore) not seen in Pancreatitis:**  
The key discriminators (present in PDAC, absent from Pancreatitis FDR<0.05 upregulated) are: Complement and Coagulation Cascades (NES +2.495), Post-translational Protein Phosphorylation (NES +2.456), IGF Transport/IGFBPs (NES +2.400), Platelet Degranulation (NES +2.275), Collagen Fibril Organization (NES +2.299), Negative Regulation of Blood Coagulation, ERBB/EGFR Signaling, MAPK/p130Cas-Integrin signaling. These represent the PDAC-specific coagulation/invasion/growth factor axis that pancreatitis does not share.

**Unique downregulated in PDAC (Moore) not in Pancreatitis:**  
Primary immunodeficiency (NES -2.324), TCR signaling (ZAP-70, CD3, downstream TCR; NES -2.175 to -2.301), RNA processing (mRNA capping, polymerase elongation complexes, NER pathways). The T cell suppression signature in cfRNA is unique to PDAC, consistent with the known PDAC immunosuppressive microenvironment, and could serve as a discriminatory biomarker.

**CFTR — a unique signal in Pancreatitis:**  
Defective CFTR Causes Cystic Fibrosis (NES -2.219) is one of the stronger signals unique to Pancreatitis (not significantly enriched in PDAC at equivalent magnitude). CFTR dysfunction is a known risk factor for both acute and chronic pancreatitis; altered CFTR signaling in pancreatic ductal cells causes impaired bicarbonate secretion, protein plug formation, and acinar injury. This CFTR signal in cfRNA may reflect ductal cell-derived cfRNA from injured pancreatic tissue, and is biologically specific to exocrine pancreatic dysfunction.

**Mitochondrial translation — enriched in Pancreatitis relative to PDAC:**  
Pancreatitis shows enrichment of Mitochondrial Translation (NES -2.302), Mitochondrial Translation Initiation/Elongation/Termination (all NES < -2.2), and Mitochondrial Gene Expression (NES -2.295) that is equally strong or stronger than in PDAC. This suggests profound mitochondrial dysfunction in pancreatitis — consistent with the role of mitochondrial injury in pancreatitis pathogenesis (trypsinogen activation is triggered by mitochondrial calcium release and ROS).

**Signal quality: MODERATE.** Pancreatitis shows a coherent metabolic suppression signal (OxPhos, mitochondrial translation). The upregulated pathways are sparse (only 30 terms) but include biologically plausible Hippo, RHOC, and ECM terms. The absence of classical acute-phase/NF-kB upregulation in cfRNA is noteworthy and reflects either chronic disease state or technical limitations of cfRNA in capturing short-lived inflammatory mediators. The CFTR signal is biologically specific. Overall signal density is lower (264 terms) than PDAC (554), consistent with pancreatitis being a less transcriptomically extreme disease state.

---

## 4. Lung Cancer (n=30)

### Summary Statistics
- Total significant terms (FDR < 0.05): **502**
- NES > 0 (upregulated): **44**
- NES < 0 (downregulated): **458**

### Top Upregulated Pathways (NES > 0)

| NES | FDR | Pathway |
|-----|-----|---------|
| +2.726 | 0.000 | GO: Positive Regulation of Integrin-Mediated Signaling Pathway |
| +2.701 | 0.000 | Reactome: Apoptotic Cleavage of Cell Adhesion Proteins |
| +2.596 | 0.000 | GO: Positive Regulation of Glycogen Biosynthetic Process |
| +2.507 | 0.000 | GO: Insulin Receptor Signaling Pathway |
| +2.457 | 0.000 | GO: Blood Vessel Endothelial Cell Migration |
| +2.406 | 0.000 | Reactome: Interleukin-7 Signaling |
| +2.376 | 0.000 | Reactome: GAB1 Signalosome |
| +2.372 | 0.000 | Reactome: Signaling By Hippo |
| +2.354 | 0.000 | GO: Cell Migration Involved in Sprouting Angiogenesis |
| +2.340 | 0.003 | GO: Positive Regulation of Cellular Response to Insulin Stimulus |
| +2.316 | 0.003 | Reactome: SHC1 Events in ERBB2 Signaling |
| +2.313 | 0.003 | Reactome: Molecules Associated with Elastic Fibres |
| +2.215 | 0.006 | Reactome: YAP1- and WWTR1 (TAZ)-stimulated Gene Expression |
| +2.213 | 0.006 | GO: Positive Regulation of Vascular Permeability |
| +2.199 | 0.007 | KEGG: Renal cell carcinoma |
| +2.189 | 0.007 | Reactome: RET Signaling |
| +2.090 | 0.029 | GO: EGFR Signaling Pathway |
| +2.050 | — | GO: TGF-beta Receptor Signaling Pathway |
| +2.037 | — | GO: Hippo Signaling |
| +2.029 | — | Reactome: Signaling By ERBB2 in Cancer |

### Top Downregulated Pathways (NES < 0)

| NES | FDR | Pathway |
|-----|-----|---------|
| -2.372 | 0.000 | KEGG: Ribosome |
| -2.362 | 0.000 | GO: Cytoplasmic Translation |
| -2.361 | 0.000 | Reactome: Eukaryotic Translation Elongation |
| -2.360 | 0.000 | Reactome: Major Pathway of rRNA Processing in Nucleolus and Cytosol |
| -2.357 | 0.000 | Reactome: Formation of A Pool of Free 40S Subunits |
| -2.355 | 0.000 | Reactome: rRNA Processing in Nucleus and Cytosol |
| -2.352 | 0.000 | Reactome: rRNA Processing |
| -2.346 | 0.000 | Reactome: Peptide Chain Elongation |
| -2.345 | 0.000 | Reactome: Eukaryotic Translation Termination |
| -2.342 | 0.000 | Reactome: SRP-dependent Cotranslational Protein Targeting to Membrane |
| -2.331 | 0.000 | GO: Mitochondrial Gene Expression |
| -2.327 | 0.000 | GO: Translation |
| -2.324 | 0.000 | GO: Mitochondrial Translation |
| -2.247 | 0.000 | GO: Proton Motive Force-Driven ATP Synthesis |
| -2.231 | 0.000 | GO: Ribosome Biogenesis |

### Biological Interpretation

**ERBB2 and EGFR signaling — strong upregulation:**  
SHC1 Events in ERBB2 Signaling (NES +2.316), Signaling by ERBB2 in Cancer (NES +2.029), ERBB2 Regulates Cell Motility (NES +2.050), and EGFR Signaling (NES +2.090) are all significantly upregulated. EGFR and HER2 (ERBB2) are among the most clinically relevant targets in non-small cell lung cancer (NSCLC). EGFR mutations are present in 10-35% of NSCLC; HER2 amplification occurs in ~2-4%. The presence of these signals in plasma cfRNA is biologically expected and represents a genuine tumor-derived transcriptomic signal. The GAB1 Signalosome (NES +2.376) — an adapter protein for EGFR and MET signaling — reinforces the EGFR-axis activation. RET Signaling (NES +2.189) is also upregulated; RET rearrangements (RET fusions) occur in ~1-2% of NSCLC and are targetable with selective RET inhibitors.

**Hippo/YAP pathway — consistent cross-cancer signal:**  
Signaling By Hippo (NES +2.372), YAP1/WWTR1 (TAZ)-stimulated Gene Expression (NES +2.215), and GO Hippo Signaling (NES +2.037) are all upregulated. YAP/TAZ activation is established in lung adenocarcinoma and squamous cell carcinoma, driving tumor cell proliferation, invasion, and therapy resistance. This is the third disease group (after PDAC Moore and Pancreatitis) to show Hippo/YAP upregulation, suggesting this is a broadly shared feature of epithelial malignancy/injury captured in plasma cfRNA.

**Angiogenesis and vascular permeability:**  
Blood Vessel Endothelial Cell Migration (NES +2.457), Cell Migration Involved in Sprouting Angiogenesis (NES +2.354), and Positive Regulation of Vascular Permeability (NES +2.213) are all upregulated. Unlike PDAC's hypovascular biology, lung cancer (particularly adenocarcinoma) can be angiogenic. VEGF-mediated angiogenesis drives lung tumor growth, and endothelial cell cfRNA shed from tumor neovasculature would contribute these transcripts to plasma.

**Insulin/IGF and glycogen metabolism:**  
Insulin Receptor Signaling (NES +2.507), Positive Regulation of Cellular Response to Insulin Stimulus (NES +2.340), Positive Regulation of Glycogen Biosynthetic Process (NES +2.596), and Positive Regulation of Glycogen Metabolic Process (NES +2.216) are among the strongest upregulated terms. This metabolic reprogramming — insulin signaling and glycogen anabolism — is consistent with the Warburg-adjacent metabolic rewiring in lung cancer. Lung cancer cells activate insulin/IGF-PI3K-AKT signaling which promotes glycogen synthesis as a metabolic reserve. The IRS1/IRS2 axis driving glycogen synthase (GSK3 inhibition) is a known downstream effector of oncogenic PI3K in NSCLC.

**Integrin signaling and cell adhesion:**  
Positive Regulation of Integrin-Mediated Signaling (NES +2.726, the top term), Molecules Associated with Elastic Fibres (NES +2.313), and Elastic Fibre Formation (NES +2.154) are upregulated. The integrin pathway drives tumor cell invasion and migration; elastic fiber/ECM components are relevant to lung tissue architecture. The Apoptotic Cleavage of Cell Adhesion Proteins (NES +2.701, second top term — including CTNNB1, CASP3, TJP1, TJP2, DSG2) represents an interesting signal: caspase-mediated cleavage of E-cadherin, beta-catenin, and tight junction proteins during tumor cell anoikis/shedding would generate these cfRNA products. CASP3 activation in circulating tumor cells fragmenting into cfRNA is a plausible biological mechanism.

**IL-7 signaling — immune context:**  
Interleukin-7 Signaling (NES +2.406) is upregulated. IL-7 is a cytokine critical for T cell survival and homeostatic proliferation. Its upregulation in lung cancer cfRNA may reflect the immune response attempting to maintain T cell homeostasis in an environment of immune exhaustion — or may represent stromal/tumor production of IL-7 as an immune evasion strategy. This is relevant to the lung cancer immunotherapy context (PD-1/PD-L1).

**Dominant downregulated signal — ribosome and translation:**  
The most striking and dominant low-Z signal in Lung Cancer is the complete suppression of ribosome biogenesis and cytoplasmic translation: KEGG Ribosome (NES -2.372) and a cascade of 35+ ribosome/translation terms (rRNA processing, eukaryotic translation elongation/termination, SRP-dependent targeting, cytoplasmic translation). This is qualitatively different from the OxPhos-dominant downregulation seen in PDAC and Pancreatitis. Ribosomal protein genes are among the most highly expressed genes in healthy cells; their suppression in cfRNA from lung cancer patients likely reflects: (1) A compositional shift — tumor-derived cfRNA has a different gene expression profile than normal epithelial cfRNA, with lower ribosomal protein mRNA abundance relative to the healthy reference; (2) Active repression of ribosome biogenesis in the disease state. The Chen et al. (eLife 2022) study noted that PI3K-AKT and longevity/metabolic pathways were associated with NSCLC in their cfRNA dataset; the ribosome suppression seen here may reflect a complementary global translation suppression in cancer.

**OxPhos suppression — present but secondary:**  
Proton Motive Force-Driven ATP Synthesis (NES -2.247) and Oxidative Phosphorylation (NES -2.124) are downregulated, similar to PDAC and Pancreatitis but with lower NES magnitude. In lung cancer, ribosome suppression dominates over OxPhos suppression, unlike pancreatic diseases. This distinction may be disease-tissue-specific.

**Hedgehog — also suppressed:**  
Multiple Hh terms (Signaling by Hedgehog NES -1.589, Hedgehog Ligand Biogenesis NES -1.930, Hh Mutants NES -2.017) are downregulated, similar to PDAC.

**cfRNA artifact signals:**  
The 44 upregulated terms in Lung Cancer include "Cardiac Muscle Cell Contraction/Action Potential" terms (absent from PDAC), and "Na+/Cl- Dependent Neurotransmitter Transporters" (NES +2.078). These are likely non-specific signals from circulating cardiomyocyte/neuronal cfRNA not specific to lung cancer biology, or may reflect cross-reactivity of gene set annotations. These should be treated as potential artifacts.

**Signal quality: MODERATE-HIGH.** Lung cancer shows coherent ERBB2/EGFR, Hippo/YAP, angiogenesis, and insulin/IGF signals in the upregulated set, all biologically expected. The dominant ribosome suppression is a characteristic feature distinguishing it from PDAC/Pancreatitis. The smaller upregulated gene set (44 terms) relative to the large downregulated set (458 terms) suggests the disease signature is primarily expressed as a suppression of normal transcriptional programs. With n=30, statistical power is moderate; replication in a larger lung cancer cohort would be needed to confirm the specificity of individual pathway terms.

---

## 5. Comparative Analysis: Pancreatic Cancer (Moore) vs. Pancreatitis

This comparison is particularly informative as both groups come from the same study (Moore et al. 2025) and the same clinical setting (pre-EUS-FNA), eliminating most technical confounders.

### Discriminating Features

**Shared (not discriminatory) — both diseases:**
- OxPhos/mitochondrial suppression: 171 shared downregulated terms dominated by electron transport chain
- Hippo/YAP upregulation
- RHOC and RHOB GTPase upregulation
- ECM-receptor interaction
- Focal Adhesion Assembly
- NF-kB suppression (all noncanonical NF-kB terms downregulated in both)
- Hedgehog pathway suppression
- Global translation reduction

This shared core suggests a pan-pancreatic disease signature driven by acinar cell loss, fibrosis, and ductal remodeling, present regardless of benign vs. malignant etiology.

**Unique to PDAC (Moore) — discriminatory upregulated signals:**
- Complement and Coagulation Cascades (NES +2.495) — ABSENT in Pancreatitis
- Post-translational Protein Phosphorylation / IGFBPs (NES +2.456, +2.400)
- Platelet Degranulation and Platelet Ca2+ response
- Collagen Fibril Organization
- Non-integrin membrane-ECM Interactions
- ERBB/EGFR Signaling axis
- Scavenging by Class A Receptors
- Rho GTPase cycle broader: RHOA, RHOV, RND3 (beyond RHOB/RHOC shared with pancreatitis)

**Unique to PDAC (Moore) — discriminatory downregulated signals:**
- TCR signaling: ZAP-70 translocation, CD3/TCR Zeta Chains phosphorylation, Downstream TCR Signaling, Primary immunodeficiency
- RNA elongation: Formation of RNA Pol II Elongation Complex, Early Elongation Complex, mRNA Capping
- DNA repair: Dual Incision in TC-NER, Incision Complex in GG-NER
- rRNA modification

**Unique to Pancreatitis — discriminatory signals:**
- Defective CFTR Causes Cystic Fibrosis (NES -2.219) — biologically specific to pancreatic exocrine dysfunction
- NOTCH signaling (NOTCH4, RUNX3/NOTCH)
- Parkinson disease KEGG pathway (likely mitochondrial complex I overlap, not disease-specific)
- Hematopoietic Stem Cell Proliferation (NES +1.974)
- Glucocorticoid stimulus response (NES +1.981)

### Discriminability Assessment

At the pathway level, PDAC and Pancreatitis are **partially overlapping but distinguishable**. The most diagnostically useful discriminatory features in cfRNA are:

1. **Complement/coagulation activation** — strongly PDAC-specific
2. **T cell suppression (TCR/ZAP-70)** — strongly PDAC-specific  
3. **CFTR pathway** — Pancreatitis-specific
4. **EGFR/ERBB amplification** — PDAC-enriched
5. **Platelet activation** — PDAC-enriched

This aligns with Moore et al.'s finding that liver function biomarkers and coagulation-related cfRNA discriminate PDAC from pancreatitis.

---

## 6. cfRNA Background Artifacts and Technical Considerations

### Platelet-derived cfRNA contamination
Platelet Degranulation (NES +2.275 in PDAC) and Platelet Ca2+ Response (NES +2.236) represent known cfRNA artifacts. Platelets are a major source of mRNA in plasma, and incomplete platelet removal during plasma preparation can lead to platelet mRNA contamination. However, in PDAC, platelet activation is also a genuine biological phenomenon (tumor-educated platelets). Careful interpretation is required — the signal is present but cannot be attributed purely to biology without matched platelet-depleted cfRNA analysis.

### Ribosomal/translation signal in Lung Cancer
The extreme downregulation of ribosomal proteins (KEGG Ribosome NES -2.372) in Lung Cancer likely includes a technical component: ribosomal protein mRNAs are among the most abundant transcripts in most cells. If the cfRNA extraction and sequencing depth for Chen et al.'s samples differs from the HC normative population, ribosomal protein depletion/enrichment steps could produce artifactual enrichment of this pathway.

### Mitochondrial RNA leakage
In all four disease groups, mitochondrial OxPhos/translation terms are strongly downregulated. While biologically interpretable, it is also possible that plasma contains a fraction of mitochondrial RNA derived from apoptotic cells; altered apoptosis rates between disease and healthy control populations could contribute to these signals. The HC normative model trained on 693 HCs establishes the baseline mitochondrial cfRNA level; any disease-induced change in apoptosis rates would register as a Z-score deviation.

### "Camera-type Eye Development" and "Neuron Migration" terms
These non-disease-relevant terms appear with high NES in PDAC (Camera-type Eye Development NES +2.282) and are shared with Pancreatitis. These are likely driven by shared gene sets (laminin, nidogen, integrin components) with ECM pathways, where the GO term annotation captures the same molecules in a developmental context. They are not biologically meaningful in the pancreatic cancer context and represent annotation overlap artifacts.

---

## 7. Overall Signal Quality Assessment per Disease

| Disease | n | Sig. terms | Up/Down | Quality | Key signals | Caveats |
|---------|---|------------|---------|---------|-------------|---------|
| Pancreatic Cancer (Moore) | 72 | 554 | 240/314 | **HIGH** | ECM/desmoplasia, Hippo/YAP, RHOC invasion, Complement/coagulation, ERBB/EGFR; OxPhos suppression, T cell exhaustion | Platelet contamination possible; Hh/Wnt suppression complex |
| Pancreatic Cancer (Reggiardo) | 6 | 558 | 117/441 | **LOW** | OxPhos upregulation (REVERSED), Chromatin/nucleosome/SUMO suppression, Hematopoietic suppression | n=6 — unreliable; COMPLETE-seq platform incompatibility; multiple biology inversions vs Moore |
| Pancreatitis (Moore) | 79 | 264 | 30/234 | **MODERATE** | Hippo/YAP, RHOC, ECM; OxPhos/mitochondrial suppression; CFTR-specific signal | Sparse upregulation; no classic acute-phase signal in cfRNA; highly overlapping with PDAC |
| Lung Cancer | 30 | 502 | 44/458 | **MODERATE-HIGH** | ERBB2/EGFR, Hippo/YAP, Angiogenesis, Insulin/glycogen metabolism; Ribosome suppression dominant | n=30 moderate power; ribosome suppression may have technical component; cardiac artifact terms |

---

## 8. Key Conclusions

1. **PDAC Moore (n=72) shows the most biologically coherent cfRNA signature** among the four groups, with canonical PDAC hallmarks (desmoplasia, complement, Hippo, RHOC invasion, EGFR) all represented. This is a high-quality signal with clear mechanistic interpretability.

2. **Pancreatitis and PDAC share the majority of their pathway signatures** at the cfRNA level, consistent with the well-known clinical challenge of distinguishing these conditions. The discriminatory features (complement/coagulation, TCR suppression for PDAC; CFTR for pancreatitis) align with previously reported protein and cfDNA biomarker literature.

3. **Lung Cancer shows a distinct dominant signal: ribosome/translation suppression** rather than OxPhos suppression. The upregulated set is coherent (ERBB2, Hippo, angiogenesis, insulin signaling) and tissue-appropriate. This disease group is most distinguishable from the pancreatic group at the pathway level.

4. **Reggiardo PDAC (n=6) results should not be interpreted biologically.** The platform-specific (COMPLETE-seq) and n=6 instability produce a pathway profile that inverts several well-established PDAC biological findings relative to Moore's 72-patient cohort from the same disease. The chromatin/SUMO/nucleosome suppression dominant signal likely reflects cohort-specific characteristics not generalizable to PDAC.

5. **OxPhos suppression is a shared pan-disease signature** in plasma cfRNA across pancreatic cancer, pancreatitis, and lung cancer, limiting its diagnostic specificity. Its use as a disease biomarker requires comparison against disease-specific upregulated terms.

6. **Hippo/YAP upregulation is a consistently detected signal** across all three disease groups with coherent GSEA results (PDAC, Pancreatitis, Lung Cancer), suggesting it represents a broadly shared feature of epithelial injury/malignancy detectable in plasma cfRNA. It should not be used as a disease-specific marker.

---

*Evaluation written based on GSEA-prerank results from the GAMLSS normative model (HC n=693). Literature references: Moore et al. Nat Commun 2025 (PMID: 40783559); Reggiardo et al. Nat Biomed Eng 2023 (PMID: 37652985); Chen et al. eLife 2022 (DOI: 10.7554/eLife.75181).*
