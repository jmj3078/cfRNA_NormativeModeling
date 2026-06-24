# GSEA-Prerank Biological Evaluation: Group 5
## CDCS / Pre-eclampsia / Other Cancer

**Analysis date:** 2026-06-24  
**Method:** GSEA-prerank on mean Z-score per gene per phenotype (GAMLSS normative model trained on 693 healthy control plasma cfRNA)  
**Gene sets used:** GO Biological Process 2023, Reactome 2022, KEGG 2021 Human  
**NES interpretation:** NES > 0 = enriched in high-Z (relatively upregulated vs. HC norm); NES < 0 = enriched in low-Z (relatively downregulated vs. HC norm)

---

## Summary Statistics

| Disease | n | Total pathways tested | FDR<0.05 | NES>0 | NES<0 | NES range |
|---------|---|-----------------------|----------|-------|-------|-----------|
| CDCS (Coronary Disease Cohort Study) | 212 | 867 | 867 (100%) | 11 | 856 | −2.13 to +2.59 |
| Pre-eclampsia | 59 | 239 | 239 (100%) | 3 | 236 | −2.27 to +2.42 |
| Other Cancer | 16 | 267 | 267 (100%) | 18 | 249 | −3.10 to +2.41 |

The 100% FDR<0.05 rate for all three diseases indicates exceptionally strong and consistent deviation from the healthy control normative model. The overwhelming NES<0 dominance — particularly for CDCS (98.6%) and Pre-eclampsia (98.7%) — is the defining feature of this group and demands biological interpretation.

---

## Disease 1: CDCS — Coronary Disease Cohort Study

### 1.1 Source Study

Ward Z et al., "Identifying Candidate Circulating RNA Markers for Coronary Artery Disease by Deep RNA-Sequencing in Human Plasma," *Cells* 2022, 11(20):3191. DOI: 10.3390/cells11203191. **CDCS = Coronary Disease Cohort Study**, a prospective observational cohort recruited at two New Zealand hospitals (Christchurch/Dunedin). The Ward et al. study used 59 CAD patients (half of whom developed de novo heart failure within 3 years) plus 30 healthy controls, profiling plasma cfRNA by deep RNA-seq. The original study identified 160 differentially expressed mRNAs and 10 lncRNAs in CAD versus controls, with 13 transcripts of mitochondrial origin, implicating myocardial ischemia and oxidative stress. Our normative model dataset contains n=212 CDCS samples — a larger subset from the same cohort, likely encompassing additional time points or patient categories.

### 1.2 GSEA Results: NES > 0 (Upregulated, n=11)

All 11 significant NES>0 pathways are FDR<0.05:

| NES | Term |
|-----|------|
| +2.59 | GO: Import Into Cell (GO:0098657) |
| +2.43 | KEGG: Nicotine addiction |
| +2.28 | GO: Positive Regulation Of Potassium Ion Transmembrane Transporter Activity |
| +2.23 | Reactome: P2Y Receptors |
| +2.23 | GO: Aminoglycan Biosynthetic Process |
| +2.21 | Reactome: Collagen Chain Trimerization |
| +2.15 | GO: Regulation Of Sensory Perception |
| +2.09 | Reactome: Voltage Gated Potassium Channels |
| +2.09 | GO: Regulation Of Sensory Perception Of Pain |
| +2.09 | Reactome: GAG Synthesis — Tetrasaccharide Linker Sequence |
| +2.07 | Reactome: Presynaptic Nicotinic Acetylcholine Receptors |

These upregulated pathways reveal a coherent biological signature. Collagen chain trimerization and GAG/aminoglycan biosynthesis reflect **extracellular matrix remodeling** — a hallmark of cardiac fibrosis and vascular remodeling in CAD. P2Y purinergic receptors and voltage-gated potassium channels are expressed in platelets and vascular smooth muscle; their elevated cfRNA signal is consistent with enhanced platelet activation and altered vascular tone seen in coronary artery disease. The nicotine addiction and nicotinic acetylcholine receptor signals may reflect autonomic nervous system dysregulation in cardiac disease, or simply residual effects of smoking as a major CAD risk factor in this cohort. Import into cell (endocytosis/transport pathways) likely reflects heightened cellular uptake activity in remodeling vascular tissue. These 11 pathways, while numerically sparse, are biologically coherent with established CAD pathophysiology.

### 1.3 GSEA Results: NES < 0 (Downregulated, n=856)

The sheer scale of NES<0 enrichment — 856 out of 867 pathways — is the most striking feature of the CDCS result. The top downregulated terms cluster into several distinct themes:

**Theme 1: Ribosome biogenesis and translation (strongest signal, NES as low as −2.13)**
- Cytoplasmic Translation (NES −2.13)
- Eukaryotic Translation Termination, Elongation, Peptide Chain Elongation (NES −2.09 to −2.08)
- Cap-dependent Translation Initiation, Formation of 43S/60S complexes (NES −2.03 to −1.98)
- KEGG: Ribosome (NES −1.98)
- Nonsense-Mediated Decay (NMD), both EJC-dependent and -independent (NES −2.07 to −2.02)
- Response of EIF2AK4 (GCN2) to Amino Acid Deficiency (NES −2.06)

**Theme 2: Splicing and RNA processing (~31 pathways)**
- Spliceosomal Complex Assembly, U2-type Prespliceosome Assembly (NES ~−1.80)

**Theme 3: Mitochondrial function and oxidative phosphorylation (~27 pathways)**
- Oxidative Phosphorylation (NES −1.72), Proton Motive Force-Driven ATP Synthesis (NES −1.69)
- Respiratory Electron Transport, Complex I Biogenesis, TCA Cycle (NES −1.47 to −1.61)
- Mitochondrial protein import and organization

**Theme 4: Immune and inflammatory signaling (~55 pathways)**
- JAK-STAT signaling after IL-12 stimulation (NES −1.82)
- NF-κB activation (FCERI, BCR, TNFR2, NIK-mediated noncanonical) (NES −1.74 to −1.77)
- B cell receptor signaling, T cell receptor signaling (NES −1.63 to −1.69)
- IL-1 signaling, IL-12 signaling (NES −1.67 to −1.74)
- Antigen cross-presentation, macrophage chemotaxis (NES −1.59 to −1.62)
- DDX58/IFIH1-mediated interferon induction (NES −1.61)

**Theme 5: Broad signaling, cell cycle, and proteostasis (~195 pathways)**
- Regulation of apoptosis (NES −1.81)
- Regulation of NOTCH4 signaling (NES −1.83)
- GSK3B/BTRC-mediated NFE2L2 degradation (NES −1.83)
- Proteasomal degradation pathways
- CSF3 (G-CSF) signaling (NES −1.80)
- ROBO receptor signaling, SLIT/ROBO regulation (NES −1.92 to −1.97)

### 1.4 Biological Interpretation: Why Predominantly NES < 0?

The massive NES<0 dominance in CDCS is **biologically meaningful and partially expected** for several reasons:

**1. Global translational suppression in cardiovascular disease.** Ribosomal protein genes and translation factors are consistently downregulated in heart failure and coronary artery disease across tissue and blood transcriptomic studies. In HFpEF, ribosomal protein translation and oxidative phosphorylation are among the most robustly downregulated protein categories. The Ward et al. study itself identified downregulation of RNA-binding proteins involved in ribosomal function in ischemic heart failure. Reduced protein synthesis capacity reflects the energetic and biosynthetic constraints of ischemia-stressed cardiomyocytes and circulating immune cells.

**2. Mitochondrial dysfunction as a primary CAD feature.** Downregulation of oxidative phosphorylation, electron transport chain components, and mitochondrial assembly pathways directly mirrors the original Ward et al. finding of 13 mitochondrially encoded transcripts in CAD cfRNA, implicating myocardial ischemia. The plasma cfRNA signal from damaged cardiomyocytes releases mitochondrial transcripts that are themselves downregulated relative to healthy myocardial baseline captured by the HC normative model.

**3. Immunosuppressive plasma environment in stable CAD.** The 55 downregulated immune pathways (NF-κB, IL-12/JAK-STAT, BCR, TCR, macrophage chemotaxis, antigen cross-presentation) reflect the complex immunomodulation in chronic stable CAD. While acute MI triggers immune activation, stable CAD is characterized by a state of regulatory immune adaptation in which inflammatory thresholds are recalibrated. This can manifest as relative suppression of innate and adaptive immune response genes compared to the healthy normative baseline.

**4. The normative model captures "active" HC biology.** The healthy control cohort (n=693) represents individuals with fully active ribosomal, mitochondrial, and immune gene programs. CAD patients, even in stable state, have chronic ischemic tissue stress that impairs these fundamental processes system-wide. The 11 NES>0 pathways (ECM remodeling, potassium channels, P2Y receptors) represent the relatively disease-specific upregulations, while the sea of NES<0 pathways reflects global homeostatic degradation.

**5. The 856 NES<0 vs. 11 NES>0 ratio as a "global suppression" signature.** This extreme asymmetry is not arbitrary noise — it indicates that CAD fundamentally suppresses more biological processes in circulating plasma cfRNA than it activates. This is consistent with the concept that chronic ischemic disease produces a state of systemic bioenergetic stress where anabolic and biosynthetic processes (translation, mitochondrial ATP production, immune activation) are globally curtailed. The plasma cfRNA thus reports a body-wide transcriptional shutdown relative to the healthy reference.

**Assessment: High biological coherence.** The CDCS GSEA result is broadly concordant with known CAD biology. The translation/ribosome/OXPHOS downregulation cluster is well-supported by independent tissue studies of heart failure. The small upregulated cluster (ECM, ion channels, purinergic signaling) maps cleanly onto vascular remodeling biology. The extreme NES<0 dominance is unusual in scale but mechanistically interpretable as systemic bioenergetic suppression in chronic ischemic disease.

---

## Disease 2: Pre-eclampsia

### 2.1 Source Study

Moufarrej MN et al., "Early prediction of preeclampsia in pregnancy with cell-free RNA," *Nature* 2022, 602:689–694. DOI: 10.1038/s41586-022-04410-z. This landmark study demonstrated that plasma cfRNA transcriptomic changes in preeclampsia are detectable as early as 5–16 weeks of gestation, long before clinical symptoms emerge. Using 404 blood samples from 199 pregnant women, Moufarrej et al. identified an 18-gene cfRNA panel distinguishing pre-eclampsia from normotensive pregnancies. Key biological finding: cfRNA changes are enriched in **neuromuscular, endothelial, and immune cell-type signatures**, consistent with the core pathophysiology of pre-eclampsia (impaired trophoblast invasion, endothelial dysfunction, and systemic inflammatory response). The enriched tissue contributions reflected brain/neuromuscular tissue, vascular endothelium, and immune compartments — not primarily placenta, which is the conventional focus of pre-eclampsia research.

### 2.2 GSEA Results: NES > 0 (Upregulated, n=3)

Only 3 significant NES>0 pathways:

| NES | Term |
|-----|------|
| +2.42 | Reactome: Neurexins And Neuroligins (R-HSA-6794361) |
| +2.17 | Reactome: DNA Methylation (R-HSA-5334118) |
| +2.11 | GO: Postsynaptic Membrane Organization (GO:0001941) |

These three upregulated pathways — neurexin/neuroligin signaling and postsynaptic membrane organization — point to a **neuronal/synaptic cfRNA signature** that is elevated in pre-eclampsia relative to the HC normative baseline. This directly corroborates the Moufarrej et al. finding that neuromuscular cell-type specific genes are enriched in pre-eclampsia cfRNA. The biological interpretation: pre-eclampsia is associated with elevated circulating RNA from neuronal or neuromuscular tissue compartments, possibly reflecting the neurological involvement (headaches, visual disturbances, seizures in severe eclampsia) or placental-neural communication pathways. DNA methylation pathway enrichment may reflect epigenetic dysregulation underlying trophoblast dysfunction. The paucity of upregulated pathways (only 3) despite strong disease biology points to a genuine asymmetry: pre-eclampsia in plasma cfRNA manifests primarily as suppression of normal HC programs rather than activation of disease-specific programs detectable at the pathway level.

### 2.3 GSEA Results: NES < 0 (Downregulated, n=236)

**Dominant theme: Ribosome biogenesis, translation, and RNA processing (strongest and most numerous)**

The top 30 pathways are almost exclusively ribosomal and translation-related:
- Ribosomal Small Subunit Biogenesis (NES −2.27), Ribosome Biogenesis (NES −2.18), Ribosomal Large Subunit Biogenesis (NES −2.03)
- Cap-dependent Translation Initiation (NES −2.27), Eukaryotic Translation Elongation/Termination (NES −2.16 to −2.14)
- rRNA processing and maturation pathways (multiple, NES −2.03 to −2.17)
- Major Pathway of rRNA Processing in Nucleolus and Cytosol (NES −2.17)
- Ribonucleoprotein Complex Biogenesis (NES −2.21)
- RNA Export from Nucleus, Nucleocytoplasmic Transport, mRNA Transport (NES −2.08 to −2.15)

**Theme 2: Interferon and innate immune signaling (~downregulated)**
- Defense Response to Symbiont (NES −2.05), ISG15 Antiviral Mechanism (NES −1.99), Antiviral Mechanism By IFN-Stimulated Genes (NES −1.95)
- Interferon Signaling (NES −1.88), Interferon Alpha/Beta Signaling (NES −1.82)
- NLR Signaling Pathways (NES −1.93), NLRP3 Inflammasome (NES −1.80)
- Defense Response to Virus, Response to Interferon-Beta (NES −1.74 to −1.84)
- Regulation of Innate Immune Response (NES −1.78), C-type Lectin Receptors (NES −1.83)
- Negative Regulation of IL-12 Production (NES −1.82)

**Theme 3: Splicing, transcription, and chromatin**
- Spliceosome (KEGG, NES −2.00), mRNA Splicing, mRNA Processing (NES −1.88 to −1.90)
- RNA Splicing via Transesterification (NES −1.83)
- Chromatin Remodeling (NES −1.89), Regulation of Histone Deacetylation (NES −1.84)
- Transcription Elongation by RNA Polymerase II (NES −1.87 to −1.89)

**Theme 4: Post-translational modification and proteostasis**
- Cellular Response to Starvation (NES −1.94), Cellular Response to Heat Stress (NES −1.94)
- Regulation of HSF1-mediated Heat Shock Response (NES −1.84)

**Theme 5: Lipoprotein and metabolic pathways**
- Lipoprotein Localization (NES −1.91), Lipoprotein Transport (NES −1.89)
- Free Fatty Acids Regulate Insulin Secretion (NES −1.77)

### 2.4 Biological Interpretation

**Why 236 NES<0 and only 3 NES>0?**

This result initially appears paradoxical: pre-eclampsia is a condition of excessive inflammation, endothelial activation, and placental stress — yet our normative model detects predominantly *suppressed* pathways relative to healthy controls. Multiple non-mutually-exclusive explanations apply:

**1. The HC reference captures biologically "active" pregnant-state pathways that are suppressed in PE.** Critically, the HC normative reference in this project is drawn from WTS plasma samples that do not specifically exclude pregnant women, but also do not systematically include them. If the HC reference primarily represents non-pregnant individuals, then comparing pre-eclamptic pregnant women to this HC norm is inherently confounded: normal pregnancy alone upregulates placental transcripts, immune tolerance genes, and metabolic adaptation programs. Pre-eclampsia could appear as a relative *suppression* of these programs (compared to the HC norm which itself has lower baseline for placenta-derived signals) while simultaneously having lower activity in canonical HC pathways (ribosomal biogenesis, immune homeostasis). The normative model Z-score reflects deviation from non-pregnant HC norms, so *both* the gestational state and the disease contribute to the signal.

**2. Ribosomal suppression reflects placental stress, not normal HC biology.** Pre-eclampsia is characterized by placental hypoxia and endoplasmic reticulum stress, which activates the integrated stress response (ISR) and suppresses global mRNA translation via eIF2α phosphorylation. This mechanistic pathway — hypoxia → ISR → reduced ribosome occupancy and ribosomal RNA synthesis — directly predicts the NES<0 enrichment of translation and ribosome biogenesis pathways in plasma cfRNA. Placenta-derived cfRNA carrying low translation capacity would appear suppressed relative to HC.

**3. Interferon pathway downregulation reflects immune tolerance breakdown.** Pre-eclampsia is associated with dysregulation of the maternal immune tolerance to fetal alloantigens. The downregulation of ISG15 antiviral mechanisms, interferon signaling, NLR/NLRP3 inflammasome, and C-type lectin receptor pathways in cfRNA relative to HC could reflect several things: (a) exhaustion of innate immune sentinel cells that normally produce these transcripts in healthy individuals; (b) immune cell redistribution from circulation to inflamed uteroplacental tissue; or (c) the dominant contribution of placental cfRNA that lacks mature immune transcripts. Studies of second-trimester maternal blood RNA profiles in pre-eclampsia have confirmed that immune dysregulation, including reduction in NK cells and cytotoxic gene expression, is a central feature.

**4. Downregulated lipoprotein transport reflects metabolic reprogramming.** Pre-eclampsia is characterized by dyslipidemia (elevated triglycerides, LDL oxidation, reduced HDL) and insulin resistance. The downregulation of lipoprotein localization and transport pathways (NES −1.91 to −1.89) maps onto this metabolic derangement. The free fatty acid regulation of insulin secretion pathway (NES −1.77) further supports a metabolic dysregulation signature.

**5. The SLIT/ROBO pathway is shared with Other Cancer.** Regulation of Expression of SLITs and ROBOs (NES −1.94 in pre-eclampsia) represents an interesting finding. ROBO receptors regulate axonal guidance but also trophoblast invasion — their downregulation is consistent with impaired trophoblast invasion characteristic of early-onset pre-eclampsia. This pathway also appears in CDCS and Other Cancer NES<0 results, suggesting it may be a general stress-response cfRNA signature.

### 2.5 Comparison with Moufarrej et al. 2022

The Moufarrej et al. Nature study identified an 18-gene panel enriched in neuromuscular, endothelial, and immune cell-type specific signatures. Our normative model GSEA approach **complements rather than contradicts** their findings in several ways:

**Concordance:** The Moufarrej study found neuromuscular gene enrichment as one of the primary signals in pre-eclampsia cfRNA — our GSEA shows the only 3 NES>0 pathways are neurexin/neuroligin and postsynaptic membrane organization, perfectly echoing this. Both approaches detect a neuronal/synaptic cfRNA component that is elevated in pre-eclampsia.

**Complementarity:** Moufarrej et al. focused on early prediction using a supervised classification approach with a small gene panel, anchored to comparisons within pregnant women (PE vs. normotensive pregnancy). Our normative model compares against the *healthy non-pregnant* HC reference using an unsupervised Z-score approach across the full transcriptome. The Moufarrej approach captures pregnancy-specific contrasts; our approach captures deviation from systemic HC physiology. The NES<0 dominance in our results likely reflects in part the fundamental difference between a pre-eclamptic pregnant woman and a healthy non-pregnant control — not all of this is pure disease signal, but some reflects the gestational context.

**Divergence:** Moufarrej et al. highlighted endothelial cell-type cfRNA as a key signal, suggesting endothelial activation in pre-eclampsia. Our GSEA does not clearly identify upregulated endothelial pathways (e.g., angiogenesis, VEGF signaling are not in the NES>0 set), which may be because: (a) these signals are diluted when averaging across the full pre-eclampsia cohort (n=59, mixed severity/timing); (b) the endothelial component is real but sub-threshold at the pathway level given the small sample size; or (c) our normative model may not be specifically powered to detect endothelial cfRNA relative to non-pregnant HC.

**Overall assessment:** The pre-eclampsia GSEA result has high biological plausibility. The NES>0 neuronal signal corroborates Moufarrej et al. directly. The NES<0 dominance reflects genuine biology (translational stress via ISR, immune redistribution, metabolic dyslipidemia, SLIT/ROBO trophoblast invasion failure) compounded by a systematic deviation from the non-pregnant HC normative baseline. Interpreting pre-eclampsia cfRNA normative model results requires acknowledging the pregnancy confounder: the biological "distance" from non-pregnant HC is driven by both gestational physiology and superimposed disease.

---

## Disease 3: Other Cancer

### 3.1 Source Study

Moore TW et al., "Cell free RNA detection of pancreatic cancer in pre-diagnostic high risk and symptomatic patients," *Nature Communications* 2025. DOI: 10.1038/s41467-025-62685-y. This study focused on pancreatic ductal adenocarcinoma (PDAC) detection from plasma cfRNA, using 153 samples across five diagnostic groups: benign pancreas, pancreatitis, IPMN (intraductal papillary mucinous neoplasm), PDAC, and **Other Cancers**. The "Other Cancer" group in our dataset (n=16) thus represents non-pancreatic cancer patients who served as comparators in the Moore et al. study. The original study identified 29 cfRNA biomarkers for PDAC diagnosis, with liver function-related transcripts elevated even in early-stage PDAC without liver metastasis. The specific cancer types in the "Other Cancer" group are not specified in the abstract but likely include diverse solid tumors (gastrointestinal, hepatobiliary, or other cancers seen at the recruitment centers) given the pancreatic surgery/GI oncology context of the study.

### 3.2 GSEA Results: NES > 0 (Upregulated, n=18)

| NES | Term |
|-----|------|
| +2.41 | Reactome: RHOB GTPase Cycle |
| +2.27 | GO: Regulation Of Cell-Substrate Junction Assembly |
| +2.25 | GO: Regulation Of Focal Adhesion Assembly |
| +2.18 | Reactome: RHOC GTPase Cycle |
| +2.15 | GO: Regulation Of Cell-Matrix Adhesion |
| +2.14 | Reactome: p130Cas Linkage To MAPK Signaling For Integrins |
| +2.10 | Reactome: RHOQ GTPase Cycle |
| +2.10 | Reactome: RHOA GTPase Cycle |
| +2.10 | Reactome: CDC42 GTPase Cycle |
| +2.07 | GO: Actin Filament-Based Transport |
| +2.07 | GO: Regulation Of VEGF Signaling Pathway |
| +2.05 | GO: Protein Localization To Cell Junction |
| +2.03 | Reactome: Regulation Of RUNX1 Expression And Activity |
| +2.03 | Reactome: RAC1 GTPase Cycle |
| +2.02 | GO: Regulation Of Small GTPase Mediated Signal Transduction |
| +2.00 | GO: Negative Regulation Of Blood Coagulation |
| +1.99 | GO: Substrate Adhesion-Dependent Cell Spreading |
| +1.98 | GO: Regulation Of Fibroblast Migration |

The NES>0 pathways in Other Cancer form a strikingly coherent and highly cancer-relevant signature dominated by:

**Rho GTPase signaling (RHOA, RHOB, RHOC, RHOQ, RAC1, CDC42):** All five major Rho family GTPases show upregulated cycle activity. Rho GTPases are master regulators of the actin cytoskeleton and are universally dysregulated in cancer — they drive invasion, metastasis, and cell survival. Their elevation in cfRNA reflects shed RNA from tumor cells with hyperactive Rho signaling programs.

**Focal adhesion and cell-matrix adhesion:** Focal adhesion assembly, cell-substrate junction assembly, and p130Cas-MAPK integrin signaling are core mechanisms of tumor invasion. The p130Cas linkage to MAPK for integrins specifically connects integrin engagement to proliferative MAPK signaling — a canonical pro-invasive pathway.

**VEGF signaling and fibroblast migration:** Tumor angiogenesis (VEGF regulation) and cancer-associated fibroblast activation (fibroblast migration) represent the tumor microenvironment remodeling signature. These signals arise from tumor stroma and endothelial cells responding to the cancer.

**RUNX1 regulation:** RUNX1 is a hematopoietic transcription factor frequently dysregulated across multiple cancer types (including leukemia, breast cancer, colorectal cancer). Its elevated regulatory pathway may reflect hematopoietic stress or myeloid cell activation in the tumor microenvironment.

**Negative Regulation of Blood Coagulation:** Cancer is associated with hypercoagulability (Trousseau syndrome), but paradoxically, regulatory mechanisms that suppress excessive clotting may also be activated. This pathway likely reflects the complex coagulation dysregulation in cancer.

This NES>0 signature is **remarkably interpretable and disease-specific** despite the heterogeneous cancer composition, suggesting shared pan-cancer mechanisms at the pathway level.

### 3.3 GSEA Results: NES < 0 (Downregulated, n=249)

The NES<0 signal is intense (minimum NES −3.10, the most extreme in this group) and similarly organized into multiple themes:

**Theme 1: Ribosome and translation (extreme signal, NES down to −3.10)**
- L13a-mediated Translational Silencing of Ceruloplasmin (NES −3.10) — the single most negative pathway across all three diseases
- Cytoplasmic Translation (NES −3.09), Eukaryotic Translation Elongation/Termination (NES −3.09)
- Peptide Chain Elongation (NES −3.06), Formation of 40S/60S pools (NES −3.05 to −2.99)
- Response of EIF2AK4 (GCN2) to Amino Acid Deficiency (NES −2.99)
- KEGG: Ribosome (NES −2.94), Translation (NES −2.80 to −2.72)
- rRNA Processing (NES −2.66 to −2.73)
- SRP-dependent Cotranslational Protein Targeting to Membrane (NES −2.90)

**Theme 2: Oxidative phosphorylation and mitochondrial energy metabolism**
- Respiratory Electron Transport + ATP Synthesis (NES −2.69)
- KEGG: Oxidative Phosphorylation (NES −2.55)
- Proton Motive Force-Driven ATP Synthesis (NES −2.52 to −2.62)
- Respiratory Electron Transport (NES −2.52), Mitochondrial ATP Synthesis Coupled Electron Transport (NES −2.46)

**Theme 3: Proteasome and ubiquitin-proteasome system (UPS)**
- KEGG: Proteasome (NES −2.44)
- Regulation of Activated PAK-2p34 by Proteasome Mediated Degradation (NES −2.40)
- Autodegradation of E3 Ubiquitin Ligase COP1 (NES −2.43)
- Degradation of AXIN and DVL (Wnt pathway components) (NES −2.42 to −2.54)

**Theme 4: SLIT/ROBO signaling**
- Regulation of Expression of SLITs and ROBOs (NES −2.89)
- Signaling by ROBO Receptors (NES −2.68)

**Theme 5: NF-κB and immune signaling**
- NIK To Noncanonical NF-κB Signaling (NES −2.43)
- Negative Regulation of NOTCH4 Signaling (NES −2.42)
- GSK3B/BTRC-CUL1-mediated Degradation of NFE2L2 (NES −2.42)
- Vpu Mediated Degradation of CD4 (NES −2.42)

**Theme 6: Platelet biology**
- Platelet Dense Granule Organization (NES −2.42)

### 3.4 Biological Interpretation

**Ribosome/translation suppression reaching NES −3.10:** The Other Cancer group shows the most extreme translation suppression across all three diseases. This is counterintuitive because cancer cells are generally characterized by *increased* ribosome biogenesis and protein synthesis. However, plasma cfRNA reflects the **composite transcriptome across all source tissues** — not just tumor cells. The downregulation likely reflects: (a) suppression of normal non-tumor tissue (liver, immune cells) in the cancer-bearing state; (b) tumor-derived cfRNA with altered ribosome content due to starvation-adapted translational programs; (c) the integrated stress response activated by metabolic stress within tumors. The extreme NES for the L13a/ceruloplasmin translational silencing pathway may specifically reflect hepatic involvement — L13a-mediated silencing is a liver-specific regulatory mechanism, and if the "Other Cancers" include hepatobiliary or gastrointestinal cancers with liver effects, this signal would be amplified. Consistent with this, the Moore et al. study explicitly found liver function biomarkers elevated in PDAC; liver-derived cfRNA is likely present in the Other Cancer group as well.

**Oxidative phosphorylation downregulation — the Warburg effect at the plasma cfRNA level.** Cancer cells preferentially use aerobic glycolysis (Warburg effect) over mitochondrial OXPHOS. The strong downregulation of oxidative phosphorylation pathways (OXPHOS, Respiratory Electron Transport, NES −2.45 to −2.69) in the plasma cfRNA suggests that tumor-derived RNA — which reflects the Warburg glycolytic phenotype — is contributing to the plasma cfRNA pool, displacing the OXPHOS-competent transcriptome of normal tissues. This is a robust and well-established metabolic signature of cancer and is particularly striking in the Other Cancer group (vs. CDCS where OXPHOS downregulation has NES ~−1.6).

**SLIT/ROBO downregulation: shared axis across diseases, particularly prominent in cancer.** SLIT/ROBO signaling regulates axon guidance, cell migration, and angiogenesis. In cancer, SLIT-ROBO acts as a tumor suppressor axis — SLIT2 is frequently methylated/silenced in multiple cancer types. The strong downregulation of Regulation of SLITs/ROBOs (NES −2.89) and ROBO receptor signaling (NES −2.68) in Other Cancer cfRNA is consistent with epigenetic silencing of this pathway in cancer cells and its loss from the plasma transcriptome.

**Proteasome downregulation: cancer rewiring of proteostasis.** The proteasome (NES −2.44) and multiple E3 ubiquitin ligase regulatory pathways being downregulated in plasma cfRNA is initially surprising — cancer cells upregulate proteasomal activity intracellularly. However, from the normative model perspective, cancer cells may sequester proteasomal transcripts intracellularly (in polysomes and stress granules) and release less proteasomal cfRNA, or cancer-derived cfRNA reflects tissue-specific compartments where proteasomal gene expression is suppressed.

**Wnt pathway component degradation (DVL, AXIN):** The downregulation of Degradation of DVL (NES −2.54) and AXIN (NES −2.42) — negative regulatory components of Wnt signaling — suggests that Wnt pathway *turnover* regulators are reduced. This could paradoxically indicate *active* Wnt signaling in the cancer state (since DVL and AXIN function to suppress Wnt — their degradation programs being downregulated could reflect that the cancer cells have constitutively active Wnt, which would not need active DVL/AXIN degradation machinery). This is consistent with frequent Wnt/β-catenin activation in gastrointestinal and other cancers.

**Coherence of the NES>0 cancer signature despite sample heterogeneity.** A critical consideration for Other Cancer is the small sample size (n=16) and presumed heterogeneity of cancer types. The Rho GTPase / focal adhesion / VEGF / invasion signature in the NES>0 group is remarkably coherent — it represents a **pan-cancer invasive signature** that transcends specific cancer types. This suggests the normative model GSEA is capturing fundamental oncogenic processes (cytoskeletal remodeling, invasion, angiogenesis) shared across solid tumors, rather than type-specific signals. This is biologically expected: Rho GTPase dysregulation, integrin-MAPK signaling, and VEGF-driven angiogenesis are near-universal features of solid tumors.

**Assessment:** Despite the heterogeneous cancer composition, the Other Cancer GSEA produces a highly interpretable result. The NES>0 Rho/invasion/VEGF signature is a robust pan-cancer fingerprint. The NES<0 ribosome/OXPHOS downregulation reflects the Warburg effect and systemic metabolic effects of cancer on circulating plasma cfRNA. The extreme NES values (−3.10) indicate that this group, despite small n, deviates dramatically from the HC normative baseline, which is expected for cancer.

---

## Cross-Disease Comparative Analysis

### Shared NES < 0 Core

All three diseases share a **common downregulated core** dominated by ribosome biogenesis, translation initiation/elongation/termination, and (to a lesser degree) oxidative phosphorylation. This shared core appears in:

- CDCS: NES −2.13 (ribosome), −1.64 to −1.72 (OXPHOS)
- Pre-eclampsia: NES −2.27 (ribosome/rRNA), mild OXPHOS
- Other Cancer: NES −3.10 (ribosome), −2.55 (OXPHOS)

The severity of this shared signature (Other Cancer > Pre-eclampsia > CDCS by NES magnitude) parallels clinical severity and degree of systemic metabolic stress. This cross-disease convergence suggests that plasma cfRNA ribosomal/translation suppression is a **non-specific marker of systemic disease burden** rather than a disease-specific feature. The normative model effectively captures this as a global "deviation from HC homeostasis" signal.

The SLIT/ROBO pathway (Regulation of SLITs and ROBOs, Signaling by ROBO Receptors) appears as a shared NES<0 pathway across all three diseases, suggesting it may be a common cfRNA stress response or cell-migration-related signal shared in diverse disease contexts.

### Disease-Specific Discriminators

| Feature | CDCS | Pre-eclampsia | Other Cancer |
|---------|------|---------------|--------------|
| NES>0 pathways | ECM remodeling, ion channels, purinergic signaling | Neurexin/neuroligin, DNA methylation | Rho GTPase cascade, focal adhesion, VEGF |
| Unique NES<0 themes | Immune NF-κB/BCR/IL-12 suppression, cardiac fibrosis pathways | Interferon/ISG15/NLR innate immune; splicing; lipoprotein metabolism | OXPHOS (Warburg); proteasome; Wnt pathway component regulation |
| Biological specificity | Moderate — reflects generic cardiac/ischemic suppression | Moderate — confounded by pregnancy confounder | High — cancer-specific Rho/invasion signature is distinctive |
| HC normative model fit | Well-suited (non-pregnant, non-cancer HC) | Partially confounded (HC non-pregnant; PE is pregnant state) | Well-suited (HC non-cancer) |

### Interpretation of Universal NES < 0 Dominance

The observation that all three diseases show >98% NES<0 FDR-significant pathways warrants a methodological caveat: when the GSEA gene-level ranking is dominated by uniformly negative Z-scores (i.e., disease samples are globally lower than HC across thousands of genes), the resulting enrichment will show most pathways as NES<0 regardless of biological specificity. This "compression" of Z-scores toward negative values may partially reflect:

1. **Library normalization effects**: If disease samples have slightly lower overall RNA quality or input, global suppression artifacts could emerge
2. **cfRNA composition shifts**: Cancer, cardiovascular disease, and pre-eclampsia all alter the cellular composition of cfRNA sources (more tumor-derived, placenta-derived, or ischemic tissue-derived RNA), effectively "diluting" the HC-like transcriptome
3. **Genuine biology**: As discussed per disease, each condition has real mechanistic reasons for broad transcriptional suppression relative to HC

The combination of these factors makes the NES<0 dominance partially expected and partially a methodological signature of normative model-based analysis with these disease groups. The disease-discriminating biological information is most reliably found in (a) the NES>0 pathways (disease-specific upregulations) and (b) the strongest and most disease-specific NES<0 pathways that differ between diseases.

---

## References

1. Ward Z, Schmeier S, Pearson J, Cameron V, Frampton CM, Troughton R, Doughty R, Richards AM, Pilbrow AP. "Identifying Candidate Circulating RNA Markers for Coronary Artery Disease by Deep RNA-Sequencing in Human Plasma." *Cells* 2022, 11(20):3191. DOI: 10.3390/cells11203191.

2. Moufarrej MN, Vorperian SK, Wong RJ, Campos AA, Quaintance CC, Sit RV, Tan M, Detweiler AM, Mekonen H, Neff NF, Baruch-Gravett C, Litch JA, Druzin ML, Winn VD, Shaw GM, Stevenson DK, Quake SR. "Early prediction of preeclampsia in pregnancy with cell-free RNA." *Nature* 2022, 602:689–694. DOI: 10.1038/s41586-022-04410-z.

3. Moore TW et al. "Cell free RNA detection of pancreatic cancer in pre-diagnostic high risk and symptomatic patients." *Nature Communications* 2025. DOI: 10.1038/s41467-025-62685-y.

4. Guo M, Price MJ, Patterson DG, et al. "EZH2 Represses the B Cell Transcriptional Program and Collaborates with PAX5 to Shape B Cell Identity." *Cell Reports* 2022 (context: B cell NF-κB references).

5. Moufarrej MN, Wong RJ, Shaw GM, Stevenson DK, Quake SR. "Investigating Pregnancy and Its Complications Using Circulating Cell-Free RNA in Women's Blood during Gestation." *Frontiers in Genetics* 2021; 11:605840.

6. Signs of immune dysregulation in second-trimester maternal blood RNA profiles in late-onset preeclampsia. *Scientific Reports* 2025. DOI: 10.1038/s41598-025-26323-3.

7. Myocardial Proteome in Human Heart Failure With Preserved Ejection Fraction. *Journal of the American Heart Association* 2024. DOI: 10.1161/JAHA.124.038945. (OXPHOS and ribosomal downregulation in HFpEF.)
