# GSEA Biological Evaluation — Group 2: Liver/GI Cancers

**Disease group:** Liver Cancer | Liver Cirrhosis | Colorectal Cancer | Esophagus Cancer | Stomach Cancer  
**Method:** GSEA-prerank on mean Z-score per phenotype (GAMLSS normative model, 693 HC reference)  
**NES interpretation:** NES > 0 = pathway enriched in high-Z genes (anomalously upregulated vs. HC); NES < 0 = enriched in low-Z genes (anomalously downregulated vs. HC)  
**Significance threshold:** FDR q-val < 0.05  
**Gene set libraries:** Reactome 2022, GO Biological Process 2023, KEGG 2021 Human  

**Source papers:**
- Chen et al. (eLife 2022) — DOI: 10.7554/elife.75181 — cfRNA-seq of ~295 plasma samples across 5 GI/liver cancer types vs. healthy donors
- Roskams-Hieter et al. (npj Precision Oncology 2022) — DOI: 10.1038/s41698-022-00270-y — cfRNA profiling of HCC (n=8 discovery, n=20 validation), liver cirrhosis (n=4 premalignant), and non-cancer controls
- Block et al. (Frontiers in Oncology 2022) — DOI: 10.3389/fonc.2022.963641 — cfRNA variant calling in 14 HCC / 8 cirrhosis / 6 healthy controls
- Peddu et al. (bioRxiv 2025) — DOI: 10.1101/2025.07.02.662774 — LOCATE-seq nanopore cfRNA for esophageal adenocarcinoma and Barrett's esophagus

---

## Summary Statistics

| Disease | Total terms tested | FDR<0.05 significant | NES>0 | NES<0 |
|---|---|---|---|---|
| Liver Cancer | 343 | 343 (100%) | 319 | 24 |
| Liver Cirrhosis | 41 | 41 (100%) | 39 | 2 |
| Colorectal Cancer | 383 | 383 (100%) | 21 | 362 |
| Esophagus Cancer | 681 | 681 (100%) | 142 | 539 |
| Stomach Cancer | 520 | 520 (100%) | 17 | 503 |

> Note: 100% significance rates across all diseases indicate a globally shifted cfRNA transcriptome relative to the HC normative model. The directional split (NES>0 vs. NES<0) is the key biological discriminant.

---

## 1. Liver Cancer

**Source datasets:** Chen et al. (eLife 2022) — liver cancer as one of five GI cancer types; Roskams-Hieter et al. (2022) — 8 HCC discovery + 20 HCC validation; Block et al. (2022) — 14 HCC patients (11 with concurrent cirrhosis)  
**cfRNA context:** Roskams-Hieter identified liver-specific secreted proteins (APOE, C3, ceruloplasmin/CP, fibrinogen chains FGA/FGB/FGG) as HCC discriminators in plasma cfRNA. Block et al. found metabolism-related gene alterations (FASN, SREBF1, PKLR, PCK1) predominantly in HCC cfRNA. Chen et al. found ECM-receptor interactions and neutrophil extracellular trap pathways upregulated, and ribosome/immune pathways downregulated in cancer cfRNA broadly.

### 1.1 NES>0 Pathways (Anomalously Upregulated in HCC, n=319 significant)

**Top 30 by NES:**

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | 2.192 | <0.001 | KEGG: Complement and coagulation cascades |
| 2 | 2.118 | 0.005 | GOBP: Negative Regulation Of Blood Coagulation |
| 3 | 2.115 | 0.003 | GOBP: Reverse Cholesterol Transport |
| 4 | 2.101 | 0.003 | Reactome: Plasma Lipoprotein Remodeling |
| 5 | 2.089 | 0.002 | GOBP: High-Density Lipoprotein Particle Remodeling |
| 6 | 2.083 | 0.002 | Reactome: Post-translational Protein Phosphorylation |
| 7 | 2.075 | 0.002 | Reactome: Regulation Of IGF Transport And Uptake By IGFBPs |
| 8 | 2.073 | 0.003 | GOBP: Regulation Of Endothelial Cell Migration |
| 9 | 2.067 | 0.003 | GOBP: Phospholipid Efflux |
| 10 | 2.060 | 0.003 | Reactome: Regulation Of Complement Cascade |
| 11 | 2.022 | 0.007 | Reactome: Metabolism Of Porphyrins |
| 12 | 2.020 | 0.007 | GOBP: Cholesterol Transport |
| 13 | 2.018 | 0.006 | KEGG: Cholesterol metabolism |
| 14 | 2.017 | 0.006 | Reactome: Platelet Degranulation |
| 15 | 2.009 | 0.006 | Reactome: Binding And Uptake Of Ligands By Scavenger Receptors |
| 16 | 2.009 | 0.006 | Reactome: Plasma Lipoprotein Assembly, Remodeling, And Clearance |
| 17 | 2.005 | 0.005 | Reactome: Scavenging Of Heme From Plasma |
| 18 | 1.989 | 0.005 | Reactome: RHOV GTPase Cycle |
| 19 | 1.987 | 0.005 | Reactome: Retinoid Metabolism And Transport |
| 20 | 1.985 | 0.005 | Reactome: Aspirin ADME |
| 21 | 1.984 | 0.005 | Reactome: Metabolism Of Fat-Soluble Vitamins |
| 22 | 1.983 | 0.005 | Reactome: Response To Elevated Platelet Cytosolic Ca2+ |
| 23 | 1.973 | 0.005 | Reactome: Signaling By High-Kinase Activity BRAF Mutants |
| 24 | 1.972 | 0.005 | Reactome: p130Cas Linkage To MAPK Signaling For Integrins |
| 25 | 1.971 | 0.005 | GOBP: Artery Development |
| 26 | 1.970 | 0.005 | GOBP: Aortic Valve Development |
| 27 | 1.969 | 0.005 | Reactome: Plasma Lipoprotein Assembly |
| 28 | 1.968 | 0.005 | Reactome: HDL Remodeling |
| 29 | 1.958 | 0.006 | GOBP: Positive Regulation Of Substrate Adhesion-Dependent Cell Spreading |
| 30 | 1.954 | 0.006 | KEGG: Adherens junction |

**Key lead genes identified:**
- Complement/coagulation: FGB, FGA, PLG, KNG1, FGG, CFB, A2M, C3, C9, CFHR1, C1S, CLU, CFH, SERPINC1, SERPING1, F13A1, C1QB, F2, F8, F9 (classical and alternative complement + fibrinogen cascade)
- Lipoproteins: ALB, APOB, APOE, APOC1, APOC2, APOC3, APOA1, APOA4, APOA2, ANGPTL3, LPL, LIPG, LDLR
- Heme/porphyrin: ALB, FABP1, SLCO1B1, BLVRB, HMBS, ALAS2, FECH, CPOX, SLCO1B3, ABCC2 (liver transporters prominent)
- IGF pathway: ALB, APOB, FGA, FGG, TF, CP, IGFBP4, IGFBP5, IGFBP7, IGF1, IGF2
- Heme scavenging: ALB, HPX (hemopexin), AMBP, APOA1 (classical plasma heme-binding proteins)
- Platelet degranulation: FGB, ALB, APOH, FGA, PLG, KNG1, FGG, TF, ORM1, SERPINA3, ITIH3, SELENOP, IGF1, IGF2, TTR

**Biological interpretation — NES>0:**

The dominant upregulated signal in HCC cfRNA is a striking enrichment of **hepatocyte-secreted plasma proteins** spanning complement, coagulation, and lipoprotein metabolism. This is biologically coherent on two levels:

1. **Acute-phase response and liver stress:** HCC is typically superimposed on chronic liver disease with systemic inflammation. The complement cascade lead genes (FGA/FGB/FGG fibrinogen, C3, C1QB/C1S, SERPING1, A2M) are all canonical acute-phase proteins synthesized by hepatocytes and released into blood. Their cfRNA appearance in plasma may reflect increased shedding of hepatocyte-origin transcripts during liver disease. Published literature confirms complement activation is prominent in HCC: CFHR4 is a reported poor-prognostic complement biomarker, and C1QB/C4BPA are elevated in HCC-derived extracellular vesicles (PMC9724420).

2. **Liver-specific gene enrichment as a normative model artifact/signal:** The HC normative model was trained on 693 healthy controls. The HCC samples upregulate hepatocyte-secreted transcripts relative to HC baseline, which is **counterintuitive** at first — one would expect HCC to *reduce* normal hepatocyte function. However, multiple cfRNA studies (Roskams-Hieter 2022) found that HCC plasma cfRNA is enriched for liver-specific secreted proteins (APOE, C3, CP, FGA/FGB/FGG) as discriminating biomarkers. This likely reflects: (a) increased hepatocyte turnover and cell death releasing transcripts into plasma, (b) inflammatory amplification of acute-phase gene expression, and (c) tumor-derived shedding of hepatocyte-lineage transcripts. The Z-score model captures these as higher-than-expected signals relative to HC.

3. **Lipoprotein metabolism (ALB, APOB, APOE, APOCs, APOA1):** Hepatocytes are the primary producers of apolipoproteins. Elevated cfRNA for these in HCC aligns with chronic liver stress and the observation that HCC disrupts normal lipid metabolism while still releasing hepatocyte-derived RNA. Reverse cholesterol transport and HDL remodeling are liver-centric processes; their enrichment here likely tracks the same hepatocyte-origin transcript pool.

4. **Porphyrin/heme metabolism (ALAS2, HMBS, FECH, CPOX; liver transporters SLCO1B1/SLCO1B3):** These are hallmark liver-specific genes. Their upregulation in the Z-score ranking reflects the hepatic origin of the cfRNA signal in HCC. ALB leading multiple terms is expected (albumin is the most abundant secreted protein in plasma, produced exclusively by hepatocytes).

5. **IGF transport (IGFBPs, ALB, APOB, TF, CP):** IGFBP signaling is altered in HCC. IGF2 is a known HCC driver gene. IGF pathway enrichment in cfRNA is consistent with reports of elevated IGF2/IGFBP expression in HCC tumor tissue and plasma.

6. **Platelet degranulation:** Many platelet granule proteins (fibrinogen chains, APOH, VTN, SERPINC1, F13A1, TTR, SPARC) are shared with hepatocyte-secreted products; this term likely co-clusters with the coagulation/complement signal rather than reflecting true platelet activation per se.

7. **Retinoid metabolism, fat-soluble vitamins, aspirin ADME:** All highly liver-specific processes (hepatocytes are the primary site of vitamin A/D/E/K storage and metabolism). Their enrichment further confirms the liver-origin nature of the upregulated cfRNA pool in HCC.

### 1.2 NES<0 Pathways (Anomalously Downregulated in HCC, n=24 significant)

**Top 24 (all significant):**

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | -2.919 | <0.001 | Reactome: Respiratory Electron Transport, ATP Synthesis (full OXPHOS) |
| 2 | -2.833 | <0.001 | KEGG: Oxidative phosphorylation |
| 3 | -2.822 | <0.001 | Reactome: Respiratory Electron Transport |
| 4 | -2.733 | <0.001 | GOBP: Proton Motive Force-Driven ATP Synthesis |
| 5 | -2.675 | <0.001 | GOBP: Aerobic Electron Transport Chain |
| 6 | -2.610 | <0.001 | GOBP: Mitochondrial ATP Synthesis Coupled Electron Transport |
| 7 | -2.581 | <0.001 | GOBP: Mitochondrial Respiratory Chain Complex Assembly |
| 8 | -2.573 | <0.001 | Reactome: Complex I Biogenesis |
| 9 | -2.535 | <0.001 | GOBP: Cellular Respiration |
| 10 | -2.359 | <0.001 | GOBP: Mitochondrial Electron Transport, Cytochrome C To Oxygen |
| 11 | -2.322 | 0.002 | GOBP: Proton Motive Force-Driven Mitochondrial ATP Synthesis |
| 12 | -2.316 | 0.001 | GOBP: Mitochondrial Respiratory Chain Complex I Assembly |
| 13 | -2.316 | 0.001 | GOBP: NADH Dehydrogenase Complex Assembly |
| 14 | -2.305 | 0.003 | GOBP: Oxidative Phosphorylation |
| 15 | -2.111 | 0.015 | GOBP: MHC Class II Protein Complex Assembly |
| 16 | -2.111 | 0.015 | GOBP: Peptide Antigen Assembly With MHC Class II Protein Complex |
| 17 | -2.105 | 0.015 | Reactome: Citric Acid (TCA) Cycle And Respiratory Electron Transport |
| 18 | -2.048 | 0.025 | GOBP: Mitochondrial Electron Transport, NADH To Ubiquinone |
| 19 | -2.037 | 0.025 | GOBP: Late Endosome To Vacuole Transport Via MVB Sorting Pathway |
| 20 | -2.000 | 0.035 | GOBP: Regulation Of Endosome Size |
| 21 | -1.994 | 0.037 | GOBP: Mitochondrial Cytochrome C Oxidase Assembly |
| 22 | -1.980 | 0.042 | GOBP: Peroxisome Organization |
| 23 | -1.961 | 0.047 | KEGG: Asthma |
| 24 | -1.957 | 0.048 | GOBP: ATP Biosynthetic Process |

**Key lead genes:** NDUFA1/2/5/6/11/13, NDUFB1/2/3/4/9/11, NDUFC1, NDUFS2/5/7, NDUFV1, UQCR11, UQCRH, UQCRQ, ATP5MG, ATP5F1E, ATP5PF, ATP5ME, COX4I1, COX5B, COX6A1, COX6C, COX7A2L, COX7C, SURF1, CYCS, ETFA, ETFB; MHC-II: HLA-DRB1, HLA-DMB, HLA-DRB5, HLA-DPA1, HLA-DQA1, HLA-DMA, HLA-DRA

**Biological interpretation — NES<0:**

1. **Oxidative phosphorylation (OXPHOS) / mitochondrial respiratory chain — dominant downregulated signal (NES up to -2.92):** This is the most strongly downregulated pathway cluster in HCC, spanning 14 of 24 significant negative terms. The lead genes are almost exclusively NADH dehydrogenase complex I subunits (NDUF family), ubiquinol-cytochrome c reductase subunits (UQCR), cytochrome c oxidase subunits (COX), and ATP synthase subunits (ATP5). This is biologically validated: HCC is well-known to reprogram metabolism toward glycolysis (Warburg effect), with OXPHOS subunit downregulation documented at both tumor tissue and plasma cfRNA levels (PMC9315135). Mitochondrial dysfunction with Complex I defects has been confirmed as a major pathologic feature of HCC. The normative model correctly captures this as a significant reduction relative to HC, suggesting that OXPHOS cfRNA transcripts (likely of immune/erythroid cell origin that are reduced, or reflecting hepatocellular metabolic reprogramming) are depleted in HCC plasma.

2. **MHC Class II complex assembly (HLA-DRB1, HLA-DRB5, HLA-DPA1, HLA-DQA1, HLA-DMA):** Downregulation of MHC-II gene transcripts in HCC cfRNA is consistent with known immunosuppression in HCC. HCC creates an immunosuppressive tumor microenvironment, downregulating antigen presentation to evade immune surveillance. The reduction of MHC-II cfRNA likely reflects depleted or dysregulated antigen-presenting cells (dendritic cells, B cells) in the circulation of HCC patients.

3. **Peroxisome organization, MVB/endosome trafficking:** Less characterized in cfRNA, but could reflect reduced normal hepatocyte lipid handling (peroxisomes are essential for fatty acid beta-oxidation in healthy hepatocytes, a process disrupted in HCC).

4. **Asthma pathway (NES=-1.96, FDR=0.047):** Likely reflects the downregulation of several immune effector genes (IL-4/IL-5/IL-13 signaling components, IgE-related) that share membership in this KEGG pathway. Not biologically specific to asthma in this context.

### 1.3 Special Assessment: Liver-Specific Proteins

- **Albumin (ALB):** Appears as a lead gene in NES>0 terms: Plasma Lipoprotein Remodeling, Scavenging Of Heme From Plasma, Metabolism Of Porphyrins, IGF Transport, Platelet Degranulation. ALB cfRNA is thus **anomalously elevated** in HCC relative to HC normative model. This contrasts with the serum protein result (albumin decreases in advanced liver disease), suggesting the cfRNA model captures a different signal — possibly increased hepatocyte death releasing ALB transcripts, or tumor-derived ALB expression in well-differentiated HCC.
- **Fibrinogen chains (FGA, FGB, FGG):** Lead genes for complement/coagulation (NES=2.192), platelet degranulation (NES=2.017), and IGF transport terms. Upregulated in HCC cfRNA — consistent with published reports of elevated fibrinogen alpha/beta chains in HCC plasma proteomics (Nature Scientific Reports 2025: fibrinogen alpha chain significantly elevated in HCC plasma).
- **Complement components (C3, C9, CFB, CFHR1, CFH):** Upregulated — consistent with complement activation in HCC inflammation and known HCC prognostic relevance of complement factor H-related proteins (CFHR4 PMC9632743).
- **Apolipoproteins (APOE, APOB, APOA1/2/4, APOC1/2/3):** All upregulated — Roskams-Hieter 2022 explicitly identified APOE and fibrinogen chains as top HCC biomarker genes in cfRNA, validating this finding independently.

### 1.4 Signal Quality Assessment

**SIGNAL QUALITY: HIGH**

Rationale: The NES>0 cluster (complement, coagulation, lipoprotein/apolipoprotein, heme metabolism, IGF) is highly coherent and fully consistent with published HCC cfRNA biomarker literature (Roskams-Hieter 2022 identifies APOE, C3, CP, fibrinogens as HCC discriminators). The NES<0 OXPHOS signal (NES up to -2.92) is the strongest negative signal in the entire dataset and is validated by multiple HCC transcriptomic studies documenting mitochondrial OXPHOS defects. The Z-score model produces a biologically interpretable liver-cancer-specific profile with minimal noise. The 319:24 NES>0:NES<0 imbalance reflects that HCC cfRNA is globally elevated (positive acute-phase/secreted protein response) with a narrow but very strong downregulation of mitochondrial energy metabolism.

---

## 2. Liver Cirrhosis

**Source datasets:** Roskams-Hieter et al. (2022) — 4 liver cirrhosis samples (premalignant cohort in HCC study); Block et al. (2022) — 8 cirrhosis patients without HCC  
**Sample size caveat:** Only n=4 cirrhosis samples remained after OOD filtering in this normative modeling pipeline. All statistical conclusions should be treated with caution given this extremely small sample size. The 41 tested gene sets with 100% significance likely reflects low power to differentiate signal from noise with n=4.

### 2.1 NES>0 Pathways (n=39 significant)

**Top 30 by NES (all 41 terms listed in full):**

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | 2.363 | <0.001 | Reactome: HDACs Deacetylate Histones |
| 2 | 2.203 | 0.005 | GOBP: Regulation Of Endothelial Cell Migration |
| 3 | 2.165 | 0.004 | KEGG: Systemic lupus erythematosus |
| 4 | 2.157 | 0.003 | Reactome: Amyloid Fiber Formation |
| 5 | 2.100 | 0.004 | GOBP: Positive Regulation Of Vasoconstriction |
| 6 | 2.041 | 0.018 | KEGG: Complement and coagulation cascades |
| 7 | 2.007 | 0.032 | GOBP: Negative Regulation Of Epithelial Cell Migration |
| 8 | 1.982 | 0.039 | GOBP: Negative Regulation Of Blood Coagulation |
| 9 | 1.982 | 0.035 | Reactome: Metabolic Disorders Of Biological Oxidation Enzymes |
| 10 | 1.980 | 0.033 | GOBP: Positive Regulation Of Ubiquitin-Protein Transferase Activity |
| 11 | 1.978 | 0.033 | Reactome: Basigin Interactions |
| 12 | 1.976 | 0.030 | Reactome: Pre-NOTCH Transcription And Translation |
| 13 | 1.968 | 0.029 | Reactome: RNA Polymerase I Promoter Opening |
| 14 | 1.966 | 0.029 | Reactome: Metalloprotease DUBs |
| 15 | 1.960 | 0.029 | GOBP: Regulation Of Vasoconstriction |
| 16 | 1.954 | 0.032 | KEGG: Neutrophil extracellular trap formation (NETs) |
| 17 | 1.952 | 0.031 | Reactome: Oxidative Stress Induced Senescence |
| 18 | 1.950 | 0.030 | Reactome: Activated PKN1 Stimulates Transcription Of AR-Regulated KLK2/KLK3 |
| 19 | 1.947 | 0.031 | Reactome: Regulation Of TLR By Endogenous Ligand |
| 20 | 1.945 | 0.030 | GOBP: Regulation Of Membrane Protein Ectodomain Proteolysis |
| 21 | 1.941 | 0.029 | Reactome: Condensation Of Prophase Chromosomes |
| 22 | 1.939 | 0.028 | Reactome: Binding And Uptake Of Ligands By Scavenger Receptors |
| 23 | 1.935 | 0.030 | GOBP: Phosphatidylcholine Metabolic Process |
| 24 | 1.934 | 0.030 | GOBP: Negative Regulation Of Extrinsic Apoptotic Signaling Via Death Domain Receptors |
| 25 | 1.932 | 0.030 | Reactome: Defective Pyroptosis |
| 26 | 1.929 | 0.030 | Reactome: ERCC6/EHMT2 Positively Regulate rRNA Expression |
| 27 | 1.928 | 0.030 | GOBP: Positive Regulation Of Signaling Receptor Activity |
| 28 | 1.927 | 0.029 | GOBP: HDL Particle Remodeling |
| 29 | 1.924 | 0.029 | GOBP: Phospholipid Catabolic Process |
| 30 | 1.918 | 0.031 | Reactome: Scavenging Of Heme From Plasma |
| 31 | 1.916 | 0.030 | Reactome: PRC2 Methylates Histones And DNA |
| 32 | 1.915 | 0.030 | GOBP: Positive Regulation Of Long-Term Synaptic Potentiation |
| 33 | 1.915 | 0.030 | Reactome: Packaging Of Telomere Ends |
| 34 | 1.913 | 0.030 | GOBP: Sterol Metabolic Process |
| 35 | 1.912 | 0.030 | Reactome: DNA Methylation |
| 36 | 1.905 | 0.030 | Reactome: SIRT1 Negatively Regulates rRNA Expression |
| 37 | 1.894 | 0.030 | Reactome: Pre-NOTCH Expression And Processing |
| 38 | 1.889 | 0.030 | GOBP: RHO GTPases Activate PKNs |
| 39 | 1.881 | 0.030 | GOBP: Negative Regulation Of Fat Cell Differentiation |

### 2.2 NES<0 Pathways (n=2 significant)

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | -2.280 | 0.024 | GOBP: Proton Motive Force-Driven Mitochondrial ATP Synthesis |
| 2 | -2.203 | 0.040 | GOBP: Oxidative Phosphorylation |

### 2.3 Biological Interpretation

**NES>0 — Histone and chromatin signal (top hit: HDAC Deacetylate Histones, NES=2.363):**

The top enriched pathway is dominated by histone protein transcripts (H2AC4, H2AC7, H2AC11/12/14/19/20/21/25, H2BC11/12/13/14/17/18/21/26, H3C7/13, H2AX). This is a classical cfRNA background signal: linker and core histones (H2A, H2B, H3, H4) are among the most highly expressed genes in proliferating cells, and their cfRNA is abundant in plasma from circulating blood cells (neutrophils, platelets, erythroid precursors). In cirrhosis, there are two biological amplifiers:

1. **cfRNA background enrichment:** Cirrhosis is associated with systemic inflammation and neutrophil activation. Neutrophil extracellular traps (NETs) are explicitly enriched (rank 16, NES=1.954; lead genes include H2AC/H2BC histones + FGG, FGB, FCGR2A, ITGB3, PRKCA, C3). NETs release massive quantities of nucleosome-associated histone proteins and DNA into the bloodstream, explaining the histone cfRNA signal.

2. **Chromatin remodeling in liver fibrosis:** BRG1 (SMARCA4) chromatin remodeling drives hepatic stellate cell (HSC) activation and fibrosis. HDAC and PRC2 (rank 31, NES=1.916) activity is elevated in cirrhotic liver. It is plausible that fibrosis-related chromatin remodeling transcripts are shed into plasma, but the dominant contribution is likely the NET/neutrophil background.

The **Systemic lupus erythematosus** pathway (rank 3, NES=2.165) shares its lead genes almost entirely with the histone/complement signature: H2AC/H2BC histones + C3, C4B, C1QA, C9, FCGR2A. This term captures the nucleosome/complement-mediated immune activation seen in cirrhosis, where chronic inflammatory signaling activates complement through pattern recognition of released chromatin.

**NES>0 — Endothelial cell migration and vascular remodeling:**

Rank 2 (GOBP: Regulation Of Endothelial Cell Migration, NES=2.203) includes genes: GFUS, APOH, NR2F2, SASH1, PROX1, CALR, PIK3CB, SPARC, ITGB3, GPI, FGF1, EPHA2, PRKCA, EGF, APOE, RHOB, DCN, WNT5A, FLT4 (VEGFR-3), KDR (VEGFR-2), DLL4, FGF2, EDN1. This is biologically compelling for cirrhosis: portal hypertension drives extensive angiogenesis and sinusoidal endothelial dysfunction. LSEC (liver sinusoidal endothelial cell) dysfunction is an early event in fibrosis-to-cirrhosis progression. The cfRNA enrichment of endothelial remodeling genes (KDR/VEGFR-2, FLT4/VEGFR-3, DLL4/Notch ligand, EDN1/endothelin) likely reflects the vascular remodeling signal from activated LSECs, consistent with published LSEC transcriptome studies showing these as stage-specific cirrhosis markers.

**NES>0 — Amyloid fiber formation, Oxidative stress-induced senescence, Defective Pyroptosis:**

- **Amyloid fiber formation (NES=2.157):** Lead genes include APOA4, histones, UBB, APOE, FGA, APP (amyloid precursor protein), SNCA (alpha-synuclein), IAPP — again a mixed signal of liver-secreted apolipoproteins and histones. Cirrhosis is associated with protein aggregation and ER stress in hepatocytes. APP and SNCA elevation may reflect hepatocyte stress.
- **Oxidative Stress-Induced Senescence (NES=1.952):** Cirrhosis is a classic context for cellular senescence of hepatocytes and activated stellate cells. This biological signal is expected.
- **Defective Pyroptosis (NES=1.932):** Inflammasome/pyroptosis pathways are activated in cirrhosis; this enrichment aligns with the inflammatory milieu.

**NES>0 — Complement and coagulation (rank 6, NES=2.041):**

Shared with HCC (complement+coagulation upregulated). In cirrhosis, coagulopathy is a key clinical complication; upregulation of these cfRNA signals is consistent with hepatocyte stress and activation of the coagulation cascade.

**NES>0 — Epigenetic regulators (DNA Methylation, SIRT1 rRNA, PRC2, ERCC6/EHMT2):**

Multiple epigenetic regulation pathways upregulated. These are consistent with the known role of epigenetic reprogramming in liver fibrosis: DNMT3A/3B-driven methylation changes, SIRT1 deacetylase activity, and PRC2-mediated H3K27me3 in fibrogenic gene regulation. These signals in cirrhosis cfRNA suggest shedding of epigenetically active transcripts from stressed hepatocytes or activated stellate cells.

**NES<0 — OXPHOS (NES=-2.280, -2.203):**

Same as HCC — oxidative phosphorylation genes are downregulated in cirrhosis cfRNA relative to HC. Mitochondrial dysfunction is a recognized early event in liver fibrosis. The cirrhosis OXPHOS signal (NES=-2.28) is weaker than HCC (NES=-2.92), suggesting progressive mitochondrial impairment from cirrhosis to HCC.

**Note on hdac-histone artifact:** The HDACs Deacetylate Histones pathway achieving the top NES>0 position in cirrhosis (vs. rank well outside top 30 in HCC) is noteworthy. Liver Cancer samples may have a more diverse upregulated signal (complement/lipoprotein dominates), while cirrhosis shows a relatively more histone-enriched plasma cfRNA background, possibly due to active NETosis in portal hypertension.

### 2.4 Signal Quality Assessment

**SIGNAL QUALITY: LOW-TO-MODERATE (caveat: n=4)**

Rationale: With only n=4 samples post-OOD filtering, statistical conclusions are unreliable. The histone/NET signal is a known cfRNA background present across many conditions, not cirrhosis-specific. The endothelial migration and complement signals are biologically coherent but require validation with larger n. The OXPHOS downregulation mirrors HCC and has some biological plausibility for cirrhosis. Overall, the signal is biologically plausible in direction but severely underpowered.

---

## 3. Colorectal Cancer

**Source dataset:** Chen et al. (eLife 2022) — colorectal cancer as one of five GI cancer types in the ~295-sample cfRNA-seq study  
**Published cfRNA findings:** Downregulated cfRNAs in colorectal cancer patients were enriched in ribosome biogenesis pathways; downregulation of translation-related pathways in tumor-educated platelets and plasma of CRC patients has been independently reported. Published plos ONE 2024 cfRNA study for CRC confirms translation/ribosome pathway suppression in plasma.

### 3.1 NES>0 Pathways (Anomalously Upregulated, n=21 significant)

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | 2.542 | <0.001 | Reactome: Signaling By Hippo |
| 2 | 2.542 | <0.001 | Reactome: YAP1- And WWTR1(TAZ)-stimulated Gene Expression |
| 3 | 2.488 | <0.001 | GOBP: Positive Regulation Of Integrin-Mediated Signaling Pathway |
| 4 | 2.410 | <0.001 | GOBP: Positive Regulation Of Insulin Receptor Signaling Pathway |
| 5 | 2.408 | <0.001 | GOBP: Regulation Of Platelet Aggregation |
| 6 | 2.406 | <0.001 | GOBP: Positive Regulation Of Cellular Response To Insulin Stimulus |
| 7 | 2.291 | 0.009 | GOBP: Maintenance Of Blood-Brain Barrier |
| 8 | 2.251 | 0.008 | Reactome: Regulation Of KIT Signaling |
| 9 | 2.232 | 0.007 | GOBP: Protein Localization To Cell-Cell Junction |
| 10 | 2.228 | 0.006 | Reactome: Constitutive Signaling By Overexpressed ERBB2 |
| 11 | 2.213 | 0.005 | Reactome: Nephrin Family Interactions |
| 12 | 2.160 | 0.008 | Reactome: Signaling By KIT In Disease |
| 13 | 2.137 | 0.014 | GOBP: Detection Of Calcium Ion |
| 14 | 2.129 | 0.017 | GOBP: Protein Localization To Cell Junction |
| 15 | 2.126 | 0.016 | Reactome: Signaling By ERBB2 In Cancer |
| 16 | 2.123 | 0.017 | GOBP: Regulation Of Receptor Signaling Pathway Via STAT |
| 17 | 2.084 | 0.034 | Reactome: Constitutive Signaling By EGFRvIII |
| 18 | 2.052 | 0.048 | GOBP: Positive Regulation Of Glucose Import |
| 19 | 2.041 | 0.046 | GOBP: Cellular Response To Growth Hormone Stimulus |
| 20 | 2.010 | 0.050 | GOBP: Embryonic Appendage Morphogenesis |
| 21 | 2.000 | 0.050 | Reactome: Apoptotic Cleavage Of Cell Adhesion Proteins |

**Key lead genes:**
- Hippo/YAP: TJP1, STK24, YAP1, WWTR1, TJP2, AMOT, MOB1A, LATS2, AMOTL2; HIPK2, KAT2B, TEAD4, HIPK1
- Integrin signaling: LIMS1, NID1, DMTN, LAMB2, LAMC1, FLNA, LOXL3, LIMS2, LAMB1
- KIT signaling: SH2B2, GRB2, KITLG, LYN, SRC, PRKCA, STAT5B, PIK3R3, STAT1, FYN, JAK2
- ERBB2: HSP90AA1, PTPN12, GRB2, CDC37, SOS1, SHC1, PIK3R1, EGF
- STAT: EGF, TGFB1, PPARG, JAK3, JAK2, CRLF2

### 3.2 NES<0 Pathways (Anomalously Downregulated, n=362 significant)

**Top 30 by NES:**

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | -2.432 | <0.001 | Reactome: Eukaryotic Translation Elongation |
| 2 | -2.427 | <0.001 | Reactome: Peptide Chain Elongation |
| 3 | -2.418 | <0.001 | Reactome: Viral mRNA Translation |
| 4 | -2.403 | <0.001 | Reactome: NMD Independent Of EJC |
| 5 | -2.397 | <0.001 | Reactome: Eukaryotic Translation Termination |
| 6 | -2.386 | <0.001 | Reactome: Selenocysteine Synthesis |
| 7 | -2.371 | <0.001 | Reactome: SRP-dependent Cotranslational Protein Targeting To Membrane |
| 8 | -2.350 | <0.001 | GOBP: Cytoplasmic Translation |
| 9 | -2.350 | <0.001 | Reactome: NMD Enhanced By EJC |
| 10 | -2.350 | <0.001 | Reactome: Formation Of A Pool Of Free 40S Subunits |
| 11 | -2.339 | <0.001 | KEGG: Ribosome |
| 12 | -2.339 | <0.001 | Reactome: Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency |
| 13 | -2.327 | <0.001 | Reactome: L13a-mediated Translational Silencing Of Ceruloplasmin Expression |
| 14 | -2.319 | <0.001 | Reactome: GTP Hydrolysis And Joining Of 60S Ribosomal Subunit |
| 15 | -2.311 | <0.001 | Reactome: Cap-dependent Translation Initiation |
| 16 | -2.252 | <0.001 | GOBP: Peptide Biosynthetic Process |
| 17 | -2.246 | <0.001 | GOBP: Ribosomal Small Subunit Biogenesis |
| 18 | -2.246 | <0.001 | Reactome: Selenoamino Acid Metabolism |
| 19 | -2.226 | <0.001 | Reactome: Regulation Of Expression Of SLITs And ROBOs |
| 20 | -2.220 | <0.001 | Reactome: Influenza Viral RNA Transcription And Replication |
| 21 | -2.205 | <0.001 | Reactome: Cellular Response To Starvation |
| 22 | -2.194 | <0.001 | GOBP: Translation |
| 23 | -2.187 | <0.001 | GOBP: Macromolecule Biosynthetic Process |
| 24 | -2.180 | <0.001 | Reactome: Translation |
| 25 | -2.179 | <0.001 | Reactome: SARS-CoV-2 Modulates Host Translation Machinery |
| 26 | -2.179 | <0.001 | Reactome: Major Pathway Of rRNA Processing In Nucleolus And Cytosol |
| 27 | -2.174 | <0.001 | Reactome: rRNA Processing In Nucleus And Cytosol |
| 28 | -2.171 | <0.001 | Reactome: rRNA Processing |
| 29 | -2.159 | <0.001 | GOBP: Ribonucleoprotein Complex Biogenesis |
| 30 | -2.159 | <0.001 | Reactome: mRNA Activation Upon Binding Of Cap-Binding Complex And eIFs |

### 3.3 Biological Interpretation

**NES>0 — Hippo/YAP signaling (top signal, NES=2.542):**

YAP1 (Yes-associated protein) and WWTR1/TAZ are the terminal effectors of the Hippo tumor suppressor pathway. When Hippo kinases (LATS1/2, MOB1A) are inactivated, YAP/TAZ translocate to the nucleus and drive TEAD-mediated transcription of pro-proliferative and anti-apoptotic genes. YAP1 is upregulated in CRC tumor tissue and drives invasion; LATS2 (a Hippo kinase) appears in the lead genes alongside YAP1 and WWTR1, indicating pathway deregulation. Importantly, YAP/TAZ-stimulated gene expression (rank 2) includes HIPK2, KAT2B, TEAD4 — direct transcriptional targets. The cfRNA upregulation of Hippo pathway components in CRC plasma represents a **validated cancer hallmark signal**: Hippo-YAP dysregulation is a recognized oncogenic driver in CRC (PMC12221920), and circulating transcripts of YAP1/WWTR1 likely originate from tumor-educated platelets or shed tumor cells.

**NES>0 — Integrin signaling (NES=2.488):**

Lead genes LIMS1/2 (PINCH proteins), NID1 (nidogen), LAMB1/LAMB2 (laminin beta), LAMC1 (laminin gamma), FLNA (filamin A), LOXL3 (lysyl oxidase-like). Integrin-ECM interactions are fundamental to CRC invasion and metastasis. LAMC1 and LAMB genes are components of basement membrane laminins; their cfRNA upregulation in CRC plasma may reflect shed ECM remodeling transcripts from tumor stroma.

**NES>0 — ERBB2/KIT/EGFR signaling (ranks 8, 10, 12, 15, 17):**

GRB2, SRC, LYN, FYN, HSP90AA1, SOS1, SHC1, PIK3R1, EGF, STAT1/5B, JAK2 appear across multiple receptor tyrosine kinase (RTK) pathway terms. ERBB2 amplification occurs in ~5-8% of CRC. KIT signaling is relevant in GI stromal tumors but may reflect stromal signaling in CRC microenvironment. The coherence of RTK signaling enrichment suggests cfRNA capture of tumor and tumor microenvironment (TME) signaling transcripts.

**NES>0 — Insulin receptor/glucose signaling (ranks 4, 6, 18):**

Metabolic reprogramming in CRC includes enhanced glucose uptake (Warburg effect). GLUT1/SLC2A1 upregulation is common in CRC. The insulin receptor/glucose import signaling enrichment may reflect the metabolic shift in CRC cells.

**NES<0 — Translation/ribosome biogenesis (dominant signal, NES up to -2.43):**

The most strongly downregulated pathways are global translation machinery: eukaryotic translation elongation/termination, cap-dependent initiation, 40S/60S ribosome assembly, rRNA processing, NMD, and SRP-dependent co-translational protein targeting. This is an extremely coherent ribosome/translation suppression signature. This mirrors exactly what Chen et al. (eLife 2022) reported: "ribosome biogenesis pathways were substantially downregulated" in cancer patients' cfRNA. The plos ONE 2024 CRC cfRNA study (PMC11326608) independently found "downregulated cfRNAs are highly enriched in pathways mainly related to ribosome biogenesis" and cited translation-related pathway suppression in tumor-educated platelets as a convergent finding. This likely reflects global suppression of ribosomal protein gene (RPG) expression in cfRNA derived from blood cells in the CRC milieu — a systemic effect of cancer on circulating cells, not direct tumor shedding.

**NES<0 — GCN2/amino acid deficiency, starvation response (ranks 12, 21):**

EIF2AK4/GCN2 is activated under amino acid deprivation and triggers the integrated stress response (ISR), globally suppressing translation. Its enrichment in the downregulated set is consistent with systemic nutrient stress in CRC patients.

**NES<0 — SLIT/ROBO regulation (rank 19):**

Unexplained by translation biology; likely co-ranked due to ribosomal protein gene overlaps within the SLIT/ROBO regulatory gene set.

**cfRNA background note:** "Viral mRNA Translation" and "Influenza Viral RNA Transcription" appearing in downregulated terms (ranks 3, 20) is a database artifact — these Reactome gene sets contain ribosomal subunit genes (RPS/RPL family) as shared members, not actual viral signatures.

### 3.4 Signal Quality Assessment

**SIGNAL QUALITY: HIGH**

Rationale: The NES<0 translation/ribosome suppression is a dominant, highly coherent signal (362/383 terms negative) validated independently by Chen et al. (eLife 2022) and published CRC cfRNA studies. The NES>0 Hippo/YAP signal is validated by CRC molecular biology. The 21:362 NES>0:NES<0 asymmetry indicates a fundamentally suppressed cfRNA transcriptome in CRC plasma relative to HC, with a narrow set of oncogenic pathway transcripts elevated — an interpretable, biologically meaningful pattern.

---

## 4. Esophagus Cancer

**Source datasets:** Chen et al. (eLife 2022) — esophageal cancer included in 5-cancer cfRNA study; Peddu et al. (bioRxiv 2025) — LOCATE-seq nanopore cfRNA for esophageal adenocarcinoma and Barrett's esophagus, identifying metabolic, signaling, and immune checkpoint pathways upregulated in precancer and cancer  
**Published findings:** Peddu et al. (2025) found metabolic, signaling, and immune checkpoint pathway upregulation in esophageal cancer/precancer cfRNA. FGFR2 signaling drives gastric/esophageal cancer progression and activates YAP1 downstream.

### 4.1 NES>0 Pathways (Anomalously Upregulated, n=142 significant)

**Top 30 by NES:**

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | 2.673 | <0.001 | GOBP: Embryonic Organ Development |
| 2 | 2.573 | <0.001 | Reactome: YAP1- And WWTR1(TAZ)-stimulated Gene Expression |
| 3 | 2.519 | <0.001 | GOBP: Adherens Junction Organization |
| 4 | 2.433 | 0.006 | GOBP: Gamma-Aminobutyric Acid Signaling Pathway |
| 5 | 2.426 | 0.004 | Reactome: Signaling By Type 1 IGF1R |
| 6 | 2.415 | 0.004 | KEGG: Nicotine addiction |
| 7 | 2.382 | 0.003 | Reactome: IGF1R Signaling Cascade |
| 8 | 2.371 | 0.003 | GOBP: Regulation Of Postsynaptic Membrane Potential |
| 9 | 2.332 | 0.002 | GOBP: Mammary Gland Development |
| 10 | 2.305 | 0.002 | Reactome: IRS-related Events Triggered By IGF1R |
| 11 | 2.300 | 0.002 | Reactome: PI-3K cascade:FGFR1 |
| 12 | 2.281 | 0.002 | GOBP: SMAD Protein Signal Transduction |
| 13 | 2.270 | 0.002 | GOBP: Sodium Ion Homeostasis |
| 14 | 2.259 | 0.005 | GOBP: Excitatory Postsynaptic Potential |
| 15 | 2.255 | 0.004 | GOBP: Negative Regulation Of Epithelial Cell Differentiation |
| 16 | 2.248 | 0.004 | GOBP: Positive Regulation Of Insulin Receptor Signaling Pathway |
| 17 | 2.245 | 0.004 | GOBP: Positive Regulation Of Glucose Import |
| 18 | 2.234 | 0.005 | Reactome: SHC-mediated cascade:FGFR2 |
| 19 | 2.231 | 0.005 | GOBP: ERK1 And ERK2 Cascade |
| 20 | 2.219 | 0.004 | GOBP: Protein Localization To Cell Junction |
| 21 | 2.216 | 0.004 | Reactome: Activated Point Mutants Of FGFR2 |
| 22 | 2.210 | 0.004 | GOBP: Regulation Of Monoatomic Ion Transport |
| 23 | 2.206 | 0.004 | Reactome: SHC-mediated cascade:FGFR1 |
| 24 | 2.198 | 0.006 | Reactome: Signal Attenuation |
| 25 | 2.191 | 0.005 | Reactome: FGFR2 Ligand Binding And Activation |
| 26 | 2.158 | 0.010 | GOBP: Negative Regulation Of miRNA Transcription |
| 27 | 2.158 | 0.010 | GOBP: Regulation Of Glucose Import |
| 28 | 2.151 | 0.010 | Reactome: PI-3K cascade:FGFR3 |
| 29 | 2.142 | 0.010 | GOBP: Sodium Ion Transmembrane Transport |
| 30 | 2.141 | 0.010 | GOBP: Positive Regulation Of Cellular Response To Insulin Stimulus |

**Key lead genes:**
- Embryonic organ development: IGF2, NES, TGFB1, TGFB2, SOX17/SOX18, HEY1/HEY2, KDR, FGFR2, TEAD1/TEAD2/TEAD4, GLI2, GATA3, FOXE1, WNT7B
- YAP/TAZ: HIPK2, WWTR1, HIPK1, YAP1, TEAD1/TEAD2/TEAD4, KAT2B
- Adherens junction: FERMT2, TJP1, CDH2/CDH5/CDH8/CDH13/CDH15, CTNNB1, VCL, SMAD7, JAM3
- FGFR2 (highly represented): FGF3/5/6/7/8/9/10/16/17/22/23, FGFBP1/2/3, FGFR2, PTPN11, IRS1/2/4, PI3K/AKT
- IGF1R: IGF2, FGFR2/3, IRS1/2/4, AKT2, PIK3R1, GRB2
- GABA/nicotine/ion channel: GABRA3/4/5, GABRG1/2/3, GRIN1/2A/2B/2C, CHRNA4/7, CACNA1A/1B

### 4.2 NES<0 Pathways (Anomalously Downregulated, n=539 significant)

**Top 30 by NES:**

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | -2.387 | <0.001 | GOBP: Mitochondrial Respiratory Chain Complex Assembly |
| 2 | -2.354 | <0.001 | Reactome: SRP-dependent Cotranslational Protein Targeting To Membrane |
| 3 | -2.351 | <0.001 | GOBP: Mitochondrial ATP Synthesis Coupled Electron Transport |
| 4 | -2.339 | <0.001 | GOBP: Cellular Respiration |
| 5 | -2.327 | <0.001 | Reactome: Translation |
| 6 | -2.326 | <0.001 | Reactome: Eukaryotic Translation Termination |
| 7 | -2.317 | <0.001 | GOBP: Cytoplasmic Translation |
| 8 | -2.313 | <0.001 | GOBP: Aerobic Electron Transport Chain |
| 9 | -2.313 | <0.001 | GOBP: Mitochondrial Gene Expression |
| 10 | -2.311 | <0.001 | GOBP: Ribosomal Small Subunit Biogenesis |
| 11 | -2.310 | <0.001 | Reactome: Respiratory Electron Transport |
| 12 | -2.307 | <0.001 | Reactome: Viral mRNA Translation |
| 13 | -2.300 | <0.001 | Reactome: Respiratory Electron Transport + ATP Synthesis (OXPHOS) |
| 14 | -2.298 | <0.001 | GOBP: Mitochondrial Translation |
| 15 | -2.297 | <0.001 | Reactome: rRNA Processing In Nucleus And Cytosol |
| 16 | -2.296 | <0.001 | KEGG: Oxidative phosphorylation |
| 17 | -2.295 | <0.001 | GOBP: Translation |
| 18 | -2.295 | <0.001 | Reactome: rRNA Processing |
| 19 | -2.294 | <0.001 | GOBP: Ribonucleoprotein Complex Biogenesis |
| 20 | -2.293 | <0.001 | Reactome: Formation Of A Pool Of Free 40S Subunits |
| 21 | -2.293 | <0.001 | Reactome: Eukaryotic Translation Elongation |
| 22 | -2.281 | <0.001 | GOBP: NADH Dehydrogenase Complex Assembly |
| 23 | -2.281 | <0.001 | GOBP: Mitochondrial Respiratory Chain Complex I Assembly |
| 24 | -2.280 | <0.001 | Reactome: Major Pathway Of rRNA Processing |
| 25 | -2.270 | <0.001 | Reactome: GTP Hydrolysis And Joining Of 60S Ribosomal Subunit |
| 26 | -2.266 | <0.001 | Reactome: Selenoamino Acid Metabolism |
| 27 | -2.262 | <0.001 | Reactome: Peptide Chain Elongation |
| 28 | -2.260 | <0.001 | Reactome: NMD Independent Of EJC |
| 29 | -2.258 | <0.001 | GOBP: Ribosome Biogenesis |
| 30 | -2.257 | <0.001 | Reactome: Cap-dependent Translation Initiation |

### 4.3 Biological Interpretation

**NES>0 — Embryonic organ development (top hit, NES=2.673):**

Lead genes include SOX17, SOX18 (endodermal transcription factors), TEAD1/2/4 (YAP effectors), WNT7B, FGFR2, GATA3, GLI2, HEY1/HEY2 (Notch targets), TGFB1/TGFB2, FOXE1, SALL1. This gene set captures a **developmental/stem-like reactivation program** in esophageal cancer. Esophageal cancer (both squamous cell and adenocarcinoma) is well-known to reactivate embryonic gene expression programs as part of epithelial-to-mesenchymal transition (EMT) and cancer stem cell formation. The presence of SOX17 (required for endodermal identity), FOXE1 (thyroid/foregut transcription factor mutated in cancer), and GATA3 (foregut and mammary lineage TF) suggests cfRNA shedding of dedifferentiated, embryonic-like tumor cells or their exosomes. TEAD1/2/4 connecting this term to YAP signaling (rank 2) is mechanistically significant.

**NES>0 — YAP1/TAZ signaling (rank 2, NES=2.573):**

Same YAP1, WWTR1, TEAD family lead genes as CRC. YAP/TAZ activation is a cancer hallmark across GI cancers. In esophageal cancer, YAP nuclear localization is a poor prognostic marker. Importantly, FGFR2 signaling has been shown to activate YAP1 downstream in gastric/esophageal cancer (PMC7581496 for gastric cancer; the FGF18-FGFR2-YAP1 axis). The co-enrichment of FGFR2 signaling (ranks 18, 21, 25) with YAP/TAZ (rank 2) suggests mechanistic connectivity: FGFR2-driven esophageal cancer activates YAP1, and both signatures are captured in cfRNA.

**NES>0 — FGFR2 signaling (highly represented, ranks 18, 21, 23, 25, 28):**

FGFR2 is amplified/overexpressed in gastric cancer (~10%) and esophageal adenocarcinoma. The FGF ligands leading these terms (FGF3/5/6/7/8/9/10/16/17/22/23) and FGFR2 itself, plus FGFBP1/2/3 (FGF-binding proteins that facilitate FGF action), represent a coherent FGFR2-centric oncogenic signal in cfRNA. Lead genes overlap strongly with IGF1R terms (IRS1/2, AKT2, PIK3R1, GRB2), indicating co-activation of growth factor receptor signaling networks. This is expected in esophageal adenocarcinoma where RTK signaling crosstalk is a driver.

**NES>0 — IGF1R signaling (ranks 5, 7, 10):**

IGF2 (imprinted oncogene) and the full IGF1R signaling cascade (IRS1/2/4, AKT2, PI3K, GRB2, SHC1) are upregulated. IGF2/IGF1R axis is relevant in Barrett's esophagus progression to adenocarcinoma. IGF2 is also a key embryonic growth factor explaining its appearance in the top Embryonic Organ Development term.

**NES>0 — GABA/nicotine/neurotransmitter signaling (ranks 4, 6, 8, 13, 14):**

GABRA3, GABRG1/2/3, GRIN1/2A/2B/2C, CHRNA4/7, CACNA1A/1B — ion channels and neurotransmitter receptors. This unexpected signature in esophageal cancer cfRNA may reflect: (a) ectopic expression of neuroendocrine/enteric nervous system genes in esophageal tumor (neuroendocrine differentiation), (b) enteric nervous system cell-derived cfRNA shed near the tumor, or (c) nicotine addiction pathway relevance given tobacco as a major esophageal cancer risk factor (KEGG Nicotine addiction, rank 6). Esophageal squamous cell carcinoma has strong tobacco association and evidence of nicotinic receptor expression. Sodium ion homeostasis (rank 13) links through the same ion channel gene pool.

**NES>0 — SMAD signaling (rank 12, NES=2.281):**

SMAD5, SMAD6, SMAD9, TGFB1/TGFB2, BMP2/5/10, GDF11, INHBB, NODAL — TGF-beta/BMP superfamily signaling. This is a recognized EMT driver in esophageal cancer. TGF-beta/SMAD signaling drives invasion. TGFBR2 loss-of-function mutations are found in microsatellite-unstable CRC and esophageal cancers. The cfRNA enrichment of SMAD pathway genes may reflect tumor stroma and EMT signaling.

**NES>0 — Adherens junction organization (rank 3, NES=2.519):**

CDH2 (N-cadherin), CDH5 (VE-cadherin), CDH8/13/15, CTNNB1 (beta-catenin), FERMT2 (kindlin-2), VCL (vinculin), SMAD7, JAM3. The presence of N-cadherin (CDH2) alongside E-cadherin-related junctional genes suggests the classic cadherin switch of EMT. CDH5/VE-cadherin indicates endothelial cell involvement (angiogenesis at tumor site). CTNNB1/beta-catenin connects to Wnt signaling. This term reflects ECM/EMT remodeling in esophageal cancer.

**NES<0 — OXPHOS + Translation (combined, n=539 terms downregulated):**

Esophageal cancer shows the **combined OXPHOS + ribosome/translation suppression** signature seen in HCC (OXPHOS) and CRC (translation) respectively. The NES<0 signal is split between:
- OXPHOS/mitochondrial respiration (ranks 1, 3, 4, 8, 9, 11, 13, 16, 22, 23): NDUF/COX/UQCR/ATP5 subunit genes — same as HCC
- Translation/ribosome (ranks 2, 5, 6, 7, 10, 12, 15, 17, 18, 19, 20, 21): RPS/RPL/EIF genes — same as CRC

This dual suppression in esophageal cancer is biologically coherent: esophageal cancer, like CRC, suppresses systemic translation in cfRNA, while also exhibiting the mitochondrial OXPHOS downregulation seen in HCC. The co-occurrence suggests these are general features of GI/upper GI cancers captured by the normative model.

**GOBP: Mitochondrial Gene Expression (rank 9, NES=-2.313):** Specifically covers mitochondrial ribosomal genes (MRPs) and mitochondrial tRNA/rRNA processing, bridging the OXPHOS and translation suppression themes.

### 4.4 Signal Quality Assessment

**SIGNAL QUALITY: HIGH**

Rationale: The NES>0 signals are mechanistically coherent (FGFR2-YAP1 axis, embryonic/stem-like reactivation, IGF2/IGF1R, SMAD/TGF-beta, EMT) and consistent with published esophageal cancer molecular pathology. The FGFR2 signal is particularly strong and specific. The dual OXPHOS+translation suppression in NES<0 is validated by parallel findings in HCC (OXPHOS) and CRC (translation). The largest dataset in this group (681 terms tested, 142 positive) provides good coverage.

---

## 5. Stomach Cancer

**Source dataset:** Chen et al. (eLife 2022) — stomach cancer included in 5-cancer cfRNA study  
**Published findings:** Chen et al. found ECM-receptor interactions and neutrophil extracellular traps upregulated, and ribosome pathways downregulated in gastric cancer cfRNA. FGFR2 amplification is a recognized driver in ~10% of gastric cancer (PMC7581496); YAP1 is downstream of FGFR2 in gastric cancer.

### 5.1 NES>0 Pathways (Anomalously Upregulated, n=17 significant)

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | 2.643 | <0.001 | GOBP: Cell-Cell Junction Organization |
| 2 | 2.583 | <0.001 | Reactome: Formation Of Cornified Envelope |
| 3 | 2.477 | <0.001 | GOBP: Positive Regulation Of Integrin-Mediated Signaling Pathway |
| 4 | 2.390 | <0.001 | Reactome: RET Signaling |
| 5 | 2.328 | 0.005 | GOBP: Extracellular Matrix Assembly |
| 6 | 2.317 | 0.004 | GOBP: Regulation Of Branching Involved In Ureteric Bud Morphogenesis |
| 7 | 2.312 | 0.007 | GOBP: Aortic Valve Morphogenesis |
| 8 | 2.281 | 0.010 | GOBP: Regulation Of Monoatomic Ion Transport |
| 9 | 2.267 | 0.012 | GOBP: Mammary Gland Development |
| 10 | 2.199 | 0.036 | GOBP: Positive Regulation Of Morphogenesis Of An Epithelium |
| 11 | 2.178 | 0.040 | GOBP: Intermediate Filament Organization |
| 12 | 2.162 | 0.037 | Reactome: Cell-Extracellular Matrix Interactions |
| 13 | 2.125 | 0.043 | GOBP: Body Fluid Secretion |
| 14 | 2.120 | 0.039 | GOBP: Positive Regulation Of Vascular Permeability |
| 15 | 2.105 | 0.043 | GOBP: Male Gonad Development |
| 16 | 2.078 | 0.047 | GOBP: Protein Activation Cascade |
| 17 | 2.078 | 0.047 | GOBP: Blood Coagulation, Fibrin Clot Formation |

**Key lead genes:**
- Cell-cell junction: PKP4, PKP3, PKP1 (desmosomal plakophilins), TJP1/TJP2, CDH3/5/6/9, CLDN1/CLDN5, DSG2, OCLN, CTNND1, MARVELD2, FSCN1, TGFB1/2, GJA4, GJC1
- Cornified envelope: PKP4, PKP1, PKP3, FURIN, CAPN1/CAPNS1, PCSK6, TCHH (trichohyalin), DSC3 (desmocollin), IVL (involucrin), PPL (periplakin), LCE2A/2D/1F/1B/6A (late cornified envelope), SPRR1B, TGM5, KLK12/14, DSG4, SPINK5, CASP14, EVPL (envoplakin), CDSN (corneodesmosin)
- Integrin: FLNA, DMTN, NID1, LIMS2, LAMB2
- RET: DOK2, PTPN11, SHANK3, PRKCA, RET, SRC, PIK3R3, GFRA3/4, GRB2, SHC1, IRS2

### 5.2 NES<0 Pathways (Anomalously Downregulated, n=503 significant)

**Top 30 by NES:** (Translation/ribosome + OXPHOS suppression, same structure as CRC and Esophagus)

| Rank | NES | FDR | Term |
|---|---|---|---|
| 1 | -2.284 | <0.001 | Reactome: Formation Of A Pool Of Free 40S Subunits |
| 2 | -2.281 | <0.001 | Reactome: Eukaryotic Translation Elongation |
| 3 | -2.280 | <0.001 | Reactome: Peptide Chain Elongation |
| 4 | -2.274 | <0.001 | Reactome: Viral mRNA Translation |
| 5 | -2.271 | <0.001 | Reactome: L13a-mediated Translational Silencing Of Ceruloplasmin Expression |
| 6 | -2.257 | <0.001 | Reactome: Selenocysteine Synthesis |
| 7 | -2.250 | <0.001 | Reactome: NMD Independent Of EJC |
| 8 | -2.248 | <0.001 | Reactome: Influenza Viral RNA Transcription And Replication |
| 9 | -2.246 | <0.001 | Reactome: GTP Hydrolysis And Joining Of 60S Ribosomal Subunit |
| 10 | -2.241 | <0.001 | GOBP: Cytoplasmic Translation |
| 11 | -2.236 | <0.001 | Reactome: Eukaryotic Translation Termination |
| 12 | -2.221 | <0.001 | GOBP: Ribosomal Small Subunit Biogenesis |
| 13 | -2.212 | <0.001 | Reactome: Selenoamino Acid Metabolism |
| 14 | -2.208 | <0.001 | Reactome: NMD Enhanced By EJC |
| 15 | -2.206 | <0.001 | Reactome: rRNA Processing |
| 16 | -2.204 | <0.001 | Reactome: Major Pathway Of rRNA Processing |
| 17 | -2.199 | <0.001 | Reactome: Cap-dependent Translation Initiation |
| 18 | -2.193 | <0.001 | Reactome: SRP-dependent Cotranslational Protein Targeting |
| 19 | -2.193 | <0.001 | Reactome: rRNA Processing In Nucleus And Cytosol |
| 20 | -2.189 | <0.001 | GOBP: Ribonucleoprotein Complex Biogenesis |
| 21 | -2.184 | <0.001 | KEGG: Ribosome |
| 22 | -2.171 | <0.001 | Reactome: Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency |
| 23 | -2.160 | <0.001 | GOBP: Anterograde Synaptic Vesicle Transport |
| 24 | -2.160 | <0.001 | GOBP: Synaptic Vesicle Transport Along Microtubule |
| 25 | -2.136 | <0.001 | Reactome: Translation Initiation Complex Formation |
| 26 | -2.134 | <0.001 | KEGG: Ribosome biogenesis in eukaryotes |
| 27 | -2.130 | <0.001 | GOBP: Aerobic Electron Transport Chain |
| 28 | -2.127 | <0.001 | GOBP: Proton Motive Force-Driven ATP Synthesis |
| 29 | -2.126 | <0.001 | GOBP: Oxidative Phosphorylation |
| 30 | -2.125 | <0.001 | Reactome: Cellular Response To Starvation |

### 5.3 Biological Interpretation

**NES>0 — Cornified envelope (rank 2, NES=2.583) — the most distinctive stomach cancer signal:**

The Cornified Envelope (CE) is the structural scaffold of the terminal differentiation program of keratinocytes/squamous epithelium. Lead genes: involucrin (IVL), periplakin (PPL), envoplakin (EVPL), SPRR1B (small proline-rich protein), corneodesmosin (CDSN), late cornified envelope proteins (LCE2A/2D/1B/1F/6A), transglutaminase 5 (TGM5), calpain 1/S1 (CAPN1/CAPNS1), desmosomal proteins (DSC3, DSG4, PKP1/3/4).

This is a striking and unexpected signal in gastric cancer cfRNA. Three explanations:

1. **Squamous cell histology co-represented:** The gastric cancer samples in Chen et al. (eLife 2022) likely include both adenocarcinoma and squamous-cell-like subtypes. Esophagogastric junction tumors (Siewert classification) frequently exhibit squamous differentiation. Squamous cell carcinoma of the stomach/cardia expresses cornification genes.

2. **Tumor microenvironment — squamous metaplasia:** Squamous metaplasia, where gastric or esophageal epithelium undergoes squamous differentiation under chronic irritation (H. pylori, acid reflux), produces cornification transcripts. These would be shed into plasma cfRNA.

3. **cfRNA from overlapping anatomical regions:** Gastric cancer patients may have co-occurring esophageal or skin/mucosal inflammatory cfRNA contributions. The overlapping signal with esophageal (where squamous epithelium is normal tissue) is plausible.

The plakophilins (PKP1/3/4) are desmosomal structural proteins; PKP4 (also known as p0071) is associated with epithelial junction formation. DSG2 (desmoglein) is a component of desmosomes found in gastric epithelium and frequently overexpressed in gastric cancer. The co-enrichment of Cell-Cell Junction Organization (rank 1, NES=2.643) with the same PKP/DSG/CDH lead genes confirms this is a coherent desmosomal/epithelial junction signal.

**NES>0 — Cell-cell junction organization (rank 1, NES=2.643):**

PKP4/3/1, TJP1/2, CDH3/5/6/9, CLDN1/5, DSG2, OCLN, CTNND1 (p120-catenin), MARVELD2 (tricellulin). This is a tight junction + desmosomal gene cluster. In gastric cancer, disruption of cell-cell junctions is paradoxically accompanied by upregulation of some junction components in tumor and stromal cells. CDH3 (P-cadherin) is overexpressed in gastric cancer and associated with poor prognosis. CDH5 (VE-cadherin) reflects tumor angiogenesis. The cfRNA from these junction proteins likely reflects tumor cell shedding and stromal/endothelial activation.

**NES>0 — RET signaling (rank 4, NES=2.390):**

RET (rearranged during transfection) proto-oncogene; lead genes: RET, GFRA3/GFRA4 (GFR-alpha co-receptors), DOK2/6, PTPN11, PRKCA, SRC, IRS2, GRB2, SHC1. RET is the receptor tyrosine kinase activated by GDNF family ligands. While classically associated with thyroid cancer and MEN2, RET rearrangements/fusions are found in ~0.5-2% of gastric cancers, and RET-positive tumors are responsive to targeted therapy. GFRA3/4 co-receptor upregulation suggests enteric nervous system signaling (GFRA3 is the co-receptor for Artemin/ARTN, which promotes gastric cancer invasion). The cfRNA enrichment of RET signaling is likely driven by both tumor cells expressing RET and enteric neurons proximate to the tumor.

**NES>0 — ECM assembly, integrin signaling (ranks 5, 3):**

FLNA, DMTN, NID1, LIMS2, LAMB2 for integrin; ECM assembly includes additional laminin/collagen genes. ECM remodeling is a hallmark of gastric cancer invasion. This mirrors the CRC integrin signal (shared lead genes: FLNA, DMTN, NID1, LIMS2, LAMB2 — nearly identical).

**NES>0 — Fibrin clot formation (rank 17, NES=2.078):**

Lead genes: GP1BB, FLNA. Only two lead genes; this FDR=0.047 hit is marginal. The coagulation activation in gastric cancer cfRNA is expected (cancer-associated coagulopathy), but the signal here is weak.

**NES>0 — Mammary gland development, Aortic valve morphogenesis, Ureteric bud morphogenesis, Male gonad development (ranks 9, 7, 6, 15):**

These developmental pathway terms at lower NES (2.10-2.27, FDR 0.01-0.05) represent ectopic expression of developmental transcription factors in gastric cancer. Ion transport (rank 8) shares lead genes with the nervous system/enteric signaling cluster. These developmental signals may reflect cancer stem cell dedifferentiation programs.

**NES<0 — Translation/ribosome suppression (dominant, ~480/503 negative terms):**

Identical pattern to CRC and esophageal cancer: complete suppression of cytoplasmic translation, rRNA processing, ribosome biogenesis in stomach cancer cfRNA. NES values (-2.28 to -2.16) are slightly weaker than CRC (-2.43 to -2.19) but the pattern is consistent. This appears to be a shared feature of GI cancers — a systemic suppression of ribosomal protein gene cfRNA in the plasma milieu of cancer patients.

**NES<0 — OXPHOS (ranks 27-29):**

Weaker OXPHOS signal in gastric cancer (NES=-2.13 to -2.13) compared to HCC (NES=-2.92) and esophageal cancer (NES=-2.39). OXPHOS suppression is present but less dominant than translation suppression, suggesting gastric cancer has a relatively more translation-centric cfRNA reduction phenotype.

**NES<0 — Synaptic vesicle transport (ranks 23, 24, NES=-2.16):**

Unusual in gastric cancer. Lead genes include KIF1A, KRAS-related vesicle trafficking genes. May reflect enteric nervous system signaling disruption in the gastric tumor environment. Alternatively this is a noise term at this NES level.

### 5.4 Signal Quality Assessment

**SIGNAL QUALITY: HIGH**

Rationale: The NES>0 cornified envelope/desmosomal junction signal is distinctive and biologically interpretable (squamous differentiation in gastric/cardia cancer). RET signaling is a known gastric cancer oncogenic pathway. ECM/integrin signals are shared cancer hallmarks. The NES<0 translation suppression is consistent with CRC findings. The tight 17:503 NES>0:NES<0 ratio (even more extreme than CRC) indicates a similarly global cfRNA translation suppression in gastric cancer, with a small set of tumor-specific signals elevated.

---

## Cross-Disease Comparison

### Shared NES<0 Signatures (Translation + OXPHOS suppression)

| Signal | Liver Cancer | Liver Cirrhosis | Colorectal Cancer | Esophagus Cancer | Stomach Cancer |
|---|---|---|---|---|---|
| OXPHOS/mitochondrial | **Primary** (NES -2.92) | Present (NES -2.28) | Absent | **Primary** (NES -2.39) | Secondary (NES -2.13) |
| Translation/ribosome | Absent | Absent | **Primary** (NES -2.43) | **Primary** (NES -2.39) | **Primary** (NES -2.28) |

- Liver cancers (HCC, cirrhosis) preferentially suppress **OXPHOS** cfRNA
- GI cancers (CRC, stomach) preferentially suppress **translation/ribosome** cfRNA
- Esophageal cancer suppresses **both equally** (intermediate anatomical position between liver and lower GI)

### Shared NES>0 Signatures

| Signal | Liver Cancer | Liver Cirrhosis | Colorectal Cancer | Esophagus Cancer | Stomach Cancer |
|---|---|---|---|---|---|
| Complement/coagulation | **Primary** | Present | Absent | Absent | Minor (fibrin) |
| Lipoprotein/apolipoprotein | **Primary** | Present (HDL) | Absent | Absent | Absent |
| YAP/Hippo | Absent | Absent | **Top** | **Top** | Absent |
| FGFR2/RTK signaling | Absent | Absent | ERBB2/KIT | **FGFR2-dominant** | RET |
| ECM/integrin | Minor | Absent | Present | Present | Present |
| Desmosomal/junction | Absent | Absent | Present | Present | **Primary** |
| Histone/chromatin | Minor | **Primary** | Absent | Absent | Absent |

### cfRNA Background Artifacts

1. **Histone transcripts (H2AC, H2BC, H3C series):** Predominantly appear in Liver Cirrhosis (top hit: HDAC pathway, NES=2.363). Also embedded in Systemic Lupus Erythematosus and Amyloid Fiber Formation terms in cirrhosis. Present at lower levels in other diseases as minor contributors to non-specific terms. These are products of nucleosome-associated cfRNA from dying cells and NETosis, not disease-specific signals. Their prominence in cirrhosis (but not liver cancer) suggests the histone cfRNA background is amplified in the inflammatory neutrophil activation seen in portal hypertension rather than being tumor-specific.

2. **Viral pathway terms as ribosomal surrogates:** "Viral mRNA Translation," "Influenza Viral RNA Transcription and Replication," "SARS-CoV-2 Modulates Host Translation Machinery" appear in NES<0 terms for CRC, esophageal, and stomach cancer. These Reactome/KEGG gene sets contain large overlapping pools of ribosomal subunit genes (RPS/RPL), eIF genes, and RNA polymerase components. They are ranking negatively because of shared gene membership with ribosome/translation gene sets, not because of true viral signatures in cancer patients.

3. **Developmental/morphogenesis terms at moderate NES (aortic valve, mammary gland, ureteric bud):** These terms (present in esophageal and stomach cancer NES>0, FDR 0.01-0.05) reflect ectopic expression of developmental transcription factors in cancer. While not pure artifacts, they should be interpreted as reflecting cancer stem cell/EMT programs rather than literal developmental biology.

---

## Summary Assessment Table

| Disease | Signal Quality | Primary NES>0 Biological Meaning | Primary NES<0 Biological Meaning | Key Caveat |
|---|---|---|---|---|
| **Liver Cancer** | HIGH | Hepatocyte acute-phase/secreted protein release (complement, coagulation, apolipoproteins, heme); liver-specific cancer signal validated by Roskams-Hieter 2022 | OXPHOS/mitochondrial suppression (Warburg effect); MHC-II immune suppression | 100% significant terms; OXPHOS NES strongest in entire dataset (-2.92) |
| **Liver Cirrhosis** | LOW-MODERATE | Histone/NET background + endothelial migration (portal hypertension); complement shared with HCC | OXPHOS (early mitochondrial dysfunction) | **n=4 samples — severely underpowered; interpret with caution** |
| **Colorectal Cancer** | HIGH | YAP/Hippo oncogenic signaling; ERBB2/KIT/RTK; ECM/integrin; insulin signaling | Global translation/ribosome suppression (validated by Chen 2022 and independent CRC cfRNA studies) | 362/383 terms negative; near-complete translational suppression signal |
| **Esophagus Cancer** | HIGH | FGFR2-YAP1 axis; embryonic/stem-like reactivation (SOX17, TEAD, GLI2); IGF1R; SMAD/TGF-beta; adherens junction/EMT | Dual OXPHOS + translation suppression (intermediate between HCC and CRC) | Largest dataset (681 terms); FGFR2 signal highly represented across multiple terms |
| **Stomach Cancer** | HIGH | Cornified envelope/desmosomal junctions (squamous differentiation); RET signaling; ECM; cell-cell junction | Global translation/ribosome suppression (stronger than OXPHOS; mirrors CRC) | Cornified envelope top signal likely reflects squamous-type gastric/cardia tumors or GEJ overlap |

---

*Evaluation completed: 2026-06-24*  
*Analysis: GAMLSS normative model GSEA-prerank on mean Z-score; FDR<0.05; Reactome 2022 + GO BP 2023 + KEGG 2021 Human*
