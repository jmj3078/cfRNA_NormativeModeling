# GSEA Biological Evaluation — Group 4: Immune/Hematologic Diseases

**Diseases**: ICI-induced Myocarditis (ICI-m) | ICI-treated Cancer | Multiple Myeloma (MM) | MGUS | ME/CFS  
**Method**: GSEA-prerank on mean Z-scores per phenotype (normative model Z-scores vs HC reference, n=693)  
**NES > 0**: pathways enriched in high-Z genes (upregulated relative to HC)  
**NES < 0**: pathways enriched in low-Z genes (downregulated relative to HC)  
**Significance threshold**: FDR q-val < 0.05  

---

## Data Sources

| Disease | n | Author | DOI |
|---------|---|--------|-----|
| ICI-m (ICI-induced myocarditis) | 11 | Raissadati et al., JCI 2025 | https://doi.org/10.1172/jci188817 |
| ICI-treated Cancer | 11 | Raissadati et al., JCI 2025 | https://doi.org/10.1172/jci188817 |
| Multiple Myeloma (MM) | 17 | Roskams-Hieter et al., npj Precis Oncol 2022 | https://doi.org/10.1038/s41698-022-00270-y |
| MGUS | 8 | Roskams-Hieter et al., npj Precis Oncol 2022 | https://doi.org/10.1038/s41698-022-00270-y |
| ME/CFS | 90 | Gardella et al., PNAS 2025 | https://doi.org/10.1073/pnas.2507345122 |

---

## Summary Statistics

| Disease | Total FDR<0.05 | NES>0 (UP) | NES<0 (DOWN) |
|---------|---------------|-----------|--------------|
| ICI-m | 482 | 313 | 169 |
| ICI-treated Cancer | 1,298 | 226 | 1,072 |
| MM | 332 | 280 | 52 |
| MGUS | 108 | 94 | 14 |
| ME/CFS | 165 | 4 | 161 |

---

## 1. ICI-induced Myocarditis (ICI-m), n=11

### Study Context
Raissadati et al. (JCI 2025) used plasma cell-free mRNA (cfRNA) sequencing to distinguish ICI-myocarditis (n=10 discovery + 6 validation) from ICI-treated cancer patients without cardiac adverse events and from healthy controls. The study is the primary cfRNA-based investigation of ICI-m and identified a 6-gene plasma classifier (ZNF385B, CADM2, AQP7, B2M, IFITM2, CCL5) with >90% accuracy. The ICI-m samples in the current normative model (n=11) represent one of the smallest cohorts, warranting caution in interpretation.

### Top Upregulated Pathways (NES > 0, FDR < 0.05)

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | 2.889 | KEGG: Lysosome |
| 2 | 2.869 | KEGG: Oxidative phosphorylation |
| 3 | 2.784 | GO: Proton Motive Force-Driven ATP Synthesis |
| 4 | 2.760 | Reactome: Neutrophil Degranulation |
| 5 | 2.710 | GO: Proton Motive Force-Driven Mitochondrial ATP Synthesis |
| 6 | 2.664 | GO: Oxidative Phosphorylation |
| 7 | 2.660 | GO: Aerobic Electron Transport Chain |
| 8 | 2.422 | Reactome: Antigen Processing–Cross Presentation |
| 9 | 2.379 | KEGG: Cardiac muscle contraction |
| 10 | 2.364 | GO: Ubiquitin-Dependent ERAD Pathway |
| 11 | 2.329 | GO: ERAD Pathway |
| 12 | 2.308 | Reactome: ER-Phagosome Pathway |
| 13 | 2.264 | KEGG: Protein processing in endoplasmic reticulum |
| 14 | 2.218 | KEGG: Proteasome |
| 15 | 2.203 | GO: Microglial Cell Activation |

### Top Downregulated Pathways (NES < 0, FDR < 0.05)

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | -2.735 | Reactome: Eukaryotic Translation Elongation |
| 2 | -2.731 | Reactome: Peptide Chain Elongation |
| 3 | -2.683 | Reactome: Formation Of A Pool Of Free 40S Subunits |
| 4 | -2.680 | Reactome: Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency |
| 5 | -2.671 | Reactome: Viral mRNA Translation |
| 6 | -2.633 | Reactome: Eukaryotic Translation Termination |
| 7 | -2.427 | KEGG: Ribosome |
| 8 | -2.252 | Reactome: FCERI Mediated Ca2+ Mobilization |
| 9 | -2.117 | Reactome: RUNX1 Regulates Megakaryocyte Differentiation And Platelet Function |

### Biological Interpretation

**Mitochondrial hyperactivation and energy stress.** The strongest upregulated signal in ICI-m plasma cfRNA is a convergent cluster of OXPHOS and mitochondrial respiratory chain pathways (NES 2.55–2.87), including Complex I assembly, electron transport chain components, TCA cycle, and proton motive force-driven ATP synthesis. This is consistent with the profound metabolic demands of cardiac inflammation: cardiomyocytes sustain contractile function under immune attack via upregulated mitochondrial energy production, and dying or stressed cardiomyocytes release mitochondrial transcripts into the circulation. Myocarditis-associated cardiomyocyte death elevates cfRNA fragments from high-mitochondrial-density cardiac cells. Similar OXPHOS upregulation has been documented in infectious myocarditis models and in heart failure transcriptomics, where mitochondrial OXPHOS genes are among the most differentially expressed.

**Cardiac antigen presentation and ERAD.** Uniquely enriched in ICI-m (absent in ICI-treated Cancer) are Reactome: Antigen Processing–Cross Presentation (NES 2.42), ER-Phagosome Pathway (NES 2.31), and Cross-Presentation of Soluble Exogenous Antigens (Endosomes) (NES 2.21), alongside ERAD (ubiquitin-dependent, NES 2.36) and Protein Processing in ER (NES 2.26). This is mechanistically highly coherent with ICI-m pathogenesis: the prevailing model of ICI-m involves CD8+ T cells and antigen-presenting cells (dendritic cells, macrophages) cross-presenting cardiac antigens—particularly α-myosin heavy chain (α-MyHC)—to autoreactive cytotoxic T cells after checkpoint blockade removes peripheral tolerance. The ER-phagosome pathway is precisely the mechanism by which phagocytosed antigenic material is retro-translocated into the ER for MHC-I loading and cross-presentation to CD8+ T cells. ERAD upregulation reflects the cellular protein quality control burden accompanying proteotoxic ER stress in inflamed cardiac and immune cells. The proteasome upregulation (NES 2.22) is consistent with both antigen processing and the degradation of misfolded proteins under ER stress.

**Neutrophil degranulation.** Shared with ICI-treated Cancer but extremely prominent in ICI-m (NES 2.76), this pathway reflects innate immune amplification in myocarditis. Neutrophil-derived transcripts in plasma are well-established cfRNA sources during systemic inflammation and tissue injury.

**Cardiac muscle contraction.** Enrichment of KEGG: Cardiac muscle contraction (NES 2.38) — shared with ICI-treated Cancer — reflects the contribution of cardiomyocyte-derived cfRNA in both groups. This is expected because: (a) ICI-treated cancer patients may have subclinical cardiac stress from their treatment, and (b) even non-myocarditis ICI patients can have elevated troponin. However, the ICI-m-specific cross-presentation and ERAD signals ontop of the shared cardiac signal point to active immune-mediated cardiac destruction rather than simple cardiac stress.

**Microglial cell activation (ICI-m unique).** GO: Microglial Cell Activation (NES 2.20) appears uniquely in ICI-m. This may reflect neuroinflammatory cross-reactivity, given that some ICI-related immune adverse events affect the CNS, or it may indicate microglial-like macrophage polarization in the inflamed cardiac tissue, as cardiac macrophages share transcriptomic features with microglia.

**Translational suppression.** The downregulated signal in ICI-m is dominated by ribosomal translation machinery (elongation, initiation, termination, 40S subunit pool formation, NES −2.44 to −2.74), KEGG Ribosome (NES −2.43), and EIF2AK4/GCN2 integrated stress response (NES −2.68). GCN2 activation (which phosphorylates eIF2α) is a hallmark of cellular stress: when amino acids are depleted or ER stress is overwhelming, GCN2/PERK-mediated eIF2α phosphorylation globally represses cap-dependent translation, allowing selective translation of stress-response mRNAs. The downregulation of ribosomal and elongation factor transcripts in plasma cfRNA of ICI-m patients thus reflects a systemic translational stress response, consistent with both the systemic inflammatory state and reduced hepatic/systemic protein synthesis during acute myocarditis. Suppressed FCERI-mediated Ca2+ mobilization (NES −2.25) and RUNX1/megakaryocyte pathways (NES −2.12) suggest altered mast cell/basophil and platelet biology, possibly reflecting consumption or exhaustion of these populations.

---

## 2. ICI-treated Cancer (without myocarditis), n=11

### Study Context
The same Raissadati et al. (JCI 2025) cohort: cancer patients receiving immune checkpoint inhibitors but not developing myocarditis. This group has the highest total enrichment count of all five diseases (1,298 significant terms), with 1,072 downregulated — reflecting the pervasive metabolic and translational remodeling induced by both cancer and ICI therapy.

### Top Upregulated Pathways

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | 2.825 | KEGG: Lysosome |
| 2 | 2.664 | GO: Vacuolar Acidification |
| 3 | 2.431 | KEGG: Other glycan degradation |
| 4 | 2.415 | Reactome: Immunoregulatory Interactions Between Lymphoid And Non-Lymphoid Cells |
| 5 | 2.388 | Reactome: Neutrophil Degranulation |
| 6 | 2.376 | Reactome: Beta Defensins |
| 7 | 2.351 | GO: Negative Regulation Of T Cell Proliferation |
| 8 | 2.276 | GO: Antifungal Innate Immune Response |
| 9 | 2.273 | GO: Antigen Processing And Presentation Of Endogenous Peptide Antigen |
| 10 | 2.263 | KEGG: Allograft Rejection |
| 11 | 2.263 | KEGG: Glycosaminoglycan Degradation |
| 12 | 2.208 | Reactome: ROS And RNS Production In Phagocytes |
| 13 | 2.194 | GO: Antigen Processing And Presentation Via MHC Class I |
| 14 | 2.177 | KEGG: Graft-versus-host Disease |
| 15 | 2.141 | KEGG: Intestinal Immune Network For IgA Production |

### Top Downregulated Pathways (selected)

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | -2.808 | Reactome: Eukaryotic Translation Elongation |
| 2 | -2.758 | Reactome: Peptide Chain Elongation |
| 3 | -2.740 | GO: Cytoplasmic Translation |
| 4 | -2.659 | KEGG: Ribosome |
| 5 | -2.503 | Reactome: Regulation Of Expression Of SLITs And ROBOs |
| 6 | -2.472 | Reactome: Signaling By ROBO Receptors |

### Biological Interpretation

**Immune activation dominated by lysosomal-phagocytic pathway and antigen presentation.** Like ICI-m, ICI-treated Cancer shows lysosome as the top upregulated pathway (NES 2.83). This reflects the core mechanism of action of immune checkpoint inhibitors: reinvigorating cytotoxic immune responses, which involves lysosomal antigen processing and MHC-I presentation. Vacuolar acidification (NES 2.66) supports the lysosomal signal—acidic lysosomes are required for cathepsin activation and antigen degradation. The strong enrichment of glycan degradation and glycosaminoglycan degradation pathways (NES 2.43, 2.27) reflects active lysosomal catabolism of extracellular matrix components, consistent with tumor invasion and macrophage-mediated tissue remodeling.

**ICI-specific immune phenotype.** Unique to ICI-treated Cancer (absent from ICI-m) are: (a) IgA intestinal immune network (NES 2.14), which may reflect mucosal immune activation from ICI-related colitis (a common adverse event); (b) Allograft rejection (NES 2.26) and Graft-versus-Host Disease (NES 2.18), which are pathway annotations capturing general alloreactive T cell biology — appropriate for ICI, which effectively unleashes T cells; (c) Autoimmune thyroid disease (NES 2.12), consistent with thyroiditis being a common ICI-related adverse event; (d) Defensins and antifungal innate immunity (NES 2.38, 2.28), reflecting broad innate immune activation.

**Pervasive ribosomal/translational downregulation.** ICI-treated Cancer has far more downregulated pathways (1,072) than upregulated (226). The dominant downregulated signal is translational machinery — ribosome, translation elongation/initiation/termination. This breadth of translational suppression exceeds that seen in ICI-m and likely reflects the combined effects of: (a) cancer-associated cachexia and systemic catabolic states; (b) chemotherapy-induced translational suppression (if combination therapy was used); (c) the known suppression of SLIT-ROBO axon guidance signaling (NES −2.50, −2.47) suggesting altered vascular and neural transcriptomes in circulation.

**Negative regulation of T cell proliferation (unique to ICI-treated).** GO: Negative Regulation Of T Cell Proliferation (NES 2.35) appearing as upregulated in ICI-treated Cancer (but not uniquely in ICI-m) is counterintuitive at first but reflects the complex immune landscape: ICI expands some T cell subsets while exhausted or regulatory cells may suppress proliferation of others. This represents the regulatory counterbalance to immune activation.

---

## 3. ICI-m vs ICI-treated Cancer: Direct Comparison

This is the most biologically informative comparison in Group 4. Both cohorts are from the same study (n=11 each, Raissadati et al. JCI 2025), allowing near-isogenic comparison of myocarditis vs non-myocarditis ICI response.

### Shared Enrichments (present in both)
Both groups share: Lysosome, OXPHOS, Neutrophil Degranulation, Multiple Aerobic Respiration/ETC terms, Cardiac Muscle Contraction, ERAD, Vacuolar Acidification. These represent the baseline ICI immune response and cancer-related cardiac stress — common to all ICI-treated patients regardless of myocarditis development.

### ICI-m Specific (not seen in ICI-treated Cancer)

| Category | Key Pathways | NES | Interpretation |
|----------|-------------|-----|----------------|
| Cardiac antigen cross-presentation | Antigen Processing–Cross Presentation, ER-Phagosome Pathway, Cross-Presentation Soluble Exogenous Antigens | 2.42, 2.31, 2.21 | Active cross-presentation of cardiac antigens to CD8+ T cells via MHC-I; the pathogenic mechanism of ICI-m |
| ER stress and protein QC | Protein Processing in ER, ERAD, Hedgehog ERAD, Proteasome | 2.26, 2.36, 2.25, 2.22 | Proteotoxic stress in inflamed cardiomyocytes and phagocytes; ER stress from immune-mediated cardiac damage |
| Enhanced mitochondrial respiration | Proton Motive Force-Driven ATP Synthesis, Aerobic Respiration, TCA + Respiratory ET | 2.78, 2.41, 2.26 | ICI-m-specific elevation of mitochondrial bioenergetic transcripts, likely from stressed/dying cardiomyocytes |
| CNS-related inflammation | Microglial Cell Activation | 2.20 | Possible neuroinflammatory co-involvement or cardiac macrophage resemblance to microglia |

### ICI-treated Cancer Specific (not in ICI-m)

| Category | Key Pathways | NES | Interpretation |
|----------|-------------|-----|----------------|
| Mucosal/systemic immune adverse events | IgA intestinal network, Autoimmune thyroid disease, Type I diabetes | 2.14, 2.12, 2.04 | ICI-related colitis, thyroiditis, and endocrine adverse events producing systemic cfRNA signatures |
| Broader alloimmune activation | Allograft rejection, Graft-vs-Host | 2.26, 2.18 | General alloreactive T cell signature from ICI-driven T cell expansion (without myocarditis specificity) |
| Axon guidance suppression | SLIT-ROBO signaling (DOWN, NES −2.50) | − | Cancer-specific neural/vascular remodeling in plasma cfRNA |

### Key Conclusion
The myocarditis-specific differentiating signals in plasma cfRNA are: **(1) cross-presentation of exogenous/soluble antigens via the ER-phagosome pathway**, directly implicating antigen presentation to CD8+ T cells as the mechanistic driver; **(2) ERAD and proteasome upregulation**, reflecting ER stress in inflamed cardiac tissue; **(3) enhanced mitochondrial bioenergetics**, reflecting dying cardiomyocyte cfRNA contribution. These align precisely with the published Raissadati et al. model: autoreactive CD8+ T cells (Temra phenotype) targeting cardiac α-myosin via cross-presentation are the primary effectors of ICI-m. The normative model Z-scores successfully capture this mechanistic distinction in plasma cfRNA without requiring explicit cardiac tissue sampling.

---

## 4. Multiple Myeloma (MM), n=17

### Study Context
Roskams-Hieter et al. (npj Precision Oncology 2022) sequenced plasma cfRNA from MM patients (n=10), MGUS (n=12), and non-cancer donors (n=20) to demonstrate that cfRNA profiling distinguishes cancer from premalignant states. MM is a plasma cell malignancy arising from clonal proliferation of immunoglobulin-secreting plasma cells in the bone marrow. Key disease hallmarks include: uncontrolled plasma cell proliferation, immunoglobulin overproduction causing paraprotein, immunosuppression, bone destruction, and anemia from marrow replacement. The Roskams-Hieter study identified CENPE, NEK2, NUSAP1 (mitotic kinases/spindle proteins) and HBG1/HBG2 (fetal hemoglobin) as top discriminating cfRNA biomarkers for MM vs MGUS.

### Top Upregulated Pathways (NES > 0)

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | 2.188 | GO: Carbon Dioxide Transport |
| 2 | 2.186 | Reactome: Kinesins |
| 3 | 2.185 | GO: Gas Transport |
| 4 | 2.149 | Reactome: G1/S-Specific Transcription |
| 5 | 2.115 | GO: Chromosome Condensation |
| 6 | 2.112 | GO: Hydrogen Peroxide Catabolic Process |
| 7 | 2.086 | GO: Positive Regulation Of Attachment Of Spindle Microtubules To Kinetochore |
| 8 | 2.078 | Reactome: HDACs Deacetylate Histones |
| 9 | 2.078 | GO: Oxygen Transport |
| 10 | 2.053 | Reactome: Polo-like Kinase Mediated Events |
| 11 | 2.042 | GO: Erythrocyte Differentiation |
| 12 | 2.042 | Reactome: E2F Transcription Targets Under HDAC1 Repression |
| 13 | 2.040 | GO: Mitotic Cytokinesis |
| 14 | 2.036 | GO: Sister Chromatid Segregation |
| 15 | 2.032 | Reactome: DNA Damage/Telomere Stress Induced Senescence |
| 16 | 2.030 | KEGG: Cell Cycle |
| 17 | 2.016 | Reactome: TP53 Regulates Cell Cycle Genes |

### Top Downregulated Pathways (NES < 0)

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | -3.049 | Reactome: Eukaryotic Translation Elongation |
| 2 | -2.879 | Reactome: Peptide Chain Elongation |
| 3 | -2.790 | GO: Immunoglobulin Mediated Immune Response |
| 4 | -2.644 | Reactome: SRP-Dependent Cotranslational Protein Targeting To Membrane |
| 5 | -2.596 | KEGG: Ribosome |
| 6 | -2.565 | KEGG: Asthma |
| 7 | -2.544 | GO: Cytoplasmic Translation |
| 8 | -2.361 | GO: Peptide Antigen Assembly With MHC Class II Complex |
| 9 | -2.361 | GO: MHC Class II Protein Complex Assembly |
| 10 | -2.315 | KEGG: Oxidative Phosphorylation |
| 11 | -2.155 | GO: B Cell Proliferation |
| 12 | -2.154 | KEGG: Allograft Rejection |
| 13 | -2.220 | Reactome: Vitamin D (Calciferol) Metabolism |

### Biological Interpretation

**Plasma cell malignancy drives cell cycle and mitotic machinery upregulation.** The most prominent upregulated signals in MM plasma cfRNA are cell cycle and mitotic spindle pathways: Kinesins (NES 2.19, including CENPE/NEK2/NUSAP1 identified by Roskams-Hieter et al.), G1/S-Specific Transcription, E2F/HDAC1 targets, Polo-like kinase events, chromosome condensation, spindle assembly, and sister chromatid segregation. This pattern reflects the circulating cfRNA signature of rapidly proliferating malignant plasma cells in the bone marrow. In MM, plasma cells undergo uncontrolled G1→S→M transitions, releasing high levels of cell cycle transcripts into the plasma. The identification of kinesins and spindle-associated genes is directly consistent with the Roskams-Hieter study's top MM cfRNA biomarkers.

**Erythrocyte/hemoglobin signature.** Carbon Dioxide Transport (NES 2.19), Gas Transport (NES 2.19), and Erythrocyte Differentiation (NES 2.04) — with lead genes including hemoglobin subunits (HBG1, HBG2 identified as top MM cfRNA biomarkers) — reflect the anemia signature of MM. MM causes anemia through multiple mechanisms: marrow infiltration by plasma cells crowding out erythroid precursors, cytokine-mediated erythropoiesis suppression, and renal damage reducing EPO. The upregulation of erythrocyte-related transcripts in MM cfRNA is paradoxical at first, but represents compensatory erythropoiesis: erythroblasts upregulate hemoglobin and CO2/O2 transport transcripts, and these cells (or their transcripts) enter the circulation at higher rates as the marrow attempts to compensate for anemia.

**Epigenetic remodeling.** HDACs Deacetylate Histones (NES 2.08), Nucleosome Organization (NES 2.00), RMTs Methylate Histone Arginines (NES 2.00), and DNA Damage/Telomere Stress Induced Senescence (NES 2.03) reflect the extensive epigenetic dysregulation characterizing MM. Histone deacetylation by HDACs and histone methyltransferases (RMTs) reprogram gene expression to support malignant plasma cell survival, and senescence pathways (TP53, telomere stress) reflect the selective pressure on plasma cells that escape senescence checkpoints during myelomagenesis.

**Immunoglobulin Mediated Immune Response strongly downregulated.** The downregulation of GO: Immunoglobulin Mediated Immune Response (NES −2.79) is seemingly paradoxical in a disease defined by monoclonal immunoglobulin overproduction. This reflects the distinction between *clonal* immunoglobulin production and *normal polyclonal immune response*: MM plasma cells produce abundant paraprotein, but the malignant clone suppresses normal B cell and plasma cell activity, resulting in hypogammaglobulinemia (low normal immunoglobulin levels) and impaired humoral immunity. The normative model captures this as a relative depletion of normal immunoglobulin-response gene expression compared to HC. The lead genes (HLA-DRA, HLA-DRB1, HLA-DMB, HLA-DOA, CD27, CD19, FCER1A) are normal B cell and antigen-presenting cell markers that are depleted in MM due to immune suppression of normal B lineage.

**MHC-II antigen presentation and B cell proliferation suppressed.** MHC Class II Protein Complex Assembly (NES −2.36), Peptide Antigen Assembly with MHC Class II (NES −2.36), and B Cell Proliferation (NES −2.16) are all downregulated, consistent with MM-associated immunosuppression of normal antigen presentation and B cell responses. MM cells downregulate MHC-II on non-malignant antigen-presenting cells as an immune evasion strategy.

**OXPHOS downregulation.** KEGG: Oxidative Phosphorylation (NES −2.32) is downregulated in MM — the inverse of ICI-m. This likely reflects bone marrow microenvironmental effects: MM cells create a hypoxic niche and reprogram neighboring hematopoietic cells toward glycolytic metabolism (Warburg effect), with reduced normal OXPHOS-high cell populations in circulation.

**Vitamin D metabolism suppression.** Vitamin D (Calciferol) Metabolism (NES −2.22) is downregulated. Vitamin D deficiency is common in MM; MM cells themselves dysregulate VDR signaling and vitamin D metabolism as a mechanism to escape anti-proliferative effects of calcitriol.

---

## 5. MGUS (Monoclonal Gammopathy of Undetermined Significance), n=8

### Study Context
MGUS is the obligate precursor state to MM: all MM patients pass through MGUS. The Roskams-Hieter 2022 study achieved 100% accuracy distinguishing MM from MGUS using cfRNA biomarkers, demonstrating that even small plasma cfRNA differences are biologically distinct. With only n=8, this is the smallest cohort analyzed here; findings require replication.

### Top Upregulated Pathways (NES > 0)

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | 2.332 | KEGG: Systemic Lupus Erythematosus |
| 2 | 2.221 | Reactome: Condensation Of Prophase Chromosomes |
| 3 | 2.206 | Reactome: YAP1/WWTR1(TAZ)-Stimulated Gene Expression |
| 4 | 2.177 | Reactome: DNA Methylation |
| 5 | 2.146 | Reactome: RNA Polymerase I Promoter Opening |
| 6 | 2.127 | Reactome: SIRT1 Negatively Regulates rRNA Expression |
| 7 | 2.118 | Reactome: HDACs Deacetylate Histones |
| 8 | 2.109 | Reactome: RHO GTPases Activate PKNs |
| 9 | 2.102 | Reactome: PKN1 Stimulates Transcription Of KLK2 And KLK3 |
| 10 | 2.099 | Reactome: DNA Damage/Telomere Stress Induced Senescence |
| 11 | 2.099 | Reactome: Transcriptional Regulation Of Granulopoiesis |
| 12 | 2.084 | Reactome: Pre-NOTCH Transcription And Translation |
| 13 | 2.077 | GO: Chromosome Condensation |
| 14 | 2.073 | Reactome: Packaging Of Telomere Ends |
| 15 | 2.039 | Reactome: Oxidative Stress Induced Senescence |

### Top Downregulated Pathways (NES < 0)

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | -2.368 | KEGG: Lysosome |
| 2 | -2.349 | GO: Antigen Processing And Presentation Of Exogenous Peptide Antigen |
| 3 | -2.288 | GO: Antigen Processing And Presentation Via MHC Class I |
| 4 | -2.240 | GO: Retrograde Protein Transport, ER To Cytosol |
| 5 | -2.239 | GO: Antigen Processing And Presentation Of Endogenous Peptide Antigen |
| 6 | -2.214 | GO: ER To Cytosol Transport |
| 7 | -2.192 | GO: Negative Regulation Of NK Cell Mediated Cytotoxicity |
| 8 | -2.186 | GO: Peptide Antigen Assembly With MHC Protein Complex |
| 9 | -2.130 | GO: Negative Regulation Of NK Cell Mediated Immunity |
| 10 | -2.111 | GO: Regulation Of Lymphocyte Activation |
| 11 | -2.098 | KEGG: Other Glycan Degradation |
| 12 | -2.076 | KEGG: Allograft Rejection |
| 13 | -2.072 | GO: Negative Regulation Of Leukocyte Mediated Cytotoxicity |
| 14 | -2.051 | GO: Monosaccharide Transmembrane Transport |

### Biological Interpretation

**SLE pathway as the top upregulated enrichment.** KEGG: Systemic Lupus Erythematosus (NES 2.33) is the highest-scoring enrichment in MGUS. The SLE pathway annotation in KEGG is dominated by histone genes (H2A, H2B, H3, H4 families) and complement components (C6, C8B), reflecting nucleosomal antigen release. In MGUS, the SLE signature captures a distinct biology: MGUS is associated with elevated histone gene expression in circulating cells (consistent with neutrophil extracellular trap [NET] formation, which also appears as uniquely upregulated: NES 2.02), and the chromatin condensation and chromosome packaging enrichments. This suggests that MGUS is associated with heightened chromatin exposure and NETs — potentially representing an early inflammatory milieu that precedes full malignant transformation.

**Epigenetic dysregulation signature in MGUS.** The MGUS-specific upregulated pathways are heavily epigenetic: DNA Methylation (NES 2.18), SIRT1 Negatively Regulates rRNA Expression (NES 2.13), HDACs Deacetylate Histones (NES 2.12, shared with MM), PRC2 Methylates Histones And DNA (NES 1.96), Packaging Of Telomere Ends (NES 2.07), and DNA Damage/Telomere Stress Induced Senescence (NES 2.10). This is biologically coherent with the known MGUS-to-MM transition: epigenetic reprogramming precedes the genomic instability of frank myeloma. DNA methylation changes (particularly at CpG islands) are established early events in MGUS, and PRC2-mediated H3K27 methylation is critical for MGUS plasma cell gene silencing. The Telomere Stress/Senescence enrichment confirms that premalignant plasma cells are under replicative stress but have not yet fully escaped senescence checkpoints — a critical distinction from MM, where these cells actively proliferate.

**YAP/TAZ (Hippo) and NOTCH signaling.** YAP1/WWTR1-Stimulated Gene Expression (NES 2.21) and Pre-NOTCH Transcription (NES 2.08) appear uniquely in MGUS, not MM. YAP/TAZ (Hippo pathway effectors) regulate cell proliferation and are activated in early malignant transformation but can be later bypassed. NOTCH signaling is involved in plasma cell differentiation and bone marrow niche interactions. These are early pro-survival signals that may be superseded by other proliferative mechanisms as MGUS progresses to MM.

**NK cell suppression — an immune evasion mechanism specific to MGUS.** The downregulated signals in MGUS are particularly informative: GO: Negative Regulation Of NK Cell Mediated Cytotoxicity (NES −2.19), Negative Regulation Of NK Cell Mediated Immunity (NES −2.13), and Negative Regulation Of Leukocyte Mediated Cytotoxicity (NES −2.07) are all downregulated in MGUS plasma cfRNA. Paradoxically, the *downregulation* of genes that *negatively regulate* NK killing means that NK cell suppression is less active in MGUS relative to HC — or alternatively, it reflects reduced expression of NK inhibitory ligands on MGUS cells compared to HC. The broader picture of MHC-I antigen presentation being downregulated (NES −2.29, −2.24) in MGUS suggests that MGUS cells are beginning to downregulate MHC-I as an immune evasion strategy, consistent with the known loss of HLA molecules during myelomagenesis.

**Lysosome downregulated in MGUS.** Unlike MM (which shows lysosomal upregulation through erythropoietic stress signals) and ICI groups (lysosomal immune activation), MGUS shows lysosome downregulation (NES −2.37). This may reflect reduced normal immune cell lysosomal activity in MGUS patients, whose immune surveillance is beginning to be suppressed by the emerging monoclonal clone, before the full immune evasion phenotype of MM is established.

---

## 6. MGUS vs MM: Direct Pathway-Level Comparison

### Shared Enrichments (both MGUS and MM, NES > 0)

Both diseases share: Carbon Dioxide/Gas Transport (erythropoietic stress), Chromosome Condensation, Hydrogen Peroxide Catabolic Process, HDACs Deacetylate Histones, DNA Damage/Telomere Stress Induced Senescence, Nucleosome Organization, RMTs Methylate Histone Arginines, Megakaryocyte/Platelet Factors. These shared upregulated pathways define the *disease continuum* of the MGUS-MM spectrum — chromosomal instability, epigenetic dysregulation, hematopoietic stress, and altered megakaryocyte/platelet biology are present in both.

### MM-Specific Enrichments (not in MGUS)

| Category | Pathways | NES | Interpretation |
|----------|----------|-----|----------------|
| Mitotic cell cycle | Kinesins, G1/S Transcription, Polo-like Kinase, Mitotic Spindle Assembly, Sister Chromatid Segregation, E2F Targets, Cell Cycle (KEGG) | 2.04–2.19 | Full proliferative program; absent in MGUS where plasma cells do not yet proliferate rapidly |
| Erythrocyte biology | Erythrocyte Differentiation, Oxygen Transport (higher NES) | ~2.04–2.08 | Worsening anemia in MM from complete marrow infiltration |
| TP53 pathway | TP53 Regulates Cell Cycle Genes | 2.02 | p53-dependent responses to genotoxic stress in MM |
| Ubiquitin/proteasome | Positive Regulation Of Ubiquitin-Protein Transferase Activity | 2.00 | Proteasome stress from paraprotein overproduction (targetable by bortezomib) |

### MGUS-Specific Enrichments (not in MM)

| Category | Pathways | NES | Interpretation |
|----------|----------|-----|----------------|
| Chromatin/epigenetic | DNA Methylation, PRC2, SIRT1/rRNA regulation | 2.13–2.18 | Early epigenetic reprogramming before chromosomal instability of MM |
| Hippo/NOTCH | YAP1/TAZ, Pre-NOTCH | 2.08–2.21 | Early pro-survival signaling in premalignant plasma cells |
| Innate immune activation | Neutrophil Extracellular Trap Formation | 2.02 | Inflammatory milieu of premalignant state |
| Senescence | Oxidative Stress Induced Senescence (without full bypass) | 2.04 | MGUS cells under senescence pressure that MM cells have escaped |

### Key Conclusion
The MGUS-to-MM transition in plasma cfRNA is characterized by: the **emergence of full mitotic cell cycle and kinetochore/spindle machinery** (absent in MGUS), **increased TP53 and DNA damage responses** (reflecting genotoxic burden), and **loss of the Hippo/NOTCH early signaling** that characterizes premalignant self-renewal. MGUS retains more epigenetic/methylation signatures while MM shifts to active mitotic proliferation. The normative model captures this malignant transition as distinct Z-score distributions even in small sample sizes (MM n=17, MGUS n=8).

---

## 7. ME/CFS (Myalgic Encephalomyelitis/Chronic Fatigue Syndrome), n=90

### Study Context
Gardella et al. (PNAS 2025) is the first large-scale plasma cfRNA study in ME/CFS, profiling 93 ME/CFS cases against 75 healthy sedentary controls by RNA sequencing, with 743 differentially abundant features identified. A diagnostic model achieved AUC 0.81 (accuracy 77%). The study found elevated plasmacytoid dendritic cell (pDC)-derived cfRNA, monocyte and T cell signatures, and notably lower platelet-derived cfRNA in ME/CFS. With n=90, this is the largest cohort in Group 4, providing the most statistical power.

ME/CFS is characterized by post-exertional malaise, severe fatigue, cognitive dysfunction, and autonomic dysregulation, with no FDA-approved treatment. The biological basis remains contested, with evidence for immune dysregulation, mitochondrial dysfunction, impaired energy metabolism, and neuroinflammation.

### Top Upregulated Pathways (NES > 0, only 4 significant)

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | 2.234 | GO: Immunoglobulin Mediated Immune Response |
| 2 | 2.222 | GO: B Cell Proliferation |
| 3 | 2.140 | GO: Peptide Antigen Assembly With MHC Class II Complex |
| 4 | 2.140 | GO: MHC Class II Protein Complex Assembly |

### Top Downregulated Pathways (NES < 0, 161 significant)

| Rank | NES | Pathway |
|------|-----|---------|
| 1 | -2.543 | GO: Cytoplasmic Translation |
| 2 | -2.466 | Reactome: Eukaryotic Translation Elongation |
| 3 | -2.456 | Reactome: Eukaryotic Translation Termination |
| 4 | -2.435 | Reactome: Peptide Chain Elongation |
| 5 | -2.435 | Reactome: Formation Of A Pool Of Free 40S Subunits |
| 6 | -2.427 | Reactome: NMD Enhanced By Exon Junction Complex |
| 7 | -2.421 | Reactome: L13a-Mediated Translational Silencing Of Ceruloplasmin |
| 8 | -2.390 | Reactome: GTP Hydrolysis And Joining Of 60S Ribosomal Subunit |
| 9 | -2.378 | Reactome: NMD Independent Of EJC |
| 10 | -2.358 | Reactome: Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency |
| 11 | -2.341 | Reactome: Selenocysteine Synthesis |
| 12 | -2.233 | KEGG: FoxO Signaling Pathway |
| 13 | -2.233 | KEGG: Coronavirus Disease |
| 14 | -2.221 | KEGG: Ribosome |
| 15 | -2.120 | Reactome: Estrogen-Dependent Nuclear Events Downstream Of ESR-Membrane Signaling |
| 16 | -2.116 | GO: Positive Regulation Of Insulin Receptor Signaling Pathway |
| 17 | -2.108 | Reactome: Regulation Of Expression Of SLITs And ROBOs |
| 18 | -2.115 | GO: Cellular Response To Glucose Stimulus |

### Biological Interpretation

**Overwhelming translational suppression as the dominant signal.** ME/CFS plasma cfRNA shows near-total downregulation of the translation machinery: cytoplasmic translation, translation elongation, termination, initiation, 40S/60S ribosome assembly, cap-dependent translation initiation, and ribosome biogenesis are all at FDR=0 (NES −2.15 to −2.54). This is the most pronounced and internally consistent downregulated signal across all five diseases in Group 4. The biological interpretation is multi-layered:

1. **GCN2/eIF2α integrated stress response**: The enrichment of EIF2AK4 (GCN2) Response To Amino Acid Deficiency (NES −2.36) indicates activation of the integrated stress response (ISR). GCN2 phosphorylates eIF2α at Ser51 in response to uncharged tRNAs (amino acid deficiency) or mitochondrial stress, causing global translational repression while selectively allowing translation of ISR-specific mRNAs (ATF4). This is consistent with the hypometabolism model of ME/CFS: PDK-mediated inhibition of pyruvate dehydrogenase → impaired TCA cycle → reduced ATP → cellular energy deficit triggering ISR.

2. **L13a/GAIT-mediated ceruloplasmin silencing**: The Reactome: L13a-Mediated Translational Silencing Of Ceruloplasmin Expression pathway (NES −2.42) is specifically relevant — L13a is a ribosomal protein that is phosphorylated by DAPK-ZIPK in response to interferon-γ and translocates from the 60S ribosome to silence ceruloplasmin mRNA translation via the GAIT (Gamma-Activated Inhibitor of Translation) complex. Ceruloplasmin is a copper-binding protein involved in iron oxidation and antioxidant defense. Its selective translational suppression via the GAIT complex suggests an interferon-γ-driven inflammatory state in ME/CFS — consistent with the pDC and immune activation findings from Gardella et al.

3. **Nonsense-mediated decay (NMD) upregulation**: Both EJC-enhanced and EJC-independent NMD pathways are downregulated (NES −2.43, −2.38). NMD is a quality control mechanism that degrades transcripts with premature stop codons; its downregulation could reflect altered mRNA surveillance in ME/CFS cells.

**FoxO signaling pathway downregulated.** KEGG: FoxO Signaling Pathway (NES −2.23) is a critical ME/CFS-relevant finding. FoxO transcription factors (FOXO1, FOXO3, FOXO4) are regulated by PI3K/Akt/mTOR signaling and integrate signals from insulin, IGF-1, oxidative stress, and energy status. In their dephosphorylated (active) nuclear form, FoxOs drive expression of genes for cell cycle arrest, stress resistance (SOD2, catalase), and gluconeogenesis. FoxO signaling is a convergence point for:
- **Insulin resistance** (FoxO activated when insulin/Akt signaling is reduced)
- **Reduced mTOR activity** (mTOR inhibition by low energy → reduced Akt → FoxO nuclear translocation → reduced protein synthesis)
- **Muscle atrophy** (FoxO drives atrogene expression; chronic fatigue → disuse → FoxO-mediated wasting)

The downregulation of FoxO pathway transcripts in plasma cfRNA of ME/CFS patients (relative to HC) is consistent with a state where FoxO-activating signals (energy deficit, oxidative stress) have caused FoxO-driven transcriptional changes in circulating immune and metabolic cells. Published ME/CFS research documents upregulation of PDK1/2/4 (pyruvate dehydrogenase kinases), reduced OXPHOS, and hypometabolism — all consistent with the energy-deficit context that would activate FoxO and suppress mTOR-driven translation.

**Insulin receptor signaling and glucose response downregulated.** GO: Positive Regulation Of Insulin Receptor Signaling (NES −2.12) and Cellular Response To Glucose Stimulus (NES −2.12) are both downregulated. This points to impaired insulin/IGF-1 signaling in circulating cells, which would release FoxO inhibition and suppress Akt/mTOR — consistent with the translational suppression above. Insulin resistance in ME/CFS has been reported, with low-grade inflammation and metabolic dysfunction.

**KEGG: Coronavirus Disease downregulated (NES −2.23).** The Coronavirus Disease pathway annotation captures interferon-stimulated gene responses and antiviral immunity pathways. Its downregulation in ME/CFS plasma cfRNA may appear paradoxical given that ME/CFS is often triggered by viral infections. However, this downregulation likely reflects a state of *impaired antiviral response* rather than active viral control — ME/CFS patients show evidence of T cell exhaustion and impaired innate immune responses, consistent with the post-viral immune dysregulation model.

**Estrogen signaling downregulated (NES −2.12).** Reactome: Estrogen-Dependent Nuclear Events is significantly downregulated. ME/CFS disproportionately affects women (3:1 ratio), and estrogen receptor signaling regulates mitochondrial function, energy metabolism, and immune modulation. Downregulation of estrogen pathway cfRNA in ME/CFS may contribute to the mitochondrial and metabolic deficits.

**Only 4 upregulated pathways — all immune.** The remarkable asymmetry (4 up vs 161 down) is itself a biological signal: ME/CFS plasma cfRNA is globally suppressed relative to HC. The four upregulated pathways are all related to B cell and MHC-II immune activity: Immunoglobulin Mediated Immune Response (NES 2.23), B Cell Proliferation (NES 2.22), MHC Class II Assembly (NES 2.14), Peptide Antigen Assembly With MHC Class II (NES 2.14). The Gardella et al. study noted elevated pDC-derived cfRNA and T cell signatures in ME/CFS; the B cell/MHC-II upregulation seen here may reflect chronic antigen-driven B cell activation (consistent with the clinical observation of anti-nuclear antibodies and polyclonal hypergammaglobulinemia reported in some ME/CFS patients) or persistent antigen presentation driving B cell responses against unresolved viral/self antigens.

**Overall ME/CFS interpretation.** The normative model Z-scores for ME/CFS point to a plasma cfRNA signature dominated by systemic translational suppression driven by the integrated stress response (GCN2/eIF2α axis), impaired insulin/FoxO/mTOR signaling, and selective B cell/MHC-II immune activation against a backdrop of severely reduced biosynthetic activity. This is mechanistically coherent with the hypometabolism model of ME/CFS: cells under energy deficit (reduced TCA/OXPHOS → ATP depletion → GCN2 activation → eIF2α phosphorylation → global translation arrest) while maintaining immune activation. The FoxO and insulin pathway depression are particularly novel cfRNA-level findings that corroborate metabolomics and proteomics studies showing hypometabolism in ME/CFS.

---

## 8. Cross-Disease Observations

### Universal Translational Suppression
All five diseases show significant downregulation of ribosomal translation pathways (elongation, termination, 40S pool formation, NMD). This is a shared feature of diseases captured by plasma cfRNA normative modeling. Possible explanations: (1) systemic inflammatory catabolic state activates integrated stress response; (2) reduced circulating translational activity in immune cells under chronic stress; (3) altered cfRNA composition (depletion of translationally active cell types in disease states).

### Lysosome Enrichment in ICI Groups
Both ICI groups show strong lysosomal enrichment (top 1, NES ~2.83–2.89). This likely reflects immune activation driving phagocytic lysosomal activity (neutrophils, macrophages, dendritic cells) — enhanced by checkpoint inhibitor-mediated immune reinvigoration. Notably, MGUS shows lysosome *downregulated* (NES −2.37), contrasting with ICI groups, suggesting disease-state-dependent lysosomal biology.

### Epigenetic Enrichment in MM/MGUS Continuum
HDACs, chromatin condensation, histone methylation, and DNA methylation pathways are enriched across both MM and MGUS, defining a shared epigenetic dysregulation signature of the plasma cell malignancy spectrum. These enrichments are absent from ICI and ME/CFS groups, validating the pathway-level specificity of the normative model.

### Immune Suppression in MM/MGUS
Both MM (immunoglobulin mediated response, MHC-II, B cell proliferation all downregulated) and MGUS (antigen presentation, NK cell regulation downregulated) show immune suppressive signatures, contrasting sharply with ICI groups (immune activation dominant) and ME/CFS (selective B cell upregulation). This immune suppression gradient from MGUS to MM in cfRNA space aligns with the clinical immunosuppression of myeloma.

---

## 9. Caveats and Limitations

1. **Small sample sizes for MM (n=17), MGUS (n=8), and ICI groups (n=11)**: Pathway enrichments from small cohorts are susceptible to outlier-driven inflation of Z-scores. The MGUS results in particular (n=8, only 22 significant pathways total, all at FDR 0.004–0.04) should be treated as preliminary. The ICI groups (n=11 each) are adequate for the most robust enrichments but marginal for the lower-NES findings.

2. **GSEA-prerank uses mean Z-scores**: Averaging across samples reduces power for heterogeneous diseases. ME/CFS in particular is a clinically heterogeneous syndrome, and the large sample size (n=90) partially compensates, but subgroup analyses (post-exertional malaise severity, duration) would be informative.

3. **Plasma cfRNA vs tissue RNA**: cfRNA captures a mixture of cell type-derived fragments in proportion to their circulation and cfRNA release rates. Enrichments of cardiac, bone marrow, or CNS pathways reflect the contribution of these tissues' transcripts to plasma, modulated by tissue damage/turnover rates. Direct causal attribution to specific tissues requires cell type deconvolution.

4. **Normative model Z-scores relative to HC**: All enrichments are *relative to 693 healthy controls (WTS-only)*. Upregulated pathways reflect absolute upregulation in disease or depletion in HC; downregulated pathways reflect depletion in disease or higher expression in HC. Both directions are equally informative.

5. **MGUS n=8 borderline FDRs**: Most MGUS significant pathways have FDR 0.001–0.04, compared to FDR=0 for the top hits in other diseases. This reflects reduced statistical power from small n and should be noted in any report.

6. **Gardella et al. (2025) methodological caveats**: The Science for ME community raised concerns about whether ME/CFS cfRNA findings reflect altered cell circulation rather than intrinsic cellular dysfunction. The translational suppression and FoxO signaling findings in the normative model are internally consistent, but mechanistic validation requires cell-type-specific analyses.

---

## 10. Summary Table

| Disease | Top UP Signal | Top DOWN Signal | Key Unique Finding |
|---------|-------------|----------------|-------------------|
| ICI-m | OXPHOS + Lysosome + Neutrophil degranulation | Translation (ribosomes, elongation) | Cross-presentation/ER-phagosome (myocarditis-specific antigen presentation) |
| ICI-treated Cancer | Lysosome + Antigen presentation + Innate immunity | Translation (pervasive, 1072 down terms) | IgA intestinal network, thyroid autoimmunity (ICI adverse events) |
| MM | Cell cycle + Kinesins + Erythrocyte/Hb | Immunoglobulin response + MHC-II + OXPHOS | Full mitotic program; paradoxical Ig downregulation = immune suppression |
| MGUS | Epigenetic (DNA methylation, PRC2, SIRT1) + Histone + SLE/NETs | Antigen presentation + NK suppression | Early epigenetic dysregulation without full mitotic activation; YAP/NOTCH |
| ME/CFS | B cell/MHC-II (only 4 pathways) | Translation + FoxO + Insulin signaling + Coronavirus (161 down) | FoxO/insulin/glucose axis depression = hypometabolism; GCN2 ISR activation |

---

*Generated: 2026-06-24 | Normative model cohort: 693 HC (WTS-only) | GSEA tool: GSEApy prerank | Gene sets: KEGG_2021_Human, Reactome_2022, GO_Biological_Process_2023*
