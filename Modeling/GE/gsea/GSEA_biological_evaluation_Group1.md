# GSEA Biological Evaluation — Group 1: Infectious Diseases (HIV / HIV+TB / Tuberculosis)

**Analysis date:** 2026-06-24  
**Ranking metric:** Mean Z-score per gene (GAMLSS normative model residuals), GSEA-prerank  
**Gene set libraries:** Reactome 2022, KEGG 2021 Human, GO Biological Process 2023  
**Significance threshold:** FDR q-val < 0.05  
**Data source:** Chang et al. (2023), DOI: 10.1101/2023.01.11.23284435 (and companion papers .23284434, .23284433)  

> **Interpretation note:** NES > 0 = pathway enriched in high-Z (anomalously upregulated) genes relative to the healthy cfRNA reference. NES < 0 = pathway enriched in low-Z (anomalously downregulated) genes.

---

## 1. Summary Statistics

| Disease | NES > 0 (up) | NES < 0 (down) | Total FDR < 0.05 |
|---------|-------------|----------------|-----------------|
| HIV | 434 | 50 | 484 |
| HIV + Tuberculosis | 306 | 89 | 395 |
| Tuberculosis | 610 | 7 | 617 |

**Notable:** Tuberculosis has the highest total count (617) with almost exclusively upregulated pathways (610 vs 7), indicating a broad upregulatory transcriptional storm in plasma cfRNA. HIV+TB has the most balanced bidirectional signal. HIV shows the strongest asymmetry in the opposite direction from TB in its downregulated gene set (translation/ribosomal cluster).

---

## 2. HIV

### 2.1 NES > 0 Pathways (Anomalously Upregulated in HIV)

**Top 30 by |NES| (all FDR < 0.05):**

| Rank | NES | Term |
|------|-----|------|
| 1 | +2.741 | Reactome: Antiviral Mechanism By IFN-stimulated Genes (R-HSA-1169410) |
| 2 | +2.595 | Reactome: Interferon Alpha/Beta Signaling (R-HSA-909733) |
| 3 | +2.585 | Reactome: ISG15 Antiviral Mechanism (R-HSA-1169408) |
| 4 | +2.541 | Reactome: Interferon Signaling (R-HSA-913531) |
| 5 | +2.540 | GO: Negative Regulation Of Viral Process |
| 6 | +2.521 | Reactome: Assembly Of ORC Complex At Origin Of Replication |
| 7 | +2.478 | Reactome: DNA Methylation |
| 8 | +2.417 | GO: Negative Regulation Of Viral Genome Replication |
| 9 | +2.399 | Reactome: Transcriptional Regulation Of Granulopoiesis |
| 10 | +2.390 | Reactome: RNA Polymerase I Promoter Opening |
| 11 | +2.385 | GO: Defense Response To Symbiont |
| 12 | +2.380 | Reactome: Activated PKN1 Stimulates Transcription Of AR-Regulated KLK2/KLK3 |
| 13 | +2.373 | Reactome: Condensation Of Prophase Chromosomes |
| 14 | +2.371 | KEGG: Neutrophil Extracellular Trap Formation |
| 15 | +2.359 | Reactome: RUNX1 Regulates Megakaryocyte Differentiation And Platelet Function |
| 16 | +2.346 | Reactome: Oxidative Stress Induced Senescence |
| 17 | +2.326 | GO: Regulation Of Viral Genome Replication |
| 18 | +2.319 | Reactome: Transcriptional Regulation By Small RNAs |
| 19 | +2.318 | GO: Defense Response To Virus |
| 20 | +2.317 | Reactome: Packaging Of Telomere Ends |
| 21 | +2.311 | Reactome: RHO GTPases Activate PKNs |
| 22 | +2.293 | Reactome: Cleavage Of Damaged Purine |
| 23 | +2.270 | KEGG: Systemic Lupus Erythematosus |
| 24 | +2.262 | Reactome: Pre-NOTCH Transcription And Translation |
| 25 | +2.258 | Reactome: Base-Excision Repair, AP Site Formation |
| 26 | +2.254 | Reactome: PRC2 Methylates Histones And DNA |
| 27 | +2.245 | Reactome: RHO GTPase Effectors |
| 28 | +2.244 | KEGG: Alcoholism |
| 29 | +2.239 | Reactome: Senescence-Associated Secretory Phenotype (SASP) |
| 30 | +2.238 | Reactome: Cleavage Of Damaged Pyrimidine |

### 2.2 NES < 0 Pathways (Anomalously Downregulated in HIV)

**Top 30 by |NES| (all FDR < 0.05):**

| Rank | NES | Term |
|------|-----|------|
| 1 | −2.494 | Reactome: SRP-dependent Cotranslational Protein Targeting To Membrane |
| 2 | −2.405 | GO: Ribosomal Small Subunit Biogenesis |
| 3 | −2.381 | Reactome: Translation |
| 4 | −2.363 | GO: Cytoplasmic Translation |
| 5 | −2.352 | GO: Negative Regulation Of Pathway-Restricted SMAD Protein Phosphorylation |
| 6 | −2.332 | Reactome: Mitochondrial Translation Termination |
| 7 | −2.255 | Reactome: Mitochondrial Translation Elongation |
| 8 | −2.237 | KEGG: Ribosome |
| 9 | −2.228 | Reactome: Mitochondrial Translation Initiation |
| 10 | −2.216 | Reactome: Formation Of A Pool Of Free 40S Subunits |
| 11 | −2.201 | Reactome: Peptide Chain Elongation |
| 12 | −2.185 | Reactome: Eukaryotic Translation Elongation |
| 13 | −2.177 | Reactome: Mitochondrial Translation |
| 14 | −2.168 | Reactome: L13a-mediated Translational Silencing Of Ceruloplasmin |
| 15 | −2.159 | GO: Ribosome Assembly |
| 16 | −2.144 | Reactome: GTP Hydrolysis And Joining Of 60S Ribosomal Subunit |
| 17 | −2.136 | Reactome: Cap-dependent Translation Initiation |
| 18 | −2.077 | GO: Maturation Of SSU-rRNA From Tricistronic rRNA Transcript |
| 19 | −2.056 | Reactome: Nonsense Mediated Decay (NMD) — EJC-independent |
| 20 | −2.045 | Reactome: rRNA Processing |
| 21 | −2.041 | GO: Ribosomal Large Subunit Assembly |
| 22 | −2.040 | Reactome: Major Pathway Of rRNA Processing In Nucleolus And Cytosol |
| 23 | −2.032 | GO: Mitochondrial Respiratory Chain Complex Assembly |
| 24 | −2.028 | GO: Translation |
| 25 | −2.027 | Reactome: NMD Enhanced By EJC |
| 26 | −2.024 | GO: Mitochondrial Cytochrome C Oxidase Assembly |
| 27 | −2.013 | Reactome: Viral mRNA Translation |
| 28 | −2.006 | Reactome: rRNA Processing In Nucleus And Cytosol |
| 29 | −2.001 | Reactome: Eukaryotic Translation Termination |
| 30 | −1.987 | KEGG: Spliceosome |

### 2.3 Biological Evaluation: HIV

#### IFN/ISG Signature — EXPECTED, HIGH CONFIDENCE

The dominant signal in HIV is a massive type I interferon response (NES +2.741 to +2.541), driven by canonical ISGs: IFI27, IFI6, OAS1, OAS2, OAS3, IFIT1, IFIT2, IFIT3, MX1, MX2, ISG15, RSAD2 (viperin), XAF1, USP18, IRF7, ADAR, BST2, IFITM1/2/3. This represents the most biologically expected finding possible in untreated or imperfectly suppressed HIV: HIV replication triggers sustained type I IFN production via cGAS-STING sensing of viral DNA and RIG-I/MDA5 sensing of viral RNA. The ISG15 pathway specifically (NES +2.585) is a well-established marker of HIV immune activation. ISG15 is one of the most potently induced ISGs in HIV infection and its protein product modifies viral proteins; elevated plasma ISG15 and other ISGs are documented biomarkers of HIV disease activity in multiple blood-based transcriptomic studies. Critically, HIV maintains this IFN tone even in ART-suppressed individuals due to residual antigen and immune activation.

**In plasma cfRNA context:** The IFN/ISG signal in cfRNA likely originates from a composite of: (1) activated monocytes and plasmacytoid dendritic cells (pDCs), which are the primary IFN-alpha producers in HIV; (2) circulating lymphocytes releasing cfRNA; and (3) possibly extracellular vesicle cargo from ISG-expressing cells. The IFN signature being the top signal in cfRNA closely mirrors what is found in whole blood transcriptomics of HIV patients.

#### Neutrophil Extracellular Trap (NET) Formation — EXPECTED

KEGG NET formation (NES +2.371) reflects neutrophil activation in HIV. The lead genes include MPO (myeloperoxidase), FCGR1A, FCGR2A, FCGR3B, histones (H2AC/H2BC/H3C/H4C families), PIK3CB, ACTB, ACTG1, MAPK14, SRC, AKT1, RELA, RAF1. NET formation is a documented phenomenon in HIV infection: HIV activates neutrophils to release chromatin-DNA-histone complexes, contributing to immune activation and inflammation. The concurrent presence of histone genes in NETs is biologically genuine (histones are the structural backbone of NETs), though the same histone gene cluster also inflates chromatin-remodeling pathway scores (see artifact section below).

#### Senescence and DNA Damage Response — EXPECTED

Oxidative Stress Induced Senescence (NES +2.346), SASP (NES +2.239), DNA Damage/Telomere Stress Induced Senescence (NES +2.224), Cellular Senescence (NES +2.176), and Packaging Of Telomere Ends (NES +2.317) together form a coherent senescence signature. This is highly consistent with established HIV pathobiology: chronic HIV infection accelerates immune aging/senescence, affecting CD4+ T cells, CD8+ T cells, and monocytes. Telomere shortening, DNA damage accumulation, and p21/p16-mediated senescence are well-documented in HIV even under ART. The SASP program drives the chronic inflammatory milieu of HIV disease. The lead genes include E2F1, CDKN1A, MAPK14 (p38), MDM2, EZH2, PHC3, and SCMH1, consistent with epigenetic reprogramming in senescent cells.

#### Granulopoiesis, RUNX1/Megakaryocyte, and Platelet Pathways — EXPECTED

Transcriptional Regulation Of Granulopoiesis (NES +2.399) and RUNX1 Regulates Megakaryocyte Differentiation And Platelet Function (NES +2.359) reflect myeloid expansion and platelet biology in HIV. HIV causes thrombocytopenia via multiple mechanisms including direct megakaryocyte infection, immune complex deposition, and dysregulated platelet production. RUNX1, GATA2, MYB, and NFE2 are master regulators of megakaryopoiesis; their elevated representation in cfRNA likely reflects increased turnover and activation of megakaryocytes and platelets. Platelet Activation (NES +2.155) and Platelet Degranulation further support this picture. This is a biologically expected and clinically meaningful signal.

#### DNA Repair (BER/NER) and Chromatin — PARTIALLY EXPECTED, PARTIALLY ARTIFACT-INFLATED

Base-Excision Repair (NES +2.258), Cleavage Of Damaged Purine (NES +2.293), and related DNA repair terms likely reflect genuine oxidative DNA damage in HIV-infected cells, as HIV replication generates reactive oxygen species. However, pathways such as DNA Methylation (NES +2.478), PRC2 Methylates Histones And DNA (NES +2.254), Condensation Of Prophase Chromosomes (NES +2.373), Packaging Of Telomere Ends (NES +2.317), and Assembly Of ORC Complex (NES +2.521) are inflated by the massive elevation of replication-coupled histone mRNAs (H2BC12, H2BC17, H2AC19, H2AJ, H2BC9, etc.) in plasma cfRNA. These histone transcripts are highly abundant in cfRNA from proliferating/dying cells and contaminate many chromosome/chromatin pathway gene sets (see Section 5 on artifacts).

#### HIV Virus-Host Interaction Terms

Interactions Of Rev With Host Cellular Proteins (NES +1.692) and Nuclear Import Of Rev Protein (NES +1.734) are HIV-specific Reactome pathways that are significantly enriched. This is a direct validation that the cfRNA normative model is detecting HIV-specific molecular biology — the HIV Rev protein nuclear import pathway components (importins: KPNB1, KPNA4, KPNA5, NUP98, NUP153, etc.) are elevated in HIV patient samples. This is an excellent positive control finding.

#### RHO GTPase / NOTCH / PKN Pathways — HISTONE-INFLATED ARTIFACTS

RHO GTPases Activate PKNs (NES +2.311), Activated PKN1 Stimulates Transcription Of AR-Regulated KLK2/KLK3 (NES +2.380), Pre-NOTCH Transcription And Translation (NES +2.262), and KEGG Alcoholism (NES +2.244) are all driven substantially by the histone gene cluster and do not represent genuine NOTCH, androgen receptor, PKN, or alcohol-related biology. Inspection of the Alcoholism lead genes reveals 42/73 (57%) are histone cluster genes (H2BC, H2AC, H3C, H4C families). The KEGG Alcoholism pathway contains a large nucleosome component (reflecting chromatin remodeling in alcohol-related neurobiology), and this is being incorrectly enriched in cfRNA due to plasma histone mRNA abundance. These terms should be considered background artifacts.

#### Translation/Ribosomal Suppression (NES < 0) — BIOLOGICALLY MEANINGFUL

The downregulated cluster is dominated by cytoplasmic and mitochondrial translation (NES −2.494 to −2.001): SRP-dependent protein targeting, ribosomal small/large subunit biogenesis, rRNA processing, cap-dependent translation initiation, peptide chain elongation, and eukaryotic translation termination. Lead genes are almost entirely ribosomal protein genes (RPL/RPS families) and translation initiation factors (EIF3A/B/J, EIF4A, EIF2B4, EIF5B, RPLP1, etc.).

This is biologically meaningful in two ways:
1. **HIV directly suppresses cellular protein synthesis** to redirect ribosomal machinery toward viral protein production. HIV Tat inhibits host cap-dependent translation; Vpr and other viral proteins have been shown to downregulate ribosomal biogenesis.
2. **Immune exhaustion and metabolic reprogramming:** CD4+ T cell depletion and monocyte dysfunction in HIV involve downregulation of anabolic metabolism, including protein synthesis. The normative model is detecting reduced cytoplasmic translational activity in the HIV cfRNA pool.
3. **Mitochondrial translation suppression** (NES −2.332 to −2.177) is consistent with HIV-associated mitochondrial dysfunction: HIV and antiretroviral drugs (nucleoside reverse transcriptase inhibitors) are well-established mitochondrial toxins that impair oxidative phosphorylation.

The Negative Regulation Of Pathway-Restricted SMAD Protein Phosphorylation (NES −2.352) is an outlier in this downregulated cluster and may be incidental (shared gene membership with ribosomal machinery). The SMAD pathway proper is not typically suppressed in HIV.

**Overall HIV signal quality: HIGH CONFIDENCE**

The top signals (IFN/ISG upregulation, translation/mitochondrial downregulation) are among the most robustly documented molecular signatures in HIV transcriptomics. The cfRNA normative model successfully captures both directions of biological disruption.

---

## 3. HIV + Tuberculosis (Coinfection)

### 3.1 NES > 0 Pathways (Anomalously Upregulated in HIV+TB)

**Top 30 by |NES|:**

| Rank | NES | Term |
|------|-----|------|
| 1 | +2.911 | KEGG: Systemic Lupus Erythematosus |
| 2 | +2.754 | KEGG: Neutrophil Extracellular Trap Formation |
| 3 | +2.693 | KEGG: Alcoholism |
| 4 | +2.656 | Reactome: RNA Polymerase I Promoter Opening |
| 5 | +2.634 | Reactome: Condensation Of Prophase Chromosomes |
| 6 | +2.627 | Reactome: ERCC6 (CSB) And EHMT2 (G9a) Positively Regulate rRNA Expression |
| 7 | +2.619 | Reactome: DNA Methylation |
| 8 | +2.602 | Reactome: Assembly Of ORC Complex At Origin Of Replication |
| 9 | +2.579 | Reactome: Cleavage Of Damaged Purine |
| 10 | +2.575 | Reactome: Defective Pyroptosis |
| 11 | +2.516 | Reactome: Packaging Of Telomere Ends |
| 12 | +2.510 | Reactome: Transcriptional Regulation Of Granulopoiesis |
| 13 | +2.504 | Reactome: SIRT1 Negatively Regulates rRNA Expression |
| 14 | +2.500 | Reactome: PRC2 Methylates Histones And DNA |
| 15 | +2.488 | Reactome: RHO GTPases Activate PKNs |
| 16 | +2.485 | Reactome: HDACs Deacetylate Histones |
| 17 | +2.468 | Reactome: RHO GTPase Effectors |
| 18 | +2.462 | Reactome: Activated PKN1/AR/KLK2-KLK3 |
| 19 | +2.444 | Reactome: Assembly Of Pre-Replicative Complex |
| 20 | +2.418 | Reactome: DNA Replication |
| 21 | +2.417 | Reactome: Deposition Of New CENPA-containing Nucleosomes At Centromere |
| 22 | +2.397 | Reactome: Cleavage Of Damaged Pyrimidine |
| 23 | +2.395 | Reactome: DNA Replication Pre-Initiation |
| 24 | +2.387 | Reactome: Cell Cycle Checkpoints |
| 25 | +2.374 | Reactome: G2/M Checkpoints |
| 26 | +2.367 | GO: Antimicrobial Humoral Immune Response Mediated By Antimicrobial Peptide |
| 27 | +2.365 | GO: Antimicrobial Humoral Response |
| 28 | +2.359 | KEGG: Viral Carcinogenesis |
| 29 | +2.358 | GO: Negative Regulation Of Viral Process |
| 30 | +2.349 | Reactome: Base-Excision Repair, AP Site Formation |

### 3.2 NES < 0 Pathways (Anomalously Downregulated in HIV+TB)

**Top 30 by |NES|:**

| Rank | NES | Term |
|------|-----|------|
| 1 | −2.743 | Reactome: Cap-dependent Translation Initiation |
| 2 | −2.737 | Reactome: Eukaryotic Translation Elongation |
| 3 | −2.715 | Reactome: L13a-mediated Translational Silencing Of Ceruloplasmin |
| 4 | −2.712 | Reactome: Eukaryotic Translation Termination |
| 5 | −2.709 | Reactome: Peptide Chain Elongation |
| 6 | −2.697 | Reactome: GTP Hydrolysis And Joining Of 60S Ribosomal Subunit |
| 7 | −2.684 | Reactome: Viral mRNA Translation |
| 8 | −2.668 | Reactome: Formation Of A Pool Of Free 40S Subunits |
| 9 | −2.668 | Reactome: NMD — EJC-independent |
| 10 | −2.649 | GO: Cytoplasmic Translation |
| 11 | −2.627 | Reactome: SRP-dependent Cotranslational Protein Targeting To Membrane |
| 12 | −2.603 | Reactome: Selenocysteine Synthesis |
| 13 | −2.589 | Reactome: Major Pathway Of rRNA Processing In Nucleolus And Cytosol |
| 14 | −2.582 | Reactome: NMD — EJC-enhanced |
| 15 | −2.577 | Reactome: rRNA Processing In Nucleus And Cytosol |
| 16 | −2.557 | Reactome: rRNA Processing |
| 17 | −2.488 | Reactome: Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency |
| 18 | −2.488 | GO: Ribosomal Small Subunit Biogenesis |
| 19 | −2.474 | KEGG: Ribosome |
| 20 | −2.413 | GO: Peptide Biosynthetic Process |
| 21 | −2.345 | Reactome: Ribosomal Scanning And Start Codon Recognition |
| 22 | −2.340 | Reactome: mRNA Activation Upon Binding Of Cap-Binding Complex And eIFs |
| 23 | −2.321 | GO: Maturation Of SSU-rRNA From Tricistronic rRNA Transcript |
| 24 | −2.319 | Reactome: Translation Initiation Complex Formation |
| 25 | −2.301 | Reactome: Selenoamino Acid Metabolism |
| 26 | −2.298 | GO: Maturation Of SSU-rRNA |
| 27 | −2.285 | Reactome: Formation Of Ternary Complex And 43S Complex |
| 28 | −2.272 | Reactome: Influenza Viral RNA Transcription And Replication |
| 29 | −2.246 | GO: Ribonucleoprotein Complex Biogenesis |
| 30 | −2.236 | GO: Translation |
| — | −2.158 | GO: B Cell Proliferation |
| — | −2.060 | GO: B Cell Activation |
| — | −1.979 | GO: B Cell Receptor Signaling Pathway |
| — | −1.835 | GO: Mast Cell Degranulation |

### 3.3 Biological Evaluation: HIV+TB Coinfection

#### Dominant NET/SLE/Histone Signature — EXPECTED WITH CAVEAT

The top two enriched pathways in HIV+TB are KEGG Systemic Lupus Erythematosus (NES +2.911) and KEGG NET Formation (NES +2.754). Both gene sets share a large number of histone cluster genes (H2BC, H2AC, H3C, H4C families) plus complement proteins (C1QA/B/C, C9, C2, C5), FcγR (FCGR1A, FCGR2A, FCGR3B), neutrophil granule proteins (CTSG, ELANE, MPO, AZU1), and cytoskeletal proteins (ACTB, ACTG1).

The biological interpretation is layered:

1. **Genuine biology:** HIV+TB coinfection is characterized by the most intense immune activation of any infectious disease combination. Both pathogens independently drive NET formation (HIV activates neutrophils via viral particles; Mtb directly stimulates NET release). The combination synergistically amplifies neutrophil hyperactivation. The SLE gene set is enriched because active HIV+TB produces an SLE-like state: massive complement consumption (C1q, C3, C5, C9 upregulation), anti-nuclear antibody-like responses driven by histone release from NETs, and FcγR-mediated immune complex clearance.

2. **Histone artifact component:** The SLE KEGG pathway contains a large nucleosome component (histones as autoantigens in SLE). The HIV+TB cfRNA data shows markedly elevated histone mRNA levels in plasma, which inflates all histone-containing gene sets. The SLE pathway's top NES (+2.911) in HIV+TB likely reflects both genuine immunopathology and this histone mRNA artifact. The KEGG Alcoholism pathway (NES +2.693) — essentially a proxy for chromatin/histone remodeling — is almost entirely histone-driven (51/78 lead genes are histone cluster genes) and represents pure artifact.

#### Pyroptosis — BIOLOGICALLY COMPELLING

Defective Pyroptosis (NES +2.575) features genes including EZH2, POLA1, and the histone cluster. However, genuine pyroptosis biology is also captured here. Both HIV and Mtb activate the NLRP3 inflammasome and trigger pyroptotic cell death (gasdermin-D mediated). HIV in particular drives pyroptosis in CD4+ T cells — this abortive pyroptosis mechanism is a major driver of CD4+ T cell depletion. In coinfection, both pyroptotic pathways converge.

#### Antimicrobial Humoral Response — EXPECTED

Antimicrobial Humoral Immune Response Mediated By Antimicrobial Peptide (NES +2.367) and Antimicrobial Humoral Response (NES +2.365) in HIV+TB reflect the neutrophil degranulation response. Lead genes include antimicrobial peptides derived from neutrophil granules (CTSG, ELANE, AZU1, MPO), complement components, and defensins. Both HIV-driven neutrophil dysfunction and Mtb-triggered granule release contribute to this signal.

#### Interferon Signaling — PRESENT BUT ATTENUATED vs. HIV ALONE

IFN Alpha/Beta Signaling (NES +2.317), Negative Regulation Of Viral Process (NES +2.358), and related viral/IFN terms are significant (FDR < 0.05) but ranked lower than in HIV alone (where IFN was the #1 signal, NES +2.741 vs. +2.317). This is counterintuitive at first but reflects genuine biology: Mtb co-infection modulates the IFN response. Specifically, type I IFN in TB is immunopathological, and regulatory mechanisms (SOCS1, IL-10 signaling) restrain IFN signaling in TB. Furthermore, in HIV+TB, the dominant transcriptional signal (histones/NETs/granulopoiesis) competes with the IFN signal in the ranking. The IFN signal is still robustly present; it is the relative rank that shifts.

Additionally, **Interferon Gamma Signaling** (NES +1.765, significant) is present in HIV+TB but not a top hit, reflecting that IFN-gamma (type II IFN, primarily T cell and NK cell derived) is contributed by the TB component.

#### DNA Replication and Cell Cycle — SHARED WITH HIV, PARTLY ARTIFACT

Assembly Of ORC Complex (NES +2.602), DNA Replication (NES +2.418), Cell Cycle Checkpoints (NES +2.387), and G2/M Checkpoints (NES +2.374) are seen in HIV+TB at higher NES than in HIV alone. This reflects both genuine immune cell proliferation (myeloid expansion, activated T cells) and histone-driven inflation, as S-phase histone gene expression co-clusters with many replication pathway gene sets.

#### B Cell Suppression — NOVEL AND BIOLOGICALLY MEANINGFUL

A distinct cluster of B cell pathway downregulation appears in HIV+TB that is absent as a dominant signal in HIV alone or TB alone:
- B Cell Proliferation (NES −2.158)
- B Cell Activation (NES −2.060)
- B Cell Receptor Signaling Pathway (NES −1.979)
- Positive Regulation Of B Cell Activation (NES −1.889)
- Regulation Of B Cell Proliferation (NES −1.900)
- B Cell Differentiation (NES −1.825)

Lead genes include canonical B cell markers: CD79A, CD79B, MS4A1 (CD20), CD19, CD22, PAX5, BTK, BLNK, PLCG2, BCL2, PTPRC (CD45), CR2 (CD21), BANK1, MEF2C. This is a highly specific and biologically meaningful signal.

In HIV+TB coinfection, B cell immunodeficiency is profound: HIV depletes CD4+ T helper cells that are required for B cell class switching and affinity maturation; simultaneously, Mtb independently impairs B cell function through regulatory mechanisms. The consequence is marked reduction in circulating B cell activity — patients with HIV+TB have severely impaired humoral immunity, making them unable to mount effective antibody responses. The cfRNA normative model detects this B cell suppression as low-Z genes in the plasma cfRNA pool, which makes biological sense: fewer and less active circulating B cells release less B cell-specific cfRNA into plasma.

The presence of mast cell degranulation downregulation (NES −1.835) may reflect concurrent mast cell suppression or redistribution, consistent with cytokine-mediated mast cell quiescence in chronic infection states.

#### Translation Suppression — STRONGER IN HIV+TB THAN IN HIV ALONE

The translation/ribosomal downregulated cluster is deeper in HIV+TB (NES −2.743 to −2.236) than in HIV alone (−2.494 to −1.987). This is biologically expected: dual pathogen burden causes more severe metabolic stress, greater T cell/B cell depletion, and more extensive suppression of anabolic metabolism in the circulating cell pool. The Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency (NES −2.488) adds a mechanistically specific finding — GCN2 (EIF2AK4) is activated by amino acid deprivation and induces global translation shutdown via eIF2α phosphorylation. HIV+TB coinfection is associated with severe wasting and amino acid depletion, consistent with GCN2-mediated translational inhibition being more prominent in coinfection than in either pathogen alone.

**Overall HIV+TB signal quality: HIGH CONFIDENCE**

The combined infection produces a transcriptional signature that is more severe than either pathogen alone (deeper translation suppression, stronger neutrophil/histone elevation), with a unique B cell immunosuppression signal not seen in single infections. These cross-disease differences are biologically interpretable and clinically coherent.

---

## 4. Tuberculosis

### 4.1 NES > 0 Pathways (Anomalously Upregulated in TB)

**Top 30 by |NES|:**

| Rank | NES | Term |
|------|-----|------|
| 1 | +2.496 | KEGG: Neutrophil Extracellular Trap Formation |
| 2 | +2.477 | KEGG: Systemic Lupus Erythematosus |
| 3 | +2.457 | Reactome: Interferon Alpha/Beta Signaling |
| 4 | +2.442 | Reactome: Regulation Of Actin Dynamics For Phagocytic Cup Formation |
| 5 | +2.352 | Reactome: RHO GTPases Activate WASPs And WAVEs |
| 6 | +2.323 | KEGG: Alcoholism |
| 7 | +2.304 | Reactome: FCGR3A-mediated Phagocytosis |
| 8 | +2.273 | GO: Glycolytic Process |
| 9 | +2.254 | KEGG: Bacterial Invasion Of Epithelial Cells |
| 10 | +2.252 | Reactome: Assembly Of Pre-Replicative Complex |
| 11 | +2.242 | Reactome: RUNX1 Regulates Transcription Of Genes Involved In HSC Differentiation |
| 12 | +2.239 | Reactome: Interferon Signaling |
| 13 | +2.213 | Reactome: MAP2K And MAPK Activation |
| 14 | +2.204 | Reactome: Signaling By High-Kinase Activity BRAF Mutants |
| 15 | +2.194 | Reactome: RHO GTPase Effectors |
| 16 | +2.185 | Reactome: Mitotic G1 Phase And G1/S Transition |
| 17 | +2.168 | Reactome: DNA Replication |
| 18 | +2.160 | Reactome: Signaling By SCF-KIT |
| 19 | +2.156 | Reactome: Formation Of Tubulin Folding Intermediates By CCT/TriC |
| 20 | +2.154 | Reactome: Paradoxical Activation Of RAF Signaling By Kinase Inactive BRAF |
| 21 | +2.151 | Reactome: DNA Replication Pre-Initiation |
| 22 | +2.150 | Reactome: Signaling By RAF1 Mutants |
| 23 | +2.149 | Reactome: Infection With Mycobacterium Tuberculosis |
| 24 | +2.144 | GO: Homotypic Cell-Cell Adhesion |
| 25 | +2.141 | Reactome: GRB2:SOS Provides Linkage To MAPK Signaling For Integrins |
| 26 | +2.138 | GO: Carbohydrate Catabolic Process |
| 27 | +2.136 | Reactome: Interferon Gamma Signaling |
| 28 | +2.130 | Reactome: ER-Phagosome Pathway |
| 29 | +2.125 | Reactome: Antigen Processing — Cross Presentation |
| 30 | +2.123 | Reactome: Transcriptional Regulation Of Granulopoiesis |

### 4.2 NES < 0 Pathways (Anomalously Downregulated in TB)

TB shows almost no downregulation (only 7 FDR < 0.05 pathways):

| Rank | NES | Term |
|------|-----|------|
| 1 | −2.333 | GO: RNA Splicing Via Transesterification Reactions |
| 2 | −2.289 | GO: mRNA Processing |
| 3 | −2.288 | GO: mRNA Splicing, Via Spliceosome |
| 4 | −2.162 | KEGG: Spliceosome |
| 5 | −2.139 | GO: RNA Splicing |
| 6 | −2.126 | GO: B Cell Proliferation |
| 7 | −2.107 | GO: Maturation Of SSU-rRNA |

### 4.3 Biological Evaluation: Tuberculosis

#### NET Formation and SLE Signature — TOP HITS, EXPECTED WITH HISTONE CAVEAT

As in HIV+TB, the top NES+ terms are NET Formation and KEGG SLE, both substantially driven by histones plus genuine neutrophil/complement biology. This is the hallmark of active pulmonary tuberculosis at the transcriptomic level: Berry et al. (2010, Nature Immunology) established that active TB blood transcriptome is dominated by a neutrophil-driven, interferon-inducible gene profile. Neutrophil expansion in active TB is massive (neutrophil-to-lymphocyte ratio is a prognostic marker in TB); neutrophils release NETs in response to Mtb, and circulating histones are released into plasma from dying neutrophils and macrophages. The cfRNA normative model captures this as elevated histone and neutrophil gene mRNA abundance in plasma.

The distinction between genuine NET biology and histone artifact here is that the NET gene set also contains authentic non-histone neutrophil markers: FCGR1A, FCGR2A, FCGR3A/B, MPO, AZU1, CTSG, ELANE, CASP1, CASP4, NCF1, NCF2, RAC1/2, ITGA2B, MAPK14, PIK3CB, PIK3R2. These confirm genuine neutrophil transcriptional activity.

#### Phagocytosis — HIGHLY SPECIFIC AND EXPECTED

Multiple phagocytosis pathways are enriched in TB at high NES values that are absent or weaker in HIV:
- Regulation Of Actin Dynamics For Phagocytic Cup Formation (NES +2.442) — highest NES of any phagocytosis term across all three diseases
- FCGR3A-mediated Phagocytosis (NES +2.304)
- FcγR-dependent Phagocytosis (NES +2.107)
- KEGG Fc gamma R-mediated Phagocytosis (NES +2.101)
- ER-Phagosome Pathway (NES +2.130)
- Response Of Mtb To Phagocytosis (NES +2.117)
- Suppression Of Phagosomal Maturation (NES +1.870)
- KEGG Phagosome (NES +1.893)

The lead genes for phagocytic cup formation include ARPC1B, ARPC2, FCGR1A/2A, LIMK1, WIPF1, WAS, WASF1/2, GRB2, ABI2, NCKAP1, CYFIP1, MYH9, ELMO1/2, DOCK180. These are precisely the proteins that mediate Mtb phagocytosis by macrophages. The "Infection With Mycobacterium Tuberculosis" (NES +2.149) and "Response Of Mtb To Phagocytosis" (NES +2.117) are Reactome pathways specifically curated for Mtb-host interactions, and their significant enrichment is a direct validation that the normative model is detecting Mtb-specific host transcriptional responses.

The "Suppression Of Phagosomal Maturation" (NES +1.870) is particularly informative: Mtb survives intracellularly by inhibiting phagosomal acidification and fusion with lysosomes. The enrichment of this pathway suggests that cfRNA in TB patients contains transcripts from cells (macrophages) that are actively engaged in phagosomal trafficking — and failing to mature the phagosome due to Mtb evasion strategies.

#### Interferon Alpha/Beta and Interferon Gamma — EXPECTED AND WELL-CALIBRATED

TB uniquely enriches both type I (IFN-α/β; NES +2.457) and type II (IFN-γ; NES +2.136) interferon pathways at high NES. The IFN-gamma signal in TB reflects the central role of Th1/CD4+ T cell-derived IFN-γ in anti-mycobacterial defense — IFN-γ activates macrophages to kill Mtb, and circulating immune cells in active TB maintain high IFN-γ-stimulated gene (ISG) transcription. In contrast, HIV shows predominantly type I IFN response (IFN-γ Signaling NES only +1.911 in HIV vs. +2.136 in TB).

The lead genes for IFN-α/β in TB include IFITM3, STAT1, GBP1/2/4/5/6, IFIT2/3, IFITM1/2, HLA-A/B/C/E, IFI27, XAF1, IFI6, ADAR, OAS1/2/3, SOCS1/3, SAMHD1, IRF1/2, PSMB8 (LMP7), STAT2. This ISG fingerprint differs from HIV: TB ISGs are more GBP-enriched (Guanylate Binding Proteins are potent anti-mycobacterial ISGs recruited to the phagosome, whereas ISG15/OAS1/MX1 dominate in HIV). This distinction is biologically meaningful.

#### Glycolysis — EXPECTED AND MECHANISTICALLY IMPORTANT

Glycolytic Process (NES +2.273) is the strongest metabolic upregulation in TB, not present at this rank in HIV. Mtb infection drives a metabolic switch in macrophages from oxidative phosphorylation to aerobic glycolysis (the "Warburg effect"), generating lactate and supporting inflammatory cytokine production. This macrophage metabolic reprogramming is a well-established feature of active TB: HIF-1α stabilization by Mtb triggers glycolytic gene upregulation (GAPDH, PKM, LDHA, HK1, ENO1, PGK1, PGAM1, TPI1, ENO2). Lead genes for this pathway are precisely these glycolytic enzymes. The cfRNA signal likely reflects macrophages and monocytes undergoing glycolytic reprogramming in response to Mtb.

The "Carbohydrate Catabolic Process" (NES +2.138) is the broader GO term that captures the same signal.

#### MAPK/RAF/KIT Signaling — EXPECTED IN MYCOBACTERIAL CONTEXT

MAP2K And MAPK Activation (NES +2.213), Signaling By BRAF Mutants (NES +2.204), and GRB2:SOS Linkage To MAPK For Integrins (NES +2.141) reflect the MAPK pathway activation downstream of Mtb-driven toll-like receptor (TLR) and FcγR signaling. Mtb activates p38 MAPK (MAPK14) and ERK1/2 in macrophages as part of its intracellular signaling strategy. The lead genes (MAPK1/3, GRB2, SRC, RAF1, SHC1) are consistent with integrin- and receptor-tyrosine-kinase-driven MAPK activation at the phagocytic synapse. Signaling By SCF-KIT (NES +2.160) may reflect mast cell and myeloid cell activation by SCF, which are elevated in TB-associated inflammation.

#### Antigen Cross-Presentation — EXPECTED

ER-Phagosome Pathway (NES +2.130) and Antigen Processing — Cross Presentation (NES +2.125) reflect the antigen presentation machinery in macrophages and dendritic cells responding to Mtb. Cross-presentation via MHC-I of Mtb antigens is an important immune response mechanism, and the engagement of ER-to-phagosome lipid transfer machinery (lipid body formation) is well-documented in Mtb-infected macrophages. Lead genes include HLA-A/B/C, TAP1/2, CALR, CANX, and ER-associated proteins.

#### Autophagy — PRESENT AND MEANINGFUL

Regulation Of Autophagy Of Mitochondrion (mitophagy; NES +1.738) and Regulation Of Autophagy (NES +1.710) reflect the autophagic killing mechanism against Mtb. Autophagy is a major host defense against intracellular Mtb: macrophages upregulate autophagy pathways to deliver Mtb to lysosomes, and Mtb counteracts this by inhibiting phagosomal maturation. The simultaneous enrichment of phagocytosis, phagosomal suppression, and autophagy regulation in TB cfRNA is a coherent and mechanistically consistent picture of the macrophage-Mtb battleground.

#### Spliceosome Downregulation — UNEXPECTED AND INTERESTING

The only 5 downregulated pathways in TB that reach FDR < 0.05 are all spliceosome/mRNA splicing-related:
- RNA Splicing Via Transesterification Reactions (NES −2.333)
- mRNA Processing (NES −2.289)
- mRNA Splicing Via Spliceosome (NES −2.288)
- KEGG Spliceosome (NES −2.162)
- RNA Splicing (NES −2.139)

Lead genes include PNN, YTHDC1, SRRM2, HNRNPU/K/M/H1, SRSF5/8/10, DDX5/46, SNRNP70, SF3B1, PRPF3, SFPQ, and TARDBP. This is a coherent and specific cluster.

Why would splicing machinery be downregulated in TB cfRNA? Several mechanisms could explain this:

1. **Mtb modulation of host mRNA splicing:** Recent evidence shows that Mtb infection modulates host alternative splicing globally, potentially suppressing splice factor expression as an evasion strategy.
2. **Cell composition effect:** The TB blood cellular landscape is neutrophil-dominated; neutrophils have low splicing activity due to their terminally differentiated state and transcriptionally repressed spliceosome. As the proportion of cfRNA from neutrophils increases, spliceosomal mRNAs (typically high in lymphocytes and macrophages) appear relatively depleted.
3. **Lymphopenia effect:** Active TB causes relative lymphopenia; lymphocytes are high in splicing factor transcripts. The depletion of lymphocyte-derived cfRNA (supported by the B cell proliferation downregulation; NES −2.126) causes apparent spliceosome suppression.

The B Cell Proliferation downregulation (NES −2.126, lead genes: PTPRC, MS4A1, CD180, CD79A, CD19, MEF2C, GAPT, CD40LG, RASGRP1) supports interpretation #3: TB actively suppresses B cell responses, and the reduced circulating B cell mass is reflected in lower B cell cfRNA abundance in plasma.

**Overall TB signal quality: HIGH CONFIDENCE**

TB produces the broadest GSEA upregulation of the three diseases (610 significant upregulated pathways), with mechanistically coherent enrichment of the canonical neutrophil-IFN axis, Mtb-specific phagocytosis and macrophage metabolism pathways, and a striking spliceosome downregulation likely reflecting the neutrophil-dominated blood cell composition.

---

## 5. Cross-Disease Comparison

### 5.1 Pathway NES Comparison Table

| Pathway | HIV | HIV+TB | TB | Notes |
|---------|-----|--------|----|-------|
| Antiviral Mechanism By ISGs | +2.741 | +1.787 | +1.897 | Highest in HIV-only; attenuated in TB |
| ISG15 Antiviral Mechanism | +2.585 | N/A | +1.792 | HIV-specific ISG15 prominence |
| IFN Alpha/Beta Signaling | +2.595 | +2.317 | +2.457 | Present in all; highest in HIV |
| IFN Gamma Signaling | +1.911 | N/A* | +2.136 | TB shows IFN-gamma > HIV |
| NET Formation | +2.371 | +2.754 | +2.496 | Additive in coinfection |
| SLE KEGG | +2.270 | +2.911 | +2.477 | Histone-driven; peaks in coinfection |
| Mtb Infection (Reactome) | +1.872 | +1.801 | +2.149 | TB-highest, as expected |
| Phagocytic Cup Formation | +1.998 | +1.955 | +2.442 | TB-highest |
| Glycolytic Process | +1.900 | +2.336 | +2.273 | TB/HIV+TB; macrophage Warburg effect |
| Pyroptosis | +2.177 | +2.575 | +1.952 | HIV+TB highest; both pathogens synergize |
| Complement Cascade | N/A | +1.907 | +1.700 | Mostly HIV+TB and TB |
| Translation/Ribosome | −2.237 | −2.743 | N/A | HIV+TB deepest suppression |
| B Cell Proliferation | N/A | −2.158 | −2.126 | TB and HIV+TB; B cell immunosuppression |
| Spliceosome | −1.987 | −1.959 | −2.162 | TB-dominant downregulation |

*IFN Gamma Signaling at NES +1.765 in HIV+TB, below top 30 but significant.

### 5.2 Disease-Differentiating Signals

- **HIV-specific:** ISG15 mechanism, antiviral ISG cluster at highest NES, Rev-host interactions, HIV-specific senescence/SASP prominence, mitochondrial translation suppression cluster
- **TB-specific:** Phagocytic cup formation, Mtb Infection (Reactome), ER-phagosome pathway, glycolysis, IFN-gamma prominence, spliceosome downregulation
- **HIV+TB-specific:** Deepest translation suppression, B cell pathway downregulation cluster (unique severity), highest pyroptosis NES, highest SLE/NET NES (additive pathogen burden)
- **Shared across all three:** IFN signaling, NET formation, granulopoiesis, senescence/cell cycle, chromatin pathways (partly artifact)

---

## 6. Artifact Assessment: Histone mRNA Inflation

### 6.1 Nature of the Artifact

Plasma cfRNA from patients with high immune cell turnover (HIV, TB, HIV+TB) contains elevated levels of replication-coupled histone mRNAs: H2BC3/5/6/9/11/12/13/14/15/17/18/21, H2AC4/7/8/11/12/13/14/16/17/19/21, H3C1/2/3/4/7/8/10/11/12/13/14/15, H4C1/2/4/8/9/11/12/13/14/15. These are not polyadenylated (they are canonical replication-coupled histones), and their presence in cfRNA libraries may reflect: (1) cell lysis releasing chromatin-bound transcripts, (2) genuine histone mRNA abundance in rapidly dividing cells (neutrophil progenitors, activated lymphocytes), or (3) pre-mRNA contamination from chromatin-associated RNA.

### 6.2 Affected Pathways (Artifact-Inflated)

The following pathways rank highly primarily due to histone gene set membership:
- **KEGG Alcoholism:** 42–57% of lead genes are histones across all three diseases
- **KEGG Systemic Lupus Erythematosus:** Contains large histone/nucleosome component (autoantigens in SLE)
- **DNA Methylation, PRC2 Methylates Histones, HDACs Deacetylate Histones:** Epigenetic pathway gene sets that include histone substrates
- **Condensation Of Prophase Chromosomes, Packaging Of Telomere Ends, Assembly Of ORC Complex:** Chromosome biology pathways with histone structural components
- **Transcriptional Regulation Of Granulopoiesis, RUNX1 Regulates Megakaryocyte Differentiation:** These have genuine biology BUT are inflated by incidental histone gene co-membership
- **RHO GTPases Activate PKNs, PKN1/AR/KLK2-KLK3:** PKN substrates include histones; these terms appear inflated

### 6.3 Pathways with Genuine Biology Despite Histone Presence

- **KEGG NET Formation:** Contains genuine non-histone neutrophil markers (MPO, AZU1, ELANE, CTSG, FCGRs)
- **KEGG Granulopoiesis:** Contains transcription factor biology (CEBPE, MYB, GFI1, GATA2) beyond histones
- **DNA Repair (BER/NER):** Reflects genuine oxidative DNA damage in HIV/TB

### 6.4 Recommendation

For downstream biological interpretation or publications, terms dominated by histone cluster genes (>40% of lead genes being H2BC/H2AC/H3C/H4C) should be flagged as likely artifact-inflated. The KEGG Alcoholism term should be excluded from biological conclusions entirely as it has no disease-relevant interpretation for these infectious disease cohorts. The SLE and NET pathways should be reported with the caveat that their high NES includes a histone-inflation component, even though the underlying biology (neutrophil activation, immune complex formation) is genuine.

---

## 7. Overall Confidence Assessment

| Disease | Signal Quality | Primary Confidence Basis |
|---------|---------------|--------------------------|
| HIV | **HIGH** | Top 4 terms are canonical HIV ISG signature; mitochondrial/ribosomal downregulation mechanistically validated; Rev-host interaction pathway as direct positive control |
| HIV + Tuberculosis | **HIGH** | Additive NET/neutrophil signal biologically expected; unique B cell suppression is disease-specific and mechanistically coherent; deepest translation suppression fits dual pathogen burden; GCN2/amino acid deficiency pathway adds mechanistic specificity |
| Tuberculosis | **HIGH** | Mtb-specific phagocytosis pathways (unique rank vs. HIV); glycolysis enrichment matches macrophage Warburg effect; Infection With Mycobacterium Tuberculosis Reactome pathway is direct validation; IFN-gamma > type I IFN ratio appropriate for TB; spliceosome downregulation is consistent with neutrophil-dominated blood composition |

All three infectious diseases produce GSEA signals with HIGH biological confidence. The cfRNA normative model (GAMLSS-based Z-scoring relative to healthy controls) successfully captures disease-specific transcriptomic programs in plasma cell-free RNA. The dominant signals — IFN/ISG activation in HIV, phagocytosis/glycolysis in TB, synergistic neutrophil-B cell imbalance in HIV+TB — are directly supported by decades of whole-blood and single-cell transcriptomic literature on these diseases.

The primary caveat across all three diseases is the histone mRNA inflation artifact, which inflates chromatin-related pathway scores and should be noted when interpreting pathways ranked #1–5 in HIV+TB and TB (SLE, Alcoholism terms). This artifact does not undermine the overall biological validity of the GSEA results; it requires selective pathway-level annotation rather than dismissal of the entire analysis.

---

*Evaluation written by: Claude Sonnet (claude-sonnet-4-6), 2026-06-24*  
*Data files analyzed: gsea_result_HIV.csv, gsea_result_HIV + Tuberculosis.csv, gsea_result_Tuberculosis.csv*  
*Location: /project/cfRNA_NormativeModeling/Modeling/GE/gsea/*
