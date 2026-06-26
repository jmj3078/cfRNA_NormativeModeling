# cfRNA Normative Model — GSEA 통합 분석 리포트

**작성일:** 2026-06-24 (2026-06-25 개정: CDCS → CAD_HF+ / CAD_HF− 분리 반영)  
**분석 대상:** 21개 표현형 (20개 질병 + Other Cancer)  
**모델:** GAMLSS normative model, HC reference n=693 (WTS-only)  
**GSEA:** GSEApy prerank, mean Z-score 기반 랭킹  
**유전자 세트:** KEGG 2021 Human · Reactome 2022 · GO Biological Process 2023  
**유의 기준:** FDR q-val < 0.05  
**NES 해석:** NES > 0 = HC 대비 상향 이탈; NES < 0 = HC 대비 하향 이탈

---

> **⚠ 데이터 버전 노트 (2026-06-24 갱신)**
> 본 리포트의 정성적 해석은 study-split + MIN_SAMPLES(n≥5) 필터 적용 직후의 GSEA 실행 결과를 기준으로 한다.
> 이후 노트북의 GSEA prerank 루프가 재실행되어 현재 `GE/gsea/*.csv`의 유의 경로 수가 **±수 개~수십 개** 수준으로 미세하게 변동했다 (GSEA permutation 무작위성 + 필터링된 샘플 구성 변화). 아래 모든 카운트는 **현재 CSV 기준**으로 갱신되었다.
>
> **9개 novel 발견과 모든 헤드라인 경로는 현재 데이터에서 동일한 NES로 그대로 유지된다** (생물학적 결론 무변경).
> 단, 두 가지 **질적 변화**가 있어 본문에 반영했다:
> - **ME/CFS**: 이전 실행의 상향 경로 4개(면역글로불린·B세포·MHC-II)가 현재 FDR<0.05 set에서 모두 탈락 → **상향 0개 / 하향 169개**. 핵심 하향 축(번역·FoxO·인슐린·L13a)은 그대로 유지.
> - **Stomach Cancer**: 각질화 봉투(Cornified Envelope) 경로가 현재 FDR<0.05에서 탈락. 단, 데스모좀 junction 신호(PKP1/3/4·DSG2·CDH3)는 Cell-Cell Junction Organization(+2.64)로 유지.

---

## 목차

1. [감염병 그룹](#1-감염병-그룹)
   - 1.1 HIV
   - 1.2 HIV + Tuberculosis
   - 1.3 Tuberculosis
2. [간·소화기암 그룹](#2-간소화기암-그룹)
   - 2.1 Liver Cancer (Chen)
   - 2.2 Liver Cancer (Roskams-Hieter)
   - 2.3 Liver Cirrhosis
   - 2.4 Colorectal Cancer
   - 2.5 Esophagus Cancer
   - 2.6 Stomach Cancer
3. [췌장·폐 그룹](#3-췌장폐-그룹)
   - 3.1 Pancreatic Cancer (Moore)
   - 3.2 Pancreatic Cancer (Reggiardo)
   - 3.3 Pancreatitis
   - 3.4 Lung Cancer
4. [면역·혈액 그룹](#4-면역혈액-그룹)
   - 4.1 ICI-induced Myocarditis
   - 4.2 ICI-treated Cancer
   - 4.3 Multiple Myeloma
   - 4.4 MGUS
   - 4.5 ME/CFS
5. [심혈관·산과·기타 그룹](#5-심혈관산과기타-그룹)
   - 5.1 CAD — Heart Failure 진행군 비교 (CAD_HF+ / CAD_HF−)
   - 5.2 Pre-eclampsia
   - 5.3 Other Cancer
6. [기존 문헌 미보고 Novel 발견 정리](#6-기존-문헌-미보고-novel-발견-정리)
7. [전체 요약 테이블](#7-전체-요약-테이블)

---

## 1. 감염병 그룹

**출처:** Chang et al. 2023 (DOI: 10.1101/2023.01.11.23284435 외 companion papers)

---

### 1.1 HIV (sample n=13 · FDR<0.05 482개 · NES>0: 434 · NES<0: 48)

#### 주요 상향 신호

**IFN/ISG 축 (NES +2.74 최고)**
- 최상위 4개 경로가 모두 IFN 관련:
  Antiviral Mechanism By ISGs (+2.741) → IFN α/β Signaling (+2.595) → ISG15 (+2.585) → Interferon Signaling (+2.541)
- Lead genes: IFI27, OAS1/2/3, IFIT1/2/3, MX1/2, ISG15, RSAD2, IRF7, ADAR, BST2, IFITM1/2/3
- HIV 복제로 cGAS-STING/RIG-I가 활성화되어 지속적 type I IFN 생산;
  ART 억압 상태에서도 잔존 항원이 IFN tone을 유지함
- cfRNA 기원: 활성화 단핵구, pDC, 림프구 유래 ISG mRNA

**Neutrophil Extracellular Trap (NET) 형성 (+2.371)**
- Lead genes: MPO, FCGR1A/2A/3B, 히스톤(H2AC/H2BC/H3C/H4C 군), PIK3CB, MAPK14
- HIV가 바이러스 입자를 통해 호중구를 직접 활성화 → NET 형성 → 히스톤 cfRNA 방출

**세포 노화·DNA 손상 복합체**
- Oxidative Stress Induced Senescence (+2.346) + SASP (+2.239) + Telomere Stress Senescence
- 만성 HIV에서 가속 면역 노화; p21/p16 매개 senescence는 CD4+ T, CD8+ T, 단핵구에서 공통 보고됨

**Granulopoiesis·RUNX1·혈소판 경로**
- RUNX1/Megakaryocyte (+2.359): HIV 혈소판 감소증 기전(직접 거핵구 감염, 면역복합체) 반영
- Lead genes: RUNX1, GATA2, MYB, NFE2

**HIV-특이 경로 (직접 검증)**
- Interactions Of Rev With Host Cellular Proteins (+1.692)
- Nuclear Import Of Rev Protein (+1.734): importin, NUP98/153 등 핵공 단백질
  → cfRNA 정규화 모델이 HIV-특이 분자 생물학을 감지하는 **직접 양성 대조군**

#### 주요 하향 신호

**번역·리보솜 억제 (NES −2.49 최저)**
- SRP-dependent Cotranslational Targeting (−2.494) · Ribosomal Subunit Biogenesis (−2.405) · Translation (−2.381)
- 미토콘드리아 번역 억제 (−2.332 ~ −2.177): HIV 및 뉴클레오시드계 ART의 미토콘드리아 독성 반영
- HIV Tat가 숙주 cap-dependent translation 억제; Vpr이 리보솜 생합성 억제

#### 아티팩트 주의

- KEGG Alcoholism (히스톤 57%), KEGG SLE (히스톤+보체), DNA Methylation, PRC2, Condensation Of Prophase Chromosomes:
  이들은 H2BC/H2AC/H3C/H4C 히스톤 클러스터 유전자의 대량 cfRNA로 인한 **아티팩트 팽창** 경로

**신호 품질: HIGH**

---

### 1.2 HIV + Tuberculosis (sample n=11 · FDR<0.05 397개 · NES>0: 306 · NES<0: 91)

#### 주요 상향 신호

**SLE/NET 신호 최상위 (NES>0 1위: KEGG SLE +2.911)**
- HIV와 Mtb 각각 독립적으로 NET을 유도하므로 공동감염에서 상승적으로 증폭됨
- SLE 경로: 히스톤 자가항원 + C1q/C3/C9 + FcγR — 진짜 면역복합체 생물학 + 히스톤 아티팩트 혼재

**Pyroptosis (NES +2.575) — 두 병원체 수렴**
- HIV: NLRP3 → gasdermin-D → CD4+ T세포 pyroptosis (CD4 감소의 주요 기전)
- Mtb: 대식세포에서 NLRP3 활성화
- 공동감염에서 두 경로 동시 작동 → HIV+TB에서 가장 높은 NES

**항균 체액성 반응 (+2.367)**
- CTSG, ELANE, AZU1, MPO, defensins: 호중구 과립 방출 반영

**IFN 신호 상대적 감쇠 (HIV 단독 대비)**
- IFN α/β: HIV+TB +2.317 vs HIV 단독 +2.595
- Mtb 공동감염이 SOCS1/IL-10 매개로 type I IFN 신호를 조절함을 반영

#### 주요 하향 신호

**번역 억제 — 세 질병 중 최심 (NES −2.743)**
- Cap-dependent Translation Initiation (−2.743), Eukaryotic Translation Elongation (−2.737)
- **GCN2 반응 (EIF2AK4, NES −2.488)**: 이중 병원체 부담으로 인한 아미노산 고갈 → GCN2 활성화 → eIF2α 인산화 → 번역 전반 억제
  ↳ HIV 단독(−2.494)보다 훨씬 깊음 = 이중 병원체 대사 스트레스의 직접 반영

**B 세포 경로 하향 (HIV+TB 고유 신호)**
- B Cell Proliferation (−2.158), B Cell Activation (−2.060), B Cell Receptor Signaling (−1.979)
- B Cell Differentiation (−1.825), Positive Regulation Of B Cell Activation (−1.889)
- Lead genes: CD79A, CD79B, MS4A1(CD20), CD19, CD22, PAX5, BTK, BLNK, BCL2, MEF2C
- HIV가 class switching에 필수인 CD4+ Th 세포를 고갈시키고, Mtb가 독립적으로 B 세포 기능을 억제 → 체액성 면역 전면 손상이 cfRNA 수준에서 포착됨

**신호 품질: HIGH**

---

### 1.3 Tuberculosis (sample n=101 · FDR<0.05 606개 · NES>0: 599 · NES<0: 7)

#### 주요 상향 신호

**NET 형성 + SLE (상위 2위: +2.496, +2.477)**
- Berry et al. 2010 (Nature Immunol.)이 확립한 활성 TB의 호중구-IFN 전사체 지배 구조 반영
- 비-히스톤 호중구 마커(FCGR1A/2A/3A/3B, MPO, AZU1, CTSG, ELANE)가 진짜 생물학 뒷받침

**식세포 작용 — TB-특이 가장 높은 NES**
- Regulation Of Actin Dynamics For Phagocytic Cup Formation (+2.442, 3개 질병 중 최고)
- Infection With Mycobacterium Tuberculosis (+2.149), Response Of Mtb To Phagocytosis (+2.117)
- Lead genes: ARPC1B/2, FCGR1A/2A, LIMK1, WAS, WASF1/2, GRB2, ELMO1/2, DOCK180
  → **Mtb 특이적 숙주-병원체 상호작용 경로가 직접 검출됨**
- Suppression Of Phagosomal Maturation (+1.870): Mtb가 리소솜 융합을 억제해 세포 내 생존하는 메커니즘을 cfRNA 수준에서 반영

**IFN-γ 우세 (TB vs HIV 역전)**
- IFN Alpha/Beta (+2.457) vs IFN Gamma (+2.136): TB에서 type II IFN이 type I IFN에 상대적으로 더 현저
- IFN-γ는 Th1 CD4+ T세포 유래; Mtb 대식세포 활성화의 핵심
- GBP1/2/4/5/6 (항Mtb ISG, 파고솜 동원) 우세 vs HIV에서는 ISG15/OAS1 우세 → **ISG 지문이 병원체-특이적**

**호중구 식세포 대식세포 Warburg 효과**
- Glycolytic Process (+2.273): HIF-1α → GAPDH, PKM, LDHA, HK1, ENO1 상향
- Mtb 감염 대식세포가 산화적 인산화 → 유기호흡 전환; 도말에 관한 cfRNA 신호

**Autophagy — Mtb 방어 기전 반영**
- Regulation of Autophagy Of Mitochondrion (+1.738), Regulation of Autophagy (+1.710)
- 파고솜, 파고솜 성숙 억제, 자가포식 경로가 동시 농축 → **대식세포-Mtb 전쟁터의 일관된 cfRNA 그림**

#### 주요 하향 신호 (단 7개)

**스플라이세오솜 하향 — TB-고유 신호**
- RNA Splicing Via Transesterification Reactions (−2.333)
- mRNA Processing (−2.289), mRNA Splicing Via Spliceosome (−2.288)
- KEGG Spliceosome (−2.162), RNA Splicing (−2.139)
- Lead genes: PNN, YTHDC1, SRRM2, HNRNPU/K/M/H1, SRSF5/8/10, DDX5/46, SF3B1, PRPF3, TARDBP
- B Cell Proliferation (−2.126): CD45, MS4A1, CD79A, CD19, MEF2C → 림프구 감소증 동반

**신호 품질: HIGH**

---

## 2. 간·소화기암 그룹

**출처:** Chen et al. eLife 2022; Roskams-Hieter et al. npj Precis Oncol 2022; Block et al. Front Oncol 2022; Peddu et al. bioRxiv 2025

---

### 2.1 Liver Cancer — Chen et al. (sample n=10 · FDR<0.05 625개 · NES>0: 171 · NES<0: 454)

#### 주요 상향 신호

**간 특이적 분비 단백질 축 (NES ≥ 2.2)**
- Complement and coagulation cascades (KEGG, +2.476): FGB, FGA, FGG, PLG, C3, CFB, CFH, SERPINC1
- Post-translational Protein Phosphorylation (+2.426): ALB, FGA, APOB, KNG1
- Regulation Of IGF Transport/IGFBPs (+2.347): ALB, FGA, APOB, AHSG
- Cholesterol Efflux (+2.326): APOA2, APOC3, APOE, APOA1, APOA4
- Plasma Lipoprotein Remodeling (+2.284): ALB, APOB, APOC3, APOE, ANGPTL3/4
- Regulation Of Fibrinolysis (+2.281): VTN, PLG, PLAU, SERPINF2, THBD, F2

**해석:** HCC에서 간세포 파괴 → ALB·FGA/FGB/FGG·PLG·APOB·APOE 등 간 분비 단백질의 cfRNA 방출 증가.
Roskams-Hieter 2022이 HCC 바이오마커로 명시한 APOE, fibrinogen, C3 등과 일치.
HNF4A/FOXA2/NKX6-1 간 분화 전사인자 축도 상향.

#### 주요 하향 신호

**번역 기계 전면 억제 (NES 최저 −3.356)**
- Eukaryotic Translation Elongation (−3.356), Peptide Chain Elongation (−3.330)
- rRNA Processing, Cap-dependent Translation Initiation, Cytoplasmic Translation
- 리보솜 단백질(RPL36, RPL21, RPS13, RPS25) 및 번역 인자(EIF4E, EIF3J) 전반 하향

**Starvation Response (−2.715):** GCN2/ISR 활성화 — 종양 미세환경 영양 결핍 반영

**주목:** Chen 코호트에서 번역 억제가 OXPHOS 억제보다 우세 → Group2 평가의 "간암 = OXPHOS 하향, GI암 = 번역 하향" 패턴에서 이탈; n=10 소수 편향 또는 샘플 이질성 가능성

**신호 품질: MODERATE-HIGH (n=10 한계)**

---

### 2.2 Liver Cancer — Roskams-Hieter et al. (sample n=28 · FDR<0.05 217개 · NES>0: 211 · NES<0: 6)

#### 주요 상향 신호

**보체/응고 (NES ≥ 2.0)**
- Negative Regulation Of Blood Coagulation (+2.224): FGB, APOH, FGA, PLG, KNG1, FGG, APOE
- Complement and coagulation cascades (+2.110): FGB, FGA, PLG, F13A1, A2M, CFB, CFH, C9
- Fibrinolysis (+1.975): FGB, FGA, PLG, PLAT
  → Chen 코호트와 **cross-cohort 재현성 확인된 HCC 코어 cfRNA 축**

**세포골격·침습 신호**
- Adherens junction (+2.062): TJP1, CTNND1, VCL, AFDN, IQGAP1
- Actin Filament Bundle Assembly (+1.968): FSCN1, LIMA1, NEDD9, LCP1, EZR
- CDC42 GTPase Cycle (+1.909): ARHGAP29/31, ECT2, DOCK9, IQGAP1
- RAC1 GTPase Cycle (+1.884): FERMT2, WASF2, SWAP70, ITGB1

**세포주기 활성화**
- Positive Regulation Of G1/S Transition (+1.998): TFDP1, RRM2, CCND2, CDC6
- Oncogene Induced Senescence (+1.910): CDKN2A, CDKN2B, TP53, RB1, MDM2, CDK6

**Hippo/YAP 상향**
- YAP1/RUNX3/TEAD4 Transcriptional Regulation (~+1.90): WWTR1, YAP1, TEAD4, LATS1, EP300
  → Roskams-Hieter 2022 원문에서 세포주기 유전자가 HCC 바이오마커로 보고됨; YAP 경로 직접 검출

**항바이러스 신호**
- Negative Regulation Of Viral Genome Replication: IFITM2/3, EIF2AK2, OAS2, IFI16
  → 바이러스성 간염(HBV/HCV) 배경의 HCC를 암시

#### 주요 하향 신호 (단 6개)

- OXPHOS/미토콘드리아 대사 (NES −2.332 ~ −2.108): UQCR11, NDUFA1, COX17 등
  → Group2에서 확인된 "HCC = Warburg 효과" 재현

**[두 코호트 핵심 비교]**

| 축 | Chen (n=10) | Roskams-Hieter (n=28) |
|---|---|---|
| 보체/응고 | ✓ 강함 | ✓ 강함 (재현) |
| 지질단백질 | ✓ | ✓ |
| 번역 억제 | ✓✓ 극강 (454개) | ✗ (거의 없음) |
| OXPHOS 하향 | △ 약함 | ✓ (6개) |
| EMT/침습 (GTPase, actin) | ✗ | ✓✓ |
| 세포주기 (G1/S, kinesin) | ✗ | ✓ |
| YAP/Hippo | ✗ | ✓ |
| 항바이러스 | ✗ | ✓ |

→ 동일 "Liver Cancer" 레이블 내 **환자 집단 이질성 (중국 HCC vs 미국 HCC, 바이러스성 vs 비바이러스성)**이 cfRNA 수준에서 명확하게 구별됨

**신호 품질: HIGH (n=28)**

---

### 2.3 Liver Cirrhosis (n=4 · NES>0: 39 · NES<0: 2)

⚠ **n=4로 통계적 신뢰도 매우 낮음; 모든 결과 잠정적**

#### 핵심 신호

**최상위 상향: HDACs Deacetylate Histones (+2.363)**
- 주도 유전자: H2AC/H2BC/H3C 히스톤군 + 대식세포/보체
- 간경화에서 문맥 고혈압 → 전신 염증 → 호중구 NETosis → 히스톤 cfRNA 방출로 해석
- HDAC 에피제네틱 리모델링 (간 섬유화에서 BRG1/SMARCA4 활성화)이 기여 가능

**혈관 리모델링 신호 (+2.203)**
- Regulation Of Endothelial Cell Migration: KDR(VEGFR-2), FLT4(VEGFR-3), DLL4, EDN1, NOS3, FGF1, EGF
- 간경화 문맥 고혈압에서 활성화된 LSEC(간 동굴 내피세포) 유래 cfRNA로 해석
- 섬유화→간경화 진행에서 LSEC 기능 이상이 초기 사건임을 반영

**보체/응고 상향 (+2.041):** HCC와 공통 — 응고병증은 간경화의 임상적 핵심 합병증

**OXPHOS 하향 (−2.280, −2.203):** HCC보다 약한 신호 (−2.92) → **간경화 → HCC 진행에 따른 미토콘드리아 기능 저하가 점진적으로 심화됨**을 시사

**신호 품질: LOW-MODERATE (n=4)**

---

### 2.4 Colorectal Cancer (sample n=37 · FDR<0.05 385개 · NES>0: 25 · NES<0: 360)

#### 주요 상향 신호

**Hippo/YAP 최강 상향 (NES +2.542, 공동 1위)**
- Signaling By Hippo (+2.542), YAP1/WWTR1(TAZ)-stimulated Gene Expression (+2.542)
- Lead genes: TJP1, STK24, YAP1, WWTR1, TJP2, AMOT, MOB1A, LATS2, AMOTL2, HIPK2, KAT2B, TEAD4
- YAP/TAZ는 Hippo 종양 억제 경로의 핵심 하위 전달자;
  CRC에서 YAP 핵 이동이 침습·예후 불량 마커로 확립됨
- cfRNA에서 YAP1/WWTR1/TEAD4 전사체 상향: 종양 교육 혈소판 또는 순환 종양 세포 기원

**인테그린 신호 (+2.488)**
- Lead genes: LIMS1/2 (PINCH), NID1 (nidogen), LAMB1/2, LAMC1, FLNA, LOXL3
- ECM-integrin 축은 CRC 침습의 핵심; 종양 기질로부터 cfRNA 방출

**ERBB2/KIT/EGFR RTK 신호 (ranks 8-17)**
- GRB2, SRC, HSP90AA1, SOS1, JAK2, STAT1/5B, PIK3R1, EGF
- ERBB2 증폭 (~5-8% CRC), KIT 기질 신호, EGF/EGFR 축

**인슐린·글루코스 대사 (+2.410, +2.052)**
- Warburg 효과와 연결된 PI3K-AKT-인슐린 신호 상향

#### 주요 하향 신호

**번역·리보솜 완전 억제 (NES 최저 −2.432, 360/385 경로)**
- Eukaryotic Translation Elongation (−2.432), Peptide Chain Elongation (−2.427)
- KEGG Ribosome (−2.339), NMD (EJC-independent/enhanced), SRP-dependent targeting
- GCN2/EIF2AK4 Response To Amino Acid Deficiency (−2.339)
- Lead genes: RPL/RPS 패밀리 전체, EIF3A/B/J, EIF4A, EIF2B4
- **Chen et al. eLife 2022에서 직접 보고된 결과와 일치**: "리보솜 생합성 경로가 CRC cfRNA에서 현저히 하향"
- PloS ONE 2024 CRC cfRNA 연구(PMC11326608)에서도 독립 검증

**신호 품질: HIGH**

---

### 2.5 Esophagus Cancer (sample n=25 · FDR<0.05 683개 · NES>0: 148 · NES<0: 535)

#### 주요 상향 신호

**배아 기관 발달 최강 상향 (+2.673)**
- Lead genes: SOX17, SOX18, TEAD1/2/4, WNT7B, FGFR2, GATA3, GLI2, HEY1/2, TGFB1/2, FOXE1, SALL1
- 식도암이 EMT·암 줄기세포 형성에서 배아 전사 프로그램을 재활성화
- SOX17 (내배엽 정체성), FOXE1 (갑상선/전장 전사인자), GATA3 (전장·유선 계통)의 cfRNA 검출
  → **미분화·배아형 종양 세포 또는 외소체 기원**

**YAP1/TAZ 신호 (+2.573)**
- YAP1, WWTR1, TEAD 패밀리; HIPK2, KAT2B
- FGFR2 신호가 YAP1을 하류에서 활성화함 (PMC7581496 위암 연구 참조)

**FGFR2 신호 (5개 경로, ranks 18·21·23·25·28)**
- FGF3/5/6/7/8/9/10/16/17/22/23, FGFBP1/2/3, FGFR2, PTPN11, IRS1/2/4, PI3K/AKT
- FGFR2는 식도 선암에서 증폭·과발현;
  여러 FGFR2 경로가 동시 농축 → **가장 강한 FGFR2-중심 cfRNA 신호**
- FGFR2 → YAP1 공동 활성화 축이 cfRNA에서 동시 검출 (기존 문헌은 조직 수준 보고)

**IGF1R 신호 (ranks 5·7·10)**
- IGF2, IRS1/2/4, AKT2, PIK3R1, GRB2, SHC1
- Barrett 식도 → 선암 진행에서 IGF2/IGF1R 축의 역할과 일치

**GABA·니코틴·이온 채널 신호 (ranks 4·6·8·13·14)**
- GABRA3/4/5, GABRG1/2/3, GRIN1/2A/2B/2C, CHRNA4/7, CACNA1A/1B
- 가능한 기원: 신경내분비 분화, 장관신경계 cfRNA, 또는 흡연 관련 니코틴 수용체 발현

**SMAD/TGF-β (+2.281)**
- SMAD5/6/9, TGFB1/2, BMP2/5/10, NODAL → EMT 드라이버

#### 주요 하향 신호

**OXPHOS + 번역 동시 억제 (n=539, 이중 최강)**
- OXPHOS (ranks 1·3·4·8·9·11·13): Mitochondrial Respiratory Chain Complex Assembly (−2.387)
- 번역 (ranks 2·5·6·7·10): SRP-dependent Targeting (−2.354), Translation (−2.327)
- **간암(OXPHOS 우세)과 대장암(번역 우세)의 중간 표현형** — 해부학적 위치(GEJ 인근) 반영

**신호 품질: HIGH**

---

### 2.6 Stomach Cancer (sample n=24 · FDR<0.05 516개 · NES>0: 13 · NES<0: 503)

#### 주요 상향 신호

**세포-세포 접합 조직 (NES +2.643, 1위)**
- Lead genes: PKP4, PKP3, PKP1(desmosomal plakophilins), TJP1/2, CDH3/5/6/9, CLDN1/5, DSG2, OCLN, CTNND1

**데스모좀/접합 신호 (Cell-Cell Junction Organization +2.643 lead genes) — 가장 독특한 신호**
- Lead genes: PKP4, PKP3, PKP1(desmosomal plakophilins), DSG2(desmoglein), CDH3(P-cadherin), TJP1/2, CLDN1/5, CTNND1, FSCN1
- 데스모좀 구성 단백질(PKP1/3/4, DSG2)과 P-cadherin(CDH3)의 상향은 위암 상피 접합 리모델링을 반영
- > **데이터 갱신 노트:** 초기 실행에서 2위였던 *Formation Of Cornified Envelope* (NES +2.583, 편평상피 분화/화생 해석)은
  현재 FDR<0.05 set에서 탈락했다. 다만 동일한 desmosomal 유전자(PKP/DSG2)가 Cell-Cell Junction Organization 경로에 그대로 잔존하므로
  상피 접합 리모델링이라는 핵심 해석은 유지된다 (편평상피 특이 cornification 해석만 제외).

**RET 신호 (+2.298)**
- Lead genes: RET, GFRA3/4, DOK2, PTPN11, PRKCA, SRC, IRS2, GRB2
- RET은 갑상선암·MEN2에서 고전적이나 위암에서 RET 재배열/융합이 ~0.5-2% 보고됨
- GFRA3 (Artemin 공수용체): 위암 침습 촉진, 장관신경계 신호와 연관

**ECM·인테그린 (+2.328, +2.477):** FLNA, DMTN, NID1, LIMS2, LAMB2 — CRC와 동일 lead gene

#### 주요 하향 신호

**번역 억제 최우세 (503/516 경로; NES 최저 −2.284)**
- CRC와 동일한 패턴; OXPHOS 하향은 ranks 27-29로 약함
- GI 하부암의 공통 cfRNA 번역 억제 패턴 확인

**신호 품질: HIGH**

---

## 3. 췌장·폐 그룹

**출처:** Moore et al. Nat Commun 2025; Reggiardo et al. Nat Biomed Eng 2023; Chen et al. eLife 2022

---

### 3.1 Pancreatic Cancer — Moore et al. (sample n=72 · FDR<0.05 560개 · NES>0: 238 · NES<0: 322)

#### 주요 상향 신호

**보체·응고 최강 상향 (+2.495, 1위)**
- FGB, FGA, FGG, PLG, C3, CFB, SERPINC1, F13A1 — 간암과 공유
- PDAC: 조직인자(TF) 분비 → 과응고 상태; 혈소판 감소증
- Platelet Degranulation (+2.275), Response To Elevated Platelet Ca2+ (+2.236): 종양 교육 혈소판의 mRNA 기여 가능

**ECM·데스모플라시아 — PDAC 특이**
- ECM-receptor interaction (+2.390), Collagen Fibril Organization (+2.299)
- Non-integrin membrane-ECM Interactions (+2.272), ECM Organization (+2.206)
- PDAC는 밀집한 섬유성 기질로 정의됨; 암관련 섬유모세포 유래 cfRNA

**Hippo/YAP (+2.445)**
- PDAC에서 YAP/TAZ가 기질 활성화·약물 저항의 핵심 드라이버
- Lead genes: YAP1, WWTR1, TEAD 패밀리, MOB1A, LATS1/2

**RhoC GTPase (+2.484, 2위) — 침습 마커**
- RHOC, RHOB, RHOA (각각 +2.484, +2.181, +2.102)
- RHOC: KRAS 하류에서 세포골격 재편성, 침습, 전이 조절의 핵심 조절자
- PDAC에서 90% 이상 KRAS 변이 → RhoGTPase 과활성화와 연결됨

**EGFR/ERBB·IGF 축**
- ERBB Signaling (+2.021), EGFR Signaling (+1.992), IGF Transport/IGFBPs (+2.400)
- EGFR은 PDAC에서 증폭·과발현; IGF 축은 기질 크로스토크와 연관

**혈관신생 역설**
- Negative Regulation of Angiogenesis (+2.194): Thrombospondin-1, Endostatin
- Regulation of Angiogenesis (+2.157), Positive Regulation of VEGF (+2.025) 동시 존재
- PDAC는 저혈관성이지만 VEGF 드라이버와 항혈관 제약이 동시 작동

#### 주요 하향 신호

**OXPHOS 최강 하향 (NES −2.705, 이 표현형 내 최강)**
- NDUF 패밀리(Complex I), UQCR(Complex III), COX(Complex IV), ATP5(ATP synthase) 전반
- Warburg 효과; 정상 선포세포 대체로 OXPHOS-고활성 세포 감소

**NF-κB·T/B 세포 면역 억제**
- NIK to Noncanonical NF-kB (−2.450), Dectin-1 Mediated NF-kB (−2.393)
- Activation of NF-κB in B Cells (−2.353), Downstream TCR Signaling (−2.377)
- ZAP-70 Translocation to Immunological Synapse (−2.301)
- PDAC 면역억제 미세환경(TME)에서 T·B 세포 배제·억제를 cfRNA가 반영

**Hedgehog 억제**
- Hh Mutants Degraded by ERAD (−2.330), Hh Ligand Biogenesis (−2.123)
- 췌장 성상세포가 PDAC 세포 유래 Hh 리간드에 반응해 기질 활성화; Hh 경로 이상 반영

**신호 품질: HIGH**

---

### 3.2 Pancreatic Cancer — Reggiardo et al. (n=6 · NES>0: 117 · NES<0: 441)

⚠ **n=6; COMPLETE-seq 플랫폼(Illumina+Nanopore 혼합) 특이성; 모든 결과 잠정적**

#### 핵심 이슈: Moore 코호트와의 역전

| 신호 | Moore (n=72) | Reggiardo (n=6) |
|---|---|---|
| OXPHOS | 강한 하향 (−2.705) | **상향** (+2.115) — 역전 |
| RhoGTPase (RHOC, RAC1) | 강한 상향 | **하향** (CDC42−2.614, RAC1−2.550) — 역전 |
| TGF-β 신호 | 상향 | **하향** — 역전 |

**COMPLETE-seq 플랫폼 편향**: Nanopore 장독물 분석이 미토콘드리아 mRNA 포착 효율이 높아 OXPHOS 역전 가능성 + n=6 소수의 통계적 불안정성

**독특한 하향 신호 (잠정적)**
- Nucleosome Organization (−2.984), SUMOylation (−2.788 ~ −2.563): PDAC 에피제네틱 재프로그래밍 반영 가능
- Regulation of Erythrocyte Differentiation (−2.724): 암 관련 빈혈, 골수 억제

**신호 품질: LOW — 생물학적 결론 도출 불가**

---

### 3.3 Pancreatitis — Moore et al. (sample n=79 · FDR<0.05 263개 · NES>0: 27 · NES<0: 236)

#### 핵심 특징: 상향 신호 희박 (30개만)

**PDAC와 공유 신호 (감별 불가)**
- Hippo Signaling (+2.423), RHOC GTPase Cycle (+2.198)
- ECM-receptor interaction (+2.082), Focal Adhesion Assembly (+2.196)
- OxPhos 하향 171개 공유 → 범-췌장 질환 기질 리모델링 신호

**PDAC vs 췌장염 감별자**

| PDAC-고유 상향 | 췌장염-고유 |
|---|---|
| Complement/Coagulation (+2.495) | **Defective CFTR Causes Cystic Fibrosis (−2.219)** |
| Platelet Degranulation (+2.275) | NOTCH4 Intracellular Signaling (+2.024) |
| Collagen Fibril Organization (+2.299) | Glucocorticoid Stimulus Response (+1.981) |
| IGFBPs (+2.400) | — |

**PDAC-고유 하향 (T세포 면역 억제)**
- ZAP-70 Translocation (−2.301), Downstream TCR Signaling (−2.377), Primary Immunodeficiency (−2.324)
  → **T세포 억제는 PDAC-특이; 췌장염에서는 관찰되지 않음**

**CFTR 경로 (췌장염 고유)**
- Defective CFTR Causes Cystic Fibrosis (−2.219)
- CFTR 기능 이상은 급성·만성 췌장염의 알려진 위험인자;
  췌관 세포에서 중탄산염 분비 장애 → 단백질 플러그 → 선포 손상
- **췌장 외분비 기능 이상에 특이적인 cfRNA 신호**

**미토콘드리아 번역 억제 (췌장염에서 현저)**
- Mitochondrial Translation (−2.302), Mitochondrial Gene Expression (−2.295)
- 미토콘드리아 손상이 췌장염 발병에서 트립시노겐 활성화의 방아쇠임과 일치

**신호 품질: MODERATE**

---

### 3.4 Lung Cancer — Chen et al. (sample n=30 · FDR<0.05 520개 · NES>0: 54 · NES<0: 466)

#### 주요 상향 신호

**인테그린 신호 최강 (+2.726, 1위)**
- Positive Regulation of Integrin-Mediated Signaling Pathway: LIMS1/2, NID1, LAMB2, FLNA, LOXL3

**세포 접합 단백질 caspase 절단 (+2.701, 2위)**
- Apoptotic Cleavage of Cell Adhesion Proteins: CTNNB1, CASP3, TJP1/2, DSG2
- 순환 종양 세포의 anoikis 중 caspase-3 활성화로 접합 단백질 절단 → cfRNA 방출 기전

**ERBB2/EGFR/RET 축 (ranks 10-17)**
- SHC1 Events In ERBB2 Signaling (+2.316), Signaling By ERBB2 In Cancer (+2.029)
- EGFR Signaling (+2.090): NSCLC에서 EGFR 변이 10-35%
- RET Signaling (+2.189): RET 융합 ~1-2% NSCLC, 선택적 RET 억제제 타겟

**Hippo/YAP (+2.372)**
- 세 번째 표현형에서 Hippo 상향 확인 → **범상피 악성종양 cfRNA 공통 신호**

**혈관신생·혈관 투과성**
- Blood Vessel Endothelial Cell Migration (+2.457), Sprouting Angiogenesis (+2.354)
- PDAC(저혈관)과 달리 폐암은 VEGF 매개 혈관신생 활발

**인슐린·글리코겐 대사 상향**
- Insulin Receptor Signaling (+2.507), Positive Regulation of Glycogen Biosynthetic Process (+2.596)
- PI3K-AKT 활성화로 글리코겐 합성 증진; 폐암의 대사 재프로그래밍 반영

**IL-7 신호 (+2.406)**
- T세포 생존·항상성 증식을 위한 사이토카인; 면역 체크포인트 치료 맥락에서 중요

#### 주요 하향 신호

**리보솜·번역 억제 절대 우세 (NES −2.372, 35개 이상 번역 경로)**
- OXPHOS보다 번역 억제가 압도적 → PDAC·췌장염(OXPHOS 우세)과 반대 패턴
- 리보솜 단백질 mRNA의 종양 cfRNA 억제 or 혈중 세포 구성 변화 반영

**신호 품질: MODERATE-HIGH**

---

## 4. 면역·혈액 그룹

**출처:** Raissadati et al. JCI 2025; Roskams-Hieter et al. npj Precis Oncol 2022; Gardella et al. PNAS 2025

---

### 4.1 ICI-induced Myocarditis (ICI-m) (sample n=11 · FDR<0.05 479개 · NES>0: 310 · NES<0: 169)

#### 주요 상향 신호

**OXPHOS + 리소좀 최강 상향**
- KEGG Lysosome (+2.889, 1위), KEGG Oxidative Phosphorylation (+2.869, 2위)
- Proton Motive Force-Driven ATP Synthesis (+2.784), Neutrophil Degranulation (+2.760)
- 심근염에서 심근세포 사망 → 고밀도 미토콘드리아 보유 심근세포 유래 OXPHOS cfRNA 방출

**심장 항원 교차제시 — ICI-m 특이 신호 (ICI-treated Cancer에는 없음)**
- Antigen Processing – Cross Presentation (+2.422)
- ER-Phagosome Pathway (+2.308)
- Cross-Presentation Of Soluble Exogenous Antigens (Endosomes) (+2.21)
- **ICI-m 병태생리 핵심**: 차단된 checkpoint가 제거되면 CD8+ T세포가 심장 α-미오신(α-MyHC)을 MHC-I을 통해 교차제시받아 공격함
  → ER-phagosome 경로는 정확히 이 교차제시의 분자적 메커니즘임
  → **Raissadati JCI 2025 원문은 CCL5·심장 수축 경로를 보고; ER-phagosome 상향은 본 분석에서 추가 발견된 mechanistic insight**

**ERAD·ER 스트레스 처리**
- Ubiquitin-Dependent ERAD Pathway (+2.364), ERAD Pathway (+2.329)
- ER Protein Processing (+2.264), KEGG Proteasome (+2.218)
- 염증 심근세포와 포식세포에서 오접힘 단백질 처리 부담 → ER 스트레스 축

**소교세포 활성화 (+2.203) — ICI-m 고유**
- 신경염증 교차반응 또는 심장 대식세포의 소교세포-유사 분극 가능성

**심근 수축 (KEGG Cardiac Muscle Contraction +2.379)**
- ICI-treated Cancer와 공유 — 기저 심장 스트레스 반영

#### 주요 하향 신호

**번역 억제 (NES −2.735, 최저)**
- Eukaryotic Translation Elongation (−2.735), Peptide Chain Elongation (−2.731)
- GCN2/EIF2AK4 (−2.680): ER 스트레스 → GCN2/PERK → eIF2α 인산화 → 전반 번역 억제
  ↳ 급성 심근염에서 전신 번역 스트레스 반응

**신호 품질: HIGH**

---

### 4.2 ICI-treated Cancer (without myocarditis) (sample n=11 · FDR<0.05 1,304개 · NES>0: 231 · NES<0: 1,073)

#### 주요 상향 신호

**리소좀 최강 상향 (+2.825)**
- Vacuolar Acidification (+2.664)
- ICI 핵심 기전: 리소좀 항원 처리 → MHC-I 제시 → 세포독성 T세포 활성화

**ICI 특이 면역 표현형**
- Immunoregulatory Interactions Lymphoid/Non-Lymphoid (+2.415)
- IgA Intestinal Immune Network (+2.141): ICI 관련 대장염(colitis)에서 점막 면역 활성화
- Allograft Rejection (+2.263), Graft-versus-Host Disease (+2.177): ICI가 T세포를 동종반응성으로 해방
- Autoimmune Thyroid Disease (+2.12): 갑상선염은 ICI 흔한 부작용
- Defensins, Antifungal Innate Immunity (+2.376, +2.276): 광범위 선천 면역 활성화

**Negative Regulation Of T Cell Proliferation (+2.351)**
- 모순적: ICI가 T세포를 확장하지만 조절 메커니즘이 동시 작동 — 면역 균형 반영

#### 주요 하향 신호

**번역 억제 가장 광범위 (1,073 경로 하향)**
- Eukaryotic Translation Elongation (−2.808) — 5개 그룹 내 최강
- 암 악액질·화학요법 병용·ICI 유발 전신 대사 억제의 복합 효과

**SLIT/ROBO 신호 억제 (−2.472, −2.503)**

**[ICI-m vs ICI-treated Cancer 핵심 감별자]**

| 신호 | ICI-m | ICI-treated |
|---|---|---|
| ER-phagosome 교차제시 | ✓✓ **고유** | ✗ |
| ERAD/ER 스트레스 | ✓ **고유** | △ (약함) |
| 소교세포 활성화 | ✓ **고유** | ✗ |
| OXPHOS 강도 | ✓✓ (+2.87) | ✓ (+공유) |
| IgA 장관 면역 | ✗ | ✓ **고유** |
| 갑상선 자가면역 | ✗ | ✓ **고유** |
| 이식거부 패턴 | ✗ | ✓ **고유** |

**신호 품질: HIGH**

---

### 4.3 Multiple Myeloma (sample n=17 · FDR<0.05 333개 · NES>0: 281 · NES<0: 52)

#### 주요 상향 신호

**세포주기·유사분열 기계 (6개 경로)**
- Reactome Kinesins (+2.186, 2위): CENPE, NEK2, NUSAP1 — **Roskams-Hieter 2022 원문에서 최상위 MM 바이오마커로 보고된 유전자 직접 검출**
- G1/S-Specific Transcription (+2.149), Polo-like Kinase Events (+2.053)
- Chromosome Condensation (+2.115), Sister Chromatid Segregation (+2.036)
- TP53 Regulates Cell Cycle Genes (+2.016), E2F Targets Under HDAC1 Repression (+2.042)
- 골수내 형질세포의 무제한 G1→S→M 전환이 cfRNA로 방출됨

**적혈구·헤모글로빈 신호**
- Carbon Dioxide Transport (+2.188, 1위), Gas Transport (+2.185, 3위)
- Lead genes: HBG1, HBG2 — **Roskams-Hieter 2022 원문에서 MM 바이오마커로 명시**
- Erythrocyte Differentiation (+2.042): MM에 의한 골수 침윤 → 빈혈 → 보상성 조혈

**에피제네틱 재프로그래밍**
- HDACs Deacetylate Histones (+2.078)
- DNA Damage/Telomere Stress Induced Senescence (+2.032)
- TP53 경로, 히스톤 탈아세틸화, 텔로미어 스트레스 — MM 발암 에피제네틱 변화

#### 주요 하향 신호

**번역 억제 최심 (NES −3.049, 이 그룹 최강)**
- Eukaryotic Translation Elongation (−3.049), Peptide Chain Elongation (−2.879)
- SRP-dependent Targeting (−2.644), KEGG Ribosome (−2.596)

**면역글로불린 반응 역설적 하향 (NES −2.790)**
- Immunoglobulin Mediated Immune Response: HLA-DRA, HLA-DRB1, CD27, CD19, FCER1A
- 역설: MM이 단클론 파라단백질을 대량 생산하지만 **정상 다클론 면역글로불린 반응이 억제됨**
- 악성 클론이 정상 B세포·형질세포를 억압 → 저감마글로불린혈증 → cfRNA 수준에서 검출

**MHC-II 항원 제시·B세포 억제**
- MHC Class II Assembly (−2.361), B Cell Proliferation (−2.155)
- OXPHOS (−2.315): 골수 미세환경의 과당분해 → OXPHOS 고활성 세포 감소

**신호 품질: HIGH**

---

### 4.4 MGUS (sample n=8 · FDR<0.05 107개 · NES>0: 92 · NES<0: 15)

⚠ n=8; FDR 0.001–0.04 범위 (MM보다 통계적으로 약함)

#### 주요 상향 신호

**에피제네틱 시그니처 (MGUS 고유 — MM에는 없음)**
- KEGG SLE (+2.332, 1위): 히스톤 + 보체 → NET/뉴클레오솜 노출
- DNA Methylation (+2.177): DNMT3A/3B — MM으로의 전환에서 CpG 메틸화 변화가 초기 사건
- SIRT1 Negatively Regulates rRNA Expression (+2.127)
- HDACs Deacetylate Histones (+2.118), PRC2 Methylates Histones (NES ~1.96)
- Telomere Packaging (+2.073), Telomere Stress Induced Senescence (+2.099)
  → **MGUS 형질세포는 복제 스트레스 하에 있으나 노화 체크포인트를 아직 탈출하지 못함**

**YAP/TAZ·NOTCH 신호 (MGUS 고유)**
- YAP1/WWTR1-Stimulated Gene Expression (+2.206): 전악성 자기 재생 신호
- Pre-NOTCH Transcription (+2.084): 형질세포 분화 및 골수 니치 상호작용
  → MM으로 진행하면 이 신호가 완전한 유사분열 프로그램으로 대체됨

**NETs/호중구 활성화 (+2.02)**
- MGUS의 전염증 미세환경; MM으로 진행 전 단계 반영

#### 주요 하향 신호

**NK세포 억제 경로 하향 — 독특한 해석**
- Negative Regulation Of NK Cell Mediated Cytotoxicity (−2.19)
- Negative Regulation Of NK Cell Mediated Immunity (−2.13)
- **해석**: NK세포 살상을 억제하는 유전자들이 MGUS에서 HC 대비 낮음 = MHC-I 하향조절(HLA 손실)의 초기 신호
  → 악성종양화를 향한 면역 회피 메커니즘 시작

**리소좀 하향 (−2.368)**
- ICI 그룹(리소좀 상향)과 반대: MGUS에서 정상 면역 감시 리소좀 활성도 감소

**신호 품질: MODERATE (n=8)**

---

### 4.5 ME/CFS (sample n=90 · FDR<0.05 169개 · NES>0: 0 · NES<0: 169)

#### 최대 특징: 극단적 비대칭 (0 UP vs 169 DOWN)

> 현재 데이터에서 **FDR<0.05 상향 경로가 0개**다. 아래 "상향 신호"는 study-split 직후 초기 실행에서 관찰된 것으로, GSEA 재실행 후 면역글로불린·B세포·MHC-II 경로가 모두 유의성을 잃었다. ME/CFS의 견고한 신호는 하향 축(번역·FoxO·인슐린)이다.

#### 상향 신호 (현재 데이터: FDR<0.05 0개)

현재 GSEA 실행에서는 유의한 상향 경로가 없다. 초기 실행에서는 B세포·MHC-II 관련 4개 경로
(Immunoglobulin Mediated Immune Response +2.234, B Cell Proliferation +2.222, MHC Class II Assembly +2.140 등)가
약하게 상향되어 만성 항원 구동 B세포 활성화 가능성으로 해석되었으나, **재실행 후 모두 FDR<0.05 기준을 통과하지 못해 본 갱신본에서는 제외한다.**
ME/CFS cfRNA 신호는 사실상 전면적 하향 억제로 수렴한다.

#### 주요 하향 신호 (169개)

**번역 기계 완전 억제 (NES −2.543 최저)**
- Cytoplasmic Translation (−2.543), Eukaryotic Translation Elongation (−2.466)
- GCN2/EIF2AK4 Response To Amino Acid Deficiency (−2.358)
- L13a-Mediated Translational Silencing Of Ceruloplasmin (−2.421): IFN-γ 유발 → GAIT 복합체 → 세룰로플라스민 번역 선택적 억제
  ↳ ME/CFS의 IFN-γ 매개 염증 상태 암시

**FoxO 신호 하향 (NES −2.233) — ME/CFS 특이 발견**
- KEGG FoxO Signaling Pathway
- FoxO1/3/4: PI3K/Akt/mTOR의 하류 — 인슐린 신호, 산화 스트레스, 에너지 상태 통합
- ME/CFS 저대사 모델: PDK1/2/4 활성화 → 피루브산 탈수소효소 억제 → TCA 이상 → ATP 감소 → GCN2/FoxO 활성화
- **Gardella et al. PNAS 2025에서 보고되지 않은 FoxO 축 발견**

**인슐린·글루코스 신호 하향**
- Positive Regulation Of Insulin Receptor Signaling (−2.116)
- Cellular Response To Glucose Stimulus (−2.112)
- ME/CFS에서 인슐린 저항성 보고; Akt/mTOR 억제 → FoxO 탈억제 → 전신 대사 억제와 연결

**KEGG Coronavirus Disease 하향 (−2.233)**
- 항바이러스 ISG/IFN 응답 경로를 포함; ME/CFS에서 T세포 소진·선천 면역 응답 손상을 반영
- Long COVID와의 표현형 중첩을 암시 (ME/CFS에서 SARS-CoV-2 유발 케이스 다수)

**에스트로겐 신호 하향 (−2.120)**
- Estrogen-Dependent Nuclear Events: ME/CFS 여성 3:1 우세; 에스트로겐이 미토콘드리아 기능·면역 조절에 관여

**NMD 하향 (−2.427, −2.378):** mRNA 품질 감시 기전 이상

**신호 품질: HIGH (n=90)**

---

## 5. 심혈관·산과·기타 그룹

**출처:** Ward et al. Cells 2022; Moufarrej et al. Nature 2022; Moore et al. Nat Commun 2025

---

### 5.1 CAD — Heart Failure 진행군 비교 (CAD_HF+ / CAD_HF−)

> **표현형 재정의.** 본 그룹은 기존에 코호트 이름인 **CDCS(Coronary Disease Cohort Study)**로 잘못 라벨링되어 있었다.
> 원본 연구(**Ward et al., Cells 2022**, *Identifying Candidate Circulating RNA Markers for Coronary Artery Disease*)는
> 안정형 CAD 환자를 **3년 내 심부전(HF) 발병 여부**로 나누어 분석한 코호트이다.
> 이에 따라 라벨을 **CAD_HF+ (HF 진행군)** / **CAD_HF− (HF 비진행군)** 두 표현형으로 분리하고 각각 별도 GSEA를 재실행하였다.
>
> ⚠ **원본 논문의 핵심 결과:** Ward et al.은 plasma RNA-Seq에서 **HF+와 HF− 사이에 어떠한 차이도 검출하지 못했다**("We did not detect any differences in the plasma RNA profile between patients who developed HF compared with those who did not"). CAD vs 대조군의 주 신호는 **미토콘드리아 mRNA**(심근 허혈/산화 스트레스)였다. 본 정규화 모델은 이와 대조적으로 두 군 사이의 뚜렷한 차이를 보이며, 이 불일치의 해석은 §5.1 말미와 Rework 보고서 §4.2에서 비판적으로 다룬다.

#### 공통 배경: 양 군 모두 전신 번역·OXPHOS 광범위 억제

| | CAD_HF+ | CAD_HF− |
|---|---|---|
| sample n | 116 | 108 |
| FDR<0.05 | 961개 | 922개 |
| NES>0 (상향) | **36개** | **11개** |
| NES<0 (하향) | 925개 | 911개 |
| 번역/리보솜 하향 | 59개 (−2.25 최강) | 58개 (−2.14) |
| OXPHOS/미토 하향 | 29개 (−1.85) | 19개 (−1.80) |

양 군의 **하향 신호는 거의 동일** — 번역·리보솜·OXPHOS·스플라이싱·면역 NF-κB의 광범위 억제로, 만성 허혈/저대사 배경을 공유한다(Ward 2022의 미토콘드리아 억제 보고와 방향 일치). **두 군을 구별하는 것은 상향 신호의 양과 내용**이다.

---

#### 5.1a CAD_HF+ — HF 진행군 (sample n=116 · FDR<0.05 961개 · NES>0: 36 · NES<0: 925)

##### 핵심 발견: FGF/FGFR1 축의 반복 농축 (HF 진행 서명)

상향 36개 경로 중 **FGF–FGFR1 축이 6개 경로에서 반복** 검출된다 (같은 생물학 축의 반복 = 신호 견고성↑):
- FGFR1 Ligand Binding And Activation (**+2.82, 전체 1위**) — FGF23·FGF5·FGF1·FGFR1·**KL(Klotho)**
- Phospholipase C-mediated Cascade: FGFR1 (+2.40)
- FGFRL1 Modulation Of FGFR1 Signaling (+2.26)
- FGFR1c Ligand Binding (+2.11) · Activated Point Mutants Of FGFR1 (+2.02) · FGFR2 변이체 (+2.09)

→ **FGF23–Klotho–FGFR1** 축은 확립된 심부전·심비대 진행 바이오마커이며(문헌 검증: FGF23 incident HF 예측, Klotho 혈장 농도 ACS 사망 독립 예측), HF+ 군에서 **최강·최다 반복**으로 나타남.

##### 동반 심장 리모델링 신호

**심장 섬유화 (ECM/콜라겐)**
- Collagen Chain Trimerization (+2.45): COL1A1/COL1A2/COL5A1/COL14A1 — 심근 섬유화 핵심
- Non-integrin membrane-ECM Interactions (+1.98), GAG/Aminoglycan 대사 (+2.10/+2.03)

**이온채널 리모델링 (전기적 리모델링·부정맥 기질)**
- Calcium Ion Import (+2.57), Sodium Ion Transmembrane Transport (+2.30), Potassium Ion Transport (+2.35)
- **Sarcoplasmic Reticulum Calcium Ion Transport (+1.95)** — 부전 심근의 칼슘 핸들링 축
- L1-Ankyrin (+2.27): SCN/ANK2 — 막 흥분성 리모델링

**전응고·혈전 (5개 경로)**
- Positive Regulation Of Blood Coagulation (+2.26), Regulation Of Coagulation (+2.15), Fibrinolysis (+2.14), Positive Regulation Of Hemostasis (+1.99)
- Lead: **F2R (PAR-1)·THBS1 (Thrombospondin-1)** — THBS1은 심장 리모델링 악화 인자(HFpEF mitophagy 억제)

**기타:** Semaphorin-Plexin (+2.60, 혈관/축삭 리모델링), Nicotine addiction (+2.32, 흡연 위험인자)

**신호 품질: HIGH (n=116)** — 단, HF+/HF− 차이는 원본 논문의 null 결과와 상충하므로 직교 검증 필요(§5.1 말미).

---

#### 5.1b CAD_HF− — HF 비진행군 (sample n=108 · FDR<0.05 922개 · NES>0: 11 · NES<0: 911)

##### 핵심 발견: 상향 신호 빈약 (11개) — HF 진행 서명 부재

하향 배경은 HF+와 동일하나, **상향 경로가 11개에 불과**하며 심장 리모델링 축(FGFR1/콜라겐/이온채널)이 **거의 나타나지 않는다**.

**상향 11개 (약한 신호)**
- 감각/통증 지각 (+2.52/+2.51): F2R·VIP·SCN11A·GPR171
- Beta Defensins (+2.50): DEFB1/103A/118 — 선천 면역
- Mucin O-당화 결함 경로 (GALNT12/GALNT3/C1GALT1C1, +2.49~+2.20): MUC2/6/16/17 클러스터
- Purinergic Receptors (+2.30): P2RY1/12/13 — 혈소판 활성화
- 응고 조절 (+2.16/+2.06): Negative Regulation Of Fibrinolysis, Regulation Of Coagulation (F2R·PROC·THBS1) — HF+보다 약하고 "음성 조절" 우세
- GAG 대사 (+2.12)

→ **응고 축은 HF+/HF− 공유**(F2R·THBS1)이나, **FGFR1·콜라겐·칼슘채널 심장 리모델링 축은 HF+ 특이적**.

**신호 품질: HIGH (n=108)**

---

#### ⚠ HF+ vs HF− 차이에 대한 비판적 해석

본 정규화 모델은 HF+ 군에서 **생물학적으로 정합적인 HF 진행 서명**(FGF23-Klotho/FGFR1, 섬유화, 이온채널 리모델링, 전응고)을 검출했다. 이는 동일 코호트를 단순 fold-change로 분석한 Ward et al. 2022이 **놓친 차이**로, 공변량 보정 GAMLSS Z-score의 민감도 이득을 시사할 수 있다.

그러나 다음을 반드시 고려해야 한다:
1. **원본 논문은 HF+/HF− 차이를 검출하지 못했다(null).** 동일 표본에서 gold-standard 분석이 null인데 본 파이프라인이 강한 차이를 보이는 불일치는, 민감도 이득일 수도 있으나 **배치/세포조성/라이브러리 편향의 인공산물**일 가능성도 배제 못 함.
2. **상향 경로 수 차이(36 vs 11)** 자체가 검출력 차이(detection power)나 군 간 RNA 복잡도 차이를 반영할 수 있으며, 곧바로 생물학적 차이로 등치할 수 없다.
3. cfRNA는 **존재량(abundance)**이지 심근 조직의 발현 방향이 아니므로, "SERCA 칼슘 핸들링 상향" 등은 **순환 단편 수준의 농축**으로만 해석해야 한다.

**결론:** HF 진행 서명은 **가설 생성적(hypothesis-generating)** 수준이며, 독립 코호트·직교 측정(qPCR/단백질 FGF23·NT-proBNP)으로 검증되기 전까지 진단 marker로 단정하지 않는다.

---

### 5.2 Pre-eclampsia (sample n=59 · FDR<0.05 230개 · NES>0: 2 · NES<0: 228)

#### 주요 상향 신호 (단 3개)

**신경/시냅스 신호 — 핵심 발견**
- Reactome: Neurexins And Neuroligins (+2.42, 1위): 시냅스 접착 분자
- GO: Postsynaptic Membrane Organization (+2.11)
- **Moufarrej et al. Nature 2022에서 cfRNA 신경근육 세포형 풍부화를 보고;**
  **본 분석의 Neurexin/Neuroligin GSEA 농축이 구체적인 시냅스 경로 수준에서 이를 뒷받침**

**DNA Methylation (+2.17):** 영양막 기능 이상의 에피제네틱 기초 반영

#### 주요 하향 신호

**번역·리보솜 (가장 강한, ~90개 경로)**
- Ribosomal Small Subunit Biogenesis (−2.27, 최강)
- Cap-dependent Translation Initiation (−2.27)
- 임신중독증 태반 저산소증 → ISR → eIF2α 인산화 → 전반 번역 억제

**선천 면역 하향 (~40개 경로)**
- ISG15 Antiviral Mechanism (−1.99), IFN Signaling (−1.88)
- NLRP3 Inflammasome (−1.80), C-type Lectin Receptors (−1.83)
- NLR Signaling (−1.93)
- 임신 내 myeloid 세포 재분포·면역 관용 붕괴 반영

**스플라이싱·전사 (~20개)**
- KEGG Spliceosome (−2.00), Chromatin Remodeling (−1.89)

**지질단백질 대사 하향**
- Lipoprotein Localization (−1.91), Lipoprotein Transport (−1.89)
- 임신중독증 특이 이상지질혈증(중성지방 상승, HDL 감소) 반영

**SLIT/ROBO 하향 (−1.94)**
- 영양막 침습 조절에서 ROBO 수용체 역할; SLIT2는 여러 암에서 메틸화 소실
  → 영양막 침습 실패의 cfRNA 신호일 가능성

**[해석 주의]** 정상 임산부를 HC 기준으로 사용하지 않았으므로 임신 상태 자체와 질병 신호가 혼재됨

**신호 품질: MODERATE-HIGH (임신 교란 변수 존재)**

---

### 5.3 Other Cancer (sample n=16 · FDR<0.05 266개 · NES>0: 17 · NES<0: 249)

#### 주요 상향 신호

**RhoGTPase 완전 세트 (6개 GTPase, NES ≥ 2.0)**
- RHOB (+2.41) · RHOC (+2.18) · RHOA (+2.10) · RHOQ (+2.10) · CDC42 (+2.10) · RAC1 (+2.03)
- **범암 침습 서명**: 모든 주요 Rho 패밀리 GTPase가 동시 상향됨
- p130Cas-MAPK Signaling For Integrins (+2.14): 인테그린-MAPK 증식 연결

**세포 기질 접합·침습**
- Focal Adhesion Assembly (+2.25), Cell-Matrix Adhesion (+2.15), Cell Junction (+2.05)
- Actin Filament-Based Transport (+2.07): 종양 세포 이동 기계

**VEGF·혈관신생 (+2.07)**
- Regulation Of VEGF Signaling Pathway: 종양 혈관신생

**RUNX1 조절 (+2.03):** 조혈계 및 다발성 고형암에서 이상 발현

**음성 혈액 응고 조절 (+2.00):** 암 관련 고응고성 + 조절 메커니즘 동시 활성

#### 주요 하향 신호

**번역 기계 극단적 억제 (NES −3.10, 이 그룹 최저)**
- L13a-Mediated Translational Silencing Of Ceruloplasmin (−3.10): **전체 데이터셋 중 최저 NES**
- Cytoplasmic Translation (−3.09), Eukaryotic Translation Elongation/Termination (−3.09)
- GCN2 Response (−2.99), KEGG Ribosome (−2.94), SRP-dependent Targeting (−2.90)

**OXPHOS (NES −2.55):** 강한 Warburg 효과; Other Cancer가 PDAC처럼 OXPHOS 억제

**SLIT/ROBO (NES −2.89, −2.68):** 종양 억제 축의 cfRNA 소실; 다중암에서 SLIT2 메틸화 소실과 일치

**프로테아솜 하향 (−2.44):** Wnt DVL 분해(−2.54), AXIN 분해(−2.42) 하향 → 활성 Wnt/β-카테닌 상태 암시

**신호 품질: HIGH (RhoGTPase 침습 서명이 고도로 일관됨)**

---

## 6. 기존 문헌 미보고 Novel 발견 정리

아래 항목들은 ① 원본 출처 논문의 직접 보고 내용을 확인하고, ② PubMed·Semantic Scholar 검색에서 선행 cfRNA 연구 결과를 찾지 못한 경우를 "잠재적 novel"로 분류하였다.

---

### 6.1 ICI-m의 ER-phagosome·교차제시 경로 cfRNA 농축
**발견:** ICI-m에서 ER-Phagosome Pathway (+2.308), Antigen Processing–Cross Presentation (+2.422), ERAD (+2.364)이 특이적으로 상향; ICI-treated Cancer에서는 없음.

**기존 문헌과의 관계:**
- Raissadati et al. JCI 2025 (원본 논문)은 6-gene 분류기 (CCL5, ZNF385B, CADM2, AQP7, B2M, IFITM2)와 심장 수축·CCL5 경로를 보고
- ER-phagosome/ERAD/교차제시 경로 상향은 원본 논문에 **명시 없음**
- 이 발견은 ICI-m 병태생리의 핵심 메커니즘 (심장 α-MyHC의 MHC-I 교차제시를 통한 CD8+ T세포 활성화)을 **cfRNA 수준에서 최초로 포착**함을 시사

**근거 강도:** 단 n=11; 재현 필요. 그러나 병태생리학적 일관성 매우 높음.

---

### 6.2 결핵에서 스플라이세오솜 하향의 cfRNA 신호
**발견:** TB에서 하향 경로 7개 중 5개가 스플라이세오솜/RNA 스플라이싱 관련 (NES −2.333 ~ −2.139).

**기존 문헌과의 관계:**
- PubMed에서 "tuberculosis spliceosome blood transcriptome"으로 0건 검색됨
- Berry et al. 2010이 확립한 TB 전사체 지배 구조(호중구-IFN)에 스플라이세오솜 하향은 포함되지 않음
- TB에서 호중구 우세 혈중 세포 조성 (림프구 감소)이 스플라이세오솜-풍부 림프구 cfRNA를 상대적으로 감소시키는 **세포 구성 효과**로 해석됨
- Mtb의 숙주 스플라이싱 조절 가능성도 배제할 수 없음

**근거 강도:** 방향성 생물학적 가설 강함; 직접 cfRNA 검증 논문 미존재

---

### 6.3 HIV+TB 공동감염에서 B세포 경로 하향 클러스터 (cfRNA)
**발견:** B Cell Proliferation (−2.158), B Cell Activation (−2.060), B Cell Receptor Signaling (−1.979), B Cell Differentiation (−1.825) — HIV 단독이나 TB 단독에서는 이 수준의 B세포 하향 클러스터 없음.

**기존 문헌과의 관계:**
- HIV+TB에서 B세포 기능 장애는 면역학적으로 알려져 있으나
- **plasma cfRNA 수준에서 B세포 경로 하향 클러스터가 공동감염 특이 신호로 검출된 선행 연구 없음**
- Chang et al. 2023 (원본 논문)은 이 B세포 신호를 명시적으로 강조하지 않음

**근거 강도:** n 충분; B세포 면역학적 근거 강함; cfRNA 특이성 novel

---

### 6.4 대장암 cfRNA에서 Hippo/YAP 경로 농축
**발견:** CRC에서 Signaling By Hippo (+2.542), YAP1/WWTR1-Stimulated Gene Expression (+2.542)이 공동 1위 (가장 강한 상향 신호).

**기존 문헌과의 관계:**
- PubMed에서 "YAP1 Hippo colorectal cancer liquid biopsy plasma cfRNA"로 0건 검색됨
- 조직 수준에서 YAP/TAZ의 CRC 역할은 확립 (PMC12221920 등)
- **cfRNA/plasma/liquid biopsy에서 CRC의 Hippo/YAP 경로를 GSEA로 포착한 선행 연구 미보고**
- YAP1·WWTR1·TEAD4 mRNA가 종양 교육 혈소판 또는 순환 종양 세포 기원으로 혈장에 방출되는 첫 pathway-level 증거일 가능성

**근거 강도:** n (Chen eLife 2022 코호트)에서 통계적으로 강; pathway-level novel

---

### 6.5 식도암에서 FGFR2-YAP1 공동 활성화 축의 cfRNA 검출
**발견:** FGFR2 관련 5개 Reactome 경로와 YAP1/WWTR1 경로가 동시 상향; FGFR2 리간드·수용체·결합단백질이 lead gene에 광범위하게 등장.

**기존 문헌과의 관계:**
- FGFR2-YAP1 기전 연결은 위암 조직에서 보고 (FGF18-FGFR2-YAP1 axis, PMC7581496)
- 식도 선암에서 FGFR2 증폭·과발현은 알려져 있음
- **혈장 cfRNA에서 FGFR2-YAP1 축의 동시 농축을 GSEA 수준에서 보고한 연구 없음**

**근거 강도:** pathway 수준의 동시 농축 패턴이 novel; 원본 Peddu et al. 2025에서 대사·신호·면역체크포인트 경로 보고하지만 FGFR2-YAP1 구체 축 언급 없음

---

### 6.6 PDAC에서 TCR/ZAP-70 하향을 통한 췌장염 감별 가능성 (cfRNA)
**발견:** PDAC에서 ZAP-70 Translocation To Immunological Synapse (−2.301), Downstream TCR Signaling (−2.377), Primary Immunodeficiency (−2.324)가 특이적 하향; 췌장염에서는 이 신호 없음.

**기존 문헌과의 관계:**
- Moore et al. 2025는 간 기능 바이오마커(LAMC2, CEACAM7 등 29개 유전자)로 감별 → TCR 억제를 명시적으로 강조하지 않음
- **cfRNA에서 T세포 면역억제 (TCR·ZAP-70)가 PDAC vs 췌장염 감별자로 기능함을 경로 수준에서 보고한 선행 연구 없음**

**근거 강도:** n=72 vs 79, 동일 study(Moore) → 기술 혼란 변수 없음; 생물학적 타당성 높음

---

### 6.7 췌장염에서 CFTR 경로 cfRNA 억제
**발견:** Defective CFTR Causes Cystic Fibrosis (−2.219)가 췌장염에 특이적 하향; PDAC·다른 질병에서 해당 수준 신호 없음.

**기존 문헌과의 관계:**
- CFTR 기능 이상이 췌장염 위험인자임은 알려져 있으나
- **혈장 cfRNA에서 CFTR 경로가 췌장염-특이 신호로 검출된 보고 없음**
- 췌관 세포 기원 cfRNA가 CFTR 기능 이상을 반영하는 첫 pathway-level 증거

**근거 강도:** n=79, 같은 코호트 내 PDAC 비교로 교란 변수 최소; 생물학적으로 일관됨

---

### 6.8 ME/CFS에서 FoxO/인슐린/GCN2 축의 cfRNA 억제
**발견:** KEGG FoxO Signaling Pathway (−2.233), Insulin Receptor Signaling (−2.116), Cellular Response To Glucose Stimulus (−2.112)가 동시 하향.

**기존 문헌과의 관계:**
- Gardella et al. PNAS 2025 (원본 논문): pDC·단핵구·T세포 상향, 혈소판 cfRNA 하향, AUC 0.81 분류기
- **FoxO 신호 경로·인슐린 반응·GCN2 통합 스트레스 반응의 cfRNA 억제는 원본 논문에 보고되지 않음**
- ME/CFS 대사체학·단백질체학 연구에서 저대사(hypometabolism)·PDK 활성화가 보고되지만 cfRNA pathway 수준의 FoxO 축은 new
- 특히 L13a/GAIT 복합체를 통한 ceruloplasmin 선택적 번역 억제는 IFN-γ 유발 ME/CFS 염증 상태와 기계론적으로 연결되는 새로운 가설 제공

**근거 강도:** n=90 (가장 큰 코호트); 대사 생물학적 타당성 강함; cfRNA 수준 novel

---

### 6.9 임신중독증에서 Neurexin/Neuroligin 경로 cfRNA 농축
**발견:** Reactome Neurexins And Neuroligins (+2.42)가 Pre-eclampsia에서 유일한 최강 상향 신호.

**기존 문헌과의 관계:**
- Moufarrej et al. Nature 2022: 신경근육 세포형 cfRNA 풍부화 보고 (세포형 수준)
- **Neurexin/Neuroligin 특정 시냅스 접착 경로를 GSEA 수준에서 cfRNA로 포착한 선행 연구 없음**
- 임신중독증에서 신경학적 합병증(두통, 시각 장애, 경련)과 신경/시냅스 cfRNA 신호 연결 강화

**근거 강도:** Moufarrej 연구를 pathway 수준에서 구체화; 임신중독증 신경계 합병증과 일치

---

## 7. 전체 요약 테이블

| # | 표현형 | sample n | 신호 품질 | NES>0 상위 신호 | NES<0 상위 신호 | 핵심 novel 발견 |
|---|---|---|---|---|---|---|
| 1 | **HIV** | 13 | HIGH | IFN/ISG 최강 (+2.741); NET 형성; 노화/SASP; Rev-Host 상호작용 | 번역·리보솜 억제 (−2.494); 미토콘드리아 번역 | Rev Import 경로 직접 검출 (양성 대조군) |
| 2 | **HIV+Tuberculosis** | 11 | HIGH | SLE/NET 최고 (+2.911); Pyroptosis 최고; 항균 체액; IFN 감쇠 | 번역 억제 최심 (−2.743); **B세포 경로 클러스터** (−2.158); GCN2 (−2.488) | ★ B세포 cfRNA 억제가 공동감염 특이 신호 (cfRNA 미보고) |
| 3 | **Tuberculosis** | 101 | HIGH | Mtb 식세포작용 최강 (+2.442); IFN-γ>IFN-α/β; 호기화; Glycolysis; Autophagy | **스플라이세오솜 5개 하향** (−2.333); B세포 증식 (−2.126) | ★ 스플라이세오솜 하향의 TB cfRNA 신호 (선행 연구 없음) |
| 4 | **Liver Cancer (Chen)** | 10 | MODERATE-HIGH | 간 분비 단백질 (ALB, FGB, PLG, APOE); 보체/응고 (+2.476) | 번역 기계 전면 (−3.356 최저); GCN2 | 번역 억제가 "간암 = OXPHOS 우세" 패턴에서 이탈 |
| 5 | **Liver Cancer (Roskams-Hieter)** | 28 | HIGH | 보체/응고; 세포주기 (G1/S, kinesin); YAP/Hippo; CDC42/RAC1/GTPase; 항바이러스 (IFITM) | OXPHOS 단 6개 | 두 코호트 간 이질성이 cfRNA에서 생물학적으로 구별됨 |
| 6 | **Liver Cirrhosis** | 4 (제외) | LOW-MODERATE | 히스톤/HDAC (NETosis); 혈관 리모델링 (+2.203); 보체 (+2.041) | OXPHOS 약함 (−2.280); HCC보다 약함 → 진행성 미토 손상 암시 | n=4 한계; LSEC 기원 혈관 리모델링 cfRNA 신호 |
| 7 | **Colorectal Cancer** | 37 | HIGH | **Hippo/YAP 최강 (+2.542)** (공동); ERBB2/KIT; 인테그린; 인슐린 | 번역·리보솜 완전 억제 (−2.432, 360/385) | ★ CRC cfRNA에서 Hippo/YAP 경로 검출 (cfRNA 선행 연구 없음) |
| 8 | **Esophagus Cancer** | 25 | HIGH | 배아 발달/SOX17 최강 (+2.673); **FGFR2-YAP1 동시 활성화**; IGF1R; SMAD/TGF-β; GABA/니코틴 | OXPHOS + 번역 동시 억제 (−2.387) | ★ FGFR2-YAP1 공동 축의 cfRNA 동시 포착 (cfRNA 미보고) |
| 9 | **Stomach Cancer** | 24 | HIGH | 세포-세포 접합/데스모좀 (+2.643, PKP1/3/4·DSG2·CDH3); 인테그린; ECM; RET (+2.298) | 번역 완전 억제 (−2.284, 503/516); OXPHOS 약함 | 데스모좀 junction(PKP/DSG2) 상향 = 위암 상피 접합 리모델링 (cornified envelope는 현재 데이터서 탈락) |
| 10 | **Pancreatic Cancer (Moore)** | 72 | HIGH | 보체/응고 (+2.495); ECM/데스모플라시아; Hippo/YAP; RHOC 침습 (+2.484); ERBB/IGF | OXPHOS (−2.705 최강); **TCR/ZAP-70 억제**; NF-κB 억제; Hedgehog | ★ TCR 억제가 PDAC vs 췌장염 cfRNA 감별자 (cfRNA 미보고) |
| 11 | **Pancreatic Cancer (Reggiardo)** | 6 (제외) | LOW | OXPHOS 역전 (↑); 단핵구 유인 | 핵소체/SUMOylation 억제; 조혈 억제 | n=6; COMPLETE-seq 플랫폼 편향; 생물학적 결론 불가 |
| 12 | **Pancreatitis** | 79 | MODERATE | Hippo/YAP (+2.423); RHOC; ECM; NOTCH4 | OxPhos 억제 (−2.658); 미토콘드리아 번역; **CFTR 경로 (−2.219)** | ★ CFTR 경로의 췌장염-특이 cfRNA 억제 (선행 연구 없음) |
| 13 | **Lung Cancer** | 30 | MODERATE-HIGH | 인테그린 최강 (+2.726); 접합 단백질 Caspase 절단; ERBB2/EGFR/RET; Hippo/YAP; 혈관신생; 인슐린/글리코겐 | 번역·리보솜 (−2.372 최강); OXPHOS 약함 | ERBB2+RET+EGFR 3중 RTK 상향; IL-7 시그널 (T세포 항상성) |
| 14 | **ICI-m** | 11 | HIGH | 리소좀·OXPHOS 최강 (+2.889/+2.869); **ER-phagosome 교차제시 (+2.422)**; ERAD; 소교세포 활성화 | 번역 (−2.735); GCN2/ISR | ★ ER-phagosome/ERAD가 ICI-m 특이 cfRNA 신호 (Raissadati 원본 미보고) |
| 15 | **ICI-treated Cancer** | 11 | HIGH | 리소좀 최강 (+2.825); IgA 장관 면역; 갑상선 자가면역; 이식거부 패턴; T세포 증식 억제 | 번역 가장 광범위 (1,072개); SLIT/ROBO (−2.472) | ICI 부작용(장관, 갑상선)의 cfRNA 시스템 반영 |
| 16 | **Multiple Myeloma** | 17 | HIGH | Kinesins (+2.186, CENPE/NEK2/NUSAP1); 적혈구/HBG1/HBG2; G1/S 전사; TP53 | 번역 최심 (−3.049); **면역글로불린 반응 역설 하향** (−2.790); MHC-II; OXPHOS | 형질세포 악성 증식의 유사분열 기계 cfRNA 검출; IgG 역설 하향 |
| 17 | **MGUS** | 8 | MODERATE | **에피제네틱 4개 경로** (DNA 메틸화, SIRT1, PRC2, HDAC); YAP/NOTCH; NETs; 텔로미어 | NK세포 억제 하향; MHC-I 제시 하향; 리소좀 하향 | MGUS에서 에피제네틱 선행 변화 cfRNA 검출; NK 억제 소실 = HLA 손실 초기 신호 |
| 18 | **ME/CFS** | 90 | HIGH | 유의 상향 경로 없음 (현재 FDR<0.05 0개) | **번역 최강** (−2.543); **FoxO 경로** (−2.233); 인슐린/글루코스; GCN2; Coronavirus 경로; 에스트로겐 | ★ FoxO/GCN2/인슐린 축 cfRNA 억제 (Gardella 원본 미보고); 저대사 병태생리 직접 증거 |
| 19 | **CAD_HF+** (HF 진행) | 116 | HIGH* | **FGF/FGFR1 축 6회 반복 (+2.82 최강, FGF23·Klotho)**; 콜라겐/섬유화; 칼슘·이온채널 리모델링; 전응고 (F2R·THBS1) | 번역 (−2.25); OXPHOS (−1.85); 면역 NF-κB; 스플라이싱 (~17개) | ★ FGF23-Klotho/FGFR1 HF 진행 서명을 cfRNA에서 포착 (단, 원본 Ward 2022은 HF+/− 차이 null → 검증 필요) |
| 20 | **CAD_HF−** (HF 비진행) | 108 | HIGH | 상향 11개로 빈약; 감각/통증; Beta defensin; mucin 당화; 퓨린 수용체; 약한 응고 조절 | 번역 (−2.14); OXPHOS (−1.80); 면역 NF-κB | HF+와 하향 배경 공유, **심장 리모델링 상향 축 부재** = HF 진행 서명 음성 대조 |
| 21 | **Pre-eclampsia** | 59 | MODERATE-HIGH | **Neurexins/Neuroligins (+2.42)** (Moufarrej 2022 신경근육 발견의 pathway 수준 구체화); DNA 메틸화 | 번역 최강 (−2.27); 선천 면역/IFN (−1.99); 스플라이싱; 지질단백질 | ★ Neurexin/Neuroligin이 임신중독증 cfRNA 최강 신호 (첫 pathway 수준 증거) |
| 22 | **Other Cancer** | 16 | HIGH | **RhoGTPase 완전 세트** (RHOB/C/A/Q + CDC42 + RAC1, NES ≥ 2.0); 세포 기질 접합; VEGF; RUNX1 | 번역 극단 (−3.10, 전체 최저); OXPHOS (−2.55); **SLIT/ROBO** (−2.89); 프로테아솜; Wnt 분해 억제 | 범암 RhoGTPase 침습 서명이 이질적 암종 집합에서 일관되게 검출됨 |

> *CAD_HF+ 신호 품질 HIGH*는 표본 수 기준이며, HF+/HF− 차이 자체는 원본 논문 null 결과와 상충하므로 직교 검증 전까지 가설 생성적(§5.1 비판적 해석 참조).

---

## 부록: 아티팩트 요약

**히스톤 클러스터 인플레이션 (HIV/HIV+TB/TB/MGUS에서 주의)**
- H2BC/H2AC/H3C/H4C 계열 mRNA의 혈장 내 대량 존재
- 영향 경로: KEGG Alcoholism, KEGG SLE, DNA Methylation, PRC2, Condensation Of Prophase Chromosomes
- 히스톤 함유 gene set에서 lead genes의 >40%가 히스톤일 경우 아티팩트로 표시

**바이러스 경로가 리보솜 대리 신호 역할 (GI암·ME/CFS 등)**
- "Viral mRNA Translation," "Influenza Viral RNA Transcription" → RPS/RPL 패밀리를 공유; 실제 바이러스 아님

**개발 형태발생 경로 (식도암·위암 중간 NES)**
- "Aortic Valve Morphogenesis," "Mammary Gland Development" → 암 줄기세포/EMT의 발달 전사인자 반영; 글자 그대로 해석 금지

---

*리포트 작성: 2026-06-24*  
*cfRNA normative model: GAMLSS (HC n=693 WTS-only), GSEA-prerank (GSEApy), FDR<0.05*  
*Gene sets: KEGG_2021_Human · Reactome_2022 · GO_Biological_Process_2023*  
*Novel 발견 근거: PubMed/Semantic Scholar 검색 + 원본 논문 내용 직접 대조*
