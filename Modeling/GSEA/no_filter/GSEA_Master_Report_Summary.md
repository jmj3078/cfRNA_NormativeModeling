# GSEA Master Report — 분류별 요약

**기준:** GAMLSS normative model (HC n=693, WTS-only) · 21개 표현형 · KEGG 2021/Reactome 2022/GO BP 2023 · FDR q<0.05  
**Novel 발견 9건 확인** · 번역·리보솜 하향 및 OXPHOS 하향은 21/22개 표현형 공통 cross-disease confound — 질병 특이 해석 보류  
**제외:** Liver Cirrhosis (n=4), Pancreatic Cancer–Reggiardo (n=6, COMPLETE-seq 플랫폼 편향)

---

## Tier 1 — 보고할 가치가 매우 높음 (11개)

> 기존 signature의 고충실도 재현, 또는 cfRNA에서 최초로 포착된 Novel pathway, 또는 두 질병 상태의 차이를 정규화 모델이 명확히 분리한 케이스

---

### 1. HIV (n=13 · FDR 482개 · ↑434 ↓48) — HIGH

**핵심 재현 신호**
- IFN/ISG 축 4개 연속 1–4위: Antiviral Mechanism By ISGs (+2.741) → IFN α/β (+2.595) → ISG15 (+2.585) → Interferon Signaling (+2.541). Lead: IFI27, OAS1/2/3, IFIT1/2/3, MX1/2, ISG15, RSAD2, IRF7, ADAR, BST2, IFITM1/2/3. 만성 HIV ART 상태의 잔존 type I IFN tone을 완벽히 재현.
- NET 형성 (+2.371): MPO, FCGR1A/2A/3B, 히스톤(H2AC/H2BC/H3C/H4C 군), PIK3CB, MAPK14.
- 세포 노화·SASP (+2.346/+2.239): p21/p16 매개 면역 조기 노화.
- RUNX1/거핵구 (+2.359): RUNX1, GATA2, MYB, NFE2 — HIV 혈소판 감소증 기전 반영.

**Novel / 특이 발견**
- Interactions Of Rev With Host Cellular Proteins (+1.692), Nuclear Import Of Rev Protein (+1.734): importin, NUP98/153 → **HIV 복제 분자 생물학을 cfRNA 수준에서 직접 포착한 양성 대조군**. 이 두 경로는 비감염 대조군에서 절대 나타나지 않는 HIV-specific 경로.

**기대했으나 미검출 신호**
- CXCR4/CCR5 공수용체 신호 경로 (HIV 세포 진입 핵심): 수용체 자체는 리간드 결합 후 내재화되므로 cfRNA 상향이 예상되었으나 별도 GSEA 경로로 농축되지 않음.
- NLRP3 인플라마솜/Pyroptosis 경로 (CD4+ T세포 pyroptosis — HIV CD4 감소의 주요 기전 중 하나): HIV 단독에서는 weak, HIV+TB 공동감염에서만 NES +2.575로 현저. 단독 감염에서는 신호 강도가 FDR 기준 미달.
- Tat 매개 전사 조절 경로: Rev 경로는 검출되었으나 HIV Tat에 의한 Cap-dependent translation 억제 경로는 별도 식별되지 않음 (번역 하향 신호에 confound로 혼재).

**아티팩트 주의:** KEGG Alcoholism, SLE, DNA Methylation, PRC2 — 히스톤 클러스터 cfRNA 과다로 인한 팽창.

---

### 2. Tuberculosis (n=101 · FDR 606개 · ↑599 ↓7) — HIGH

**핵심 재현 신호**
- Berry et al. 2010 (Nature Immunol.) 호중구-IFN 전사체 지배 구조 완벽 재현: NET (+2.496), SLE/히스톤 (+2.477).
- **Mtb 특이 숙주-병원체 경로 직접 검출**: Regulation Of Actin Dynamics For Phagocytic Cup Formation (+2.442), Infection With Mycobacterium Tuberculosis (+2.149), Response Of Mtb To Phagocytosis (+2.117). Lead: ARPC1B/2, FCGR1A/2A, LIMK1, WAS, WASF1/2, GRB2, ELMO1/2, DOCK180.
- Suppression Of Phagosomal Maturation (+1.870): Mtb가 리소솜 융합을 억제해 세포 내 생존하는 기전이 cfRNA에서 포착.
- ISG 지문 병원체-특이적: GBP1/2/4/5/6 (항Mtb ISG, 파고솜 동원) 우세 vs HIV의 ISG15/OAS1 우세.
- Glycolysis (+2.273): HIF-1α → GAPDH, PKM, LDHA, HK1, ENO1 — 대식세포 Warburg 전환.
- Autophagy 경로 동시 농축: Regulation of Autophagy Of Mitochondrion (+1.738), Regulation of Autophagy (+1.710).

**Novel 발견**
- **스플라이세오솜 하향 클러스터 (TB-고유)**: RNA Splicing Via Transesterification Reactions (−2.333), mRNA Processing (−2.289), KEGG Spliceosome (−2.162). Lead: PNN, YTHDC1, SRRM2, HNRNPU/K/M/H1, SRSF5/8/10, DDX5/46, SF3B1. PubMed "tuberculosis spliceosome blood transcriptome" 검색 결과 0건 — cfRNA 미보고 신호. 호중구 우세 혈액 세포 조성으로 스플라이세오솜-풍부 림프구 cfRNA 상대적 감소 해석.

**기대했으나 미검출 신호**
- IL-12/IL-23 → Th1 분화 경로 (Mtb 면역의 핵심): IFN-γ 상향은 포착되었으나 Th1 세포 분극 경로(IL12 signaling, STAT4, T-bet) 자체가 별도 top-ranking 경로로 나타나지 않음.
- TNF-α 신호 경로 (대식세포 활성화·육아종 형성의 주요 사이토카인): TNF 관련 경로가 상위 경로 목록에서 명시적으로 언급되지 않음. 대식세포 활성화는 식세포 경로로 간접 반영.
- NOD2/PRR Innate Recognition 경로 (Mtb 세포벽 성분 인식): Dectin-1/CLR 경로가 PDAC 하향에서는 등장하나 TB 상향에서는 주요 경로로 식별되지 않음.

---

### 3. HIV+Tuberculosis (n=11 · FDR 397개 · ↑306 ↓91) — HIGH

**핵심 재현 신호**
- KEGG SLE/NET 공동감염 최고 NES (+2.911): 두 병원체 NET 유도 상승 효과.
- **Pyroptosis (+2.575) — 공동감염 수렴**: HIV(CD4 T세포) + Mtb(대식세포) NLRP3 경로 동시 작동. HIV 단독(NES 약함)보다 현저히 높음.
- 항균 체액성 반응 (+2.367): CTSG, ELANE, AZU1, MPO, defensins — 호중구 과립 방출.
- IFN α/β 신호 감쇠 (+2.317 vs HIV 단독 +2.595): Mtb의 SOCS1/IL-10 매개 type I IFN 억제 반영.

**Novel 발견**
- **B세포 경로 하향 클러스터 (공동감염 특이)**: B Cell Proliferation (−2.158), B Cell Activation (−2.060), B Cell Receptor Signaling (−1.979), B Cell Differentiation (−1.825). Lead: CD79A, CD79B, MS4A1(CD20), CD19, CD22, PAX5, BTK, BLNK, BCL2, MEF2C. HIV 단독·TB 단독에서는 이 수준의 B세포 하향 클러스터 없음. Chang et al. 2023 원본 미보고. **cfRNA 수준에서 공동감염 특이 B세포 면역 전면 손상 최초 포착**.
- GCN2/EIF2AK4 (−2.488): 이중 병원체 아미노산 고갈 → 대사 스트레스 정량화.

**기대했으나 미검출 신호**
- Mtb-특이 식세포/파고솜 경로: TB 단독(+2.442)에서 최강으로 나타났으나 HIV+TB에서는 HIV 면역 억제로 대식세포 기능이 훼손되어 해당 신호가 희석됨 — 공동감염에서의 Mtb 식세포 신호 약화 자체가 생물학적으로 의미 있지만 별도 언급은 없음.
- BAFF/APRIL 매개 B세포 생존 경로의 직접 억제 경로: B세포 기능 하향은 포착되었으나 그 기전 경로(BAFF/APRIL 결핍이 class switching에 필수인 Tfh 고갈 경유)까지는 경로 수준에서 분해되지 않음.

---

### 4. ICI-induced Myocarditis (n=11 · FDR 479개 · ↑310 ↓169) — HIGH

**핵심 재현 신호**
- KEGG Lysosome (+2.889, 1위), KEGG OXPHOS (+2.869, 2위): 심근세포 손상·사멸 → 고밀도 미토콘드리아 보유 심근세포 유래 cfRNA. OXPHOS 상향은 ICI-treated Cancer와 달리 ICI-m에서 극도로 강함.
- Neutrophil Degranulation (+2.760), ERAD (+2.364), Cardiac Muscle Contraction (+2.379).

**Novel 발견**
- **ER-phagosome 교차제시 축 (ICI-m 고유)**: Antigen Processing–Cross Presentation (+2.422), ER-Phagosome Pathway (+2.308), Cross-Presentation Of Soluble Exogenous Antigens (+2.21). ICI-m의 핵심 병태생리(차단된 체크포인트 해제 → CD8+ T세포가 심장 α-MyHC를 MHC-I 교차제시로 공격)가 cfRNA에서 기계론적으로 포착. Raissadati JCI 2025 원본(6-gene 분류기, CCL5·심장 수축 경로 보고)에 명시 없는 mechanistic insight. ICI-treated Cancer에서는 이 신호 완전 부재 → **ICI-m 특이 cfRNA 바이오마커 후보**.
- Microglia-Like Activation (+2.203): 심장 대식세포의 소교세포-유사 분극 (ICI-treated에 없음).
- GCN2/EIF2AK4 (−2.680): ER 스트레스 → PERK/GCN2 → eIF2α 인산화 연결 — disease-specific ISR.

**기대했으나 미검출 신호**
- 심장-특이 자가항원(α-MyHC, TNNI3) 관련 adaptive immune activation 경로: 교차제시 경로는 포착되었으나 심장 자가항원에 특이적인 CD8 클론 팽창 경로(TCR activation specific to cardiac antigen)는 식별 불가.
- TNF-α/IL-6 사이토카인 폭풍 경로: 임상 ICI-m에서 관찰되는 염증 사이토카인 패턴이 cfRNA GSEA 상위 경로에서 별도 농축되지 않음.
- CCL5 경로: Raissadati 2025 원본이 핵심 분류기로 보고했으나 GSEA prerank 상위 경로로 부상하지 않음 (개별 유전자 수준 신호 vs 경로 수준 신호의 차이).

---

### 5. Multiple Myeloma (n=17 · FDR 333개 · ↑281 ↓52) — HIGH

**핵심 재현 신호**
- **Roskams-Hieter 2022 원문 바이오마커 직접 재현**: Reactome Kinesins (+2.186) — CENPE, NEK2, NUSAP1이 원문 MM 최상위 바이오마커로 보고된 유전자와 정확히 일치.
- CO₂/Gas Transport (+2.188/+2.185): HBG1, HBG2 — 원문에서 MM 바이오마커로 명시된 태아 헤모글로빈 유전자. MM 골수 침윤 → 빈혈 → 보상성 조혈 반영.
- G1/S-Specific Transcription (+2.149), Polo-like Kinase Events (+2.053), Sister Chromatid Segregation (+2.036), TP53 Regulates Cell Cycle Genes (+2.016): 골수내 형질세포 무제한 증식 기계 cfRNA.
- HDACs Deacetylate Histones (+2.078), Telomere Stress Induced Senescence (+2.032): MM 에피제네틱 변화.

**Novel / 특이 발견**
- **면역글로불린 반응 역설적 하향 (−2.790)**: HLA-DRA, HLA-DRB1, CD27, CD19, FCER1A — 악성 형질세포 클론이 정상 B세포·형질세포 억압 → 저감마글로불린혈증이 cfRNA 수준에서 검출. 단클론 파라단백질 과다 생산과 정상 다클론 Ig 억제가 동시에 cfRNA에서 포착.
- MHC-II 하향 (−2.361), OXPHOS (−2.315): 골수 미세환경 면역 억제 반영.

**기대했으나 미검출 신호**
- 프로테아솜 경로 상향: 보르테조밉(bortezomib) 등 프로테아솜 억제제가 MM 1차 치료임 → 프로테아솜 과활성이 cfRNA 상향으로 나타날 것으로 예상되었으나 확인되지 않음.
- BCMA(TNFRSF17)/SLAMF7 신호 경로: MM 면역치료의 핵심 타겟 → cfRNA에서 해당 신호 검출 실패.
- IRF4/MYC/XBP1 형질세포 정체성 전사인자 경로: MM 발암에서 핵심이나 GSEA 상위 경로로 식별되지 않음 (개별 유전자 발현은 있을 수 있으나 경로 수준 농축 미확인).
- 단클론 면역글로불린 생산 경로 (IgG/IgA 중쇄 생합성): 역설적으로 생산 경로가 상향 농축되지 않음 — 단클론 과발현이 다클론 억제에 의한 net cfRNA 신호로 나타남.

---

### 6. ME/CFS (n=90 · FDR 169개 · ↑0 ↓169) — HIGH

**핵심 특징: 극단적 비대칭 (FDR<0.05 상향 경로 0개)**

**Novel 발견 — 3가지 동시**
- **FoxO 신호 경로 하향 (−2.233)**: KEGG FoxO Signaling Pathway — FoxO1/3/4: PI3K/Akt/mTOR 하류 조절자. ME/CFS 저대사 모델(PDK1/2/4 → 피루브산 탈수소효소 억제 → TCA 이상 → ATP 감소 → GCN2/FoxO 활성화)과 직접 연결. Gardella et al. PNAS 2025 미보고.
- **GCN2/EIF2AK4 (−2.358) + L13a/GAIT (−2.421)**: GCN2는 ME/CFS 저대사 기전, L13a-mediated translational silencing은 IFN-γ 유발 ME/CFS 염증 특이 성분 → cross-disease confound와 분리 가능한 disease-specific 성분.
- **Insulin Receptor Signaling (−2.116), Cellular Response To Glucose Stimulus (−2.112)**: FoxO → 인슐린 저항성 → 전신 대사 억제 연결.
- KEGG Coronavirus Disease (−2.233): Long COVID-ME/CFS 표현형 중첩 암시.
- Estrogen-Dependent Nuclear Events (−2.120): ME/CFS 여성 3:1 우세와 일치.
- NMD 하향 (−2.427, −2.378): mRNA 품질 감시 기전 이상.

**기대했으나 미검출 신호**
- NK세포 기능 이상 경로: NK세포 수·기능 저하가 ME/CFS의 핵심 면역 이상으로 널리 보고되나 GSEA에서 NK-specific 경로 농축 미확인.
- 미토콘드리아 복합체 I/III 직접 억제 경로: 대사 이상의 근원으로 예상되었으나 PDK1/2/4 경로 자체보다는 FoxO·GCN2 경로가 검출됨 (상위 조절 경로 수준의 신호).
- HPA 축 (시상하부-뇌하수체-부신) 이상 경로: 코티솔 조절 이상이 ME/CFS에서 일관되게 보고되나 glucocorticoid signaling 경로가 상위에 농축되지 않음.
- 상향 신호의 완전 부재: 만성 항원 구동 B세포 활성화 경로 (초기 실행에서 약하게 나타났으나 재실행 후 FDR<0.05 탈락) — 재현성 부족으로 해석 불가.

---

### 7. Colorectal Cancer (n=37 · FDR 385개 · ↑25 ↓360) — HIGH

**핵심 재현 신호**
- ERBB2/KIT/EGFR RTK 신호 (ranks 8-17): GRB2, SRC, HSP90AA1, SOS1, JAK2, STAT1/5B, PIK3R1, EGF.
- Integrin 신호 (+2.488): LIMS1/2, NID1, LAMB1/2, LAMC1, FLNA, LOXL3.
- Insulin/Glucose (+2.410, +2.052): Warburg 효과 연계 PI3K-AKT.

**Novel 발견**
- **Hippo/YAP 공동 1위 (+2.542)**: Signaling By Hippo (+2.542) = YAP1/WWTR1-Stimulated Gene Expression (+2.542). Lead: TJP1, STK24, YAP1, WWTR1, TJP2, AMOT, MOB1A, LATS2, AMOTL2, HIPK2, KAT2B, TEAD4. 조직 수준에서 YAP/TAZ CRC 역할은 확립되어 있으나 **plasma cfRNA/liquid biopsy에서 CRC의 Hippo/YAP 경로를 GSEA 수준에서 포착한 선행 연구 없음**. YAP1/WWTR1/TEAD4 mRNA 종양 교육 혈소판 또는 CTC 기원.

**기대했으나 미검출 신호**
- **KRAS/MAPK/ERK 경로**: CRC의 ~40-45%에서 KRAS 변이. KRAS 하류 MAPK/ERK 신호 경로가 GSEA 상위 경로에서 독립적으로 식별되지 않음 — ERBB2/EGFR 신호와 혼재 또는 cfRNA 수준의 신호 미약.
- **Wnt/APC/β-카테닌 경로**: APC 변이가 CRC의 ~80%에서 발생, CRC 발암의 가장 초기·흔한 사건. DKK1, AXIN2, TCF4(TCF7L2) 등 Wnt 표적 유전자가 상위 경로로 농축되지 않음 — 가장 중요한 미검출 드라이버.
- MMR/MSI 신호 경로: MSI-H CRC (~15%)에서 neoantigen 부하 증가·면역 활성화가 예상되나 별도 경로로 식별 안됨.
- TP53 직접 하류 경로: TP53 변이 ~50% CRC. Roskams-Hieter에서는 TP53/Oncogene Induced Senescence 경로가 상향되었으나 CRC에서는 뚜렷하지 않음.

---

### 8. Pancreatic Cancer — Moore (n=72 · FDR 560개 · ↑238 ↓322) — HIGH

**핵심 재현 신호**
- Complement/Coagulation (+2.495, 1위): FGB, FGA, FGG, PLG, C3, CFB, SERPINC1 — 과응고 상태.
- Platelet Degranulation (+2.275): 종양 교육 혈소판 cfRNA.
- ECM/Desmoplasia (4개 경로): ECM-receptor (+2.390), Collagen Fibril (+2.299), Non-integrin ECM (+2.272), ECM Organization (+2.206) — PDAC 밀집 섬유성 기질.
- RHOC GTPase (+2.484, 2위), RHOB (+2.181), RHOA (+2.102): KRAS 하류 세포골격·침습.
- Hippo/YAP (+2.445), ERBB Signaling (+2.021), IGFBPs (+2.400).
- 혈관신생 역설: 항혈관신생(+2.194) + 촉진(+2.157) 동시 — PDAC 저혈관성과 VEGF 드라이버 동시 반영.

**Novel 발견**
- **TCR/ZAP-70 하향 — PDAC vs 췌장염 cfRNA 감별자**: ZAP-70 Translocation To Immunological Synapse (−2.301), Downstream TCR Signaling (−2.377), Primary Immunodeficiency (−2.324). 췌장염에서는 이 신호 없음. Moore et al. 2025 원본이 간 기능 바이오마커(LAMC2, CEACAM7 등)로 감별함 → T세포 면역억제가 cfRNA 경로 수준의 PDAC 특이 감별자로 기능함은 미보고.
- NF-κB 비정규 경로 하향 (−2.450), Hedgehog 억제 (−2.330): PDAC 면역억제 TME의 다층적 cfRNA 신호.

**기대했으나 미검출 신호**
- **KRAS/MAPK/ERK 경로**: PDAC의 90%이상 KRAS G12D/V 변이. RhoGTPase(KRAS 하류)는 검출되었으나 KRAS → RAF → MEK → ERK 직접 경로가 상위에 명시되지 않음. cfRNA에서 활성 KRAS 돌연변이 발현 세포의 신호가 기질·면역 세포 cfRNA에 희석되는 현상.
- SMAD4/TGF-β 신호 소실: SMAD4 deletion이 PDAC ~55%에서 발생, 예후 불량 마커. SMAD/TGF-β 상향은 식도암에서 강하나 PDAC에서는 별도 농축되지 않음.
- 직접 TP53 gain-of-function 신호: PDAC TP53 변이 ~75%이나 GOF TP53 하류 경로 특이 농축 없음.

---

### 9. Esophagus Cancer (n=25 · FDR 683개 · ↑148 ↓535) — HIGH

**핵심 신호**
- 배아 기관 발달 최강 (+2.673): SOX17, SOX18, TEAD1/2/4, WNT7B, FGFR2, GATA3, GLI2, HEY1/2, TGFB1/2, FOXE1, SALL1 — 배아 전사 프로그램 재활성화.
- SMAD/TGF-β (+2.281): SMAD5/6/9, TGFB1/2, BMP2/5/10, NODAL → EMT.
- IGF1R 축 (ranks 5·7·10): IGF2, IRS1/2/4, AKT2 — Barrett 식도→선암 진행.
- GABA/니코틴/이온채널 (ranks 4·6·8·13·14): GABRA3/4/5, GRIN1/2A/2B/2C, CHRNA4/7, CACNA1A/1B — 신경내분비 분화 또는 흡연 관련 수용체.

**Novel 발견**
- **FGFR2-YAP1 공동 활성화 축의 cfRNA 동시 검출**: FGFR2 관련 5개 Reactome 경로 동시 상향 + YAP1/WWTR1 (+2.573). FGF3/5/6/7/8/9/10/16/17/22/23, FGFBP1/2/3, FGFR2, IRS1/2/4. FGFR2 → YAP1 공동 활성화 축은 위암 조직(PMC7581496)에서 보고되었으나 혈장 cfRNA에서 동시 농축을 GSEA 수준에서 보고한 연구 없음. 식도 선암의 FGFR2 증폭/과발현 배경과 일치.

**기대했으나 미검출 신호**
- **TP53 경로 (식도 편평세포암 ~100% TP53 변이, 선암 ~70%)**: TP53 GOF 또는 LOF 하류 경로가 상위에 명시되지 않음 — 가장 빈번한 드라이버이나 cfRNA 경로 수준 미검출.
- HER2/ERBB2 증폭 경로 (선암 15-32%에서 HER2 증폭·과발현): ERBB2 신호 경로가 식도암 상위에서 별도 강하게 나타나지 않음 (FGFR2 신호에 묻힘).
- PI3KCA/mTOR 경로 (식도 편평세포암 빈발 변이): mTOR 경로가 상위 경로로 부상하지 않음.
- TP63/SOX2 (편평상피 분화 마커, ESCC에서 증폭): 편평상피-특이 경로가 배아 발달 경로에 혼재되어 별도 식별 어려움.

---

### 10. CAD_HF+ vs CAD_HF− (HF+ n=116 / HF− n=108) — HIGH*

**HF+ 핵심 신호**
- **FGF23–Klotho–FGFR1 축 6회 반복 (HF 진행 서명)**: FGFR1 Ligand Binding (+2.82, 전체 1위) — FGF23, FGF5, FGF1, FGFR1, KL(Klotho). Phospholipase C-FGFR1 (+2.40), FGFRL1 Modulation Of FGFR1 (+2.26), FGFR1c Ligand Binding (+2.11), Activated Point Mutants Of FGFR1 (+2.02), FGFR2 변이체 (+2.09). 동일 생물학 축 6중 반복 = 신호 견고성 고도.
- 심장 섬유화: Collagen Chain Trimerization (+2.45) — COL1A1/COL1A2/COL5A1/COL14A1.
- 이온채널 리모델링: Calcium Ion Import (+2.57), Sarcoplasmic Reticulum Calcium Ion Transport (+1.95), Sodium/Potassium Ion Transport (+2.30/+2.35) — 부전 심근 전기적 리모델링.
- 전응고: F2R(PAR-1), THBS1 (Thrombospondin-1) — THBS1은 HFpEF mitophagy 억제 인자.
- Semaphorin-Plexin (+2.60), Nicotine addiction (+2.32).

**HF− 비교 (상향 단 11개, 심장 리모델링 축 부재)**
- 감각/통증 지각 (+2.52/+2.51), Beta Defensins (+2.50), Mucin 당화 결함, Purinergic Receptors. FGFR1·콜라겐·칼슘채널 심장 리모델링 축 완전 부재 → **HF+ 신호의 음성 대조군 확립**.

**Novel / 특이 발견**
- Ward et al. 2022이 동일 코호트 fold-change 분석에서 HF+/HF− 차이를 검출하지 못한 것(null result)과 대조적으로 정규화 모델이 **생물학적으로 정합한 HF 진행 서명을 검출**. FGF23은 incident HF 독립 예측 인자, Klotho는 ACS 사망 독립 예측 인자로 문헌에서 확립됨.

**⚠ 중요 주의사항**: 원본 논문의 null result와의 불일치로 배치 효과·세포 조성 차이 인공산물 가능성 배제 불가. **가설 생성적 수준 — 독립 코호트 및 직교 측정(qPCR, 혈장 단백질 FGF23·NT-proBNP)으로 검증 필수.**

**기대했으나 미검출 신호 (HF+)**
- BNP/NT-proBNP 전사 경로 (심부전 gold-standard 바이오마커): 분비 단백질의 cfRNA는 검출될 수 있으나 natriuretic peptide 신호 경로가 상위에 명시되지 않음.
- 심근세포 apoptosis 경로 (Caspase/Bcl-2): HF+ 군에서 직접적 심근세포 사멸 신호가 포착될 것으로 예상되었으나 경로 수준에서 미식별.
- GDF-15 pathway: 최근 emerging HF biomarker로 각광받으나 GSEA 상위 경로 미포함.
- 직접 아테로마 불안정화 경로 (Matrix metalloproteinase/TIMP): CAD에서 플라크 불안정화 신호가 cfRNA에서 검출될 것으로 예상되었으나 명시적 농축 없음.

---

### 11. Pre-eclampsia (n=59 · FDR 230개 · ↑2 ↓228) — MODERATE-HIGH

**Novel 발견**
- **Neurexins And Neuroligins (+2.42, 1위)**: 시냅스 접착 분자 경로 — Postsynaptic Membrane Organization (+2.11). Moufarrej et al. Nature 2022이 세포형 수준에서 신경근육 cfRNA 풍부화를 보고했으나 **Neurexin/Neuroligin 특정 시냅스 경로를 GSEA 수준에서 포착한 선행 연구 없음**. 임신중독증 신경학적 합병증(두통, 시각 장애, 경련)과 직접 연결.
- DNA Methylation (+2.17): 영양막 기능 이상의 에피제네틱 기초.
- 선천 면역 하향 (~40개): ISG15 (−1.99), IFN Signaling (−1.88), NLRP3 (−1.80) — 면역 관용 붕괴.
- SLIT/ROBO 하향 (−1.94): 영양막 침습 실패 신호 가능성.
- 지질단백질 하향: Lipoprotein Localization (−1.91) — 임신중독증 특이 이상지질혈증.

**기대했으나 미검출 신호**
- **sFlt-1/PlGF 항혈관신생 축**: 임신중독증의 hallmark인 sFlt-1 과발현 → VEGF/PlGF 포획 → 항혈관신생 상태. VEGF 신호 경로가 상위 상향 경로에서 명시적으로 식별되지 않음 — 단 2개의 상향 경로가 한계.
- HIF-1α/저산소 반응 경로 (태반 ischemia의 촉발 기전): 태반 저산소증이 시발점임에도 HIF-1α 경로 농축 미확인.
- 보체 시스템 말단 활성화 경로: 임신중독증에서 C3a/C5a 과활성화 보고가 많으나 cfRNA 상향 신호 미포착.
- Renin-Angiotensin 경로 (혈압 조절 이상): 혈압 관련 신호 직접 미검출.

**⚠ 해석 주의**: 정상 임산부를 HC 기준으로 사용하지 않아 임신 상태 자체와 질병 신호 혼재. FDR<0.05 상향 경로 단 2개로 신호 풍부도 매우 제한적.

---

## Tier 2 — 질병 관련 신호 확인 케이스 (7개)

> 질병 관련 signature가 잘 뽑혔으나 주요 driver 경로 미식별 또는 특이한 novel finding 없음

---

### 12. Liver Cancer — Roskams-Hieter (n=28 · FDR 217개 · ↑211 ↓6) — HIGH

**재현된 신호**: 보체/응고(+2.224, +2.110) — Chen 코호트와 cross-cohort 재현성 확인된 HCC 코어 cfRNA 축. 세포주기 (G1/S +1.998, kinesin +1.992), Oncogene Induced Senescence (+1.910: CDKN2A/B, TP53, RB1, MDM2, CDK6), YAP/Hippo (+1.90), CDC42/RAC1 GTPase, 항바이러스(IFITM2/3, OAS2: HBV/HCV 배경 암시). OXPHOS 하향 (6개 경로: UQCR11, NDUFA1, COX17 — 두 코호트에서 재현된 Warburg 효과).

**특이 발견**: 두 HCC 코호트(Chen n=10 vs Roskams-Hieter n=28) 간 이질성이 cfRNA에서 명확히 구별됨 — 번역 억제 여부, EMT/침습 신호, YAP 경로 강도, 항바이러스 신호 모두 차이. 동일 HCC 진단 내 환자 이질성 반영.

**기대했으나 미검출 신호**: Wnt/CTNNB1 경로 (HCC ~30% CTNNB1 변이). HIF-1α/VEGF 혈관신생 (고혈관성 HCC의 특징이나 VEGF 경로 상위 미포함). AFP(alpha-fetoprotein) 생산 관련 전사 프로그램. DKK1 (HCC cfRNA 바이오마커로 보고된 Wnt 억제인자).

---

### 13. Liver Cancer — Chen (n=10 · FDR 625개 · ↑171 ↓454) — MODERATE-HIGH

**재현된 신호**: 간 분비 단백질 완전 세트 (ALB, FGA/B/G, PLG, C3, CFB, APOB, APOE, APOA1-4, SERPINC1, KNG1). Complement/Coagulation (+2.476), Fibrinolysis (+2.281). GCN2/Starvation Response (−2.715): 종양 미세환경 영양 결핍 ISR.

**기대했으나 미검출 신호**: Roskams-Hieter에서 나타난 세포주기/YAP 신호가 Chen에서는 거의 없음 — n=10 한계 또는 코호트 이질성. HIF-1α/VEGF. TP53/RB1 경로.

**n=10 한계로 신호 과잉 해석 주의.**

---

### 14. ICI-treated Cancer (n=11 · FDR 1,304개 · ↑231 ↓1,073) — HIGH

**재현된 신호**: Lysosome (+2.825), IgA Intestinal Immune Network (+2.141: ICI 관련 대장염), Autoimmune Thyroid Disease (+2.12: 갑상선염 ICI 부작용), Allograft/GvHD Rejection Pattern (+2.263/+2.177), Defensins (+2.376). SLIT/ROBO (−2.472/−2.503: 종양 미세환경 변화).

**특이 발견**: ICI 치료 관련 다발성 부작용(장관·갑상선·이식거부 패턴)이 cfRNA에서 동시 시스템적으로 포착. ICI-m과 감별자 명확(ER-phagosome/ERAD/소교세포 활성화 = ICI-m 고유).

**기대했으나 미검출 신호**: PD-1/PD-L1(CD274/PDCD1) 체크포인트 경로 직접 신호. 세포독성 T세포 팽창 경로 (GZMB, PRF1, NKG7 연계 경로 상위 미포함). 종양 cell death/apoptosis 신호.

---

### 15. Lung Cancer — Chen (n=30 · FDR 520개 · ↑54 ↓466) — MODERATE-HIGH

**재현된 신호**: Integrin (+2.726, 1위: LIMS1/2, NID1, LAMB2, FLNA, LOXL3), Caspase Cleavage of Cell Adhesion Proteins (+2.701: CTNNB1, CASP3, TJP1/2, DSG2 — anoikis CTC 신호), ERBB2 (+2.316/+2.029), EGFR (+2.090), RET (+2.189), Hippo/YAP (+2.372), Angiogenesis (+2.457/+2.354), Glycogen Biosynthesis (+2.596).

**기대했으나 미검출 신호**
- **KRAS G12C 경로 (NSCLC ~13%, 최근 sotorasib 타겟)**: KRAS 직접 신호 경로 미포착. RhoGTPase(KRAS 하류)는 있으나 KRAS → RAF → MEK → ERK 경로가 상위에 없음.
- ALK/ROS1 재배열 경로 (ALK+ NSCLC ~3-5%): 샘플 수에 ALK+ 환자가 충분하지 않거나 cfRNA 신호 희박.
- MET 증폭 경로 (MET exon 14 skipping ~3% NSCLC).
- SMAD/TGF-β EMT 경로 (식도암에서는 강했으나 폐암에서 약함): 폐암 EMT 기전이 integrin/caspase 경로로 반영되어 있으나 TGF-β 경로 자체는 상위 미포함.
- PD-L1(CD274) 발현 경로.

---

### 16. Stomach Cancer (n=24 · FDR 516개 · ↑13 ↓503) — HIGH

**재현된 신호**: Cell-Cell Junction Organization (+2.643): PKP4/3/1, DSG2, CDH3(P-cadherin), TJP1/2, CLDN1/5, CTNND1, FSCN1 — 위암 상피 접합 리모델링. ECM/인테그린 (+2.477/+2.328: FLNA, DMTN, NID1, LIMS2, LAMB2).

**주의**: Formation Of Cornified Envelope (NES +2.583, 편평상피 화생 해석)이 현재 데이터에서 FDR<0.05 탈락 — 해당 해석 제외.

**기대했으나 미검출 신호**
- **HER2/ERBB2** (위암 15-20%, 트라스투주맙 1차 타겟): ERBB2 경로가 다른 암(CRC, 폐암)에서는 검출되었으나 위암 상위 경로에서 별도 강하게 나타나지 않음.
- **CDH1(E-cadherin) loss 경로** (미만형 위암·반지세포암의 핵심 변이): E-cadherin 손실 → N-cadherin switching 경로가 상위에 없음. 접합 단백질 상향은 CDH1 손실과 상반되는 방향으로도 해석 가능.
- EBV-related signal (EBV+ 위암 ~10%): 항바이러스 경로 미포함.
- MET amplification, PIK3CA, ARID1A 경로.
- RET 경로 (+2.298)는 검출되었으나 위암에서 RET 재배열 빈도 ~0.5-2%로 소수 — 주요 driver 신호로 해석하기 어려움.

---

### 17. MGUS (n=8 · FDR 107개 · ↑92 ↓15) — MODERATE

**재현된 신호**: 에피제네틱 선행 변화 — DNA Methylation (+2.177: DNMT3A/3B), SIRT1 rRNA Repression (+2.127), HDACs (+2.118), PRC2 (~+1.96), Telomere Packaging/Senescence (+2.073/+2.099). YAP/NOTCH 전악성 신호. NK세포 살상 억제 경로 하향 (−2.19/−2.13: HLA 손실 초기 신호). 리소좀 하향 (−2.368).

**기대했으나 미검출 신호**: BLIMP1/PRDM1, IRF4 형질세포 정체성 경로. MGUS vs Smoldering MM 구별 경로. MAPK/MEK 경로 (MGUS 진행 관련 변이). **n=8으로 모든 신호 재현성 제한.**

---

## Tier 3 — 애매한 케이스 (3개)

> 질병 관련 신호가 나오긴 했으나 간접 정보 위주이거나 PDAC와 중복·배경 신호에 묻힘

---

### 18. Pancreatitis (n=79 · FDR 263개 · ↑27 ↓236) — MODERATE

상향 경로 27개 중 Hippo/YAP (+2.423), RHOC (+2.198), ECM (+2.082) 등이 **PDAC와 공유** — 단독으로는 췌장염-특이 신호 구별 어려움. CFTR 경로 하향 (−2.219)이 유일한 췌장염-특이 Novel 신호(췌관 기능 이상 cfRNA 미보고)이나 하향 방향으로 대부분 OXPHOS(−2.658) 등 confound에 묻힘.

**기대했으나 미검출 신호**: 급성 췌장염 트립시노겐 활성화 경로 (PRSS1). NF-κB 염증 경로 상향 (PDAC에서는 하향으로 나타남). 직접적인 선포 세포 손상 마커.

---

### 19. CAD_HF− (n=108 · FDR 922개 · ↑11 ↓911) — HIGH

상향 11개가 감각/통증(+2.52), Beta Defensin(+2.50), Mucin O-당화 결함(+2.49), Purinergic Receptors(+2.30) 등 심혈관과 직접 관련 없는 신호 위주. 하향 배경(번역·OXPHOS)은 HF+와 동일. **FGFR1·콜라겐·칼슘채널 심장 리모델링 축 완전 부재** = HF 비진행의 음성 대조군으로 기능.

**기대했으나 미검출 신호**: CAD 자체의 cfRNA 신호(아테로마·내피 기능 이상·지질단백질 대사). 이 코호트가 "안정형 CAD" 기저 환자이므로 CAD 활성 신호가 약한 것으로 해석되나, 정상 HC 대비 어떤 질병-특이 상향 신호도 식별되지 않음.

---

### 20. Other Cancer (n=16 · FDR 266개 · ↑17 ↓249) — HIGH (신호 일관성)

RhoGTPase 완전 세트(RHOB +2.41, RHOC +2.18, RHOA +2.10, RHOQ +2.10, CDC42 +2.10, RAC1 +2.03) 동시 상향 — 이질적 암종 집합에서도 일관된 침습 서명. SLIT/ROBO 하향 (−2.89/−2.68), Wnt/β-카테닌 간접 활성화 (프로테아솜 하향 → DVL/AXIN 분해 억제). 그러나 혼합 암종 구성으로 **어떤 암종의 어떤 특이 신호인지 귀속 불가** → 임상 활용 제한.

**기대했으나 미검출 신호**: MYC/MYCN 전사 프로그램, PI3K/AKT/mTOR (범암 드라이버), TP53 LOF 공통 신호. 암종-특이 신호 없음이 이 그룹의 본질적 한계.

---

## 종합 Novel 발견 9건 위치 정리

| # | Novel 발견 | 해당 질환 | Tier |
|---|---|---|---|
| 1 | ER-phagosome/교차제시 cfRNA 농축 | ICI-m | 1 |
| 2 | TB 스플라이세오솜 하향 cfRNA 신호 | Tuberculosis | 1 |
| 3 | HIV+TB B세포 경로 하향 클러스터 | HIV+TB | 1 |
| 4 | CRC Hippo/YAP cfRNA liquid biopsy 검출 | CRC | 1 |
| 5 | ESCC/EAC FGFR2-YAP1 공동 활성화 축 cfRNA | Esophagus | 1 |
| 6 | PDAC vs 췌장염 TCR/ZAP-70 cfRNA 감별자 | Pancreatic Cancer | 1 |
| 7 | 췌장염 CFTR 경로 cfRNA 억제 | Pancreatitis | 3 |
| 8 | ME/CFS FoxO/GCN2/인슐린 축 cfRNA 억제 | ME/CFS | 1 |
| 9 | 임신중독증 Neurexin/Neuroligin cfRNA 경로 | Pre-eclampsia | 1 |

---

*GSEA Master Report Summary | 작성: 2026-06-29*  
*원본: GSEA_Master_Report.md (2026-06-24, 2026-06-25 개정)*
