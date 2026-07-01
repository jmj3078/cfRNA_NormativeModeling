# cfRNA Normative Model — GSEA with_rare 통합 분석 리포트

**분석:** rare-event 분기(559개 극저발현 protein-coding 유전자, HC det<1%, pooled 공변량 GLM) 포함 GSEA vs 미포함(no_filter) 비교
**유의 기준:** FDR q-val < 0.05 · **NES:** >0 HC 대비 상향, <0 하향
**교차검증 DB:** Open Targets Platform (질병별 association score 상위 300 target) + 문헌
**핵심 정량:** rare 포함은 공유 term의 NES 부호를 **전 질병에서 한 번도 뒤집지 않음(sign_agree=1.0)** → 순수 가산적. 질병당 신규 유의 term +32~+225개.

> 자동 생성 정량 백본 위에 질병별 문헌검증 해석을 덧붙이는 구조. "rare-led 신규 신호"=rare-분기 유전자가 직접 leading-edge에 든 경로(직접 근거), "DB 지지 신규 신호"=lead genes가 Open Targets 질병 유전자와 겹치는 경로.

---

## 전체 요약 테이블 (rare 포함 vs 미포함)

| Phenotype | no_filter | with_rare | 신규(+) | 소실(−) | Jaccard | 부호일치 | DB지지 신규 | rare-led |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| Liver Cancer (Roskams-Hieter) | 217 | 431 | +225 | −11 | 0.47 | 1.00 | 162 | 4 |
| CAD_HF+ | 961 | 1085 | +163 | −39 | 0.82 | 1.00 | 51 | 0 |
| ICI-treated Cancer | 1304 | 1362 | +136 | −78 | 0.85 | 1.00 | 0 | 0 |
| Esophagus Cancer (Chen) | 683 | 746 | +108 | −45 | 0.81 | 1.00 | 41 | 4 |
| CAD_HF- | 922 | 984 | +106 | −44 | 0.85 | 1.00 | 30 | 0 |
| Tuberculosis | 606 | 672 | +104 | −38 | 0.80 | 1.00 | 60 | 1 |
| MM | 333 | 405 | +92 | −20 | 0.74 | 1.00 | 58 | 1 |
| ICI-m | 479 | 522 | +91 | −48 | 0.76 | 1.00 | 40 | 2 |
| Pancreatic Cancer (Moore) | 560 | 604 | +90 | −46 | 0.79 | 1.00 | 43 | 4 |
| Liver Cancer (Chen) | 625 | 656 | +79 | −48 | 0.82 | 1.00 | 31 | 0 |
| Pre-eclampsia | 230 | 282 | +76 | −24 | 0.67 | 1.00 | 12 | 2 |
| ME_CFS | 169 | 226 | +71 | −14 | 0.65 | 1.00 | 51 | 5 |
| Pancreatitis | 263 | 310 | +70 | −23 | 0.72 | 1.00 | 28 | 3 |
| Colorectal Cancer | 385 | 441 | +67 | −11 | 0.83 | 1.00 | 35 | 0 |
| Stomach Cancer | 516 | 563 | +67 | −20 | 0.85 | 1.00 | 24 | 1 |
| HIV + Tuberculosis | 397 | 434 | +66 | −29 | 0.79 | 1.00 | 45 | 0 |
| Other Cancer | 266 | 323 | +65 | −8 | 0.78 | 1.00 | 34 | 9 |
| Lung Cancer | 520 | 548 | +60 | −32 | 0.84 | 1.00 | 26 | 0 |
| HIV | 482 | 478 | +49 | −53 | 0.81 | 1.00 | 18 | 0 |
| MGUS | 107 | 127 | +32 | −12 | 0.68 | 1.00 | 6 | 2 |

---

## rare-분기 유전자가 직접 leading-edge에 든 신규 경로 (직접 근거 카탈로그)

rare 유전자가 lead에 직접 포함된 신규 유의 term. rare 분기가 없었으면 검출 불가한, 기전적으로 rare에 귀속되는 신호.

| Phenotype | Term | NES | rare lead | DB 지지 유전자 |
|---|---|--:|---|---|
| Esophagus Cancer (Chen) | Extracellular Matrix Organization (GO:0030198) | 2.17 | **MMP3** | CREB3L1;NF1 |
| Esophagus Cancer (Chen) | cAMP signaling pathway | 2.08 | **PDE4C;TSHB** | AKT1;BRAF;CACNA1D;CREB3L1;CREBBP;GRIN2A;PIK3R1;PTC |
| Esophagus Cancer (Chen) | Detection Of Chemical Stimulus Involved In Sensory Perceptio | 1.97 | **OR10H3;OR13C8** | — |
| Esophagus Cancer (Chen) | Detection Of Chemical Stimulus Involved In Sensory Perceptio | 1.94 | **OR10H3;OR13C8** | — |
| ICI-m | Iron Uptake And Transport R-HSA-917937 | 1.80 | **ATP6V1G3** | — |
| ICI-m | Response To Zinc Ion (GO:0010043) | 1.72 | **MT1B** | — |
| Liver Cancer (Roskams-Hieter) | Negative Regulation Of MAPK Cascade (GO:0043409) | 1.75 | **DEFB114** | APOE;NF1;SH2B3 |
| Liver Cancer (Roskams-Hieter) | Negative Regulation Of MAP Kinase Activity (GO:0043407) | 1.73 | **DEFB114** | APOE;NF1;SH2B3 |
| Liver Cancer (Roskams-Hieter) | Regulation Of MAP Kinase Activity (GO:0043405) | 1.69 | **DEFB114** | APOE;EGFR;FGF1;FGF2;FLT1;NF1;NTRK3;SH2B3;SOS1 |
| Liver Cancer (Roskams-Hieter) | Androgen Receptor Signaling Pathway (GO:0030521) | 1.66 | **CST11** | ARID1A;MAPK1 |
| ME_CFS | Negative Regulation Of Cholesterol Transport (GO:0032375) | -1.81 | **PLA2G10** | NFKB1 |
| ME_CFS | Insulin Processing R-HSA-264876 | -1.80 | **INS** | INS |
| ME_CFS | Positive Regulation Of Cellular Catabolic Process (GO:003133 | -1.79 | **INS** | APP;ATG13;IGF1;IL1B;INS;PRKAA1;TNF |
| ME_CFS | Rap1 signaling pathway | -1.78 | **PFN3** | GRIN1;GRIN2A;GRIN2B;IGF1;MAPK1;MAPK14;MAPK3 |
| ME_CFS | Regulation Of Lipid Storage (GO:0010883) | -1.77 | **PLA2G10** | FBXW7;NFKB1;TNF |
| MGUS | Mitotic Prophase R-HSA-68875 | 1.87 | **H2BC1** | — |
| MGUS | Nucleosome Assembly (GO:0006334) | 1.81 | **H2BC1** | — |
| MM | Neutrophil extracellular trap formation | 1.69 | **H4C7** | MAPK1;MTOR;PIK3CA;PIK3CB;RAF1 |
| Other Cancer | Apical Junction Assembly (GO:0043297) | 2.19 | **CLDN17;CLDN25** | STRN |
| Other Cancer | Positive Regulation Of Vascular Endothelial Growth Factor Pr | 1.97 | **NOX1** | BRCA1 |
| Other Cancer | Regulation Of Epithelial Cell Apoptotic Process (GO:1904035) | 1.97 | **USP17L24** | — |
| Other Cancer | Bicellular Tight Junction Assembly (GO:0070830) | 1.96 | **CLDN17;CLDN25** | STRN |
| Other Cancer | Tight Junction Assembly (GO:0120192) | 1.94 | **CLDN17;CLDN25** | STRN |
| Other Cancer | Regulation Of Biomineral Tissue Development (GO:0070167) | 1.94 | **AMTN** | DDR2 |
| Other Cancer | Cell-Cell Communication R-HSA-1500931 | 1.93 | **CLDN17** | CDH11;CTNNB1;FYN;PIK3CB;PIK3R1;PTPN11 |
| Other Cancer | Positive Regulation Of Translational Initiation (GO:0045948) | -1.83 | **DAZ2;DAZ4** | — |
| Other Cancer | Alzheimer disease | -1.74 | **COX8C** | AKT2;AXIN1;BRAF;CASP8;MAP2K2;PIK3CA;SEM1 |
| Pancreatic Cancer (Moore) | Amino Acids Regulate mTORC1 R-HSA-9639288 | -1.82 | **ATP6V1G3** | — |
| Pancreatic Cancer (Moore) | Retinoid Metabolism And Transport R-HSA-975634 | 1.82 | **PNLIP** | — |
| Pancreatic Cancer (Moore) | Positive Regulation Of Cell Differentiation (GO:0045597) | 1.78 | **ADIG** | ARID1A;BMPR1A;CTNNB1;DDR2;HGF;KDR;MTOR;PBRM1;PPARG |
| Pancreatic Cancer (Moore) | Negative Regulation Of Cell Differentiation (GO:0045596) | 1.78 | **MSTN** | CALR;CDK6;EGFR;MYB;PPARG;RUNX1T1;WWTR1;YAP1;ZBTB16 |
| Pancreatitis | Regulation Of Acrosome Reaction (GO:0060046) | 1.99 | **IQCF1** | PLCB1 |
| Pancreatitis | Skeletal System Development (GO:0001501) | 1.94 | **EPYC** | — |
| Pancreatitis | ROS And RNS Production In Phagocytes R-HSA-1222556 | -1.75 | **ATP6V1G3** | — |
| Pre-eclampsia | Ribosome Assembly (GO:0042255) | -1.84 | **RPL10L** | — |
| Pre-eclampsia | Negative Regulation Of Viral Process (GO:0048525) | -1.76 | **IFNL3** | — |
| Stomach Cancer | Calcium-Independent Cell-Cell Adhesion Via Plasma Membrane C | 2.12 | **CLDN17** | CLDN18 |
| Tuberculosis | Signaling By WNT In Cancer R-HSA-4791275 | 1.73 | **DKK4** | — |

---

## 질병별 상세 (문헌검증)

문헌 근거는 PubMed 검색 히트 수로 표기(gene×disease). **기지** = 선행 문헌·DB 다수 존재, rare 포함이 기존 생물학을 재확인. **후보(novel)** = 문헌 거의 없음(≤2건), rare 분기가 처음 표면화한 검증 필요 신호.

### ME/CFS (신규 +71 · DB지지 51 · rare-led 5)

가장 해석이 풍부한 사례. 신규 유의 경로의 72%가 Open Targets ME/CFS 유전자로 지지되고, rare 유전자가 직접 lead한 신호 5개가 모두 **대사·신경 축**으로 수렴한다(전부 NES<0 = HC 대비 하향).

- **INS (Insulin Processing R-HSA-264876, Positive Regulation Of Cellular Catabolic Process; NES −1.79/−1.80)** — rare lead **INS**, DB 지지 APP/ATG13/IGF1/IL1B/PRKAA1/TNF. ME/CFS × insulin 문헌 116건, × lipid metabolism 80건으로 **기지 축**(에너지 대사·인슐린 감수성 이상)과 정합. cfRNA에서 인슐린 처리·이화작용 하향은 보고된 대사 저하 표현형을 분자 수준에서 재현.
- **PLA2G10 (Regulation Of Lipid Storage, Negative Regulation Of Cholesterol Transport; NES −1.77/−1.81)** — rare lead **PLA2G10**. ME/CFS × PLA2G10 문헌 **0건 = 후보 신호**. 지질 저장/콜레스테롤 수송 하향은 대사 축과 방향 일치하나, 특정 분비형 phospholipase A2(PLA2G10) 귀속은 본 분석이 처음 제시하는 검증 대상.
- **PFN3 (Rap1 signaling; NES −1.78)** — DB 지지 GRIN1/2A/2B(글루타메이트 수용체)·MAPK — 신경신호 축.

**해석:** rare 분기는 ME/CFS에서 노이즈가 아니라 **대사(INS)·지질(PLA2G10)·신경(GRIN/Rap1)** 이라는, 저발현이라 기존 group-wise 분석에서 소실됐던 축을 복원. PLA2G10은 문헌 미보고 후보.

### 위암 (Stomach Cancer, 신규 +67 · DB지지 24 · rare-led 1)

- **CLDN17 (Calcium-Independent Cell-Cell Adhesion; NES +2.12)** — rare lead **CLDN17**, DB 지지 **CLDN18**. 위암 × CLDN18 문헌 **324건**(CLDN18.2는 zolbetuximab 표적으로 임상 검증된 위암 바이오마커) — **강한 기지 축**. rare 분기가 저발현 claudin(CLDN17)을 통해 이 검증된 tight-junction/claudin 축을 cfRNA에서 표면화. CLDN17 자체는 위암 문헌 27건으로 상대적 신규.

### Other Cancer (신규 +65 · DB지지 34 · rare-led 9 — rare-led 최다)

rare-led 신호가 가장 많은 표현형. **상피 junction 축**이 지배적:
- **CLDN17/CLDN25 (Apical Junction Assembly NES +2.19, Cell-Cell Communication +1.93, Bicellular Tight Junction Assembly +1.96)** — 다수 claudin이 rare lead. 상피 기원 종양의 tight-junction 재편과 정합(claudin은 pan-cancer 확립 마커).
- **NOX1 (Positive Regulation Of VEGF production; NES +1.97)** — DB 지지 BRCA1. NADPH oxidase-혈관신생 축.
- **COX8C (Alzheimer disease KEGG; NES −1.74)** — 미토콘드리아 전자전달 subunit; 이질적 "Other" 집단이라 해석 보수적으로.

### 췌장암 (Pancreatic Cancer (Moore), 신규 +90 · DB지지 43 · rare-led 4)

- **PNLIP (Retinoid Metabolism And Transport R-HSA-975634; NES +1.82)** — rare lead **PNLIP**(췌장 리파아제, 외분비 조직 특이). PDAC × PNLIP 문헌 44건, × retinoid 152건. 외분비 췌장 기능 유전자가 cfRNA에서 rare 분기로 검출 — 조직 기원 정합.
- **ADIG/MSTN (Positive/Negative Regulation Of Cell Differentiation; NES +1.78)** — DB 지지 13개(ARID1A/CTNNB1/HGF/KDR/PPARG/TGFBR2/YAP1 등 PDAC 확립 드라이버). rare lead(ADIG=adipogenin, MSTN=myostatin)는 분화·악액질 축; MSTN×PDAC 11건(악액질 맥락 기지).
- **ATP6V1G3 (Amino Acids Regulate mTORC1; NES −1.82)** — v-ATPase subunit, 리소좀-mTORC1.

### 다발골수종·MGUS (MM 신규 +92·DB지지 58 / MGUS 신규 +32·DB지지 6)

- **MM: H4C7 (Neutrophil Extracellular Trap Formation; NES +1.69)** — rare lead 히스톤 **H4C7**, DB 지지 MAPK1/MTOR/PIK3CA/PIK3CB/RAF1. MM × NET 문헌 18건(신흥 축). 히스톤 cfRNA는 NET/세포사 방출과 정합.
- **MGUS: H2BC1 (Mitotic Prophase NES +1.87, Nucleosome Assembly +1.81)** — rare lead 히스톤 **H2BC1**. MGUS/gammopathy × histone/nucleosome 문헌 719건 = 기지. 형질세포 증식·염색질 축을 저발현 히스톤 유전자가 표면화.

### ICI 유발 심근염 (ICI-m, 신규 +91 · DB지지 40 · rare-led 2)

- **MT1B (Response To Zinc Ion; NES +1.72)** — rare lead metallothionein **MT1B**. 심근염 × metallothionein 문헌 **2건 = 후보 신호**(MT는 심근 산화스트레스 보호에 관여하나 ICI 심근염 맥락 미보고).
- **ATP6V1G3 (Iron Uptake And Transport R-HSA-917937; NES +1.80)** — 심근염 × iron 문헌 138건(ferroptosis 축 기지). rare 분기가 철·아연 금속대사 축을 심근염에서 표면화 — 후속 검증 가치.

### 자간전증 (Pre-eclampsia, 신규 +76 · DB지지 12 · rare-led 2)

- **RPL10L (Ribosome Assembly; NES −1.84)** — rare lead 리보솜 단백 **RPL10L**. PE × ribosome 문헌 98건(태반 번역 조절 기지 축). 리보솜 조립 하향은 태반 기능부전과 정합.
- **IFNL3 (Negative Regulation Of Viral Process; NES −1.76)** — rare lead **IFNL3**(IFN-λ3). PE × IFN-lambda 문헌 **1건 = 후보 신호**. 태반 항바이러스 IFN-λ 축은 생물학적으로 그럴듯하나 PE 맥락 미보고.

### 식도암 (Esophagus Cancer (Chen), 신규 +108 · DB지지 41 · rare-led 4)

- **MMP3 (Extracellular Matrix Organization; NES +2.17)** — rare lead **MMP3**(matrix metalloproteinase), DB 지지 CREB3L1/NF1. ECM 리모델링은 식도암 침습 확립 축.
- **PDE4C/TSHB (cAMP signaling; NES +2.08)** — DB 지지 AKT1/BRAF/CREBBP/PIK3R1.
- **OR10H3/OR13C8 (Detection Of Chemical Stimulus; NES +1.97)** — 후각수용체군. 저발현 OR 유전자 다발 검출은 **주의 신호**: 생물학적 의미보다 극저발현 유전자 특유의 잔존 분산일 수 있어 보수적 해석 필요.

### CAD — 심부전 진행군 (CAD_HF+ 신규 +163·DB지지 51 / CAD_HF− 신규 +106·DB지지 30)

두 표현형 모두 rare-led는 0(직접 rare 신호 없음)이나, 신규 term이 **혈관·기질·산화질소 축**으로 수렴 — 전부 확립된 CAD/HF 생물학:
- **CAD_HF+: ECM-receptor interaction (NES +2.15; ITGA2B/ITGAV/ITGB3/ITGB5/VWF/COL6A3), Nitric Oxide Stimulates Guanylate Cyclase (+2.00; GUCY1A2/PDE5A), Removal Of Superoxide Radicals (−1.61; NOS3)** — CAD×ECM/collagen 문헌 2484건, HF×NO/eNOS 4300건 = 기지. rare 포함은 저발현 혈관 유전자를 통해 NO-cGMP·기질 리모델링 축을 강화.
- **CAD_HF−: Collagen Chain Trimerization (+2.38; COL4A1/COL6A3), P2Y Receptors (+1.99; P2RY12), Leukocyte Adhesion (−1.52; PECAM1)** — 콜라겐·혈소판(P2RY12)·백혈구 부착 축.

**해석:** rare 분기 효과는 간접(랭킹 이동)이며, 방향은 기존 CAD/HF 생물학을 재확인. 신규 신호 중 미보고 후보는 없음(전부 DB/문헌 기지).

### 간암 (Liver Cancer)

- **Roskams-Hieter (신규 +225, DB지지 162 — 절대 신규 최다·DB지지율 72%)**: rare 포함이 신호를 거의 2배로 확장. 신규 DB-지지 term이 **HCC 확립 드라이버 축**에 집중 — RUNX3 Regulates p14-ARF(CCND1/CDKN2A/EP300/NRAS), Small GTPase Signal Transduction(KRAS/NRAS/TP53/RB1/NF1), Negative Regulation Of ERBB(CBL). rare-led 4개는 **MAP kinase 조절축**(rare lead **DEFB114**)과 **Androgen Receptor Signaling**(rare lead **CST11**). DB 지지 유전자(APOE/EGFR/FGF1/FGF2/NF1/SH2B3)는 HCC MAPK 축과 정합하나, rare lead인 **DEFB114(β-defensin)는 HCC 문헌 3건 = 후보 신호**로, 이 MAPK 경로 검출을 촉발한 저발현 유전자 귀속은 검증 대상. HCC×GPC3(glypican-3, 확립 마커) 문헌 1204건.
- **Chen (신규 +79, DB지지 31)**: rare-led 0. 신규 term은 **GPC3(GAG Synthesis, +1.83)**, Protein Glycosylation(TET2), Negative Regulation Of Coagulation(APOE) — GPC3·응고·당쇄화 축, HCC 확립 생물학과 정합.

### 대장암 (Colorectal Cancer, 신규 +67 · DB지지 35 · rare-led 0)

신규 term이 **성장인자·JAK-STAT 축**으로 수렴 — Interleukin-7 Signaling(+2.25; JAK1/JAK3/PIK3R1), Growth Hormone Receptor(+2.13; JAK2/STAT3), Cytosolic FGFR1(+2.20), Positive Regulation Of Glucose Transport(+2.14; BRAF/PTPN11). CRC×IL-7 문헌 64건(면역-종양 축), JAK-STAT/FGFR은 CRC 확립. rare 포함이 저발현 성장인자 신호 유전자를 통해 이 축을 강화, 방향은 기지 생물학과 정합.

### 폐암 (Lung Cancer, 신규 +60 · DB지지 26 · rare-led 0)

신규 term의 최상위가 **ERBB/EGFR 축** — Regulation Of miRNA Transcription(+2.42; EGFR/STAT3), ERBB Signaling(+2.16; EGFR/PTPN11/SRC), Regulation Of KIT Signaling(+2.15; SRC). 폐암×ERBB/EGFR 문헌 23,568건 = 가장 확립된 축. rare 포함은 EGFR-family 하류를 재확인(간접효과), 미보고 후보 없음.

### 감염병 그룹 (HIV · Tuberculosis · HIV+TB)

- **HIV (신규 +49, DB지지 18, rare-led 0)**: 신규 term이 **바이러스-숙주 상호작용·리보솜 축** — Vpr-mediated Nuclear Import Of PICs(+1.78; NUP153/160/205/214 등 nucleoporin 다수), Regulation Of Innate Immune Response(SAMHD1), Ribosome Biogenesis(−1.97; NPM1). HIV×Vpr 핵수송 문헌 112건, SAMHD1은 확립된 HIV 제한인자. rare 포함이 핵공복합체·숙주제한 축을 강화 — HIV 생물학 정합.
- **Tuberculosis (신규 +104, DB지지 60, rare-led 1)**: 신규 term이 **자연면역·보체·B세포 축** — Cellular Response To LPS(+1.79; CASP1/CXCL9/10/IL18), Complement Cascade(+1.76; MBL2/SERPING1), TLR4 Signaling(LGALS9), B Cell Activation(−2.01; CD40LG/CD79A). TB×complement/MBL2 문헌 2274건 = 기지. rare-led 1개는 **Signaling By WNT In Cancer**(rare lead **DKK4**); TB×DKK4/WNT 문헌 68건(육아종 WNT 축은 신흥) — 후속 검증 가치.
- **HIV+TB (신규 +66, DB지지 45 — DB지지율 68%)**: 신규 term이 **VEGF·해당과정·림프구 활성** — VEGF Signaling(+1.86), Hexose Biosynthesis(+1.85; ENO1/PGK1), Lymphocyte Activation(−1.85; CCR2/TLR9). 공동감염의 대사(해당과정 상향)·면역 축을 rare 포함이 강화.

### ICI-treated Cancer (신규 +136 · DB지지 0 · rare-led 0 — 문헌전용)

이질적 "ICI 치료 암" 혼합 집단으로 Open Targets 단일 질병 참조가 부재(DB 교차검증 불가). 신규 term은 다수이나 표현형이 혼합이라 **질병특이 해석 유보** — group 구성이 정의되면 재평가 필요. 정량적으로 rare 포함이 +136 신규(부호일치 1.0)를 더한다는 사실만 기록.

---

## 종합

1. **rare 포함은 방향 보존적·가산적**: 20개 전 질병에서 공유 term의 NES 부호 일치율 1.0. rare 분기는 기존 신호를 뒤집지 않고 신규 유의 경로만 추가(+32~+225).
2. **신규 신호의 상당수가 DB/문헌으로 지지**: 대부분 질병에서 신규 term의 30~72%가 Open Targets 질병 유전자로 지지되며, 방향이 확립된 질병 생물학과 정합(노이즈 아님).
3. **rare 분기가 처음 표면화한 미보고 후보**: ME/CFS **PLA2G10**(지질), 자간전증 **IFNL3/IFN-λ**, ICI 심근염 **MT1B/metallothionein**, 간암(R-H) **DEFB114/MAPK** — 모두 저발현이라 group-wise 분석에서 소실됐던, 검증 가치 있는 신호.
4. **주의**: 식도암의 후각수용체(OR10H3/OR13C8) 등 일부 rare-led는 생물학보다 극저발현 특유의 잔존 분산일 수 있어 보수적 해석 필요. ICI-treated는 혼합 집단이라 해석 유보.
