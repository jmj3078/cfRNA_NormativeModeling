# cfRNA ORA 결과 — 질병별 생물학적 적절성 평가 (전체 20개 완료)

각 질병의 Z-score 기반 선정 유전자에 대한 Enrichr ORA(GO-BP / KEGG / Reactome, adj-P<0.05) 결과를,
출처 논문 및 문헌과 대조해 **"질병 특이적 생물학 vs 일반적 혈액/혈소판/리보솜 배경"**으로 분류 평가.

> **공통 관찰**: 다수 질병에서 ORA가 리보솜/번역(RPL·RPS), 미토콘드리아 OXPHOS, 프로테아좀/유비퀴틴(PSM·UBC·UBA52),
> 히스톤 클러스터(H2AC·H2BC·H3C·H4C), 혈소판 활성, 적혈구/헤모글로빈 등 **혈장 cfRNA 일반 배경 신호**에 지배됨.
> 이들 허브 유전자는 ORA에서 무관한 term(신경퇴행 KEGG, 바이러스 번역, SLE/NET, "Alcoholism" 등)으로 잘못 매핑되는 아티팩트를 만듦.

---

## 종합 판정표

| 질병 (n) | 출처 | 판정 | 핵심 |
|---|---|---|---|
| **Tuberculosis** (103) | Chang | **HIGH** | IFN-γ + 호중구/NET + GBP + BATF2 + 인플라마좀 = 교과서적 활성 TB |
| **HIV** (13) | Chang | **HIGH** | 순수 제1형 IFN/ISG (MX1, OAS, ISG15, IFI27, IFI44L) |
| **Liver Cancer** (48) | Roskams-Hieter+ | **HIGH** | 간세포 분비단백 (ALB, FGA/B/G, 보체, 아포지단백, 트랜스페린) = 조직 기원 |
| **HIV + Tuberculosis** (11) | Chang | **MOD–HIGH** | TB 면역+잔여 HIV ISG, 히스톤 세포주기 배경이 일부 희석 |
| **ICI-m** (11) | Raissadati | **MODERATE** | TCR 신호 + PD-1/PD-L1 체크포인트 (CD8 T세포 기전 일치); 심근 마커는 부재 |
| **Pancreatitis** (81) | Moore | **MODERATE** | 세포사멸(BAD/BID) + HLA/면역 = 선세포 손상·염증 일관 |
| **CDCS** (224) | Ward Z | **LOW–MOD** | 미토콘드리아+혈소판+죽상경화 신호 존재하나 리보솜 잡음에 묻힘 |
| **Lung Cancer** (33) | Chen | **LOW–MOD** | 세포주기/저산소 회수, 조직 마커(SFTPC/NKX2-1) 없음 |
| **Pancreatic Cancer** (74) | Moore | **LOW–MOD** | 피브리노겐-FN1 응고병증 + 간 분비축, 선세포 효소 마커 부재 |
| **Colorectal Cancer** (41) | Chen | **LOW** | 프로테아좀/리보솜 배경 지배, 대장 기원 신호 없음 |
| **Stomach Cancer** (29) | Chen | **LOW** | 소포수송/Rho-GTPase 약간 구별, 위점막 마커 없음 |
| **Liver Cirrhosis** (5) | Roskams-Hieter+ | **LOW** | 적혈구/헤모글로빈 배경, 간 신호 미약(GPX3/APOA4), n=5 취약 |
| **Pre-eclampsia** (62) | Moufarrej | **LOW** | ACE2/JAK2(RAS) 흥미롭지만 고전적 태반 sFlt1/PAPPA2 시그니처 부재 |
| **ME/CFS** (90) | Gardella | **LOW** | 혈소판+에스트로겐, 일반 cfRNA 배경과 구별 곤란 |
| **MM** (18) | Roskams-Hieter | **LOW** | 적혈구/세포주기 배경, Ig/형질세포/ER 신호 부재 |
| **Pancreatic Cancer** (6) | Reggiardo | **LOW/NONE** | 리보솜+혈소판; 논문 실제 신호는 Alu/반복서열(ORA 포착 불가)+n=6 |
| **Esophagus Cancer** (27) | Chen | **NONE/LOW** | OXPHOS·혈소판 배경만, 편평상피 신호 없음 |
| **Other Cancer** (18) | Moore | **NONE/LOW** | 순수 리보솜/번역 배경 (혼합 암종, 소표본) |
| **ICI-treated Cancer** (11) | Raissadati | **NONE/LOW** | 리보솜/Rho-GTPase/일반 배경, 체크포인트 신호 없음 (baseline) |
| **MGUS** (8) | Roskams-Hieter | **NONE/LOW** | 히스톤 클러스터 아티팩트(NET/SLE/Alcoholism), 형질세포 신호 전무 |

---

## 핵심 결론 (파이프라인 평가)

**1. 정규화 Z-score ORA가 진짜 생물학을 회수하는 두 가지 조건:**
- **(a) 강하고 일관된 전신 면역/혈액 전사 프로그램** — 감염질환(TB, HIV)·ICI 심근염의 면역 축. 가장 깨끗한 신호.
- **(b) 병변 조직이 혈장으로 대량 분비하는 조직 특이 전사체** — 간암(HCC)이 유일한 고형암 성공 사례(간세포 분비단백 ALB·피브리노겐·아포지단백).

**2. 대부분의 고형암·형질세포질환은 LOW/NONE.**
- Chen et al. 4개 GI암(대장·식도·폐·위)은 **서로 거의 구별되지 않고** 공통 리보솜/OXPHOS/프로테아좀 배경에 지배 → 암종 특이 조직 생물학 미회수. 원 논문의 "혈소판+리보솜 감소" 배경 소견과 일치.
- MGUS·MM은 면역글로불린/형질세포/ER-UPR 신호 전무 — 형질세포질환의 핵심 생물학 미포착.

**3. 주목할 "예상된 신호의 부재" (파이프라인 점검 필요):**
- **Pre-eclampsia**: cfRNA의 교과서적 성공 사례(Moufarrej, Nature)임에도 태반 항혈관신생 시그니처(sFlt1/FLT1, PAPPA2, ENG)가 거의 부재. gene selection 또는 정규화 모델이 태반 분획을 놓치고 있을 가능성 → 조사 가치 높음.

**4. 방법론적 경고 (ORA 아티팩트):**
- 리보솜(RPL/RPS), 히스톤 클러스터, 프로테아좀/UBC 허브 유전자가 무관 term으로 promiscuous 매핑(신경퇴행 KEGG, 바이러스 번역, SLE/NET, Alcoholism). → **ORA 전 이들 유전자 마스킹**, 또는 thresholded 리스트 ORA 대신 **전체 Z-ranking 기반 GSEA-prerank** 권장.
- 소표본(n≤8: MGUS·간경변·Reggiardo췌장암)은 통계적으로 취약, 1-2개 샘플의 고발현 전사체가 ORA 전체를 견인할 수 있음.
- gene-symbol ORA는 반복서열(Reggiardo Alu) 신호를 원천적으로 포착 불가.

**5. 총평:** 강한 결과(TB·HIV·Liver Cancer·ICI-m)는 정규화 Z-score 접근이 **실제 질병 생물학을 회수할 수 있음을 검증**. 실패 사례 대부분은 정규화 모델 자체보다 **cfRNA 조성(혈액세포 지배) + ORA 허브 유전자 아티팩트**의 산물. → 모델은 유효하되, 후처리(허브 유전자 제거, GSEA 전환)와 조직 분획 회수 개선이 다음 과제.

---

## 질병별 상세 평가

### 감염질환 (Chang et al., medRxiv 2023; 10.1101/2023.01.11.23284433~435)

**Tuberculosis (n=103) — HIGH.** 가장 깨끗한 활성 TB 시그니처: IFN-γ 우세(GBP1/2/4/5, STAT1, BATF2) + 호중구 탈과립/NET(DEFA4, AZU1, PRTN3, MPO, S100A8/9; OR~5) + 인플라마좀/pyroptosis(CASP1/4/5, GBP5) + 해당과정/PPP 대사전환(GPI, PKM, ENO1, G6PD). GBP·BATF2는 확립된 TB 진단 마커. 배경 희석 적음. (PMC3492754; 10.1038/s41467-018-04579-w)

**HIV (n=13) — HIGH.** 거의 순수한 제1형 IFN/ISG 프로그램(IFN-α/β, MX1/2, OAS1/2/3, ISG15, RSAD2, IFIT1/2/3, IFI27, IFI44L; adj-P~1e-16, OR>30) — 만성 HIV 바이러스혈증의 교과서적 혈액 시그니처. IFI27은 중증도 마커. 바이러스감염 KEGG도 같은 ISG 코어. (10.1172/JCI89488; PMC3577907)

**HIV + Tuberculosis (n=11) — MOD–HIGH.** 활성 TB 면역(IFN-γ + 호중구/NET + 보체 + 인플라마좀) + 잔여 HIV ISG의 가산적 중간 표현형. 단 가장 강하고 수많은 항목이 히스톤 클러스터 세포주기/유사분열 배경(H2BC/H4C 아티팩트)이라 신호가 부분 희석.

### 관상동맥질환 (Ward Z et al., Cells 2022; 10.3390/cells11203191)

**CDCS = Coronary Disease Cohort Study (관상동맥질환/ACS, n=224) — LOW–MOD.** 질병 관련: 미토콘드리아 OXPHOS(ATP5*, NDUFB10, COX4I1 — 원 논문의 미토콘드리아 mRNA 소견 일치), 혈소판 활성/지혈(RAP1B, MYL9, GNAS — ACS 핵심), Lipid & atherosclerosis(STAT3, CCL5, SOD2). 그러나 상위 ~20개가 전부 리보솜/번역 블록(adj-P 1e-25~1e-33)에 지배 + UBC/UBA52 허브 아티팩트 + SARS/Influenza 번역 + "Axon Guidance"(리보솜 아티팩트). CAD 생물학은 잡음 아래 소수 유전자가 견인. (PMID 31226691; 21048155; 19112112)

### GI 고형암 (Chen et al., eLife 2022; 10.7554/eLife.75181)

**Lung Cancer (n=33) — LOW–MOD.** 4개 암 중 가장 깊은 세포주기/증식(APC-C, G2/M, DNA복제) + HIF/저산소(EPAS1/HIF-2α; PMID 27932066). NOTCH/Hedgehog/WNT는 PSMB 견인 아티팩트. 폐 마커(SFTPC/NKX2-1) 없음, 프로테아좀 팽창이 신뢰도 저하.

**Colorectal Cancer (n=41) — LOW.** 리보솜+OXPHOS+프로테아좀 배경 지배. 세포주기/저산소(HIF-1α/EPAS1; PMID 31412254) 약하게 존재하나 공유 프로테아좀 견인, 대장 마커(CDX2) 없음.

**Stomach Cancer (n=29) — LOW.** 소포/Golgi 수송·Rho-GTPase 강조로 다소 구별, WNT/TCF(CDC73, ASH2L; PMID 32393771) 약간. 위점막 마커(MUC/TFF/GKN/PGC) 없음, 유의성 약함(상위 adj-P~1e-9).

**Esophagus Cancer (n=27) — NONE/LOW.** 55개 항목, OXPHOS+신경퇴행 KEGG(아티팩트)·혈소판·멜라노좀 배경 지배. 편평상피/케라틴/EMT/ESCC 마커 없음.

> **Cross-cancer**: 4개 프로파일 대체로 공유(구별 안 됨). 공통 코어=리보솜+OXPHOS+프로테아좀, "암 신호전달" term들은 공유 PSMB/RBX1의 ORA 다중매핑 아티팩트. 고전적 조직 기원 상피 마커·EMT/ECM 모듈 전부 미회수.

### 간/혈액암 (Roskams-Hieter B et al., npj Precis Oncol 2022; 10.1038/s41698-022-00270-y)

**Liver Cancer (HCC, n=48) — HIGH.** 간세포 분비 혈장단백이 압도: 알부민(ALB), 피브리노겐(FGA/B/G), 응고/항응고(PLG, KNG1, APOH, A2M), 보체(C3, CFB), 아포지단백(APOA2/B/C1/E, ANGPTL3), 운반단백(TF, CP, HPX, RBP4, AMBP). 원 논문의 "간 특이 유전자" HCC 바이오마커 및 독립 cfRNA deconvolution과 일치. 혈소판/ECM은 배경이나 견인 유전자가 간 기원 → 진짜 조직 기원 신호. (PMID 35468987; 10.1186/s40246-023-00537-w)

**MM (다발골수종, n=18) — LOW.** 적혈구/헤모글로빈(HBA1/HBD/HBG, ALAS2, FECH)+세포주기(TOP2A, MKI67, CENPE/F)+비특이 유비퀴틴/자가포식 배경. **Ig/형질세포/ER-UPR(XBP1/HSPA5) 신호 전무**. 세포주기만 원 논문의 "골수·세포주기" MM 마커와 약하게 일치.

**Liver Cirrhosis (n=5) — LOW.** 적혈구/헤모글로빈(HBB/HBA1/HBD, SLC4A1, CA1, ALAS2, KLF1) 산소운반 배경 + 히스톤. 간 특이 신호는 GPX3·APOA4 정도로 미약. **n=5 극도로 취약.**

**MGUS (n=8) — NONE/LOW.** 4개 모두 히스톤 클러스터(H2AC/H2BC/H3C/H4C) 견인 — NET formation/SLE/Alcoholism/Viral carcinogenesis는 전형적 히스톤 아티팩트. 형질세포/면역글로불린 신호 전무. **n=8 취약.**

### 췌장 (Moore et al., Nat Commun 2025; 10.1038/s41467-025-62685-y)

**Pancreatitis (n=81) — MODERATE.** 세포사멸/내인성 BH3(BAD, BID, APIP — 선세포 사멸/괴사는 췌장염 핵심; PMID 19560599) + HLA class II(HLA-DPB1/DRB1)/IL 신호(염증). 가장 일관된 질병 신호. 단 선세포 효소 마커 부재, RAS-MAPK 증식 블록(NRAS/CCND1)이 다소 암-유사(배경/교차오염 가능).

**Pancreatic Cancer (PDAC, n=74) — LOW–MOD.** 피브리노겐-FN1-인테그린 응고병증(FGA/FGB/VWF/FN1; PDAC 전혈전성, PMID 28765272) + ALB/APOB/IGFBP 간 분비축(원 논문 SERPINA1/CRP 간 우세 PDAC 신호와 일치). **선세포 효소 마커(PRSS/CELA/CTRB/CPA/AMY/PNLIP) 전부 부재.** 강한 항목은 혈소판/응고·수송·OXPHOS 배경.

**Other Cancer (n=18, 혼합) — NONE/LOW.** 상위 ~30개가 전부 리보솜/번역 배경. KLK3/GSTP1/ERG 산발적이나 일관된 조직 기원 없음. 소표본+혼합 암종.

> **PDAC vs 췌장염**: 둘 다 선세포 효소 신호 부재(공유 췌장 기원 없음). 방향성은 기대대로 분기 — 췌장염=세포사멸/HLA 염증, PDAC=응고병증/간 분비축. 단 췌장염의 RAS-MAPK 증식 블록이 구분을 흐림.

### ICI 심근염 (Raissadati et al., JCI 2024/25; 10.1172/jci188817)

**ICI-m = ICI 유발 심근염 (n=11) — MODERATE.** TCR 신호(NFATC2, PTPRC, CARD11) + **PD-L1/PD-1 체크포인트 경로** + BCR/MHC-II — 논문의 CD8+ Temra T세포 매개 기전과 일치, irAE에 기전적으로 on-target. 단 심근 손상 마커(TNNT2/TNNI3/MYH6/전도) 부재 → 심장 손상 축은 미표현. 호중구/혈소판/ER/리소좀 배경 큼. (PMID 27806233)

**ICI-treated Cancer (n=11) — NONE/LOW.** 리보솜/번역 + Rho-GTPase 세포골격 + 일반 백혈구/혈소판/적혈구 배경. 체크포인트/심장 신호 없음 — "ICI 투여, irAE 없음" baseline과 부합.

### ME/CFS (Gardella et al., PNAS; 10.1073/pnas.2507345122)

**ME/CFS (n=90) — LOW.** 8개 모두 Reactome: 혈소판 활성/탈과립/지혈(LYN, ITGB3, ITGA2B, PLEK, SOD1) + 에스트로겐/핵수용체(ESR, IGF1R, MMP9). 혈소판은 ME/CFS microclot 문헌 일부 있으나 일반 cfRNA 배경과 구별 불가. 에스트로겐은 여성 우세와 연결 가능하나 추정적. ME/CFS는 검증된 분자 마커 부재로 해석성 낮음.

### Pre-eclampsia (Moufarrej et al., Nature 2022; 10.1038/s41586-022-04410-z)

**Pre-eclampsia (n=62) — LOW.** 3개만 유의: Cell Cycle(일반 증식) + ACE2 포함 "SARS therapeutics"/"Infectious Disease". **ACE2·JAK2는 PE의 RAS/혈압 기전과 진짜 관련**이나, 고전적 태반 항혈관신생 시그니처(sFlt1/FLT1, PAPPA2, ENG, INHBA) 부재 — cfRNA 교과서적 성공 사례임을 감안하면 주목할 누락.

### Pancreatic Cancer (Reggiardo et al., Nat Biomed Eng 2023; 10.1038/s41551-023-01081-7)

**Pancreatic Cancer (Reggiardo, n=6) — LOW/NONE.** 리보솜/번역 + 혈소판/세포골격 배경 지배, 췌장 조직 마커 부재. 핵심: 논문 실제 신호는 **Alu/반복서열 RNA**라 gene-symbol ORA로 원천 포착 불가 + **n=6 취약**. 가설 생성용으로만.
