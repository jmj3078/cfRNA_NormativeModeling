# GSEA Master Report — Rework
## 기존 검증 signature의 cfRNA 재현성 평가 (Reproduction-of-validated-signatures view)

**작성일:** 2026-06-25
**관점 전환:** 기존 [GSEA_Master_Report.md](GSEA_Master_Report.md)는 "cfRNA에서 무엇이 보이는가"에 집중했다.
본 Rework는 질문을 뒤집는다 →

> **기존 데이터베이스·문헌에서 독립적으로 검출되고 재현·검증된 질병 signature가, cfRNA normative model에서 다시 나타나는가?**

재현(reproduction)이 잘 될수록 모델이 **진짜 질병 생물학**을 포착한다는 강한 증거다.
Novelty는 부차적이며, 여기서는 다음과 같이 **재정의**한다 →

> **Novelty(재정의)** = 기존 **비-cfRNA** 연구(조직·전혈·유전체)에서 잘 검증되었으나,
> **cfRNA(무세포 혈장) 구획에서는 처음으로 포착**된 검증 marker / pathway.

---

## 0. 방법론·해석 프레임 (critical-thinking 적용)

### 0.1 검증 기준 데이터원
- **Open Targets Platform** — 질병별 target-disease association (유전체·문헌·약물 통합 점수). 암종 driver 검증에 사용.
- **PubMed 검증 signature** — 분야 표준으로 다중 재현된 전사체/바이오마커 패널:
  - TB: Berry et al. *Nature* 2010 (PMID 20725040, IFN-inducible 호중구 전혈 signature); Sweeney et al. *Lancet Respir Med* 2016 (PMID 26907218, 3-gene GBP5/DUSP3/KLF2)
  - MM: Shaughnessy et al. *Blood* 2007 (PMID 17105813, GEP70 고위험 증식 signature)
  - Pre-eclampsia: Maynard et al. *J Clin Invest* 2003 (PMID 12618519, sFLT1 항혈관신생 축); 임상 표준 sFLT1/PGF ratio
  - HIV: 다수 검증된 type-I IFN/ISG signature
- 암 driver: COSMIC/Open Targets 정설 driver 유전자.

### 0.2 cfRNA에 대한 필수 해석 원칙 (confounder 통제)
1. **cfRNA = 전사체 abundance, NOT 변이.**
   따라서 *기능 상실/돌연변이 driver*(KRAS, TP53, APC, CDH1)는 cfRNA에서 직접 재현될 대상이 **아니다**.
   재현 대상은 **전사적으로 활성화되는 축** — RTK/신호전달(EGFR, MET, FGFR2), 분비 단백질(APO/FG/보체), 증식 기계(kinesin/E2F), 면역 신호(ISG/GBP)다.
2. **Leading-edge 부호 구분이 결정적.**
   - **NES > 0** 에 검증 marker가 있으면 → **질병 양성 신호로 재현됨**(진짜 신호).
   - **NES < 0** 에 있으면 → 대부분 **전신 번역억제 배경**에 섞여 들어간 것 → **특이적 재현 아님**(아래 한계 참조).
   많은 암에서 driver 유전자(TP53, APC, MSH6)가 NES<0로 잡히는데, 이는 그 유전자가 "질병 신호"여서가 아니라
   리보솜/번역 억제 leading-edge에 우연히 포함된 것이다. 본 보고서는 이를 **재현으로 인정하지 않는다.**
3. **검증 한계:** leading-edge 포함 ≠ 그 유전자 자체의 차등 풍부도. 경로 수준 신호의 proxy로 해석한다.

### 0.3 재현성 등급 (GRADE-style)
| 등급 | 정의 |
|---|---|
| **STRONG** | 검증된 핵심 signature 다수가 NES>0 양성 신호로 재현 |
| **PARTIAL** | RTK/신호 축 등 일부 검증 marker는 양성 재현, 나머지(특히 mutational driver)는 미재현 |
| **WEAK / DIVERGENT** | 임상·문헌 표준 marker가 부재하거나 배경 음성으로만 출현 → 모델 신호가 검증 축과 어긋남 |
| **UNVERIFIABLE** | 분야에 합의된 재현 gold-standard signature 자체가 없음 |

---

## 1. 재현성 종합 표 (전 표현형)

| 표현형 | 검증 signature (출처) | cfRNA 재현 | 등급 | redefined novelty (조직→cfRNA 확장) |
|---|---|---|---|---|
| **HIV** | type-I IFN/ISG (정설) | ISG15·OAS1-3·MX1·IFI27·IFIT1·RSAD2·USP18 등 **15/15** NES≈+2.7 | **STRONG** | 검증된 ISG 축을 무세포 혈장에서 완전 재현 |
| **Tuberculosis** | IFN-inducible 호중구 signature (Berry 2010; Sweeney3 2016) | GBP1/2/5/6·STAT1·IFIT1/3·IFITM3·FCGR1A·BATF2·DUSP3·SOCS3 등 **14/18** NES≈+2.2~2.5 | **STRONG** | Berry/Sweeney 전혈 signature를 cfRNA에서 재현; Sweeney3의 GBP5·DUSP3 모두 양성 |
| **HIV + TB** | ISG(HIV)+GBP/IFN(TB) 합 | 양 축 모두 양성 + B세포 억제 | **STRONG** | 두 검증 축의 중첩을 단일 cfRNA에서 분해 |
| **Liver Cancer (Chen)** | 간 분비단백·HNF4A 축 (HCC cfRNA 문헌) | APOB·APOE·FGA·FGB·ALB·SERPINA1·**HNF4A**·FABP1 NES≈+2.4 | **STRONG** | HNF4A 마스터 TF·간 분비체를 cfRNA에서 재현 |
| **Liver Cancer (Roskams-Hieter)** | HCC driver (OT: CTNNB1/MET/TP53/VEGFR) | **CTNNB1·MET·TP53·APC·FLT4·KDR·FLT1·BRAF·PIK3CA·RB1·CDKN2A 등 16/30** NES>0 | **STRONG** | Wnt(CTNNB1)·MET·VEGFR 조직 driver를 cfRNA에서 양성 포착 |
| **Esophagus Cancer (Chen)** | RTK/TGF-β driver (OT) | **FGFR2(+2.67)·TGFBR2(+2.67)·FGFR1·ERBB4·CREBBP·PTPRT** 양성 | **STRONG** | FGFR2/TGFBR2 축을 cfRNA에서 재현 |
| **Multiple Myeloma** | GEP70 증식 signature (Shaughnessy 2007) | CCND1·CCND2·MAF·NSD2·CKS1B·NEK2·CENPE·NUSAP1·TOP2A·MKI67·BIRC5·KIF11 **12/14** NES>0 | **STRONG** | 골수 GEP70 고위험 증식 signature를 **비침습 혈장**에서 재현 |
| **Lung Cancer** | NSCLC RTK driver (OT) | **EGFR(+2.32)·MET(+2.12)·KDR(+2.46)** 양성; TP53·KEAP1·EML4는 음성 배경 | **PARTIAL** | EGFR/MET/VEGFR 축 cfRNA 재현; ALK/KRAS/RET 미포착 |
| **Pancreatic Cancer (Moore)** | PDAC driver (OT) | **SMAD4·TGFBR2·BRCA1/2·GNAS** 양성; **KRAS 부재·TP53 부재·CDKN2A 음성** | **PARTIAL** | TGF-β/SMAD4·DDR 축 재현; 핵심 mutational driver(KRAS/TP53) 미재현 |
| **MGUS** | t(11;14) CCND1 / t(4;14) NSD2 | CCND1(+2.07)·NSD2(+1.89) 양성 | **PARTIAL** | 전위 표적을 전악성 단계 cfRNA에서 포착 |
| **ICI-m (myocarditis)** | 심근 자가항원 + 항원제시 (자가면역 심근염 정설) | **MYH7(+2.36, 자가항원)·HLA-A·TAP1·PRF1·IFITM2·NPPA** 양성; CD8A·GZMB·TNNT2 미포착 | **PARTIAL** | α/β-myosin 자가항원·MHC-I 제시를 혈장 cfRNA에서 포착 |
| **ICI-treated Cancer** | ICI 면역활성 (리소좀/항원제시/이식거부 패턴) | 리소좀·MHC-I·alloimmune 양성 | **PARTIAL** | ICI 전신 면역활성의 cfRNA 반영 |
| **Colorectal Cancer** | MMR/APC/KRAS driver (OT) | driver 전부 NES<0 배경; **양성 신호는 Hippo/YAP·ERBB2 경로** | **PARTIAL** | mutational driver 미재현; YAP/RTK 경로 수준만 재현 |
| **Stomach Cancer** | CDH1/RHOA/ERBB2/RET (OT) | driver NES<0 배경; **양성은 데스모좀 junction·RET** | **PARTIAL** | CDH1 미재현; junction/RET 경로 재현 |
| **Other Cancer** | 범암 침습 (Rho GTPase) | RHO/RAC1/CDC42 양성 | **PARTIAL** | 이질적 암 집합의 공통 침습 축 |
| **Pancreatitis** | 췌장염 위험유전자 PRSS1/SPINK1/CFTR/CTRC/CPA1 | 유전자 수준 **전부 부재**; 단 "Defective CFTR" **경로**는 음성 농축 | **WEAK** | 위험 유전자 미재현; CFTR 경로만 간접 시사 |
| **Pre-eclampsia** | **sFLT1(FLT1)·PGF·ENG** 항혈관신생 (Maynard 2003; 임상 sFLT1/PGF) | **FLT1·PGF·ENG·LEP 전부 부재**; INHBA·INHA·VEGFA는 음성 | **DIVERGENT** | ⚠ 임상 gold-standard 미재현 — 아래 critical note |
| **CAD_HF+** (HF 진행) | HF 진행 바이오마커 (FGF23·Klotho·THBS1·심장 섬유화·SERCA 칼슘) | **FGF23·KL·FGFR1 6회 반복(+2.82)**; COL1A1/A2 섬유화; F2R·THBS1 응고; SR 칼슘 | **PARTIAL** | HF 진행 축 cfRNA 재현 ✦ 단 고전 지질/염증 CAD 축(LPA·PCSK9·IL6) 부재 + 원본 Ward 2022 HF+/− null → 검증 필요 |
| **CAD_HF−** (HF 비진행) | CAD 유전·염증 (LPA·PCSK9·LDLR·IL6·CRP) | LPA·PCSK9·IL6·CRP 부재; FGFR1/섬유화 축도 부재 | **DIVERGENT** | 검증 CAD 축 미재현 + HF 진행 축 음성 = 비특이 전신억제 (HF+ 음성 대조군 역할) |
| **ME/CFS** | 합의된 재현 signature **부재** | 후보 대사 유전자(PDK4·FOXO3·MTOR·TNF) 음성 배경만 | **UNVERIFIABLE** | 검증 기준 자체 부재 → 재현성 판정 불가 |

> n<5 (Liver Cirrhosis n=4, Pancreatic Cancer Reggiardo n=6)는 분석 제외.

---

## 2. STRONG 재현 — 검증 signature가 cfRNA에서 살아남은 케이스

### 2.1 HIV — type-I IFN/ISG signature (15/15)
**검증 출처:** 다수 전혈·PBMC 연구에서 재현된 HIV ISG 활성 (IFI27, ISG15, OAS, MX1, RSAD2).
**cfRNA 재현:** 패널 15개 전부 NES +1.9~+2.74 양성. ISG15·OAS1-3·MX1·MX2·IFIT1·USP18·HERC5 모두 최상위.
**판정:** 검증된 핵심 면역 signature를 무세포 혈장에서 **완전 재현**. 모델 타당성의 강력한 양성 대조군.
**redefined novelty:** 세포 기반에서 검증된 ISG 축이 cfRNA 구획에서도 동일하게 포착됨.

### 2.2 Tuberculosis — Berry 2010 / Sweeney3 IFN-호중구 signature (14/18)
**검증 출처:** Berry et al. *Nature* 2010 (전혈 IFN-inducible 호중구 driven signature, 다국가 재현); Sweeney et al. 2016 3-gene(GBP5·DUSP3·KLF2).
**cfRNA 재현:** GBP1/2/5/6·STAT1·IFIT1/3·IFITM3·FCGR1A·BATF2·SERPING1·DUSP3·SEPTIN4·SOCS3 = NES +1.8~+2.5 양성.
**Sweeney3 중 GBP5·DUSP3 모두 양성** 재현(KLF2만 미검출).
**판정:** TB 분야 최강 검증 signature의 거의 전부를 cfRNA에서 재현 → **STRONG**.
**critical note:** IFI44L/ANKRD22 미검출은 cfRNA 검출 한계(저발현 ISG) 가능성.
**redefined novelty:** 전혈 진단 signature(Berry/Sweeney)를 **무세포 혈장**으로 확장 검출.

### 2.3 Multiple Myeloma — GEP70 증식 signature (12/14)
**검증 출처:** Shaughnessy et al. *Blood* 2007 — 골수 형질세포 GEP70 고위험 signature(증식/유사분열 유전자 핵심).
**cfRNA 재현:** CCND1·CCND2·MAF·NSD2·CKS1B·NEK2·CENPE·NUSAP1·TOP2A·MKI67·BIRC5·KIF11 = NES +1.7~+2.2 양성.
**판정:** 침습적 **골수 생검** 기반 검증 signature를 **혈장 cfRNA**에서 재현 → **STRONG** + 임상적 가치(비침습화).
**redefined novelty:** GEP70을 비침습 혈장으로 옮긴 첫 신호.

### 2.4 Liver Cancer (Roskams-Hieter) — HCC driver 16/30 (Open Targets)
**검증 출처:** Open Targets HCC top target — CTNNB1(Wnt), MET, TP53, VEGFR(FLT4/KDR/FLT1), BRAF, PIK3CA, RB1, CDKN2A.
**cfRNA 재현:** 위 16개 NES +1.7~+2.1 **양성**. 특히 **CTNNB1(+2.06)·MET(+2.06)** = HCC 정설 driver.
**판정:** 조직 driver 절반 이상을 양성 재현(이 코호트는 NES<0가 6개뿐이라 driver가 배경에 묻히지 않음) → **STRONG**.

### 2.5 Liver Cancer (Chen) — 간 분비체 + HNF4A
**cfRNA 재현:** APOB·APOE·FGA·FGB·ALB·SERPINA1·HNF4A·FABP1 NES≈+2.4 양성. **HNF4A**(간세포 정체성 마스터 TF) 재현이 핵심.
**판정:** 검증된 간기능/HCC cfRNA 바이오마커 축을 재현 → **STRONG**. (단 mutational driver는 부재)

### 2.6 Esophagus Cancer (Chen) — FGFR2/TGFBR2 RTK 축
**cfRNA 재현:** FGFR2(+2.67)·TGFBR2(+2.67)·FGFR1(+2.21)·ERBB4(+2.32)·CREBBP·PTPRT 양성.
**판정:** Open Targets 식도암 RTK/TGF-β driver를 양성 재현 → **STRONG**.

---

## 3. PARTIAL — RTK/신호 축은 재현, mutational driver는 미재현

### 핵심 패턴 (cancer 공통)
cfRNA는 **전사적으로 활성·분비되는 driver(EGFR, MET, FGFR2, SMAD4, VEGFR)** 는 양성 재현하지만,
**기능상실/돌연변이 driver(KRAS, TP53, APC, CDH1, CDKN2A)** 는 부재하거나 번역억제 배경(NES<0)에만 출현한다.
→ 이는 **결함이 아니라 cfRNA의 본질**(abundance ≠ mutation)이다. 변이 검출은 cfDNA의 영역.

- **Lung Cancer:** EGFR(+2.32)·MET(+2.12)·KDR(+2.46) 재현; KRAS·ALK·RET·TP53 미재현/배경. → NSCLC RTK 축 재현.
- **Pancreatic Cancer (Moore):** SMAD4·TGFBR2·BRCA1/2·GNAS 양성; **KRAS·TP53 미재현**(PDAC의 >90% KRAS 변이는 cfRNA 전사체로 잡히지 않음). TCR/T세포 억제(novel)와 함께 해석.
- **Colorectal Cancer:** APC·TP53·MMR(MSH6/MLH1) 전부 NES<0 배경. 양성 재현은 **Hippo/YAP·ERBB2 경로 수준**. → driver 유전자 미재현, 경로 재현.
- **Stomach Cancer:** CDH1(검증 위암 driver) 미재현; 양성은 데스모좀 junction(PKP/DSG2)·RET. → 경로 수준 재현.
- **MGUS:** CCND1·NSD2(전위 표적) 양성 → 전악성 단계 신호 포착.
- **ICI-m:** **MYH7(+2.36)** = 자가면역 심근염의 검증 자가항원(α/β-myosin), HLA-A·TAP1(MHC-I 제시)·PRF1(perforin)·NPPA(심근 stress) 양성. CD8A·GZMB·troponin(TNNT2/TNNI3) 미포착 → T세포 effector는 cfRNA에 덜 반영. → 항원제시 축 재현.

---

## 4. WEAK / DIVERGENT — 검증 축과 어긋나는 케이스 (⚠ 가장 중요한 비판적 발견)

### 4.1 Pre-eclampsia — 임상 gold-standard 미재현 ⚠
**검증 출처:** Maynard 2003 (sFLT1 항혈관신생); **임상에서 실제 사용되는 sFLT1/PGF ratio**, 그리고 sENG(ENG).
**cfRNA 결과:** **FLT1·PGF·ENG·LEP·PAPPA·ADAM12 전부 부재**. INHBA·INHA·VEGFA는 오히려 **음성(NES<0)**.
모델의 PE 양성 신호는 Neurexin/synapse(neuro 축, Moufarrej 2022와 일치)이나, **임상 검증 항혈관신생 축과는 어긋남**.
**판정:** **DIVERGENT.** cfRNA normative model이 PE의 **임상적으로 검증된 핵심 바이오마커를 재현하지 못함.**
**critical-thinking:**
- (a) HC reference가 **비임신** 대조군 → 정상 임신조차 항혈관신생/태반 유전자가 상향이므로, PE를 비임신 HC와 비교하면 이 축이 상쇄/희석될 수 있음(**construct validity 위협, confounding by pregnancy state**).
- (b) 태반 유래 cfRNA의 검출 효율·정규화 문제 가능.
- 결론: PE에 대해서는 본 모델 신호를 임상 바이오마커 대체물로 해석하면 안 됨. neuro 신호의 의미는 추가 검증 필요.

### 4.2 CAD (CAD_HF+ / CAD_HF−) — HF 진행 축은 재현, 고전 CAD 축은 미재현 ✦⚠

> **표현형 재정의.** 기존 "CDCS"는 코호트 이름(Coronary Disease Cohort Study)으로 잘못 라벨링된 것이었다.
> 원본 연구 **Ward et al., Cells 2022**는 안정형 CAD 환자를 **3년 내 HF 발병 여부**로 나눈 코호트이므로,
> **CAD_HF+ (HF 진행, n=116)** / **CAD_HF− (HF 비진행, n=108)**로 분리해 재해석한다.

**검증 출처 (두 축):**
1. *고전 CAD 축* — 유전(LPA, PCSK9, LDLR, SORT1) + 염증(IL6, CRP). GWAS·임상 다중 검증.
2. *HF 진행 축* — **FGF23–Klotho–FGFR1**(incident HF·심비대 예측, Klotho 혈장농도 ACS 사망 독립 예측), **THBS1**(HFpEF 심장 리모델링 악화), 심장 섬유화 콜라겐, 칼슘 핸들링.

**cfRNA 결과:**
- *고전 CAD 축:* 양 군 모두 LPA·PCSK9·IL6·CRP·APOB **부재**, LDLR 음성 → **미재현**.
- *HF 진행 축:* **CAD_HF+에서 FGF/FGFR1 축이 6개 경로 반복**(FGFR1 +2.82 전체 1위, lead FGF23·KL), 콜라겐 섬유화(COL1A1/A2 +2.45), 칼슘·이온채널 리모델링(SR Ca²⁺ +1.95), 전응고(F2R·THBS1) → **HF+ 특이 양성**.
- *CAD_HF−:* 상향 11개로 빈약, FGFR1/섬유화 축 **부재** → HF 진행 서명의 **자연적 음성 대조군**.

**판정:**
- **CAD_HF+ → PARTIAL.** 고전 지질/염증 CAD 축은 미재현이나, **독립 검증된 HF 진행 바이오마커(FGF23-Klotho/FGFR1·섬유화·이온채널)가 cfRNA에서 신규 포착**됨 — 사용자 정의의 "novelty"(비-cfRNA에서 검증, cfRNA에서 신규)에 부합.
- **CAD_HF− → DIVERGENT.** 검증 CAD 축 미재현 + HF 진행 축 음성 = 신호가 비특이 전신억제.

**critical-thinking (⚠ 결정적 단서):**
- **원본 Ward et al. 2022은 plasma RNA-Seq에서 HF+/HF− 차이를 전혀 검출하지 못했다(null).** 동일 코호트에서 gold-standard 분석이 null인데 본 정규화 모델이 강한 차이를 보이는 불일치는 두 가지로 갈린다:
  - (i) 공변량 보정 GAMLSS Z-score의 **민감도 이득** — 단순 fold-change가 놓친 약한 HF 진행 신호 포착.
  - (ii) **인공산물** — 배치/세포조성/라이브러리 복잡도 차이가 군 간 상향 경로 수(36 vs 11)를 만들었을 가능성.
- cfRNA는 **존재량(abundance)**이지 심근 조직 발현 방향이 아니므로, "SERCA 칼슘 상향" 등은 순환 단편 수준 농축으로만 해석.
- **결론:** HF 진행 축 재현은 생물학적으로 정합적이나 **가설 생성적**. 독립 코호트·직교 측정(FGF23/NT-proBNP 단백질, qPCR)으로 검증되기 전까지 진단 marker로 단정 불가. PARTIAL 판정은 **잠정적**이며 §6의 indirectness/imprecision 다운그레이드가 강하게 적용됨.

### 4.3 Pancreatitis — 위험 유전자 미재현 (경로만 간접)
**검증 출처:** 췌장염 위험 유전자 PRSS1·SPINK1·CFTR·CTRC·CPA1 (유전체 정설).
**cfRNA 결과:** 유전자 수준 **전부 부재**. 단 "Defective CFTR Causes Cystic Fibrosis" **경로**가 음성 농축(NES −2.22)으로 CFTR 축을 간접 시사.
**판정:** **WEAK.** 위험 유전자 자체는 미재현; 경로 수준 신호만 존재.

### 4.4 ME/CFS — 검증 gold-standard 부재 (UNVERIFIABLE)
**critical-thinking:** ME/CFS는 분야에 **합의된 재현 가능한 분자 signature가 없다**(진단은 임상 기준). Gardella 2025 cfRNA 연구 자체가 후보 단계.
후보 대사 축(PDK4·FOXO3·MTOR·TNF)은 cfRNA에서 음성 배경으로만 출현.
→ "재현 여부"를 **판정할 기준이 없음.** 모델의 번역·FoxO 억제 신호는 흥미롭되, 검증된 외부 signature로 뒷받침되지 않으므로 **확증 아닌 탐색적(exploratory)** 으로만 해석해야 함.

---

## 5. 종합 — 무엇이 재현되고 무엇이 안 되는가

**잘 재현됨 (모델이 진짜 질병 생물학 포착):**
- 면역 활성 signature — HIV ISG, TB IFN-GBP (분야 표준 전혈 signature를 cfRNA에서 재현) ★ 가장 강력
- 증식/유사분열 signature — MM GEP70 (골수→혈장 비침습 재현)
- **전사적으로 활성인** 암 driver 축 — HCC(CTNNB1/MET/VEGFR), 식도암(FGFR2/TGFBR2), 폐암(EGFR/MET), 간 분비체(HNF4A/APO/FG)
- 자가면역 자가항원 — ICI-m의 MYH7 + MHC-I 제시

**재현 안 됨 (해석상 한계 또는 모델 한계):**
- **mutational/LoF driver** — KRAS, TP53, APC, CDH1, CDKN2A → cfRNA(abundance)의 본질적 사각지대 (cfDNA 영역)
- **혈관신생 바이오마커(PE)** — sFLT1/PGF/ENG → 비임신 HC reference confounding 강력 의심 ⚠
- **고전 지질/염증 CAD 축** — LPA/PCSK9/IL6 → 양 군 모두 미재현 (단 CAD_HF+의 HF 진행 축 FGF23/FGFR1은 부분 재현, §4.2)
- **췌장염 위험 유전자** — PRSS1/SPINK1 등 → 유전자 수준 미재현

**핵심 메시지:**
모델은 **면역·증식·RTK/분비 축**에서 검증 signature를 강하게 재현하나,
**변이 기반·혈관신생·지질대사 축**에서는 재현이 약하다.
전자는 모델 타당성의 강한 증거, 후자는 (i) cfRNA의 본질적 한계 또는 (ii) reference/confounding 문제를 구분해 해석해야 한다.

---

## 6. Critical-thinking 한계 (GRADE downgrade 요인)

1. **번역억제 배경의 교란.** 다수 질병에서 NES<0 leading-edge가 거대(예: CAD_HF+ 925개, ME/CFS 169개). 검증 driver가 이 배경에 우연히 포함되면 "재현"으로 오인될 수 있어, 본 보고서는 **NES>0만 재현으로 인정**했다. 그래도 일부 양성 재현이 단일 광역 신호의 부산물일 가능성 잔존(**indirectness**).
2. **Leading-edge ≠ 차등 풍부도.** GSEA leading-edge 멤버십은 경로 수준 proxy다. 개별 유전자 검증은 per-gene Z 분포로 별도 확인 필요(시각화의 우측 패널이 일부 보강).
3. **비임신·비암 HC reference.** PE(임신)·일부 암에서 reference mismatch가 construct validity를 위협(특히 PE에서 결정적).
4. **소표본.** ICI-m/MGUS(n≤11)·HIV(n=13)는 양성 재현도 outlier 구동 가능(**imprecision**). HIV의 ISG 재현은 효과크기가 커서 견고하나, MGUS 등은 잠정.
5. **검증 패널 선정 편향.** 본 cross-ref 패널은 대표 검증 marker를 선별한 것으로 완전하지 않다(**selection**). 누락 marker가 결론을 바꿀 수 있다.
6. **다중비교.** GSEA permutation(n=100) + FDR<0.05; 일부 phenotype의 100% 유의는 광역 분포 이동을 반영하므로 개별 경로 특이성은 보수적으로 해석.
7. **재현 ≠ 인과·진단 성능.** 재현은 생물학적 타당성의 증거일 뿐, 진단 AUC/임상 유용성은 별도 평가 필요(특히 DIVERGENT 케이스).

---

## 7. 출처
- Open Targets Platform — target-disease associations (api.platform.opentargets.org, 2026-06 조회)
- Berry MPR et al. *Nature* 2010;466:973–977. PMID 20725040 (TB 전혈 IFN signature)
- Sweeney TE et al. *Lancet Respir Med* 2016;4:213–224. PMID 26907218 (TB 3-gene GBP5/DUSP3/KLF2)
- Shaughnessy JD et al. *Blood* 2007;109:2276–2284. PMID 17105813 (MM GEP70)
- Maynard SE et al. *J Clin Invest* 2003;111:649–658. PMID 12618519 (sFLT1 preeclampsia)
- 1차 cfRNA 데이터 출처는 [GSEA_Master_Report.md](GSEA_Master_Report.md) §출처 참조.

---

*본 Rework는 cfRNA normative model GSEA 결과를 기존 검증 signature 재현성 관점에서 재해석한 것이며, 등급은 NES 부호·검증 패널 일치도·표본 크기를 종합한 정성 판정이다.*
