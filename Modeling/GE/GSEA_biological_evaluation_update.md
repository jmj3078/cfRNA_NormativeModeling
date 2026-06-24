# GSEA Biological Evaluation — Study-Split Update
## Liver Cancer (Chen) vs Liver Cancer (Roskams-Hieter)

**Context:** 기존 Group2 평가(`GSEA_biological_evaluation_Group2.md`)는 세 study (Chen, Roskams-Hieter, Block)를 합산한 결과를 분석했다. Study-level split 이후 n≥5 기준을 충족하는 두 그룹에 대해 추가 평가를 수행한다.

| Study | n (OOD 후) | DOI |
|---|---|---|
| Chen et al. (eLife 2022) | 10 | https://doi.org/10.7554/elife.75181 |
| Roskams-Hieter et al. (npj Precis Oncol 2022) | 28 | https://doi.org/10.1038/s41698-022-00270-y |

---

## 1. Liver Cancer (Chen) — n=10

**코호트 특성:** Chen et al. eLife 2022는 중국 환자 코호트 기반의 다중 암종 cfRNA 연구. Liver Cancer 샘플은 HCC (간세포암) 환자들로 구성.

### NES>0 (고Z 유전자에 농축 — 정상 대비 상향 경로) : 171개

**최상위 신호 (NES ≥ 2.2):**

| Term | NES | DB | 핵심 lead genes |
|---|---|---|---|
| Complement and coagulation cascades | 2.476 | KEGG | FGB, FGA, FGG, PLG, C3, CFB, CFH, SERPINC1 |
| Post-translational Protein Phosphorylation | 2.426 | Reactome | ALB, FGA, APOB, KNG1, FGG, TF, AHSG |
| Regulation Of IGF Transport / IGFBPs | 2.347 | Reactome | ALB, FGA, APOB, KNG1, FGG, TF, AHSG |
| Glycosaminoglycan biosynthesis | 2.331 | KEGG | B4GALT1, GLCE, CHST11, CHPF, EXTL1 |
| Cholesterol Efflux | 2.326 | GO-BP | APOA2, NPC2, APOC3, APOE, APOA1, APOA4 |
| Plasma Lipoprotein Remodeling | 2.284 | Reactome | ALB, APOB, APOC3, APOE, APOA1, ANGPTL3/4 |
| Regulation Of Fibrinolysis | 2.281 | GO-BP | VTN, PLG, PLAU, SERPINF2, THBD, F2 |

**생물학적 해석:**

- **간 특이적 분비 단백질 대량 상향**: ALB(알부민), FGB/FGA/FGG(피브리노겐 삼량체), PLG(플라스미노겐), APOB/APOE/APOC1-3(지질단백질), C3/CFB/CFH(보체) — 이들은 간세포가 합성·분비하는 혈장 단백질들. HCC에서 정상 간세포 파괴와 염증 반응으로 이 단백질들의 cfRNA가 증가하는 것은 병태생리학적으로 타당.
- **보체/응고 cascade 상향**: 간은 응고인자와 보체의 주요 생산 기관. HCC에서 급성기 반응(acute phase response)이 활성화되어 FGB/FGA/FGG, PLG, SERPINC1 등이 증가. Roskams-Hieter 2022에서도 fibrinogen과 complement가 HCC 바이오마커로 보고됨.
- **콜레스테롤/지질단백질 대사 상향**: APOA1/2, APOE, APOC1-3 — HCC에서 지질 대사 재프로그래밍이 특징적으로 나타나며, APOE와 APOC3는 HCC 고위험 마커로 문헌에 보고됨.
- **IGF transport**: IGFBP 조절이 상향됨 — HCC에서 IGF-2 과발현 및 IGFBP 축 활성화가 알려진 기전.
- **Fibrinolysis/응고 양방향 조절**: PLG(플라스미노겐), PLAU(uPA), SERPINF2(α2-antiplasmin), THBD(트롬보모듈린) — HCC에서 응고 항진 및 혈전성 변화가 일반적.
- **HNF4A/Maturity onset diabetes 경로**: HNF4A, FOXA2, NKX6-1 등 — HCC 발생의 주요 전사인자 축. 간 분화 경로와 직결.

**cfRNA 아티팩트 여부:** 히스톤(H2BC/H2AC) 관련 term 부재. 주요 lead gene들이 모두 간 특이적 단백질 → **높은 질병 특이성**.

### NES<0 (저Z 유전자에 농축 — 정상 대비 하향 경로) : 454개

**최상위 신호 (NES ≤ -3.0):**
- 리보솜/번역 기계 전체 (NES 최저 -3.356): Eukaryotic Translation Elongation, Peptide Chain Elongation, rRNA Processing, Cap-dependent Translation Initiation, Cytoplasmic Translation — RPL36, RPL21, RPS13, RPS25, EIF4E, EIF3J 등 거의 전체 번역 기계.
- **Cellular Response To Starvation** (NES=-2.715): GCN2/ISR(통합 스트레스 반응) 활성화 시 나타나는 번역 억제 패턴.
- **SLIT/ROBO Expression Regulation** (NES=-2.693): PSMC5/6, 프로테아좀 서브유닛들.

**생물학적 해석:**

- **번역 기계 전면 하향**: Chen 코호트에서 가장 강력한 신호. 리보솜 단백질과 번역 인자들이 정상 HC 대비 일관되게 낮은 Z-score. 이는 HCC 세포의 Warburg 효과 및 번역 수준 재프로그래밍과 일치하며, Group2 평가의 "간암 = OXPHOS 하향, GI암 = 번역 하향" 패턴에서 Chen 코호트가 오히려 GI암과 유사한 번역 억제를 보이는 점은 흥미로운 발견. 이는 Chen 코호트의 샘플이 HCC 외에 다른 소화기계 암을 포함할 수 있음을 시사하거나, 혹은 소수 샘플(n=10)로 인한 편향일 수 있음.
- **Starvation response**: GCN2 경로 활성화는 종양 미세환경의 영양 결핍(아미노산 고갈)을 반영.

**신호 품질 평가:** MODERATE-HIGH. 상향 경로는 간 특이적 분비 단백질로 매우 명확하나, n=10의 소수 샘플로 인해 하향 패턴의 해석에 주의 필요.

---

## 2. Liver Cancer (Roskams-Hieter) — n=28

**코호트 특성:** Roskams-Hieter et al. npj Precision Oncology 2022는 미국 환자 코호트 기반의 HCC/간경화/MM/MGUS cfRNA 연구. Liver Cancer 샘플은 HCC 환자 중심.

### NES>0 (고Z 유전자에 농축) : 211개

**최상위 신호:**

| Term | NES | DB | 핵심 lead genes |
|---|---|---|---|
| Negative Regulation Of Blood Coagulation | 2.224 | GO-BP | FGB, APOH, FGA, PLG, KNG1, FGG, APOE, VTN |
| Complement and coagulation cascades | 2.110 | KEGG | FGB, FGA, PLG, KNG1, FGG, F13A1, A2M, CFB, CFH, C9 |
| Adherens junction | 2.062 | KEGG | TJP1, CTNND1, WASF2, RAC2, AFDN, VCL, IQGAP1 |
| Positive Regulation Of G1/S Transition | 1.998 | GO-BP | TFDP1, RRM2, CCND2, CDC6, RRM1 |
| Fibrinolysis | 1.975 | GO-BP | FGB, FGA, PLG, FGG, PLAT, KRT1 |
| Actin Filament Bundle Assembly | 1.968 | GO-BP | FSCN1, LIMA1, NEDD9, LCP1, EZR |
| Oncogene Induced Senescence | 1.910 | Reactome | MAPK1, CDKN2A, CDK6, RB1, TP53, MDM2, CDKN2B |
| CDC42 GTPase Cycle | 1.909 | Reactome | ARHGAP29/31, ECT2, DOCK9, DLC1, PREX2, IQGAP1 |
| RAC1 GTPase Cycle | 1.884 | Reactome | FERMT2, FGD5, WASF2, SWAP70, ECT2, ITGB1, ARHGAP21 |
| YAP1/RUNX3/TEAD4 Transcriptional Regulation | ~1.90 | Reactome | WWTR1, YAP1, TEAD4, LATS1, EP300, CDKN2A |

**생물학적 해석:**

- **보체/응고 상향**: Chen 코호트와 공통되는 핵심 간 기능 신호. FGB/FGA/FGG, PLG, KNG1, APOH, A2M 등 간 분비 단백질들. 이는 study 간 재현성이 높은 **HCC의 강력한 cfRNA 바이오마커 축**임을 시사.
- **세포골격/침습 신호**: Adherens junction (TJP1, CTNND1, VCL, AFDN), Actin bundle assembly (FSCN1, EZR, LCP1, NEDD9), CDC42/RAC1 GTPase cycles — HCC의 EMT(상피중간엽이행)와 침습성 증가를 반영하는 경로. Chen 코호트에서는 미미했던 신호가 Roskams-Hieter에서 더 뚜렷하게 나타남.
- **세포주기 활성화**: G1/S transition (TFDP1, RRM2, CCND2, CDC6), Oncogene-Induced Senescence (CDKN2A, CDKN2B, TP53, RB1, MDM2, CDK6) — Roskams-Hieter 코호트 HCC 샘플의 종양 세포 분열 활성을 반영. Roskams-Hieter 2022 원문에서 HCC와 건강인 대비 세포주기 유전자들이 주요 바이오마커로 보고됨.
- **Hippo/YAP**: YAP1, WWTR1(TAZ), TEAD4, LATS1, NF2 — HCC에서 Hippo 경로 비활성화와 YAP 핵 이동은 확립된 발암 기전. cfRNA에서 검출된다는 점이 중요.
- **내피세포 조절 경로**: Artery Development (NOS3, ENG, APLNR), Negative Regulation Of Endothelial Cell Migration (VASH1, THBS1, SLIT2, ADAMTS9), Endothelial Cell Development — HCC에서 종양 혈관신생(angiogenesis)의 반영 가능성.
- **항바이러스 신호**: Negative Regulation Of Viral Genome Replication (IFITM2/3, EIF2AK2, OAS2, IFI16) — HCC 환자 다수가 B/C형 간염 바이러스 감염 관련임을 반영. Roskams-Hieter 코호트의 바이러스성 간염 유병률을 암시.

### NES<0 (저Z 유전자에 농축) : 단 6개

| Term | NES | DB |
|---|---|---|
| Protein Insertion Into ER Membrane | -2.332 | GO-BP |
| Oxidative phosphorylation | -2.261 | KEGG |
| Quinone Biosynthetic Process | -2.236 | GO-BP |
| Mitochondrial Electron Transport | -2.227 | GO-BP |
| Lysosomal Lumen Acidification | -2.181 | GO-BP |
| Mitochondrial Respiratory Chain Complex Assembly | -2.108 | GO-BP |

**생물학적 해석:**

- **OXPHOS/미토콘드리아 대사 하향**: Group2 평가에서 확인된 "HCC = Warburg 효과, OXPHOS cfRNA 하향"이 여기서도 재현. UQCR11, NDUFA1, COX17, NDUFS5, NDUFB4 등 미토콘드리아 전자전달계 서브유닛들.
- **하향 경로 수 자체(6개)가 Chen(454개)과 극명한 대조**: Roskams-Hieter 코호트는 NES<0 경로가 거의 없음 → 이 코호트의 HCC 샘플들이 정상 HC 대비 주로 특이적 상향 신호를 보이며, 광범위한 하향 억제 패턴은 나타나지 않음. 이는 Roskams-Hieter 코호트의 HCC 샘플들이 더 초기 단계거나 간 기능이 비교적 보존된 환자들일 가능성을 시사.

**신호 품질 평가:** HIGH. n=28로 Chen보다 통계적으로 신뢰도가 높으며, 상향/하향 신호 모두 HCC 기존 문헌과 일관됨. 세포주기, YAP, GTPase, 보체/응고의 4가지 축이 HCC-specific 바이오마커 후보군으로 강력히 지지됨.

---

## 3. 두 Study 비교 분석

### 공통 신호 (양 코호트 재현)
- 보체/응고 cascade 상향 (FGB, FGA, FGG, PLG, KNG1 core)
- 지질단백질 재조합 / 콜레스테롤 대사 (APOE, APOB, APOA1, APOC1-3)
- Fibrinolysis (PLG, PLAU/PLAT, SERPINF2)
- OXPHOS 하향 (Warburg effect)

→ 이 4가지 축은 **study에 관계없이 재현되는 HCC cfRNA 코어 신호**로 볼 수 있다.

### Chen 특이적 신호
- 리보솜/번역 기계 대규모 하향 (454개 NES<0, NES 최저 -3.356) — n=10 소수 샘플의 편향 가능성 있으나, 번역 억제 패턴은 첸 코호트에서 독특하게 강함
- 간 분화 전사인자 축 상향 (HNF4A, FOXA2, NKX6-1): 간세포 정체성 관련
- Ion channel / neuron differentiation 경로 상향: 간암에서 일부 신경내분비 분화 특성 가능성

### Roskams-Hieter 특이적 신호
- 세포골격/침습 (Adherens junction, Actin bundle, FSCN1): EMT 및 침습 반영
- CDC42/RAC1/RHOU/RHOV GTPase cycles: 세포 이동 및 침습
- Oncogene-Induced Senescence (CDKN2A, TP53, MDM2): 종양 억제 회피 기전
- YAP/Hippo 상향 (YAP1, TEAD4, LATS1): HCC 발암 경로
- 항바이러스 신호 (IFITM2/3, OAS2): 바이러스성 간염 배경 암시

### 해석적 결론
두 코호트는 동일한 "Liver Cancer" 레이블이지만 상이한 cfRNA 프로파일을 보이며, 이는 **환자 모집단의 차이(인종, 바이러스성 vs 비바이러스성 HCC, 종양 스테이지)를 반영**할 가능성이 높다. Study-level 분리의 타당성이 생물학적으로도 지지된다.

---

## 결론 요약

| 항목 | Liver Cancer (Chen) | Liver Cancer (Roskams-Hieter) |
|---|---|---|
| n (OOD 후) | 10 | 28 |
| 신호 품질 | MODERATE-HIGH | HIGH |
| 핵심 상향 경로 | 간 분비 단백질, 보체, 지질단백질 | 보체, 세포주기, YAP, GTPase, 혈관신생 |
| 핵심 하향 경로 | 리보솜/번역 전면 억제 | OXPHOS (소수) |
| 고유 특성 | 번역 기계 하향 패턴이 독특 | 침습/세포주기/YAP 신호가 풍부 |
| 코호트 추정 특성 | 진행성 HCC, 번역 억제 심각 | 초기-중기 HCC, 바이러스성 간염 배경 가능 |
| cfRNA 아티팩트 | 없음 | 없음 |
| Study-split 타당성 | ✓ 생물학적으로 구별되는 프로파일 | ✓ |
