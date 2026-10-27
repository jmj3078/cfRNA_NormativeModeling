# Disease tissues of interest for NM tissue-signal check (literature survey, 2026-10-27)

Purpose: for each scored phenotype, list the non-blood tissues whose cfRNA could appear as NM deviations, and whether existing cfRNA literature / source-paper DEGs already report them. Next step (pending approval): `_temp_` notebook in Notebook_p0 testing per-patient BH up-hits against these tissue panels (nb6 method).

## Scope

| Item | Value |
|---|---|
| Included | Pancreatic Cancer, Pancreatitis, Colorectal, Esophagus, Stomach, Lung, Liver Cancer, Liver Cirrhosis, MM, MGUS, Tuberculosis, HIV, HIV + Tuberculosis, CAD_HF+, CAD_HF-, Pre-eclampsia, ME/CFS |
| Excluded (unclear phenotype) | ICI-treated Cancer, ICI-m, Other Cancer |
| Evidence tags | [FT] full text read, [Abs] abstract only, [recall] model memory not verified, [PREPRINT], [proj] project memory |
| Marker genes | all [recall]; must be re-derived from HPA single-cell / PanglaoDB files before use |

## Inclusion criteria (revised 2026-10-27)

| Criterion | Rule |
|---|---|
| C1 (required) | In this disease, cells of the tissue die (necrosis, apoptosis, injury) or actively secrete/release RNA-bearing material |
| C2 (at least one) | cfDNA tissue-of-origin evidence for that tissue in this disease, or an established clinical blood marker of that tissue's injury/function |
| Not required | A prior cfRNA report |
| Dropped by rule | Metastatic sites as host tissue (lesion sheds tumour-derived RNA), cachexia/sarcopenia muscle and adipose (atrophy, not death) |
| Testability | Clinical subgroups exist only for Moore pancreatic (Met_Site, Tumor_Stage, Pancreatitis_Type: chronic 65 / acute 14). Chen cancers have no histology or metastasis data, so their conditional rows are not testable |

## Overview (after re-evaluation; details in Part R)

| Disease | Keep | Conditional (condition) | Dropped | Main caveat |
|---|---|---|---|---|
| Pancreatic Cancer | Tumour ductal, acinar, liver | Biliary (obstruction); islet (islet cell tumour only) | Muscle, adipose, metastatic sites, CAF stroma (C2 fails, exploratory only) | Ductal markers overlap cholangiocyte; 71 Moore benign in HC pool |
| Pancreatitis | Acinar | Liver/biliary (gallstone); islet (chronic with diabetes); lung, kidney (severe acute) | Muscle, adipose, peripancreatic fat necrosis (no C2), stroma | Aetiology and severity not in metadata |
| Colorectal | Colon epithelium, liver | - | Peritoneum, muscle, adipose, stroma | |
| Esophagus | - | Squamous (ESCC) vs columnar (EAC) by histology; liver (liver mets) | Muscle, adipose | Histology unknown: not testable |
| Stomach | Gastric epithelium | Liver (liver mets) | Peritoneum, muscle, adipose, stroma | |
| Lung | AT2 | AT1/club/basal (squamous or lung injury); lung endothelium; liver, brain (mets) | Bone, adrenal, muscle, adipose | Smoking/COPD raise lung cfDNA; HC smoking unknown |
| Liver Cancer | Hepatocyte, tumour state (AFP/GPC3) | Kidney tubule (AKI) | Lung, lymph node, adrenal, bone (mets) | Positive control |
| Liver Cirrhosis | Hepatocyte, stellate/scar mesenchyme (ELF test) | Kidney tubule (AKI); gut (decompensated/alcohol); Kupffer (activation only) | Spleen, muscle | n=4 scored; stellate mRNA release evidence in vitro only |
| MM | Blood-lineage reference only (plasma cell, erythroid) | Kidney tubule (renal impairment); osteoclast; heart (AL amyloid) | Osteoblast, marrow stroma | Clinical subgroups absent |
| MGUS | Blood-lineage reference only | - | All bone and kidney | |
| Tuberculosis | AT2, alveolar/granuloma macrophage | Liver (on-treatment DILI); lymph node, pleura (extrapulmonary TB) | Endothelium | Chang TB-negative "HC" are symptomatic coughers, lung deviation underestimated |
| HIV | Enterocyte (I-FABP, sCD14) | - | Paneth/goblet, HIVAN kidney | ART/CD4 unknown |
| HIV + TB | AT2, macrophage, enterocyte | Liver (DILI or disseminated) | Spleen | |
| CAD_HF+/- | Endothelium (low prior) | Cardiomyocyte (score against hs-TnI) | Kidney, liver, cardiac fibroblast | HF+ = future HF, groups matched on hs-TnI/NT-proBNP |
| Pre-eclampsia | Syncytiotrophoblast (EVT merged), endothelium | Liver (HELLP), brain (late gestation), cardiomyocyte | Kidney, smooth muscle | 194 pregnant samples in HC pool; most PE samples pre-symptomatic |
| ME/CFS | - | Skeletal muscle (post-exertion sampling only) | Brain, gut epithelium, stroma, endothelium | Gardella sampled at rest: effectively no testable tissue |

---

# Part R. Re-evaluation against C1/C2

## Re-evaluation of candidate tissues: pancreas / GI / lung (literature only, 2026-10-27)

Criteria: C1 (required) host cells of that tissue die or actively release RNA-bearing material at plausible scale. C2 (>=1) cfDNA tissue-of-origin evidence in this disease, or established clinical blood marker of that tissue's injury/function (tumour markers CEA/CA19-9/SCC-Ag do not count).
Principle: a metastatic site sheds tumour (primary-identity) RNA; the host tissue only qualifies if host cells are damaged (C1) and that damage is measurable (C2). Cachexia = catabolic atrophy, not mass cell death.
Tags: [FT] full text read (here or in source survey), [Abs] abstract checked via PubMed MCP, [recall] not verified.

### References (all checked in PubMed MCP unless [recall])

| Key | Citation | DOI / PMID | Tag |
|---|---|---|---|
| LW16 | Lehmann-Werman, PNAS 2016 | 10.1073/pnas.1519286113 / 26976580 | [Abs]: exocrine pancreas cfDNA in pancreatic cancer and pancreatitis; beta-cell cfDNA in new T1D and islet-graft recipients; neuron/glia cfDNA after TBI |
| Lub22 | Lubotzky, JCI Insight 2022 | 10.1172/jci.insight.153559 / 35076021 | [FT]+[Abs]: hepatocyte cfDNA up with liver mets (various primaries) and in localized pancreatic and colon cancer; neuron/oligodendrocyte/astrocyte cfDNA up in a subset with brain mets |
| Moss18 | Moss, Nat Commun 2018 | 10.1038/s41467-018-07466-6 / 30498206 | [FT]: colon strongest in 3/4 met CRC; lung strongest in 2/4 lung cancer |
| Liu20 | Liu, Ann Oncol 2020 (CCGA) | 10.1016/j.annonc.2020.02.011 / 33506766 | [Abs]: TOO panel includes colorectal, oesophagus, stomach, lung, pancreas (tumour-derived methylation) |
| Per23 | Peretz, Med 2023 | 10.1016/j.medj.2023.03.006 / 37060900 | [Abs]: lung-specific endothelial cfDNA selectively up in COPD and lung cancer |
| Nes25 | Nesselbush, Nature 2025 | 10.1038/s41586-025-08834-1 / 40240612 | [FT]: AT2-enriched cfRNA in LUAD; normal lung cfRNA up in COPD/injury/smoking |
| Ten13 | Tenner, Am J Gastroenterol 2013 (ACG guideline AP) | 10.1038/ajg.2013.218 / 23896955 | [Abs]; lipase/amylase >3x ULN diagnostic and ALT >150 suggesting gallstone aetiology are [recall] of guideline body |
| Banks13 | Banks, Gut 2013 (revised Atlanta) | 10.1136/gutjnl-2012-302779 / 23100216 | [Abs]: pancreatic and peripancreatic necrosis; organ failure defines severity. Marshall score (creatinine, PaO2/FiO2) [recall] |
| Pan09 | Pannala, Lancet Oncol 2009 | 10.1016/S1470-2045(08)70337-1 / 19111249 | [Abs]: up to 80% hyperglycaemic/diabetic; improves after resection |
| Fea11 | Fearon, Lancet Oncol 2011 | 10.1016/S1470-2045(10)70218-7 / 21296615 | [Abs]: cachexia = ongoing muscle mass loss, negative protein/energy balance (catabolism, not necrosis) |
| WL20 | Winther-Larsen, Cancers 2020 | 10.3390/cancers12102852 / 33023150 | [Abs]: serum NfL higher in lung cancer brain mets (AUC 0.77); low IF |
| Kim24 | Kim, Int J Mol Sci 2024 | 10.3390/ijms25126397 / 38928104 | [Abs]: sNfL and sGFAP up in lung cancer brain mets (n=700); low IF |
| Ish12 | Ishikawa, Respir Investig 2012 | 10.1016/j.resinv.2012.02.001 / 22554854 | [Abs]: KL-6 from injured/regenerating AT2, serum ILD marker |
| Miki11 | Miki, Proc Jpn Acad B 2011 | 10.2183/pjab.87.405 / 21785258 | [Abs]: serum pepsinogen I/II used to stage atrophic gastritis in gastric cancer risk screening |
| Cives18 | Cives & Strosberg, CA Cancer J Clin 2018 | 10.3322/caac.21493 / 30295930 | from source survey [Abs]; chromogranin A as PanNET marker [recall] |
| Tao23 | Tao, Cell Rep Med 2023 | 10.1016/j.xcrm.2023.101281 / 37992683 | [FT]: CAF/endothelial cfRNA deconvolution signal (cfRNA, not C2) |

### 1. Pancreatic cancer (PDAC / IPMN / islet cell tumour)

| Tissue | C1 mechanism | C2 cfDNA | C2 clinical marker | Verdict | Reason |
|---|---|---|---|---|---|
| Neoplastic ductal epithelium | Tumour cell turnover/death (Kleeff 2016 [Abs]) | Tumour TOO methylation, Liu20 [Abs] | None host (CA19-9 is tumour) | keep | Primary source; cfDNA TOO supports |
| Exocrine acinar | Upstream duct obstruction, obstructive pancreatitis/atrophy [recall mechanism] | Exocrine cfDNA in pancreatic cancer, LW16 [Abs] | Lipase/amylase (raised only in a minority of PDAC [recall]) | keep | C1+C2(cfDNA) met; NM already recovered CTRB2/CPB1/PNLIPRP1 [proj] |
| Liver (hepatocyte) | Biliary obstruction, liver mets, micrometastasis | Hepatocyte cfDNA up in localized PDAC and liver mets, Lub22 [FT] | ALT/AST/bilirubin [recall] | keep (strongest in head tumours / liver-met) | Best-supported host tissue; remove acute-phase genes |
| Biliary (cholangiocyte) | Obstruction, cholangitis, stenting injury [recall] | Not found | Bilirubin, ALP, GGT [recall] | conditional: biliary-obstruction/stent patients | Marker overlap with ductal tumour genes |
| Endocrine islet | PDAC: beta-cell dysfunction, death not shown [recall]. Islet tumour: tumour itself | Beta-cell cfDNA shown only in T1D/islet graft, LW16 [Abs]; not in PDAC | Glucose/HbA1c (Pan09 [Abs]); CgA for PanNET (Cives18, [recall]) | conditional: islet cell tumour patients only | For PDAC-DM, C1 (cell death) unsupported |
| Stroma (CAF/stellate) | EV secretion plausible; scale unknown | Not found | None | drop (C2 fails); exploratory only | Only cfRNA evidence (Tao23, CRC/STAD, not PDAC) |
| Skeletal muscle | Cachexia = catabolic atrophy (Fea11 [Abs]) | Not found | CK not used/raised in cachexia [recall] | drop | C1 fails |
| Adipose | Lipolysis, not death (Danai 2018 [Abs]) | Not found | None | drop | C1 fails |
| Lung | Mets shed tumour RNA | Not found for PDAC | None | drop | Host injury unsupported |
| Peritoneum | Mets shed tumour RNA | Not found | CA125 partly mesothelial [recall, unsure] | drop | C1/C2 unsupported |

### 2. Pancreatitis (acute / chronic)

| Tissue | C1 mechanism | C2 cfDNA | C2 clinical marker | Verdict | Reason |
|---|---|---|---|---|---|
| Exocrine acinar | Acinar injury/necrosis (AP), fibrotic replacement (CP) (Banks13, Boxhoorn 2020 [Abs]) | Exocrine cfDNA in pancreatitis, LW16 [Abs] | Lipase/amylase >3x ULN (Ten13/Banks13, criterion [recall]) | keep | Strongest case of all diseases here |
| Liver / biliary | Gallstone passage, cholestasis (Mederos 2021 [Abs], Beyer 2020 [Abs]) | Not found for pancreatitis | ALT >150 for gallstone AP (Ten13 [recall]); bilirubin/ALP | conditional: biliary/gallstone aetiology | Aetiology unknown in our metadata |
| Endocrine islet | CP fibrosis destroys islets; diabetes 38-40% (Singh 2019 [Abs]); mechanism of loss [recall] | Beta-cell cfDNA in islet-graft (TPIAT, includes CP) LW16 [Abs], not native CP | Glucose/HbA1c | conditional: CP with diabetes | C1 plausible but chronic, low rate |
| Stroma (stellate) | Activation/proliferation, not death; secretion unknown | Not found | None | drop (C2 fails) | Fibrosis is gain, not release |
| Peripancreatic adipose (added) | Fat necrosis in necrotizing AP (Banks13 [Abs]) | Not found | None | drop (C2 fails); revisit if necrotizing AP present | Real necrosis but no C2; outpatient cohort |
| Lung (AT1/2) | ARDS in severe AP [recall] | Not found | PaO2/FiO2 (Marshall, [recall]) | conditional: severe AP with respiratory failure | Moore pancreatitis = EUS outpatients, likely none |
| Kidney (added) | AKI in severe AP [recall] | Not found | Creatinine (Marshall, [recall]) | conditional: severe AP with AKI | Same caveat |
| Muscle / adipose (malnutrition) | Atrophy | Not found | None | drop | C1 fails |

### 3. Colorectal cancer

| Tissue | C1 mechanism | C2 cfDNA | C2 clinical marker | Verdict | Reason |
|---|---|---|---|---|---|
| Colon epithelium (tumour) | Tumour cell death | Colon cfDNA 3/4 met CRC (Moss18 [FT]); TOO (Liu20 [Abs]) | None host (CEA is tumour) | keep | Primary source |
| Liver (hepatocyte) | Liver mets; localized colon cancer also up | Hepatocyte cfDNA up in localized colon cancer, Lub22 [Abs] | ALT/AST/ALP [recall] | keep (strongest if liver-met) | Direct cfDNA evidence even without mets |
| Stroma (CAF/endothelium) | EV secretion; scale unknown | Not found | None | drop (C2 fails); exploratory only | Tao23 cfRNA only, overlapping cohort |
| Peritoneum | Mets shed tumour RNA | Not found | None (CA125 [recall, unsure]) | drop | Unsupported |
| Lung | Mets shed tumour RNA | Not found | None | drop | Unsupported |
| Muscle / adipose | Atrophy (Fea11) | Not found | None | drop | C1 fails |

### 4. Oesophageal cancer (ESCC / EAC)

| Tissue | C1 mechanism | C2 cfDNA | C2 clinical marker | Verdict | Reason |
|---|---|---|---|---|---|
| Squamous epithelium (tumour) | Tumour cell death | Oesophagus in TOO (Liu20 [Abs]) | None host (SCC-Ag tumour) | conditional: ESCC histology | Histology unknown in our data |
| Columnar / Barrett (tumour) | Tumour cell death | Oesophagus TOO (Liu20 [Abs]) | None | conditional: EAC histology | Same |
| Liver | Liver mets (Ai 2017 [Abs]) | Hepatocyte cfDNA up with liver mets of various primaries (Lub22 [Abs]); oesophageal not specifically shown [unsure] | ALT/AST/ALP [recall] | conditional: liver-met patients | Chen cohort 65-86% early stage |
| Lung | Mets shed tumour RNA | Not found | None | drop | Unsupported |
| Bone | Osteolysis by osteoclasts [recall] | Not found | ALP/bone turnover [recall] | drop | C1 scale doubtful; no cfDNA |
| Muscle / adipose | Atrophy (Fea11) | Not found | None | drop | C1 fails |

### 5. Stomach cancer

| Tissue | C1 mechanism | C2 cfDNA | C2 clinical marker | Verdict | Reason |
|---|---|---|---|---|---|
| Gastric epithelium (tumour + atrophic mucosa) | Tumour cell death; chronic atrophic gastritis background [recall] | Stomach in TOO (Liu20 [Abs]) | Pepsinogen I/II (chief-cell mass/atrophy; Miki11 [Abs]) | keep | Only GI tissue with a host function marker |
| Liver | Liver mets (48% of metastatic, Riihimaki 2016 [Abs]) | Lub22 liver-met effect (pan-primary) [Abs] | ALT/AST/ALP [recall] | conditional: liver-met patients | Early-stage cohort |
| Peritoneum | Mets shed tumour RNA | Not found | None | drop | Unsupported |
| Lung, bone | Mets shed tumour RNA | Not found | None | drop | Unsupported |
| Stroma (CAF) | EV secretion; scale unknown | Not found | None | drop (C2 fails); exploratory only | Tao23 cfRNA only |
| Muscle / adipose | Atrophy (Fea11) | Not found | None | drop | C1 fails |

### 6. Lung cancer

| Tissue | C1 mechanism | C2 cfDNA | C2 clinical marker | Verdict | Reason |
|---|---|---|---|---|---|
| AT2 (tumour + injured lung) | Tumour death (LUAD); smoking/COPD injury (Nes25 [FT]) | Lung cfDNA 2/4 (Moss18 [FT]); TOO (Liu20 [Abs]) | KL-6, SP-D from injured AT2 (Ish12 [Abs]) | keep | C1+both C2; not tumour-specific (smoking) |
| AT1 / club / basal | Basal = LUSC tumour identity; AT1 injury | Same lung cfDNA | KL-6 (AT2-biased) | conditional: LUSC (basal) or injury context | Histology/smoking unknown |
| Lung endothelium (added) | Vascular turnover in tumour lung | Lung-specific endothelial cfDNA up in lung cancer (Per23 [Abs]) | None | conditional: low specificity (also COPD) | New C2-supported candidate |
| Brain (neuron/glia) | Neural damage by brain mets | Neuron/glia cfDNA up in subset with brain mets (Lub22 [Abs]) | NfL, GFAP up with brain mets (WL20, Kim24 [Abs]) | conditional: brain-met patients only | Passes C1/C2 but BBB limits plasma RNA [recall]; cohort 86% early stage |
| Liver | Liver mets (SCLC 35%, Riihimaki 2014 [Abs]) | Lub22 liver-met effect [Abs] | ALT/AST/ALP [recall] | conditional: liver-met patients | Early-stage cohort |
| Bone | Osteolysis [recall] | Not found | ALP/bone turnover [recall] | drop | C1 scale doubtful; no cfDNA |
| Adrenal | Mets shed tumour RNA; adrenal insufficiency rare [recall] | Not found | None | drop | Unsupported |
| Muscle / adipose | Atrophy (Fea11) | Not found | None | drop | C1 fails |

### Dropped tissues

| Tissue | Diseases | One-line reason |
|---|---|---|
| Skeletal muscle | all | Cachexia is catabolic atrophy, not cell death (Fea11); no cfDNA or injury marker |
| Adipose | all (except necrotizing AP, also dropped) | Lipolysis, not death; no C2 |
| Lung (as metastatic host) | PDAC, CRC, ESCA, STAD | Metastases shed tumour RNA; no host-injury C2 |
| Peritoneum | PDAC, CRC, STAD | No host-injury C2; CA125 origin uncertain |
| Bone | ESCA, STAD, lung | No cfDNA; marrow/bone RNA release scale doubtful |
| Adrenal | lung | No C1/C2 |
| Stroma (CAF/stellate) | PDAC, CRC, STAD, pancreatitis | Fails C2; keep only as exploratory given Tao23 cfRNA |
| Peripancreatic adipose | pancreatitis | Real necrosis but no C2; outpatient cohort |

### Uncertainties

| Point | Uncertainty | Verify manually |
|---|---|---|
| Clinical cut-offs (lipase 3x ULN, ALT >150, Marshall) | Abstract-only; values [recall] | Ten13 and Banks13 full text |
| Lub22 per-primary liver-met effect for gastric/oesophageal | Abstract says "different primary tumors" | Lub22 Fig. on liver-met by primary |
| Brain RNA crossing BBB into plasma | [recall] | Lub22 brain-met fraction; cfRNA brain-met reports |
| Met/stage status, histology, AP severity in our cohorts | Not in metadata reviewed | obs metadata |

---
## Re-evaluation of candidate tissues: HCC, Cirrhosis, MM, MGUS, TB, HIV, HIV+TB (literature only, 2026-10-27)

Criteria: C1 (required) = cells of the tissue die or actively release RNA-bearing material at plausible scale in this disease. C2 (>=1) = cfDNA tissue-of-origin evidence OR established clinical blood marker of that tissue's injury/function.
Tags: [FT] full text, [Abs] abstract checked via PubMed MCP this session or in the input file, [recall] not verified. Keys from input file (RH22, CH22, SA20, SU15, GE21, LW18, LU22, RA19, GI21, IK19, NG19, TE21, Chang24, Hunt14, Brenchley06, Gounder17, Alemu20, Gupta15) are re-used; DOI/PMID are in the input section 3.0 / 4.

### 0. New references (verified this session via PubMed MCP)

| Key | Paper | DOI | PMID | Tag |
|---|---|---|---|---|
| JO15 | Johnson, J Clin Oncol 2015, ALBI grade (albumin+bilirubin = liver function in HCC) | 10.1200/JCO.2014.57.9151 | 25512453 | [Abs] |
| RO04 | Rosenberg, Gastroenterology 2004, ELF panel (HA, PIIINP, TIMP1) detects fibrosis, n=1021 | 10.1053/j.gastro.2004.08.052 | 15578508 | [Abs] |
| HF11 | Holland-Fischer, Gut 2011, sCD163 (Kupffer activation) 3x higher in cirrhosis, tracks portal pressure | 10.1136/gut.2010.234542 | 21572121 | [Abs] |
| CH14 | Charrier, Surgery 2014, activated HSC release exosomes carrying CCN2 mRNA (mouse/LX-2, in vitro) | 10.1016/j.surg.2014.04.014 | 24882759 | [Abs] |
| ZH23 | Zhang, Hepatol Commun 2023, plasma I-FABP, sCD14, LBP associated with severity in alcohol-related cirrhosis | 10.1097/HC9.0000000000000159 | 37204414 | [Abs] |
| TE11 | Terpos, Int J Cancer 2011, MM: sclerostin (osteocyte) up, correlates with CTX (+) and BALP (-); MGUS not elevated vs controls | 10.1002/ijc.27342 | 22052418 | [Abs] |
| DI10 | Dimopoulos, J Clin Oncol 2010, IMWG renal impairment consensus (creatinine/eGFR; cast nephropathy) | 10.1200/JCO.2010.30.8791 | 20956629 | [Abs] |
| SE20 | Seefried, J Bone Oncol 2020, MGUS = risk factor for osteoporosis/fracture; measured CTX/NTX, TRAP5b, P1NP (exercise trial, n=15) | 10.1016/j.jbo.2020.100323 | 33083217 | [Abs] |
| EN14 | Enomoto, Kekkaku 2014, serum SP-A, SP-D higher with cavities, fall after TB treatment (n=27, Japanese, low IF) | DOI not listed in PubMed | 25195297 | [Abs] |
| CM21 | Cheng, Med 2021, cfDNA methylation tissue-of-origin: lung and liver injury in severe COVID-19 (analogy, not TB) | 10.1016/j.medj.2021.01.001 | 33521749 | [Abs] |
| RK12 | Ramakrishnan, Nat Rev Immunol 2012, granuloma biology (macrophage necrosis drives dissemination) | 10.1038/nri3211 | 22517424 | [Abs]; necrosis detail [recall] |
| SK06 | Saukkonen, AJRCCM 2006, ATS statement on anti-TB DILI; ALT monitoring, HIV coinfection raises risk | 10.1164/rccm.200510-1666ST | 17021358 | [Abs] |
| SA11 | Sandler, J Infect Dis 2011, I-FABP, LPS, sCD14 raised in HIV; sCD14 predicts mortality | 10.1093/infdis/jiq118 | 21252259 | [Abs] |
| MO18b | Mosepele, Cardiovasc J Afr 2018, sVCAM-1/sICAM-1 raised in ART-suppressed HIV (low IF) | 10.5830/CVJA-2018-003 | 29771268 | [Abs] |

Not found by PubMed (this session): Moss-type cfDNA tissue-of-origin study in TB, HIV, MM, or cirrhosis; serum KL-6 in TB as an established marker (only EN14).

### 1. Liver Cancer (HCC)

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| Liver: hepatocytes (normal + malignant) | Tumour and peritumoral hepatocyte turnover; liver mRNA dominates HCC cfRNA (RH22, CH22 [FT]; SA20 [Abs]) | Elevated liver cfDNA fraction in HCC (SU15 [Abs]); liver main source of tumour-variant cfDNA (GE21 [Abs]) | ALT/AST, bilirubin, albumin; ALBI grade (JO15 [Abs]) | keep | C1+both C2 arms; positive control |
| Liver: oncofetal tumour hepatocyte program (AFP, GPC3) | Same death/release as above, tumour cells | Tumour cfDNA (GE21 [Abs]) | AFP (tumour marker; [recall] for guideline) | keep | Same cells as row 1, tumour state; markers [recall] |
| Liver: hepatocyte CYP/zonation genes (down) | Loss of functional mass (synthetic failure), not release | As row 1 | Albumin/bilirubin (JO15) | conditional: test as down-deviation only, within hepatocyte set | Down-direction is a function readout, not new release |
| Liver: stellate / scar mesenchyme | Activated HSC release EVs with CCN2 mRNA (CH14 [Abs], in vitro); scale in plasma unknown | None found | ELF: HA, PIIINP, TIMP1 are ECM/HSC products (RO04 [Abs]) | conditional: only if background cirrhosis is present (11/14 in BL22 [FT]) | C1 weak (in vitro only); relevant via cirrhosis, not tumour |
| Liver: LSEC / scar endothelium | LSEC capillarisation [recall]; death/release scale not shown | None found | HA (LSEC-cleared) in ELF (RO04) is an indirect function marker [recall for LSEC attribution] | conditional: as cirrhosis background, exploratory | C1 unverified |
| Liver: Kupffer / scar macrophages | Activation, not death; monocyte-derived (RA19 [Abs]) | None found | sCD163 (HF11 [Abs], cirrhosis) | conditional: separate from blood monocyte signal | Overlaps blood myeloid; C1 is activation/shedding, not death |
| Liver: cholangiocyte / ductular reaction | Ductular reaction with fibrosis [recall]; cholangiocyte injury in cholestasis [recall] | None found | ALP, GGT (cholestatic; [recall]) | conditional: mixed HCC-CCA or cholestasis | C1 [recall] only |
| Lung, lymph node, adrenal, bone (metastasis) | Metastatic deposits shed tumour (hepatocyte-identity) RNA, not host-organ RNA | None | None specific | drop | Fails C1 by user rule |
| Spleen / gut / kidney (portal hypertension) | See cirrhosis | See cirrhosis | See cirrhosis | conditional: only via background cirrhosis rows | Not HCC-specific |

### 2. Liver Cirrhosis

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| Liver: hepatocytes | Ongoing hepatocyte death drives fibrosis (GI21 [Abs]); hepatocyte cfDNA rises with acute death (LW18 [Abs]) | Hepatocyte cfDNA shown for acute injury, not tested in cirrhosis (LW18 [Abs]); unsure | ALT/AST, bilirubin, albumin (JO15 [Abs] validated in non-HCC liver disease too) | keep | Direction may be up (death) or down (synthetic failure) |
| Liver: activated stellate / scar mesenchyme | HSC EV-mRNA release (CH14 [Abs], in vitro); fibrotic niche expansion (RA19 [Abs]) | None found | ELF HA/PIIINP/TIMP1 (RO04 [Abs]) | keep | Disease-defining cell; C2 strong; C1 moderate (in vitro) |
| Liver: scar endothelium / LSEC | PLVAP+/ACKR1+ expansion (RA19 [Abs]); death/release unverified | None found | HA clearance [recall attribution] | conditional: exploratory, check HC detection rate | C1 unverified |
| Liver: Kupffer / scar macrophages | Activation (HF11 [Abs]); not death | None found | sCD163 3x up, tracks portal pressure (HF11 [Abs]) | conditional: needs non-monocyte-specific markers (CD5L, VSIG4, MARCO) | Blood myeloid overlap |
| Liver: cholangiocytes / ductular reaction | Ductular reaction scales with fibrosis [recall] | None found | ALP, GGT [recall] | conditional: cholestatic aetiology (PBC/PSC) | C1 [recall] |
| Gut epithelium (enterocyte) | Enterocyte damage + translocation (ZH23 [Abs], alcohol-related cirrhosis) | None found | I-FABP, sCD14, LBP (ZH23 [Abs]); LBP 70% up in cirrhosis (HF11 [Abs]) | conditional: decompensated / alcohol-related; low effect expected | C1+C2 met but evidence from one aetiology |
| Kidney: tubule (AKI/ATN) | ATN in decompensated cirrhosis [recall]; HRS is functional, little death | None found | Creatinine (guideline [recall]); NGAL for ATN vs HRS [recall; PubMed hits not opened] | conditional: only with documented AKI | Most cirrhosis patients lack tubular death |
| Spleen (red pulp) | Congestion/hypersplenism; cells mostly blood lineage; no death mechanism | None | None specific | drop | Fails C1 |
| Skeletal muscle (sarcopenia, missed in input) | Atrophy, not mass death | None | None | drop | Fails C1 (user rule) |

### 3. Multiple Myeloma

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| BM clonal plasma cells (blood lineage) | Disease cell, high turnover (TE21, EH25 [Abs]) | Tumour cfDNA in MM [recall] | M-protein, FLC (IMWG RJ14 [recall]) | keep | Disease cell; reference set |
| BM erythroid precursors (blood lineage) | Marrow; already dominant reported cfRNA signal (RH22 [FT]; IB20 [Abs]) | Not checked | Hemoglobin (CRAB "A", RJ14 [recall]) | keep | Positive control |
| Bone: osteoclasts | Hyperactive resorption in ~80% MM (TE21, MBD17 [Abs]); osteoclasts secrete TRAP5b; RNA release scale unknown | None found | CTX up, tracks sclerostin (TE11 [Abs]) | conditional: lytic disease present; check HC detection of CALCR/OSCAR/DCSTAMP | C1 is activity, not proven RNA release |
| Bone: osteocytes (sclerostin) (missed in input) | Osteocyte secretion of sclerostin up; osteocyte death in MM [recall] | None found | Serum sclerostin up, highest with fractures (TE11 [Abs]) | conditional: exploratory, SOST likely undetectable | Bone-encased cells; release to plasma doubtful |
| Bone: osteoblasts | Suppressed (fewer cells) (TE21 [Abs]) | None found | BALP inverse to sclerostin (TE11 [Abs]) | drop | Fails C1 (suppression, no death/release) |
| BM stroma / niche (MSC, adipocyte) | Niche support [recall]; no death/release mechanism | None found | None (DKK1 is myeloma-cell secreted) | drop | Fails C1 |
| Kidney: proximal tubule | Light-chain proximal tubulopathy, tubular injury (IK19 [Abs]; DI10 [Abs]) | None found | Creatinine/eGFR (DI10 [Abs]) | conditional: renal impairment (CRAB "R") documented | Only subset; check RH22 metadata |
| Kidney: distal tubule / TAL (cast nephropathy) | Cast nephropathy with tubular injury (DI10 [Abs]); UMOD role [recall] | None found | Creatinine/eGFR (DI10 [Abs]) | conditional: same as above | Subset only |
| Heart / other (AL amyloid, EMD) | Cardiomyocyte injury in AL [recall] | None | Troponin/NT-proBNP in AL staging [recall] | drop (unless AL amyloid known) | Rare; not verified |

### 4. MGUS

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| BM clonal plasma cells (blood lineage) | Low-burden clone (EH25 [Abs]) | Not checked | M-protein [recall] | keep | Disease cell; weak expected signal |
| BM erythroid (blood lineage) | No anaemia by definition [recall]; RH22 gradient NC to MGUS to MM [FT] | Not checked | Hemoglobin normal by definition [recall] | conditional: only as gradient comparator to MM | Little biology expected |
| Bone: osteoclasts | Fracture risk and microstructure change; DKK1/MIP-1a up (NG19 [Abs]); osteoporosis risk (SE20 [Abs]) | None found | CTX/NTX, TRAP5b measured (SE20 [Abs]); sclerostin not raised vs controls (TE11 [Abs]) | drop (or conditional: exploratory only) | C1 weak at MGUS scale; no marker shown elevated vs controls |
| Bone: osteoblasts / osteocytes | No clear change (TE11 [Abs]) | None | Sclerostin not raised (TE11) | drop | Fails C1 |
| Kidney (MGRS) | Rare; biopsy-defined (IK19 [Abs]) | None | Creatinine [recall] | drop | Rare; not expected in this cohort |

### 5. Tuberculosis

Note: Chang "TB-negative" participants are symptomatic coughers (Chang24 [FT]) and sit in our HC pool. Lung/airway/inflammation deviation of TB vs this pool will be attenuated, and the HC reference itself is contaminated.

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| Lung: alveolar AT2 | Caseous necrosis, cavitation destroys parenchyma; AT2 genes up in TB cfRNA (Chang24 [FT]) | No TB study; lung cfDNA injury shown in COVID-19 (CM21 [Abs], analogy) | SP-A/SP-D up with cavities, fall on treatment (EN14 [Abs], low IF) | keep | C1 strong; C2 marker weak-moderate |
| Lung: alveolar AT1 | Same necrosis; AT1 cover most alveolar surface [recall] | As AT2 (CM21 analogy) | KL-6 (AT2/AT1 regenerating epithelium) not changed with treatment in EN14 | conditional: if AT2 positive; AT1 markers (AGER) check HC detection | C2 weak |
| Airway epithelium (club, ciliated) | Endobronchial/cavity erosion [recall] | None | CC16 [recall] | conditional: exploratory; HC pool coughers share airway irritation | C1 [recall] only |
| Alveolar / granuloma macrophages | Infected macrophage death and necrotic core (RK12 [Abs]; necrosis detail [recall]); MARCO, C1QB up (Chang24 [FT]) | None | None specific (ADA is pleural fluid) | keep | C1 strong; C2 gap noted |
| Lymph node (germinal centre, FDC) | Caseating TB lymphadenitis in EPTB (Alemu20 [Abs]) | None | None | conditional: EPTB documented; mostly HIV+TB | C2 missing; lymphoid cells overlap blood |
| Pleura (mesothelium) | Pleural TB, 15% of EPTB in HIV (Alemu20 [Abs]); mesothelial injury [recall] | None | None in blood (ADA is pleural fluid) | conditional: pleural TB documented | Fails C2 |
| Liver (DILI) | Hepatocyte injury from isoniazid/rifampicin/pyrazinamide (SK06 [Abs]) | Liver cfDNA up in acute injury (LW18 [Abs]) | ALT (SK06 [Abs]) | conditional: samples taken during treatment | Chang sampling pre-treatment unverified |
| Endothelium | Endothelial injury in TB not verified | None | None verified for TB | drop | Fails C1 (not verified) |
| Neutrophil granule (blood reference) | NETs / neutrophil death (Cell Death Dis 2024, PMID 39085192 [Abs]) | NET-DNA, cfDNA up (same) | MPO-DNA (same) | keep (reference, not a tissue) | Positive control for already-reported |

### 6. HIV (TB-negative)

Note: Chang HIV+/TB- are also symptomatic coughers; lung signal may be non-TB pneumonia. No Chang DE for HIV vs healthy exists.

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| Small-intestinal enterocytes | Enterocyte damage, barrier loss, translocation (Brenchley06, Hunt14 [Abs]) | None found | I-FABP, LPS, sCD14 raised (SA11 [Abs]); I-FABP predicts mortality (Hunt14 [Abs]) | keep | C1+C2 strong |
| Paneth / goblet | Mucosal antimicrobial loss [recall] | None | None specific | drop | C1 [recall], no C2 |
| GALT / lymph node (CD4 T-cell depletion) | CD4 death in GALT (Mehandru04 [Abs]) | None | CD4 count (blood lineage) | conditional: treat as immune (blood) set, not tissue | Lymphoid = blood lineage |
| Endothelium | Endothelial injury/activation [recall for death] | None | sVCAM-1, sICAM-1 up even on ART (MO18b [Abs], low IF) | conditional: exploratory; C1 unverified | Activation markers, not death |
| Liver | ART hepatotoxicity, HBV/HCV coinfection [recall] | HCC cfDNA methylation in HIV cohorts (input section 4.2 [Abs]) | ALT/AST [recall] | conditional: coinfection or ART hepatotoxicity documented | Not HIV per se |
| Lung (AT2) | Only via intercurrent pneumonia in coughers | CM21 analogy | SP-D [recall] | conditional: intercurrent pneumonia | Not HIV-specific |
| Kidney (HIVAN) | Podocyte/tubule injury [recall], rare on ART | None | Creatinine [recall] | drop | Rare; unverified |

### 7. HIV+TB

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| Lung: AT2 / AT1 | As TB (Chang24 [FT]); cavitation often less in advanced HIV [recall] | CM21 analogy | SP-A/SP-D (EN14 [Abs]) | keep (AT2), conditional (AT1) | As TB |
| Alveolar / granuloma macrophages | As TB (RK12 [Abs]) | None | None | keep | As TB |
| Enterocytes | As HIV (Brenchley06, Hunt14 [Abs]) | None | I-FABP, sCD14 (SA11 [Abs]) | keep | As HIV |
| Liver (hepatic TB, DILI) | Hepatic TB with high ALP/GGT in HIV (Gounder17 [Abs]); DILI higher with HIV (SK06 [Abs]) | LW18 (acute hepatocyte death) [Abs] | ALT, ALP, GGT (SK06, Gounder17 [Abs]) | conditional: disseminated TB or on treatment | Stronger prior than TB-only |
| Lymph node (EPTB) | 56% of EPTB in HIV (Alemu20 [Abs]); caseation | None | None | conditional: EPTB documented | No C2 |
| Disseminated multi-organ (liver, spleen) | Disseminated in 87.9% of autopsy TB in HIV (Gupta15 [Abs]) | None | Liver enzymes [recall] | conditional: covered by liver row; spleen dropped | Spleen = blood lineage |
| Pleura | As TB | None | None | conditional: pleural TB documented | No C2 |
| Endothelium | As HIV | None | sVCAM-1 (MO18b) | conditional: exploratory | C1 unverified |

### 8. Dropped tissues

| Disease | Tissue | One-line reason |
|---|---|---|
| HCC | Lung, LN, adrenal, bone (metastasis) | Metastases shed tumour hepatocyte-identity RNA, not host cells |
| Cirrhosis | Spleen | No death mechanism; cells mostly blood lineage |
| Cirrhosis | Skeletal muscle (sarcopenia) | Atrophy, not cell death |
| MM | Osteoblasts | Suppressed, not dying/releasing |
| MM | BM stroma / adipocyte | No death or release mechanism cited |
| MM | Heart / AL amyloid / EMD | Rare; unverified; keep only if AL documented |
| MGUS | Osteoclasts / osteoblasts / osteocytes | No marker shown raised vs controls (TE11); C1 weak |
| MGUS | Kidney (MGRS) | Rare, biopsy-defined |
| TB | Endothelium | TB endothelial injury not verified |
| HIV | Paneth / goblet | C1 [recall] only, no C2 |
| HIV | Kidney (HIVAN) | Rare on ART, unverified |
| HIV+TB | Spleen | Blood lineage |

### 9. Caveats

| Evidence source | Evidence strength | Main uncertainty | What to verify manually |
|---|---|---|---|
| Clinical markers (JO15, RO04, HF11, SA11, TE11, DI10, SK06) | Moderate-strong (JCO, Gastroenterology, Gut, JID, AJRCCM) | Protein marker does not prove mRNA release to plasma | HC detection rate of each tissue's marker genes before testing |
| HSC EV-mRNA (CH14) | Weak (in vitro, mouse/LX-2) | Plasma scale unknown | Check stellate genes' detection in HC |
| TB lung markers (EN14) | Weak (n=27, low IF, Japanese) | No larger SP-D/KL-6 TB study verified | Search Kekkaku/Int J Tuberc Lung Dis for larger cohorts |
| cfDNA tissue-of-origin for TB/HIV/MM/cirrhosis | Absent | May exist outside PubMed | bioRxiv/medRxiv search |
| Kidney/cholangiocyte/DILI rows | Conditional on clinical metadata | Metadata (AKI, eGFR, treatment timing) may be missing | Check RH22 and Chang metadata |

---
## Re-evaluation: CAD_HF+/HF-, Pre-eclampsia, ME/CFS (C1 required, C2 >=1)

Source section: `_temp_disease_tissue_lit.md` "Tissues of interest: CAD (HF+/HF-), pre-eclampsia, ME/CFS". All PMIDs below were checked in PubMed (2026-10-27). Tags: [FT] full text, [Abs] abstract only, [recall] not verified this pass.

Sampling-time context:
- Ward 2022: CAD patients sampled about 4 months after ACS while stable. HF+ was sampled **before** HF onset and matched to HF- on NT-proBNP, hs-TnI, LVEF and creatinine ([FT], 10.3390/cells11203191 / 36291058).
- Moufarrej 2022: PE samples were mostly taken before clinical onset (the panel uses 5-16 weeks) ([Abs], 10.1038/s41586-022-04410-z / 35140405). End-organ injury (HELLP, AKI, eclampsia) is a feature of clinical disease, so these samples are expected to show little of it. Check sample timing in our metadata.
- Gardella 2025: ME/CFS samples taken at rest. Tabula Sapiens deconvolution found minimal solid-tissue cfRNA ([FT], 10.1073/pnas.2507345122 / 40789036).

### CAD_HF+ / CAD_HF-

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| Cardiomyocyte | Chronic low-grade myocyte injury: hs-cTnT detectable in 97.7% of stable CAD and predicts HF (Omland, NEJM 2009, 10.1056/NEJMoa0805299 / 19940289 [Abs]) | Cardiomyocyte cfDNA is about zero in healthy plasma and high only in acute STEMI or sepsis (Zemmour, Nat Commun 2018, 10.1038/s41467-018-03961-y / 29691397 [Abs]). No stable-CAD data | hs-cTnT/hs-TnI and NT-proBNP (Omland 2009). Ward used both for matching [FT] | conditional: only as a per-sample continuous score checked against hs-TnI; do not expect a group effect, and drop for HF+ vs HF- | C1 and C2 are met, but release is ng/L-level. Ward matched HF+/HF- on hs-TnI and NT-proBNP, so no contrast between them is expected |
| Vascular endothelium | Endothelial apoptosis: CD31+/annexin V+ apoptotic microparticles rise in CAD and track coronary endothelial dysfunction (Werner, ATVB 2005, 10.1161/01.ATV.0000191634.13057.15 / 16239600 [Abs]). EMPs remain raised 6 months after ACS (Koganti, Thromb Res 2021, 10.1016/j.thromres.2021.07.018 / 34418680 [Abs]; small n, low IF) | Endothelium is about 10% of healthy cfDNA (Moss, Nat Commun 2018, 10.1038/s41467-018-07466-6 / 30498206 [Abs]). Endothelial cfDNA falls after LVAD in advanced HF (Park, Circ Heart Fail 2025, 10.1161/CIRCHEARTFAILURE.124.013302 / 41307146 [Abs]). No stable-CAD cfDNA data | No routine clinical endothelial-injury marker. CEC/EMP counts are research-only | keep (low prior) | C1 is met, and C2 cfDNA exists but only for HF. The endothelium has a real baseline share in plasma, so a shift is detectable in principle. Ward's post-ACS timing fits the persistent-EMP data |
| Kidney (tubule/podocyte) | Cardiorenal syndrome needs HF or a low-output state (Rangaswami, Circulation 2019, 10.1161/CIR.0000000000000664 / 30852913 [Abs]). No evidence that renal cells release RNA into blood in stable CAD | None found for plasma kidney cfDNA in CAD/HF | Creatinine/eGFR exists, but Ward matched HF+/HF- on creatinine [FT] | drop | C1 not met. Kidney injury products go mostly to urine, and HF had not yet occurred at sampling |
| Liver (hepatocyte) | Congestive/ischaemic hepatopathy happens in established chronic HF (Sessa, JHEP Rep 2021, 10.1016/j.jhepr.2021.100249 / 33665589 [Abs]) | Liver cfDNA falls after LVAD, in advanced HF only (Park 2025 [Abs]) | ALT/AST/bilirubin in HF (Sessa 2021) | drop (conditional only if a cohort has established, congested HF) | C1 depends on HF that was absent at sampling. Hepatocytes are already a large baseline cfRNA source, so a small shift would be swamped |
| Cardiac fibroblast | No death or RNA-release mechanism. Fibroblasts activate and deposit ECM. FGF23 is not fibroblast-specific [recall] | None | None | drop | C1 not met |
| Platelet / megakaryocyte (missed; blood) | Platelet activation, and PMPs persist after ACS despite DAPT (Koganti 2021 [Abs]) | n/a (anucleate) | Platelet count; antiplatelet therapy is universal in CAD [recall] | conditional: handle as a blood/treatment confound, not as a tissue of interest | Drug exposure (aspirin/P2Y12 inhibitors) differs from controls. A platelet shift would not be disease-specific |

### Pre-eclampsia

Shared caveat: the HC pool contains 194 normotensive pregnant Moufarrej samples. A trophoblast Z is therefore a PE-vs-normal-pregnancy contrast only if gestational age, or at least pregnancy, enters as a covariate. Otherwise it mixes the effect of pregnancy with the effect of PE.

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| Placenta: syncytiotrophoblast | STB stress, with increased shedding of STB-derived EVs and soluble factors (Burton, BMJ 2019, 10.1136/bmj.l2381 / 31307997 [Abs]; Erlandsson, J Intern Med 2021, 10.1111/joim.13349 / 34146434 [Abs]). Placental CRH mRNA in maternal plasma is about 10x higher in PE (Ng, Clin Chem 2003, 10.1373/49.5.727 / 12709362 [Abs]) | Fetal (placental) cfDNA is 5x higher in PE (Lo, Clin Chem 1999, PMID 9931039, no DOI in PubMed [Abs]), partly from slower clearance (Lau, Clin Chem 2002, PMID 12446469 [Abs]). In EOPE after diagnosis, the relative placental fraction is lower while total cfDNA is 20x higher (van Vliet, Placenta 2025, 10.1016/j.placenta.2025.12.006 / 41477915 [Abs]) | sFlt-1:PlGF ratio (Zeisler, NEJM 2016, 10.1056/NEJMoa1414838 / 26735990 [Abs]); sFlt1 is of placental origin (Maynard, JCI 2003, 10.1172/JCI17189 / 12618519 [Abs]) | keep, conditional on a gestational-age/pregnancy covariate | Strongest C1+C2 in this set. Direction is unclear: Moufarrej saw lower placental scores early, while CRH mRNA is higher at diagnosis. Expect direction to depend on timing |
| Placenta: EVT | Defective invasion in early-onset PE is a failure to invade, not cell death. EVT has little contact with maternal blood compared with STB (Burton 2019 [Abs]; Tsang, PNAS 2017, 10.1073/pnas.1710470114 / 28830992 [Abs], which reports an EVT cfRNA signal) | None tissue-specific | None EVT-specific. PAPP-A/PAPPA2 are mainly made by STB/decidua [recall, unsure] | conditional: merge into one trophoblast signature, and report EVT as a separate score only for early-onset PE | C1 is weak on its own. Tsang 2017 gives cfRNA support, but EVT cannot be separated from STB by C2 |
| Systemic endothelium | Endothelial damage: circulating endothelial cells 88/mL in PE vs 16/mL in normal pregnancy, falling after delivery (Grundmann, AJOG 2008, 10.1016/j.ajog.2007.09.049 / 18068139 [Abs]). sFlt1 causes endothelial dysfunction (Maynard 2003) | Endothelium-specific cfDNA rises in EOPE after diagnosis (van Vliet 2025 [Abs]) | sFlt-1:PlGF (anti-angiogenic state; Zeisler 2016). CEC count is research-only | keep | C1 shows direct cell detachment/death. C2 cfDNA is direct. Moufarrej also found an endothelial enrichment [FT, input] |
| Kidney: glomerular endothelium / podocyte | Podocytes are shed into **urine** (podocyturia), and this precedes proteinuria (Craici, Hypertension 2013, 10.1161/HYPERTENSIONAHA.113.01115 / 23529165 [Abs]). Glomerular endotheliosis (Maynard 2003). No route from kidney into blood is shown | Kidney fraction is relatively **lower** in EOPE after diagnosis (van Vliet 2025 [Abs]) | Proteinuria; creatinine in severe PE [recall] | drop (glomerular endothelial injury is covered by "endothelium") | C1 release goes to urine, not plasma. The cfDNA direction is down, and Moufarrej's kidney score also fell, which fits expression change rather than release |
| Liver: hepatocyte | Hepatocyte injury in HELLP (AST >=70, LDH >600), in 10-20% of severe PE (Haram, BMC Pregnancy Childbirth 2009, 10.1186/1471-2393-9-8 / 19245695 [Abs]; review, moderate IF) | Liver fraction relatively lower in EOPE (van Vliet 2025 [Abs]) | ALT/AST/LDH/platelets (Haram 2009) | conditional: only in samples with HELLP or raised transaminases at draw | C1 holds only in a minority with clinical disease. Pre-symptomatic Moufarrej samples are unlikely to qualify. The reported down-shift is not release |
| Brain: neurons/oligodendrocytes | Neuroaxonal injury: plasma NfL higher at weeks 33-37 and tau at week 37 in women who later develop PE (Bergman, PLoS One 2018, 10.1371/journal.pone.0196025 / 29719006 [Abs]; n=16 cases) | None | NfL, tau, S100B, NSE (Bergman 2018) | conditional: late-gestation or eclampsia samples only; treat a neuronal hit as unconfirmed until it is checked against non-CNS expression of the same genes | C1 and C2 are met only in late gestation and with small n. CNS RNA must cross the BBB. The Moufarrej "neuronal" genes at 5-16 weeks are unlikely to be brain-derived [inference] |
| Smooth muscle | No death or release mechanism. Vascular tone change is functional | None | None | drop | C1 not met. Moufarrej direction conflicts [FT, input] |
| Cardiomyocyte (missed) | Subclinical myocardial stress/injury in PE [recall on mechanism] | None | cTnI and NT-proBNP differ between PE and normotensive pregnancy (meta-analysis: Alma, ESC Heart Fail 2017, 10.1002/ehf2.12129 / 28451444 [Abs]) | conditional: severe/clinical PE only; low prior | C2 is met, but C1 is ng/L-level release as in CAD. Moufarrej's "dilated cardiomyopathy" term is weak support |
| Platelet (missed; blood) | Platelet consumption in HELLP (Haram 2009). Moufarrej saw a platelet signal before 12 weeks [FT, input] | n/a | Platelet count (Haram 2009) | conditional: handle as a blood compartment, not a tissue | It is a real PE-related axis, but it belongs to blood deconvolution |

### ME/CFS

| Tissue | C1 mechanism (cite) | C2 cfDNA evidence (cite) | C2 clinical blood marker (cite) | Verdict | Reason |
|---|---|---|---|---|---|
| Skeletal muscle | Necrosis plus regeneration in biopsies taken **1 day after exercise** in post-COVID ME/CFS; mitochondrial damage (Scheibenbogen, J Cachexia Sarcopenia Muscle 2025, 10.1002/jcsm.13669 / 39727052 [Abs], review). WASF3/ER stress, which is dysfunction, not death (Wang, PNAS 2023, 10.1073/pnas.2302738120 / 37579159 [Abs]) | None found | CK is **lower**, not higher: 54 U/L in severe ME/CFS vs 101.5 in HC (Nacul, Diagnostics 2019, 10.3390/diagnostics9020041 / 30974900 [Abs]; low IF, large n=408) | conditional: only for post-exertion samples (not available in Gardella's design) | At rest, the evidence argues against muscle release: CK is low and Gardella's Tabula Sapiens deconvolution found no difference. Necrosis is shown only after exertion |
| Brain/CNS | No neuronal death shown. Walitt 2024 is functional/neurochemical [input, Abs] | None | Plasma NfL higher in ME/CFS (Biomedicines 2024, 10.3390/biomedicines12071539 / 39062112 [Abs]; single study, low IF) | drop | C1 not met. C2 is one low-IF study. A BBB-crossing RNA signal at this scale is implausible |
| Gut epithelium | Microbial translocation: LPS, LBP and sCD14 higher in ME/CFS (Giloteaux, Microbiome 2016, 10.1186/s40168-016-0171-4 / 27338587 [Abs]). I-FABP was measured, but the abstract does not report it as raised [unsure] | None | LPS/LBP/sCD14 are translocation markers, not markers of epithelial death | drop | C1 not met: barrier leak is not shown to involve enterocyte death or release into blood. The marker of enterocyte injury (I-FABP) was not reported raised |
| Fibroblast/stroma | No death mechanism. Gardella genes (COL1A1, MMP2, ADAM12, TGFB2) are **lower** in cases [FT, input] | None | None | drop as a tissue; keep as a gene-level observation | C1 not met. The direction is opposite to release |
| Endothelium (missed) | Peripheral endothelial **dysfunction** by EndoPAT (Scherbakov, ESC Heart Fail 2020, 10.1002/ehf2.12633 / 32154656 [Abs]). This is dysfunction, not shown to involve death | None | None established | drop | C1 not met |

### Dropped tissues

| Disease | Tissue | One-line reason |
|---|---|---|
| CAD_HF± | Kidney | No release mechanism without HF; output goes to urine; creatinine matched |
| CAD_HF± | Liver | Congestion needs established HF, absent at sampling; large hepatocyte baseline |
| CAD_HF± | Cardiac fibroblast | Activation and ECM deposition, no death or RNA release |
| PE | Kidney | Podocytes shed into urine; kidney cfDNA and cfRNA scores go down, not up |
| PE | Smooth muscle | Functional tone change only; direction conflicting |
| ME/CFS | Brain/CNS | No neuronal death; single low-IF NfL study |
| ME/CFS | Gut epithelium | Translocation markers, not epithelial injury; I-FABP not reported raised |
| ME/CFS | Fibroblast/stroma | No release; genes lower in cases |
| ME/CFS | Endothelium | Dysfunction only |

### Key uncertainties

| Evidence Source | Evidence Strength | Main Uncertainty | What to Verify Manually |
|---|---|---|---|
| Omland 2009, Zemmour 2018 | High (NEJM, Nat Commun) | Whether ng/L-level troponin release yields detectable cardiomyocyte cfRNA | Ward per-sample hs-TnI vs a cardiomyocyte score |
| Grundmann 2008, van Vliet 2025 | Moderate (n=10; n=21) | Small samples; cfDNA reflects death, cfRNA can also reflect secretion | Endothelial direction in Moufarrej per-sample scores |
| Ng 2003, Lo 1999 | Moderate (Clin Chem) | Measured at diagnosis; Moufarrej samples are pre-symptomatic | Gestational age at draw for each of our PE samples |
| Nacul 2019, Scheibenbogen 2025 | Low-moderate | Muscle necrosis shown only after exertion; review evidence | Whether any ME/CFS sample was drawn after exertion |
| Moufarrej sampling timing | [recall]/[Abs] | Share of PE samples drawn after clinical onset | Metadata: collection week vs diagnosis week |

---

# Part S. Original survey (pre-criteria, kept for reference)

## 1. Pancreas: tissues of interest (PDAC / IPMN / islet cell tumour / pancreatitis)

Survey date: 2026-10-27. All papers checked with PubMed MCP (title, journal, DOI, PMID).
Tags: [FT] = PMC full text read. [Abs] = abstract only. [recall] = from memory, not verified. [proj] = this project's memory notes, not literature.
Moore 2025 and the cfRNA paradigm papers are already covered in `lit/A_cfrna_paradigms.md` §3. Here I only add tissue-specific details.

### A. Pancreatic cancer (PDAC, IPMN, islet cell tumour)

#### A1. Source papers of our data

| Paper | DOI / PMID | Method | Blood/immune signals | Non-blood tissue signals | Supplementary data |
|---|---|---|---|---|---|
| Moore et al., Nat Commun 2025 | 10.1038/s41467-025-62685-y / 40783559 | cf-normalisation (NMF, platelet factor) [FT]. Group DE: Wilcoxon, PDAC vs all other diagnoses, FDR ≤ 0.05, 29 genes [FT]. HPA tissue z-score enrichment [FT]. nu-SVR deconvolution over 256 HPA tissues [FT]. RF classifier [FT] | Genes tied to immune tissues (lymph node, bone marrow, tonsil, spleen) all down in PDAC [FT]. Gene names were stripped from the PMC text. From [proj]: RPS4X, RPS13, RPS14, CD37, MS4A1, TNFRSF13C, PAX5, FAM117A (a B-cell axis) | Liver: 4 liver-unique DE genes up [FT]. From [proj]: SERPINA1, ACOT12, CRP, AFMID. Liver deconvolution share ">3 times" higher in PDAC, including early-stage patients and those without liver metastasis [FT]. "Pancreas" is listed among the tissues where the DE genes are highly expressed, but which genes drive it is not stated [FT]. No acinar, islet, muscle or adipose result | Supplementary Dataset 1 only (pdac_meta.csv per [proj]). No full DE table was released. Code: github ohsu-cedar-comp-hub/2024_ngolab_cfRNA [proj]. The project has already reproduced the DE (config.MOORE2025_DE_DIR) [proj] |
| Reggiardo et al., Nat Biomed Eng 2023 | 10.1038/s41551-023-01081-7 / 37652985 | COMPLETE-seq: repeat-aware Salmon quantification. DESeq2 (~age + gender + input_volume + condition), padj < 0.01. Elastic-net LR classifier [FT]. ExoRNeasy (EV-enriched) cfRNA [FT] | None reported | No tissue-of-origin analysis. Signal is Alu subfamilies (AluY, AluSc, AluSg7, AluSc8, AluSx3, AluSg) up in pancreatic cancer [FT] | Supplementary Information plus cohort metadata (32 members) [FT]. A gene-level DE table was not found in the text |

#### A2. Other plasma cfRNA and cfDNA-methylation studies reporting tissue signals in PDAC

| Paper | DOI / PMID | Analyte | Tissue signal | Tag |
|---|---|---|---|---|
| Lubotzky et al., JCI Insight 2022 | 10.1172/jci.insight.153559 / 35076021 | cfDNA methylation (hepatocyte markers) | Pancreatic cancer had "the highest concentrations of hepatocyte cfDNA" among primaries, in both local and metastatic disease. Local head/neck tumours trended higher than tail tumours (p = 0.06). Bilirubin and ALP did not differ between head and tail, "arguing against prominent bile duct damage". Cause unresolved: obstruction vs micrometastasis vs systemic factors | [FT] |
| Lehmann-Werman et al., PNAS 2016 | 10.1073/pnas.1519286113 / 26976580 | cfDNA methylation | Exocrine pancreas DNA detected in pancreatic cancer and pancreatitis | [Abs] (PMC text was empty) |
| Yu et al., Gut 2020 (epub 2019) | 10.1136/gutjnl-2019-318860 / 31562239 | Plasma EV long RNA (exLR-seq), SVM | 8-gene signature: FGA, KRT19, HIST1H2BK, ITIH2, MARCH2, CLDN1, MAL2, TIMP1. FGA and ITIH2 are hepatocyte-expressed and KRT19, CLDN1 and MAL2 are epithelial/ductal [recall for the tissue labels]. The abstract does not state whether the authors did a tissue attribution | [Abs] |
| Li et al., Mol Ther Nucleic Acids 2021 | 10.1016/j.omtn.2021.08.017 / 34631279 | EV long RNA, subtyping plus deconvolution | Immune-cell deconvolution and immunogenic subtypes (n = 345 PDAC) | [Abs] |
| Moss et al., Nat Commun 2018 | 10.1038/s41467-018-07466-6 / 30498206 | cfDNA methylation atlas | Healthy baseline: WBC 55%, erythroid progenitors 30%, endothelium 10%, hepatocytes 1%. No PDAC cohort in the abstract | [Abs] |
| Loyfer et al., Nature 2023 | 10.1038/s41586-022-05580-6 / 36599988 | WGBS atlas of 39 sorted cell types | Reference only. Supplies cell-type methylation markers for liquid biopsy. Whether acinar, duct and islet cells are among the 39 types is [recall], so check | [Abs] |

#### A3. Pathophysiology-based candidate tissues (ranked)

| Rank | Tissue | Mechanism | Disease subtype | Evidence |
|---|---|---|---|---|
| 1 | Pancreas: neoplastic ductal epithelium and stroma (CAF/stellate) | Tumour and desmoplastic stroma shed RNA directly | PDAC, IPMN | Kleeff, Nat Rev Dis Primers 2016, 10.1038/nrdp.2016.22 / 27158978 [Abs]. Stroma "central" (Roth, Nat Rev Dis Primers 2026, 10.1038/s41572-026-00699-6 / 42098138) [Abs] |
| 2 | Liver: hepatocytes, cholangiocytes | Biliary obstruction from head tumours, liver metastasis (the most common site), portal-circulation factors | PDAC | Lubotzky 2022 [FT]. Moore 2025 [FT]. Liver vs lung vs distant-node metastasis in 13,233 stage IV cases: Oweira, World J Gastroenterol 2017, 10.3748/wjg.v23.i10.1872 / 28348494 [Abs]. Jaundice as a presenting feature [recall], from Kleeff 2016 MeSH "Jaundice" only |
| 3 | Pancreas: acinar (exocrine) | Duct obstruction, parenchymal atrophy and exocrine insufficiency. Tumour-adjacent acinar damage | PDAC, IPMN (via AP) | Lehmann-Werman 2016 [Abs]. Danai, Nature 2018, 10.1038/s41586-018-0235-7 / 29925948: decreased exocrine function drives adipose loss [Abs] |
| 4 | Pancreas: islet endocrine | Up to 80% of patients are hyperglycaemic or diabetic. Diabetes improves after resection | PDAC. Islet cell tumour (as the tumour itself) | Pannala, Lancet Oncol 2009, 10.1016/S1470-2045(08)70337-1 / 19111249 [Abs] |
| 5 | Skeletal muscle, adipose | Cachexia: adipose and muscle loss early and common at diagnosis | PDAC | Danai 2018 [Abs] |
| 6 | Lung, distant lymph node, peritoneum | Metastatic sites | PDAC | Oweira 2017 [Abs] (lung and distant node). Peritoneum [recall] |
| 7 | Neuroendocrine tumour cells, plus liver | PanNET hormone and granin secretion. Liver-directed therapy is common, implying hepatic metastasis | Islet cell tumour | Cives & Strosberg, CA Cancer J Clin 2018, 10.3322/caac.21493 / 30295930 [Abs]. Liver as the dominant metastatic site [recall] |

#### A4. Final tissue-of-interest list (PDAC / IPMN / islet cell tumour)

Cell-type names follow HPA single cell type (HPA) and PanglaoDB (PDB). The exact term strings are [recall], so check them against the database before use. Marker genes are canonical [recall].

| Tissue | Cell types (HPA / PDB) | Markers | Rationale | Reported in cfRNA literature? |
|---|---|---|---|---|
| Liver | HPA: Hepatocytes. PDB: Hepatocytes | ALB, APOA1, APOH, AHSG, FGB, TF, HP, SERPINC1, CYP2E1, F2 | Obstruction, metastasis or systemic effect (A3 rank 2) | Yes: Moore 2025 [FT]. FGA and ITIH2 in Yu 2020 [Abs]. cfDNA: Lubotzky 2022 [FT] |
| Biliary | HPA: Cholangiocytes. PDB: Cholangiocytes | KRT7, KRT19, SOX9, CFTR, SPP1, ANXA4, EPCAM, SCTR | Duct obstruction and cholestasis. Overlaps with ductal tumour markers, so interpret with care | No dedicated report found. KRT19 appears in Yu 2020 [Abs], not attributed to tissue |
| Exocrine pancreas | HPA: Exocrine glandular cells. PDB: Acinar cells | PRSS1, CPA1, CPB1, CELA3A, CTRB2, CTRC, PNLIP, PNLIPRP1, CEL, SPINK1 | Parenchymal damage or atrophy (rank 3) | Not found in cfRNA. cfDNA: yes, Lehmann-Werman 2016 [Abs]. [proj]: CTRB2, CPB1 and PNLIPRP1 were recovered per-patient by NM |
| Pancreatic duct (normal and neoplastic) | HPA: Ductal cells. PDB: Ductal cells | KRT19, MUC1, TFF1, AGR2, S100P, CEACAM6, MSLN, CLDN18, LGALS4, KRT17 | Tumour-derived (rank 1). IPMN mucin phenotype: add MUC5AC, MUC2 | Partly: Yu 2020 KRT19, CLDN1, MAL2 [Abs]. Moore 2025: none [FT/proj] |
| Pancreatic stroma | HPA: Fibroblasts. PDB: Pancreatic stellate cells | COL1A1, COL1A2, COL3A1, FAP, POSTN, LUM, DCN, SPARC, THBS2, MMP11 | Desmoplasia (rank 1) | Partly: TIMP1 in Yu 2020 [Abs]. No stromal tissue attribution found |
| Endocrine pancreas | HPA: Pancreatic endocrine cells. PDB: Beta, Alpha, Delta, Gamma (PP) cells | INS, IAPP, GCG, SST, PPY, CHGA, CHGB, PCSK1, MAFA, SCG2 | Diabetes in PDAC (rank 4). Islet cell tumour is the tumour itself (rank 7) | Not found in cfRNA. [proj]: MAFA, PPY and PYY were recovered per-patient |
| Skeletal muscle | HPA: Skeletal myocytes. PDB: Skeletal muscle cells | ACTA1, CKM, MB, MYH1, MYH2, MYH7, MYL1, TNNT1, TNNT3, DES | Cachexia (rank 5) | Not found |
| Adipose | HPA: Adipocytes. PDB: Adipocytes | ADIPOQ, LEP, PLIN1, PLIN4, CIDEA, CIDEC, LPL, FABP4, GPD1 | Cachexia (rank 5) | Not found |
| Lung | HPA: Alveolar cells type 2/1. PDB: Pulmonary alveolar type II cells | SFTPC, SFTPB, SFTPA1, SFTPA2, NAPSA, LAMP3, AGER | Metastasis (rank 6). Low prior | Not found |

### B. Pancreatitis (acute and chronic)

#### B1. Source paper of our data

| Paper | DOI / PMID | Role of pancreatitis | Reported signals |
|---|---|---|---|
| Moore et al., Nat Commun 2025 | 10.1038/s41467-025-62685-y / 40783559 | Control class only: CEDAR has 59 acute or chronic pancreatitis cases, BCC has 23 [FT] | No pancreatitis-vs-healthy or pancreatitis-specific DE, and no tissue analysis [FT]. Reggiardo 2023 has no pancreatitis samples [FT] |

#### B2. Other plasma studies with tissue signals

| Paper | DOI / PMID | Analyte | Tissue signal | Tag |
|---|---|---|---|---|
| Lehmann-Werman 2016 | 10.1073/pnas.1519286113 / 26976580 | cfDNA methylation | Exocrine pancreas DNA in pancreatitis | [Abs] |
| Yu 2020 | 10.1136/gutjnl-2019-318860 / 31562239 | EV long RNA | Chronic pancreatitis used only as a comparator (n = 100). No chronic-pancreatitis tissue signal in the abstract | [Abs] |
| Al-Mashat, Int J Mol Sci 2025 | 10.3390/ijms262110792 / 41226826 | cfDNA promoter hypermethylation | Hypermethylation in acute pancreatitis over time vs healthy. Not tissue-of-origin | [Abs] (low IF) |
| Vasu, Clin Transl Med 2023 | 10.1002/ctm2.1434 / 37846205 | Plasma miRNA (total pancreatectomy with islet autotransplant) | miRNAs of pancreatic damage in chronic pancreatitis and at islet infusion | [Abs] |

No plasma cfRNA tissue-of-origin study of pancreatitis was found in PubMed. Coverage of this search is uncertain.

#### B3. Pathophysiology-based candidate tissues (ranked)

| Rank | Tissue | Mechanism | Evidence |
|---|---|---|---|
| 1 | Pancreas: acinar | Acinar injury and necrosis in AP. Fibrotic replacement in CP | Boxhoorn, Lancet 2020, 10.1016/S0140-6736(20)31310-6 / 32891214 [Abs]. Beyer, Lancet 2020, 10.1016/S0140-6736(20)31318-0 / 32798493 [Abs] |
| 2 | Pancreatic stroma (stellate cells, fibrosis) | CP is "fibroinflammatory" with "extensive fibrotic tissue replacement" | Beyer 2020 [Abs]. Thierens, Lancet 2024, 10.1016/S0140-6736(24)02187-1 / 39647500 [Abs] |
| 3 | Pancreas: islet endocrine | Diabetes in 38-40% of CP patients | Singh, JAMA 2019, 10.1001/jama.2019.19411 / 31860051 [Abs] |
| 4 | Liver, biliary | Gallstone aetiology of AP. Cholestasis needing intervention in CP | Mederos, JAMA 2021, 10.1001/jama.2020.20317 / 33496779 [Abs] (cholecystectomy). Beyer 2020 [Abs] ("cholestasis") |
| 5 | Muscle, adipose (malnutrition) | Exocrine insufficiency in 30-48% of CP patients, and malnutrition | Singh 2019 [Abs]. Beyer 2020 [Abs] |
| 6 | Lung, kidney (distant organ failure) | Prognosis "mainly depends on the development of organ failure" | Boxhoorn 2020 [Abs]. That lung and kidney are the main failing organs is [recall] |

#### B4. Final tissue-of-interest list (pancreatitis)

| Tissue | Cell types (HPA / PDB) | Markers | Rationale | Reported in cfRNA literature? |
|---|---|---|---|---|
| Exocrine pancreas | Exocrine glandular cells / Acinar cells | PRSS1, CPA1, CPB1, CELA3A, CTRB2, CTRC, PNLIP, CEL, SPINK1, REG1A | Primary injury | cfRNA: not found. cfDNA: yes, Lehmann-Werman 2016 [Abs] |
| Pancreatic stroma | Fibroblasts / Pancreatic stellate cells | COL1A1, COL3A1, FAP, POSTN, LUM, DCN, SPARC, PDGFRB, ACTA2 | Chronic pancreatitis fibrosis | Not found |
| Endocrine pancreas | Pancreatic endocrine cells / Beta, Alpha cells | INS, IAPP, GCG, SST, PPY, CHGA, PCSK1, MAFA | Diabetes in chronic pancreatitis | Not found (plasma miRNA only, Vasu 2023 [Abs]) |
| Liver, biliary | Hepatocytes, Cholangiocytes | ALB, APOA1, FGB, SERPINC1, CYP2E1; KRT7, KRT19, SOX9, CFTR | Gallstone acute pancreatitis, cholestasis | Not found for pancreatitis |
| Lung | Alveolar cells type 2 / Pulmonary alveolar type II cells | SFTPC, SFTPB, SFTPA1, NAPSA, AGER | Organ failure in severe acute pancreatitis. Only severe cases; our cohort is EUS outpatients, so a low prior | Not found |

### C. Caveats

| Issue | Note |
|---|---|
| Ductal and cholangiocyte markers overlap (KRT19, SOX9, CFTR, SPP1) | Tumour-derived and biliary origin cannot be separated by these genes. Report them as a joint "ductal/biliary" axis |
| Hepatocyte genes respond to acute-phase signalling (CRP, SAA1, SERPINA1, HP) | Liver cfRNA up may reflect hepatocyte transcription rather than tissue damage. Remove acute-phase genes from the "liver damage" marker set |
| Moore cohort controls | Moore's control class contains pancreatitis, IPMN and benign EUS patients (all pancreatic disease). Our HC pool includes 71 Moore benign patients [proj], so pancreatic tissue signal in HC may be inflated |
| Islet cell tumour n | 3 in CEDAR and 6 in BCC [FT]. Only case-level description is feasible |
| Unverified | HPA/PDB term strings, the marker lists, PanNET liver metastasis frequency and the peritoneum as a PDAC site are all [recall] |

---

## 2. Tissues of interest: CRC, esophageal, stomach and lung cancer (survey, 2026-10-27)

**Tags:** [FT] = PMC full text read. [Abs] = abstract only. [recall] = from memory, not verified. [PREPRINT] = not peer-reviewed.
**Sources:** PubMed MCP metadata and full text, bioRxiv MCP, and the bioRxiv JATS XML (read with WebFetch, so it is a model summary).
**Overlap:** A_cfrna_paradigms.md already covers Larson 2021, Nesselbush 2025, Wang 2024, Tao 2023 and Chen 2022. Here I only add the tissue-level details.

### 0. Shared source paper (all four diseases)

| Item | Content | Tag |
|---|---|---|
| Paper | Chen S et al. eLife 2022;11:e75181. DOI 10.7554/eLife.75181, PMID 35816095, PMC9273212 | [FT] |
| Cohort | 295 plasma samples sequenced (54 CRC, 37 STAD, 35 lung, 31 ESCA newly sequenced, plus liver cancer and HD); 263 passed QC. Pre-treatment samples from 6 Chinese centres, collected in EDTA tubes and spun within 2 h. Norgen extraction, SMARTer Stranded Total RNA-Seq Pico (DASH rRNA depletion) | [FT] |
| Stage | Early stage (I-II) ranged from 65% (stomach) to 86% (lung) | [FT] |
| DE method | edgeR quasi-likelihood, abs(log2FC) > 1, FDR < 0.05, cancer vs HD and one-vs-rest. KEGG enrichment with clusterProfiler. RUVg (25% least-variable genes as empirical controls) used only for the ML input | [FT] |
| Histology | Not stated in the text: lung (LUAD vs LUSC) and esophageal (ESCC vs EAC) subtypes are not given, and the abstract lists "LUAD, ESCA" (A_cfrna). Check the Supplementary file or GEO metadata | [FT] |
| Suppl. / data | Figure 5 (recurrent one-vs-rest features) and Figure 3 supplement 1 (per-cancer KEGG) hold the per-cancer results. I did not open the supplementary file names or the DEG tables. Data are probably GEO GSE174302 [recall, unsure] | [FT] / [recall] |

### 1. Colorectal cancer (CRC)

#### 1.1 Source paper signals (Chen 2022)

| Class | Reported | Tag |
|---|---|---|
| Blood/immune (pan-cancer vs HD) | Up: platelet activation, NET formation, ECM-receptor, focal adhesion. Down: ribosome (RPL8, RPS8, RPL10A), PD-1 checkpoint, Th1/Th2/Th17, NOD-like, NF-kB | [FT] |
| Pan-cancer features | ADAM10 and TMEM165 up (invasion-related, not tissue-specific) | [FT] |
| CRC-specific human | "Top features for colorectal cancer and stomach cancer are mainly circRNAs". GSEA on tumour-upregulated circRNAs from MiOncoCirc showed "mild but significant" enrichment | [FT] |
| Non-blood tissue | No colon-epithelium, liver or other tissue gene was reported for CRC | [FT] |
| Microbial | CRC-specific genera were reported, but the genus names were stripped from the PMC text | [FT] |

#### 1.2 Other plasma studies

| Study | DOI / PMID | Tissue signal | Tag |
|---|---|---|---|
| Tao et al. Cell Rep Med 2023 (same lab; reuses "our previously published data", n=143 CRC+STAD, likely overlapping Chen) | 10.1016/j.xcrm.2023.101281 / 37992683 | EPIC deconvolution: lymphocytes, macrophages, **endothelial cells and CAFs**. CAF signal was "nearly undetectable" in PBMC. CAF score rises with stage; γδ-T score falls | [FT] |
| PLoS One 2024, CRC plasma total RNA-seq, then capture, then RT-qPCR | 10.1371/journal.pone.0308711 / 39146279 | Validated markers: FGA, FGB, GPR107, CDH3, RP23AP7. FGA and FGB are hepatocyte genes; the liver reading is mine, not the authors' | [Abs] |
| Wang et al. Nat Commun 2024 (37 CRC, 55 GC) | 10.1038/s41467-023-44461-y / 38168054 | Only checked which HPA tissue-elevated genes (intestine, stomach, liver, ...) are detectable. No per-cancer tissue-of-origin result found | [FT] |
| Moss et al. Nat Commun 2018 (cfDNA methylation) | 10.1038/s41467-018-07466-6 / 30498206 | Metastatic colon cancer (n=4): colon was the strongest signal in 3/4. The atlas merges colon and rectum. Healthy gut epithelium probably sheds into the lumen rather than blood | [FT] |
| Liu et al. Ann Oncol 2020 (CCGA, cfDNA methylation) | 10.1016/j.annonc.2020.02.011 / 33506766 | Colorectal, oesophagus, stomach and lung are in the TOO panel; 93% TOO accuracy overall | [Abs] |
| Nat Biotechnol 2025 (cfRNA modifications) | 10.1038/s41587-025-02731-8 / 40629040 | Microbiome-derived cfRNA, not host tissue | [Abs] |

#### 1.3 Pathophysiology candidates (ranked)

| Rank | Tissue | Basis | Citation | Tag |
|---|---|---|---|---|
| 1 | Colon epithelium (primary) | Tumour origin. Colon was detected in cfDNA | Moss 2018 | [FT] |
| 2 | Liver (hepatocytes) | Most common metastatic site. Hepatocyte-like FGA/FGB found in CRC cfRNA | Riihimaki Sci Rep 2016, 10.1038/srep29765 / 27416752; PLoS One 2024 | [Abs] |
| 3 | Peritoneum / mesothelium | Commoner in mucinous and signet-ring tumours (OR 3.8) | Riihimaki 2016 | [Abs] |
| 4 | Lung | Thoracic metastases, especially rectal (OR 2.4) | Riihimaki 2016 | [Abs] |
| 5 | Stroma (CAF, endothelium) | Tumour microenvironment | Tao 2023 | [FT] |
| 6 | Skeletal muscle / adipose | Cachexia occurs in "bowel" cancers | Baracos Nat Rev Dis Primers 2018, 10.1038/nrdp.2017.105 / 29345251 | [Abs] |

### 2. Esophageal cancer (ESCC and EAC)

#### 2.1 Source paper signals

| Paper | Class | Reported | Tag |
|---|---|---|---|
| Chen 2022 | Blood/immune | The pan-cancer pattern in 1.1 | [FT] |
| Chen 2022 | Non-blood tissue | No ESCA-specific human tissue gene is named in the text. One ESCA-enriched microbial genus (name stripped) is "related to several human intestinal diseases" | [FT] |
| Peddu et al. bioRxiv 2025 (nanopore cfRNA; 16 HC, 12 Barrett's with HGD, 19 EAC, mostly stage I) | 10.1101/2025.07.02.662774 (no PMID; not published as of the bioRxiv record) | Up in HGD: CTLA4, LAG3 (immune checkpoint). Pathways in HGD/EAC: MYC targets, OXPHOS, PI3K-AKT-mTOR, mTORC1. HER2, PD-L1, TIGIT, MDM2 about flat. **No tissue-of-origin or deconvolution analysis, and no ESCC**. No suppl. DE table named. GEO GSE271389 (nanopore), GSE271390 (Illumina) | [PREPRINT] (JATS via WebFetch summary) |

#### 2.2 Other plasma studies

| Study | Tissue signal | Tag |
|---|---|---|
| Moss 2018 | Atlas merges stomach and oesophagus. No oesophageal cancer plasma samples | [FT] |
| Liu 2020 CCGA | Oesophagus is in the TOO panel | [Abs] |
| ESCC/gastric cfRNA tissue-of-origin study | **None found** in a PubMed search (2019+). Hits were cfDNA/ctDNA only | search, not exhaustive |

#### 2.3 Pathophysiology candidates

| Rank | Tissue | Basis | Citation | Tag |
|---|---|---|---|---|
| 1a | Squamous oesophageal epithelium (ESCC) | Tumour origin | Textbook | [recall] |
| 1b | Columnar / intestinal-metaplastic epithelium (EAC, Barrett's) | EAC arises from Barrett's metaplasia | Peddu cohort design | [PREPRINT]; biology [recall] |
| 2 | Liver | Most common metastatic site, then lung, bone, brain (SEER n=9,934) | Ai et al. J Thorac Dis 2017, 10.21037/jtd.2017.08.72 / 29221275 | [Abs] |
| 3 | Lung | Second most common metastatic site | Ai 2017 | [Abs] |
| 4 | Bone | Third most common metastatic site | Ai 2017 | [Abs] |
| 5 | Skeletal muscle / adipose | Cachexia is strongly associated with oesophageal cancer | Baracos 2018 | [Abs] |

### 3. Stomach cancer (STAD)

#### 3.1 Source paper signals

| Class | Reported | Tag |
|---|---|---|
| Blood/immune | Pan-cancer pattern (1.1) | [FT] |
| Non-blood tissue | Top features mostly circRNAs, mildly enriched for tumour-upregulated circRNAs. No gastric-epithelium gene named | [FT] |
| Microbial | One STAD-specific genus (name stripped), "enriched in oral cancer" | [FT] |
| Tao 2023 (overlapping cohort) | CAF and endothelial signals. CAF score correlates with stage; low CAF score predicts better TCGA survival | [FT] |

#### 3.2 Pathophysiology candidates

| Rank | Tissue | Basis | Citation | Tag |
|---|---|---|---|---|
| 1 | Gastric epithelium (pit/foveolar, chief, parietal, neck cells) | Tumour origin | Textbook | [recall] |
| 2 | Liver | 48% of metastatic patients | Riihimaki Oncotarget 2016, 10.18632/oncotarget.10740 / 27447571 | [Abs] |
| 3 | Peritoneum / mesothelium | 32%, mainly non-cardia and signet-ring | Riihimaki 2016 | [Abs] |
| 4 | Lung | 15% | Riihimaki 2016 | [Abs] |
| 5 | Bone | 12%. Also ovary in signet-ring tumours | Riihimaki 2016 | [Abs] |
| 6 | Stroma (CAF) | Shown in plasma | Tao 2023 | [FT] |
| 7 | Skeletal muscle / adipose | Cachexia in stomach cancer | Baracos 2018 | [Abs] |

### 4. Lung cancer

#### 4.1 Source paper signals

| Class | Reported | Tag |
|---|---|---|
| Blood/immune | Lung-specific features IL1R2 and CLEC4E, both immune-regulatory (myeloid). Plus the pan-cancer pattern | [FT] |
| Non-blood tissue | No lung-epithelium gene reported | [FT] |

#### 4.2 Other plasma studies (cfRNA and cfDNA)

| Study | DOI / PMID | Tissue signal | Tag |
|---|---|---|---|
| Larson et al. Nat Commun 2021 (stage III lung, n=30) | 10.1038/s41467-021-22444-1 / 33883548 | 12 lung "dark channel" genes, 6/12 HPA lung-specific. Named: **SFTA3, SFTPA2, SLC34A2, BPIFA1** (the other 8 names were stripped from the PMC text). Found in LUAD, not LUSC. Deconvolution: lung-tissue transcripts present in cancer plasma | [FT] |
| Nesselbush et al. Nature 2025 (RARE-seq; 139 LUAD) | 10.1038/s41586-025-08834-1 / 40240612 | GSEA: LUAD cfRNA "most enriched for ... pulmonary alveolar type II cells". **Normal lung RNA in COPD, COVID-19, ARDS, and tobacco-exposed controls** ("lung cfRNA is induced by ongoing injury ... tobacco smoke"). TOO panel: LIHC, LUAD, PAAD, PRAD only | [FT] |
| Moss 2018 | 10.1038/s41467-018-07466-6 / 30498206 | Lung cancer (n=4): lung was the strongest signal in 2/4 | [FT] |
| J Immunother Cancer 2026 (EV cfRNA, advanced NSCLC) | 10.1136/jitc-2026-015873 / 42875050 | ITGB3 (platelet integrin) prognostic. ITGB3-low patients had more immune (monocyte, MAIT) fraction | [Abs] |
| Diagnostics 2022 pilot (NSCLC n=12) | 10.3390/diagnostics12122897 / 36552904 | 251 of 2357 DE cfRNAs were previously reported in primary NSCLC; no tissue named | [Abs] |

#### 4.3 Pathophysiology candidates

| Rank | Tissue | Basis | Citation | Tag |
|---|---|---|---|---|
| 1 | Lung epithelium: AT2 (LUAD), basal/squamous (LUSC), club, AT1 | Tumour origin. AT2 already reported in plasma | Larson 2021; Nesselbush 2025 | [FT] |
| 2 | Bone | 39% of metastatic adenocarcinoma | Riihimaki Lung Cancer 2014, 10.1016/j.lungcan.2014.07.020 / 25130083 | [Abs] |
| 3 | Brain / nervous system | Most frequent site overall (47% in SCLC) | Riihimaki 2014 | [Abs] |
| 4 | Liver | Common in SCLC (35%) | Riihimaki 2014 | [Abs] |
| 5 | Adrenal | Among the most frequent sites | Riihimaki 2014 | [Abs] |
| 6 | Skeletal muscle / adipose | Cachexia in lung cancer | Baracos 2018 | [Abs] |
| - | Paraneoplastic (e.g. SIADH, PTHrP hypercalcaemia) | Not verified (searches returned noise). Most of our lung cases are early-stage NSCLC, so low priority | - | [recall] |

**Lung confounder:** smoking and benign lung injury raise lung cfRNA (Nesselbush [FT]). An AT2/AT1 deviation in a lung cancer patient is therefore not tumour-specific unless the healthy-control smoking status is known.

### 5. Final tissue-of-interest list

**Marker genes:** all are [recall], canonical markers. Check each against HPA single-cell-type enrichment before use.
**Cell-type names:** "HPA" = HPA single-cell-type name; "Panglao" = PanglaoDB term. Both are [recall]; confirm the exact strings in the database.

| Disease | Tissue | Cell types (HPA / Panglao) | Markers (5-10) | Rationale | Already reported in cfRNA lit? |
|---|---|---|---|---|---|
| CRC | Colon epithelium | HPA: Distal enterocytes, Intestinal goblet cells. Panglao: Enterocytes, Goblet cells | CDX1, CDX2, CDH17, VIL1, LGALS4, KRT20, SATB2, GUCA2A, CEACAM5, MUC2 | Primary organ | **No** (Chen: circRNAs only [FT]; Wang 2024: detectability only [FT]). cfDNA yes (Moss [FT]) |
| CRC | Liver | HPA: Hepatocytes. Panglao: Hepatocytes | ALB, APOA1, APOH, FGA, FGB, FGG, HP, SERPINA1, TF, AHSG | Top metastatic site | **Partly**: FGA/FGB as markers, not framed as liver (PLoS One 2024 [Abs]) |
| CRC | Peritoneum | HPA: Mesothelial cells (unsure the type exists). Panglao: Mesothelial cells | MSLN, UPK3B, CALB2, WT1, PRG4, ITLN1 | Mucinous/signet spread | No |
| CRC, STAD | CAF / stroma | HPA: Fibroblasts. Panglao: Fibroblasts | COL1A1, COL1A2, COL3A1, DCN, LUM, FAP, PDGFRB, POSTN | TME | **Yes**: Tao 2023 [FT] (overlapping cohort) |
| ESCC | Squamous oesophageal epithelium | HPA: Squamous epithelial cells, Basal squamous epithelial cells. Panglao: Keratinocytes (closest) | KRT4, KRT13, KRT5, TP63, CRNN, SPRR3, S100A8* | Primary (ESCC) | No |
| EAC | Barrett's / columnar metaplasia | HPA: Gastric mucus-secreting cells, Intestinal goblet cells. Panglao: Foveolar cells, Goblet cells | TFF1, TFF2, TFF3, MUC5AC, AGR2, CDX2, MUC2, REG4 | Primary (EAC) | No (Peddu: no tissue analysis [PREPRINT]) |
| STAD | Gastric epithelium | HPA: Gastric mucus-secreting cells (+ chief/parietal if present). Panglao: Gastric chief cells, Parietal cells, Foveolar cells | PGC, PGA3, LIPF, GKN1, GKN2, TFF1, MUC5AC, ATP4A, ATP4B, CLDN18 | Primary organ | No |
| STAD, ESCA, CRC | Liver | as above | as above | 48% of gastric metastases; top site in oesophageal cancer | **No** for STAD and ESCA |
| STAD | Peritoneum | as above | as above | 32% of metastases | No |
| Lung | AT2 | HPA: Alveolar cells type 2. Panglao: Pulmonary alveolar type II cells | SFTPC, SFTPB, SFTPA1, SFTPA2, NAPSA, LAMP3, SLC34A2, SFTA3 | Primary (LUAD) | **Yes**: Larson [FT], Nesselbush [FT] |
| Lung | AT1 / club / basal | HPA: Alveolar cells type 1, Club cells, Respiratory basal cells. Panglao: Pulmonary alveolar type I cells, Clara cells, Basal cells | AGER, EMP2, SCGB1A1, SCGB3A2, BPIFA1, KRT5, TP63, SOX2 | Injury and LUSC | Partly: BPIFA1 (Larson); normal lung in injury (Nesselbush) |
| Lung | Bone | HPA: (no osteoblast type; unsure). Panglao: Osteoblasts, Osteoclasts | BGLAP, IBSP, SP7, RUNX2, COL1A1, ACP5, CTSK | 39% of metastatic LUAD | No |
| Lung | Brain | HPA: Astrocytes, Oligodendrocytes, Excitatory neurons. Panglao: Astrocytes, Oligodendrocytes | GFAP, AQP4, MBP, PLP1, SNAP25, NEFL | Most frequent site (nervous system) | No |
| Lung | Adrenal cortex | HPA: (tissue-level only; unsure). Panglao: Adrenocortical cells (unsure) | STAR, CYP11A1, CYP11B1, CYP17A1, CYP21A2, NR5A1, MC2R | Frequent site | No |
| All four | Skeletal muscle | HPA: Skeletal myocytes. Panglao: Myocytes | ACTA1, CKM, MB, MYH1, MYH2, MYH7, TNNT1, TNNT3, DES | Cachexia | No |
| All four | Adipose | HPA: Adipocytes. Panglao: Adipocytes | ADIPOQ, LEP, PLIN1, PLIN4, CIDEA, CIDEC, LIPE | Cachexia | No |

*S100A8 is neutrophil-high: drop it to avoid blood confounding.
**General caveat:** COL1A1, MMP9, ACP5 and SPP1 are also expressed by blood or myeloid cells, so test non-blood specificity in HC plasma first.

### 6. Key gaps / uncertainties

| Point | Evidence | Uncertainty | Verify manually |
|---|---|---|---|
| Chen 2022 reports no non-blood tissue genes for CRC, STAD, ESCA or lung | [FT], main text | Supplementary DEG tables were not read | eLife Supplementary files and Figure 5c gene list |
| Lung subtype and oesophageal histology in our samples are unknown | [FT] gives no histology | AT2 markers are expected only in LUAD | GEO / SRA metadata |
| GI liver signal is untested in cfRNA (except FGA/FGB) | [Abs] | The PLoS One authors did not interpret it as liver | PLoS One 2024 full text |
| Lung tissue cfRNA rises with smoking and injury | Nesselbush [FT] | Smoking status of the HC pool is unknown | HC metadata |
| Tao 2023 and Chen 2022 cohorts overlap | Tao [FT] says "our previously published data" | Extent of overlap unclear | Tao STAR Methods |

---

## 3. Liver Cancer (HCC), Liver Cirrhosis, Multiple Myeloma, MGUS: tissues of interest

Tags: [FT] full text read via PMC, [Abs] abstract only, [recall] unverified memory, [PREPRINT]. "A file" = scratchpad/lit/A_cfrna_paradigms.md (already covered there, not repeated).
All papers checked through PubMed MCP (2026-10-27).

### 0. Verified references

| Key | Paper | DOI | PMID | Read |
|---|---|---|---|---|
| RH22 | Roskams-Hieter, NPJ Precis Oncol 2022 (OHSU/CEDAR; HCC, Cirr, MM, MGUS, NC) | 10.1038/s41698-022-00270-y | 35468987 | [FT] main text, supplement not opened |
| CH22 | Chen, eLife 2022 (5 cancers incl. liver, n~300) | 10.7554/eLife.75181 | 35816095 | [FT] |
| BL22 | Block, Front Oncol 2022 (HCC 14, LC 8, NHC 6; SRA PRJNA907745) | 10.3389/fonc.2022.963641 | 36582804 | [FT] |
| SA20 | Sayeed, Oncotarget 2020 (same Blumberg group, HCC/LC cf-mRNA) | 10.18632/oncotarget.27617 | 32577166 | [Abs] |
| ZE25 | Zezulinski, Liver Cancer 2025 (targeted ctmutRNA, HCC 50 / LC 35) | 10.1159/000545366 | 40331063 | [Abs] |
| ZH21 | Zhu, Theranostics 2021 (Lu lab; 35 liver cancer + 30 HD cfRNA reused in CH22) | 10.7150/thno.48206 | 33391469 | [Abs] |
| KI23 | Kim, Commun Biol 2023 (OHSU; EV-fractionated cf-mRNA, liver ca, MM, lung ca) | 10.1038/s42003-023-05232-z | 37644220 | [Abs] |
| IB20 | Ibarra, Nat Commun 2020 (BM cf-mRNA; MM/AML transplant) | 10.1038/s41467-019-14253-4 | 31964864 | [Abs] |
| SU15 | Sun, PNAS 2015 (plasma DNA methylation tissue mapping; HCC) | 10.1073/pnas.1508736112 | 26392541 | [Abs] |
| GE21 | GETMap, eLife 2021 (HCC tumor cfDNA origin) | 10.7554/eLife.64356 | 33752803 | [Abs] |
| MO18 | Moss, Nat Commun 2018 (cfDNA methylation atlas) | 10.1038/s41467-018-07466-6 | 30498206 | [Abs] |
| LW18 | Lehmann-Werman, JCI Insight 2018 (hepatocyte cfDNA) | 10.1172/jci.insight.120687 | 29925683 | [Abs] |
| LU22 | Lubotzky, JCI Insight 2022 (collateral tissue damage, hepatocyte cfDNA) | 10.1172/jci.insight.153559 | 35076021 | [Abs] |
| RA19 | Ramachandran, Nature 2019 (cirrhosis scRNA fibrotic niche) | 10.1038/s41586-019-1631-3 | 31597160 | [Abs] |
| MP18 | MacParland, Nat Commun 2018 (human liver scRNA atlas) | 10.1038/s41467-018-06318-7 | 30348985 | [Abs] |
| GI21 | Gines, Lancet 2021 (cirrhosis review) | 10.1016/S0140-6736(21)01374-X | 34543610 | [Abs] |
| RJ14 | Rajkumar, Lancet Oncol 2014 (IMWG diagnostic criteria, CRAB) | 10.1016/S1470-2045(14)70442-5 | 25439696 | PMID resolved by DOI search; content [recall] |
| TE21 | Myeloma bone disease review, Int J Mol Sci 2021 | 10.3390/ijms22126208 | 34201396 | [Abs] (low-IF review) |
| MBD17 | Myeloma and bone disease, Curr Osteoporos Rep 2017 | 10.1007/s11914-017-0397-5 | 28861842 | [Abs] |
| NG19 | Blood 2011 (author not checked), MGUS bone microstructure, DKK1 and MIP-1a | 10.1182/blood-2011-04-351437 | 22042700 | [Abs] |
| IK19 | IKMG MGRS consensus, Nat Rev Nephrol 2019 | 10.1038/s41581-018-0077-4 | 30510265 | [Abs] |
| EH25 | EHA-EMN MM guidelines, Nat Rev Clin Oncol 2025 | 10.1038/s41571-025-01041-x | 40624367 | [Abs] |
| HM19 | HCC seminal vesicle metastasis case report, Medicine 2019 (lists usual extrahepatic sites) | 10.1097/MD.0000000000013770 | 30653089 | [Abs] (weak source) |

---

### A. Liver Cancer (HCC)

#### A1. Source papers of our data

| Paper | Method | Blood/immune signal | Non-blood tissue signal | Supplementary DEG tables |
|---|---|---|---|---|
| RH22 | SMARTer pico v2 total RNA; DESeq2 (Wald, padj<0.01): 12 DEGs HCC vs NC; LVQ top 10; RF/LDA; n=8 HCC pilot, 20 validation [FT] | None named for HCC | Liver: 8 of the top 10 LVQ genes are liver-specific secreted proteins (APOE, C3, CP, DHCR24, FGA, FGB, FGG, HRG). GO up: plasminogen activation [FT] | "Supplementary Table" (DESeq2 + LVQ lists) on the journal page; file names/URL not opened, unsure |
| CH22 | SMART-total, edgeR QL (abs(log2FC)>1, FDR<0.05), RUVg, balanced RF [FT] | Pan-cancer up: platelet activation, NET formation; down: ribosome, PD-1, Th1/Th2/Th17 (immune suppression) [FT] | Liver-cancer-specific: TF, HRG, CP, FGA ("specifically expressed in the liver") [FT]. Viral: Orthohepadnavirus (HBV), TTVs [FT] | Fig 3-5 source data / supplementary files on eLife; exact file names not checked |
| BL22 | Total RNA-seq, GATK/SnpEff variant calling, plasma + EV + tumor [FT] | None (no DE analysis) | Variant load in liver transcripts: FGL1, CYP2E1, SERPINA1, ALB, TF; complement (C2, C4A/B, C1R, C5, C8G, KLKB1, KNG1, CFHR5, MASP2), F2, F11, F12 [FT]. Expression-level DE not reported | Supplementary tables of variants; data at SRA PRJNA907745 [FT] |

#### A2. Other cfRNA / cfDNA tissue-of-origin studies (HCC)

| Paper | Analyte | Tissue signal reported | Tag |
|---|---|---|---|
| SA20 | Plasma cf-mRNA | Up in HCC: ALB, APOA1, APOA2, APOH, SERPINA1, SERPINE1, FTL, FGL1. Liver mRNA mostly outside EVs | [Abs] |
| ZH21 | exoRNA + cfRNA | Liver cancer noncoding RNA panel; tissue mapping not stated in the abstract | [Abs] |
| KI23 | EV-fractionated cf-mRNA | EV-enriched liver cancer genes linked to "liver function" | [Abs] |
| ZE25 | Targeted ctmutRNA | Variant panel, fusions (SENP7, HYI, SAR1A, RASA2, TUBA); no tissue mapping | [Abs] |
| SU15 | cfDNA methylation | Elevated liver contribution in HCC plasma | [Abs] |
| GE21 | cfDNA methylation + variants | Liver is the predominant source of tumor-variant cfDNA in HCC | [Abs] |
| LU22 | cfDNA methylation | Hepatocyte cfDNA up in liver metastasis; HCC itself not the focus | [Abs] |
| Moore 2025, Vorperian 2022 (NAFLD), Chalasani 2023 (PSC: hepatocyte, HSC, Kupffer) | cfRNA | See A file | A file |

#### A3. Pathophysiology-based candidate tissues (HCC); rank 1 = most plausible

| Rank | Tissue / cell | Rationale | Citation |
|---|---|---|---|
| 1 | Liver: hepatocytes (normal + malignant) | Tumor and peritumoral hepatocyte turnover; already the dominant reported signal | RH22, CH22 [FT]; SA20, SU15, GE21 [Abs] |
| 2 | Liver: tumor-specific/oncofetal hepatocyte program | Distinguishes tumor from background-cirrhosis hepatocytes (AFP, GPC3 etc.) | Markers [recall]; not reported in RH22/CH22 text [FT] |
| 3 | Liver: stellate cells / myofibroblasts, LSEC, Kupffer, cholangiocytes (ductular reaction) | Most HCC arise on cirrhosis (11/14 in BL22 [FT]); fibrotic niche cells expand (TREM2+CD9+ macrophages, PLVAP+/ACKR1+ endothelium) | RA19, MP18 [Abs]; BL22 [FT] |
| 4 | Lung, lymph node, adrenal, bone (metastasis) | "Most common sites of extrahepatic metastases are lung, regional lymph node, adrenal gland, and bone" | HM19 [Abs] (case report; frequencies not verified) |
| 5 | Spleen / gut / kidney | Via portal hypertension in background cirrhosis; see B3 | GI21 [Abs] |

### B. Liver Cirrhosis

#### B1. Source papers of our data

| Paper | Method | Blood/immune signal | Non-blood tissue signal | Supplementary |
|---|---|---|---|---|
| RH22 | LVQ top 10 Cirr vs NC; RF LOOCV 100% (n=4 Cirr pilot only) [FT] | Not named in the main text | Gene names for Cirr vs NC are in the supplement only (not opened). Main text names liver genes for HCC vs NC only [FT] | Supplementary Table, unsure of contents |
| BL22 | Variant calling; LC plasma used as the negative comparator [FT] | None | No expression signal reported for LC. LC EVs higher CD63 (n=1) [FT] | PRJNA907745 |

#### B2. Other studies (cirrhosis / chronic liver disease)

| Paper | Tissue signal | Tag |
|---|---|---|
| SA20 | Circulating mRNA "similar from person to person" across HCC and LC; LC-specific genes not stated in the abstract | [Abs] |
| Vorperian 2022 (NAFLD), Chalasani 2023 (PSC: hepatocyte, HSC, Kupffer, subset liver cluster) | cfRNA cell-type scores in chronic liver disease | A file |
| LW18 | Hepatocyte cfDNA rises with acute hepatocyte death (transplant, hepatectomy, sepsis); cirrhosis itself not tested in the abstract | [Abs] |
| Direct cirrhosis cfRNA tissue-of-origin study | None found by PubMed search (2 queries) | unsure, may exist |

#### B3. Pathophysiology-based candidate tissues (cirrhosis)

| Rank | Tissue / cell | Rationale | Citation |
|---|---|---|---|
| 1 | Liver: activated stellate cells / scar mesenchyme (PDGFRA+) | Fibrosis drives cirrhosis; mesenchymal expansion in the fibrotic niche | RA19 [Abs], GI21 [Abs] |
| 2 | Liver: scar-associated endothelium (PLVAP+, ACKR1+) and loss of LSEC phenotype | Expand in cirrhosis, restricted to the fibrotic niche | RA19 [Abs]; LSEC capillarisation [recall] |
| 3 | Liver: hepatocytes (direction ambiguous) | Inflammation and regenerative nodules raise turnover; synthetic failure may lower secreted-protein mRNA | GI21 [Abs]; direction [recall] |
| 4 | Liver: cholangiocytes / ductular reaction | Ductular reaction scales with fibrosis | [recall] |
| 5 | Liver: Kupffer and scar-associated macrophages (TREM2+CD9+, monocyte-derived) | Expand in fibrosis. Note: monocyte-derived, so overlaps with blood myeloid | RA19 [Abs] |
| 6 | Spleen (red pulp, littoral cells) | Portal hypertension causes splenomegaly / hypersplenism; most splenic cells are blood lineage, so specific markers are weak | GI21 [Abs]; splenomegaly [recall] |
| 7 | Gut epithelium | Portal hypertensive enteropathy, bacterial translocation and systemic inflammation | GI21 [Abs] ("systemic inflammation"); translocation [recall] |
| 8 | Kidney (hepatorenal syndrome) | Mainly in decompensated cirrhosis; functional, little tubular death, so the cfRNA signal is probably weak | [recall] |

### C. Multiple Myeloma and MGUS

#### C1. Source paper of our data

| Paper | Method | Blood/immune signal | Non-blood tissue signal | Supplementary |
|---|---|---|---|---|
| RH22 | DESeq2: 110 DEGs MM vs NC; GO up: oxygen transport / gas transport. LVQ top 10 MM vs NC includes CA1, EPB42, HBG1, HBG2, CENPE, CPOX, NUSAP1, NEK2. "7/10 higher expression in bone marrow" (HPA). Gradient NC to MGUS to MM [FT] | All named genes are erythroid (CA1, EPB42, HBG1/2, CPOX) or proliferation (CENPE, NUSAP1, NEK2). The authors call this "bone marrow", but it is hematopoietic, not plasma-cell specific [FT; interpretation mine] | None reported [FT] | Supplementary Table (DESeq2 110 genes, LVQ MGUS vs NC); not opened |

#### C2. Other studies (MM / MGUS)

| Paper | Tissue signal | Tag |
|---|---|---|
| IB20 | cf-mRNA is enriched for bone-marrow transcripts relative to circulating cells; tracks BM reconstitution after transplant in MM/AML (hematopoietic lineages) | [Abs] |
| KI23 | MM EV-enriched genes linked to "toxic substance regulation" | [Abs] |
| MM cfDNA-methylation tissue-of-origin (bone, kidney) | None found (PubMed query returned 0) | unsure |
| MM cfRNA with bone or kidney markers | None found | unsure |

#### C3. Pathophysiology-based candidate tissues (MM / MGUS)

Plasma cells and erythroid precursors are blood lineage. The plausible non-blood tissues are bone (osteoblast, osteocyte, osteoclast niche), BM stroma, and kidney.

| Rank | Tissue / cell | Rationale | MGUS relevance | Citation |
|---|---|---|---|---|
| 1 (blood) | BM clonal plasma cells | Disease cell; 80-90% of disease occurs in BM | Yes (lower burden) | TE21 [Abs], EH25 [Abs] |
| 2 (blood) | BM erythroid precursors | Anemia (CRAB "A"); already the dominant reported signal | Weak | RH22 [FT], RJ14 [recall] |
| 3 (non-blood) | Bone: osteoclasts (up) and osteoblasts (suppressed) | Osteolytic lesions in ~80% of MM; increased osteoclast activity with osteoblast inhibition | MGUS: higher fracture risk, altered microstructure, raised DKK1 / MIP-1a (CCL3) | TE21, MBD17, NG19 [Abs] |
| 4 (non-blood) | BM stroma: mesenchymal stromal/CAR cells, adipocytes | Myeloma niche; supports plasma cells | Possible | [recall] |
| 5 (non-blood) | Kidney: proximal tubule, distal tubule / TAL (cast nephropathy) | CRAB "R"; light-chain toxicity to tubules; casts form with uromodulin | MGUS: MGRS only (rare; needs biopsy) | IK19 [Abs]; frequency and UMOD role [recall] |
| 6 (non-blood) | Heart, liver, other organs (AL amyloid, extramedullary disease) | Rare subgroups | Rare | [recall] |

---

### D. Final tissue-of-interest list

HPA = HPA single-cell-type name; PDB = PanglaoDB cell-type name. Exact name strings are [recall], so check them against the downloaded tables. Markers are [recall] canonical markers; check specificity in HPA/GTEx before use.

#### D1. HCC and cirrhosis

| Tissue | Cell type (HPA / PDB) | Markers | Rationale | Already reported in cfRNA lit? |
|---|---|---|---|---|
| Liver: hepatocyte (secreted proteins) | Hepatocytes / Hepatocytes | ALB, APOA1, APOH, FGA, FGB, FGG, HRG, CP, TF, SERPINA1 | Main reported signal; positive control for NM | Yes: RH22, CH22 [FT]; SA20 [Abs]; Moore 2025 (PDAC) A file |
| Liver: hepatocyte, oncofetal/tumor | Hepatocytes (tumor state; no separate HPA/PDB type) | AFP, GPC3, DLK1, IGF2, H19, AKR1B10, SPINK1, MDK | Separates tumor from background liver | No in RH22/CH22 main text [FT]; SA20 supplement not checked |
| Liver: hepatocyte, CYP/zonation (loss of function) | Hepatocytes | CYP2E1, CYP3A4, CYP1A2, GLUL, OAT, HAMP | Down-direction deviation (synthetic/metabolic failure) is testable per patient | Partly: CYP2E1 only as a variant locus in BL22 [FT]; expression not reported |
| Liver: cholangiocyte / ductular reaction | Cholangiocytes / Cholangiocytes | KRT19, KRT7, SOX9, EPCAM, SPP1, CFTR, ANXA4, AQP1 | Ductular reaction in cirrhosis; mixed HCC-CCA | No (CCA only as variants in BL22 [FT]) |
| Liver: stellate / scar mesenchyme | Ito cells (HPA name, unsure) / Hepatic stellate cells | PDGFRB, PDGFRA, RGS5, LRAT, DCN, LUM, COL1A1, COL3A1, ACTA2 | Fibrosis driver | No for HCC/cirrhosis; yes for PSC (Chalasani 2023, A file) |
| Liver: LSEC and scar endothelium | Endothelial cells (no LSEC-specific HPA type, unsure) / Endothelial cells | STAB2, CLEC4G, CLEC4M, CRHBP, FCN2, FCN3, OIT3, PLVAP, ACKR1 | LSEC capillarisation; PLVAP+/ACKR1+ expansion (RA19) | No |
| Liver: Kupffer / scar-associated macrophage | Kupffer cells / Kupffer cells | CD5L, MARCO, VSIG4, TIMD4, CD163, TREM2, CD9 | Tissue-resident macrophage; partly overlaps with monocytes | No for HCC/cirrhosis; yes for PSC (A file) |
| Lung (metastasis) | Alveolar cells type 2 / Pneumocytes | SFTPC, SFTPB, SFTPA1, NAPSA, AGER | Most common extrahepatic site | No |
| Bone (HCC metastasis) | see D2 bone row | see D2 | Bone metastasis | No |
| Spleen (cirrhosis) | Red pulp: no clean non-blood HPA type, unsure | STAB2, CD8A (littoral), unsure | Portal hypertension; low specificity | No |
| Gut epithelium (cirrhosis) | Enterocytes / Enterocytes | FABP2, FABP6, APOA4, RBP2, ALPI, MUC2 | Portal hypertensive enteropathy / translocation | No |
| Kidney (hepatorenal) | Proximal tubular cells / Proximal tubule cells | LRP2, CUBN, SLC34A1, UMOD, HAVCR1 | Low prior | No |

#### D2. MM and MGUS

| Tissue | Cell type (HPA / PDB) | Markers | Rationale | Already reported in cfRNA lit? |
|---|---|---|---|---|
| BM plasma cell (blood) | Plasma cells / Plasma cells | JCHAIN, MZB1, SDC1, TNFRSF17, XBP1, PRDM1, IGKC, IGLC2, TXNDC5 | Disease cell; Ig genes are dominated by the clone (kappa/lambda skew is a per-patient signal) | Not among the top 10 named in RH22 [FT]; DEG supplement not checked |
| BM erythroid (blood) | Erythroid cells / Erythroid-like and erythroid precursor cells | HBB, HBG1, HBG2, CA1, ALAS2, EPB42, SLC4A1, CPOX | Anemia / marrow; already the main reported signal, so it is the positive control | Yes: RH22 [FT]; IB20 [Abs] |
| Bone: osteoclast | Osteoclasts (HPA type unsure) / Osteoclasts | ACP5, CTSK, MMP9, OSCAR, DCSTAMP, CALCR, NFATC1 | Lytic bone disease. Hematopoietic origin but bone-resident; ACP5/CTSK/MMP9 are also in macrophages, so CALCR, OSCAR and DCSTAMP are more specific | No |
| Bone: osteoblast / osteocyte | Osteoblasts (HPA unsure) / Osteoblasts, Osteocytes | RUNX2, SP7, BGLAP, IBSP, ALPL, SPP1, COL1A1, SOST, DMP1 | Suppressed in MM (direction: down or no change), raised DKK1 in MGUS (NG19) | No |
| BM stroma / niche | Fibroblasts? (HPA unsure) / Mesenchymal stem cells, Adipocytes | CXCL12, LEPR, KITLG, VCAM1, ADIPOQ, DKK1 | Myeloma niche. Note: DKK1 is myeloma-cell secreted | No |
| Kidney: proximal tubule | Proximal tubular cells / Proximal tubule cells | LRP2, CUBN, SLC34A1, SLC22A6, HAVCR1, ALDOB (shared with liver) | Light-chain proximal tubulopathy, CRAB "R" | No |
| Kidney: distal tubule / TAL | Distal tubular cells / Loop of Henle cells | UMOD, SLC12A1, SLC12A3, CALB1 | Cast nephropathy site | No |

---

### E. Caveats

| Evidence source | Evidence strength | Main uncertainty | What to verify manually |
|---|---|---|---|
| RH22 liver/MM gene lists [FT] | Moderate (small n: HCC 8, Cirr 4, MM 10, MGUS 8 pilot) | Cirr and MGUS gene names are only in the supplement | Download the RH22 Supplementary Table and extract Cirr vs NC and MGUS vs NC LVQ/DE genes |
| RH22 MM "bone marrow" claim | Weak as tissue claim | Genes are erythroid, so "bone marrow" is the hematopoietic compartment, not plasma cell or bone | Check whether plasma-cell (JCHAIN, MZB1) genes appear in the 110-gene DEG list |
| CH22 liver genes [FT] | Moderate (RUVg across 6 centers) | Only 4 liver genes named | Check the eLife Fig 5C source data |
| Non-blood MM tissues (bone, kidney) | Clinical only (reviews [Abs]) | No cfRNA/cfDNA study found showing a bone or kidney signal in MM | Osteoclast/osteoblast transcripts may be undetectable in plasma; check HC detection rate before testing |
| Cirrhosis secondary organs (spleen, gut, kidney) | Weak | Mostly [recall]; specific cell markers unclear | Treat as exploratory only |
| Marker lists and HPA/PDB names | [recall] | Specificity varies (COL1A1, SPP1, CD163, ACP5 are shared with other cell types) | Cross-check against the HPA single-cell-type and PanglaoDB tables |

---

## 4. TB / HIV / HIV+TB: candidate tissues of interest (literature survey, 2026-10-27)

Sources were verified through PubMed metadata. Full text came from PMC / Europe PMC XML.
Tags: [FT] full text read; [Abs] abstract only; [recall] not verified, from memory; [PREPRINT] not peer reviewed.
Papers already covered in `scratchpad/lit/A_cfrna_paradigms.md` are marked "(see A)".

### 0. Critical design note (affects every table below)

| Point | Evidence | Tag | Implication for NM |
|---|---|---|---|
| All 251 Chang participants were adults with cough of 2 weeks or more. "TB negative" means symptomatic people with other respiratory conditions, not healthy people | "251 individuals with a cough lasting at least two weeks"; the TB-negative group "reflect[s] the local distribution of other conditions that also present with a cough" | [FT] Chang 2024 | Chang's DEGs compare TB with non-TB respiratory illness. Any lung or inflammation signal shared by both groups is invisible in Chang's DE, but NM against healthy controls could detect it |
| Gaudio 2026 describes the same cohort as "active tuberculosis versus symptomatic respiratory controls" | Abstract | [Abs] Gaudio 2026, Nat Commun, 10.1038/s41467-026-74077-x / PMID 42276999 | Independent confirmation of the point above |
| If our "Chang HC" are Chang's TB-negative/HIV-negative participants, they are symptomatic, not healthy | Project pilot found a lung signal in Chang HC | Not verified | **Check manually** how the project labels Chang TB-/HIV- samples (HC pool vs a separate group). This could explain the pilot's lung signal biologically, not only as batch |
| Chang includes no HIV-vs-healthy comparison. HIV is analysed only as a modifier of classifier performance | 7 of 8 misclassified HIV+ samples were HIV+/TB-; the signature kept accuracy in HIV+/TB+ (1 of 18 misclassified) | [FT] Chang 2024 | For HIV-only (TB-) samples there is no published cfRNA DEG list, tissue or blood |

### 1. Source paper: Chang et al.

Preprint: medRxiv 2023, 10.1101/2023.01.11.23284433, PMID 36711999, n=182 [PREPRINT][Abs].
Published version: Nat Commun 2024;15:4949, 10.1038/s41467-024-49245-6, PMID 38858368, PMC11164910, n=251 [FT]. (see A)

#### 1a. Method

| Item | Detail | Tag |
|---|---|---|
| Cohorts | C1 Uganda: Streck tubes plus an extra 16,000 g spin. C2 Uganda/Vietnam/Philippines and C3 Uganda: K2EDTA tubes, single 1600 g spin | [FT] |
| Size | 142 TB+, 37 HIV+ (18 HIV+/TB+) | [FT] |
| Library | Norgen extraction, SMARTer Stranded Pico v3, NextSeq 500/2000. Reads trimmed to 61 bp. Mito, ribosomal, X and Y genes removed | [FT] |
| DE | DESeq2 TB+ vs TB- (all 251 samples), BH < 0.05, giving 1956 genes. Then IPA, plus KEGG/GO/Hallmark GSEA | [FT] |
| Deconvolution | BayesPrism with Tabula Sapiens v1 (grouped as in Vorperian 2022). Used descriptively only, with no TB vs non-TB fraction test reported | [FT] |
| Classifier | Greedy forward search giving 6 genes: GBP5, BNIP3L, KLF6, DYSF, LASP1, PCBP1 | [FT] |
| Data | GEO GSE255071, GSE255073, GSE255074 | [FT] |

#### 1b. Reported genes and tissue signals

| Compartment | Genes or signal reported (TB+ vs TB-) | Tag |
|---|---|---|
| Blood/immune: "macrophage markers" | MARCO, SOCS3, FCGR1A, MPO, C1QB | [FT] |
| Blood/immune: "neutrophil markers" | ELANE, FCGR3A, S100A8, S100A9, ERG | [FT] |
| Blood/immune: interferon | IFI27L1, IFIT2, IFIT3, IFITM3, IRF1 | [FT] |
| Blood/immune: antimicrobial | AZU1, CTSG, DEFA4, STAT1, GBP1, GBP2, GBP4, GBP5, GBP6 | [FT] |
| Blood/immune: pathways | Neutrophil degranulation, IFN-gamma signalling, antimicrobial peptides | [FT] |
| **Non-blood: lung** | "lung-specific markers (SFTPC, SFTPB, SLC34A2, SFTPA1, SFTPA2)" are elevated in TB. These are alveolar type 2 cell genes | [FT] |
| Non-blood: deconvolution | Main sources across all samples: platelets, myeloid progenitors, B cells, endothelial cells, NK cells, monocytes, neutrophils. Solid-organ fraction 0.23 / 0.39 / 0.33 by cohort (whole blood: 0.025). No organ-level TB vs non-TB result | [FT] |
| Batch | Platelet fraction 0.11 / 0.61 / 0.36 by cohort, attributed to tube type and processing | [FT] |
| Oddity | ERG is listed as a "neutrophil marker". ERG is usually an endothelial transcription factor [recall] | Verify in Supp Data 1 |

#### 1c. Supplementary files (file names from the PMC XML [FT]; contents mapped from the main text)

| File | Content according to the main text |
|---|---|
| 41467_2024_49245_MOESM4_ESM.xls ... MOESM8_ESM.xls, MOESM10_ESM.xlsx | Supplementary Data 1 = 1956 DEGs (TB+ vs TB-). Data 2 = IPA pathways. A further Data file = KEGG/GO/Hallmark GSEA. A further Data file = wbRNA vs cfRNA DEGs, plus overlap with published wbRNA signatures. The exact mapping of MOESM numbers to Data 1-5 is **unverified** |
| Base URL (pattern, [recall]) | https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-49245-6/MediaObjects/41467_2024_49245_MOESM4_ESM.xls |

### 2. Other plasma/serum cfRNA or cfDNA tissue-of-origin studies in TB/HIV

| Paper | Journal/Year | DOI / PMID | Finding relevant to tissues | Tag |
|---|---|---|---|---|
| Wu et al. | iScience 2026 | 10.1016/j.isci.2026.115178 / 41858621 | Reanalysis of the Chang GEO data. CIBERSORT LM22 (immune only). SFTPA2 among 340 robust DEGs. Perturbation score: PTB+HIV > PTB > non-PTB+HIV > controls (ordering partly inferred). No non-immune deconvolution | [FT] |
| Gaudio et al. | Nat Commun 2026 | 10.1038/s41467-026-74077-x / 42276999 | LLM benchmark on the same TB cohort. No tissue analysis | [Abs] |
| Wang et al. | Pathogens 2026 | 10.3390/pathogens15010055 / 41599039 | Microbial-read reanalysis of a TB cfRNA cohort (n=51). Mtb reads ≤0.001%. Host tissues not studied | [Abs] |
| Kowarsky et al. | PLoS Pathog 2026 | 10.1371/journal.ppat.1014066 / 41955174 | HIV: MHC-I antigen-presentation cfRNA signature (n=14). Immune, not a tissue signal (see A) | [Abs] |
| Gu et al. | Open Forum Infect Dis 2022 | 10.1093/ofid/ofac540 / 36519116 | Reference genes for qPCR of Mtb small RNAs in TB plasma. No host tissue signal | [Abs] |
| Lin et al. (first author not checked) | Diagn Microbiol Infect Dis 2024 | 10.1016/j.diagmicrobio.2024.116559 / 39476771 | Global 5mC of cfDNA by ELISA: TB 6.4% vs LTBI 4.1%. Not tissue-resolved | [Abs] |
| Cell Death Dis 2024 (neutrophil extracellular traps, NETs) | Cell Death Dis 2024 | 10.1038/s41419-024-06892-3 / 39085192 | Serum NET markers (NE-DNA, MPO-DNA) and cfDNA rise in TB, especially with extensive tissue damage. A neutrophil (blood) source | [Abs] |
| HCC in people with HIV, cfDNA methylation | Int J Cancer 2025; Hepatol Int 2024 | 10.1002/ijc.70204 / 41129333; 10.1007/s12072-024-10768-1 / 39704909 | Liver-disease cfDNA methylation within HIV+ cohorts. Not HIV vs healthy | [Abs] |
| **Gap** | No cfDNA methylation tissue-of-origin study (Moss-type) and no cfRNA non-immune deconvolution contrast for TB or HIV was found in PubMed | Searches as above | Uncertain: supplements, bioRxiv/medRxiv and non-PubMed sources were not searched exhaustively |

### 3. Pathophysiology-based candidate tissues (ranked)

| Rank | Tissue / process | Disease | Clinical evidence | DOI / PMID | Tag |
|---|---|---|---|---|---|
| 1 | Lung alveolar epithelium (AT2, AT1) | TB, TB+HIV | AT2 genes are already elevated in Chang cfRNA. Cavitary disease destroys parenchyma. Chang's TB score correlated with chest X-ray CAD scores (r = 0.50-0.62) | Chang 2024 | [FT] |
| 2 | Airway/BAL neutrophils (biologically blood-derived, but lung-recruited) | TB | Neutrophils make up 86.6% of sputum cells and 78.8% of BAL cells, and are the main infected phagocytes | 10.1378/chest.09-0903 / 19749004 (Eum, Chest 2010) | [Abs] |
| 3 | Granuloma macrophages / alveolar macrophages | TB | MARCO and C1QB are in Chang's DEGs. Caseating granulomas carry NETs | Chang 2024; 39085192 | [FT]/[Abs] |
| 4 | Lymph node (extrapulmonary TB) | HIV+TB > TB | Among EPTB in people with HIV: lymph node 56%, pleural 15%, disseminated 15%. Risk rises with CD4 < 200 | 10.1371/journal.pone.0232426 / 32374773 (Alemu, PLoS One 2020) | [Abs] |
| 5 | Disseminated TB, multi-organ (liver, spleen) | HIV+TB | At autopsy in adults with HIV, TB was present in 39.7% and disseminated in 87.9% of those cases; it was undiagnosed in 45.8% | 10.1097/QAD.0000000000000802 / 26266773 (Gupta, AIDS 2015) | [Abs] |
| 6 | Liver (hepatic TB; anti-TB drug hepatotoxicity) | HIV+TB | Hepatic TB in HIV: median CD4 47, ALP and GGT high, 75% had pulmonary TB. Drug-induced injury is a confounder if samples were taken during treatment | 10.1186/s12879-017-2222-2 / 28148232 (Gounder, BMC Infect Dis 2017) | [Abs] |
| 7 | Pleura (mesothelium) | TB, HIV+TB | Pleural TB is 15% of EPTB in HIV (as for rank 4) | 32374773 | [Abs] |
| 8 | Gut mucosa: GALT CD4 T-cell depletion, enterocyte damage, microbial translocation | HIV, HIV+TB | Preferential and persistent GI CD4 depletion despite ART | 10.1084/jem.20041196 / 15365095 (Mehandru, J Exp Med 2004) | [Abs] |
| 8 | (same) | (same) | LPS-indexed microbial translocation drives immune activation | 10.1038/nm1511 / 17115046 (Brenchley, Nat Med 2006) | [Abs] |
| 8 | (same) | (same) | Plasma I-FABP (FABP2, enterocyte) predicts mortality in treated HIV | 10.1093/infdis/jiu238 / 24755434 (Hunt, J Infect Dis 2014) | [Abs] |
| 9 | Endothelium | TB, HIV | Endothelium is a major baseline cfRNA source in Chang. No disease contrast reported. Evidence for TB endothelial injury was not verified here | Chang 2024 | [FT] / [recall] for injury |
| ref | Neutrophil granule (blood, kept as a reference set, not a tissue) | TB | ELANE, AZU1, CTSG, DEFA4, MPO are already Chang DEGs | Chang 2024 | [FT] |

Notes:
- HIV ART status and CD4 in Chang were not checked; the full text has no CD4 data [FT]. Gut and disseminated signals should scale with immunosuppression.
- The Chang HIV-only group contains symptomatic coughers, so lung signals there may come from non-TB pneumonia.

### 4. Final tissue-of-interest list

Gene markers are canonical HPA/PanglaoDB-style choices [recall]. **Verify each against HPA single-cell-type enrichment before use.** Lung AT2 genes are from Chang [FT].

| Tissue | Cell types (HPA / PanglaoDB names) | Canonical markers (5-10) | Disease | Rationale | Already reported in cfRNA lit? |
|---|---|---|---|---|---|
| Lung alveolar AT2 | HPA "Alveolar cells type 2" / Panglao "Pulmonary alveolar type II cells" | SFTPC, SFTPB, SFTPA1, SFTPA2, SLC34A2, NAPSA, LAMP3, ABCA3, SFTPD | TB, HIV+TB | Parenchymal destruction | **Yes**: Chang 2024 [FT] (5 genes), Wu 2026 [FT] (SFTPA2). Only vs symptomatic controls, never vs healthy |
| Lung alveolar AT1 | HPA "Alveolar cells type 1" / Panglao "Pulmonary alveolar type I cells" | AGER, RTKN2, CLIC5, EMP2, SCEL, HOPX | TB, HIV+TB | AT1 cells cover most of the alveolar surface; necrosis | No (not found) |
| Airway epithelium | HPA "Club cells", "Ciliated cells", "Basal respiratory cells" / Panglao "Club cells", "Airway epithelial cells" | SCGB1A1, SCGB3A1, SCGB3A2, BPIFA1, FOXJ1, CAPS, TPPP3, KRT5 | TB (cavity/bronchial) | Endobronchial disease; cough | No. Not checked in Chang Supp Data 1 |
| Alveolar / granuloma macrophage | HPA "Alveolar macrophages", "Macrophages" / Panglao "Alveolar macrophages" | MARCO, FABP4, PPARG, INHBA, C1QB, CHIT1, CHI3L1, APOC1 | TB | Granuloma core | **Partly**: MARCO and C1QB (Chang [FT]) |
| Lymph node / germinal center | HPA "B-cells", "Follicular dendritic cells" (check the name exists in HPA) / Panglao "Germinal center B cells" | CR2, FDCSP, CXCL13, MS4A1, AICDA, RGS13 | HIV+TB (EPTB), HIV | Most common EPTB site in HIV | No |
| Pleura | HPA "Mesothelial cells" / Panglao "Mesothelial cells" | MSLN, UPK3B, CALB2, WT1, ITLN1, PRG4 | TB, HIV+TB | Pleural TB | No |
| Liver | HPA "Hepatocytes" / Panglao "Hepatocytes" | ALB, APOA1, APOC3, FGB, AHSG, TF, CYP2E1, HPX | HIV+TB (dissemination, drug injury) | Hepatic TB, drug hepatotoxicity | No (not in TB/HIV cfRNA). Liver cfRNA is reported in PDAC (Moore 2025, see A) |
| Small intestine enterocyte | HPA "Proximal enterocytes", "Distal enterocytes" / Panglao "Enterocytes" | FABP2, APOA4, RBP2, FABP6, ALPI, SI, LCT, MTTP | HIV, HIV+TB | Gut barrier damage; I-FABP prognosis | No |
| Paneth / goblet (gut innate) | HPA "Paneth cells", "Intestinal goblet cells" / Panglao same | DEFA5, DEFA6, REG3A, MUC2, TFF3, CLCA1 | HIV | Mucosal antimicrobial loss | No |
| Endothelium | HPA "Endothelial cells" / Panglao "Endothelial cells" | VWF, PECAM1, CDH5, CLDN5, EMCN, ERG | TB, HIV | Major cfRNA source; ERG appears in Chang DEGs | **Ambiguous**: Chang lists ERG, but as a "neutrophil marker" [FT] |
| Neutrophil granule (reference, blood) | HPA "Neutrophils" / Panglao "Neutrophils" | ELANE, MPO, AZU1, CTSG, DEFA4, LTF, CAMP, S100A8 | TB | Control set for "already reported" | **Yes**: Chang 2024 [FT] |

### 5. What to verify manually

| Item | Why |
|---|---|
| Labelling of Chang TB-/HIV- samples in our HC pool | Symptomatic controls; possible lung or inflammation contamination of the reference |
| Chang Supp Data 1 (1956 DEGs): presence of AT1, airway, mesothelial, hepatocyte and enterocyte genes | Determines "already reported = yes/no" more strictly than the main text allows |
| Whether Chang samples were taken before treatment, and CD4/ART status | Drug hepatotoxicity confound; gut and dissemination scale with CD4 |
| HPA cell-type names and enrichment for each marker | Markers are [recall] |
| Earlier recall error: PMID 19815809 is a bronchial thermoplasty trial, not a TB paper. The correct Eum paper is PMID 19749004 | Recall unreliable |

---

## Tissues of interest: CAD (HF+/HF-), pre-eclampsia, ME/CFS

All papers were checked in PubMed. Tags: [FT] full text read (PMC); [Abs] abstract only; [recall] from memory, not verified; [verify] check HPA or PanglaoDB before use. "A_" means the paper is already covered in `lit/A_cfrna_paradigms.md`.

Marker genes in section 4 are canonical textbook markers [recall]. None was checked against HPA in this pass. Check each one in HPA single-cell type and Tabula Sapiens before building the gene sets.

---

### 1. CAD with HF (CAD_HF+) and without HF (CAD_HF-): Ward et al.

#### 1a. Source paper: Ward Z et al., Cells 2022, DOI 10.3390/cells11203191, PMID 36291058 (A_)

| Item | Content | Tag |
|---|---|---|
| Design | Patients from the CDCS cohort, sampled about 4 months after an acute coronary event while stable. HF+ means HF developed **within 3 years after sampling**, so most HF+ patients did not have HF when sampled. HF- patients were matched to HF+ on LVEF, E/e', NT-proBNP, hs-TnI, creatinine and medications. Controls: n=30 Canterbury healthy volunteers | [FT] |
| Method | DESeq2 with FDR<0.01 and \|FC\|>1.2. Genes kept only if ≥1 TPM in ≥90% of samples (4,150 genes). Storage time was about 5 years longer for patients than controls | [FT] |
| Results: HF+ vs HF- | **Zero DE transcripts** | [FT] |
| Results: patients vs controls, blood/immune | No blood or immune signal was reported explicitly | [FT] |
| Results: patients vs controls, non-blood | 170 DE transcripts. 13 of the top 20 are mitochondrial-genome mRNAs, mostly higher in patients and correlated with NT-proBNP and hsTnI. The authors read this as "myocardial ischaemia and oxidative stress". COA6 and NDUFAF2 were lower in patients; the paper calls them mitochondrial, but they are nuclear-encoded [recall]. FGF23 (cardiac hypertrophy) and STAG2 were also named. No tissue deconvolution or enrichment was done | [FT] |
| Supplementary | The full DE list (170 transcripts) is in the MDPI supplement ("Supplementary Materials", Table S-numbered). The exact table name was not captured from PMC text | [verify] |

#### 1b. Other cfDNA/cfRNA tissue-of-origin studies (cardiac)

| Study | DOI / PMID | Tissue signal | Relevance to stable CAD | Tag |
|---|---|---|---|---|
| Zemmour, Nat Commun 2018 | 10.1038/s41467-018-03961-y / 29691397 | Cardiomyocyte cfDNA (methylation): essentially zero in healthy plasma, high in STEMI | Acute only. Stable CAD is expected to be near zero | [Abs] |
| Moss, Nat Commun 2018 | 10.1038/s41467-018-07466-6 / 30498206 | Healthy cfDNA comes from WBC 55%, erythroid progenitors 30%, vascular endothelium 10%, hepatocytes 1% | Baseline that includes endothelium | [Abs] |
| Cuadrat, NAR Genom Bioinform 2023 | 10.1093/nargab/lqad061 / 37388821 | cfDNA methylation deconvolution separates ACS types | Acute | [Abs] |
| Dash, Clin Epigenetics 2024 | 10.1186/s13148-024-01755-2 / 39415189 | AMI cfDNA: about 75% neutrophil and about 10% "left atrium" | Acute. Small n=25, lower-IF journal | [Abs] |
| LVAD cfDNA, Circ Heart Fail 2025 (first author not recorded) | 10.1161/CIRCHEARTFAILURE.124.013302 / 41307146 | Advanced HF (LVAD): nuclear cfDNA 4x healthy. After LVAD, innate immune, **vascular endothelium, GI and liver** cfDNA fell | Supports end-organ (liver, gut, endothelium) shedding in advanced HF only | [Abs] |
| Nidadavolu, J Gerontol A 2025 | 10.1093/gerona/glaf081 / 40235369 | A "cardiovascular-enriched" cfDNA cluster in older adults was associated with MI, CHD and frailty | Chronic setting, but array-based | [Abs] |
| Loy, Sci Rep 2025 (Kawasaki disease) | 10.1038/s41598-025-15843-7 / 40813796 | cfRNA: hepatocyte injury and platelet differences | Not CAD. Only shows that cfRNA can detect hepatocyte signals | [Abs] |
| cfRNA in stable CAD/HF besides Ward | none found | Only Ward 2022 appeared in a PubMed search ("cell-free RNA" AND HF/CAD) | Gap | search |

#### 1c. Pathophysiology-based candidate tissues (ranked)

| Rank | Tissue / cell | Rationale | Clinical citation | Tag |
|---|---|---|---|---|
| 1 | Cardiomyocyte | Low-level chronic myocyte injury in stable CAD: hs-cTnT is detectable in 97.7% and predicts HF and CV death | Omland, NEJM 2009, 10.1056/NEJMoa0805299 / 19940289 | [Abs] |
| 2 | Vascular endothelium | Atherosclerosis; endothelium is about 10% of healthy cfDNA (Moss 2018) | Moss 2018 | [Abs] |
| 3 | Kidney (tubule/podocyte) | Cardiorenal syndrome | Rangaswami, Circulation 2019, 10.1161/CIR.0000000000000664 / 30852913 | [Abs] |
| 4 | Liver (hepatocyte) | Congestive hepatopathy and cardiogenic liver injury in chronic HF | Sessa, JHEP Rep 2021, 10.1016/j.jhepr.2021.100249 / 33665589; Brankovic, J Clin Transl Hepatol 2023, 10.14218/JCTH.2022.00294 / 37408802 | [Abs] |
| 5 | Cardiac fibroblast | Post-MI remodelling. FGF23 is in Ward's DE list | Ward 2022 | [FT] |
| Caveat | HF+ and HF- were sampled before HF onset and matched on cardiac severity. Expect **no** congestion (liver/kidney) signal and no HF+ vs HF- difference | Ward 2022 | [FT] |

---

### 2. Pre-eclampsia (PE): Moufarrej et al.

#### 2a. Source paper: Moufarrej MN et al., Nature 2022, DOI 10.1038/s41586-022-04410-z, PMID 35140405 (A_)

| Item | Content | Tag |
|---|---|---|
| Method | DE was tested for enrichment (one-sided hypergeometric) against HPA, "HPA strict" and TSP+. TSP+ is Tabula Sapiens plus per-organ atlases for placenta, liver, kidney, heart and brain. Per-sample scores were the sum of logCPM | [FT] |
| Results: blood/immune | Enrichment in the adaptive immune system (bone marrow) among increased genes. Innate cells (granulocytes, neutrophils, thymus) among genes that switch sign over gestation. Increased platelet signal before 12 weeks. 18-gene panel: 2 immune genes | [FT] |
| Results: non-blood | **Neuromuscular**: excitatory neurons, oligodendrocytes, smooth muscle, plus "muscle contraction" and "dilated cardiomyopathy" terms. **Endothelial cells**. Signature scores: **placenta and syncytiotrophoblast decreased** before 20 weeks; **hepatocyte, kidney, endothelial and smooth muscle decreased** across gestation. Authors' summary: "dysfunction in at least five organ systems (brain, liver, kidney, muscle and bone marrow)". 18-gene panel: placenta and neuromuscular genes per HPA (symbols stripped in PMC text) | [FT] |
| Inconsistency | Smooth muscle goes up in the enrichment result but down in the signature score | [FT] |
| Supplementary | 18 genes: Extended Data Table. DEGs: Supplementary Tables (exact numbers not captured) | [verify] |
| Our-data caveat | Our healthy pool contains about 194 normotensive pregnant samples from this same study. Placental (STB/EVT) genes are therefore partly "normal" in the reference. Unless pregnancy and gestational age are covariates, the placental Z of a PE sample mixes the pregnancy effect with the PE effect. Moufarrej reports placental signal **decreased** in PE, so it will be hard to separate | project memory |

#### 2b. Other cfRNA/cfDNA PE studies with tissue signals

| Study | DOI / PMID | Tissue signal | Tag |
|---|---|---|---|
| Munchel, Sci Transl Med 2020 (A_) | 10.1126/scitranslmed.aaz0131 / 32611681 | 30 DE transcripts at early-onset PE diagnosis, consistent with "placental dysfunction, fetal development, maternal immune and cardiovascular" changes. Gene list not read | [Abs] |
| Rasmussen, Nature 2022 (A_) | 10.1038/s41586-021-04249-w / 34987224 | 7-gene PE signature that includes PAPPA2 (trophoblast), CLDN7 and FABP1 (cytotrophoblast-to-STB, fetal liver); the other symbols were stripped. Pathways up: placental blood vessel development and artery morphogenesis. Pathways down: immune | [FT] |
| Tsang, PNAS 2017 | 10.1073/pnas.1710470114 / 28830992 | Placental scRNA signatures mapped to cfRNA: **EVT dysfunction** in early PE | [Abs] |
| Del Vecchio, Epigenetics 2021 | 10.1080/15592294.2020.1816774 / 33045922 | Placental cfDNA fraction did not change in first-trimester PE. cfRNA: S100A8, MS4A3 and MMP8 up (neutrophil), ALPL down | [Abs] |
| MeD-seq EOPE, Placenta 2025 | 10.1016/j.placenta.2025.12.006 / 41477915 | After diagnosis: placenta, liver and kidney cfDNA fractions decreased **relatively**; **endothelium increased**. Total cfDNA 20x higher | [Abs] |
| cfDNA fragmentomics, Nat Commun 2026 | 10.1038/s41467-026-72682-4 / 42069757 | TSS coverage differs for blood-cell and placenta genes; used for risk prediction | [Abs] |
| AIF1L, Transl Res 2026 | 10.1016/j.trsl.2026.07.014 / 42526532 | cfRNA 10-gene panel (AIF1L, FLT1, SH3PXD2A). scRNA places AIF1L in **EVT** | [Abs] |
| GW13 reanalysis, Placenta 2026 | 10.1016/j.placenta.2026.09.008 / 42762647 | Reanalysis of Moufarrej data: MYZAP, CAV2, OXTR, GBP1 | [Abs] |
| J Reprod Immunol 2026 | 10.1016/j.jri.2026.104952 / 42623724 | KRT19, CCND1, ZBTB16, TUBB6, HMGB2. Lower-IF journal | [Abs] |

#### 2c. Pathophysiology-based candidate tissues (ranked)

| Rank | Tissue / cell | Rationale | Citation | Tag |
|---|---|---|---|---|
| 1 | Placenta: syncytiotrophoblast | STB stress releases soluble factors (sFlt1) that drive the maternal syndrome | Burton, BMJ 2019, 10.1136/bmj.l2381 / 31307997; Maynard, JCI 2003, 10.1172/JCI17189 / 12618519 | [Abs] |
| 2 | Placenta: EVT | Defective placentation in early-onset PE | Burton 2019; Tsang 2017 | [Abs] |
| 3 | Systemic endothelium | sFlt1 causes endothelial dysfunction | Maynard 2003 | [Abs] |
| 4 | Kidney: glomerular endothelium/podocyte | Proteinuria; glomerular endotheliosis in the sFlt1 rat | Maynard 2003 | [Abs] |
| 5 | Liver: hepatocyte | Liver involvement is part of the PE definition (HELLP) | Moufarrej 2022 intro | [FT] |
| 6 | Brain: neurons, oligodendrocytes | Brain involvement is part of the PE definition (eclampsia); neuronal enrichment in Moufarrej | Moufarrej 2022 | [FT] |
| 7 | Smooth muscle | Reported by Moufarrej, with conflicting direction | Moufarrej 2022 | [FT] |

---

### 3. ME/CFS: Gardella et al.

#### 3a. Source paper: Gardella AE et al., PNAS 2025, DOI 10.1073/pnas.2507345122, PMID 40789036 (A_)

| Item | Content | Tag |
|---|---|---|
| Method | DESeq2 (padj<0.05, \|log2FC\|>0.5): 743 features. BayesPrism with the Vu 2024 ME/CFS PBMC atlas. A Tabula Sapiens deconvolution was also run. IPA pathways. Kraken2 for viral reads | [FT] |
| Results: blood/immune | Up in cases: pDC, monocytes, naive CD8 T cells, MAIT cells, T-cell clusters. Platelets down. Up genes: IL18R1, CCR7, FCRL5, PRF1, IFNLR1, TCR genes, HNRNPLL. Pathways: T-cell exhaustion, JAK-STAT | [FT] |
| Results: non-blood | **Tabula Sapiens deconvolution: solid-tissue cfRNA was "minimal" with no case-control difference.** Higher in **controls**: COL1A1, TGFB2, CD109, CCL3, ADAM12, MMP2. These are stromal/ECM genes, interpreted by the authors as "tissue repair". mtRNA fraction lower in cases (16.4% vs 19.4%); MTG1 up; several MT-tRNAs down | [FT] |
| Supplementary | Dataset S01 (XLSX) and Appendix 01 (PDF) on the PNAS site. Dataset S01 most likely holds the DA list | [verify] |

#### 3b. Other circulating tissue-of-origin studies in ME/CFS

| Study | DOI / PMID | Finding | Tag |
|---|---|---|---|
| None found (cfDNA or cfRNA tissue deconvolution) | PubMed search, 29 hits | Hits are miRNA/EV or PBMC studies only. Gardella is the only whole-transcriptome cfRNA study | search |
| EV miRNA, J Transl Med 2026 | 10.1186/s12967-026-08695-w / 42533331 | EV miRNA classifier; no tissue-of-origin analysis per title | title only |

#### 3c. Pathophysiology-based candidate tissues (ranked)

| Rank | Tissue / cell | Rationale | Citation | Tag |
|---|---|---|---|---|
| 1 | Skeletal muscle (myocyte) | Mitochondrial damage and necrosis/regeneration in muscle biopsies; WASF3 up and ER stress in ME/CFS muscle | Scheibenbogen, J Cachexia Sarcopenia Muscle 2025 (review), 10.1002/jcsm.13669 / 39727052; Wang, PNAS 2023, 10.1073/pnas.2302738120 / 37579159 | [Abs] |
| 2 | Brain/CNS | Dysfunction of integrative brain regions and catechol pathway changes | Walitt, Nat Commun 2024, 10.1038/s41467-024-45107-3 / 38383456 | [Abs] |
| 3 | Gut epithelium | Dysbiosis and reduced butyrate in short-duration disease. This is evidence of microbial dysbiosis, not epithelial shedding | Xiong, Cell Host Microbe 2023, 10.1016/j.chom.2023.01.001 / 36758521 | [Abs] |
| 4 | Fibroblast/stroma | COL1A1, MMP2, ADAM12 and TGFB2 lower in cases (direction: down) | Gardella 2025 | [FT] |
| Caveat | The only cfRNA study found no solid-tissue difference. Expected tissue signal is weak, and sampling was at rest (no post-exertion samples) | Gardella 2025 | [FT] |

---

### 4. Final tissue-of-interest list

HPA cell-type names and PanglaoDB terms are [recall] [verify]. Marker genes are [recall].

| Disease | Tissue | Cell types (HPA / PanglaoDB) | Markers (5-10) | Rationale | Already in cfRNA literature? |
|---|---|---|---|---|---|
| CAD_HF± | Heart, myocyte | Cardiomyocytes / Cardiomyocytes | MYH6, MYH7, TNNT2, TNNI3, MYL7, NPPA, NPPB, ACTC1, MYBPC3 | Chronic hs-TnT release (Omland 2009) | **No** for cfRNA (Ward found only mito-genome transcripts). Yes for acute cfDNA (Zemmour 2018) |
| CAD_HF± | Heart, fibroblast | Fibroblasts / Fibroblasts (cardiac subtype not separable) | POSTN, DCN, COL1A1, COL3A1, FAP, PDGFRA | Post-MI remodelling. Low specificity | Partial: FGF23 named by Ward |
| CAD_HF± | Endothelium | Endothelial cells / Endothelial cells | PECAM1, VWF, CDH5, CLDN5, EMCN, NOS3 | Atherosclerosis | **No** for CAD cfRNA. Yes for HF cfDNA (LVAD 2025) |
| CAD_HF+ | Liver | Hepatocytes / Hepatocytes | ALB, APOA1, APOC3, FGB, HP, SERPINA1, TF, CYP2E1 | Congestion. Weak expectation given sampling before HF onset | **No** for cfRNA. Yes for cfDNA in advanced HF (LVAD 2025) |
| CAD_HF+ | Kidney | Proximal tubular cells; Podocytes? / Proximal tubule cells, Podocytes | LRP2, CUBN, SLC34A1, UMOD, NPHS1, NPHS2, PODXL | Cardiorenal | **No** |
| PE | Placenta, STB | Syncytiotrophoblasts / Trophoblast cells | CGA, CGB3, PSG1, PSG4, KISS1, CYP19A1, ERVW-1, GDF15, FLT1, LGALS13 | Core of the disease. Confounded by pregnant controls in our pool | **Yes** (Moufarrej 2022; Munchel 2020) |
| PE | Placenta, EVT | Extravillous trophoblasts / Trophoblast cells | HLA-G, PAPPA2, ITGA5, MMP2, ASCL2, AIF1L | Early-onset PE | **Yes** (Tsang 2017; AIF1L 2026; PAPPA2 in Rasmussen 2022) |
| PE | Endothelium | Endothelial cells | PECAM1, VWF, CDH5, CLDN5, EMCN, ENG | Maternal syndrome | **Yes** (Moufarrej 2022; cfDNA 2025) |
| PE | Kidney | Podocytes; Proximal tubular cells | NPHS1, NPHS2, PODXL, PTPRO, LRP2, UMOD | Proteinuria | **Yes**, as a decreased kidney score (Moufarrej 2022) |
| PE | Liver | Hepatocytes | ALB, APOA2, APOH, FGB, HP, SERPINA1, AHSG | HELLP | **Yes**, decreased (Moufarrej 2022) |
| PE | Brain | Excitatory neurons; Oligodendrocytes | SNAP25, SLC17A7, NRGN, MBP, PLP1, MOG | Eclampsia | **Yes** (Moufarrej 2022) |
| PE | Smooth muscle | Smooth muscle cells | ACTA2, MYH11, CNN1, DES, TAGLN | Vascular tone | **Yes**, direction conflicting (Moufarrej 2022) |
| ME/CFS | Skeletal muscle | Skeletal myocytes / Myocytes? | ACTA1, MYH1, MYH2, MYH7, CKM, MB, TNNT1, TNNT3, DES | PEM, muscle mitochondrial damage | **No** (Gardella TS deconvolution found no difference) |
| ME/CFS | Brain/CNS | Excitatory neurons; Astrocytes; Oligodendrocytes | SNAP25, SYT1, GFAP, AQP4, MBP, PLP1 | Walitt 2024 | **No** |
| ME/CFS | Gut epithelium | Distal/Proximal enterocytes? / Enterocytes | FABP2, APOA4, VIL1, CDX2, LGALS4, MUC2 | Dysbiosis. Indirect evidence only | **No** |
| ME/CFS | Fibroblast/stroma | Fibroblasts | COL1A1, COL3A1, DCN, LUM, MMP2, ADAM12 | Lower in cases | **Yes**, as individual DE genes; not framed as a tissue (Gardella 2025) |

#### Key uncertainties

| Evidence source | Evidence strength | Main uncertainty | What to verify manually |
|---|---|---|---|
| Ward 2022 [FT] | Moderate: n=89, 1 batch | HF+ was sampled before HF; storage time is confounded with group | Supplement DE list; whether any cardiac genes pass the TPM filter |
| Moufarrej 2022 [FT] | High (Nature) | Tissue directions conflict; HC pool includes pregnant women | Extended Data 18 genes; per-sample scores |
| Gardella 2025 [FT] | High (PNAS) | Little solid-tissue signal, so a no-tissue-signal result is likely | Dataset S01; TS deconvolution SI figure |
| Marker lists [recall] | Low until checked | HPA or Tabula Sapiens specificity; whether markers are detectable in plasma | HPA single-cell type enrichment; detection rate in our HC pool (see RQR zero-count note) |

---

