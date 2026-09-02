# Bulk RNA-seq Dataset Summary

**Project:** SCA2 / SCA7 / APOE4 bulk RNA-seq analysis  
**Snapshot date:** 2026-09-02  
**Primary sources:** `config/datasets.tsv`, `metadata/analysis/*.tsv`, and `config/contrasts.tsv`

## Scope

The analysis includes 11 accessions represented by 13 analysis units, 300
biological samples, and 33 registered differential contrasts. GSE102956 is
split by sequencing layout; GSE163857 is split by species. The dataset registry
also documents seven excluded accessions, all excluded because they are single
cell or single nucleus experiments.

## Included Analysis Units

| Analysis unit | Accession | Focus | Species | Samples | Tissue / cell type | Time point / age | Sex | Groups | Registered contrasts |
| --- | --- | --- | --- | ---: | --- | --- | --- | ---: | ---: |
| E-MTAB-6293 | E-MTAB-6293 | ATXN2 | Mouse | 64 | Cerebellum | P1, 3 wk, 6 wk | Male | 6 | 3 |
| GSE138527 | GSE138527 | ATXN7 | Mouse | 6 | Cerebellum | 40 wk | Not recorded | 2 | 1 |
| GSE139090 | GSE139090 | ATXN7 | Mouse | 12 | Not recorded | 12 wk, 29 wk | Not recorded | 4 | 2 |
| GSE271392 | GSE271392 | ATXN7 | Mouse | 35 | Medulla; cervical spinal cord | P8, 5 wk, 9 wk | Not recorded | 12 | 5 |
| GSE102956_paired | GSE102956 | APOE | Human | 8 | Microglia-like cells | Not recorded | Not recorded | 2 | 1 |
| GSE102956_single | GSE102956 | APOE | Human | 18 | Neurons; astrocytes; iPSCs | Not recorded | Not recorded | 6 | 3 |
| GSE140205 | GSE140205 | APOE | Mouse | 20 | Hippocampus | 16 months | Female, male | 4 | 2 |
| GSE163857_human | GSE163857 | APOE | Human | 24 | iPSC-derived microglia | Not recorded | Female, male | 12 | 3 |
| GSE163857_mouse | GSE163857 | APOE | Mouse | 30 | MACS-sorted CD11b+ brain cells | Not recorded | Female, male | 8 | 3 |
| GSE234350 | GSE234350 | APOE | Mouse | 33 | Brain; phagocytic and non-phagocytic microglia | 8 months | Female, male | 8 | 2 |
| GSE234480 | GSE234480 | APOE | Mouse | 8 | Astrocytes | Not recorded | Male | 2 | 1 |
| GSE234488 | GSE234488 | APOE | Human | 25 | Brain whole tissue | 67-91 years | Female, male | 2 | 1 |
| GSE254205 | GSE254205 | APOE | Human | 17 | iMicroglia | Not recorded | Not recorded | 7 | 6 |

## Registered Contrasts

The reference group is the first group and the alternative group is the second
group. These names are the exact labels passed to DESeq2 and Shiba.

| Analysis unit | Contrast | Reference group | Alternative group |
| --- | --- | --- | --- |
| E-MTAB-6293 | SCA2_vs_WT_P1 | wild_type_genotype__1_day | ATXN2Q127__1_day |
| E-MTAB-6293 | SCA2_vs_WT_3wk | wild_type_genotype__3_week | ATXN2Q127__3_week |
| E-MTAB-6293 | SCA2_vs_WT_6wk | wild_type_genotype__6_week | ATXN2Q127__6_week |
| GSE138527 | SCA7_vs_WT_40wk | wildtype | SCA7_140Q_5Q_knock_in_Ki |
| GSE139090 | SCA7_vs_WT_12wk | wild_type__12_weeks | fxSCA7_92Q__12_weeks |
| GSE139090 | SCA7_vs_WT_29wk | wild_type__29_weeks | fxSCA7_92Q__29_weeks |
| GSE271392 | SCA7_vs_WT_medulla_P8 | Wildtype__Medulla__p8 | SCA7_266Q_5Q__Medulla__p8 |
| GSE271392 | SCA7_vs_WT_medulla_5wk | Wildtype__Medulla__5_weeks | SCA7_266Q_5Q__Medulla__5_weeks |
| GSE271392 | SCA7_vs_WT_medulla_9wk | Wildtype__Medulla__9_weeks | SCA7_266Q_5Q__Medulla__9_weeks |
| GSE271392 | SCA7_vs_WT_cervical_5wk | Wildtype__Cervical_spinal_cord__5_weeks | SCA7_266Q_5Q__Cervical_spinal_cord__5_weeks |
| GSE271392 | SCA7_vs_WT_cervical_9wk | Wildtype__Cervical_spinal_cord__9_weeks | SCA7_266Q_5Q__Cervical_spinal_cord__9_weeks |
| GSE102956_paired | APOE4_vs_APOE3_microglia | APOE3__microglia_like_cells | APOE4__microglia_like_cells |
| GSE102956_single | APOE4_vs_APOE3_astrocytes | APOE3__astrocytes | APOE4__astrocytes |
| GSE102956_single | APOE4_vs_APOE3_iPSCs | APOE3__iPSCs | APOE4__iPSCs |
| GSE102956_single | APOE4_vs_APOE3_neurons | APOE3__neuron | APOE4__neuron |
| GSE140205 | APOE4_vs_APOE3_female | APOE3__Female | APOE4__Female |
| GSE140205 | APOE4_vs_APOE3_male | APOE3__Male | APOE4__Male |
| GSE163857_human | E3E4_vs_E3E3_CTL_female | E3_E3__Healthy_Control__CTL__Female | E3_E4__Healthy_Control__CTL__Female |
| GSE163857_human | E4E4_vs_E3E3_CTL_female | E3_E3__Healthy_Control__CTL__Female | E4_E4__Healthy_Control__CTL__Female |
| GSE163857_human | E3E4_vs_E3E3_CTL_male | E3_E3__Healthy_Control__CTL__Male | E3_E4__Healthy_Control__CTL__Male |
| GSE163857_mouse | E4E4_vs_E3E3_FAD_female | E3_E3__FAD_transgenic__Female | E4_E4__FAD_transgenic__Female |
| GSE163857_mouse | E4E4_vs_E3E3_FAD_male | E3_E3__FAD_transgenic__Male | E4_E4__FAD_transgenic__Male |
| GSE163857_mouse | E4E4_vs_E3E3_TR_female | E3_E3__Targeted_Replacement_Control__Female | E4_E4__Targeted_Replacement_Control__Female |
| GSE234350 | APOE4_vs_APOE3_KI_phagocytic | APOE3_KI__Phagocytic_microglia | APOE4_KI__Phagocytic_microglia |
| GSE234350 | APOE4_vs_APOE3_KI_nonphagocytic | APOE3_KI__Non_phagocytic_microglia | APOE4_KI__Non_phagocytic_microglia |
| GSE234480 | APOE4_cKO_vs_APOE4 | WT_APOE4 | WT_APOE4_cKO |
| GSE234488 | APOEe3e4_vs_APOEe3e3 | AD_APOEe3_e3 | AD_APOEe3_e4 |
| GSE254205 | fAB_vs_untreated_APOE44 | APOE_4_4__untreated | APOE_4_4__fAB_treated |
| GSE254205 | GNE_fAB_vs_fAB_APOE44 | APOE_4_4__fAB_treated | APOE_4_4__GNE_317_fAB_treated |
| GSE254205 | LDhigh_vs_LDlow_APOE33 | APOE_3_3__APOE_3_3_iMG_LD_low__LD_low | APOE_3_3__APOE_3_3_iMG_LD_high__LD_high |
| GSE254205 | LDhigh_vs_LDlow_APOE44 | APOE_4_4__APOE_4_4_iMG_LD_low__LD_low | APOE_4_4__APOE_4_4_iMG_LD_high__LD_high |
| GSE254205 | APOE44_vs_APOE33_LDlow | APOE_3_3__APOE_3_3_iMG_LD_low__LD_low | APOE_4_4__APOE_4_4_iMG_LD_low__LD_low |
| GSE254205 | APOE44_vs_APOE33_LDhigh | APOE_3_3__APOE_3_3_iMG_LD_high__LD_high | APOE_4_4__APOE_4_4_iMG_LD_high__LD_high |

## Excluded Accessions

| Accession | Focus | Exclusion reason |
| --- | --- | --- |
| GSE164507 | APOE | Single-nucleus RNA-seq |
| GSE213446 | APOE | Single-nucleus RNA-seq |
| GSE223719 | APOE | Single-cell RNA-seq |
| GSE237718 | APOE | Single-nucleus RNA-seq |
| GSE242153 | APOE | Single-nucleus RNA-seq |
| GSE248020 | APOE | Single-cell RNA-seq |
| GSE295612 | APOE | Single-cell RNA-seq |

## Notes

- Only bulk RNA-seq datasets are in scope. GSE254205 retains only the 17 bulk
  iPSC-derived microglia samples from its mixed series.
- GSE271392 medulla and cervical spinal cord samples are analyzed separately;
  the cervical spinal cord P8 comparison is not registered because it has one
  wild-type sample.
- Job identifiers in `jobs/` document submissions and dependencies. They are
  not a substitute for checking current scheduler state or result completion.