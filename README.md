This repository consists of all the analysis performed for the findings in manuscript: "A distinct monocyte cellular state links systemic immune dysregulation to pulmonary impairment in long COVID".

All the analysis was performed using Seurat version 5 and Signac version 1.13. 
The scripts are organized in different folders, based on the analysis performed:
- RobjCreation_integration includes codes used to create Seurat and Signac object and integration of datasets.
- MajorCellTypes_dgeAnalysis includes Rscripts used to perform differential gene expression analysis for each major cell types and across different categories.
- pseudobulk_pathwayAndCorrelations includes R scripts for each major cell type, GSEA analysis and calcualting AUC score on pseudobulk of samples to correlate with fatigue scores.
- NeighbourhoodAnalysis includes neighbourhood analysis done using MiloR for CD14 monoctyes, CD8 T cells and NK cells. Also includes code to calculate correlation with FAS and pO2.
- ValidationCohort3_analysis includes analysis done for cohort 3 validation cohort, starting from data integration, to differences within monocytes and enrichment of MC4 state.
- OpenChromatinAnalysis includes all the code for peaks analysis using Signac and Chromvar.
- cohort5_validation_stimulation_analysis includes validation of MC4 in public dataset from PMID: 39018367 and stimulation dataset generated in-house.

Contact:
Saumya Kumar Center for Individualized Infection Medicine. email: saumya.dileepkumar@helmholtz-hzi.de.
