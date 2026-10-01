# Identifying and Intercepting Pathologic Pore-Forming Proteins in Ventricular Tachycardia
This repository contains the documentation for the project Engineering Immunotherapy to Prevent Ventricular Arrhythmia.
## General Information
* PI's Name: Noor Momin
* PI's ORCID: https://orcid.org/0000-0003-3145-6089 
* Data collected Jan 2024- March 2026
* Proposal Number (pC ID): 1269197
* Award Number: 24IPA1269197
  
## PORCUPINE 
We developed a computational pipeline, <ins>POR</ins>e-forming proteins <ins>C</ins>ompUtational <ins>P</ins>rediction pipl<ins>INE</ins> (PORCUPINE), combining transcriptomics, AlphaFold3 protein structure prediction, and coarse-grained molecular dynamics simulations. 
<img width="711" height="107" alt="Screenshot 2026-10-01 at 3 09 23 PM" src="https://github.com/user-attachments/assets/bd684704-f186-4f03-bf14-daa1b752d5f3" />
### Transcriptomic Analysis 
We used [R code file](Jafarah-et-al-analysis/Psuedobulk%20and%20DEG.R) to re-analyze data deposited in GEO.
<img width="570" height="355" alt="Screenshot 2026-10-01 at 3 17 50 PM" src="https://github.com/user-attachments/assets/4c2fdbab-9a9d-4980-96e4-c4cd8a11c0f9" />

### Molecular Dynamics Simulation
Top hits from the transcriptomic analysis was shuffled into protein structure prediction and molecular dynamic simulation. The [README_insane](Jafarah-et-al-analysis/README_insane) was used as the input parameter file for GROMACS to run MARTINI based simulations. 
