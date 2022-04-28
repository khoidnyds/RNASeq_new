# RNA-Seq analysis pipeline

## Data
Colorectal cancer (CRC) patients frequently experience disease recurrence and distant metastasis. This study aimed to identify prognostic indicators, including individual responses to chemotherapy, in CRC patients. RNA-seq data was generated using 54 samples (normal colon, primary CRC, and liver metastases) from 18 CRC patients and genes associated with CRC aggressiveness were identified. A risk score based on these genes was developed and validated in four independent CRC patient cohorts (n = 1063). Diverse statistical methods were applied to validate the risk scoring system, including a generalized linear model likelihood ratio test, Kaplan-Meier curves, a log-rank test, and the Cox model. TREM1 and CTGF were identified as two activated regulators associated with CRC aggressiveness. A risk score based on 19 genes regulated by TREM1 or CTGF activation (TCA19) was a significant prognostic indicator. In multivariate and subset analyses based on pathological staging, TCA19 was an independent risk factor (HR = 1.894, 95% CI = 1.227-2.809, P = 0.002). Subset stratification in stage III patients revealed that TCA19 had prognostic potential and identified patients who would benefit from adjuvant chemotherapy, regardless of age. The TCA19 predictor represents a novel diagnostic tool for identifying high-risk CRC patients and possibly predicting the response to adjuvant chemotherapy. 

> Ref: https://pubmed.ncbi.nlm.nih.gov/25049118/


> RNAseq data: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA218851

Instrument: Illumina HiSeq 2500 - RNA-Seq

Reference genomes: Human GRCh38.p14

> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/

## Set-up
- data/
    SRR_Acc_List.txt
- run.sh
- deseq2.r
## Tools: 
sratoolskit, fastqc, multiqc, catadapt, star, featureCounts, deseq2

## Usage:
```
chmod +x run.sh

./run.sh
```