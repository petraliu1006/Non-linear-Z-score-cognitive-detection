# Nonlinear Z-score for Abnormal Cognitive Detection

## Overview
This project implements nonlinear Z-score calculations for detecting abnormal cognitive performance, with a focus on demographic adjustments for age, education, sex, and race.

## Key Features
- Box-Cox transformation for normalization
- Nonlinear age-adjusted means and standard deviations
- Race-specific cognitive trajectories
- Automated Z-score lookup table generation

## Data Sources
- NACC (National Alzheimer's Coordinating Center) dataset
- Trail Making Test normative data

## Methods
1. Data preprocessing and filtering
2. Box-Cox transformation
3. Nonlinear modeling using SCAM (Shape-Constrained Additive Models)
4. Z-score calculation with demographic adjustments

## Installation
```bash
git clone https://github.com/username/nonlinear-zscore-cognitive
cd nonlinear-zscore-cognitive
Rscript requirements.txt
