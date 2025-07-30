# 🌿 carbon-optimisation-png

A MATLAB-based modelling project to identify and prioritise forest conservation sites in Papua New Guinea using static and dynamic optimisation for carbon credit investment.

## 📌 Project Summary

This project was developed to support The Nature Conservancy (TNC) in identifying and prioritising Intact Forest Landscapes (IFLs) in Papua New Guinea for conservation investment through carbon credit generation. The goal is to balance ecological, economic, and social considerations to determine the most viable locations for resource allocation using both static and dynamic optimisation models.

## 🧠 Key Features

- **Real-world Data Integration**: Utilises datasets from NASA, IUCN, and The Nature Conservancy.
- **Robust Prioritisation Framework**: Considers cost, risk, ecological value, and feasibility.
- **Dynamic Resource Allocation**: Solves time-dependent conservation investments using ODEs and Pontryagin’s Maximum Principle.
- **Sensitivity and Uncertainty Analysis**: Validates model reliability under data variability.

## 🗂️ MATLAB File Descriptions

| File Name         | Description |
|------------------|-------------|
| `carbon reader.m` | Extracts and processes carbon density data from NASA's Earth Data Carbon Map. Calculates average carbon density for each IFL. |
| `data collection.m` | Integrates and cleans raw datasets, maps them onto IFL polygons, and aligns resolutions for model input. |
| `static model.m` | Implements the static prioritisation model. Computes expected cost-effectiveness and performs both full and partial allocation using Lagrange multipliers. |
| `dynamic model.m` | Solves the system of coupled ODEs that represent conservation impact over time. Uses Pontryagin’s Maximum Principle to determine optimal investment strategies. |
| `sort.m` | Utility script for ranking IFLs and visualising results. |
| `filter clip.m` | Pre-processes spatial data, filtering and clipping geographic datasets to focus on the PNG region. |

## 📈 Results Overview

- SEA_9 and SEA_21 were among the top-ranked IFLs based on their ecological value and accessibility.
- Partial allocation models demonstrated higher biodiversity returns by distributing investment across several IFLs.
- Dynamic modelling revealed how conservation strategies should shift over time as conditions and land availability change.

## 🧪 Limitations

- **Outdated Data**: Carbon and roads datasets from 2010 or earlier may not reflect current conditions.
- **Simplified Logistics**: Travel times and accessibility approximated using static geodesic distances and assumptions.
- **Single IFL Allocation (Dynamic Model)**: Current model assumes only one IFL is funded at a time. Real-world conservation may involve simultaneous investment in multiple sites.

## 📚 Full Report

Please see the [Identification of Areas for Investment in PNG IFLs.pdf](./Identification%20of%20Areas%20for%20Investment%20in%20PNG%20IFLs.pdf) for complete methodology, results, discussion, and model derivations.

## 🤝 Acknowledgements

This project was completed as part of a capstone unit at Queensland University of Technology.  
Special thanks to The Nature Conservancy (TNC) for their support and guidance.

---

*Created by Zane Bowser | Faculty of Science, QUT | 2024*
