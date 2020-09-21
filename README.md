# Coral microbiome transplants

This repository contains the scripts used to analyze data and create figures for the manuscript “Surface topography, bacterial carrying capacity, and the prospect of microbiome transplants in the sea anemone coral model Aiptasia”.  

## Summary
the development of coral model systems is imperative. In this paper we set the foundation for Aiptasia as a model for the study of coral-bacterial interactions by 1) compared ultrastructural surface architecture of Aiptasia and stony corals, 2) determined Aiptasia bacterial carrying capacity and 3) established a protocol for bacterial depletion of Aiptasia polyps to facilitate microbial transplants from two species of stony corals

## Workflow for 16S rDNA amplicon data analysis

1. Mothur pipeline is described in `Transplants_mothur.sh`
2. The script `Contaminant_removal.R` identifies and remove bacterial contamination from negative controls and exports a new OTU table 
3. The script `Transplants_ordination.R` plots PCoAs for each symbiotic state colored by treatment.
4. The script `Transplants_VennDiagram.R` plots Venn diagrams to compare overlapping OTUs between treatments.
5. The script `Transplants_adonis.R` tests bacterial community composition between treatments and symbiotic states.
6. The script `Transplants_TaxaComposition.R` makes bar plots to visualize most abundant OTUs in eahc treatment.

