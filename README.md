# XAI_boron_LA

## Project Description
This repository provides resources and tools for applying Explainable AI (XAI) techniques to the study and design of Boron Lewis acids. It accompanies the article:
"Constructing and explaining machine learning models for chemistry: example of the exploration and design of boron derivatives with targeted Lewis acidity."

This project demonstrates the use of machine learning models to predict and interprete Fluoride Ion Affinities (FIA) of various molecular scaffolds, showcasing how explainability can aid in understanding and guiding molecular design.

#### Reference:

## Installation
After cloning the repository, install the necessary dependencies running in the requirements

## Repository Structure
- scripts/: Contains the main scripts used in the project to perform the analyses and optimize ML models. The optimized ML models can be found in the scripts/models.py file.
- data/: Includes datasets with Fluoride Ion Affinities (FIA) values (computed and predicted) for molecules, represented as SMILES strings. Quantum descriptors featurization of the database is also provided.
- results/: Stores plots and tables published in the main text and supplementary materials of the article.
- notebooks/: Contains Jupyter notebooks for reproducing results. These notebooks can be run independently and provide an interactive way to explore the data and models.


## Usage
Use the Jupyter notebooks in the notebooks/ folder to reproduce analyses or customize the workflows for new data. For example, this material can be used to unravel structure-property relationships for another chemical problem (other chemical property, reaction selectivities, yields...).


