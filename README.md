# Heredity Genetics Probability Model (CS50 AI)

## What it does
* This is a Python program that assesses the probability of a genetic trait (such as gene counts and expressed traits) being passed down through a family tree, using joint probability and Bayesian inference.
* It parses family data from CSV files, computes joint probabilities for gene inheritance and trait expression, and normalizes the final probability distributions for each person.

## Technologies Used
* Python
* Probability and Statistics / Bayesian Inference
* CSV data processing and iterative combinatorics (`itertools`)
* Git / GitHub

## How to Run It
* Make sure you have your family data CSV file (e.g., `family0.csv`) in your directory.
* Run the main script with the data file path:
* "```bash"
* python heredity.py data/family0.csv
* View the calculated probability distributions for genes and traits for every individual in the family tree.
