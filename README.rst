=============================
AMRtime
=============================

.. image:: https://badge.fury.io/py/AMRtime.png
    :target: http://badge.fury.io/py/AMRtime

.. image:: https://travis-ci.org/fmaguire/AMRtime.png?branch=master
    :target: https://travis-ci.org/fmaguire/AMRtime

Metagenomic AMR gene detection using hierarchical machine learning models based on either
AMR curated gene families or sequence identity based clusters.


1. Initial read analysis using heuristically-accelerated homology searches (DIAMOND)

2. Xgboost based classification of metagenomics reads to AMR gene families or sequence clusters

3. Xgboost based classification of reads within family/clusters to specific genes

4. Localised assembly of classified reads

5. Attempted extension of family/cluster member-specific contigs 


Installation
--------

External dependencies:

- DIAMOND

- vsearch (cd-hit actually used)

- art

- scikit-learn

- biopython

- tqdm

