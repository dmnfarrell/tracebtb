#!/usr/bin/env python

"""
    Source attribution module for tracebtb.
    Created Nov 2025
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys,os,subprocess,glob,shutil,re
import pickle
import random,time
import json
import platform
import numpy as np
import pandas as pd
import pylab as plt
import matplotlib as mpl
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO, Align
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from . import core, tools

module_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(module_path,'data')

SNP_BINS = {
    "snp3":  (0, 3),
    "snp5":  (0, 5),
    "snp12": (0, 12),
}
SPATIAL_RADIUS = 5000
MOVEMENT_WINDOW_YEARS = 2

from typing import Dict, Any, List, Set

def score_within_herd(animal_row, context):
    """
    Scores the likelihood of within-herd transmission (0-10) using SNP5 cluster size
    within the current herd breakdown.
    Args:
        animal_row: The specific animal being scored, pandas Series
        context: The pre-computed dictionary for the animal's herd.
    Returns:
        An integer score (0-10).
    """

    my_cluster = animal_row['snp5']
    is_homebred = animal_row['Homebred']
    metrics = context['metrics']
    data = context['data']
    if metrics['herd_isolates'] == 1:
        # Only one isolate in the herd. Cannot confirm or rule out within-herd spread.
        return 1, "Singleton"

    # Find all animals that share the same SNP5 and SNP12 cluster IDs as the current animal.
    herd_sub_df = data['herd_isolates']
    cluster_size5 = len(herd_sub_df[herd_sub_df['snp5'] == my_cluster])
    cluster_size12 = len(herd_sub_df[herd_sub_df['snp12'] == my_cluster])
    # Apply Scoring Logic for clusters
    if cluster_size5 > 1:
        # Strong evidence: The animal shares a recent genetic ancestor
        # (SNP5 cluster) with another animal currently in the herd.
        return 10, f"SNP5 Cluster match, size={cluster_size5}"
    elif cluster_size12 > 1:
        # Medium evidence, shares less recent ancestor
        return 5, f"SNP12 Cluster match, size={cluster_size12}"
    return 0, "Multi isolate herd with no strain match"

def score_residual(animal_row, context):
    """
    Scores the likelihood of residual within-herd transmission (0-10) by
    checking if the current SNP5 cluster has been previously isolated in this herd
    in past breakdowns.
    Args:
        animal_row: The specific animal being scored, pandas Series
        context: The pre-computed dictionary for the animal's herd.
    Returns:
        An integer score (0-10).
    """

    tag = animal_row['sample']
    my_herd = animal_row['HERD_NO']
    curr_year = animal_row['Year']
    snp_cluster5 = animal_row['snp5']
    snp_cluster12 = animal_row['snp12']
    data = context['data']
    sub = data['herd_isolates']
    sr = data['reactor_history']

    # Look up Historical Strains - prior to current year
    cols = ['snp5','snp12','Year']
    historic = sub[(sub.Year<curr_year)][cols]

    years = sub.Year.dropna().unique()
    print (tag,curr_year,snp_cluster5,snp_cluster12)

    # --- Step 2: Apply Scoring Logic ---
    if sr is None or sr.sum()==0:
        # No known previous breakdowns.
        return 0, 'No known previous breakdowns'

    previous_snp5 = historic.snp5
    previous_snp12 = historic.snp12
    if len(historic) == 0:
        return 1, 'No previous breakdown prior to this sample'
    elif snp_cluster5 in previous_snp5:
        # Strong evidence: The current breakdown strain is a close genetic match
        # to a strain that caused a previous breakdown
        return 9, 'Current breakdown strain same as previous at snp5 level'
    elif snp_cluster12 in previous_snp12:
        # Medium evidence: The current breakdown strain is a match on snp12 level
        return 5, 'Current breakdown strain same as previous at snp12 level'
    else:
        return 0, 'No evidence of residual infection in past breakdowns'

def score_local(animal_row, context):
    """
    Scores the likelihood of local area transmission (0-10) by checking if
    the animal's SNP5 cluster is present in a neighboring herd.
    Args:
        animal_row: The specific animal being scored.
        context: The pre-computed dictionary for the animal's herd.
    Returns:
        An integer score (0-10).
    """

    my_cluster = animal_row['snp5']
    data = context['data']
    nb = data['neighbour_parcels']
    # Check if any neighboring herds were identified in the context setup
    if len(nb) == 0:
        # Cannot attribute to local spread if no neighbors are defined or found
        return 0, "No neighbouring herds"

    # Retrieve the pre-computed set of unique SNP5 clusters found in all neighbors
    neighbour_clusters = data['neighbour_strains']

    # 1. Check for Strong Evidence: Exact SNP5 Cluster Match
    if my_cluster in neighbour_clusters:
        # High likelihood: The animal's specific, highly related cluster
        # is circulating in a neighbouring herd.
        return 9, f"Cluster match ({neighbour_clusters})"

    # 2. Low Evidence: No Match Found
    # Neighbours exist, but they do not share this specific genetic cluster.
    # We assign a baseline score of 1 to ensure this pathway is considered,
    # but highly unlikely to be the primary source.
    return 1, "Neighbours exist but no cluster shared"

def score_movement(animal_row, context, global_strain_map):
    """
    Scores the likelihood of movement-based transmission (0-10) by checking
    the herds movement history against the WGS data of its source herds.

    Args:
        animal_row: The specific animal being scored.
        context: The pre-computed dictionary for the animal's current herd.
        global_strain_map: Global lookup for all herd strains/clusters.

    Returns:
        An integer score (0-10).
    """

    data = context['data']
    sub = data['herd_isolates']
    my_animal_id = animal_row['Animal_ID']
    my_cluster = animal_row['snp5']
    current_herd = animal_row['HERD_NO']
    #print (context['moves'])
    data = context['data']
    herd_moves = data['isolate_moves']
    if herd_moves is None:
        return 0, "No herd moves"
    # 1. Homebred Check (Must be the highest priority)
    #if animal_row['Homebred'] is True:
    #    return 0
    max_score = 1 # Baseline score if movement happened but no match is found

    # Get the animal's movement history
    animal_moves = herd_moves[herd_moves.tag == my_animal_id]

    # --- Check for Scenario A (Direct Link) ---
    prior_herds = set(animal_moves.move_to.dropna())
    prior_herds.discard(current_herd)
    source_strains = []
    for source_herd in prior_herds:
        #print (source_herd)
        # Look up strains in the historical source herd
        s=global_strain_map.get(source_herd)
        if s is not None:
            source_strains.extend(s)
    #print (prior_herds, source_strains)
    if my_cluster in source_strains:
        # Strongest possible evidence: Direct move from a known herd with same strain
        return 10, "Direct move from X"

    # --- Check for Scenario C (Passive Introduction by other animals) ---

    # Use all animals that moved into the current herd n years prior


    # Look up the strains associated with these incoming animals
    # (requires linking 'tag' back to WGS metadata, or using a move-specific strain map)

    # We simplify this check by looking for ANY incoming animal that matches the cluster.
    # (Requires a slightly augmented context or a helper function)

    # Simplified Logic for C (Requires further helper implementation):
    # If the animal's cluster matches the strain of ANY other animal moved
    # into the current herd (and sequenced):
    # if my_cluster in tools.get_strains_from_incoming_moves(all_incoming_moves, metadata_df):
    #     max_score = max(max_score, 6)

    # For now, let's keep the primary focus on the strongest, direct link (Scenario A).
    # Since you will implement the logic for the complex links later,
    # we return the best score found.
    return max_score, ''

def run_pathways(herd_df, gdf, moves, testing, feedlots,
                     lpis, lpis_cent, grid, meta_cols=False):
    """Run pathway attribution on a set of herds"""

    global_strain_map = gdf.groupby('HERD_NO')['snp5'].apply(lambda x: set(x.dropna())).to_dict()
    results_list = []
    #default cols from metadata
    cols=['snp3','snp5','CFU','in_homerange','Homebred','Year']

    for i, r in herd_df.iterrows():
        herd_no = r['HERD_NO']
        # Build context for the single herd
        hc = tools.get_herd_context(herd_no, gdf, moves, testing, feedlots,
                     lpis, lpis_cent, grid=grid, dist=1000)
        sub = hc['data']['herd_isolates']
        #print (hc)
        #run scoring per animal
        for _, row in sub.iterrows():
            s_within, e_within = score_within_herd(row, hc)
            #print (herd_no, row.Animal_ID, s_within)
            s_local, e_local = score_local(row, hc)
            s_movement, e_movement = score_movement(row, hc, global_strain_map)
            s_residual, e_residual = score_residual(row, hc)

            # Package the final score results
            scores= {'within': s_within,
                    'local': s_local,
                    'movement': s_movement,
                    'residual': s_residual}
            evidence = {
                    'within': e_within,
                    'local': e_local,
                    'movement': e_movement,
                    'residual': e_residual,
                    'indeterminate': 'NA' }
            pathway = max(scores, key=scores.get)
            if max(scores.values()) == 0:
                pathway = 'indeterminate'

            result = {'HERD_NO': row.HERD_NO,
                      'Animal_ID': row.Animal_ID}
            if meta_cols == True:
                #add other cols
                scols=[c for c in cols if c in sub.columns]
                result = result.update(row[scols].to_dict())
            result.update(scores)
            result['best_pathway'] = pathway
            for key in evidence.keys():
                if key != 'indeterminate':
                    result['e_'+key] = evidence[key]
            result['evidence'] = evidence[pathway]
            #print (scores.values(),max(scores), pathway)
            results_list.append(result)
    return pd.DataFrame(results_list)

def summarize_pathways(df):
    """
    Prints a formatted summary of animal movement and infection pathways.
    """
    total_animals = len(df)
    unique_herds = df['HERD_NO'].nunique()

    print(f"SUMMARY REPORT")
    print("-" * 40)
    print(f"Herd:     {df['HERD_NO'].iloc[0]}")
    print(f"Total Animals Processed: {total_animals}")
    print("-" * 40)

    # 1. Best Pathway Breakdown
    print("\n[Best Pathway Distribution]")
    pathway_counts = df['best_pathway'].value_counts()
    for pathway, count in pathway_counts.items():
        percentage = (count / total_animals) * 100
        print(f"- {pathway:15}: {count} ({percentage:.1f}%)")

    # 2. Evidence Summary (Top explanations)
    print("\n[Primary Evidence/Explanations]")
    # We aggregate the 'evidence' column to see the most common reasons
    evidence_counts = df['evidence'].value_counts().head(5)
    for reason, count in evidence_counts.items():
        print(f"- {count} animal(s): {reason}")

    # 3. Binary Indicator Totals
    print("\n[Binary Flags (Presence)]")
    indicators = ['within', 'local', 'movement', 'residual']
    for col in indicators:
        if col in df.columns:
            sum_val = df[col].sum()
            print(f"- {col.title():10}: {sum_val} positive hits")

    print("-" * 40)

def calculate_benchmark_metrics(results_df):
    """
    Calculates agreement (Kappa) between 'True_Pathway' and 'Final_Pathway'.
    Ignores 'Reference' rows (historical anchors/sources).
    """
    # Filter out reference animals (sources/historical anchors) that aren't being scored
    test_set = results_df[results_df['true_pathway'] != 'Reference'].copy()

    y_true = test_set['true_pathway']
    y_pred = test_set['estimated_pathway']

    # Calculate Cohen's Kappa
    kappa = cohen_kappa_score(y_true, y_pred)

    print("--- Benchmarking Results ---")
    print(f"Cohen's Kappa: {kappa:.4f}")
    print("\nDetailed Classification Report:")
    print(classification_report(y_true, y_pred, zero_division=0))

    return kappa