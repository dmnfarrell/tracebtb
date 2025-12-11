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

def get_herd_context(herd_no, metadata, parcels, moves, testing, lpis_cent, dist=5000):
    """
    Computes and returns all WGS and epidemiological context for a single target herd.
    Args:
        herd_no: The HERD_NO to process
        metadata: geodataframe of all metadata
        moves: relevant movement data
        testing: per herd testing data
        lpis_cent: centroid points of all herds in lpis
        dist: distance at which to consider neighbours
    Returns:
        A dictionary containing all context elements for the target herd.
    """

    # --- A. Current Herd Data (Isolates) ---
    herd_isolates = metadata[metadata['HERD_NO'] == herd_no].copy()
    # --- B. Neighboring Herds Data ---
    if herd_isolates is not None:
        neighbour_data = tools.find_neighbouring_points(metadata, herd_isolates, dist)
        neighbour_strains = set(neighbour_data['snp5'].dropna())
    else:
        neighbour_data = None
        neighbour_strains = None
    # --- C. Movement Data ---
    # Get moves into this herd for all sampled animals
    herd_movements = tools.get_moves_bytag(herd_isolates, moves, lpis_cent)
    if herd_movements is not None:
        herd_movements = herd_movements.reset_index()
    # Get all moves into the herd in past 5 years! TODO
    all_movements = None
    # 3. Assemble context dict
    current_strains = set(herd_isolates['short_name'].dropna())
    current_snp5 = set(herd_isolates['snp5'].dropna())
    is_singleton: bool = len(herd_isolates) <= 1
    if herd_no in testing.index:
        te = testing.loc[herd_no]
    else:
        te = None

    herd_context = {
        # Core data for this breakdown
        'herd_isolates': herd_isolates,
        'is_singleton': is_singleton,
        # WGS context for Residual/Within-Herd
        'current_strains': current_strains,
        'current_snp5': current_snp5,
        # Spatial Context (Local Pathway)
        'neighbour_herds': neighbour_data,
        'neighbour_strains': neighbour_strains,
        # Movement Context (Movement Pathway)
        'moves': herd_movements,
        'breakdown_history': te
        #
        #'herd_summary': herd_summary
        # NOTE: Residual (Past Infections) requires historical data lookup here
        # E.g., 'historical_strains': tools.get_historic_strains(target_herd_no)
    }
    return herd_context

def score_within_herd(animal_row, context) -> int:
    """
    Scores the likelihood of within-herd transmission (0-10) using SNP5 cluster size
    within the current herd breakdown.
    Args:
        animal_row: The specific animal being scored, dict or Series
        context: The pre-computed dictionary for the animal's herd.
    Returns:
        An integer score (0-10).
    """

    my_cluster = animal_row['snp5']
    is_homebred = animal_row['Homebred']

    # Use pre-computed flag for singleton status
    if context['is_singleton']:
        # Only one isolate in the herd. Cannot confirm or rule out within-herd spread.
        return 0, "Singleton"

    # 1. Calculate Cluster Size within the Herd
    # Filter the subset of metadata (isolates_metadata) to find all animals
    # that share the same SNP5 cluster ID as the current animal.
    herd_sub_df = context['herd_isolates']
    cluster_size = len(herd_sub_df[herd_sub_df['snp5'] == my_cluster])
    # 2. Apply Scoring Logic (0-10)
    if cluster_size > 1:
        # Strong evidence: The animal shares a recent genetic ancestor
        # (SNP5 cluster) with another animal currently in the herd.
        return 9, "Cluster match"
    else:
        # Negative evidence: The animal is a singleton cluster in a multi-isolate herd.
        # Implies introduction, not internal spread.
        return 0, "Multi isolate herd with no strain match"
    # 3. Apply Homebred Bonus (optional)
    #if is_homebred is True and score < 10:
        # Homebred status slightly increases the plausibility of a local/within source.
    #    score += 1
    return

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
    nb = context['neighbour_herds']
    # Check if any neighboring herds were identified in the context setup
    if len(nb)==0:
        # Cannot attribute to local spread if no neighbors are defined or found
        return 0, "Singleton sample"

    # Retrieve the pre-computed set of unique SNP5 clusters found in all neighbors
    neighbor_clusters = context['neighbour_strains']

    # 1. Check for Strong Evidence: Exact SNP5 Cluster Match
    if my_cluster in neighbor_clusters:
        # High likelihood: The animal's specific, highly related cluster
        # is circulating in a neighbouring herd.
        return 9, f"Cluster Match ({neighbor_clusters})"

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

    sub = context['herd_isolates']
    my_animal_id = animal_row['Animal_ID']
    my_cluster = animal_row['snp5']
    current_herd = animal_row['HERD_NO']
    #print (context['moves'])
    herd_moves = context['moves']
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
    return max_score, "Default score"

def run_attribution(herd_df, gdf, parcels, moves, meta_cols=False):
    """Run attribution on a set of herds"""

    global_strain_map = gdf.groupby('HERD_NO')['snp5'].apply(lambda x: set(x.dropna())).to_dict()

    results_list = []
    #default cols from metadata
    cols=['snp3','snp5','CFU','in_homerange','Homebred','Year']

    for i, r in herd_df.iterrows():
        herd_no = r['HERD_NO']
        # Build context for the single herd
        hc = get_herd_context(herd_no, gdf, parcels, moves)
        sub = hc['herd_isolates']

        #run scoring per animal
        for _, row in sub.iterrows():
            s_within, e_within = score_within_herd(row, hc)
            #print (herd_no, row.Animal_ID, s_within)
            s_local, e_local = score_local(row, hc)
            s_movement, e_movement = score_movement(row, hc, global_strain_map)
            # s_residual = score_residual(animal_row, current_context)

            # Package the final score results
            scores= {'within': s_within,
                    'local': s_local,
                    'movement': s_movement}
                    #'residual': s_residual}
            evidence = {
                    'within': e_within,
                    'local': e_local,
                    'movement': e_movement,
                    'indeterminate': 'NA'
                    #'Residual': e_resid
                }
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
            result['estimated_pathway'] = pathway
            result['evidence'] = evidence[pathway]
            #print (scores.values(),max(scores), pathway)
            results_list.append(result)
    return pd.DataFrame(results_list)

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