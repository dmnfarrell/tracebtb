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
from dateutil.relativedelta import relativedelta

from . import tools, movement

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

    snp5_cluster = animal_row['snp5']
    snp12_cluster = animal_row['snp5']
    data = context['data']
    nb = data['neighbour_parcels']
    badgers = data['near_badger']

    print (len(nb),snp5_cluster,snp12_cluster,len(badgers))
    # Check if any neighboring herds were identified in the context setup
    if len(nb) == 0:
        # Cannot attribute to local spread if no neighbors are defined or found
        return 0, "No neighbouring herds"

    # unique SNP5 clusters found in all neighbors
    nb_isolates = data['neighbour_herd_isolates']
    species = set(nb_isolates.Species)
    neighbour_snp5 = list(nb_isolates.snp5)
    neighbour_snp12 = nb_isolates.snp12
    print (f'neighbour snp5={neighbour_snp5}')

    # 1. Check for Strong Evidence: Exact SNP5 Cluster Match
    if len(nb_isolates) == 0:
        # 2. Low Evidence: No Match Found
        # Neighbours exist, but they do not share this specific genetic cluster.
        # We assign a baseline score of 1 to ensure this pathway is considered,
        return 1, f"No sampled herds in {len(nb)} neighbours"
    elif snp5_cluster in neighbour_snp5:
        # High likelihood: The animal's specific, highly related cluster
        # is circulating in a neighbouring herd.
        if 'badger' in species:
            print ('badger link')
        return 10, f"Cluster match at snp5 level ({snp5_cluster})"
    elif snp12_cluster in neighbour_snp12:
        # High likelihood: The animal's specific, highly related cluster
        # is circulating in a neighbouring herd.
        return 5, f"Cluster match at snp12 level ({snp12_cluster})"
    return 0, f"{len(nb_isolates)} neighbouring isolates found with no related cluster"

def score_movement(animal_row, context):
    """
    Scores the likelihood of long distance movement-based transmission (0-10)
    by checking the herds movement history against any related strain isolates.

    Args:
        animal_row: The specific animal being scored.
        context: The pre-computed dictionary for the animal's current herd.

    Returns:
        An integer score (0-10).
    """

    data = context['data']
    sub = data['herd_isolates']
    my_animal_id = animal_row['Animal_ID']
    herd = animal_row.HERD_NO
    tag = animal_row.Animal_ID
    last_move = animal_row.last_move
    snp5_cluster = animal_row['snp5']
    snp12_cluster = animal_row['snp12']

    data = context['data']
    herd_moves = data['isolate_moves']

    #first check if any strains outside local area
    #all herds with snp5 cluster

    x = data['snp5_related']
    r5 = x[x.snp5==snp5_cluster]
    x = data['snp12_related']
    r12 = x[x.snp12==snp12_cluster]
    #remove local herds - haven't done this yet
    print (herd, snp5_cluster, len(r5),len(r12))

    target_herds = list(r5.HERD_NO)
    #print (target_herds)
    # Get the animals own movement history - do we need it now?
    animal_moves = herd_moves[herd_moves.tag == my_animal_id]
    #print (r5[['HERD_NO','snp5','snp12']])

    if len(r5) == 0 and len(r12) == 0:
        return 1, 'no related strains outside herd at snp5 or snp12'
    elif len(r5) > 1:
        #get all moves from sources to herd in 1 year prior to breakdown
        #these really only need to be found once her perd ...?
        date1 = last_move - relativedelta(years=1)
        moves = movement.query_herd_moves_all(herd, date1, last_move)
        animal_moves = moves[moves.tag==my_animal_id]
        direct_moves = moves[(moves.herd==herd) & (moves.tag!=my_animal_id)]
        rel_moves = moves[moves.herd!=herd]

        print (f'{len(animal_moves)} sampled moves;{len(direct_moves)} direct;{len(rel_moves)} related')

        sampled_from_strain_herds = animal_moves[animal_moves.herd.isin(target_herds)]
        all_from_strain_herds = direct_moves[direct_moves.herd.isin(target_herds)]
        if len(sampled_from_strain_herds)>0:
            #did animal move directly from a strain herd?
            return 10, 'direct move of sampled animal from herd with same strain'
        elif len(all_from_strain_herds)>0:
            #any direct moves of any animals from herds with related isolates
            return 7, 'move of other animals in herd from another herd with same strain'
        #elif
        #strains same as those in herd in intermediate herds of moves

        return 3, 'related strain in herd >10km but no move'

def run_pathways(herds, gdf, moves, testing, feedlots,
                     lpis, lpis_cent, grid, meta_cols=False, dist=4000):
    """Run attribution on a set of herds"""

    results_list = []
    #default cols from metadata
    cols=['snp3','snp5','CFU','in_homerange','Homebred','Year']

    for herd_no in herds:
        # Build context for the single herd
        hc = tools.get_herd_context(herd_no, gdf, moves, testing, feedlots,
                     lpis, lpis_cent, grid=grid, dist=dist)
        if hc is None:
            continue
        data = hc['data']
        sub = data['herd_isolates']
        #print (hc)
        #run scoring per animal
        for _, row in sub.iterrows():
            s_within, e_within = score_within_herd(row, hc)
            #print (herd_no, row.Animal_ID, s_within)
            s_local, e_local = score_local(row, hc)
            s_movement, e_movement = score_movement(row, hc)
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
            if len(set(scores.values()))==1 or max(scores.values()) == 0:
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