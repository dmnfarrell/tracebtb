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
import math,random,time
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

def scale_value(x, old_max=15, new_min=1, new_max=10):
    """Scale a value to a score range"""
    if x > old_max:
        return 0
    # Calculate the scale factor (9 / 15 = 0.6)
    scale_factor = (new_max - new_min) / old_max
    # Invert the logic: start at max and subtract based on x
    score = new_max - (x * scale_factor)
    # Optional: Clamp the results between 1 and 10
    return max(new_min, min(new_max, score))

def score_within_herd(animal_row, context, snp_dist=None):
    """
    Scores the likelihood of within-herd transmission (0-10) using SNP5 cluster size
    within the current herd breakdown.
    Args:
        animal_row: The specific animal being scored, pandas Series
        context: The pre-computed dictionary for the animal's herd.
    Returns:
        An integer score (0-10).
    """

    snp5_cluster = animal_row['snp5']
    snp12_cluster = animal_row['snp12']
    is_homebred = animal_row['Homebred']
    metrics = context['metrics']
    data = context['data']
    if metrics['herd_isolates'] == 1:
        # Only one isolate in the herd. Cannot confirm or rule out within-herd spread.
        return 1, "Singleton"

    # Find all animals that share the same clusters IDs as current animal.
    sub = data['herd_isolates']
    cluster5 = sub[sub['snp5'] == snp5_cluster]
    cluster12 = sub[sub['snp12'] == snp12_cluster]

    #are animals slaughtered within 60 days of each other?
    print (cluster5.groupby('last_move')['snp5'].count())
    #filter for rows within 60 days of any other

    if len(cluster5) <= 1 and len(cluster12)<=1:
        return 0, "Multi isolate herd with no strain match"
    try:
        df = cluster5.sort_values('last_move')
        df['days_diff'] = df['last_move'].diff().dt.days
        close_rows = df[df['days_diff'] <= 90]
        if len(close_rows)==0:
            #animals are >60 days apart in death - need total overlap time too?
            return 3, 'matching isolates but >90 days apart'
    except:
        pass

    #switch scoring to use scale based on snp dist?
    if snp_dist is not None:
        names = list(cluster12['sample'])
        dm = snp_dist.loc[names, names]
        md = dm.stack().min()
        #print (dm)
        print (f'snp dist={md},score={math.ceil(scale_value(md))}')

    if len(cluster5) > 1:
        # Strong evidence: The animal shares a recent genetic ancestor
        # (SNP5 cluster) with another animal currently in the herd.
        return 10, f"SNP5 Cluster match, size={len(cluster5)}"
    elif len(cluster12) > 1:
        # Medium evidence, shares less recent ancestor
        return 5, f"SNP12 Cluster match, size={len(cluster12)}"

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
    snp5_cluster = animal_row['snp5']
    snp12_cluster = animal_row['snp12']
    last_move = animal_row.last_move
    data = context['data']
    sub = data['herd_isolates']
    sr = data['std_reactors']
    #print (sr)
    #print (sr.loc[years])

    # --- Check testing first ---
    if sr is None or sr.sum()==0:
        # No known previous breakdowns.
        return 0, 'No known previous breakdowns'

    # find Historical Strains - prior to this animals breakdown
    cols = ['snp5','snp12','last_move']
    #print (sub.last_move)
    historic = sub[(sub.last_move<last_move)][cols]
    print (tag,curr_year,snp5_cluster,snp12_cluster,last_move)
    #print (historic)
    previous_snp5 = set(historic.snp5)
    previous_snp12 = set(historic.snp12)
    #print (snp5_cluster,previous_snp5)
    if len(historic) == 0:
        return 1, 'No previous breakdown prior to this sample'
    elif snp5_cluster in previous_snp5:
        return 10, 'Current breakdown strain same as previous at snp5 level'
    elif snp12_cluster in previous_snp12:
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
    nb_isolates = data['neighbour_isolates']
    species = set(nb_isolates.Species)
    neighbour_snp5 = list(nb_isolates.snp5)
    neighbour_snp12 = nb_isolates.snp12
    #print (f'neighbour snp5={neighbour_snp5}')

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

def score_movement(animal_row, context, lpis_cent, moves_in=None):
    """
    Scores the likelihood of long distance movement-based transmission (0-10)
    by checking the herds movement history against any related strain isolates.

    Args:
        animal_row: The specific animal being scored.
        context: The pre-computed dictionary for the animal's current herd.
        moves_in: All moves into herd can be provided, else we calculate

    Returns:
        An integer score (0-10).
    """

    data = context['data']
    sub = data['herd_isolates']
    herd = animal_row.HERD_NO
    tag = animal_row.Animal_ID
    last_move = animal_row.last_move
    snp5_cluster = animal_row['snp5']
    snp12_cluster = animal_row['snp12']
    data = context['data']
    #first check if any strains outside local area

    x = data['snp5_related']
    r5 = x[x.snp5==snp5_cluster]
    x = data['snp12_related']
    r12 = x[x.snp12==snp12_cluster]
    print (herd, snp5_cluster, len(r5),len(r12), last_move)
    target_herds = list(r5.HERD_NO)
    #print (target_herds)

    if len(r5) == 0:# and len(r12) == 0:
        return 1, 'No related strains outside herd at snp5 or snp12'
    elif len(r5) >= 1:
        #get all moves from sources to herd in 1 year prior to breakdown
        try:
            date1 = last_move - relativedelta(months=24)
        except:
            return 1, 'Cannot parse last_move date'
        if moves_in is None:
            moves_in = movement.query_herd_moves_all(herd, date1, last_move)
            v = movement.get_movement_vectors(moves_in, lpis_cent)
            #remove local moves
            v = v[v.dist_km>=15]

        animal_moves = moves_in[moves_in.tag==tag]
        #print (animal_moves)
        any_from_strain_herds = v[v.herd.isin(target_herds)]
        sample_from_strain_herds =  v[(v.herd.isin(target_herds)) & (v.tag==tag)]
        #print (any_from_strain_herds)
        if len(sample_from_strain_herds)>0:
            #did animal move directly from a strain herd?
            return 10, 'Direct move of sampled animal from herd with same strain'
        elif len(any_from_strain_herds)>0:
            #any direct moves of any animals from herds with related isolates
            return 7, 'Move into this herd from another herd with same strain'
        #elif
        #strains same as those in herd in intermediate herds of moves

        return 3, 'Related strain in herd >10km but no move'

def run_pathways(herds, *args):
    """Run attribution on a set of herds"""

    results = []
    for herd_no in herds:
        result = run_pathway(herd_no, *args)
        results.append(result)
    return pd.concat(results)

def run_pathway(herd_no, gdf, moves, testing, feedlots,
                 lpis, lpis_cent, grid, hc=None, moves_in=None,
                 dist=4000, meta_cols=False):
    """Run attribution on a single herd"""

    print (herd_no)
    #default cols from metadata
    cols=['snp3','snp5','CFU','in_homerange','Homebred','Year']
    # Build context for the single herd
    if hc is None:
        hc = tools.get_herd_context(herd_no, gdf, moves, testing, feedlots,
                 lpis, lpis_cent, grid=grid, dist=dist)
    if hc is None:
        return
    data = hc['data']
    sub = data['herd_isolates']
    results_list = []
    for _, row in sub.iterrows():
        s_within, e_within = score_within_herd(row, hc)
        #print (herd_no, row.Animal_ID, s_within)
        s_local, e_local = score_local(row, hc)
        s_movement, e_movement = score_movement(row, hc, lpis_cent)
        s_residual, e_residual = score_residual(row, hc)

        # Package the final score results
        scores= {'W': s_within,
                'L': s_local,
                'M': s_movement,
                'R': s_residual}
        evidence = {
                'W': e_within,
                'L': e_local,
                'M': e_movement,
                'R': e_residual,
                'I': 'NA' }
        '''pathway = max(scores, key=scores.get)
        if max(scores.values()) == 1:
            pathway = 'U'
        elif len(set(scores.values()))==1:
            #multiple non-zero ties - hand this.
            pathway = 'I'
        '''
        # Find the maximum score
        max_val = max(scores.values())
        # Identify if any pathways share the maximum
        winners = [k for k, v in scores.items() if v == max_val]
        # Apply conditional logic
        if max_val <= 1:
            pathway = 'U'
        elif len(winners) > 1:
            # Sort winners alphabetically to ensure consistency
            pathway = "".join(sorted(winners))
        else:
            # Only one clear winner
            pathway = winners[0]

        result = {'HERD_NO': row.HERD_NO,
                  'Animal_ID': row.Animal_ID,
                  'snp5':row.snp5}
        if meta_cols == True:
            #add other cols
            scols=[c for c in cols if c in sub.columns]
            result = result.update(row[scols].to_dict())
        result.update(scores)
        result['best_pathway'] = pathway
        for key in evidence.keys():
            if key != 'I':
                result['e_'+key] = evidence[key]

        #print (scores.values(),max(scores), pathway)
        results_list.append(result)
    return pd.DataFrame(results_list).set_index('Animal_ID')

def summarize_pathways(df):
    """
    Summarizes herd data focusing on pathway involvement, tie-break
    concurrence, and genomic data gaps (orphans).
    """
    total_animals = len(df)
    unique_herds = df['HERD_NO'].unique()

    # Path columns for score/orphan checking
    score_cols = ['W', 'L', 'M', 'R']

    # Mapping for clear report labels
    mapping = {
        'W': 'Within-Herd',
        'L': 'Local',
        'M': 'Movement',
        'R': 'Residual',
        'U': 'Unknown (Low Score)'
    }

    lines = []
    lines.append("=" * 60)
    lines.append(f"{'GENOMIC PATHWAY SUMMARY REPORT':^60}")
    lines.append("=" * 60)
    lines.append(f"Total Animals: {total_animals}")
    lines.append(f"Total Herds:   {len(unique_herds)}")

    # 1. Best Pathway Distribution (Assigned Categories)
    lines.append("\n[1. Top Pathway Attribution]")
    path_counts = df['best_pathway'].value_counts()
    for path, count in path_counts.items():
        # Build label for combos (e.g., 'WL' -> 'Within-Herd+Local')
        if path == 'U':
            label = mapping['U']
        else:
            label = "+".join([mapping.get(c, c) for c in path])

        pct = (count / total_animals) * 100
        lines.append(f"- {label:<35}: {count:>3} animals ({pct:>5.1f}%)")

    # 2. Total Pathway Involvement (Weighted Lead Analysis)
    # Counts how many times a pathway was a 'winner', including within ties
    lines.append("\n[2. Total Pathway Involvement]")
    lines.append(" (How often each route was a top-tier suspect")

    involvement = {k: 0 for k in score_cols}
    for path in df['best_pathway']:
        if path == 'U': continue
        for char in path:
            if char in involvement:
                involvement[char] += 1

    for char, count in involvement.items():
        pct = (count / total_animals) * 100
        lines.append(f"- {mapping[char]:<15}: {count:>3} cases ({pct:>5.1f}%)")

    # 3. Genomic Gap Analysis (Orphan Check)
    lines.append("\n[3. Genomic Gap Analysis]")
    # Identifies animals where all scored pathways resulted in exactly 0
    orphans = df[(df[score_cols] <= 1).all(axis=1)]
    orphan_count = len(orphans)
    orphan_pct = (orphan_count / total_animals) * 100

    lines.append(f"- 'Genomic Orphans': {orphan_count} ({orphan_pct:.1f}%)")
    lines.append("  (Isolates with zero genetic matches in any known route)")

    # 4. Per-Herd Breakdown
    if len(unique_herds) > 1:
        lines.append("\n" + "=" * 60)
        lines.append(f"{'PER-HERD TRENDS':^60}")
        lines.append("=" * 60)

        # Aggregate by herd to see primary driver and tie frequency
        herd_summary = df.groupby('HERD_NO').agg({
            'best_pathway': [
                ('Primary', lambda x: x.mode()[0] if not x.empty else "N/A"),
                ('Tied_Cases', lambda x: (x.str.len() > 1).sum())
            ],
            'W': 'count' # Just to get a 'Size' column
        })

        # Flatten multi-index columns for cleaner printing
        herd_summary.columns = ['Primary Pathway', 'Tied Cases', 'Herd Size']
        lines.append(herd_summary[['Herd Size', 'Primary Pathway', 'Tied Cases']].to_string())

    lines.append("\n" + "=" * 60)
    return "\n".join(lines)

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