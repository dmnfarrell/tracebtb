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

    days_diff = 150
    snp5_cluster = animal_row['snp5']
    snp10_cluster = animal_row['snp10']
    is_homebred = animal_row['Homebred']
    last_move = animal_row.last_move
    metrics = context['metrics']
    data = context['data']
    sr = data['std_reactors']
    lp = data['lesion_positives']

    if metrics['herd_isolates'] == 1:
        # Only one isolate in the herd.
        #check lesion positives in same breakdown - can't tell?
        return 1, "Singleton"

    # Find all animals that share the same clusters IDs as current animal.
    sub = data['herd_isolates']
    cluster5 = sub[sub['snp5'] == snp5_cluster]
    cluster10 = sub[sub['snp10'] == snp10_cluster]

    print (animal_row.Animal_ID,snp5_cluster,snp10_cluster,last_move)

    if len(cluster5) <= 1 and len(cluster10)<=1:
        return 0, "Multi isolate herd with no strain match"
    if len(cluster5) > 1:
        df = cluster5.sort_values('last_move')
        df['days_diff'] = df['last_move'].diff().dt.days
        close_rows = df[df['days_diff'] < days_diff]
        if len(close_rows)>0:
            return 10, f"SNP5 Cluster match, size={len(cluster5)}"
        else:
            return 1, f'Isolates >{days_diff} days apart'
    elif len(cluster10) > 1:
        df = cluster10.sort_values('last_move')
        df['days_diff'] = df['last_move'].diff().dt.days
        close_rows = df[df['days_diff'] < days_diff]
        #print(df[['Animal_ID','snp10','days_diff']])
        if len(close_rows)>0:
            return 5, f"SNP10 Cluster match, size={len(cluster10)}"
        else:
            return 1, f'Isolates >{days_diff} days apart'

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

    days_diff = 150
    tag = animal_row['sample']
    my_herd = animal_row['HERD_NO']
    curr_year = animal_row['Year']
    snp5_cluster = animal_row['snp5']
    snp10_cluster = animal_row['snp10']
    last_move = animal_row.last_move
    data = context['data']
    sub = data['herd_isolates']
    sr = data['std_reactors']
    #print (sr)
    #print (sr.loc[years])

    # --- Check testing first ---
    if sr is None or sr.sum()==0:
        # No known previous breakdowns.
        return 1, 'No known previous breakdowns found'

    # find Historical Strains - prior to this animals breakdown
    # must be >70 days prior to death date
    cols = ['snp5','snp10','last_move']
    #print (sub.last_move)
    try:
        last_date = last_move-relativedelta(days=days_diff)
    except:
        last_date = last_move
    historic = sub[(sub.last_move<last_date)][cols]
    #print (tag,curr_year,snp5_cluster,snp10_cluster,last_move)
    #print (historic)
    previous_snp5 = set(historic.snp5)
    previous_snp10 = set(historic.snp10)
    #print (snp5_cluster,previous_snp5)
    if len(historic) == 0:
        return 1, 'No previous breakdowns found'
    elif snp5_cluster in previous_snp5:
        return 10, 'Current breakdown strain same as previous at snp5 level'
    elif snp10_cluster in previous_snp10:
        return 5, 'Current breakdown strain same as previous at snp10 level'
    else:
        return 0, 'Isolates in past breakdowns present but >10 snps.'

def score_local(animal_row, context):
    """
    Scores the likelihood of local area transmission (0-10) by checking if
    the animal's SNP5 cluster is present in a neighboring herd.
    Args:
        animal_row: The specific animal being scored.
        context4: Herd context using 4km for neighbours
        context10: context using 10km for neighbours
    Returns:
        An integer score (0-10).
    """

    snp5_cluster = animal_row['snp5']
    data = context['data']
    nb = data['neighbour_parcels']
    badgers = data['near_badger']
    last_move = animal_row.last_move

    #print (len(nb),snp5_cluster,len(badgers))
    # Check if any neighboring herds were identified in the context setup
    if len(nb) == 0 and len(nb10) == 0:
        # Cannot attribute to local spread if no neighbors are defined or found
        return 1, "No neighbouring herds at 4 or 10km radius"

    nb_isolates = data['neighbour_isolates']
    nb10_isolates = data['neighbour10_isolates']
    nb_snp5 = list(nb_isolates.snp5)
    nb10_snp5 = list(nb10_isolates.snp5)
    #test data for neighbours one year prior to this breakdown
    nb_lp = data['nb_lesion_positives']
    year = last_move.year
    mask = (nb_lp.columns <= year) & (nb_lp.columns >= year - 1)
    nb_lp = nb_lp.loc[:, mask]
    #print (nb_lp.values.sum())
    #print (len(nb_isolates),len(nb10_isolates))
    #is neighbour contiguous?
    #how far is neighbour exactly?
    #badgers?
    badger = nb_isolates[nb_isolates.Species=='Badger']

    # Check for Strong Evidence: Exact SNP5 Cluster Match
    if len(nb_isolates) == 0 and len(nb10_isolates) == 0:
        # Low Evidence: No Match Found
        # Neighbours exist, but they do not share this specific genetic cluster.
        # We assign a baseline score of 1 to ensure this pathway is considered,
        return 1, f"No sampled herds in {len(nb)} neighbours"
    elif snp5_cluster in nb_snp5:
        # High likelihood: The animal's specific, highly related cluster
        # is in a neighbouring herd with 4km.
        match = nb_isolates[nb_isolates.snp5==snp5_cluster]
        return 10, f"Cluster match within 4km. {len(match.HERD_NO)} herds"
    elif snp5_cluster in nb10_snp5:
        # The animal's specific, highly related cluster
        # is in a neighbouring herd with 10km.
        match = nb10_isolates[nb10_isolates.snp5==snp5_cluster]
        return 5, f"Cluster match within 10km {len(match.HERD_NO)} herds"
    #elif nb_lp.values.sum()>0:
    #    return 3, f"Lesion positive animal within 4km 1 year prior"
    return 0, f"{len(nb_isolates)} neighbouring isolates with no related cluster"

def score_movement(animal_row, herd_context, lpis, lpis_cent, metadata, moves_in=None):
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

    months_prior = 36
    data = herd_context['data']
    sub = data['herd_isolates']
    herd = animal_row.HERD_NO
    tag = animal_row.Animal_ID
    last_move = animal_row.last_move
    snp5_cluster = animal_row['snp5']
    snp10_cluster = animal_row['snp10']
    #check if animal is homebred?
    homebred = animal_row['Homebred']
    data = herd_context['data']

    x = data['snp5_related']
    r5 = x[x.snp5==snp5_cluster]
    #x = data['snp10_related']
    #r10 = x[x.snp10==snp10_cluster]
    target_herds = list(r5.HERD_NO)
    #print (target_herds)

    if len(r5) == 0:
        return 1, 'No related strains outside herd at snp5 or snp10'
    if len(r5) >= 1:
        #get all moves from sources to herd in 1 year prior to breakdown
        try:
            date1 = last_move - relativedelta(months=months_prior)
        except:
            return 1, 'Cannot parse last_move date'
        if moves_in is None:
            moves_in = movement.query_all_herd_moves_in(herd, date1, last_move)
            msp = movement.get_moves_spans(moves_in)
            v = movement.get_movement_vectors(msp, lpis_cent)
            #remove local moves
            #print (v)
            v = v[v.dist_km>=15]

        animal_moves = moves_in[moves_in.tag==tag]
        #print (set(target_herds))
        any_from_strain_herds = v[v.move_from.isin(target_herds)]
        #print (any_from_strain_herds)
        sample_from_strain_herds =  v[(v.move_from.isin(target_herds)) & (v.tag==tag)]

        if len(sample_from_strain_herds)>0:
            #did animal move directly from a strain herd?
            return 10, 'Direct move of sampled animal from herd with same strain'
        elif len(any_from_strain_herds)>0:
            #any direct moves of any animals from herds with related isolates
            #hard to do manually
            return 7, 'Move into this herd from another herd with same strain'
        elif set(animal_moves.move_from).intersection(set(target_herds)):
            #strain seen in intermediate herds of the animal
            return 5, 'Same strain seen in intermediate herd of animal'
        else:
            #try intermediate neighbours - score 7 too high?
            #get isolates local to intermediate herds
            ac = tools.get_animal_context(tag, metadata, lpis, lpis_cent)
            if ac is not None:
                int_nb_strains = ac['metrics']['int_neighbour_strains']
                if snp5_cluster in int_nb_strains:
                    return 7, 'Related strain within 4km of an intermediate herd'
        return 5, 'Related strain in herd >10km away but no move'

def bin_confidence(df):
    """
    Bins pathway results into High, Low, and Inconclusive
    based on the winning score value.
    """
    score_cols = ['W', 'L', 'M', 'R']
    df['max_score'] = df[score_cols].max(axis=1)
    conditions = [
        (df['max_score'] >= 8),                          # High: Clear WGS evidence
        (df['max_score'] >= 2) & (df['max_score'] < 8),  # Low: Circumstantial/Weak evidence
        (df['max_score'] <= 1)                           # Inconclusive/Orphan
    ]
    choices = ['High', 'Low', 'Inconclusive']
    df['confidence'] = np.select(conditions, choices, default='Unknown')
    return df

def run_pathways(herds, *args):
    """Run attribution on a set of herds"""

    results = []
    for herd_no in herds:
        result = run_pathway(herd_no, *args)
        print (result)
        results.append(result)
        print ('-----------------')
    df = pd.concat(results)
    df = bin_confidence(df)
    return df

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
        print ('herd context is null')
        return
    data = hc['data']
    sub = data['herd_isolates']
    results_list = []
    for _, row in sub.iterrows():
        s_within, e_within = score_within_herd(row, hc)
        #print (herd_no, row.Animal_ID, s_within)
        s_local, e_local = score_local(row, hc)
        s_movement, e_movement = score_movement(row, hc, lpis, lpis_cent, gdf)
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

        # Find the maximum score
        max_val = max(scores.values())
        # Identify if any pathways share the maximum
        winners = [k for k, v in scores.items() if v == max_val]
        # Apply conditional logic
        if sum(scores.values())==0:
            pathway = 'E'
        elif max_val <= 1:
            pathway = 'U'
        elif len(winners) > 1:
            # Sort winners alphabetically to ensure consistency
            pathway = "".join(sorted(winners))
        else:
            # Only one clear winner
            pathway = winners[0]

        result = {'sample': row['sample'],
                  'HERD_NO': row.HERD_NO,
                  'Animal_ID': row.Animal_ID}
                  #'max_score':max_val}
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
    return pd.DataFrame(results_list).set_index('sample')

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
