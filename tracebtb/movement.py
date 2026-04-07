#!/usr/bin/env python

"""
    tracebtb movement functions with polars, delta-rs
    Created Oct 2024
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

import sys,os,traceback
import glob,platform,shutil
from datetime import date
import time
import math
import pandas as pd
import numpy as np
import polars as pl
from deltalake import DeltaTable

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(module_path,'data')
from .core import config
BASEDIR = config.settings['moves_lake'] #deltalake folder


def write_delta_lake(csv_file, outdir):
    """Write delta lake folder from movement csv file.
        We run this once with new csv."""

    keep_columns = ['tag', 'move_from','move_to','mart','move_date',
                    'data_type','year','dob','birth_herd']

    (pl.scan_csv(
        csv_file,
        infer_schema_length=10000,
        truncate_ragged_lines=True,
        rechunk=False,
        quote_char=None,
        schema_overrides={
            "tag": pl.Utf8,
            "move_date": pl.Utf8,
            "data_type": pl.Categorical,
            "birth_id": pl.Utf8,
            "birth_herd": pl.Utf8,
            "dob": pl.Utf8,
            "year": pl.Int64
        },
        null_values=[],
        ignore_errors=True
    )
    .select(keep_columns)
    .with_columns(pl.col(pl.String, pl.Categorical).cast(pl.String))
    #.with_columns(pl.col(pl.Utf8).str.strip_chars())
    .with_columns([
        pl.col("move_date").str.strptime(pl.Date, "%Y%m%d", strict=False),
        pl.col("dob").str.strptime(pl.Date, "%Y%m%d", strict=False)
    ])
    #.with_columns(year = pl.col("data_year").alias("year"))
    #.filter(pl.col("event_type") != 'p')
    #.sink_parquet(os.path.join(path,"cattle_moves_2023_cleaned.parquet"))
    .sink_delta(outdir,
                mode="overwrite",
                # Partitioning by year allows Polars to skip entire folders
                # when you filter by date.
                delta_write_options={"partition_by": "year"})
    )
    return

def z_ordering(basedir):
    """Perform z ordering on delta lake files"""

    dt = DeltaTable(basedir)
    # Run the Z-Order optimization
    # This reorganizes the files so herd/date searches are lightning fast
    dt.optimize.z_order(["move_to", "move_date"])
    return

def check_delta_lake(basedir):

    dt = DeltaTable(basedir)
    print(f"Approx Row Count: {dt.count()}")
    table = pl.scan_delta(BASEDIR)
    print(table.head(5).collect())
    bad_dates = (
        table
        .filter(pl.col("event_date").is_null())
        .select(pl.len())
        .collect()
    )
    print(f"Rows with unparseable dates: {bad_dates.item()}")
    return

def inject_birth_moves(df):
    """
    Creates a 'Birth' move for each animal.
    """

    df['move_date'] = pd.to_datetime(df['move_date']).dt.date
    df['dob'] = pd.to_datetime(df['dob']).dt.date
    # Extract the 'First Move' for every animal to get birth context
    birth_rows = df.groupby('tag').first().reset_index()
    # set 'move_date' to the 'dob' and 'move_to' to the birth herd
    birth_df = birth_rows.copy()
    birth_df['move_date'] = birth_df['dob']
    birth_df['move_to'] = birth_df['birth_herd']
    birth_df['move_from'] = None  # It didn't come from anywhere else
    birth_df['data_type'] = 'Birth'
    # Combine and Re-sort
    combined_df = pd.concat([birth_df, df], ignore_index=True)
    combined_df = combined_df.sort_values(['tag', 'move_date'])
    return combined_df

def query_tags(tags):
    """
    All moves for sampled tags.
    Returns: A pandas dataframe.
    """

    if BASEDIR == None:
        print ('BASEDIR is not set')
        return
    table = pl.scan_delta(BASEDIR)
    if type(tags) is str:
        tags=[tags]
    q = table.filter((pl.col("tag").is_in(tags)))
    try:
        res = q.collect(engine='streaming')
    except:
        return None
    res = res.sort(['tag','move_date']).to_pandas()
    res = inject_birth_moves(res)
    return res

def query_herd_moves_in(herds, start_date, end_date):
    """Query for direct moves into herds in time frame.
        Returns: A pandas dataframe.
        Usage:
            herd_moves = query_herds(
            ['herdA','herd'B'],
            date(2018, 1, 1),
            date(2022, 12, 30)
            )
    """

    if BASEDIR == None:
        print ('BASEDIR is not set')
        return
    if type(herds) is str:
        herds = [herds]
    return (
        pl.scan_delta(BASEDIR)
        .filter([
            pl.col("move_to").is_in(herds),
            pl.col("move_date").is_between(start_date, end_date)
        ])
        .collect(engine='streaming')
    ).to_pandas()

def query_all_herd_moves_in(herds, start_date, end_date):
    """
    Returns moves into 'herds' and moves of animals in leading up to that arrival.
    """
    if isinstance(herds, str):
        herds = [herds]

    # Get the 'arrival' moves (the moment they entered your study herds)
    arrivals = query_herd_moves_in(herds, start_date, end_date)
    # Create a mapping of Tag -> Latest allowed date (Arrival Date)
    max_dates = arrivals.groupby('tag')['move_date'].max().to_dict()
    full_history = query_tags(list(arrivals['tag'].unique()))
    # remove subsequent moves
    full_history['max_allowed'] = full_history['tag'].map(max_dates)
    filtered_history = full_history[full_history['move_date'] <= full_history['max_allowed']]
    return filtered_history.drop(columns=['max_allowed'])

def query_herd_moves_out(herds, start_date, end_date):
    """Query for herds in time frame.
        Returns: A pandas dataframe.
        Usage:
            herd_moves = query_herds(
            ['herdA','herd'B'],
            date(2018, 1, 1),
            date(2022, 12, 30)
            )
    """

    if BASEDIR == None:
        print ('BASEDIR is not set')
        return
    if type(herds) is str:
        herds = [herds]
    return (
        pl.scan_delta(BASEDIR)
        .filter([
            pl.col("move_from").is_in(herds),
            pl.col("move_date").is_between(start_date, end_date)
        ])
        .collect(engine='streaming')
    ).to_pandas()

def query_herd_moves_all(herds, start_date, end_date):
    """
    Query for herds moves in and out.
    """

    if BASEDIR == None:
        print ('BASEDIR is not set')
        return
    if type(herds) is str:
        herds = [herds]
    return (
        pl.scan_delta(BASEDIR)
        .filter([
            pl.col("move_from").is_in(herds),
            pl.col("move_to").is_in(herds),
            pl.col("move_date").is_between(start_date, end_date)
        ])
        .collect(engine='streaming')
    ).to_pandas()
    return

def get_strain_moves_in(herd, hc, moves_in=None,
                        start=date(2015, 1, 1),
                        end=date(2024, 1, 1)):
    """
    Returns moves from herds with shared strains.
    """
    r5 = hc['data']['snp5_related']
    target_herds = list(r5.HERD_NO)
    if moves_in is None:
        moves_in = query_all_herd_moves_in(herd, start, end)
    msp = get_moves_spans(moves_in)
    #moves from or to target herds
    amoves = msp[(msp['move_from'].isin(target_herds)) | msp['move_to'].isin(target_herds)]
    #get all moves for those animals
    tags = list(amoves.tag)
    moves = msp[msp.tag.isin(tags)]
    return moves

def get_moves_spans(df):
    """Get moves time spans"""

    df['move_date'] = pd.to_datetime(df['move_date'])
    # Sort by tag and date to ensure chronological order
    df = df.sort_values(by=['tag', 'move_date']).reset_index(drop=True)
    df['departure_date'] = df.groupby('tag')['move_date'].shift(-1)
    # Calculate 'days_on_herd' as the difference between departure and arrival
    df['days_on_herd'] = (df['departure_date'] - df['move_date']).dt.days
    return df

def get_movement_vectors(m, lpis_cent):
    """Get movement data locations with end/start positions"""

    coords_df = lpis_cent[['SPH_HERD_N','X_COORD','Y_COORD']]
    #msp = convert_moves_to_spans(m)
    df_vectors = m.merge(
        coords_df, left_on='move_from', right_on='SPH_HERD_N', how='left'
    ).rename(columns={'X_COORD': 'start_x', 'Y_COORD': 'start_y'})
    #Join for Destination coordinates (Where they moved next)
    df_vectors = df_vectors.merge(
        coords_df, left_on='move_to', right_on='SPH_HERD_N', how='left'
    ).rename(columns={'X_COORD': 'end_x', 'Y_COORD': 'end_y'})
    dist_m = np.sqrt(
        (df_vectors['end_x'] - df_vectors['start_x'])**2 +
        (df_vectors['end_y'] - df_vectors['start_y'])**2
    )
    df_vectors['dist_km'] = dist_m / 1000
    df_vectors = df_vectors.drop(columns=['SPH_HERD_N_x','SPH_HERD_N_y'])
    df_vectors['herd'] = df_vectors.move_from
    return df_vectors

def categorize_moves(df_vectors):

    bins = [0, 10, 50, np.inf]
    labels = ['Local', 'Regional', 'Long']
    df_vectors['move_category'] = pd.cut(df_vectors['dist_km'], bins=bins, labels=labels)
    return df_vectors

def create_herd_network(m, herds, lpis_cent):
    """
    Creates a Directed Graph and a position dictionary.
    Only includes herds that have valid X-Y coordinates.
    """

    if type(herds) is not list:
        herds = [herds]
    import networkx as nx
    cols=['SPH_HERD_N','X_COORD','Y_COORD']
    coords_df = lpis_cent[cols]
    df_vectors = get_movement_vectors(m, coords_df)
    df_vectors['link'] = np.where(
        df_vectors['herd'].isin(herds),
        'direct',
        'indirect'
    )
    # 1. Identify all active herds from the movement data
    active_herds = set(df_vectors['herd'].unique()) | \
                set(df_vectors['move_to'].dropna().unique())
    #print (active_herds)
    # 2. Efficiently build 'pos' dict using only herds present in this move set
    mask = coords_df['SPH_HERD_N'].isin(active_herds)
    relevant_coords = coords_df[mask]
    # Create the dictionary: { 'herd_id': (x, y) }
    pos = dict(zip(
        relevant_coords['SPH_HERD_N'],
        zip(relevant_coords['X_COORD'], relevant_coords['Y_COORD'])
    ))
    G = nx.DiGraph()
    # 4. Add edges ONLY if both the source and destination have a position
    valid_moves = df_vectors[
        df_vectors['move_to'].notnull() &
        df_vectors['herd'].isin(pos) &
        df_vectors['move_to'].isin(pos)
    ]
    # Vectorized edge addition: (source, target, weight)
    # We use 'days_on_herd' as the edge weight
    edges = valid_moves[['herd', 'move_to', 'days_on_herd']].values
    G.add_weighted_edges_from(edges)

    node_attr_df = valid_moves[['herd', 'link']].drop_duplicates('herd').set_index('herd')
    # Also include the destination herds that might not be in the 'herd' column
    # but still need a 'relation' attribute
    dest_nodes = [n for n in G.nodes() if n not in node_attr_df.index]
    if dest_nodes:
        # If they aren't in our target list, they are 'related'
        extra_attrs = pd.DataFrame({'link': 'indirect'}, index=dest_nodes)
        node_attr_df = pd.concat([node_attr_df, extra_attrs])

    # Convert dataframe to dict of dicts: {herd_id: {'relation': 'direct'}, ...}
    node_data = node_attr_df.to_dict('index')
    nx.set_node_attributes(G, node_data)
    return G, pos