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
from . import core

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(module_path,'data')
BASEDIR = None #deltalake folder


def write_delta_lake(csv_file, outdir):
    """Write delta lake folder from movement csv file.
    We run this once with new csv."""

    keep_columns = ["tag", "herd", "event_date", "event_type", "data_type"]
    (pl.scan_csv(
        csv_file,
        infer_schema_length=10000,
        truncate_ragged_lines=True,
        low_memory=True,
        schema_overrides={
            "tag": pl.Utf8,
            "dob": pl.Utf8,
            "event_date": pl.Utf8,
            "data_type": pl.Utf8
        },
        null_values=[],
        ignore_errors=True
    )
    .select(keep_columns)
    .with_columns([
        pl.col("event_date").str.strptime(pl.Date, "%Y%m%d", strict=False)
    ])
    .with_columns(year = pl.col("event_date").dt.year().alias("year"))
    .filter(pl.col("event_type") != 'p')
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
    dt.optimize.z_order(["herd", "event_date"])
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

def query_tags(tags):
    """All moves for sampled tags"""

    if BASEDIR == None:
        print ('BASEDIR is not set')
        return
    table = pl.scan_delta(BASEDIR)
    if type(tags) is str:
        tags=[tags]
    q = table.filter((pl.col("tag").is_in(tags)))
    res = q.collect(engine='streaming')
    return res.sort(['tag','event_date']).to_pandas()

def query_herd_moves_in(herds, start_date, end_date):
    """Query for herds in time frame.
        Returns: A polars dataframe.
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
            pl.col("herd").is_in(herds),
            pl.col("event_date").is_between(start_date, end_date)
        ])
        .select([
            "tag", "herd", "event_date", "event_type", "data_type","year"
        ])
        .collect(engine='streaming')
    ).to_pandas()

def query_herd_moves_all(herds, start_date, end_date):
    """Query for all animal moves who were in the herds in time frame"""

    if type(herds) is str:
        herds = [herds]
    moves_in = query_herd_moves_in(herds, start_date, end_date)
    res = query_tags(list(moves_in['tag']))
    return res

def convert_moves_to_spans(df):
    """Convert moves for animal to get spans"""

    df['event_date'] = pd.to_datetime(df['event_date']).dt.date
    # 2. Ignore Marts and sort
    temp_df = df[df.event_type != 'Mart'].copy()
    # 3. Calculate Departure Date
    temp_df['departure_date'] = temp_df.groupby('tag')['event_date'].shift(-1)
    temp_df['destination_herd'] = temp_df.groupby('tag')['herd'].shift(-1)
    # 4. Handle "Still Present" animals
    today = date.today()
    temp_df['departure_date'] = temp_df['departure_date'].fillna(today)
    temp_df['days_on_herd'] = (temp_df['departure_date'] - temp_df['event_date'])
    # Remove the FACT (death) rows as they mark the end of a residency
    result = temp_df[temp_df['data_type'] != 'FACT'].copy()
    return result

def get_movement_vectors(m, coords_df):
    """Get movement data locations with end/start positions"""

    msp = convert_moves_to_spans(m)
    df_vectors = msp.merge(
        coords_df, left_on='herd', right_on='SPH_HERD_N', how='left'
    ).rename(columns={'X_COORD': 'start_x', 'Y_COORD': 'start_y'})
    #Join for Destination coordinates (Where they moved next)
    df_vectors = df_vectors.merge(
        coords_df, left_on='destination_herd', right_on='SPH_HERD_N', how='left'
    ).rename(columns={'X_COORD': 'end_x', 'Y_COORD': 'end_y'})
    return df_vectors

def categorize_moves(df_vectors):

    dist_m = np.sqrt(
        (df_vectors['end_x'] - df_vectors['start_x'])**2 +
        (df_vectors['end_y'] - df_vectors['start_y'])**2
    )
    df_vectors['dist_km'] = dist_m / 1000
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
                set(df_vectors['destination_herd'].dropna().unique())
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
        df_vectors['destination_herd'].notnull() &
        df_vectors['herd'].isin(pos) &
        df_vectors['destination_herd'].isin(pos)
    ]
    # Vectorized edge addition: (source, target, weight)
    # We use 'days_on_herd' as the edge weight
    edges = valid_moves[['herd', 'destination_herd', 'days_on_herd']].values
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