#!/usr/bin/env python

"""
    tracebtb movement functions
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
import pickle
import time
import math
import pandas as pd
import numpy as np
import polars as pl
from . import core

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(module_path,'data')
BASEFILE = None

def get_default_query():
    """    
    This is a base out of memory query for the large movement table using polars.
    Set BASEFILE before running this.
    """
    if BASEFILE is None:
        print ('set BASEFILE attribute of module')
        return
    default_query = (
        pl.scan_csv(
            os.path.join(BASEFILE),
            infer_schema_length=1000,
            truncate_ragged_lines=True,
            schema_overrides={
                "tag": pl.Utf8,
                "dob": pl.Utf8,
                "event_date": pl.Utf8,
                "data_type": pl.Categorical
            },
            null_values=[],
            ignore_errors=True
        )
        .with_columns([
            pl.col("event_date").str.strptime(pl.Date, "%Y%m%d", strict=False),
            pl.col("dob").str.strptime(pl.Date, "%Y%m%d", strict=False)
        ])
        .with_columns(year = pl.col("event_date").dt.year().alias("year"))
        .filter(pl.col("event_type") != 'p')
    )
    return default_query

def get_moves_with_tags(tags):

    q = get_default_query().filter((pl.col("tag").is_in(tags)))
    print (q.explain(streaming=True))
    res = q.collect(streaming=True)
    return res

def get_related_moves(herds):
    """Get all related moves for a herd as follows:
        1. all moves for each herd
        2. moves for all the animals in step 1
    """

    default = get_default_query()
    q1 = default.filter((pl.col("herd").is_in(herds)))
    print (q1.explain(streaming=True))
    res1 = q1.collect(streaming=True)
    found_tags = res1.select(pl.col("tag").unique())
    q2 = default.filter(pl.col("tag").is_in(found_tags))
    res2 = q2.collect(streaming=True)
    return res2