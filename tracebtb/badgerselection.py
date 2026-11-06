"""
    TraceBTB badger sample selection dashboard.
    Triage of new badger samples by location, before culture/sequencing.
    Created Oct 2026
    Copyright (C) Damien Farrell

    The user supplies a location for a new sample. We look up how well that
    area is already covered by sequence data (adjusted Good's coverage on the
    hex grid, see tools.grid_metrics) and combine that with local context
    (nearby sequenced isolates, badger post mortem records) to give a
    priority score for whether the sample should be cultured/sequenced.

    Layout of this module:
        load_badger_pm_table()   - PLACEHOLDER loader for the badger PM table
        get_location_context()   - gathers grid + local context for a point
        pm_adjustment()          - PLACEHOLDER use of PM data in the score
        compute_priority()       - PLACEHOLDER scoring, swap in the real thing
        BadgerSelectionDashboard - the panel app

    Launch from weblauncher.py like the other dashboards, e.g.
        layout4 = create_bootstrap_layout_cls(
            badgerselection.BadgerSelectionDashboard, 'Badger Selection', '#8E5EA2',
            layers=layers, treefile=treefile, lpis=lpis, settings=settings, **data)
        layouts = {..., '/badger': layout4}
"""

import os
from functools import lru_cache
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.colors as mcolors
from shapely.geometry import Point
from pyproj import Transformer

from bokeh.models import (ColumnDataSource, GeoJSONDataSource, HoverTool, Range1d,
                          ColorBar, LinearColorMapper)
from bokeh.events import Tap
import panel as pn
import panel.widgets as pnw

from tracebtb import core, tools, bokeh_plot
from tracebtb.dashboards import Dashboard, stylesheet, df_stylesheet, get_icon

IRISH_GRID = 'EPSG:29902'
WEB_MERCATOR = 'EPSG:3857'

#coordinate systems the user can enter a location in:
#label: (epsg, x label, y label, step, decimals to show)
INPUT_CRS = {
    'Irish Grid (EPSG:29902)': (IRISH_GRID, 'Easting', 'Northing', 1, 0),
    'ITM (EPSG:2157)': ('EPSG:2157', 'Easting', 'Northing', 1, 0),
    'Lat/Lon (EPSG:4326)': ('EPSG:4326', 'Longitude', 'Latitude', 0.001, 6),
}
DEFAULT_INPUT_CRS = 'Irish Grid (EPSG:29902)'

#priority score thresholds (0-100), PLACEHOLDER values
TIERS = [(66, 'Tier 1: High'), (33, 'Tier 2: Medium'), (0, 'Tier 3: Low')]
#columns of the hex grid (as made by tools.grid_metrics) that we report
GRID_COLS = ['count', 'n_herds', 'spat_disp', #'shannon_diversity',
             'goods_coverage', 'adj_goods_coverage']

# PM table settings. ASSUMPTION: same coordinate columns/CRS as the meta table
PM_X_COL = 'X_COORD'
PM_Y_COL = 'Y_COORD'
PM_CRS = IRISH_GRID
PM_COLUMNS = ['pm_id', 'Year', 'HERD_NO', 'result']

#Responsive layout. Wide screens: widgets | map | results side by side.
#Narrow screens (phones): the three panels stack vertically at full width.
NARROW_PX = 900
NARROW_ROW_CSS = f"""
@media (max-width: {NARROW_PX}px) {{
  :host {{ flex-direction: column !important; flex-wrap: nowrap !important;
           width: 100% !important; max-width: 100% !important; min-width: 0 !important;
           height: auto !important; max-height: none !important; overflow: visible !important; }}
}}
"""
#each stacked panel: full width, natural height, never squeezed by the flex container
NARROW_PANEL_CSS = f"""
@media (max-width: {NARROW_PX}px) {{
  :host {{ width: 100% !important; max-width: 100% !important; min-width: 0 !important;
           height: auto !important; flex: none !important; }}
}}
"""
#the map itself needs an explicit height once stacked (60% of the screen height)
NARROW_MAP_CSS = f"""
@media (max-width: {NARROW_PX}px) {{
  :host {{ width: 100% !important; max-width: 100% !important; min-width: 0 !important;
           height: 60vh !important; min-height: 380px !important; flex: none !important; }}
}}
"""

host_colors = core.host_colors

@lru_cache(maxsize=None)
def _transformer(src, dst):
    return Transformer.from_crs(src, dst, always_xy=True)

def transform_xy(x, y, src, dst):
    """Transform a single coordinate pair between CRSs. Returns (x,y)"""

    if src == dst:
        return x, y
    return _transformer(src, dst).transform(x, y)

def empty_pm_table():
    """Empty badger PM table with the expected columns"""

    return gpd.GeoDataFrame(pd.DataFrame(columns=PM_COLUMNS),
                            geometry=gpd.GeoSeries([], crs=IRISH_GRID), crs=IRISH_GRID)

def load_badger_pm_table(path=None):
    """
    Load the table of previous badger post mortem results (where badgers were
    taken from). Expects a csv with columns:
        SAMPLE, ID, SETT, YEAR, RESULT, LAT, LONG, County (plus others, kept as is)
    LAT/LONG are WGS84 decimal degrees. Returns a GeoDataFrame in Irish Grid
    (EPSG:29902) with standard columns pm_id, Year, HERD_NO (the sett number)
    and result, which the map hover and scoring use. If no file is given an
    empty table is returned so the app still runs.
    Rows without coordinates are dropped. Rows that fall outside the Ireland
    county boundaries (bad coordinates?) are kept but flagged in 'in_ireland'.
    """

    if not path or not os.path.exists(path):
        return empty_pm_table()
    df = pd.read_csv(path, dtype={'SETT': str, 'SAMPLE': str})
    df['LAT'] = pd.to_numeric(df['LAT'], errors='coerce')
    df['LONG'] = pd.to_numeric(df['LONG'], errors='coerce')
    nocoords = df['LAT'].isnull() | df['LONG'].isnull()
    if nocoords.any():
        print(f'badger PM table: dropped {nocoords.sum()} rows with no coordinates')
    df = df[~nocoords].copy()
    df = df.rename(columns={'SAMPLE': 'pm_id', 'YEAR': 'Year', 'SETT': 'HERD_NO',
                            'RESULT': 'result'})
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['LONG'], df['LAT']),
                           crs='EPSG:4326').to_crs(IRISH_GRID)
    gdf['in_ireland'] = gdf.geometry.within(tools.ireland_mask)
    if (~gdf['in_ireland']).any():
        print(f"badger PM table: {(~gdf['in_ireland']).sum()} rows outside Ireland, check coordinates")
    return gdf.reset_index(drop=True)

def _nearby(gdf, pt, radius):
    """Subset of points gdf within radius (m) of pt, and nearest distance in km"""

    if gdf is None or len(gdf) == 0:
        return gdf, np.nan
    d = gdf.geometry.distance(pt)
    nearest = d.min() / 1000 if d.notnull().any() else np.nan
    return gdf[d <= radius], nearest

def get_location_context(pt, meta, grid, pm=None, radius=4000):
    """
    Gather everything we know about a location.
    Args:
        pt: shapely Point in Irish Grid (EPSG:29902)
        meta: geodataframe of sequenced isolates (all species)
        grid: hex grid with metrics (the 'ireland_grid' in the project file)
        pm: badger post mortem table (geodataframe)
        radius: local neighbourhood radius in metres
    Returns: dict
    """

    ctx = {'x': pt.x, 'y': pt.y, 'radius_m': radius,
           'in_ireland': bool(tools.ireland_mask.contains(pt))}

    #grid cell the point falls in
    cell = None
    if grid is not None and len(grid) > 0:
        hit = grid[grid.geometry.contains(pt)]
        if len(hit) > 0:
            cell = hit.iloc[0]
    ctx['grid_id'] = None if cell is None else cell.get('grid_id')
    ctx['County'] = None if cell is None else cell.get('County')
    for col in GRID_COLS:
        val = np.nan if cell is None else cell.get(col, np.nan)
        ctx[col] = float(val) if pd.notna(val) else np.nan

    #sequenced isolates nearby
    badger = meta[meta.Species == 'Badger']
    cattle = meta[meta.Species == 'Bovine']
    nb, ctx['nearest_badger_seq_km'] = _nearby(badger, pt, radius)
    nc, _ = _nearby(cattle, pt, radius)
    ctx['n_badger_seq'] = len(nb)
    ctx['n_cattle_seq'] = len(nc)
    local = pd.concat([nb, nc])
    ctx['n_strains_local'] = local['short_name'].nunique() if 'short_name' in local else np.nan

    #badger PM records nearby
    npm, ctx['nearest_pm_km'] = _nearby(pm, pt, radius)
    ctx['n_pm'] = 0 if npm is None else len(npm)
    return ctx

def pm_adjustment(ctx):
    """
    PLACEHOLDER. How the PM records should change the score (points added to
    the 0-100 score, can be negative). e.g. many PM badgers nearby but few
    sequenced might raise priority. Not decided yet so no effect for now.
    """

    return 0.0

def compute_priority(ctx):
    """
    PLACEHOLDER scoring. Replace with the real priority calculation.
    Currently: priority = 100 * (1 - adjusted Good's coverage of the grid cell),
    i.e. the worse sampled the area the higher the priority. Cells with too
    few isolates have undefined coverage and are treated as unsampled (0).
    Args:
        ctx: dict from get_location_context
    Returns: dict with score (0-100 or None), tier, reasons (list of str)
    """

    if ctx['grid_id'] is None:
        return {'score': None, 'tier': 'N/A',
                'reasons': ['Location is outside the sampling grid (not on land in Ireland?)']}
    reasons = []
    adj = ctx['adj_goods_coverage']
    if pd.isna(adj):
        coverage = 0.0
        reasons.append(f"Grid cell {ctx['grid_id']} has only {int(ctx['count'])} sequenced "
                       "isolates. Coverage is undefined so treated as unsampled.")
    else:
        coverage = float(adj)
        reasons.append(f"Adjusted Good's coverage of grid cell {ctx['grid_id']} is {coverage:.2f}.")
    reasons.append(f"{ctx['n_badger_seq']} sequenced badgers and {ctx['n_cattle_seq']} "
                   f"sequenced cattle within {ctx['radius_m']/1000:g} km.")
    if ctx['n_pm']:
        reasons.append(f"{ctx['n_pm']} badger PM records within {ctx['radius_m']/1000:g} km.")
    score = 100 * (1 - coverage) + pm_adjustment(ctx)
    score = float(np.clip(score, 0, 100))
    tier = next(name for threshold, name in TIERS if score >= threshold)
    return {'score': score, 'tier': tier, 'reasons': reasons}

def _txt(series):
    """Series to strings for hover text, blanks for missing"""

    return series.astype(object).where(series.notna(), '').astype(str).tolist()

def points_source(gdf, cols):
    """ColumnDataSource of x,y (web mercator) plus text columns, for scatter plots"""

    if gdf is None or len(gdf) == 0:
        return ColumnDataSource(data={c: [] for c in ['x', 'y'] + cols})
    g = gdf[gdf.geometry.notnull() & ~gdf.geometry.is_empty]
    m = g.to_crs(WEB_MERCATOR)
    data = {'x': m.geometry.x.tolist(), 'y': m.geometry.y.tolist()}
    for c in cols:
        data[c] = _txt(g[c]) if c in g.columns else [''] * len(g)
    return ColumnDataSource(data=data)

def hover_for(renderer, fields):
    """Hover tool listing (label, column) fields"""

    return HoverTool(renderers=[renderer], tooltips=[(l, f'@{c}') for l, c in fields])

class BadgerSelectionDashboard(Dashboard):
    """
    Triage new badger samples by location using sequence sampling coverage.
    Takes the same project data as the other dashboards (meta, ireland_grid, ..).
    Optional extras:
        badger_pm: badger post mortem table as a geodataframe (or in project data)
        settings['badger_pm_file']: path to PM table, read with load_badger_pm_table
    """

    def __init__(self, **kwargs):
        #the base class calls setup_widgets, so set what that needs first
        settings = kwargs.get('settings') or {}
        self.iregrid = kwargs.get('ireland_grid')
        pm = kwargs.get('badger_pm')
        if pm is None:
            pm = load_badger_pm_table(settings.get('badger_pm'))
        self.pm = pm
        self.results = []
        self.new_point = None
        super(BadgerSelectionDashboard, self).__init__(**kwargs)
        if self.meta.crs is not None and str(self.meta.crs) != IRISH_GRID:
            self.meta = self.meta.to_crs(IRISH_GRID)
        return

    def setup_widgets(self):
        """Create widgets and layout. Widgets stretch to the width of their
        container so the same code works side by side on desktop and stacked
        on a phone (see NARROW_ROW_CSS)."""

        sw = {'sizing_mode': 'stretch_width'}
        self.sample_input = pnw.TextInput(name='Sample ID (optional)', **sw)
        self.crs_select = pnw.Select(name='Coordinate system', options=list(INPUT_CRS),
                                     value=DEFAULT_INPUT_CRS, **sw)
        #widget names can't change, so labels cover all the coordinate systems
        self.x_input = pnw.FloatInput(name='X (Easting / Longitude)', value=None, step=1, **sw)
        self.y_input = pnw.FloatInput(name='Y (Northing / Latitude)', value=None, step=1, **sw)
        #self.radius_input = pnw.IntSlider(name='Local radius (m)', value=4000,
        #                                  start=500, end=10000, step=500, **sw)
        self.score_btn = pnw.Button(name='Score Sample', button_type='primary', height=42, **sw)
        self.clear_btn = pnw.Button(name='Clear', button_type='default', **sw)
        self.provider_input = pnw.Select(name='Map Provider', options=[''] + bokeh_plot.providers,
                                         value='CartoDB Positron', **sw)
        #on a touch screen a tap also means 'look at this point', so allow turning it off
        self.tap_mode = pnw.Checkbox(name='Tap map to set location', value=True, **sw)
        hint = pn.pane.Markdown('Enter a location or **tap the map** to set it. '
                                'Then press *Score Sample*.', **sw)

        self.crs_select.param.watch(self.crs_changed, 'value')
        self.x_input.param.watch(self.update_marker, 'value')
        self.y_input.param.watch(self.update_marker, 'value')
        #self.radius_input.param.watch(self.update_marker, 'value')
        self.provider_input.param.watch(self.set_provider, 'value')
        self.score_btn.on_click(self.score_sample)
        self.clear_btn.on_click(self.clear)

        #map. needs a minimum height or it collapses when the panels are stacked
        self.map_pane = pn.pane.Bokeh(self.build_map(), sizing_mode='stretch_both', min_height=420,
                                      stylesheets=[NARROW_MAP_CSS])
        legend_note = pn.pane.HTML(
            "<small>Grid: adjusted Good's coverage (red = poorly sampled, green = well sampled, "
            "grey = fewer than 5 sequenced isolates, undefined). "
            "Use the Layers buttons to hide/show layers.</small>", **sw)
        map_col = pn.Column(self.map_pane, legend_note, sizing_mode='stretch_both', min_height=460,
                            stylesheets=[NARROW_PANEL_CSS])

        #results
        self.status_pane = pn.pane.Alert('', alert_type='warning', visible=False, **sw)
        self.score_indicator = pn.indicators.Number(
            name='Priority score (0-100)', value=0, format='{value:.0f}', font_size='48pt',
            title_size='13pt', visible=False,
            colors=[(33, '#7f8c8d'), (66, '#e6a100'), (100, '#c0392b')])
        mstyles = {'font-size': '16px', 'margin-left': '10px'}
        #self.tier_pane = pn.pane.Markdown('', styles=mstyles, **sw)
        self.reasons_pane = pn.pane.Markdown('', **sw)
        self.context_table = pn.pane.DataFrame(stylesheets=[df_stylesheet], **sw)
        result = pn.Column(self.status_pane, self.score_indicator, #self.tier_pane,
                           self.reasons_pane, self.context_table, **sw)
        self.history_table = pnw.Tabulator(show_index=False, disabled=True, page_size=50,
                                           stylesheets=[stylesheet], sizing_mode='stretch_both',
                                           min_height=250)
        self.history_table.value = pd.DataFrame()
        self.about_pane = pn.pane.Markdown(self.notes(), **sw)
        right_panel = pn.Tabs(('Result', result), ('History', self.history_table),
                              ('Notes', self.about_pane), width=430, dynamic=True,
                              stylesheets=[NARROW_PANEL_CSS])

        widgets = pn.Column(self.sample_input, self.crs_select, self.x_input, self.y_input,
                            self.score_btn, self.clear_btn, hint, self.tap_mode,
                            self.provider_input, self.layer_toolbar(),
                            width=240, stylesheets=[NARROW_PANEL_CSS])
        app = pn.Row(widgets, map_col, right_panel, sizing_mode='stretch_both',
                     stylesheets=[NARROW_ROW_CSS])
        return app

    def layer_toolbar(self):
        """Row of icon toggles (with hover tooltips) to show/hide map layers"""

        #key: (tooltip, icon, shown at start). Icons are Tabler icon names, chosen to
        #echo the glyph used on the map: hexagon grid, square badgers, triangle PM, circle cattle
        layers = {'grid': ('Coverage grid', get_icon('hexbin'), True),
                  'badger': ('Badgers (sequenced)', get_icon('badger'), False),
                  'pm': ('Badger PM records', get_icon('badger_pm'), False),
                  'cattle': ('Cattle (sequenced)', get_icon('cow'), False)}
        self.layer_toggles = {}
        for key, (label, icon, shown) in layers.items():
    
            t = pnw.ToggleIcon(icon=icon, size='2.4em', active_icon=icon,
                               description=label, value=shown, width=32, height=32)
            t.param.watch(lambda event, key=key: self.toggle_layer(key, event.new), 'value')
            self.layer_toggles[key] = t
            self.toggle_layer(key, shown) #apply the start state to the map
        return pn.Column(pn.pane.Markdown('**Layers**', margin=(5, 10, 0, 10)),
                         pn.Row(*self.layer_toggles.values()), sizing_mode='stretch_width')

    def toggle_layer(self, key, visible):
        """Show or hide the map renderers for a layer"""

        for model in self.layer_renderers.get(key, []):
            model.visible = bool(visible)
        return

    def notes(self):
        """Explanatory text"""

        return ("### How the score works (placeholder)\n"
                "Priority = 100 x (1 - adjusted Good's coverage) of the hex grid cell "
                "containing the location. The Adjusted coverage is Good's coverage x herd count "
                "factor x spatial dispersion. This takes into account the previous sampling"
                "of all sequenced samples (all species) and how well distibuted they are in the cell"
                "Cells with too few isolates are treated as unsampled.\n\n"
                "Badger PM records and nearby isolates are shown as context but do not "
                "change the score yet (`pm_adjustment`).\n\n"
                "Replace `compute_priority` in `badgerselection.py` with the real calculation.")

    def build_map(self):
        """Build the map once. Layers are toggled by the legend and the sample
        marker is updated in place, so the figure and tap handler persist."""

        #renderers that the layer toolbar switches on/off: name -> list of models
        self.layer_renderers = {'grid': [], 'badger': [], 'pm': [], 'cattle': []}
        p = bokeh_plot.init_figure(provider='CartoDB Positron')
        xmin, ymin, xmax, ymax = core.counties_gdf.total_bounds
        p.x_range = Range1d(xmin, xmax)
        p.y_range = Range1d(ymin, ymax)

        #county outlines
        cs = GeoJSONDataSource(geojson=core.counties_gdf.to_json())
        p.patches('xs', 'ys', source=cs, fill_alpha=0, line_color='#555555', line_width=0.8)

        #coverage grid
        grid = self.iregrid
        if grid is not None and len(grid) > 0:
            g = grid.copy()
            cov = g['adj_goods_coverage'] if 'adj_goods_coverage' in g else pd.Series(np.nan, index=g.index)
            cmap = mpl.colormaps['RdYlBu_r']
            g['color'] = [mcolors.rgb2hex(cmap(float(np.clip(v, 0, 1)))) if pd.notna(v) else '#bbbbbb'
                          for v in cov]
            g['adj_txt'] = [f'{v:.2f}' if pd.notna(v) else 'NA' for v in cov]
            gc = g['goods_coverage'] if 'goods_coverage' in g else pd.Series(np.nan, index=g.index)
            g['good_txt'] = [f'{v:.2f}' if pd.notna(v) else 'NA' for v in gc]
            keep = [c for c in ['grid_id', 'County', 'count', 'n_herds', 'adj_txt', 'good_txt', 'color']
                    if c in g.columns]
            gsrc = GeoJSONDataSource(geojson=g[keep + ['geometry']].to_crs(WEB_MERCATOR).to_json())
            r = p.patches('xs', 'ys', source=gsrc, fill_color='color', fill_alpha=0.5,
                          line_color='gray', line_width=0.5, legend_label='Coverage grid')
            p.add_tools(HoverTool(renderers=[r], tooltips=[
                ('Grid cell', '@grid_id'), ('County', '@County'), ('Isolates', '@count'),
                ('Herds/setts', '@n_herds'), ("Good's coverage", '@good_txt'),
                ("Adj. coverage", '@adj_txt')]))
            palette = [mcolors.rgb2hex(cmap(i / 255)) for i in range(256)]
            mapper = LinearColorMapper(palette=palette, low=0, high=1)
            colorbar = ColorBar(color_mapper=mapper, title="Adj. Good's coverage",
                                orientation='horizontal', width=220, height=12,
                                title_standoff=4, label_standoff=6, location='center')
            p.add_layout(colorbar, 'below') #below rather than right: saves width on phones
            self.layer_renderers['grid'] = [r, colorbar]

        #sequenced isolates
        fields = [('Sample', 'sample'), ('Herd/Sett', 'HERD_NO'), ('Year', 'Year'),
                  ('Strain', 'short_name'), ('snp5', 'snp5')]
        cols = [c for _, c in fields]
        meta = self.meta
        csrc = points_source(meta[meta.Species == 'Bovine'], cols)
        rc = p.scatter('x', 'y', source=csrc, marker='circle', size=5, fill_alpha=0.6,
                       fill_color=host_colors['Bovine'], line_color='black', line_width=0.3,
                       legend_label='Cattle (sequenced)')
        rc.visible = False #many points, off by default (matches the toolbar default)
        self.layer_renderers['cattle'] = [rc]
        p.add_tools(hover_for(rc, fields))
        bsrc = points_source(self.badger, cols)
        rb = p.scatter('x', 'y', source=bsrc, marker='square', size=6, fill_alpha=0.8,
                       fill_color=host_colors['Badger'], line_color='black', line_width=1,
                       legend_label='Badger (sequenced)')
        self.layer_renderers['badger'] = [rb]
        p.add_tools(hover_for(rb, fields))

        #badger PM records
        pmfields = [('PM id', 'pm_id'), ('Year', 'Year'), ('Sett', 'HERD_NO'), ('Result', 'result')]
        psrc = points_source(self.pm, [c for _, c in pmfields])
        rp = p.scatter('x', 'y', source=psrc, marker='triangle', size=10, fill_alpha=0.6,
                       fill_color='#7b3fa0', line_color='black', line_width=1,
                       legend_label='Badger PM record')
        self.layer_renderers['pm'] = [rp]
        p.add_tools(hover_for(rp, pmfields))

        #new sample and its radius, updated in place
        self.radius_source = ColumnDataSource(data={'xs': [], 'ys': []})
        p.patches('xs', 'ys', source=self.radius_source, fill_color='gold', fill_alpha=0.08,
                  line_color='black', line_dash='dashed', line_width=1.5)
        self.sample_source = ColumnDataSource(data={'x': [], 'y': []})
        p.scatter('x', 'y', source=self.sample_source, marker='star', size=26, fill_color='gold',
                  line_color='black', line_width=1.5, legend_label='New sample')

        p.legend.location = 'top_left'
        p.legend.click_policy = 'none' #layers are switched from the toolbar instead
        p.legend.background_fill_alpha = 0.8
        p.legend.label_text_font_size = '9pt'
        p.legend.glyph_height = p.legend.glyph_width = 14
        p.legend.spacing = 1
        p.legend.padding = 4
        p.on_event(Tap, self.map_tapped)
        self.fig = p
        return p

    #location input

    def input_xy(self):
        """Current (x, y) in the chosen CRS, or None if incomplete"""

        x, y = self.x_input.value, self.y_input.value
        if x is None or y is None or not np.isfinite(x) or not np.isfinite(y):
            return None
        return x, y

    def get_input_point(self):
        """Entered location as a shapely Point in Irish Grid, or None"""

        xy = self.input_xy()
        if xy is None:
            return None
        epsg = INPUT_CRS[self.crs_select.value][0]
        try:
            x, y = transform_xy(xy[0], xy[1], epsg, IRISH_GRID)
        except Exception:
            return None
        if not (np.isfinite(x) and np.isfinite(y)):
            return None
        return Point(x, y)

    def crs_changed(self, event):
        """Change coordinate system. Re-express any entered location in it"""

        old = INPUT_CRS[event.old][0]
        epsg, step = INPUT_CRS[event.new][0], INPUT_CRS[event.new][3]
        xy = self.input_xy() #values still in the old CRS
        self.x_input.step = self.y_input.step = step
        if xy is not None:
            x, y = transform_xy(xy[0], xy[1], old, epsg)
            self.set_xy(x, y, epsg)
        return

    def set_xy(self, x, y, epsg):
        """Set inputs from a coordinate in epsg, rounded for the current CRS"""

        dec = INPUT_CRS[self.crs_select.value][4]
        self.x_input.value = round(float(x), dec)
        self.y_input.value = round(float(y), dec)
        return

    def map_tapped(self, event):
        """Map click sets the location"""

        if not self.tap_mode.value:
            return
        epsg = INPUT_CRS[self.crs_select.value][0]
        x, y = transform_xy(event.x, event.y, WEB_MERCATOR, epsg)
        self.set_xy(x, y, epsg)
        return

    def update_marker(self, event=None, zoom=False):
        """Show the new sample marker and radius on the map"""

        pt = self.get_input_point()
        if pt is None:
            self.sample_source.data = {'x': [], 'y': []}
            self.radius_source.data = {'xs': [], 'ys': []}
            return
        mx, my = transform_xy(pt.x, pt.y, IRISH_GRID, WEB_MERCATOR)
        self.sample_source.data = {'x': [mx], 'y': [my]}
        #buffer in metres in the grid CRS then reproject for display
        circle = gpd.GeoSeries([pt.buffer(4000)], crs=IRISH_GRID).to_crs(WEB_MERCATOR)
        xs, ys = circle.iloc[0].exterior.coords.xy
        self.radius_source.data = {'xs': [list(xs)], 'ys': [list(ys)]}
        if zoom:
            x0, y0, x1, y1 = circle.iloc[0].bounds
            padx, pady = (x1 - x0) * 0.6, (y1 - y0) * 0.6
            self.fig.x_range.start, self.fig.x_range.end = x0 - padx, x1 + padx
            self.fig.y_range.start, self.fig.y_range.end = y0 - pady, y1 + pady
        return

    #scoring

    def score_sample(self, event=None):
        """Score the entered location and show the result"""

        pt = self.get_input_point()
        if pt is None:
            self.show_status('Enter a valid location (or click the map) first.')
            return
        try:
            ctx = get_location_context(pt, self.meta, self.iregrid, self.pm,
                                       radius=4000)
            res = compute_priority(ctx)
        except Exception as e:
            self.show_status(f'Scoring failed: {type(e).__name__}: {e}')
            raise
        self.update_marker(zoom=True)
        self.show_result(ctx, res)
        return

    def show_status(self, msg=''):
        self.status_pane.object = msg
        self.status_pane.visible = bool(msg)
        return

    def show_result(self, ctx, res):
        """Update result panes and history"""

        msg = '' if ctx['in_ireland'] else 'Location is outside the Ireland county boundaries.'
        self.show_status(msg)
        score = res['score']
        self.score_indicator.visible = score is not None
        if score is not None:
            self.score_indicator.value = round(score)
        #self.tier_pane.object = f"**{res['tier']}**"
        self.reasons_pane.object = '\n'.join(f'* {r}' for r in res['reasons'])
        self.context_table.object = self.context_frame(ctx)

        self.results.append({
            'sample': self.sample_input.value or f'sample {len(self.results)+1}',
            'X (IG)': round(ctx['x']), 'Y (IG)': round(ctx['y']),
            'grid_id': ctx['grid_id'], 'County': ctx['County'],
            'adj_coverage': None if pd.isna(ctx['adj_goods_coverage']) else round(ctx['adj_goods_coverage'], 2),
            'score': None if score is None else round(score)})
        self.history_table.value = pd.DataFrame(self.results)
        return

    def context_frame(self, ctx):
        """Context dict as a display table"""

        def f(v, d=2):
            if v is None or (isinstance(v, float) and np.isnan(v)):
                return 'NA'
            return f'{v:.{d}f}' if isinstance(v, (float, np.floating)) else str(v)
        r = f"{ctx['radius_m']/1000:g} km"
        rows = {
            'Grid cell': f(ctx['grid_id']), 'County': f(ctx['County']),
            'Sequenced isolates in cell': f(ctx['count'], 0),
            'Herds/setts in cell': f(ctx['n_herds'], 0),
            "Good's coverage": f(ctx['goods_coverage']),
            'Spatial dispersion': f(ctx['spat_disp']),
            "Adj. Good's coverage": f(ctx['adj_goods_coverage']),
            #'Shannon diversity': f(ctx['shannon_diversity']),
            f'Sequenced badgers within {r}': f(ctx['n_badger_seq']),
            'Nearest sequenced badger (km)': f(ctx['nearest_badger_seq_km'], 1),
            f'Sequenced cattle within {r}': f(ctx['n_cattle_seq']),
            f'Distinct strains within {r}': f(ctx['n_strains_local']),
            f'Badger PM records within {r}': f(ctx['n_pm']),
            'Nearest PM record (km)': f(ctx['nearest_pm_km'], 1),
        }
        return pd.DataFrame(pd.Series(rows), columns=['value'])

    def clear(self, event=None):
        """Clear entered location and results (history is kept)"""

        self.x_input.value = None
        self.y_input.value = None
        self.sample_input.value = ''
        self.score_indicator.visible = False
        #self.tier_pane.object = ''
        self.reasons_pane.object = ''
        self.context_table.object = None
        self.show_status('')
        return