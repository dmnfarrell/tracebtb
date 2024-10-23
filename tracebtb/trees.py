"""
    Tree methods for bacterial phylogenetics, mostly using ete3.
    Created Nov 2019
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

import sys,os,subprocess,glob,shutil,re,random
import platform
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO
import numpy as np
import pandas as pd
from  . import tools

qcolors = ['blue','green','crimson','blueviolet','orange','cadetblue','chartreuse','chocolate',
            'coral','gold','cornflowerblue','palegreen','khaki','orange','pink','burlywood',
            'red','lime','mediumvioletred','navy','teal','darkblue','purple','orange',
            'salmon','maroon']

def set_tiplabels(t, labelmap):
    for l in t.iter_leaves():
        #print (l.name)
        if l.name in labelmap:
            l.name = labelmap[l.name]
    return

def get_colormap(values):

    import pylab as plt
    labels = values.unique()
    cmap = plt.cm.get_cmap('Set1')
    colors = [cmap(i) for i in range(len(labels))]
    #colors=qcolors
    #clrs = {labels[i]:cmap(float(i)/(len(labels))) for i in range(len(labels))}
    clrs = dict(list(zip(labels,colors)))
    return clrs

def run_fasttree(infile, outpath='', bootstraps=100):
    """Run fasttree on fasta alignment"""

    fc = tools.get_cmd('fasttreeMP')
    out = os.path.join(outpath,'tree.newick')
    cmd = '{fc} -nt {i} > {o}'.format(fc=fc,b=bootstraps,i=infile,o=out)
    try:
        tmp = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except Exception as e:
        print(e)
    return out

def run_RAXML(infile, name='variants', threads=8, bootstraps=100, outpath='.'):
    """Run Raxml pthreads.
        Returns:
            name of .tree file.
    """

    outpath = os.path.abspath(outpath)
    if not os.path.exists(outpath):
        os.makedirs(outpath, exist_ok=True)

    model = 'GTRCAT'
    s1 = random.randint(0,1e8)
    s2 = random.randint(0,1e8)

    files = glob.glob(os.path.join(outpath,'RAxML_*'))
    for f in files:
        os.remove(f)
    if platform.system() == 'Windows':
        cmd = tools.get_cmd('RAxML')
    else:
        cmd = 'raxmlHPC-PTHREADS'
    cmd = '{c} -f a -N {nb} -T {t} -m {m} -V -p {s1} -x {s2} -n {n} -w {w} -s {i}'\
            .format(c=cmd,t=threads,nb=bootstraps,n=name,i=infile,s1=s1,s2=s2,m=model,w=outpath)
    print (cmd)
    try:
        tmp = subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print ('Error building tree. Is RAxML installed?')
        return None
    out = os.path.join(outpath,'RAxML_bipartitions.variants')
    return out

def convert_branch_lengths(treefile, outfile, snps):

    tree = Phylo.read(treefile, "newick")
    for parent in tree.find_clades(terminal=False, order="level"):
            for child in parent.clades:
                if child.branch_length:
                    child.branch_length *= snps
    #Phylo.draw(tree)
    Phylo.write(tree, outfile, "newick")
    return

def tree_from_aln(aln):
    """Make tree from core snp matrix"""

    if len(aln) == 0:
        return
    AlignIO.write(aln, 'temp.fa', 'fasta')
    treefile = run_fasttree('temp.fa')
    ls = len(aln[0])
    convert_branch_lengths(treefile,treefile, ls)
    return treefile

def tree_from_snps(snpmat):
    """Make tree from core snp matrix"""

    aln = tools.alignment_from_snps(snpmat)
    treefile = tree_from_aln(aln)
    return treefile

def njtree_from_snps(df):
    """NJ tree from core SNP alignment"""

    aln = tools.alignment_from_snps(df)
    # Calculate the pairwise distances
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(aln)
    # Build the Neighbor-Joining tree
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(dm)
    # Plot and display the tree
    ax=Phylo.draw(nj_tree)
    return

def tree_from_distmatrix(df, treefile=None):
    """Make nj tree from snp dist matrix"""

    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    #trim zeroes out
    names = list(df.index)
    matrix = np.tril(df.values.tolist())
    M=[]
    for i in range(len(matrix)):
        M.append(list(matrix[i][:i+1]))
    #print (M)
    # Convert the distance matrix to a BioPython DistanceMatrix object
    dm = DistanceMatrix(names, M)
    # Build the tree using the Neighbor-Joining algorithm
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    def remove_inner_labels(node):
        for child in node.clades:
            remove_inner_labels(child)
        if not node.is_terminal():
            node.name = None
    remove_inner_labels(tree.root)
    #Phylo.draw(tree)
    if treefile != None:
        tree = Phylo.write(tree, treefile, 'newick')
    return tree

def biopython_draw_tree(filename):

    from Bio import Phylo
    tree = Phylo.read(filename,'newick')
    Phylo.draw(tree)
    return

def draw_tree(filename, df=None, col=None, cmap=None, tiplabelcol=None, markercol=None,
              width=500, height=500, **kwargs):
    """Draw newick tree with toytree"""

    import toytree
    tre = toytree.tree(filename)
    idx = tre.get_tip_labels()
    if len(idx) > 80:
        tiplabelcol = None
    if df is None:
        tip_colors = None
        node_colors = None
        node_sizes = None
        tip_labels = False
    else:
        if col not in [None, '']:
            #get color from col
            df = df.fillna('')
            labels = df[col].unique()
            if cmap == None:
                cmap = ({c:tools.random_hex_color() if c in labels else 'black' for c in labels})
            else:
                c,cmap = tools.get_color_mapping(df, col, cmap)
            df['color'] = df[col].apply(lambda x: cmap[x])
        #if color column already present will use that or return
        if 'color' not in df.columns:
            return
        df = df.loc[idx]
        tip_colors = list(df.color)
        node_colors = [df.loc[n]['color'] if n in df.index else 'black' for n in tre.get_node_values('name', True, True)]
        node_sizes=[2 if i else 8 for i in tre.get_node_values(None, 1, 0)]

    if tiplabelcol not in [None, '']:
        tip_labels = [df.loc[n][tiplabelcol] if n in df.index else n for n in idx]
    else:
        tip_labels = False
    if markercol == 'Species':
        smap = {'Bovine':'o','Badger':'s','Deer':'x'}
        markers = []
        for n in tre.get_node_values('name', True, True):
            if n in df.index:
                k = df.loc[n][markercol]
                if k in smap:
                    markers.append(smap[k])
                else:
                    markers.append(None)
            else:
                markers.append(None)
    else:
        markers = None
    canvas,axes,mark = tre.draw(scalebar=True,edge_widths=1,height=height,width=width,
                                tip_labels_colors=tip_colors,node_colors=node_colors,
                                tip_labels=tip_labels,node_markers=markers,
                                node_sizes=node_sizes,**kwargs)
    return canvas

def run_treecluster(f, threshold, method='max_clade'):
    """Run treecluster on a newick tree.
       Clustering Method (options: avg_clade, length,
            length_clade, max, max_clade, med_clade, root_dist,
            single_linkage_clade) (default: max_clade)
        see https://github.com/niemasd/TreeCluster
    """

    import io
    cmd = 'TreeCluster.py  -i {f} -t {t} -m {m}'.format(f=f,t=threshold,m=method)
    #print (cmd)
    cl=subprocess.check_output(cmd, shell=True)
    cl=pd.read_csv(io.BytesIO(cl),sep='\t')
    return cl

def get_clusters(tree):
    """Get snp clusters from newick tree using TreeCluster.py"""

    dists = [3,5,7,10,12,20,50,100]
    c=[]
    for d in dists:
        clust = run_treecluster(tree, threshold=d, method='max_clade')
        #print (clust.ClusterNumber.value_counts()[:10])
        clust['d']='snp'+str(d)
        c.append(clust)

    clusts = pd.pivot_table(pd.concat(c),index='SequenceName',columns='d',values='ClusterNumber').reset_index()
    return clusts

def biopython_tree(tree, tip_labels=True, tip_colors=None, color='blue', s=10, ax=None):
    """Plot tree with biopython Phylo"""

    from Bio import Phylo
    if ax==None:
        fig,ax=plt.subplots(1,1,figsize=(8,4))

    if tip_labels==True:
        func = lambda x: x.name if x.is_terminal() else None
    else:
        func = lambda x: None
    Phylo.draw(tree, label_func=func, axes=ax,
            show_confidence=False, do_show=False)#, fontsize={'fontsize':'5'})

    for line in ax.lines:
        line.set_linewidth(.1)
    tip_positions = calculate_tip_positions(tree)
    tip_labels = list(tip_positions.keys())

    # Create a list of colors for the tips
    #tip_colors = [color_mapping.get(label, 'black') for label in tip_labels]

    if tip_colors != None:
        # Color the tips according to the mapping
        for label, color in zip(tip_labels, tip_colors):
            x, y = tip_positions[label]
            ax.scatter(x, y, color=color, s=s)
    else:
        for label in tip_labels:
            x, y = tip_positions[label]
            #print (x,y)
            ax.scatter(x, y,color=color, s=s)
    ax.axis('off')
    return

def phylocanvas_tree(treefile, df=None, col=None):
    """Draw newick tree"""

    with open(treefile, 'r') as file:
        newick_data = file.read()
    #color_mapping = df[col].to_dict()
    #color_json = json.dumps(color_mapping)

    html = """

        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Phylocanvas Tree Example</title>
            <script src="https://unpkg.com/@phylocanvas/phylocanvas.gl@latest/dist/bundle.min.js"></script>

        </head>
        <body>
            <h1> TEST </h1>
            <div id="demo" style="border: 1px solid lightgray"></div>

            <script>
            const phylotree = new phylocanvas.PhylocanvasGL(
            document.querySelector("#demo"),
            {{
                showLabels: true,
                showLeafLabels: true,
                size: {{ width: 400, height: 500 }},
                source: `{n}`,
            }},
            );
            </script>
        </body>
        </html>
    """.format(n=newick_data)
    return html