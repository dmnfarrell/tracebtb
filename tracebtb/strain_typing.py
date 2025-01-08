"""
    Strain typing methods for tracebtb.
    Created Dec 2024
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

import sys,os,subprocess,glob,re,random
import platform
import numpy as np
import pandas as pd
from  . import tools, trees

def generate_strain_names(df):
    """Generate strain names from fastbaps cluster levels
    Args:
        cl (dataframe): cluster level labels
    Returns:
        new strain names in a dataframe
    """
    
    new = []
    #print (df.columns)
    for idx, sample in df.iterrows():
        #short_name = f"IE{sample['Level.1']}"
        strain_name = f"MB{sample['Level.1']}.{sample['Level.2']}.{sample['Level.3']}"       
        df.loc[idx,'IE_clade'] = f"MB{sample['Level.1']}-{sample['Level.2']}"        
        df.loc[idx,'strain_name'] = strain_name
    return df

def identify_clade_defining_snps(snp_df, threshold=0.55, min_size=None):
    """
    Identify clade-defining SNPs from a SNP table.

    Parameters:
        snp_table (pd.DataFrame): A DataFrame where rows are SNP positions and columns are samples.
                                 The first row contains reference nucleotides.
                                 The last column indicates the clade assignment for each sample.
        threshold (float): The minimum frequency of the ALT allele within a clade to consider it clade-defining. 
                            0.55 is the lowest value that seems to be useful before clade assignment will start to fail

    Returns:
        pd.DataFrame: A DataFrame listing clade-defining SNPs with columns for reference, SNP, clade, position, and frequency.
    """
    if min_size != None:
        c=snp_df['clade'].value_counts()
        omit = c[c<=min_size].index
        snp_df = snp_df[~snp_df['clade'].isin(omit)]
    
    # Extract the reference row and clade information
    ref_row = snp_df.iloc[0, :-1]  # First row, excluding clade column
    clade_column = snp_df.iloc[1:, -1]  # Clade information, excluding the reference row
    
    # List to store clade-defining SNPs
    unique_snps = []

    # Iterate through unique clades
    unique_clades = clade_column.unique()
    for clade in unique_clades:
        # Identify samples belonging to the current clade
        clade_samples = snp_df.iloc[1:, :-1][clade_column == clade]
        non_clade_samples = snp_df.iloc[1:, :-1][clade_column != clade]

        # Find SNPs that are nearly fixed in the clade (ALT allele) and REF in others
        for position in clade_samples.columns:
            clade_alleles = clade_samples[position]
            non_clade_alleles = non_clade_samples[position]

            # Calculate allele frequencies
            clade_alt_freq = (clade_alleles != ref_row[position]).mean()
            non_clade_ref_freq = (non_clade_alleles == ref_row[position]).mean()

            # Check if the SNP meets the criteria
            if clade_alt_freq >= threshold and non_clade_ref_freq == 1.0:
                unique_snps.append({
                    "ref": ref_row[position],
                    "snp": clade_alleles.mode()[0],
                    "clade": clade,
                    "pos": position,
                    "freq": clade_alt_freq
                })

    # Convert results to a DataFrame
    res = pd.DataFrame(unique_snps)
    return res

def assign_clade(new_sample, unique_snps_df):
    """
    Assigns a clade to a new sample based on unique SNPs.

    Args:
        new_sample (pd.Series): SNP data for a new sample, with SNP positions as the index.
        unique_snps_df (pd.DataFrame): Contains strictly unique SNPs with their associated clades.

    Returns:
        str: The assigned clade for the new sample.
    """
    # Extract SNP positions, clades, and alleles from the unique SNPs DataFrame
    snp_positions = unique_snps_df["pos"].values
    snp_alleles = unique_snps_df["snp"].values
    clade_labels = unique_snps_df["clade"].values

    # Filter the new sample to include only SNP positions of interest
    relevant_snps = new_sample[new_sample.index.isin(snp_positions)]

    # Match SNPs and count occurrences for each clade
    clade_matches = pd.Series(0, index=unique_snps_df["clade"].unique(), dtype=int)
    for idx, pos in enumerate(snp_positions):
        if pos in relevant_snps.index and relevant_snps[pos] == snp_alleles[idx]:
            clade_matches[clade_labels[idx]] += 1

    # Assign clade based on the highest match count
    assigned_clade = clade_matches.idxmax() if clade_matches.max() > 0 else "Unclassified"
    return assigned_clade

# Function to generate a random reference genome
def generate_reference_genome(length=1000):
    return ''.join(random.choices(['A', 'C', 'G', 'T'], k=length))

def introduce_mutations(genome, mutation_rate):
    genome_list = list(genome)
    for i in range(len(genome_list)):
        if random.random() < mutation_rate:
            current_base = genome_list[i]
            # Create a list of nucleotides excluding the current base
            nucleotides = [nuc for nuc in ['A', 'C', 'G', 'T'] if nuc != current_base]
            genome_list[i] = random.choice(nucleotides)  # Choose a random nucleotide different from the current base
    return ''.join(genome_list)

def mutate_genome(genome, mutation_rate):
    """
    Introduce random mutations into a genome.
    
    Parameters:
    genome (str): The original genome sequence.
    mutation_rate (float): Probability of mutation at each site.
    
    Returns:
    str: A mutated genome sequence.
    """
    genome_list = list(genome)
    for i in range(len(genome_list)):
        if random.random() < mutation_rate:
            genome_list[i] = random.choice([base for base in 'ACGT' if base != genome_list[i]])
    return ''.join(genome_list)
    
def simulate_hierarchical_clades(ref_genome, num_clades, num_subclades, samples_per_subclade, mutation_rate):
    """
    Simulate hierarchical clades with a nested structure.
    
    Parameters:
    ref_genome (str): Reference genome sequence.
    num_clades (int): Number of main clades.
    num_subclades (int): Number of subclades per main clade.
    samples_per_subclade (int): Number of samples per subclade.
    mutation_rate (float): Mutation rate for introducing changes to genomes.
    
    Returns:
    dict: A dictionary where keys are clade/subclade identifiers and values are lists of genome sequences.
    """
    clades = {}
    
    for clade_idx in range(num_clades):
        # Generate a main clade genome by mutating the reference genome
        clade_genome = mutate_genome(ref_genome, mutation_rate)
        clade_name = f"Clade{clade_idx+1}"
        
        for subclade_idx in range(num_subclades):
          
            # Generate a subclade genome by further mutating the clade genome
            # Subclades diverge less than main clades
            subclade_genome = mutate_genome(clade_genome, mutation_rate / 2)  
                       
            # Generate individual genomes for the subclade
            subclade_sequences = [
                # Individual variation within the subclade
                mutate_genome(subclade_genome, mutation_rate / 4) 
                for _ in range(samples_per_subclade)
            ]
            clades[clade_name] = subclade_sequences
    
    return clades
    
def generate_snp_table(clades, reference):
    snp_data = {}
    clade_info = {}
    
    # Iterate through each clade and its samples
    for clade, genomes in clades.items():
        for sample_id, genome in enumerate(genomes):
            sample_name = f"{clade}_S{sample_id + 1}"
            snp_data[sample_name] = list(genome)
            clade_info[sample_name] = clade  # Store the clade for each sample
    
    # Create a DataFrame with SNP data
    snp_df = pd.DataFrame(snp_data, index=[f"{i + 1}" for i in range(len(reference))])
    snp_df.insert(0, "ref", list(reference))  # Add reference genome as the first column
    
    # Add clade information as a new row
    clade_row = pd.Series(clade_info)
    snp_df.loc['clade'] = clade_row    
    return snp_df
    
def snp_table_to_fasta(snp_table, output_fasta="output.fasta"):
    """
    Converts SNP table to FASTA format.

    Parameters:
    snp_table (pd.DataFrame): The SNP table with a reference row and individual genomes as columns.
    output_fasta (str): The output file path for the FASTA file.
    """
    with open(output_fasta, "w") as fasta_file:
        # Iterate over each genome (column) in the SNP table, excluding the reference column
        for sample in snp_table.columns[1:]:  # Skip the first column ("ref")
            fasta_file.write(f">{sample}\n")  # Write the header for each genome
            genome_sequence = ''.join(snp_table[sample])  # Join the SNPs into a full sequence
            fasta_file.write(f"{genome_sequence}\n")  # Write the sequence to the FASTA file

def test_simulated():

    ref_genome = generate_reference_genome(length=50)
    clades = simulate_hierarchical_clades(ref_genome, num_clades=5, num_subclades=5, 
                                          samples_per_subclade=3, mutation_rate=0.1)
    sim_table = generate_snp_table(clades, ref_genome)
    snp_table_to_fasta(sim_table[:-1], "simulated.fa")
    trees.run_fasttree('simulated.fa','simulated.newick')
    sim_df = sim_table.T
    clade_snps = identify_clade_defining_snps(sim_df)
    new = sim_table.T.sample(1).squeeze()
    assign_clade(new, clade_snps)
    test = sim_df.apply(lambda x: assign_clade(x, clade_snps),1)
    test
    return