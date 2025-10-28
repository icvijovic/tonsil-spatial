import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from collections import Counter
from scipy.optimize import curve_fit
from scipy.stats import poisson
from scipy.special import zeta
import seaborn as sb
import matplotlib.cm as cm
from scipy.stats import binom
import copy
from Bio import Phylo
from ete3 import Tree, TreeStyle, NodeStyle
import os
import warnings
from scipy.interpolate import interp1d
from scipy.optimize import fsolve



# Error bar adjustment for log plots
def log_error_bar(y, y_err):
    return [y*(1- np.exp(-y_err/y)), y*(np.exp(y_err/y) -1)]

def survival_func(data, xvals=None):
    if xvals is None:
        xvals = np.sort(data)
    yvals = np.zeros(len(xvals))
    for i in range(len(xvals)):
        yvals[i] = np.sum(data >= xvals[i]) / len(data)
    return xvals, yvals

# Divergence (1-ANI) between two sequences, or between one sequence and root
def pairwise_div(row1, row2 = None):
    if row2 is None:
        muts1 = set() if isinstance(row1['v_mutations'], float) else row1['v_mutations'].split(',')
        return len(muts1)/len(row1['v_sequence_no_trunc'])

    tot_length = (len(row1['v_sequence_no_trunc']) + len(row2['v_sequence_no_trunc']))/2

    muts1 = set() if isinstance(row1['v_mutations'], float) else row1['v_mutations'].split(',')
    muts2 = set() if isinstance(row2['v_mutations'], float) else row2['v_mutations'].split(',')
    shared = len(set(muts1) & set(muts2))

    return (len(muts1) + len(muts2) - 2*shared) / tot_length

# Calculate empirical PDF of data
def pdf_histogram(data, edges):

    x = np.zeros(len(edges)-1)
    y = np.zeros(len(x))
    y_err = np.zeros(len(x))

    for i in range(len(x)):
        mask = np.logical_and(data >= edges[i], data < edges[i+1])
        if sum(mask) == 0: continue

        x[i] = np.mean(data[mask])
        y[i] = np.sum(mask)/( len(data) * (edges[i+1]-edges[i]) )
        y_err[i] = np.sqrt(np.sum(mask))/( len(data) * (edges[i+1]-edges[i]) )

    return x[y > 0], y[y > 0], y_err[y > 0]

def geometric_pdf(x, p, logY = True):
    ans = (x-1)*np.log(p) + np.log(1-p) if logY else (1-p) * p**(x-1)
    return ans

def expon_pdf(x, X0, logY = True):
    ans = -x/X0 - np.log(X0) if logY else np.exp(-x/X0)/X0
    return ans

# Integer power law PDF conditioned on x > N
def power_law_cond_pdf(x, a, N, logY = True, cond_adj = True):
    ans = -a * np.log(x) - np.log(zeta(a)) if logY else x**(-a) / zeta(a)
    adj = 1
    if cond_adj:
        for n in range(1,N+1): adj -= power_law_cond_pdf(n, a, 0, False)
    ans = ans - np.log(adj) if logY else ans / adj
    return ans

def power_law_combined_pdf(x, a, logY = True, cond_adj = True):
    ans = np.zeros(len(x))
    len_x = int(len(x)/3) + 1
    ans[:len_x] = power_law_cond_pdf(x[:len_x], a, 0, logY, cond_adj)
    ans[len_x:2*len_x-1] = power_law_cond_pdf(x[len_x:2*len_x-1], a, 1, logY, cond_adj)
    ans[2*len_x-1:] = power_law_cond_pdf(x[2*len_x-1:], a, 2, logY, cond_adj)
    return ans

def single_power_law_pdf(x, a, logY = True):
    if logY:
        return -a * np.log(x) - np.log(zeta(a))
    else:
        return x**(-a) / zeta(a)

# Euclidean distance in the XY plane
def euclidean_distance(row1, row2):
    sq_dist = (float(row1['x'])-float(row2['x']))**2 + (float(row1['y'])-float(row2['y']))**2
    return np.sqrt(sq_dist)

def pairwise_helper(mut1, mut2, skip_same = True, tot_len = 0):
    muts1 = set() if isinstance(mut1, float) else mut1.split(',')
    muts2 = set() if isinstance(mut2, float) else mut2.split(',')
    shared = len(set(muts1) & set(muts2))

    if skip_same and shared == len(muts1) and shared == len(muts2):
        return -1, -1, -1
    
    if tot_len == 0:
        return shared, len(muts1) - shared, len(muts2) - shared
    else:
        return shared/tot_len, (len(muts1) - shared)/tot_len, (len(muts2) - shared)/tot_len

def binom_diff_prediction(n_vals, give_prop = False):
    unique_n_vals, counts = np.unique(n_vals, return_counts=True)

    if give_prop:
        edges = np.linspace(0,1,21)
        x_vals = edges[:-1]
        ans = np.zeros(len(x_vals))
        for n_ind in range(len(unique_n_vals)):
            n = unique_n_vals[n_ind]
            for i in range(0,n+1):
                index = len(edges) - 2 if i == n else np.argmax(edges > i/n) - 1
                prob = binom.pmf(i, n, 0.5)
                #x_vals[index] += counts[n_ind]*(i/n)*prob/len(n_vals)
                ans[index] += counts[n_ind]*prob/len(n_vals)
    else:
        x_vals = np.arange(0, max(n_vals)+1)
        ans = np.zeros(len(x_vals))
        for n_ind in range(len(unique_n_vals)):
            n = unique_n_vals[n_ind]
            for i, x in enumerate(x_vals):
                if (n+x) % 2 == 0 and x <= n: ans[i] += counts[n_ind]*binom.pmf(int((x+n)/2), n, 0.5)/len(n_vals)
                if (n-x) % 2 == 0 and x <= n: ans[i] += counts[n_ind]*binom.pmf(int((n-x)/2), n, 0.5)/len(n_vals)
    return x_vals, ans

def get_internal_ypos(clade, terminal_ypos):
    if clade.is_terminal(): return terminal_ypos[clade]
    else: return np.mean([get_internal_ypos(child, terminal_ypos) for child in clade.clades])

def visualize_tree(tree, df_IGH, lineage_id, ignore_EF = False, use_internal_nodes = False, savefig = 0):
    foll_values = df_IGH[df_IGH['lineage_id'] == lineage_id]['follicle'].unique()

    foll_colors = {foll: np.minimum(1, np.random.rand(3,) + 0.1) for foll in foll_values}
    foll_colors['EF'] = np.array([0, 0, 0])

    clade_depths = tree.depths()
    y_positions = {clade: j for j, clade in enumerate(tree.get_terminals())}


    Phylo.draw(tree, do_show=False, label_func=lambda x: None, branch_labels=lambda c: None)

    # Plot points at the tips (leaves) of the tree
    for clade in tree.get_terminals():
        if clade.name.startswith("ROOT"): continue

        muts = df_IGH.loc[int(clade.name[6:])]['v_mutations']
        matching_rows = df_IGH[(df_IGH['lineage_id'] == lineage_id) & (df_IGH['v_mutations'] == muts)]

        inc = 0
        for i, foll in enumerate(matching_rows['follicle'].sort_values()):
            if ignore_EF and foll == 'EF': continue
            plt.scatter(clade_depths[clade]+0.001*inc, y_positions[clade]+1, color=foll_colors[foll],s=30)
            inc += 1
    
    if use_internal_nodes:
        for clade in tree.get_nonterminals():
            if ignore_EF and clade.metadata['foll'] == 'EF': continue
            x = clade_depths[clade]
            y = get_internal_ypos(clade, y_positions)
            plt.scatter(x, y+1, color=foll_colors[clade.metadata['foll']], s=30)

    plt.xlabel("Divergence")
    plt.ylabel(None)
    plt.yticks([])
    #plt.xlim([0,0.065])

    if savefig > 0: plt.savefig(r'Figures\Tree_' + str(savefig) + '.pdf', format='pdf')

    plt.show()

def flag_insertions_deletions(row):
    seq = list(row['v_sequence_no_trunc'])
    if not isinstance(row['v_mutations'], str): return False
    
    mut_list = row['v_mutations'].split(',')
    for mut in mut_list:
        if mut[-4] != ':': return True
        
    return False

def get_root_seq(row, print_info = False):
    seq = list(row['v_sequence_no_trunc'])

    if not isinstance(row['v_mutations'], str): return ''.join(seq)
    
    mut_list = row['v_mutations'].split(',')
    for mut in mut_list:
        if mut[-4] != ':':
            print(mut)
            print(seq)
            continue
        ind = int(mut[:-4]) - 1
        if seq[ind] != 'N': seq[ind] = mut[-3]

    if print_info:
        print("Inferring germline from mutations")
        print("Original:", row['v_sequence_no_trunc'])
        print("Mutations:", row['v_mutations'])
        print("Germline:", ''.join(seq))
        
    return ''.join(seq)

def get_muts_from_consensus_seq(seq,consensus, print_info = False):
    muts = ''
    diff_flag = False

    for i in range(len(seq)):
        if seq[i] != consensus[i] and seq[i] != 'N' and consensus[i] != 'N':
            diff_flag = True
            muts += str(i+1) + ':' + consensus[i] + '>' + seq[i] + ','

    if not diff_flag: return float('nan')

    muts = muts[:-1]
    if print_info:
        print("Inferring mutations from germline")
        print("Original:", seq)
        print("Germline:", consensus)
        print("Mutations:", muts)
    return muts

def survival_func(data, xvals = None):
    if xvals is None: xvals = np.sort(data)
    yvals = np.zeros(len(xvals))
    for i in range(len(xvals)):
        yvals[i] = np.sum(data >= xvals[i]) / len(data)
    return xvals, yvals

# Labels the terminals of the tree with follicular location
# Ignores extrafollicular reads if ignore_EF is True
# If prioritize_IF is true, only labels as EF if sequence has no IF reads
def label_tips(tree, df_IGH, lineage_id, ignore_EF = False, prioritize_IF = True):
    for clade in tree.get_terminals():
        
        muts = df_IGH.loc[int(clade.name[6:])]['v_mutations']
        matching_rows = df_IGH[(df_IGH['lineage_id'] == lineage_id) & (df_IGH['v_mutations'] == muts)]

        if not isinstance(muts, str):
            matching_rows = df_IGH[(df_IGH['lineage_id'] == lineage_id) & [not isinstance(mut, str) for mut in df_IGH['v_mutations']]]

        foll_dict = {}
        for i, foll in enumerate(matching_rows['follicle'].sort_values()):
            if ignore_EF and foll == 'EF': continue
            if foll not in foll_dict: foll_dict[foll] = 0
            foll_dict[foll] += 1

        EF_reads = foll_dict['EF'] if 'EF' in foll_dict else 0
        if prioritize_IF and 'EF' in foll_dict: foll_dict['EF'] = 0

        max_key = []
        max_val = 0
        for key, val in foll_dict.items():
            if val > max_val:
                max_key = [key]
                max_val = val
            elif val == max_val: max_key.append(key)

        if prioritize_IF:
            if max_val == 0 and EF_reads > 0:
                max_key = ['EF']
                max_val = EF_reads
            if 'EF' in foll_dict: foll_dict['EF'] = EF_reads

        clade.metadata = {'foll': max_key[np.random.randint(len(max_key))]}
        if len(foll_dict) > 1: clade.metadata['foll_dict'] = foll_dict

    for clade in tree.get_nonterminals():
        clade.metadata = {}
        
# Labels internal nodes of the tree based on maximum parsimony
# track_dominant sets the algorithm to resolve arbitrary states to the most common state in the clade
# disallow_EF_to_IF prevents the algorithm from inferring EF to IF transitions
def fill_internal_nodes(tree, lineage_id = None, df_IGH = None, track_dominant = True, disallow_EF_to_IF = True):

    # Function to find the follicle in foll_set that is most common in the clade
    def find_majority(clade, foll_set):
        foll_list = list(foll_set)
        ct = np.zeros(len(foll_set))
        freq = np.zeros(len(foll_set))

        for term in clade.get_terminals():
            muts = df_IGH.loc[int(term.name[6:])]['v_mutations']
            matching_rows = df_IGH[(df_IGH['lineage_id'] == lineage_id) & (df_IGH['v_mutations'] == muts)]

            if not isinstance(muts, str):
                matching_rows = df_IGH[(df_IGH['lineage_id'] == lineage_id) & [not isinstance(mut, str) for mut in df_IGH['v_mutations']]]
            
            for _, row in matching_rows.iterrows():
                if row['follicle'] in foll_list:
                    ct[foll_list.index(row['follicle'])] += 1
                    freq[foll_list.index(row['follicle'])] += 1/sum(df_IGH['follicle'] == row['follicle'])

        mask = ct == max(ct)
        freq[np.logical_not(mask)] = -1
        mask = np.logical_and(mask, freq == max(freq))
        choice = np.random.choice(np.array(foll_list)[mask])
        return set([choice])

    # Bottom-up step of maximum parsimony, filling in each state with possibilities
    def fitch_bottom_up(clade, disallow_EF_to_IF):
        if clade.is_terminal(): return clade.metadata['foll']
        else:
            child_states = [fitch_bottom_up(child, disallow_EF_to_IF) for child in clade.clades]
            union = set.union(*child_states)
            if 'EF' in union and len(union) > 1 and disallow_EF_to_IF:
                child_states = [state - {'EF'} if 'EF' in state else state for state in child_states]
            intersection = set.intersection(*child_states)
            if intersection: clade.metadata['foll'] = intersection
            else: clade.metadata['foll'] = set.union(*child_states)
            return clade.metadata['foll']
        
    # Top-down step of maximum parsimony, choosing from possibilities
    def fitch_top_down(clade, track_dom, parent_state=None):
        if parent_state and (clade.metadata['foll'] & parent_state): # If parent state is resolved, choose that
            clade.metadata['foll'] = clade.metadata['foll'] & parent_state
        else:
            if track_dom and len(list(clade.metadata['foll'])) > 1: # Otherwise, choose from possibilities randomly or using find_majority
                clade.metadata['foll'] = find_majority(clade, clade.metadata['foll'])
            clade.metadata['foll'] = set([np.random.choice(list(clade.metadata['foll']))])

        for child in clade.clades:
            fitch_top_down(child, track_dom, clade.metadata['foll'])

    for clade in tree.get_terminals():
        clade.metadata['foll'] = set([clade.metadata['foll']])


    fitch_bottom_up(tree.root, disallow_EF_to_IF)
    fitch_top_down(tree.root, track_dominant)

    for clade in tree.find_clades():
        clade.metadata['foll'] = next(iter(clade.metadata['foll']))


def get_branch_annotations(tree):
    branches = []
    for clade in tree.find_clades(order='level'):
        if clade.is_terminal():
            continue
        for child in clade.clades:
            branches.append({
                'branch_length': child.branch_length,
                'parent_foll': clade.metadata.get('foll', None),
                'child_foll': child.metadata.get('foll', None),
                'parent_depth': tree.depths().get(clade, None)
            })
    return branches

# Calculate phylogenetic diversity below clade (total branch length), with option to stop at migration events
def phylo_diversity(clade, stop_at_migration = True, stop_at_EF = True, start_from_div = 0, tree = None):
    total_length = 0
    if clade.is_terminal(): return 0
    for child in clade.clades:
        if (child.metadata.get('foll', None) == clade.metadata.get('foll', None)) or (not stop_at_migration):
            if child.metadata.get('foll', None) != 'EF' or (not stop_at_EF):
                if start_from_div > 0:
                    if tree.depths().get(clade) >= start_from_div:
                        total_length += child.branch_length + phylo_diversity(child, stop_at_migration, stop_at_EF, start_from_div, tree)
                    else:
                        total_length += phylo_diversity(child, stop_at_migration, stop_at_EF, start_from_div, tree)
                else:
                    total_length += child.branch_length + phylo_diversity(child, stop_at_migration, stop_at_EF, start_from_div, tree)
    return total_length

def n_descendants(clade, stop_at_migration = True):
    total_descendants = 0
    if clade.is_terminal(): return 0
    for child in clade.clades:
        if (child.metadata.get('foll', None) == clade.metadata.get('foll', None)) or (not stop_at_migration):
            total_descendants += (child.branch_length > 0) + n_descendants(child)
    return total_descendants

def total_depth(clade, stop_at_migration = True):
    depth = 0
    if clade.is_terminal(): return 0
    for child in clade.clades:
        if (child.metadata.get('foll', None) == clade.metadata.get('foll', None)) or (not stop_at_migration):
            depth = max(depth, child.branch_length + total_depth(child))
    return depth

def find_first_entry(tree, foll):
    if tree.root.metadata['foll'] == foll: return 0

    depth_list = []
    for clade in tree.find_clades(order='level'):
        if clade.is_terminal(): continue
        for child in clade.clades:
            if child.metadata['foll'] == foll and clade.metadata['foll'] != foll:
                depth_list.append((tree.depths()[clade] + tree.depths()[child]) / 2)
    if len(depth_list) == 0: return float('nan')
    return min(depth_list)

# Get statistics for subtrees, starting at migration events
# Ignore migrations that occur before a certain depth, div_cutoff
def get_subtree_data(tree, cutoff = 0, div_cutoff = 0):
    subtree_divs = []
    subtree_descendants = []
    subtree_depths = []

    for clade in tree.find_clades(order='level'):
        if clade.is_terminal(): continue
        for child in clade.clades:
            # Check for branches where migration has occurred
            if clade.metadata.get('foll', None) != child.metadata.get('foll', None) and tree.depths().get(clade, None) >= div_cutoff:
                child_terminals = child.get_terminals()
                foll_matches = [term.metadata.get('foll',None) == child.metadata.get('foll',None) for term in child_terminals]
                if np.mean(foll_matches) >= cutoff:
                    subtree_divs.append(phylo_diversity(child))
                    subtree_descendants.append(n_descendants(child))
                    subtree_depths.append(total_depth(child))

    return subtree_divs, subtree_descendants, subtree_depths

# Generate a tree with the same structure as underlying tree, but with random migration locations
def randomize_migration_locations(input_tree, stop_at_EF = True, div_cutoff = 0):

    # Step down tree, incrementing running_branch_length over each step until all migrations are assigned
    def step_down_tree(clade, running_branch_length, mig_locs, stop_at_EF, div_cutoff, tree):
        if clade.is_terminal(): return running_branch_length, mig_locs
        for child in clade.clades:
            if stop_at_EF and child.metadata['foll'] == 'EF': continue

            if tree.depths().get(clade) >= div_cutoff:
                running_branch_length += child.branch_length
                if len(mig_locs) > 0:
                    if running_branch_length >= mig_locs[0]:
                        child.metadata['foll'] = clade.metadata['foll'] + 1
                        mig_locs = mig_locs[1:]
                    else: child.metadata['foll'] = clade.metadata['foll']
                else: child.metadata['foll'] = clade.metadata['foll']
                running_branch_length, mig_locs = step_down_tree(child, running_branch_length, mig_locs, stop_at_EF, div_cutoff, tree)
            else:
                child.metadata['foll'] = clade.metadata['foll']
                running_branch_length, mig_locs = step_down_tree(child, running_branch_length, mig_locs, stop_at_EF, div_cutoff, tree)

        return running_branch_length, mig_locs

    tree = copy.deepcopy(input_tree)

    n_migrations = 0
    total_div = phylo_diversity(tree.root, stop_at_migration = False, stop_at_EF = True, start_from_div = div_cutoff, tree = tree)

    for clade in tree.find_clades():
        if clade.is_terminal(): continue
        for child in clade.clades:
            if clade.metadata.get('foll', None) != child.metadata.get('foll', None) and child.metadata.get('foll', None) != 'EF':
                if tree.depths().get(clade) >= div_cutoff: n_migrations += 1

    # Generate locations where migration should occur
    migration_locations = np.random.rand(n_migrations) * total_div
    migration_locations.sort()

    if len(migration_locations) == 0: return tree

    tree.root.metadata['foll'] = 0
    bl, ml = step_down_tree(tree.root, 0, migration_locations, stop_at_EF, div_cutoff, tree)

    #for clade in tree.find_clades(order='level'):
    #    if clade.is_terminal(): continue
    #    for child in clade.clades:
    #        running_branch_length += child.branch_length
    #        if len(migration_locations) > 0:
    #            if running_branch_length >= migration_locations[0]:
    #                child.metadata['foll'] = clade.metadata['foll'] + 1
    #                migration_locations = migration_locations[1:]
    #            else: child.metadata['foll'] = clade.metadata['foll']
    #        else: child.metadata['foll'] = clade.metadata['foll']

    return tree

def get_ef_transition_info(tree):
    results = []

    def traverse(clade, branch_len, tot_len, current_migrations, current_foll, first_migration_time, has_EF_sibling):

        dict_flag = True if clade.metadata.get('foll_dict') is not None else False
        if dict_flag:
            if 'EF' not in clade.metadata.get('foll_dict'): dict_flag = False

        if clade.metadata.get('foll') == 'EF':
            if has_EF_sibling == -1: has_EF_sibling = 0
            results.append((1,branch_len, tot_len, current_migrations, first_migration_time, has_EF_sibling))
        elif dict_flag:
            if has_EF_sibling == -1: has_EF_sibling = 0
            results.append((1,branch_len, tot_len, current_migrations, first_migration_time, has_EF_sibling))
            results.append((0,branch_len, tot_len, current_migrations, first_migration_time, 1))
        elif clade.is_terminal():
            if has_EF_sibling == -1: has_EF_sibling = 1
            results.append((0,branch_len, tot_len, current_migrations, first_migration_time, has_EF_sibling))
        else:
            #results.append((0,branch_len, tot_len, current_migrations))
            EF_flag = False
            for child in clade.clades:
                if child.metadata.get('foll') == 'EF':
                    EF_flag = True
                    break
            for child in clade.clades:
                new_length = tot_len + child.branch_length
                is_migrant = (1 if child.metadata.get('foll') != current_foll else 0) - (1 if child.metadata.get('foll') == 'EF' else 0)
                new_migrations = current_migrations + is_migrant
                new_migration_time = new_length if is_migrant == 1 and first_migration_time == -1 else first_migration_time

                if EF_flag and has_EF_sibling == 0: has_EF_sibling = -1
                elif has_EF_sibling == -1: has_EF_sibling = 1

                traverse(child, child.branch_length, new_length, new_migrations, child.metadata.get('foll'), new_migration_time, has_EF_sibling)

    traverse(tree.root, 0, 0, 0, tree.root.metadata.get('foll'), -1, 0)
    return results

# Identifies lineages with at least two distinct IF sequences, and prepares sequences for tree analysis
def prepare_sequences_for_trees(df_IGH, df_IGH_IF, df_combined, use_EF=False, write_files=False):
    min_lineage_size = 2
    lineage_lens = []
    lineage_ids = []
    N_lineages_tot = 0
    N_trees = 0

    # Return labels of lineages used in tree analysis
    for lineage_id, group in df_IGH.groupby('lineage_id'):
        N_lineages_tot += 1
        group_trimmed = group.drop_duplicates(subset='v_mutations')
        group_trimmed = group_trimmed[group_trimmed['follicle'] != 'EF']
        if len(group_trimmed) >= min_lineage_size:
            lineage_lens.append(len(group_trimmed))
            lineage_ids.append(lineage_id)
            N_trees += 1

    lineage_ids = np.array(lineage_ids)
    lineage_lens = np.array(lineage_lens)

    sorted_indices = np.flip(np.argsort(lineage_lens))
    lineage_lens = lineage_lens[sorted_indices]
    lineage_ids = lineage_ids[sorted_indices]

    # Write sequences to files for each lineage for use in FastTree
    if write_files:
        if not use_EF:
            output_dir = "Data\Lineage_Seqs"
            os.makedirs(output_dir, exist_ok=True)

            for lineage_id, group in df_IGH_IF.groupby('lineage_id'):
                root_flag = False
                group_trimmed = group.drop_duplicates(subset='v_mutations')
                if len(group_trimmed) >= min_lineage_size:
                    filename = os.path.join(output_dir, f"lin_{lineage_id}_no_ef.fasta")
                    with open(filename, "w") as fasta_file:
                        seq_len = None
                        for idx, row in group_trimmed.iterrows():
                            if not isinstance(row['v_mutations'], str):
                                root_flag = True
                                header = f">ROOT__{idx}"
                            else:
                                header = f">CLONE_{idx}"
                            sequence = row['v_sequence_no_trunc']
                            if seq_len is None:
                                seq_len = len(sequence)
                            elif seq_len != len(sequence):
                                print(lineage_id)
                            fasta_file.write(f"{header}\n{sequence}\n")

                        if not root_flag:
                            root_seq = get_root_seq(group_trimmed.iloc[0])
                            fasta_file.write(f">ROOT__INFERRED\n{root_seq}\n")

        else:
            output_dir = "Data\Lineage_Seqs_EF"
            os.makedirs(output_dir, exist_ok=True)

            for lineage_id, group in df_combined.groupby('lineage_id'):
                root_flag = False
                group_IF_trimmed = group[group['follicle'] != 'EF'].drop_duplicates(subset='v_mutations')
                group_trimmed = group.drop_duplicates(subset='v_mutations')
                if len(group_IF_trimmed) >= min_lineage_size:
                    filename = os.path.join(output_dir, f"lin_{int(float(lineage_id))}_inc_ef.fasta")
                    with open(filename, "w") as fasta_file:
                        seq_len = None
                        for idx, row in group_trimmed.iterrows():
                            if not isinstance(row['v_mutations'], str):
                                root_flag = True
                                header = f">ROOT__{idx}"
                            else:
                                header = f">CLONE_{idx}"
                            sequence = row['v_sequence_no_trunc']
                            if seq_len is None:
                                seq_len = len(sequence)
                            elif seq_len != len(sequence):
                                print(lineage_id)
                            fasta_file.write(f"{header}\n{sequence}\n")

                        if not root_flag:
                            root_seq = get_root_seq(group_trimmed.iloc[0])
                            fasta_file.write(f">ROOT__INFERRED\n{root_seq}\n")

    return lineage_ids, lineage_lens, N_trees, N_lineages_tot

def data_readin(tonsil_vdjs_path, metadata_path, UMI_collapsed_path, clone_list_path):
    # Read in processed sequence data
    with open(tonsil_vdjs_path, 'r') as file:
        df = pd.read_csv(file, sep='\t', dtype=str)
    df['lineage_id'] = df['lineage_id'].astype(float).astype(int)

    # Read in spatial barcode metadata
    with open(metadata_path, 'r') as file:
        metadata_df = pd.read_csv(file, sep=',', dtype=str)
    metadata_df['plasma_frac'] = metadata_df['Plasmablast'].astype(float) / (metadata_df['Plasmablast'].astype(float) + metadata_df['B.Cells'].astype(float))

    # Read in T cell receptor data
    df_collapsed = pd.read_csv(UMI_collapsed_path)
    df_clones = pd.read_csv(clone_list_path)
    cloneid_to_type = df_clones.set_index('cloneId')['type']
    df_collapsed['type'] = df_collapsed['cloneId'].map(cloneid_to_type)
    df_TRB = df_collapsed[df_collapsed['type'] == 'TRB']
    df_TRB = df_TRB.drop_duplicates(subset=['cloneId', 'st_umi'])
    df_TRB = df_TRB.merge(metadata_df[['spatial_bc','x','y','Follicles_seurat']], left_on='st_barcode', right_on='spatial_bc', how='left')
    df_TRB.drop(columns=['spatial_bc'], inplace=True)
    df_TRB['x'] = df_TRB['x'].astype(int)
    df_TRB['y'] = df_TRB['y'].astype(int)

    # Focus on heavy chain reads
    df_IGH = df[df['locus'] == 'IGH'].copy()

    # Set follicular labels
    df_IGH = df_IGH.merge(metadata_df[['spatial_bc', 'section', 'Follicles_seurat','x','y','plasma_frac']], left_on='st_barcode', right_on='spatial_bc', how='left')

    df_IGH['follicle'] = df_IGH['Follicles_seurat']
    df_IGH.drop(columns=['spatial_bc', 'Follicles_seurat'], inplace=True)
    df_IGH['follicle'].replace('nonFoll', 'EF', inplace=True)

    mask = [type(df_IGH['follicle'][i]) != type('a') for i in range(len(df_IGH))].copy()
    df_IGH['follicle'].replace(df_IGH[mask]['follicle'].values[0], 'EF', inplace=True)

    # Prune lineages with multiple v calls
    max_id = 0
    n_new_lineage = 0

    for lineage_id, group in df_IGH.groupby('lineage_id'):
        max_id = max(max_id, lineage_id)
    for lineage_id, group in df_IGH.groupby('lineage_id'):
        v_counts = group['v_db_call'].value_counts()

        if len(v_counts) > 1:
            majority_call = v_counts.idxmax()
            for v_call in v_counts.index:
                if v_call != majority_call:
                    df_IGH.loc[(df_IGH['lineage_id'] == lineage_id) & (df_IGH['v_db_call'] == v_call), 'lineage_id'] = max_id + 1
                    max_id += 1
                    n_new_lineage += 1

    print(f"Added {n_new_lineage} new lineages after pruning inconsistent V gene calls.")

    # Pad V gene sequences to ensure consistent lengths
    mask = (df_IGH['lineage_id'] == 410) & (df_IGH['v_sequence_no_trunc'].apply(len) == 288)
    df_IGH.loc[mask, 'v_sequence_no_trunc'] = 'N' + df_IGH.loc[mask, 'v_sequence_no_trunc']
    mask = (df_IGH['lineage_id'] == 785) & (df_IGH['v_sequence_no_trunc'].apply(len) == 283)
    df_IGH.loc[mask, 'v_sequence_no_trunc'] = 'NNNNNN' + df_IGH.loc[mask, 'v_sequence_no_trunc']
    mask = (df_IGH['lineage_id'] == 2392) & (df_IGH['v_sequence_no_trunc'].apply(len) == 283)
    df_IGH.loc[mask, 'v_sequence_no_trunc'] = 'NNNNNN' + df_IGH.loc[mask, 'v_sequence_no_trunc']
    mask = (df_IGH['lineage_id'] == 4985) & (df_IGH['v_sequence_no_trunc'].apply(len) == 287)
    df_IGH.loc[mask, 'v_sequence_no_trunc'] = 'NN' + df_IGH.loc[mask, 'v_sequence_no_trunc']
    mask = (df_IGH['lineage_id'] == 5836) & (df_IGH['v_sequence_no_trunc'].apply(len) == 283)
    df_IGH.loc[mask, 'v_sequence_no_trunc'] = 'NNNNNN' + df_IGH.loc[mask, 'v_sequence_no_trunc']
    mask = (df_IGH['lineage_id'] == 13105) & (df_IGH['v_sequence_no_trunc'].apply(len) == 280)
    df_IGH.loc[mask, 'v_sequence_no_trunc'] = 'NNNNNN' + df_IGH.loc[mask, 'v_sequence_no_trunc']

    # Split lineages where V genes have different lengths and cannot be padded
    mask = (df_IGH['lineage_id'] == 2430) & (df_IGH['v_sequence_no_trunc'].apply(len) == 289)
    df_IGH.loc[mask, 'lineage_id'] = max_id + 1
    max_id += 1
    mask = (df_IGH['lineage_id'] == 2981) & (df_IGH['v_sequence_no_trunc'].apply(len) == 289)
    df_IGH.loc[mask, 'lineage_id'] = max_id + 1
    max_id += 1

    # Delete rows with insertions or deletions in mutation calls
    print(f"Removing {df_IGH.apply(flag_insertions_deletions, axis=1).sum()} reads with indels in V gene mutation calls.")
    df_IGH = df_IGH[~df_IGH.apply(flag_insertions_deletions, axis=1)].copy()

    # Determine consensus germline for each lineage
    update_ct = 0
    for lineage_id, group in df_IGH.groupby('lineage_id'):
        root_seqs = group.apply(get_root_seq, axis=1)
        unique_root_seqs = root_seqs.value_counts().to_dict()
        if len(unique_root_seqs) > 1:
            max_ct = 0
            consensus_seq = None
            for seq, count in unique_root_seqs.items():
                if count > max_ct:
                    max_ct = count
                    consensus_seq = seq
            for idx, row in group.iterrows():
                if get_root_seq(row) != consensus_seq:
                    df_IGH.at[idx, 'v_mutations'] = get_muts_from_consensus_seq(row['v_sequence_no_trunc'], consensus_seq)
                    update_ct += 1

    print(f"Updated {update_ct} mutation calls to ensure single consensus germline per lineage.")

    df_IGH_IF = df_IGH[df_IGH['follicle'] != 'EF'].copy()
    df_IGH_EF = df_IGH[df_IGH['follicle'] == 'EF'].copy()
    df_IGH_EF_PBdom = df_IGH[np.logical_and(df_IGH['follicle'] == 'EF', df_IGH['plasma_frac'] > 0.9)].copy()
    df_combined = pd.concat([df_IGH_IF, df_IGH_EF_PBdom]).copy()

    return df_IGH_IF, df_IGH_EF_PBdom, df_combined, df_TRB, metadata_df, df_IGH_EF, df_IGH

def safe_fsolve(func, x0, args):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", RuntimeWarning) 
        result = fsolve(func, x0, args=args)

        for warning in w:
            if issubclass(warning.category, RuntimeWarning):
                print("Error", args, result, func)

        return result

def poiss_CI(n_obs, div, alpha = 0.05):

   def eq1(log_tau_min, n_obs, div, alpha):
      return poisson.cdf(n_obs, div / np.exp(log_tau_min)) - alpha
   
   def eq2(log_tau_max, n_obs, div, alpha):
      return poisson.sf(n_obs-1, div / np.exp(log_tau_max)) - alpha

   if n_obs == 0:
       tau_min = safe_fsolve(eq1, x0=np.log(div/2), args=(n_obs, div, alpha))[0]
       tau_max = np.inf
   else:
      tau_min = safe_fsolve(eq1, x0=np.log(div/(2*n_obs)), args=(n_obs, div, alpha))[0]
      tau_max = safe_fsolve(eq2, x0=np.log(2*div/n_obs), args=(n_obs, div, alpha))[0]

   if n_obs == 0: return np.inf, np.exp(tau_min), np.inf
   return div/n_obs, np.exp(tau_min), np.exp(tau_max)

def surv_func(data, x_vals = None):
   if x_vals is not None:
      survival = np.zeros(len(x_vals))
      for i in range(len(x_vals)): survival[i] = np.sum(data >= x_vals[i]) / len(data)
      return x_vals, survival
   else:
      data_mat = np.copy(data)
      data_sorted = data_mat[np.argsort(data_mat)]
      survival = 1 - np.arange(len(data_sorted)) / len(data_sorted)
      return data_sorted, survival

    # Fit survival functions of inferred rates
def fit_survival_func_single_rate(
    rates_data, phylo_divs, rate_init, 
    N_samples=100, 
    test_vals=np.linspace(0.1, 2, 200),
    points_side=20,
    plot_fit=False
):
    err_vals = np.zeros(len(test_vals))

    rates_data_sorted, survival_data = surv_func(rates_data)

    for i in range(len(test_vals)):
        for j in range(N_samples):
            sampled_counts = np.random.poisson(test_vals[i] * phylo_divs * rate_init)
            filt = sampled_counts >= 1

            # If not enough samples, skip this iteration
            if sum(filt) <= 5:
                err_vals[i] += 1000
                continue

            tot_ct = np.copy(sampled_counts[filt])
            tot_denom = np.copy(phylo_divs[filt])

            rates_sampled_sorted, survival_sampled = surv_func(tot_ct / tot_denom)

            x_min = min(rates_data_sorted[0], rates_sampled_sorted[0])
            x_max = max(rates_data_sorted[-1], rates_sampled_sorted[-1])
            x_common = np.logspace(np.log10(x_min), np.log10(x_max), 500)

            interp_survival = interp1d(rates_data_sorted, survival_data, bounds_error=False, fill_value=(1, 0))
            interp_sampled = interp1d(rates_sampled_sorted, survival_sampled, bounds_error=False, fill_value=(1, 0))

            err_vals[i] += np.sum(
                (interp_survival(x_common) - interp_sampled(x_common))**2
            ) / (len(x_common) * N_samples)

    min_idx = np.argmin(err_vals)

    # Check bounds for parabola fit
    if min_idx < points_side or min_idx >= len(err_vals) - points_side:
        print("Check range of tested values")
        plt.plot(test_vals[err_vals < 1000], err_vals[err_vals < 1000])
        plt.xlabel("Test rate × rate_init")
        plt.ylabel("Error")
        plt.show()
        return test_vals[min_idx] * rate_init

    # Select points around the minimum
    idx_range = range(min_idx - points_side, min_idx + points_side + 1)
    x_fit = test_vals[list(idx_range)]
    y_fit = err_vals[list(idx_range)]

    # Fit parabola: y = ax² + bx + c
    coeffs = np.polyfit(x_fit, y_fit, 2)
    a, b, c = coeffs

    if a <= 0: rate_opt = test_vals[min_idx]
    else:
        # Vertex of parabola
        rate_opt = -b / (2 * a)

    if plot_fit:
        xs = np.linspace(min(x_fit), max(x_fit), 200)
        plt.plot(test_vals, err_vals, "o-", label="Error values")
        plt.plot(xs, np.polyval(coeffs, xs), "r--", label="Parabola fit")
        plt.axvline(rate_opt, color="k", linestyle=":", label="Parabola min")
        plt.xlabel("Test rate")
        plt.ylabel("Error")
        plt.legend()
        plt.show()

    return rate_opt * rate_init

# Get survival function at single rate for plotting
def get_survival_func_single_rate(phylo_divs, rate, x_vals, N_samples = 100):
    surv_sampled = np.zeros((N_samples, len(x_vals)))
    count_dist = np.zeros(19)
    frac_monofoll = 0
    for i in range(N_samples):
        sampled_counts = np.random.poisson(rate * phylo_divs)

        filt = sampled_counts >= 1
        _, survival_sampled = surv_func(sampled_counts[filt] / phylo_divs[filt], x_vals)
        surv_sampled[i,:] = survival_sampled

        count_dist += np.histogram(sampled_counts,bins=np.arange(0,20,1),density=True)[0]/N_samples

        frac_monofoll += np.sum(sampled_counts == 0) / (len(sampled_counts)*N_samples)
    return np.mean(surv_sampled, axis=0), np.std(surv_sampled,axis=0), frac_monofoll, count_dist