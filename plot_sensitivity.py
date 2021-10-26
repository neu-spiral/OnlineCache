import matplotlib.pyplot as plt
import argparse
import numpy as np
from matplotlib.dates import date2num
from matplotlib.transforms import Bbox
import datetime
import os
import pickle


graph = {'erdos_renyi': 100, 'balanced_tree': 341, 'hypercube': 128, 'geant': 22, 'dtelekom': 68}
Abbr = {'erdos_renyi': 'ER', 'balanced_tree': 'BT', 'hypercube': 'HC', 'geant': 'geant', 'dtelekom': 'dtelekom'}
Algorithms = ['BT','dtelekom','ER','geant','HC']
colors =['r', 'sandybrown', 'gold', 'darkseagreen', 'c', 'dodgerblue', 'm']
lin_styles = ['-.', '*:', 's-', '^-', 'v-', '--']
dir = "sen_result/"


def lin_ex2(DICS, outfile, xaxis, xaxis_label, yaxis_label, LEGEND):
    fig, ax = plt.subplots()
    #fig.set_size_inches(7,7)
    x_ax = xaxis
    ax.get_xaxis().get_major_formatter().labelOnlyBase = False
    i = 0
    for alg in DICS:
        plt.xscale('log')
        plt.plot(x_ax, np.array(DICS[alg]), lin_styles[i], markersize=15, color=colors[i], label=alg, linewidth=3)
        i += 1
    plt.xlabel(xaxis_label, fontsize=24)
    plt.ylabel(yaxis_label, fontsize=24)
    plt.xticks(fontsize=21)
    plt.yticks(fontsize=21)
    if LEGEND:
        lgd = plt.legend(loc='best', ncol=1, borderaxespad=0., fontsize=21, handlelength=1.5, handletextpad=0.5)

    plt.ylim(0.1, 1)
#    plt.ylim(50, 2000)
    plt.show()

    if LEGEND:
        fig.savefig(outfile + '.pdf', bbox_extra_artists=(lgd,), bbox_inches = 'tight')
    else:
        fig.savefig(outfile + '.pdf')


def read_file(fname):
    with open(fname, 'rb') as f:
        args, construct_stats, optimal_stats, demand_stats, node_stats, network_stats = pickle.load(f)
    time = sorted(network_stats['demand'].keys())[-1]
    tag = network_stats['fun'][time][4] - network_stats['demand'][time]['weight']
    print(network_stats['fun'][time][3]-network_stats['fun'][time][4])
    F = optimal_stats['F']
    return tag/F

def read_num_demand(fname):
    with open(fname, 'rb') as f:
        args, construct_stats, optimal_stats, demand_stats, node_stats, network_stats = pickle.load(f)

    num_demand = 0
    for d in demand_stats:
        num_demand += demand_stats[d]['queries_satisfied']
    return num_demand

def lin_ex(DICS, outfile, xaxis, xaxis_label, yaxis_label, LEGEND):
    fig, ax = plt.subplots()
    #fig.set_size_inches(7,7)
    x_ax = xaxis
    ax.get_xaxis().get_major_formatter().labelOnlyBase = False
    i = 0
    for alg in DICS:
        plt.xscale('log')
        plt.plot(x_ax[alg], np.array(DICS[alg]), lin_styles[i], markersize=15, color=colors[i], label=alg, linewidth=3)
        i += 1
    plt.xlabel(xaxis_label, fontsize=24)
    plt.ylabel(yaxis_label, fontsize=24)
    plt.xticks(fontsize=21)
    plt.yticks(fontsize=21)
    if LEGEND:
        lgd = plt.legend(loc='best', ncol=1, borderaxespad=0., fontsize=21, handlelength=1.5, handletextpad=0.5)

    plt.ylim(0.1, 1)
#    plt.ylim(50, 2000)
    plt.show()

    if LEGEND:
        fig.savefig(outfile + '.pdf', bbox_extra_artists=(lgd,), bbox_inches = 'tight')
    else:
        fig.savefig(outfile + '.pdf')




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate plots comparing different algorithms.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--lgd', action='store_true', help='Pass to make a legened.')
    args = parser.parse_args()
    #trace_location = 'None'
    #trace_location = 'traces/fixed_popularity_catalog_50.pkl'
    trace_location = 'traces/changing_popularity_catalog_50.pkl'

    # DICS = {'BT': [],'dtelekom': [],'ER': [],'geant': [],'HC': []}
    #
    # for action_selector_eta in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]:
    #     for graph_type in graph:
    #         problem_instance = "%s_TBGRD_%s_%dnodes_%fchange_%feta_%dcolors_%dfreq_%fT" % (
    #             graph_type, trace_location[7:11], graph[graph_type], 0., action_selector_eta, 100, -1, 1)
    #         DICS[Abbr[graph_type]].append(read_file(dir+problem_instance))
    #
    # outfile = dir + 'eta_' + trace_location[7:11]
    # lin_ex2(DICS, outfile, [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2],r'$\epsilon$','TACG', args.lgd)
    #
    # DICS = {'BT': [],'dtelekom': [],'ER': [],'geant': [],'HC': []}
    #
    # for color in [100, 200, 500, 1000, 2000, 5000]:
    #     for graph_type in graph:
    #         problem_instance = "%s_TBGRD_%s_%dnodes_%fchange_%feta_%dcolors_%dfreq_%fT" % (
    #             graph_type, trace_location[7:11], graph[graph_type], 0., 0.005, color, -1, 1)
    #         DICS[Abbr[graph_type]].append(read_file(dir+problem_instance))
    # outfile = dir + 'color_' + trace_location[7:11]
    # lin_ex2(DICS, outfile, [100, 200, 500, 1000, 2000, 5000],'$M$','TACG', args.lgd)
    #
    # DICS = {'BT': [],'dtelekom': [],'ER': [],'geant': [],'HC': []}
    #
    # for freq in [1, 5, 10, 50, 100, 500, 1000, -1]:
    #     for graph_type in graph:
    #         problem_instance = "%s_TBGRD_%s_%dnodes_%fchange_%feta_%dcolors_%dfreq_%fT" % (
    #             graph_type, trace_location[7:11], graph[graph_type], 0., 0.005, 100, freq, 1)
    #         DICS[Abbr[graph_type]].append(read_file(dir + problem_instance))
    # outfile = dir + 'frequency_' + trace_location[7:11]
    # lin_ex2(DICS, outfile, [1, 5, 10, 50, 100, 500, 1000, 5000], '$K$', 'TACG', args.lgd)


    DICS = {'BT': [],'dtelekom': [],'ER': [],'geant': [],'HC': []}
    xaxis = {'BT': [],'dtelekom': [],'ER': [],'geant': [],'HC': []}


    for T in [0.5, 1, 2, 5, 10, 20]:
        for graph_type in graph:
            problem_instance = "%s_TBGRD_%s_%dnodes_%fchange_%feta_%dcolors_%dfreq_%fT_4156910908seed" % (
                graph_type, trace_location[7:11], graph[graph_type], 0., 0.005, 100, -1, T)
            DICS[Abbr[graph_type]].append(read_file(dir + problem_instance))
            xaxis[Abbr[graph_type]].append(read_num_demand(dir + problem_instance)*T/5000)
    outfile = dir + 'averageR_' + trace_location[7:11]
    lin_ex(DICS, outfile, xaxis, '$\overline{|\mathcal{R}_t}|$', 'TACG', args.lgd)
