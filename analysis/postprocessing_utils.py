import sys
import glob
import os
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
import collections
from pylab import *
import pandas
import pickle


### Loading data 



def read_in_solutions(input):
    """Read in pickle solution files generated from the deimmunization solver.
 
    Parameters
    ----------
    input : list of str
        Mixed list of single *.pcl files and directories containing *.pcl files.
    zip_filenames : bool, optional (default: False)
        If true, also filenames are returned (along with their corresponding
        `Solution`).
 
    Returns
    -------
    {list of `Solution`, list of tuple (str, `Solution`)}
        If `zip_filenames` is False, then only the solutions are returned.
        Else, also the filenames are returned, zipped to the corresponding
        `Solution`.
 
    """
    res = []
    for fil_name in input:
        with open(fil_name, "r") as f:
            res.append((fil_name,pickle.load(f)))
    return res



def get_mutant_seq(sol):
    """Extract mutant sequence from a single solution.
 
    Parameters
    ----------
    sol : `Solution`
        Solution from de-immunization solver.
 
    Returns
    -------
    str
        Mutant sequence.
 
    """
    return ''.join(item[1] for item in sorted(
        map(lambda x: eval(x[1:].replace(x[-2], "'%s'" % x[-2])),
            filter(lambda x: x.startswith("x("),sol.vars.keys())),
        key=lambda y: y[0]
    ))



### Plots



def get_epitope_mesh(seq, pwms, alleles, eplen=9, binary=False):
    """
    predicts remaining epitopes in seq and returns for a dict with key=(allele_name,start) value=affinity
    
    :param seq: AA string
    :param pwd: dictionaries of PWDS for alleles in alleles
    :param alleles: dictionary with key=allele_name value=binding threshold
    :param binary: determins whether to return binary matrix or affinity (affinity=max(pred-shift, 0))
    """
    res = {}
    for allele_name, pwm in pwms.items():

        for i in range(len(seq)-eplen+1):
            ep = seq[i:i+eplen]
            score = max(sum( pwm[j,aa] for j,aa in enumerate(ep))- alleles[allele_name],0)
            
            if binary:
                score = 1.0 if score > 0 else 0
            if score > 0:
                res[(allele_name,i)]=score
    return res



def load_pwms(alleles, library_path, epi_len=9):
    """
    loads pwms for specified alleles
    
    :param alleles: pandas dataframe of allele that are supported by EVdeimunization
    :param library_path: the path to where EVdeimmunization is saved
    :return: dictionary of dictionaries {allele_name:{pos:{aa:score}}}
    """
    import sys
    sys.path.append(library_path)
    from collections import OrderedDict
    
    pssms = OrderedDict()
    for name, pssm_thresh in zip(alleles.name, alleles.pssm_thresh):
        raw_pssm = getattr(__import__("preprocessing.mhc_matrices", fromlist=[name]), name)

        # flatten PSSM
        pssm = {}
        for epitope_position, d in raw_pssm.items():
            for aa, value in d.items():
                pssm[(epitope_position, aa)] = value

        for epitope_position in range(epi_len):
            if (epitope_position, "C") not in pssm:
                pssm[(epitope_position, "C")] = 0
        pssms[name] = pssm
    return pssms



def calc_wildtype_point(wt, pwms, alleles, ev_couplings, start_pos=0):
    """Calculate immunogenicity and energy of the wildtype.
 
    Parameters
    ----------
    ev_couplings : `EVcouplings`
        Read in e_ij binaries used for preparing input files.
    config : `ConfigParser`
        Configurations.
 
    Returns
    -------
    immuno : float
        Immunogenicity of the wildtype.
    hamiltonian : float
        Hamiltonian of the wildtype.
 
    """
    # find lower-case positions and delete
    idx_map = {start_pos+i:i for i in range(len(wt))}
    idx = set(ev_couplings.index_list).intersection(set(range(start_pos, start_pos+len(wt))))
    cov_seq = "".join(wt[idx_map[i]] for i in sorted(idx)) 
    hamiltonians = ev_couplings.hamiltonians([cov_seq])[0]
    hamiltonian =  -1.0 * hamiltonians[0] 

    # calculate immunogenicity
    immuno = 0
    for i in range(len(wt)-9+1):
        ep = wt[i:i+9]
        for al,pwm in pwms.items():
            immuno += alleles.loc[al,"p"]*max(sum(pwm[j,a] for j,a in enumerate(ep))-alleles.loc[al,"pssm_thresh"] ,0)
    return immuno, hamiltonian



def plot_pareto_front(sols, wt_im, wt_en):
    import re
    
    fig = figure(figsize=(5, 5), dpi=300)
    labels = []
    wt_en = -1.0*wt_en
    x_min,x_max = float("inf"),-float("inf")
    for fn, sol in sols:
        y = []
        x = []
        m = re.search('_k(.+?).pcl', fn)
        name = m.group(1)
        for s in sorted(set(sol)):
            im = s.objs[0]/10000
            en = -1.0*s.objs[1]
            x.append( ((im - wt_im)/wt_im)*100)
            y.append( ((en - wt_en)/wt_en)*100)
        x_min = min(min(x),x_min)
        x_max = max(max(x),x_max)
        plt.plot(x,y,'*--')
        labels.append("k={}".format(name))
    plt.hlines(0.0,x_min,x_max, linestyles="dashed")
    plt.xlabel("Immunogenicity [percent change to wt]")
    plt.ylabel("Hamiltonian [percent change to wt]")
    plt.legend(labels)
    return fig



def plot_epitope_cluster(epitopes, wt_seq, mutant=None, epi_length=9, alleles=None, enrichment=None, enrichment_title="", plot_title=None, start=0, out=None, knwon_epitopes=None):
    """
        plots the ASCII output as actual plot similer to
        DeGroots plots
        :param epitopes: a dict with key: (allele,start_pos)= affinity_pred
        :param mut_seq: a string with the mutated sequence
        :param wt_seq: The wilde type sequence
        :param enrichment: some kind of enrichment scores for each position of the sequence (i.e EC, conservation) must be a dictionary with key=pos, value=score
        :param enrichment_title: string specifying what kind of enrichment it is (gets ploted)
        :param alleles: a dictionary with key: allele_name value: frequency
                        -> is used to calculated weighted immuno scores
        :param plot_title: optional title for plot
        :param start: start is added to the x-axis to match for example thr Uniprot positions
        param knwon_epitopes: is a dictionary similar to epitopes with already known ones. these will be ploted beneath the predicted ones in a different color
    """
 
 
 
    def _get_epitopes(sol, seq_len,  alleles=None, epi_length=9):
        import collections
        import pandas
 
        def __init_seq_dic(seq_len=seq_len):
            return {i:0.0 for i in xrange(seq_len)}
 
        epitopes = collections.defaultdict(__init_seq_dic)
 
        for key, val in sol.iteritems():
            #hack here for degroot plot -> but needed has to be refactored
            if len(key) <= 2:
                allele, pos =  key
            else:
                allele,pos,epi_length = key
 
            for i in xrange(epi_length):
                epitopes[allele][pos+i] += float(val) if alleles is None else float(alleles.get(allele, 1.0))*float(val)
        return pandas.DataFrame(epitopes)
 
    def _thin_lines(axis, linewidth):
        for line in ['top', 'bottom', 'right', 'left']:
            try:
                axis.spines[line].set_linewidth(linewidth)
            except Exception:
                pass
 
 
    if enrichment is not None:
        ec_enriched = enrichment
    else:
        ec_enriched = {i:0.0 for i in xrange(len(wt_seq))}
    epitopes_df = _get_epitopes(epitopes, len(wt_seq),alleles=alleles, epi_length=epi_length)
    # print epitopes_df.to_string()
 
    if knwon_epitopes is not None:
        knwon_epis_df = _get_epitopes(knwon_epitopes, len(wt_seq),alleles=alleles, epi_length=epi_length)
        knwon_epi_cmp = deepcopy(plt.cm.Greens)
        knwon_epi_cmp.set_bad(color='0.75', alpha=None)
 
    wt_tmp = wt_seq.upper()
    mut_pos = [i for i in xrange(len(wt_seq)) if mutant[i] != wt_tmp[i] ] if mutant is not None else []
 
    ##here starts the plot######
    matplotlib.rc('font', family='monospace')
 
    LABEL_X_OFFSET = 0.5
    LABEL_Y_OFFSET = 0.5
    LABEL_SIZE = 7
    LINEWIDTH = 0.05
    LINE_TRANSPARENCY = 0.7
 
    colormap = deepcopy(plt.cm.Reds)
    colormap.set_bad(color='0.75', alpha=None)
 
    matrix_width = len(wt_seq)
    print("Matrix-width", matrix_width)
    ratio = matrix_width / 20.0
    fig = figure(figsize=(ratio*2*1.2, 3), dpi=300)
    if knwon_epitopes is not None:    #
        gs = matplotlib.gridspec.GridSpec(5, 1, height_ratios=[10,10, 2, 1, 15], 
                                          width_ratios = [matrix_width, 1, 1, 1, 1], wspace=0.01, hspace=0.02)
    else:
        #[matrix_width, 1, 1, 1]
        gs = matplotlib.gridspec.GridSpec(4, 1, height_ratios=[10, 2, 1, 15], 
                                          width_ratios = [1] , wspace=0.01, hspace=0.02)
    if plot_title is not None:
        fig.suptitle(plot_title, fontsize=7, fontweight="bold")
 
    offset_axis = 0 if knwon_epitopes is None else 1
    # --------------------------------
    # main epitope matrix
    # --------------------------------
    epitope_matrix = epitopes_df.as_matrix().T*1.0
    #epitope_matrix = epitopes_df.as_matrix().T
    m_max = epitope_matrix.max()
    #ax1 = subplot(gs[0:8,:])
    ax1 = subplot(gs[0])
    #print epitope_matrix
    #print "MAx ", m_max, " MIN ", epitope_matrix.min()
 
 
    cm = pcolormesh(epitope_matrix, cmap=colormap, vmax=5, vmin=0)
 
    row, col = epitopes_df.shape
    for i in xrange(col):
        for j in xrange(row):
            if j in mut_pos:
                marker = plt.Circle((j + 0.5, i + 0.5), 0.1, fc='k')
                plt.gca().add_patch(marker)
            #else:
            #
            marker = plt.Rectangle((j, i), 1, 1, fc='None', linewidth=LINEWIDTH, alpha=LINE_TRANSPARENCY)
            plt.gca().add_patch(marker)
 
 
 
    tick_params(axis='x', labelsize=LABEL_SIZE)
    tick_params(axis='y', labelsize=LABEL_SIZE)
 
    xlim([0, matrix_width])
    xticks(np.arange(LABEL_X_OFFSET, matrix_width+LABEL_X_OFFSET), mutant if mutant is not None else wt_seq)
    ax1.xaxis.tick_top()
 
    if mutant is not None:
        ax_mut = ax1.twiny()
        tick_params(axis='x', labelsize=LABEL_SIZE)
        tick_params(axis='y', labelsize=LABEL_SIZE)
        xlim([0, matrix_width+0.1])
        xticks(np.arange(LABEL_X_OFFSET, matrix_width+LABEL_X_OFFSET), wt_seq)
        ax_mut.xaxis.tick_bottom()
        ax_mut.xaxis.tick_top()
        ax_mut.xaxis.set_ticks_position("none")
        setp(ax_mut.get_yticklabels(), visible=False)
 
    ylim([0, col])
    yticks(np.arange(LABEL_Y_OFFSET, col + LABEL_Y_OFFSET), epitopes_df.columns)
    ax1.invert_yaxis()
 
 
    ax1.xaxis.set_ticks_position("none")
    ax1.yaxis.set_ticks_position("none")
    _thin_lines(ax1, 0.5)
 
 
    # --------------------------------
    # Known Epitopes grid
    # --------------------------------
    if knwon_epitopes is not None:
        ax_known_epi = subplot(gs[1])
        knwon_epis_matrix = knwon_epis_df.as_matrix().T*1.0
        #epitope_matrix = epitopes_df.as_matrix().T
        kn_m_max = knwon_epis_matrix.max()
        knwon_epi_cm = pcolormesh(knwon_epis_matrix, cmap=knwon_epi_cmp, vmax=kn_m_max, vmin=0)
 
        row, col = epitopes_df.shape
        for i in xrange(col):
            for j in xrange(row):
                marker = plt.Rectangle((j, i), 1, 1, fc='None', linewidth=LINEWIDTH, alpha=LINE_TRANSPARENCY)
                plt.gca().add_patch(marker)
 
 
        tick_params(axis='x', labelsize=LABEL_SIZE)
        tick_params(axis='y', labelsize=LABEL_SIZE)
 
        xlim([0, matrix_width])
        xticks(np.arange(LABEL_X_OFFSET, matrix_width+LABEL_X_OFFSET), mutant)
 
 
        ylim([0, col])
        yticks(np.arange(LABEL_Y_OFFSET, col + LABEL_Y_OFFSET), knwon_epis_df.columns)
        ax_known_epi.invert_yaxis()
 
 
        ax_known_epi.xaxis.set_ticks_position("none")
        ax_known_epi.yaxis.set_ticks_position("none")
        _thin_lines(ax1, 0.5)
 
 
    # --------------------------------
    # Middel line EC distribution
    # --------------------------------
    if enrichment is not None:
        #ax2 = subplot(gs[9,:])
        ax2 = subplot(gs[2+offset_axis])
        ec=[0]*len(wt_seq)
        for pos, score in ec_enriched.iteritems():
            ec[pos] = score
        ec_max = np.max(ec)
        ec_min = np.min(ec)-ec_max
        ec_colormap = deepcopy(plt.cm.Blues)
        ec_colormap.set_bad(color='0.75', alpha=None)
        #print "EC max", np.max(ec), "EC min ", np.min(ec)
        cm_EC = pcolormesh(np.array([ec]), cmap=ec_colormap, vmax=ec_max, vmin=0)
 
        row, col = epitopes_df.shape
        for i in xrange(col):
            for j in xrange(row):
                # skip unspecified entries
                marker = plt.Rectangle((j, i), 1, 1, fc='None', linewidth=LINEWIDTH, alpha=LINE_TRANSPARENCY)
                plt.gca().add_patch(marker)
 
        tick_params(axis='x', labelsize=LABEL_SIZE)
        tick_params(axis='y', labelsize=LABEL_SIZE)
 
        ylim([0, 1])
        #yticks(range(0,1), [""])
 
        xlim([0, matrix_width])
        yticks(np.arange(LABEL_Y_OFFSET, 1 + LABEL_Y_OFFSET), [enrichment_title])
        xticks(np.arange(LABEL_X_OFFSET, matrix_width+LABEL_X_OFFSET), wt_seq)
        ax2.xaxis.tick_top()
        ax2.xaxis.set_ticks_position("none")
        ax2.yaxis.set_ticks_position("none")
        setp(ax2.get_xticklabels(), visible=False)
        _thin_lines(ax1, 0.5)
 
    # --------------------------------
    # bottom Immunogenicity distribution
    # --------------------------------
    #ax3 = subplot(gs[10:18,:])
    ax3 = subplot(gs[3+offset_axis])
 
    immuno_scores = list(epitopes_df.sum(axis=1))
    immuno_scores.append(0)
    tmp = immuno_scores[:]
    for i in xrange(len(immuno_scores)-1):
        if immuno_scores[i+1] > 0:
            tmp[i+1] = immuno_scores[i]*1.0
            #tmp[i+1] = immuno_scores[i]
        else:
            tmp[i+1] = 0
 
    tmp.append(0)
    #step(range(matrix_width+1), immuno_scores, where="post")
    fill_between(range(matrix_width+2), tmp, 0, facecolor='#F97E60')
    for i in xrange(6+1):
        for j in xrange(row):
            # skip unspecified entries
            marker = plt.Rectangle((j, i), 1, 1, fc='None', linewidth=LINEWIDTH, alpha=LINE_TRANSPARENCY)
            plt.gca().add_patch(marker)
 
    tick_params(axis='x', labelsize=LABEL_SIZE)
    tick_params(axis='y', labelsize=LABEL_SIZE)
 
    xlim([start, start+len(wt_seq)])
    #xticks(np.arange(LABEL_X_OFFSET, matrix_width+LABEL_X_OFFSET), [start+i for i in xrange(len(wt_seq))])
    ylim([0,6+0.1])
    ax3.xaxis.set_ticks_position("none")
    #ax3.yaxis.set_ticks_position("none")
    setp(ax3.get_xticklabels(), visible=False)
    _thin_lines(ax3, 0.5)
 
 
    # --------------------------------------------
    # Scales
    # --------------------------------------------
    #0.90, 0.4, 0.008, 0.25
    offset = 0.0 if knwon_epitopes is None else 0.05
    cbaxes = fig.add_axes([0.89, 0.63+offset, 0.006, 0.25])
    cb = colorbar(cm, cax = cbaxes, ticks=[0,  3])
    cb.ax.set_yticklabels(["{:+>4.3f}".format(v) for v in [0, 3]])
    cb.ax.xaxis.set_ticks_position("none")
    cb.ax.yaxis.set_ticks_position("none")
    cb.outline.set_linewidth(LINEWIDTH)
    cb.ax.tick_params(labelsize=LABEL_SIZE)
 
    if knwon_epitopes is not None:
        cbaxes = fig.add_axes([0.89, 0.395, 0.006, 0.25])
        cb = colorbar(knwon_epi_cm, cax = cbaxes, ticks=[0, kn_m_max])
        cb.ax.set_yticklabels(["{:+>4.3f}".format(v) for v in [0,kn_m_max]])
        cb.ax.xaxis.set_ticks_position("none")
        cb.ax.yaxis.set_ticks_position("none")
        cb.outline.set_linewidth(LINEWIDTH)
        cb.ax.tick_params(labelsize=LABEL_SIZE)
 
    if enrichment is not None:
        offset = 0.0 if knwon_epitopes is None else -0.23
        cbaxes = fig.add_axes([0.89, 0.33+offset, 0.006, 0.25])
        cb = colorbar(cm_EC, cax = cbaxes, ticks=[0, ec_max])
        cb.ax.set_yticklabels(["{:+>4.3f}".format(v) for v in [0,ec_max]])
        cb.ax.xaxis.set_ticks_position("none")
        cb.ax.yaxis.set_ticks_position("none")
        cb.outline.set_linewidth(LINEWIDTH)
        cb.ax.tick_params(labelsize=LABEL_SIZE)
 
    if out is None:
        plt.draw()
    else:
        savefig(out, bbox_inches='tight')