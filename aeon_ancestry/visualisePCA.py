import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mp_lines

from aeon_ancestry.genotype_from_vcf import Genotypes

def plotReferencePCA (num_loci=128097):
    # DATA
    ref_pcas = pd.read_table("refs/reference_PC3s.txt", sep=" ")

    colour_map = {
            "ACB": "#2BCE48",   "ASW": "#94FFB5",   "BEB": "#990000",    "CDX": "#C20088",
            "CEU": "#FFA405",   "CHB": "#740AFF",   "CHS": "#4C005C",    "CLM": "#5EF1F2",
            "ESN": "#005C31",   "FIN": "#FFFF80",   "GBR": "#FF0010",    "GIH": "#FFCC99",
            "GWD": "#8F7C00",   "IBS": "#FF5005",   "ITU": "#993F00",    "JPT": "#F0A3FF",
            "KHV": "#FFA8BB",   "LWK": "#E0FF66",   "MSL": "#426600",    "MXL": "#00998F",
            "PEL": "#003380",   "PJL": "#808080",   "PUR": "#0075DC",    "STU": "#191919",
            "TSI": "#FFE100",   "YRI": "#9DCC00"
            }
    pop_order = ['ESN','YRI','ACB','ASW','GWD','MSL','LWK','GIH','PJL',
                'BEB','STU','ITU','CHB','CHS','CDX','KHV','JPT','CEU',
                'GBR','IBS','TSI','FIN','PEL','MXL','CLM','PUR']

    # PLOTTING
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    fig.set_figheight(16) #16
    fig.set_figwidth(16)  #16

    #Fig1: PC1 vs PC2
    ax1.grid(color='0.9')
    ax1.set_axisbelow(True)
    ax1.scatter(ref_pcas["PC1"]/num_loci, ref_pcas["PC2"]/num_loci, c=ref_pcas["Population"].map(colour_map), s=6)
    ax1.tick_params(top=True, left=True, labelleft=False)
    ax1.xaxis.set_label_position("top")
    ax1.set_xlabel("PC1", fontsize=20)
    ax1.set_ylabel("PC2", fontsize=20)

    #Fig2: PC3 vs PC2
    ax2.grid(color='0.9')
    ax2.set_axisbelow(True)
    ax2.scatter(ref_pcas["PC3"]/num_loci, ref_pcas["PC2"]/num_loci, c=ref_pcas["Population"].map(colour_map), s=6)
    ax2.tick_params(top=True, bottom=False, labelbottom=False, left=False, labelleft=False, right=True, labelright=False)
    ax2.xaxis.set_label_position("top")
    ax2.yaxis.set_label_position("right")
    ax2.set_xlabel("PC3", fontsize=20)
    ax2.set_ylabel("PC2", fontsize=20)


    # Fig3 done slightly differently to aid legend creation
    # Adapted from https://stackoverflow.com/questions/47006268/scatter-plot-with-color-label-and-legend-specified-by-c-option
    ax3.grid(color='0.9')
    ax3.set_axisbelow(True)
    for g in pop_order:
        gpoints = ref_pcas[ref_pcas["Population"] == g]
        ax3.scatter(gpoints["PC1"]/num_loci, gpoints["PC3"]/num_loci, c=colour_map[g], s=6, label=g)

    ax3.tick_params(bottom=True, labelbottom=False, left=True, labelleft=False)
    ax3.set_xlabel("PC1", fontsize=20)
    ax3.set_ylabel("PC3", fontsize=20)
    pop_legend = ax3.legend(loc='center left', bbox_to_anchor=(1.3,0.5), title="Population")
    ax3.add_artist(pop_legend)

    fig.delaxes(ax4)
    plt.subplots_adjust(wspace=0, hspace=0)
    return fig, [ax1, ax2, ax3]

def transformToPCA (genotypes: Genotypes, subset=False):
    
    scale_f = np.loadtxt("refs/g1k_PCA_scale_factors.csv")
    centres = np.loadtxt("refs/g1k_PCA_centres.csv")
    rot = np.genfromtxt("refs/g1k_PCA_rotations.csv", delimiter=",")
    n = 128097

    if subset:
        full_set_df = pd.read_table("refs/g1k_allele_freqs.txt").filter(["VAR_ID"])
        subids = genotypes.variantList()
        subset = full_set_df[full_set_df["VAR_ID"].isin(subids)]
        n = len(subids)

        scale_f = scale_f[subset.index]
        centres = centres[subset.index]
        rot = rot[subset.index,:]

    coords = dict()
    for sample in genotypes.getSamples():
        d = genotypes.dosageForSample(sample)
        scaled_gt = ((d - centres) / scale_f)
        rotated_gt = np.matmul(scaled_gt, rot)

        coords[sample] = rotated_gt[0:3]

    return (pd.DataFrame(coords, index=["PC1", "PC2", "PC3"])/n)

def addTrioToPCAplot(coord_df, base_axs):
    # Requires an array of length 3 for base_axs, with axes:
    #       axs[0]: x -> PC1, y -> PC2
    #       axs[1]: x -> PC3, y -> PC2
    #       axs[3]: x -> PC1, y -> PC3
    # (i.e. 2nd return value of plotReferencePCA)

    samples = coord_df.columns
    colours = ['crimson','mediumblue','lime']
    markers = []

    i = 0
    for s in samples:
        base_axs[0].plot(coord_df[s]["PC1"], coord_df[s]["PC2"], marker="D", c=colours[i], mec='black', label=s)
        base_axs[1].plot(coord_df[s]["PC3"], coord_df[s]["PC2"], marker="D", c=colours[i], mec='black', label=s)
        base_axs[2].plot(coord_df[s]["PC1"], coord_df[s]["PC3"], marker="D", c=colours[i], mec='black', label=s)
        
        m = mp_lines.Line2D([],[], color=colours[i], marker="D", mec='black', linestyle='None', label=s) 
        markers.append(m)

        i += 1

    base_axs[2].legend(handles=markers, loc='center left', bbox_to_anchor=(1.5, 0.5), title="Samples")
    return base_axs

def saveTrioPCAplot(title, case_coords, verbose=True, num_loci=128097):
    fig, axs = plotReferencePCA(num_loci)
    fig.suptitle(title, size=25, y=0.95)
    updated = addTrioToPCAplot(case_coords, axs)

    plt.savefig(f"{title}_PCA_plot.png")
    if verbose: print (f"saved figure to {title}_PCA_plot.png")

def saveIndividualPCAplot (sample, all_coords, num_loci=128097):
    fig, axs = plotReferencePCA(num_loci)
    fig.suptitle(sample, size=25, y=0.95)

    s_coords = all_coords.filter(items=[sample])
    updated = addTrioToPCAplot(s_coords, axs)

    plt.savefig(f"{sample}_PCA_plot.png")
