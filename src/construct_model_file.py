import sys
import argparse
import os
import yaml
from collections import OrderedDict
import pandas as pd
import numpy as np
import pickle
from Bio import AlignIO

from preprocessing.utils import EVcouplings, frequencies

ALPHABET_PROTEIN_NOGAP = "ACDEFGHIKLMNPQRSTVWY"
MODEL_FILE_IMM = "AMPL/continuous_immuno_hemiltonian_imm.mod"
MODEL_FILE_EN = "AMPL/continuous_immuno_hemiltonian_en.mod"

# Table 2 Supertype clusters as defined by Lund et al. (25) and their HLA
# allele representatives used for binding affinity prediction.
# Supertype Representative Supertype Representative
# A01 A*01:01 B08 B*08:01
# A02 A*02:01 B27 B*27:05
# A03 A*03:01 B39 B*39:01
# A24 A*24:02 B44 B*40:01
# A26 A*26:01 B58 B*58:01
# B07 B*07:02 B62 B*15:01


def main():
    args = parse_args()

    # make relative paths absolute
    pwd = os.path.abspath(os.path.dirname(__file__))
    model_file_imm = os.path.join(pwd, MODEL_FILE_IMM)
    model_file_en = os.path.join(pwd, MODEL_FILE_EN)

    # read in config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    if config["sets"]["ignore_pos"] is None:
        config["sets"]["ignore_pos"] = []

    if config["sets"]["exclude_pos"] is None:
        config["sets"]["exclude_pos"] = []

    # read in alignment
    ali = AlignIO.read(args.alignment, "fasta")
    wild_type = str(ali[0].seq)
    offset = 0
    for a in wild_type:
        if a.islower():
            offset += 1
        else:
            break

    # read in EV couplings binary file
    ev_couplings = EVcouplings(
        args.model,
        file_format=args.ev_file_format
    )

    # read in alleles
    alleles = pd.read_csv(
        config["sets"]["allele_file"],
        names=["name", "pssm_thresh", "p"]
    )

    # read in PSSMs
    pssms = OrderedDict()
    pssm_thresh_norm = []
    for name, pssm_thresh in zip(alleles.name, alleles.pssm_thresh):
        raw_pssm = getattr(__import__("preprocessing.mhc_matrices", fromlist=[name]), name)

        # flatten PSSM
        pssm = {}
        for epitope_position, d in raw_pssm.items():
            for aa, value in d.items():
                pssm[(epitope_position, aa)] = value

        for epitope_position in range(config["parameters"]["epi_len"]):
            if (epitope_position, "C") not in pssm:
                pssm[(epitope_position, "C")] = 0

        # normalize PSSM
        pssm_mean = np.mean(list(pssm.values()))
        pssm_std = np.std(list(pssm.values()))
        pssm = {
            k: (v - pssm_mean) / pssm_std
            for k, v in pssm.items()
        }

        pssm_thresh_norm.append(
            (pssm_thresh - config["parameters"]["epi_len"] * pssm_mean) / pssm_std
        )

        pssms[name] = pssm

    alleles["pssm_thresh_norm"] = pssm_thresh_norm

    data_file = args.out_basename + ".data"
    with open(data_file, "w") as f:
        # amino acid alphabet, alleles
        f.write("set SIGMA := {};\n".format(' '.join(ALPHABET_PROTEIN_NOGAP)))
        f.write("set A := {};\n\n".format(' '.join(alleles.name)))

        # wild type length, epitope length, number of mutations
        f.write("param N := {};\n".format(len(wild_type)))
        f.write("param eN := {};\n".format(config['parameters']['epi_len']))
        f.write("param k := {};\n\n".format(config['parameters']['k']))

        f.write("param pssm_thresh :=\n")
        for name, _, _, pssm_thresh_norm in alleles.itertuples(index=False):
            f.write("\t{}\t{}\n".format(name, pssm_thresh_norm))
        f.write(";\n\n")

        f.write("param p :=\n")
        for name, _, p, _ in alleles.itertuples(index=False):
            f.write("\t{}\t{}\n".format(name, p))
        f.write(";\n\n")

        # PSSMs
        f.write("param pssm :=\n")
        for name, pssm in pssms.items():
            f.write("[{},*,*]: {} :=\n".format(name, ' '.join(map(str, range(1, config['parameters']['epi_len'] + 1)))))
            for aa in ALPHABET_PROTEIN_NOGAP:
                f.write(aa + "\t")
                for i in range(config['parameters']['epi_len']):
                    f.write(str(pssm[(i, aa)]) + " ")
                f.write("\n")
        f.write(";\n\n")

        # valid indices of the model

        model_inds = [
            i + 1 + offset for i in ev_couplings.index_map.values()
            if i + offset not in config["sets"]["ignore_pos"]
        ]

        # valid pair indices of the model
        eij_inds = sorted(set(
            [(i, j) if i < j else (j, i)
             for i in model_inds
             for j in model_inds
             if i != j and (i < j or j < i)]
        ))

        # pair indices
        eij_inds_str = "\t".join(map(lambda t: "{} {}".format(t[0], t[1]), eij_inds))
        f.write("set Eij := {};\n".format(eij_inds_str))

        # single indices
        hi_inds_str = "\t".join(map(str, model_inds))
        f.write("set E := {};\n\n".format(hi_inds_str))

        # wild type sequence
        for i, aa in enumerate(wild_type):
            f.write("set WT[{}] := {};\n".format(i+1, aa.upper()))
        f.write("\n")

        # compute mutations based on aa frequencies
        mutations = [set(res.upper()) for res in wild_type]
        matrix = np.array([list(rec) for rec in ali])

        f_i = frequencies(matrix, ev_couplings)
        print("fi", args.freq_thresh)
        excluded_pos = set(config["sets"]["ignore_pos"] + config["sets"]["exclude_pos"])
        for i, muts in enumerate(mutations):
             if i not in excluded_pos:
                above_thresh = [aa for aa in list(ALPHABET_PROTEIN_NOGAP)
                                    if f_i[i, ev_couplings.alphabet_map[aa]] >= args.freq_thresh]
                muts.update(above_thresh)

        # allowed mutations per position
        for i, muts in enumerate(mutations):
            f.write("set M[{}] := {};\n".format(i+1, ' '.join(muts)))
        f.write("\n")

        # single site EV parameters
        hi = {}
        for i in model_inds:
            for aa in list(ALPHABET_PROTEIN_NOGAP):
                hi[(i, aa)] = -1.0 * ev_couplings.h_i[
                    i-1-offset,
                    ev_couplings.alphabet_map[aa]
                ]
        f.write("param h: {} :=\n".format(' '.join(list(ALPHABET_PROTEIN_NOGAP))))
        for i in model_inds:
            f.write(str(i) + "\t" + " ".join(str(hi[i, aa]) for aa in ALPHABET_PROTEIN_NOGAP) + "\n")
        f.write(";\n\n")

        # evolutionary couplings
        eij = {}
        for i, j in eij_inds:
            for ai in mutations[i - 1]:
                for aj in mutations[j - 1]:
                    eij[(i, j, ai, aj)] = -1.0 * ev_couplings.J_ij[
                        i-1-offset,
                        j-1-offset,
                        ev_couplings.alphabet_map[ai],
                        ev_couplings.alphabet_map[aj]
                    ]
        f.write("param eij :=\n")
        for i, j in eij_inds:
            mut_list = list(mutations[j - 1])
            f.write("[{},{},*,*]: {} :=\n".format(i,j,' '.join(mut_list)))
            for ai in mutations[i - 1]:
                f.write(ai + "\t" + " ".join(str(eij[i, j, ai, aj]) for aj in mut_list) + "\n")
        f.write(";\n\n")
        f.write("end;")

    base_path, _ = os.path.splitext(args.out_basename)

    # generate immunogenicity and energy lp files
    for model_file, lp_file in [
        (model_file_imm, base_path + "_imm.lp"),
        (model_file_en, base_path + "_en.lp")
    ]:
        command = "glpsol -m {} -d {} --check --wcpxlp {}".format(model_file, data_file, lp_file)
        print("RUN COMMAND: {}".format(command))
        os.system(command)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate specification files necessary to run EVdeimmunization"
    )

    parser.add_argument("alignment", help="alignment in fasta format that served as input to EVcouplings")
    parser.add_argument("model", help="binary eij couplings file with the biotherapeutic sequence as target")
    parser.add_argument("config", help="config file in YAML format")
    parser.add_argument("out_basename", help="basename of out file (multiple files will be created "
                                             "with appropriate extensions added)")

    parser.add_argument("--freq_thresh", "-t", type=float, default=0.0,
                        help="amino acid frequency threshold used to determine allowed mutations")
    parser.add_argument("--ev_file_format", "-f", choices=["plmc_v1", "plmc_v2"], default="plmc_v2",
                        help="file format of EVcouplings model file (default: 'plmc_v2')")

    return parser.parse_args()


if __name__ == '__main__':
    main()