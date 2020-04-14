# EVdeimmunization


## Introduction:
This software de-immunizes therapeutic sequences by solving the dual optimization problem of reducing immunogenicity and preserving protein function using model output from Evolutionary Couplings.

## Requirement:
1. Python 2.7
2. Numpy 1.9.1
3. Cplex (+Python API) 12.5
4. Polygon 2.0.7
5. Psutil

## Installation:
After installing the above requirements, simply clone this repository onto your machine. Add the folder to your path in `.bashrc` to run from any location.

## Basic Usage:
### Preprocessing:
This module preprocesses the model file from EV couplings for use with EVdeimmunization.
```
usage: deimm preprocess [options] alignment model config out_basename

positional arguments:
  alignment             alignment in fasta format that served as input to
                        EVcouplings
  model                 binary eij couplings file with the biotherapeutic
                        sequence as target
  config                config file in YAML format
  out_basename          basename of out file (multiple files will be created
                        with appropriate extensions added)

optional arguments:
  -h, --help            show this help message and exit
  --freq_thresh FREQ_THRESH, -t FREQ_THRESH
                        amino acid frequency threshold used to determine
                        allowed mutations
  --ev_file_format {plmc_v1,plmc_v2}, -f {plmc_v1,plmc_v2}
                        file format of EVcouplings model file (default:
                        'plmc_v2')
```

The config file contains parameters for building the model and pointers to an allele file. See an example below:
```yaml
general:
    # amino acid frequency threshold used to determine allowed mutations
    frequency_thresh: 0.01

sets:
    # path to PSSM file
    allele_file: peptide_design_allele_file.txt
    # positions excluded from mutations
    exclude_pos: [1, 2, 3, 50, 51]
    # positions completely ignored by model, i.e. mutations excluded AND eijs do not influence decision process
    ignore_pos:

parameters:
    # epitope length
    epi_len: 9
    # number of mutations to introduce
    k: 2
```
The example allele file includes representative alleles of supertypes [1]. Along with each allele is a PSSM threshold value to select peptides in the top 1% of binders and a frequency parameter. 

TODO: Describe how PSSM threshold value is calculated
```
A0101,0.086238,1.0
A0201,0.059324,1.0
B2705,-0.006773,1.0
A1101,0.123383,1.0
B3901,-0.033835,1.0
A2402,-0.004717,1.0
B4001,0.093183,1.0
A2601,0.049978,1.0
B5801,0.087361,1.0
B3501,0.081655,1.0
B1501,0.102201,1.0
```

### Solve:
This module takes the model files from the preprocessing step and solves the bi-objective mixed integer problem using a rectangle splitting approach [2].
```
usage: deimm solve  [options] --port PORT  --output OUTPUT model_imm.pl model_en.pl

positional arguments:
  model                 model files (generated from preprocessing)

optional arguments:
  -h, --help            show this help message and exit
  --port PORT, -p PORT  Port
  --output OUTPUT, -o OUTPUT
                        Solution output file (ex: output.pcl)
  --approximate APPROXIMATE, -a APPROXIMATE
                        Bound on approximation (default:0.009)
  --relTol RELTOL, -rel RELTOL
                        The relative tolerance for floating point comparison
  --absTol ABSTOL, -abs ABSTOL
                        The absolut tolerance for floating point comparison, also used as epsilon in the model
  -t THREADS, --threads THREADS
                        Number of threads (default: based on cpu count)
  --verbose VERBOSE, -v VERBOSE
                        Verbosity (default 0)
```

## Examples:
TODO: Find fast example
### Test solver installation:

```
deimm solve ./test/test_imm.lp ./test/test_en.lp -p 6882 -o ./test/out.pcl -v 1
```

### Extended example:
First process the EV couplings output:
```
deimm preprocess ./example/PvLEA4_repeats_b0.75.fasta \
                 ./example/PvLEA4_repeats_b0.75.model \
                 ./example peptide_design_config.cfg \ 
                 ./example/PvLEA4_repeats_b0.75.k2
```
Next, solve:
```
deimm solve ./example/PvLEA4_repeats_b0.75.k2_imm.pl ./example/PvLEA4_repeats_b0.75.k2_en.pl \
            -o ./example/PvLEA4_repeats_b0.75.k2.pcl \ 
            -p 6882 \
            -v 1
```

## Advanced Usage:
TODO: Add link to FRED2 EVDeimmunization

The solver can be used in distributed systems. For that, first start the manager script, which implements the algorithm and handles the work distribution.
Make sure that the chosen IP and Port combination is reachable from all other used systems.
```
usage: deimmunization_manager.py [-h] [--grid GRID] --port PORT
                                 [--approximate APPROXIMATE] [--key KEY]
                                 --output OUTPUT [--resolve RESOLVE]

Rectangle Manager implementation

Arguments:
  -h, --help            show this help message and exit
  --grid GRID, -g GRID  Number of Epsilon grid points (default:3)
  --port PORT, -p PORT  Port
  --approximate APPROXIMATE, -a APPROXIMATE
                        Bound on approximation (default:0.009)
  --key KEY, -k KEY     Authentication key (default: rectangle)
  --output OUTPUT, -o OUTPUT
                        Solution output as pickel
  --resolve RESOLVE, -r RESOLVE
                        Reinitialize with partial solution
```

After initializing the manager process, one can start multiple worker processes on the same or distributed machines. The worker process connects to the manager process via TCP/IP at the specified port (must be the same as the one the manager listens at) and obtains single work packages from the manager process. Using the flag -r (--resolve) one can specify a intermediate solution and refine or restart the solving process from there. The specified intermediate solution must be a pickled list of Solution objects.

```
usage: deimmunization_worker.py [-h] --input INPUT1 INPUT2 --masterip MASTERIP
                                --port PORT --authkey AUTHKEY --threads
                                THREADS

Rectangle Worker Grid implementation

Arguments:
  -h, --help            show this help message and exit
  --input INPUT1 INPUT2, -i INPUT INPUT
                        model files
  --masterip MASTERIP, -m MASTERIP
                        The IP of the master node
  --port PORT, -p PORT  port to connect
  --authkey AUTHKEY, -a AUTHKEY
                        authentication key
  --threads THREADS, -t THREADS
                        nof of core
```

INPUT1 and INPUT2 are single-objective files (in CPLEX compatible formats) of the bi-objective problems, each containing one of the objective function and all constraints of the bi-objective problem. 

### Examples:
Manager:
```
python deimmunization_manager.py -p 6882 -g 3 -a 0.009 -k rectangle -o ./example/FA8_HUMAN_hmmer_plm_n5_m50_f70_t01_g_r2188-2345_e20_k1_output.pcl
```
Worker:
```
python deimmunization_worker.py -i ./example/FA8_HUMAN_hmmer_plm_n5_m50_f70_t01_g_r2188-2345_e20_k1_pssm05_f01_f_he_experiment_local_imm.lp ./example/FA8_HUMAN_hmmer_plm_n5_m50_f70_t01_g_r2188-2345_e20_k1_pssm05_f01_f_he_experiment_local_en.lp -m 127.0.0.1 -p 6882 -a rectangle -t 4
```

## References

Please cite the following reference for the EVDeimmunization software:

Schubert B, Schärfe C, Dönnes P, Hopf T, Marks D, Kohlbacher O. Population-specific design of de-immunized protein biotherapeutics. PLoS Comput Biol. 2018;14(3):e1005983. Published 2018 Mar 2. doi:10.1371/journal.pcbi.1005983


Other sources:

[1] Lund O, Nielsen M, Kesmir C, et al. Definition of supertypes for HLA molecules using clustering of specificity matrices. Immunogenetics. 2004;55(12):797–810. doi:10.1007/s00251-004-0647-4

[2] Boland, N, Charkhgard, H, and Savelsbergh, M. A criterion space search algorithm for biobjective integer programming: The balanced box method. INFORMS J. Comput. 2015; 27, 735–754.
