# The code and data
## Constructing a simulation dataset
- make_simulation_data.py
## Hinge-annotated dataset
- pdb_12_with_hinges.csv 
- coord_csv.zip
## Determining the variance $\sigma^2$
- Generate simulated structures by $\sigma = 0.5, \dots, 5.0$ and determine $\sigma$ values
  - make_simulation_sigma.py
- Find protein pairs such that each protein pair has the same amino acid sequence
  - same_amino_acid_pdb.py
- Determine the $\sigma$ value
  - rmsd_for_same_amino_acid.py
## Hinge number accuracy and computation time
- Results on the simulation dataset
  - plot_accuracy_f_measure.py
- Results on the hinge-annotated dataset
  - eval_result_shibuya.py
## Experiments of previous methods
- DynDom
  - Simulation dataset: dyndom_simulation.py
  - The Shibuya 2008 dataset: dyndom_experiment_shibuya.py
- FATCAT
  - Simulation dataset: fatcat_simulation.py
  - The Shibuya 2008 dataset: fatcat_experiment_shibuya.py
## Data of experimental results
- original data of protein structures for generating simulated structures
  - selected_files.csv
- Results on the simulation dataset
  - AIC Exact: rmsdh_result/simulation_parallel_rmsdh_aic_0.894893output_simulation_file.csv
  - AIC LH: rmsdh_result/simulation_parallel_fast_rmsdh_aic_0.894893output_simulation_file.csv
  - BIC Exact: rmsdh_result/simulation_parallel_rmsdh_bic_0.894893output_simulation_file.csv
  - BIC LH: rmsdh_result/simulation_parallel_fast_rmsdh_bic_0.894893output_simulation_file.csv
  - SE: rmsdh_result/simulation_auto_shibuya_parallel.csv
  - DynDom: dyndom_simulation_hinge_count_result.csv, dyndom_simulation_result/dyndom_simulation_execution_time.csv
  - FATCAT: fatcat_simulation_hinge_count_result.csv, fatcat_simulation_execution_time.csv
- Results on the Shibuya 2008 dataset
  - AIC Exact: rmsdh_result/shibuya_rmsdh_aic_0.894893.csv
  - AIC LH: rmsdh_result/shibuya_fast_rmsdh_aic_0.894893.csv
  - BIC Exact: rmsdh_result/shibuya_rmsdh_bic_0.894893.csv
  - BIC LH: rmsdh_result/shibuya_fast_rmsdh_bic_0.894893.csv
  - SE: rmsdh_result/shibuya_auto_rmsdhk.csv
  - DynDom: all_pdb/dyndom_execution_time_shibuya.csv
  - FATCAT: all_pdb/fatcat_dyndom_execution_time_shibuya.csv
