#!/usr/bin/python
__author__="Michael Yang"
__purpose__="Conversion of titanCNA files into curation format"
_comments_="Script to convert TitanCNA output into pngs, then compile into curation format *******"

#Load modules
import os
import sys
import argparse
import pandas as pd
import subprocess

import contextlib

file_pattern1 = ".params.txt"
clusters = ['_cluster1', '_cluster2', '_cluster3']
ploidys = ['ploidy2', 'ploidy3', 'ploidy4']

@contextlib.contextmanager
def new_cd(x):
    d = os.getcwd()
    os.chdir(x)
    try:
        yield
    finally:
        os.chdir(d)


def run_curation_conversion(input_path):
    # command1="cd " + input_path + "/"
    # os.system(command1)
    print("     ///////////     ")
    print(input_path)
    print("     ///////////     ")
    # os.chdir(input_path)
    with new_cd(input_path):
        command= "for pdf_file in *.pdf; do convert ${pdf_file} ${pdf_file%%.*}.png; done;"
        os.system(command)
    # os.system(command2)
    # subprocess.run(command2, cwd=input_path)
    return 0;


def run_curation_compilation(input_path, output_path):
    curation_directory = output_path.rsplit("/", 1)[0]
    print("     ///////////     ")
    print(input_path)
    print(curation_directory)
    print(output_path)
    print("     ///////////     ")
    if not os.path.exists(curation_directory):
        os.makedirs(curation_directory)
    # command1="cd " + input_path
    # os.chdir(input_path)
    # os.system(command1)
    with new_cd(input_path):
        command="convert *_CNA.png *_LOH.png *chr1.png *chr2.png *chr3.png *chr4.png *chr5.png *chr6.png *chr7.png *chr8.png *chr9.png *chr10.png *chr11.png *chr12.png *chr13.png *chr14.png *chr15.png *chr16.png *chr17.png *chr18.png *chr19.png *chr20.png *chr21.png *chr22.png -quality 100 " + output_path
        os.system(command)
    # subprocess.run(command2, cwd=input_path)
    return 0;


def run_params_compilation(input_sample, input_path, ploidy, output_path):
    print("     ///////////     ")
    print(input_path)
    print(ploidy)
    print(output_path)
    print("     ///////////     ")
    input_list = list(input_sample)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    ploidy_val = ploidy.split('_')[1]
    samps = [x for x in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, x))]
    samps2 = []
    for x in samps:
        samp_name = x.rsplit('_cluster')[0]
        samps2.append(samp_name)
    # samps_working = [s for s in samps2 if input_sample in samps2]
    samps_working = list(filter(lambda x: any(input_sample in x for input_sample in input_list), samps2))
    print("     -------------      ")
    print(input_sample)
    print(samps_working)
    print("     -------------       ")
    for s in samps_working:
        tissue_ploidy = []
        tissue_purity = []
        tissue_loglikelihood = []
        tissue_phi = []
        tissue_cp = []
        tissue_sdbw_log = []
        tissue_sdbw_allele = []
        tissue_sdbw = []
        tissue_samps_clusters = []
        tissue_samp_barcode = []
        tissue_norm = []
        for c in clusters:
            # Reading in sample
            sample = s + c
            tissue_samp_barcode.append(s)
            tissue_samps_clusters.append(sample)
            file = input_path + '/{0}.params.txt'.format(sample)
            text_file = open(file, "r")

            #read whole file to a string
            data = text_file.read()
            
            #close file
            text_file.close()

            strs = data.split('\t')
            # Purity
            purity_data = (1 - pd.to_numeric(strs[1].split('\n')[0]))
            tissue_norm.append(pd.to_numeric(strs[1].split('\n')[0]))
            tissue_purity.append(purity_data)
            # Ploidy
            ploid_data = strs[2].split('\n')[0]
            tissue_ploidy.append(ploid_data)

            # Phi
            phi_val = cur_ploidy[-1]
            tissue_phi.append(phi_val)
            
            for n in range(0, len(strs)):
                cur = (strs[n].split('\n'))
                if len(cur) > 1:
                    # Log Likelihood
                    if cur[1] == 'Log likelihood:':
                        log_data = strs[(n+1)].split('\n')[0]
                        tissue_loglikelihood.append(log_data)
                    ## clonal cluster cp
                    if ('Clonal cluster cellular prevalence' in (cur[1])):
                        cp = ((strs[n])[-4] + '=' + (strs[n])[-2] + ':  ' + strs[n+1].split('\n')[0])
                        tissue_cp.append(cp)
                    # s_dbw_validity index (log)
                    if cur[1] == 'S_Dbw validity index (LogRatio):':
                        sdbw_data_log = strs[(n+1)].split('\n')[0]
                        tissue_sdbw_log.append(sdbw_data_log) 
                    # s_dbw_validity_index (allele)
                    if cur[1] == 'S_Dbw validity index (AllelicRatio):':
                        sdbw_data_allele = strs[(n+1)].split('\n')[0]
                        tissue_sdbw_allele.append(sdbw_data_allele)                
                    # s_dbw validity index (both)
                    if cur[1] == 'S_Dbw validity index (Both):':
                        # sdbw_data = strs[17].split('\n')[0]
                        sdbw_data = strs[(n+1)].split('\n')[0]
                        tissue_sdbw.append(sdbw_data)            
            # Concatenate                      
            t_samps = pd.DataFrame(tissue_samps_clusters, columns=['Sample'])
            t_barcode = pd.DataFrame(tissue_samp_barcode, columns=["Barcode"])
            t_ploids = pd.DataFrame(tissue_ploidy, columns=['Ploidy'])
            t_purity = pd.DataFrame(tissue_purity, columns=['Purity'])
            t_norm = pd.DataFrame(tissue_norm, columns=["Norm"])
            t_phi = pd.DataFrame(tissue_phi, columns=["Phi"])
            t_cp = pd.DataFrame(tissue_cp, columns=["Clonal_Cluster_Cellular_Prevalence"])
            t_log = pd.DataFrame(tissue_loglikelihood, columns=['Log_Likelihood'])
            t_sdbw_log = pd.DataFrame(tissue_sdbw_log, columns=['S_Dbw_Validity_Index_(LogRatio)'])
            t_sdbw_allele = pd.DataFrame(tissue_sdbw_allele, columns=['S_Dbw_Validity_Index_(AllelicRatio)'])
            t_sdbw = pd.DataFrame(tissue_sdbw, columns=['S_Dbw_Validity_Index_(Both)'])
            t_table = pd.concat([t_phi, t_samps, t_barcode, t_ploids, t_purity, t_norm, t_cp, t_log, t_sdbw_log, t_sdbw_allele, t_sdbw], axis=1)
            
            # Save as Params.txt - specific to samples PER ploidy
            t_table.to_csv(output_path+'{0}_{1}_collected_params.txt'.format(s, ploidy_val), sep='\t', index=False)
    # Save as Params.txt - for ALL ploidys together 


def run_params_samplewide_compilation(input_sample, input_path, params_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    input_list = list(input_sample)
    # ploidy_val = ploidy.split('_')[1]
    samps = [x for x in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, x))]
    samps2 = []
    for x in samps:
        samp_name = x.rsplit('_cluster')[0]
        samps2.append(samp_name)
    # samps_working = [s for s in samps2 if input_sample in samps2]
    samps_working = list(filter(lambda x: any(input_sample in x for input_sample in input_list), samps2))
    params_path_generalized = params_path.rsplit('_', 1)[0]
    print("     -------------      ")
    print(input_path)
    print(input_sample)
    print(samps_working)
    print(params_path)
    print(params_path_generalized)
    print(output_path)
    print("     -------------       ")
    num_nulls = 2
    spacing_df = [' ']
    null_samps_clusters = [' ']
    null_ploidy = [' ']
    null_purity = [' ']
    null_phi = [' ']
    null_cp = [' ']
    null_loglikelihood = [' ']
    null_sdbw_log = [' ']
    null_sdbw_allele = [' ']
    null_sdbw = [' ']
    null_norm = [' ']
    null_barcode = [' ']

    for n in range(0, num_nulls):
        n_samps = pd.DataFrame(null_samps_clusters, columns=['Sample'])
        n_barcode = pd.DataFrame(null_barcode, columns=["Barcode"])
        n_ploids = pd.DataFrame(null_ploidy, columns=['Ploidy'])
        n_purity = pd.DataFrame(null_purity, columns=['Purity'])
        n_norm = pd.DataFrame(null_norm, columns=["Norm"])
        n_phi = pd.DataFrame(null_phi, columns=["Phi"])
        n_cp = pd.DataFrame(null_cp, columns=["Clonal_Cluster_Cellular_Prevalence"])
        n_log = pd.DataFrame(null_loglikelihood, columns=['Log_Likelihood'])
        n_sdbw_log = pd.DataFrame(null_sdbw_log, columns=['S_Dbw_Validity_Index_(LogRatio)'])
        n_sdbw_allele = pd.DataFrame(null_sdbw_allele, columns=['S_Dbw_Validity_Index_(AllelicRatio)'])
        n_sdbw = pd.DataFrame(null_sdbw, columns=['S_Dbw_Validity_Index_(Both)'])
        n_table = pd.concat([n_phi, n_samps, n_barcode, n_ploids, n_purity, n_norm, n_cp, n_log, n_sdbw_log, n_sdbw_allele, n_sdbw], axis=1)
        if n == 0:
            n_table_concat = pd.DataFrame()
            n_table_concat = n_table
        else:
            n_table_concat = pd.concat([n_table_concat, n_table], axis=0)
    ## Concatenate PER sample
    for s in samps_working:
        cur_sample_table = pd.DataFrame()
        for p in ploidys:
            print('     /-/-/-/-/-/-/-/-/       ')
            print(params_path_generalized + '_' + p + '/{0}_{1}_collected_params.txt')
            print('     /-/-/-/-/-/-/-/-/       ')   
            params = pd.read_csv(params_path_generalized + '_' + p + '/{0}_{1}_collected_params.txt'.format(s, p), sep='\t')
            cur_sample_table = pd.concat([cur_sample_table, params, n_table_concat], axis=0)
        cur_sample_table.to_csv(output_path+'{0}_all_ploidy_params.txt'.format(s), sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process files.')
    overall_wd = os.getcwd()
    parser.add_argument('--input_params_file_path',help='please enter a path to the params file from titan',required=True)
    args = parser.parse_args()
    input_params_file_path = args.input_params_file_path        # Like /results/titan/hmm/titanCNA_ploidy__/{tumor}_cluster{clustNum}.params.txt

    ## Conversion from PNG to PDF
    sample_directory = overall_wd + "/" + input_params_file_path.split(file_pattern1)[0]
    run_curation_conversion(sample_directory)

    ## Compile a PDF to a new location
    sample_cluster_names = (input_params_file_path.split('/')[-1]).split(file_pattern1)[0]
    sample_names = (sample_cluster_names).split('_cluster')[0]
    print('---')
    print(sample_names)
    print('---')
    cur_ploidy = input_params_file_path.split('/')[3]
    out_pdf_file_name = sample_cluster_names + ".pdf"
    output_directory = overall_wd + "/results/CurationFiles/" + cur_ploidy + "/" + out_pdf_file_name
    run_curation_compilation(sample_directory, output_directory)

    ## Generate the Parameter Files as "collected_params.txt"
    params_directory = overall_wd + "/" + input_params_file_path.rsplit('/',1)[0]
    output_directory_params = overall_wd + '/results/CurationFiles/' + cur_ploidy + "/"
    run_params_compilation(sample_names, params_directory, cur_ploidy, output_directory_params)

    overall_params_outdir = overall_wd + '/results/CurationFiles/OverallSolutionParams/'
    ## Compile the Parameter Files that go into 'Overall_sample_params'
    run_params_samplewide_compilation(sample_names, params_directory, output_directory_params, overall_params_outdir)
