#!/usr/bin/env python3

import pandas as pd

def read_file(filename):
    df = pd.read_csv(filename, sep = '\t', header = 0)
    return df

'''Binary pathway = Reactome + KEGG + ACSN'''
def write_binary_pathway(reactome, kegg, acsn, file_out):
    concat = reactome
    concat = concat.merge(kegg, on = 'Node')
    concat = concat.merge(acsn, on = 'Node')
    concat.to_csv(file_out, sep = '\t', index = False)
    return concat


'''Binary all = Reactome + ACSN + KEGG + GO + COSMIC + Drug target'''
def write_binary_all(l_pathways, cosmic, go, drug_target, file_out):
    concat = l_pathways
    concat = concat.merge(cosmic, on = 'Node')
    concat = concat.merge(go, on = 'Node')
    concat = concat.merge(drug_target, on = 'Node')
    concat.to_csv(file_out, sep = '\t', index = False)

'''Numeric features = graph related features of each node of the model (thus drug target is not included)'''
def write_numeric_features(numeric_f, file_out):
    numeric_f = numeric_f.drop(columns = 'Drug.target')
    numeric_f.to_csv(file_out, sep = '\t', index = False)

if __name__ == "__main__":
    # List of files #
    file_node_features = '../../../data/features/node_features.txt'
    file_reactome = '../../../data/features/pathwayAnalysis/Reactome_pathways_binary.txt'
    file_acsn = '../../../data/features/pathwayAnalysis/ACSN_pathways_binary.txt'
    file_kegg = '../../../data/features/pathwayAnalysis/KEGG_pathways_binary.txt'
    file_go = '../../../data/features/Genes2Go/20190205_Gene2GO_node_list.txt'
    file_cosmic = '../../../data/features/COSMIC/20190206_COSMIC_node_census.tier_oncogene_tumor.suppressor.txt'

    # Get dataframes from files #
    df_node_features = read_file(file_node_features)
    df_reactome = read_file(file_reactome)
    df_acsn = read_file(file_acsn)
    df_kegg = read_file(file_kegg)
    df_go = read_file(file_go)
    df_cosmic = read_file(file_cosmic)

    # Export input files for Random Forest analysis #
    file_output_pathways = './input.data/binary_all_pathways.txt'
    df_pathways = write_binary_pathway(df_reactome, df_kegg, df_acsn, file_output_pathways)

    file_output_all_binaries = './input.data/binary_all.txt'
    write_binary_all(df_pathways, df_cosmic, df_go, df_node_features[['Node', 'Drug.target']], file_output_all_binaries)

    file_output_numeric_features = './input.data/numeric_features.txt'
    write_numeric_features(df_node_features, file_output_numeric_features)
