#!/usr/bin/env python3

import pandas as pd

'''
Extraction pathways in which at least 20 entities are found in Reactome (pvalue < 0.05)
Return a dataframe with two columns: 'pathway' and 'entities' (which is a list of HGNC symbols)
'''
def extract_reactome_pathways(file):
    df = pd.read_csv(file, sep = ',', header = 0)
    # Filter for pvalue <0.05 and number of entities found in the model >20 for a specific pathway
    df = df[(df['Entities pValue'] < 0.05 ) & (df['#Entities found'] > 20)]
    # Combine pathway identifier and pathway name in a single column
    df['pathway'] = df[['Pathway identifier', 'Pathway name']].apply(lambda x: ':'.join(x), axis=1)
    # Create a list of entities within a column
    df['Submitted entities found'] = df['Submitted entities found'].str.split(';')
    # Renaming of columns to match 'get_pathways_and_nodes()' function
    df = df.rename(index=str, columns={"Submitted entities found": "entities"})
    df = df[['pathway', 'entities']]
    return df

'''
Extraction of KEGG pathways with at least 20 entities found and pvalue < 0.05
Returns a dataframe with two columns: 'pathway' and 'entities' (which is a list of HGNC symbols)
'''
def extract_KEGG_pathways(file, mapping_file):
    df = pd.read_csv(file, sep = '\t')
    df = df[(df['Category'] == 'KEGG_PATHWAY' ) & (df['PValue'] < 0.05) & (df['Count'] > 20)]
    df = df[['Term', 'Genes']]
    df['Genes'] = df['Genes'].str.split(', ')

    df_uniprot = pd.read_csv(mapping_file, sep = ';', header=0)
    df_uniprot = df_uniprot.rename(index = str, columns = {'yourlist:M201902226746803381A1F0E0DB47453E0216320D2AC3A4Q': 'HGNC.symbols'})
    df_uniprot[['HGNC.symbols','Entry']]

    hgnc = []
    for data in df.iterrows():
        # uniprot 'genes' in each row
        l_hgnc = []
        for genes in data[1]['Genes']:
            df2 = df_uniprot[df_uniprot['Entry'] == genes]
            for value in df2['HGNC.symbols'].str.split(',').values:
                for v in value:
                    l_hgnc.append(v)

        hgnc.append(l_hgnc)

    df['HGNC.symbols'] = hgnc
    df = df.rename(index = str, columns = {'Term':'pathway', 'HGNC.symbols':'entities'})
    df = df[['pathway', 'entities']]

    return df

def extract_acsn_pathways(pathway_file):
    df = pd.read_csv(pathway_file, sep = '\t', header = 0)
    # Only keep top level pathways - remove modules (+ the acsn2 master map is also removed)
    df = df[(df['Module '].str.contains(':master')) & (df['Module '] != 'acsn2:master ')]
    df = df[['Module ', 'Genes']]
    df['Genes'] = df['Genes'].str.split(' ')
    df = df.rename(index = str, columns = {'Module ':'pathway', 'Genes':'entities'})
    df = df[['pathway', 'entities']]

    return df

def write_pathways_and_nodes(df_pathways, file_hgnc, file_output):
    df_hgnc = pd.read_csv(file_hgnc, sep = '\t', header = 0)
    df_hgnc['HGNC.symbols'] = df_hgnc['HGNC.symbols'].str.split(', ')

    l_pathways = list(df_pathways[df_pathways.columns[0]])
    l_nodes = list(df_hgnc[df_hgnc.columns[0]])
    df_binaries = pd.DataFrame(columns=l_pathways, index = l_nodes)
    df_binaries = df_binaries.fillna(0)

    for data in df_hgnc.iterrows():
        if isinstance(data[1]['HGNC.symbols'], list):
            for node in data[1]['HGNC.symbols']:
                for row in df_pathways.iterrows():
                    if node in row[1]['entities']:
                        # Update values to 1 if node is present in pathway
                        df_binaries.loc[data[1]['Node'], row[1]['pathway']] = 1
    df_binaries.index.name = 'Node'

    # Write to file binary data (node present (1) or not (0) in the pathway)
    df_binaries.to_csv(file_output, sep  = "\t")


if __name__ == "__main__":
    # Mapping of HGNC ids to the nodes in the model
    file_hgnc_nodes = "../20190130_nodes_to_HGNC.txt"

    # Extraction of Reactome pathways binaries
    file_reactome = "./20190212_Reactome_result.csv"
    output_file_reactome = './Reactome_pathways_binary.txt'
    df_reactome = extract_reactome_pathways(file_reactome)
    write_pathways_and_nodes(df_reactome, file_hgnc_nodes, output_file_reactome)

    # Extraction of KEGG pathways binaries
    file_kegg = "../../../../taskNodeAssessment/results/DavidEnrichmentAnalysis/20190301_David_enrichment_analysis_result_allNodes.txt"
    # KEGG analysis with DAVID gives a mapping to uniprot ids. The following file contains a mapping between uniprot & hgnc
    file_uniprot_hgnc = "../../../../taskNodeAssessment/results/DavidEnrichmentAnalysis/HGNC_to_uniprot_accession.csv"
    df_kegg = extract_KEGG_pathways(file_kegg, file_uniprot_hgnc)
    write_pathways_and_nodes(df_kegg, file_hgnc_nodes, './KEGG_pathways_binary.txt')

    # Extraction of ACSN pathways binaries
    file_acsn = "./20190212_ACSN_results.tsv"
    df_acsn = extract_acsn_pathways(file_acsn)
    write_pathways_and_nodes(df_acsn, file_hgnc_nodes, './ACSN_pathways_binary.txt')
