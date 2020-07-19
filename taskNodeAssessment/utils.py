#!/usr/bin/env python3

import biolqm
import seaborn as sns
import csv
import pandas as pd

'''Split a list l into lists of size n '''
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


'''Get number of perturbations '''
def get_nb_perturbations(perturbations):
    return len(perturbations)


'''Get number of drugs from the list of perturbations'''
def get_nb_drugs(perturbations):
    nb_drugs = 0
    for key in perturbations.keys():
        if '-' not in key:
            nb_drugs = nb_drugs + 1
    return nb_drugs


''' Generate dict of drug with nodes ko perturbations '''
def get_perturbations(comb_drugs, drug_targets):
    perturbations = {}
    for drug in comb_drugs:
        str_drug = ''
        for i in range(0, len(drug)):
            str_drug = str_drug + '-' + drug[i]
        str_drug = str_drug[1:]
        perturbation = []
        for d in drug_targets:
            if drug[0] == d[0]:
                for i in range(2, len(d), 1):
                    perturbation.append(d[i] + '%0')
            if len(drug) == 2:
                if drug[1] == d[0]:
                    for i in range(2, len(d), 1):
                        perturbation.append(d[i] + '%0')
        perturbations[str_drug] = perturbation
    return perturbations


def get_stable_states_wt(file, cell_model, perturbations):
    file.write("node_ko\tdrugs\tgrowth\tnb_fixpoints\tstable_state\n")
    for drug, ko in perturbations.items():
        l_nko = ''
        mutated_model = cell_model
        for node_ko in ko:
            mutated_model = biolqm.perturbation(mutated_model, node_ko)
            l_nko = l_nko + ',' + node_ko
        l_nko = l_nko[1:]

        '''calculate fixpoints'''
        fixpoints = biolqm.fixpoints(mutated_model)
        if (fixpoints):
            ''''get the mean of pro/anti-survival nodes in case of multiple fixpoints'''
            pro = 0
            anti = 0
            ss = ''
            for point in fixpoints:
                for key, val in point.items():
                    ss = ss + str(val)
                    if 'Prosurvival' in key:
                        pro = pro + val
                    if 'Antisurvival' in key:
                        anti = anti + val
                if (len(fixpoints) > 1):
                    ss = ss + ';'
            growth = (pro / len(fixpoints)) - (anti / len(fixpoints))
            file.write(l_nko + '\t' + str(drug) + '\t' + str(growth) + '\t' + str(len(fixpoints)) + '\t' + ss + '\n')
        else:
            '''No fixpoints'''
            file.write(l_nko + '\t' + str(drug) + '\tn/a\t0\t-\n')
    file.close()


def get_classification_wt(prediction, classification):
    model_pred = csv.reader(prediction, delimiter='\t')
    next(model_pred)
    model_pred = list(model_pred)

    list_classification = []
    TP = 0
    TN = 0
    FP = 0
    FN = 0

    for n in model_pred:
        if n[3] == 'TP':
            TP += 1
        elif n[3] == 'TN':
            TN += 1
        elif n[3] == 'FP':
            FP += 1
        elif n[3] == 'FN':
            FN += 1

    list_classification.append(TP)
    list_classification.append(TN)
    list_classification.append(FP)
    list_classification.append(FN)

    df = pd.DataFrame(columns=['TP', 'TN', 'FP', 'FN'])
    df.loc[-1] = list_classification
    df.index = df.index + 1  # shifting index
    df = df.sort_index()
    df.to_csv(classification, sep='\t')

def get_predictions(file_fxpts, file_prediction, hsa_synergy,  nb_perturbations, nb_drugs, phenotype, cell_name):
    prediction = open(file_prediction, 'w+')
    fxpts = pd.read_csv(file_fxpts, sep = '\t', header=0, keep_default_na=False)

    # WT doesn't have a mutated_node column
    if phenotype == 'WT':
        prediction.write('drugs\tprediction\texperimental\tclassification\n')
    else:
        prediction.write('mutated_node\tdrugs\tprediction\texperimental\tclassification\n')

    # Split the models into lists: one for each mutated node to be able to compute the drug synergies
    l_models = list(chunks(fxpts, nb_perturbations))

    # Go through each occurrence result of mutated nodes in the model
    for activity in l_models:
        # Create a tuple of triplets with all info (AK-BI, AK, BI)
        single = activity[0:nb_drugs]
        double = activity[nb_drugs:nb_perturbations]
        for d in double.iterrows():
            if phenotype == 'WT':
                str_triplet = d[1]['drugs']
            else:
                str_triplet = d[1]['mutated_node'] + '\t' + d[1]['drugs']
            triplet = []
            triplet.append(d[1])
            for s in single.iterrows():
                if s[1]['drugs'] in d[1]['drugs']:
                    triplet.append(s[1])

            exp_synergy = ''
            for val in hsa_synergy.iterrows():
                if val[1]['Combination'] == (triplet[0])['drugs']:
                    exp_synergy = str(val[1][cell_name]).upper()

            # Below are data where fixpoints are not found
            if ((triplet[0])['growth'] == 'n/a') and ((triplet[1])['growth'] == 'n/a') and ((triplet[2])['growth'] == 'n/a'):
                # Single and double drugs do not lead to a fixpoint
                prediction.write(str_triplet + '\tnone\t' + exp_synergy + '\tn/a\n')

            elif ((triplet[0])['growth'] == 'n/a') and ((triplet[1])['growth'] != 'n/a') and ((triplet[2])['growth'] != 'n/a'):
                # Single drugs finds stable state but not double drugs
                prediction.write(str_triplet + '\tsingle\t' + exp_synergy + '\tn/a\n')

            elif ((triplet[0])['growth'] != 'n/a') and ((triplet[1])['growth'] == 'n/a') and ((triplet[2])['growth'] == 'n/a'):
                # No stable states for single drugs - Stable state for combination
                prediction.write(str_triplet + '\tdouble\t' + exp_synergy + '\tn/a\n')

            elif ((triplet[0])['growth'] == 'n/a') or ((triplet[1])['growth'] == 'n/a') or ((triplet[2])['growth'] == 'n/a'):
                # Either one of single drug, double drugs do not lead to a stable state
                prediction.write(str_triplet + '\teither\t' + exp_synergy + '\tn/a\n')

            else:
                # Compute synergy and compare with HSA synergy data
                if (triplet[0])['growth'] < min(((triplet[1])['growth']), (triplet[2])['growth']):
                    if exp_synergy == 'FALSE':
                        prediction.write(str_triplet + '\tsynergy\t' + exp_synergy + '\tFP\n')
                    else:
                        prediction.write(str_triplet + '\tsynergy\t' + exp_synergy + '\tTP\n')
                else:
                    if exp_synergy == 'FALSE':
                        prediction.write(str_triplet + '\tno_synergy\t' + exp_synergy + '\tTN\n')
                    else:
                        prediction.write(str_triplet + '\tno_synergy\t' + exp_synergy + '\tFN\n')

    prediction.close()


def get_stable_states_mut(model, filename, mutations, ko_nodes):
    file = open(filename, "w+")
    file.write("mutated_node\tnode_ko\tdrugs\tgrowth\tnb_fixpoints\tstable_state\n")

    ''' For each inverted node, create a new mutated model'''
    for index, row in mutations.iterrows():
        mut = str(index) + '%' + str(row.values[0])

        '''Apply drug perturbations'''
        for drug, ko in ko_nodes.items():
            '''List of nodes knocked out by the drugs'''
            l_nodeko = ''

            '''create a model that will be perturbed and apply first the node activity inversion '''
            mutated_model = model
            mutated_model = biolqm.perturbation(mutated_model, mut)

            '''knock-out  the perturbed model drug targets one by one'''
            for node_ko in ko:
                mutated_model = biolqm.perturbation(mutated_model, node_ko)
                l_nodeko = l_nodeko + ',' + node_ko
            l_nodeko = l_nodeko[1:]

            '''calculate fixpoints'''
            fixpoints = biolqm.fixpoints(mutated_model)
            if (fixpoints):
                '''ss: variable to keep the stable states'''
                ss = ''
                ''''get the mean of pro/anti-survival nodes in case of multiple fixpoints'''
                pro = 0
                anti = 0
                for point in fixpoints:
                    for key, val in point.items():
                        if 'Prosurvival' in key:
                            pro = pro + val
                        if 'Antisurvival' in key:
                            anti = anti + val
                        ss = ss + str(val)

                    if (len(fixpoints) > 1):
                        ss = ss + ';'
                growth = (pro / len(fixpoints)) - (anti / len(fixpoints))
                file.write(mut + '\t' + l_nodeko + '\t' + str(drug) + '\t' + str(growth) + '\t' + str(
                    len(fixpoints)) + '\t' + ss + '\n')
            else:
                '''No fixpoints'''
                file.write(mut + '\t' + l_nodeko + '\t' + str(drug) + '\tn/a\t0\t-\n')
    file.close()


def get_classification_gain_loss(file_mutated_pred, file_WT_pred, nb_double_comb):
    # Open and read the mutated (fixed or inverted) and WT files - omit first line (header)
    m = open(file_mutated_pred, 'r')
    m = csv.reader(m, delimiter='\t')
    next(m)
    m = list(m)
    wt_m = open(file_WT_pred, 'r')
    wt_m = csv.reader(wt_m, delimiter='\t')
    next(wt_m)
    wt_m = list(wt_m)

    # The data is split into collected lines of the nb of drug combinations tested
    chunked_m = list(chunks(m, nb_double_comb))
    l_nodes = []
    l_loss_TP = []
    l_loss_TN = []
    l_gain_TP = []
    l_gain_TN = []
    l_na = []

    for line in chunked_m:
        l_nodes.append(line[0][0])
        loss_TP = 0;
        gain_TP = 0;
        loss_TN = 0;
        gain_TN = 0;
        na = 0;

        for i in range(0, len(wt_m), 1):
            if line[i][4] == wt_m[i][3]:
                # Same predictions with WT and mutant
                '''print(line[i] , wt_m[i])'''
            else:
                if wt_m[i][3] == 'TP' and (line[i][4] == 'FN' or line[i][4] == 'n/a'):
                    loss_TP = loss_TP + 1
                else:
                    if wt_m[i][3] == 'TN' and (line[i][4] == 'FP' or line[i][4] == 'n/a'):
                        loss_TN = loss_TN + 1
                    else:
                        if (wt_m[i][3] == 'FN' or wt_m[i][3] == 'n/a') and line[i][4] == 'TP':
                            gain_TP = gain_TP + 1
                        else:
                            if (wt_m[i][3] == 'FP' or wt_m[i][3] == 'n/a') and line[i][4] == 'TN':
                                gain_TN = gain_TN + 1
                            else:
                                if (line[i][4] == 'n/a'):
                                    na = na + 1
        l_loss_TP.append(loss_TP)
        l_loss_TN.append(loss_TN)
        l_gain_TP.append(gain_TP)
        l_gain_TN.append(gain_TN)
        l_na.append(na)
    # Create and fill dataframe with computed data
    df = pd.DataFrame(index=l_nodes, columns=['loss_TP', 'loss_TN', 'gain_TP', 'gain_TN', 'na'])
    df['loss_TP'] = l_loss_TP
    df['loss_TN'] = l_loss_TN
    df['gain_TP'] = l_gain_TP
    df['gain_TN'] = l_gain_TN
    df['na'] = l_na
    # Export dataframe as CSV
    file_output = file_mutated_pred.replace('prediction', 'classification_gainloss')
    df.to_csv(file_output, sep='\t')
    return file_output


def get_classification_heatmap(file_model, cmap):
    df = pd.read_csv(file_model, sep='\t', header=0, index_col = 0)
    df = (df - df.mean()) / df.std()
    df = df.fillna(value=0)
    g = sns.clustermap(df, yticklabels=True, cmap=cmap, figsize=(20,25))

    file_output = file_model.replace('.txt', '_heatmap.png')
    file_output = file_output.replace('results/', 'results/figures/')
    g.savefig(file_output)
