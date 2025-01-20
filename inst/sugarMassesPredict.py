#!/usr/bin/env python
# import modules
import pandas as pd
import numpy as np
import itertools
import math
import re
from operator import itemgetter
import os
# suppress warnings
pd.options.mode.chained_assignment = None  # default='warn'

#possible modifications
possible_modifications = ['carboxylicacid',
                          'sialicacid',
                          'phosphate',
                          'deoxy',
                          'nacetyl',
                          'omethyl',
                          'anhydrobridge',
                          'oacetyl',
                          'unsaturated',
                          'alditol',
                          'amino',
                          'dehydrated',
                          'sulfate']
# hexose and water masses to build molecule base
hex_mass = 180.0633881
water_mass = 18.01056468
# mass differences for modifications
pent_mdiff = -30.01056468
modifications_mdiff = {
    "sulfate": 79.95681486,
    "anhydrobridge": -water_mass,
    "omethyl": 14.01565006,
    "carboxylicacid": 13.97926456,
    "sialicacid" : 129.0425931,
    "nacetyl": 41.0265491,
    "oacetyl": 42.01056468 ,
    "phosphate": 79.96633052 ,
    "deoxy": -15.99491462,
    "unsaturated": -water_mass,
    "alditol": 2.015650064,
    "amino": -0.984015588,
    "dehydrated": -water_mass}

# mass differences for labels
proca_mdiff = 219.173546
pa_mdiff = 78.058183
aba_mdiff = 121.052764
aba_nonreductive_mdiff = 119.037114
pmp_mdiff = 330.148061
ab_mdiff = 120.068748
aq_nonreductive_mdiff = 126.058183
aminopentyllinker_mdiff= 85.08914935

#possible names for labels
proa_names = {"procainamide", "proca", "procA", "ProA"}
pa_names = {"2-ap", "2-AP", "pa", "PA", "2-aminopyridine"}
aba_names = {"2-aa", "2-AA", "aba", "ABA", "2-aminobenzoic acid"}
aba_nonreductive_names = {"2-aa nonreductive", "2-AA nonreductive", "aba nonreductive", "ABA nonreductive", "2-aminobenzoic acid nonreductive"}
ab_names = {"2-ab", "2-AB", "ab", "AB", "2-aminobenzamide"}
pmp_names = {"pmp", "PMP", "1-phenyl-3-methyl-5-pyrazolone"}
aq_nonreductive_names = {'3AQ nonreductive', '3-AQ nonreductive', '3-aminoquinoline nonreductive'}
aminopentyllinker_names = {"aminopentyllinker"}


# mass differences for ions
ion_mdiff = {
    "H": 1.007825032,
    "Na": 22.98976928,
    "Cl": 34.96885268,
    "CHOO": 44.997654272,
    "NH4": 18.034374128,
    "K": 38.96370668,
    "Ca": 39.96259098
}

e_mdiff = 0.00054857990924

# formulas

#chnosp
formulas = {
    "hex": [6, 12, 0, 6, 0, 0],
    "pent": [5, 10, 0, 5, 0, 0],
    "water": [0, -2, 0, -1, 0, 0],
    "sulfate": [0, 0, 0, 3, 1, 0],
    "anhydrobridge": [0, -2, 0, -1, 0, 0],
    "omethyl": [1, 2, 0, 0, 0, 0],
    "carboxylicacid": [0, -2, 0, 1, 0, 0],
    "sialicacid": [11, 19, 1, 9, 0, 0],
    "nacetyl": [2, 3, 1, 0, 0, 0],
    "oacetyl": [2, 2, 0, 1, 0, 0],
    "phosphate": [0, 1, 0, 3, 0, 1],
    "deoxy": [6, 12, 0, 5, 0, 0],
    "proca": [13, 21, 3, 0, 0, 0],
    "unsaturated":  [0, -2, 0, -1, 0, 0],
    "alditol": [0, +2, 0, 0, 0, 0],
    "amino": [0, +1, +1, -1, 0, 0],
    "dehydrated": [0, -2, 0, -1, 0, 0],
    "pa": [5, 6, 2, -1, 0, 0] ,
    "aba": [7, 7, 1, 1, 0, 0] ,
    "aba_nonreductive": [7, 7, 1, 1, 0, 0],
    "pmp": [20, 18, 4, 1, 0, 0],
    "ab": [7, 8, 2, 0, 0, 0],
    "aq_nonreductive" : [9, 6, -1, 2, 0, 0],
    "aminopentyllinker": [5, 11, 1, 0, 0, 0]
}
# modification types
modifications_anionic = {"sulfate",
                         "phosphate",
                         "carboxylicacid",
                         "sialicacid"}
modifications_neutral = {"anhydrobridge",
                         "omethyl",
                         "nacetyl",
                         "oacetyl",
                         "deoxy",
                         "unsaturated",
                         "amino",
                         "dehydrated"}

#names
names_iupac = {
    "hex" : 'Hex',
    'pent' : 'Pen',
    "sulfate": 'Sulfate',
    "anhydrobridge": 'AnhydroBridge',
    "omethyl": 'O-Methyl',
    "carboxylicacid": 'CarboxylicAcid',
    "sialicacid": "NeuAc",
    "nacetyl": 'N-Acetyl',
    "oacetyl": 'O-Acetyl',
    "phosphate": 'Phosphate',
    "deoxy": 'DeoxyHex',
    "unsaturated": 'Unsaturated',
    "alditol": 'Alditol',
    "amino": 'Amino',
    "dehydrated": 'Dehydrated'}
names_glycoct = {
    "hex": 'HEX',
    'pent': 'PEN',
    "sulfate": 'SO4',
    "anhydrobridge": 'ANH',
    "omethyl": 'OMe',
    "carboxylicacid": 'COOH',
    "sialicacid": "SIA",
    "nacetyl": 'NAc',
    "oacetyl": 'Ac',
    "phosphate": 'PO4',
    "deoxy": 'DHEX',
    "unsaturated": 'UNS',
    "alditol": 'ALD',
    "amino": 'NH2',
    "dehydrated": 'Y'}
names_oxford = {
    "hex": 'H',
    'pent': 'P',
    "sulfate": 'S',
    "anhydrobridge": 'B',
    "omethyl": 'M',
    "carboxylicacid": 'A',
    "sialicacid": 'SA',
    "nacetyl": 'N',
    "oacetyl": 'Ac',
    "phosphate": 'P',
    "deoxy": 'D',
    "unsaturated": 'U',
    "alditol": 'o',
    "amino": 'Am',
    "dehydrated": 'Y'}

#brackets
bracket_mapping = {
    "hex": '()',
    'pent': '()',
    "sulfate": '[]',
    "anhydrobridge": '[]',
    "omethyl": '[]',
    "carboxylicacid": '[]',
    "nacetyl": '[]',
    "oacetyl": '[]',
    "phosphate": '[]',
    "deoxy": '()',
    "sialicacid": '()',
    "unsaturated": '[]',
    "alditol": '[]',
    "amino": '[]',
    "dehydrated": '[]'}

#n-glycan limits
nglycan_limits = {
    "hex": 20,
    "nacetyl": 20,
    "deoxy": 6,
    "sialicacid": 5,
    "pent": 4,
    "sulfate": 3,
    "phosphate": 2,
    "amino": 0
}
oglycan_limits = {
    "hex": 14,
    "nacetyl": 14,
    "deoxy": 6,
    "sialicacid": 7,
    "pen": 3,
    "sulfate": 6,
    "phosphate": 6,
    "amino": 2
}

def predict_sugars(dp= [1, 6], polarity='neg', scan_range=[175, 1400], pent_option=False, modifications='none', nmod_max=1, double_sulfate=False, label='none', ion_type = "ESI", format="long", adducts = ["all"], naming = "IUPAC", glycan_linkage = ["none"], modification_limits = 'none', custom_label_name = 'none', custom_label_mdiff = 'none', custom_label_formdiff = [0,0,0,0,0,0]):
    if 'all' in adducts:
        adducts=['H', 'Cl', 'CHOO', 'Na', 'NH4', 'K']
    if type(adducts)==str:
        adducts = [adducts]
    dp_range_list = list(range(dp[0], dp[1] + 1))
    if "all" in modifications:
        modifications = possible_modifications
    if "sulfate" in modifications and len(modifications) > 1:
        modifications.append(modifications.pop(modifications.index('sulfate')))
    if "alditol" in modifications:
        alditol_option = 'y'
        modifications.remove("alditol")
    elif "alditol" not in modifications:
        alditol_option = 'n'
    if "unsaturated" in modifications:
        unsaturated_option = 'y'
        modifications.remove("unsaturated")
    elif "unsaturated" not in modifications:
        unsaturated_option = 'n'
    if "dehydrated" in modifications:
        dehydrated_option = 'y'
        modifications.remove("dehydrated")
    elif "dehydrated" not in modifications:
        dehydrated_option = 'n'
    #calculate possible masses
    # build hexose molecules
    def getHexMasses(dp_range_list):
        dp = np.array(dp_range_list, dtype=int)
        mass = dp * hex_mass - (dp - 1) * water_mass
        masses = np.column_stack((dp, dp, mass))
        return masses
    masses = getHexMasses(dp_range_list)
    #define additional functions
    def dpRepeats(dp_range_list):
        for i in dp_range_list:
            index = dp_range_list.index(i)
            if index == 0 : repeats_list = np.arange(i+1)
            if index > 0: repeats_list = np.concatenate((repeats_list,np.arange(i+1)))
        return repeats_list
    def getPentMasses(masses):
        dp = np.repeat(masses[:,1], masses[:,1].astype(int)+1, axis = 0)
        pent = dpRepeats(dp_range_list)
        hex = dp - pent
        mass = np.repeat(masses[:,2], masses[:,1].astype(int)+1, axis = 0)
        mass = mass + pent * pent_mdiff
        masses = np.column_stack((dp, hex, pent, mass))
        return masses
    def getModificationNumbers(dp_range_list, m, pent_option, modifications):
        for i in dp_range_list:
            index = dp_range_list.index(i)
            a = list(range(0, i + 1))
            if pent_option == True:
                b = np.array(a)[np.rollaxis(np.indices((len(a),) * m), 0, m + 1).reshape(-1, m)]
                b = np.tile(b,(i+1,1))
                if index == 0: modification_numbers = b
                if index > 0: modification_numbers = np.concatenate((modification_numbers, b), axis=0)
            if pent_option == False:
                b = np.array(a)[np.rollaxis(np.indices((len(a),) * m), 0, m + 1).reshape(-1, m)]
                if index == 0: modification_numbers = b
                if index > 0: modification_numbers = np.concatenate((modification_numbers, b), axis = 0)
        return modification_numbers
    if pent_option == 1:
        masses = getPentMasses(masses)
    #add modifications
    if "none" not in modifications and pent_option == True and len(modifications) != 0:
        m = len(modifications)
        modification_numbers = getModificationNumbers(dp_range_list, m, pent_option, modifications)
        masses = np.repeat(masses, repeats=(masses[:,0].astype(int)+1) ** m, axis=0)
        masses = np.column_stack((masses, modification_numbers))
        nmod_sums =  np.sum(modification_numbers, axis = 1)
        nmod_avg = nmod_sums / masses[:,0]
        nmod_avg_filter = nmod_avg <= nmod_max
        masses = masses[nmod_avg_filter]
        del modification_numbers, nmod_sums, nmod_avg, nmod_avg_filter
        masses_columns = ['dp', 'hex', 'pent', 'mass'] + modifications
        #subtract deoxy from hex
        if "deoxy" in modifications:
            # remove rows with deoxypentose
            hex_index = masses_columns.index('hex')
            deoxy_index = masses_columns.index('deoxy')
            masses = masses[masses[:, deoxy_index] <= masses[:, hex_index]]
            masses[:,hex_index] = masses[:,hex_index] - masses[:,deoxy_index]
        # subtract sialicacid from hex
        if "sialicacid" in modifications:
            # remove rows with sialicacid pentose
            hex_index = masses_columns.index('hex')
            sialicacid_index = masses_columns.index('sialicacid')
            masses = masses[masses[:, sialicacid_index] <= masses[:, hex_index]]
            masses[:,hex_index] = masses[:,hex_index] - masses[:,sialicacid_index]
        masses = pd.DataFrame(masses, columns=masses_columns)
        if "nglycan" in glycan_linkage:
            masses = masses[masses['hex'] != 0]
            relevant_columns = list(set(masses.columns) & set(nglycan_limits.keys()))
            if "nacetyl" in modifications:
                masses['hex'] = masses['hex'] - masses['nacetyl']
            conditions = masses[relevant_columns].apply(lambda col: col.le(nglycan_limits.get(col.name, float('inf'))), axis=0)
            masses = masses[conditions.all(axis=1)]
            if "nacetyl" in modifications:
                masses['hex'] = masses['hex'] + masses['nacetyl']
            if "sulfate" in modifications and "phosphate" in modifications:
                masses = masses[(masses['sulfate'] == 0) | (masses['phosphate'] == 0)]
            if "deoxy" in modifications:
                masses = masses[masses['deoxy']+1 <= masses['hex']]
            if "sialicacid" in modifications and "nacetyl" in modifications:
                to_remove = masses[(masses['nacetyl'] <= 2) & ((masses['hex'] - masses['nacetyl']) > 2) & masses['sialicacid'] != 0]
                masses.drop(to_remove.index, axis=0, inplace=True)
                del to_remove
        if "oglycan" in glycan_linkage:
            relevant_columns = list(set(masses.columns) & set(oglycan_limits.keys()))
            conditions = masses[relevant_columns].apply(lambda col: col.le(oglycan_limits.get(col.name, float('inf'))),
                                                        axis=0)
            masses = masses[conditions.all(axis=1)]
        modification_masses = masses.apply(
            lambda row: sum(row[col] * modifications_mdiff[col] for col in modifications), axis=1)
        masses['mass'] += modification_masses
        del modification_masses
    if "none" not in modifications and pent_option == False and len(modifications) != 0:
        m = len(modifications)
        modification_numbers = getModificationNumbers(dp_range_list, m, pent_option, modifications)
        masses = np.repeat(masses, repeats=(masses[:,0].astype(int)+1) ** m, axis=0)
        masses = np.column_stack((masses, modification_numbers))
        nmod_sums =  np.sum(modification_numbers, axis = 1)
        nmod_avg = nmod_sums / masses[:,0]
        nmod_avg_filter = nmod_avg <= nmod_max
        masses = masses[nmod_avg_filter]
        del modification_numbers, nmod_sums, nmod_avg, nmod_avg_filter
        masses_columns = ['dp', 'hex', 'mass'] + modifications
        #subtract deoxy from hex
        if "deoxy" in modifications:
            hex_index = masses_columns.index('hex')
            deoxy_index = masses_columns.index('deoxy')
            masses[:,hex_index] = masses[:,hex_index] - masses[:,deoxy_index]
        # subtract sialicacid from hex
        if "sialicacid" in modifications:
            # remove rows with sialicacid pentose
            hex_index = masses_columns.index('hex')
            sialicacid_index = masses_columns.index('sialicacid')
            masses[:,hex_index] = masses[:,hex_index] - masses[:,sialicacid_index]
        masses = pd.DataFrame(masses, columns=masses_columns)
        if "nglycan" in glycan_linkage:
            masses = masses[masses['hex'] != 0]
            relevant_columns = list(set(masses.columns) & set(nglycan_limits.keys()))
            if "nacetyl" in modifications:
                masses['hex'] = masses['hex'] - masses['nacetyl']
            conditions = masses[relevant_columns].apply(lambda col: col.le(nglycan_limits.get(col.name, float('inf'))), axis=0)
            masses = masses[conditions.all(axis=1)]
            if "nacetyl" in modifications:
                masses['hex'] = masses['hex'] + masses['nacetyl']
            if "sulfate" in modifications and "phosphate" in modifications:
                masses = masses[(masses['sulfate'] == 0) | (masses['phosphate'] == 0)]
            if "deoxy" in modifications:
                masses = masses[masses['deoxy']+1 <= masses['hex']]
            if "sialicacid" in modifications and "nacetyl" in modifications:
                to_remove = masses[(masses['nacetyl'] <= 2) & ((masses['hex'] - masses['nacetyl']) > 2) & masses['sialicacid'] != 0]
                masses.drop(to_remove.index, axis=0, inplace=True)
                del to_remove
        if "oglycan" in glycan_linkage:
            relevant_columns = list(set(masses.columns) & set(oglycan_limits.keys()))
            conditions = masses[relevant_columns].apply(lambda col: col.le(oglycan_limits.get(col.name, float('inf'))),
                                                        axis=0)
            masses = masses[conditions.all(axis=1)]
        modification_masses = masses.apply(
            lambda row: sum(row[col] * modifications_mdiff[col] for col in modifications), axis=1)
        masses['mass'] += modification_masses
        del modification_masses
    if "none" in modifications or len(modifications) == 0:
        if pent_option == True: masses = pd.DataFrame(masses, columns=['dp', 'hex', 'pent', 'mass'])
        if pent_option == False: masses = pd.DataFrame(masses, columns=['dp', 'hex', 'mass'])
    if "sulfate" in modifications and double_sulfate == True:
        #print("--> adding extra sulfate groups")
        masses_s1 = masses.loc[masses['sulfate'] >= 1]
        masses_s2 = masses_s1
        masses_s2.sulfate = masses_s1.sulfate + masses_s1.dp
        masses_s2.mass = masses_s1.mass + masses_s1.dp * modifications_mdiff['sulfate']
        masses = pd.concat([masses, masses_s2]).reset_index()
        del masses_s1
        del masses_s2
    if label in proa_names:
        masses['mass'] = masses.mass + proca_mdiff
    if label in pa_names:
        masses['mass'] = masses.mass + pa_mdiff
    if label in aba_names:
        masses['mass'] = masses.mass + aba_mdiff
    if label in ab_names:
        masses['mass'] = masses.mass + ab_mdiff
    if label in pmp_names:
        masses['mass'] = masses.mass + pmp_mdiff
    if label in aba_nonreductive_names:
        masses['mass'] = masses.mass + aba_nonreductive_mdiff
    if label in aq_nonreductive_names:
        masses['mass'] = masses.mass + aq_nonreductive_mdiff
    if label in aminopentyllinker_names:
        masses['mass'] = masses.mass + aminopentyllinker_mdiff
    if custom_label_name != 'none':
        masses['mass'] = masses.mass + custom_label_mdiff
    if unsaturated_option == 'y':
        masses_a = masses.copy()
        masses_a['unsaturated'] = 1
        masses['unsaturated'] = 0
        masses_a.mass = masses.mass + modifications_mdiff['unsaturated']
        masses = pd.concat([masses, masses_a]).reset_index(drop=True)
        del masses_a
    if alditol_option == 'y':
        masses_a = masses.copy()
        masses_a['alditol'] = 1
        masses['alditol'] = 0
        masses_a.mass = masses_a.mass + modifications_mdiff['alditol']
        masses = pd.concat([masses, masses_a]).reset_index(drop=True)
        del masses_a
    if dehydrated_option == 'y':
        #print("--> adding dehydration to sugars")
        masses_a = masses.copy()
        masses_a['dehydrated'] = 1
        masses['dehydrated'] = 0
        masses_a.mass = masses_a.mass + modifications_mdiff['dehydrated']
        masses = pd.concat([masses, masses_a]).reset_index(drop=True)
        del masses_a
    #remove modification combinations that are not possible
    if "none" not in modifications and pent_option == True and len(modifications) != 0:
        if 'carboxylicacid' in modifications:
            if 'phosphate' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['carboxylicacid'] + masses['phosphate'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['carboxylicacid'] + masses['phosphate'] <= masses['hex'] + masses['pent'])]
            if 'anhydrobridge' in modifications:
                masses = masses[(masses['carboxylicacid'] + masses['anhydrobridge'] <= masses['hex'] + masses['pent'])]
            if 'deoxy' in modifications:
                masses = masses[(masses['carboxylicacid'] <= masses['hex'] + masses['pent'])]
        if 'anhydrobridge' in modifications:
            if 'phosphate' in modifications:
                masses = masses[(masses['anhydrobridge'] + masses['phosphate'] <= masses['hex'] )]
            if unsaturated_option == 'y':
                masses = masses[(masses['unsaturated'] + masses['anhydrobridge'] <= masses['hex'] + masses['pent'])]
            if 'oacetyl' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['anhydrobridge'] + masses['oacetyl'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['anhydrobridge'] + masses['oacetyl'] <= masses['hex'] + masses['pent'])]
            if 'deoxy' in modifications:
                masses = masses[(masses['anhydrobridge'] <= masses['hex'])]
            if 'amino' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['anhydrobridge'] + masses['amino'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['anhydrobridge'] + masses['amino'] <= masses['hex'] + masses['pent'])]
            if dehydrated_option == 'y':
                if 'deoxy' in modifications:
                    masses = masses[(masses['anhydrobridge'] + masses['dehydrated'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['anhydrobridge'] + masses['dehydrated'] <= masses['hex'] + masses['pent'])]
            if alditol_option == 'y':
                masses = masses[(masses['anhydrobridge'] + masses['alditol'] <= masses['hex'] + masses['pent'])]
        if 'oacetyl' in modifications:
            if 'amino' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['oacetyl'] + masses['amino'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['oacetyl'] + masses['amino'] <= masses['hex'] + masses['pent'])]
        if unsaturated_option == 'y':
            if 'nacetyl' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['unsaturated'] + masses['nacetyl'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['unsaturated'] + masses['nacetyl'] <= masses['hex'] + masses['pent'])]
            if 'amino' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['unsaturated'] + masses['amino'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['unsaturated'] + masses['amino'] <= masses['hex'] + masses['pent'])]
            if dehydrated_option == 'y':
                if 'deoxy' in modifications:
                    masses = masses[(masses['unsaturated'] + masses['dehydrated'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['unsaturated'] + masses['dehydrated'] <= masses['hex'] + masses['pent'])]
            if alditol_option == 'y':
                masses = masses[(masses['unsaturated'] + masses['alditol'] <= masses['hex'] + masses['pent'])]
        if dehydrated_option == 'y':
            if 'nacetyl' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['dehydrated'] + masses['nacetyl'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['dehydrated'] + masses['nacetyl'] <= masses['hex'] + masses['pent'])]
            if 'amino' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['dehydrated'] + masses['amino'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['dehydrated'] + masses['amino'] <= masses['hex'] + masses['pent'])]
        if alditol_option == 'y':
            if 'deoxy' in modifications:
                masses = masses[(masses['alditol'] <= masses['hex'] + masses['pent'])]
            if dehydrated_option == 'y':
                if 'deoxy' in modifications:
                    masses = masses[(masses['alditol'] + masses['dehydrated'] <= masses['hex'] + masses['pent'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['alditol'] + masses['dehydrated'] <= masses['hex'] + masses['pent'])]
    if "none" not in modifications and pent_option == False and len(modifications) != 0:
        if 'carboxylicacid' in modifications:
            if 'phosphate' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['carboxylicacid'] + masses['phosphate'] <= masses['hex'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['carboxylicacid'] + masses['phosphate'] <= masses['hex'])]
            if 'anhydrobridge' in modifications:
                masses = masses[(masses['carboxylicacid'] + masses['anhydrobridge'] <= masses['hex'] )]
            if 'deoxy' in modifications:
                masses = masses[(masses['carboxylicacid'] <= masses['hex'] )]
        if 'anhydrobridge' in modifications:
            if 'phosphate' in modifications:
                masses = masses[(masses['anhydrobridge'] + masses['phosphate'] <= masses['hex'] )]
            if unsaturated_option == 'y':
                masses = masses[(masses['unsaturated'] + masses['anhydrobridge'] <= masses['hex'] )]
            if 'oacetyl' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['anhydrobridge'] + masses['oacetyl'] <= masses['hex'] +  masses['deoxy'])]
                else:
                    masses = masses[(masses['anhydrobridge'] + masses['oacetyl'] <= masses['hex'] )]
            if 'deoxy' in modifications:
                masses = masses[(masses['anhydrobridge'] <= masses['hex'])]
            if 'amino' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['anhydrobridge'] + masses['amino'] <= masses['hex']  + masses['deoxy'])]
                else:
                    masses = masses[(masses['anhydrobridge'] + masses['amino'] <= masses['hex'] )]
            if dehydrated_option == 'y':
                if 'deoxy' in modifications:
                    masses = masses[(masses['anhydrobridge'] + masses['dehydrated'] <= masses['hex'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['anhydrobridge'] + masses['dehydrated'] <= masses['hex'] )]
            if alditol_option == 'y':
                masses = masses[(masses['anhydrobridge'] + masses['alditol'] <= masses['hex'] )]
        if 'oacetyl' in modifications:
            if 'amino' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['oacetyl'] + masses['amino'] <= masses['hex'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['oacetyl'] + masses['amino'] <= masses['hex'] )]
        if unsaturated_option == 'y':
            if 'nacetyl' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['unsaturated'] + masses['nacetyl'] <= masses['hex'] ++ masses['deoxy'])]
                else:
                    masses = masses[(masses['unsaturated'] + masses['nacetyl'] <= masses['hex'] )]
            if 'amino' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['unsaturated'] + masses['amino'] <= masses['hex'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['unsaturated'] + masses['amino'] <= masses['hex'] )]
            if dehydrated_option == 'y':
                if 'deoxy' in modifications:
                    masses = masses[(masses['unsaturated'] + masses['dehydrated'] <= masses['hex'] + masses['deoxy'])]
                else:
                    masses = masses[(masses['unsaturated'] + masses['dehydrated'] <= masses['hex'])]
            if alditol_option == 'y':
                masses = masses[(masses['unsaturated'] + masses['alditol'] <= masses['hex'])]
        if dehydrated_option == 'y':
            if 'nacetyl' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['dehydrated'] + masses['nacetyl'] <= masses['hex'] + + masses['deoxy'])]
                else:
                    masses = masses[(masses['dehydrated'] + masses['nacetyl'] <= masses['hex'] )]
            if 'amino' in modifications:
                if 'deoxy' in modifications:
                    masses = masses[(masses['dehydrated'] + masses['amino'] <= masses['hex']  + masses['deoxy'])]
                else:
                    masses = masses[(masses['dehydrated'] + masses['amino'] <= masses['hex'] )]
        if alditol_option == 'y':
            if 'deoxy' in modifications:
                masses = masses[(masses['alditol'] <= masses['hex'] )]
            if dehydrated_option == 'y':
                if 'deoxy' in modifications:
                    masses = masses[(masses['alditol'] + masses['dehydrated'] <= masses['hex'] +  masses['deoxy'])]
                else:
                    masses = masses[(masses['alditol'] + masses['dehydrated'] <= masses['hex'] )]
    if modification_limits != 'none':
        conditions = masses[list(modification_limits.keys())].apply(
            lambda col: col.le(modification_limits.get(col.name, float('inf'))), axis=0)
        masses = masses[conditions.all(axis=1)]
    if 'sialicacid' in modifications:
        if 'anhydrobridge' in modifications:
            masses = masses[(masses['anhydrobridge'] <= masses['hex'])]
        if pent_option == False:
            if 'carboxylicacid' in modifications:
                masses = masses[(masses['carboxylicacid'] <= masses['hex'])]
            if alditol_option == 'y':
                masses = masses[(masses['alditol'] <= masses['hex'])]
            if 'deoxy' in modifications:
                if 'sulfate' in modifications:
                    if double_sulfate == True:
                        masses = masses[(masses['sulfate'] <= (masses['hex'] + masses['deoxy']) * 2)]
                    if double_sulfate == False:
                        masses = masses[(masses['sulfate'] <= masses['hex'] + masses['deoxy'])]
                if 'phosphate' in modifications:
                    masses = masses[(masses['phosphate'] <= masses['hex'] + masses['deoxy'])]
                if 'omethyl' in modifications:
                    masses = masses[(masses['omethyl'] <= masses['hex'] + masses['deoxy'])]
                if 'oacetyl' in modifications:
                    masses = masses[(masses['oacetyl'] <= masses['hex'] + masses['deoxy'])]
                if 'nacetyl' in modifications:
                    masses = masses[(masses['nacetyl'] <= masses['hex'] + masses['deoxy'])]
                if 'amino' in modifications:
                    masses = masses[(masses['amino'] <= masses['hex'] + masses['deoxy'])]
                if unsaturated_option == 'y':
                    masses = masses[(masses['unsaturated'] <= masses['hex'] + masses['deoxy'])]
            if 'deoxy' not in modifications:
                if 'sulfate' in modifications:
                    if double_sulfate == True:
                        masses = masses[(masses['sulfate'] <= (masses['hex']) * 2)]
                    if double_sulfate == False:
                        masses = masses[(masses['sulfate'] <= masses['hex'] )]
                if 'phosphate' in modifications:
                    masses = masses[(masses['phosphate'] <= masses['hex'] )]
                if 'omethyl' in modifications:
                    masses = masses[(masses['omethyl'] <= masses['hex'] )]
                if 'oacetyl' in modifications:
                    masses = masses[(masses['oacetyl'] <= masses['hex'] )]
                if 'nacetyl' in modifications:
                    masses = masses[(masses['nacetyl'] <= masses['hex'] )]
                if 'amino' in modifications:
                    masses = masses[(masses['amino'] <= masses['hex'] )]
                if unsaturated_option == 'y':
                    masses = masses[(masses['unsaturated'] <= masses['hex'] )]
        if pent_option == True:
            if 'carboxylicacid' in modifications:
                masses = masses[(masses['carboxylicacid'] <= masses['hex'] + masses['pent'])]
            if alditol_option == 'y':
                masses = masses[(masses['alditol'] <= masses['hex'] + masses['pent'])]
            if 'deoxy' in modifications:
                if 'sulfate' in modifications:
                    if double_sulfate == True:
                        masses = masses[(masses['sulfate'] <= (masses['hex'] + masses['pent']+ masses['deoxy']) * 2)]
                    if double_sulfate == False:
                        masses = masses[(masses['sulfate'] <= masses['hex'] + masses['pent']+ masses['deoxy'])]
                if 'phosphate' in modifications:
                    masses = masses[(masses['phosphate'] <= masses['hex'] + masses['deoxy'] + masses['pent'])]
                if 'omethyl' in modifications:
                    masses = masses[(masses['omethyl'] <= masses['hex'] + masses['deoxy'] + masses['pent'])]
                if 'oacetyl' in modifications:
                    masses = masses[(masses['oacetyl'] <= masses['hex'] + masses['deoxy'] + masses['pent'])]
                if 'nacetyl' in modifications:
                    masses = masses[(masses['nacetyl'] <= masses['hex'] + masses['deoxy'] + masses['pent'])]
                if 'amino' in modifications:
                    masses = masses[(masses['amino'] <= masses['hex'] + masses['deoxy'] + masses['pent'])]
                if unsaturated_option == 'y':
                    masses = masses[(masses['unsaturated'] <= masses['hex'] + masses['deoxy'] + masses['pent'])]
            if 'deoxy' not in modifications:
                if 'sulfate' in modifications:
                    if double_sulfate == True:
                        masses = masses[(masses['sulfate'] <= (masses['hex'] + masses['pent']) * 2)]
                    if double_sulfate == False:
                        masses = masses[(masses['sulfate'] <= masses['hex'] + masses['pent'])]
                if 'phosphate' in modifications:
                    masses = masses[(masses['phosphate'] <= masses['hex'] + masses['pent'])]
                if 'omethyl' in modifications:
                    masses = masses[(masses['omethyl'] <= masses['hex'] + masses['pent'])]
                if 'oacetyl' in modifications:
                    masses = masses[(masses['oacetyl'] <= masses['hex'] + masses['pent'])]
                if 'nacetyl' in modifications:
                    masses = masses[(masses['nacetyl'] <= masses['hex'] + masses['pent'])]
                if 'amino' in modifications:
                    masses = masses[(masses['amino'] <= masses['hex'] + masses['pent'])]
                if unsaturated_option == 'y':
                    masses = masses[(masses['unsaturated'] <= masses['hex'] + masses['pent'])]
    molecules = list(masses.drop(['dp', "mass"], axis=1).columns)
    if "index" in molecules:
        molecules.remove("index")
        masses = masses.drop(columns=['index'])
    molecule_numbers = np.array(masses[molecules])
    atom_names = ["C", "H", "N", "O", "S", "P"]
    for i in range(len(atom_names)):
        atom = atom_names[i]
        tmp = np.multiply(np.array(molecule_numbers), np.array([list(map(itemgetter(i), [formulas.get(key) for key in molecules]))]).repeat(len(masses), axis=0)).sum(axis=1)
        water_loss = np.multiply(np.array(masses.dp)-1, formulas['water'][i])
        tmp = tmp + water_loss
        if label in proa_names:
            tmp = tmp + [formulas['proca'][i]]
        if label in pa_names:
            tmp = tmp + [formulas['pa'][i]]
        if label in aba_names:
            tmp = tmp + [formulas['aba'][i]]
        if label in ab_names:
            tmp = tmp + [formulas['ab'][i]]
        if label in pmp_names:
                tmp = tmp + [formulas['ab'][i]]
        if label in aba_nonreductive_names:
                tmp = tmp + [formulas['aba_nonreductive'][i]]
        if label in aq_nonreductive_names:
            tmp = tmp + [formulas['aq_nonreductive'][i]]
        if label in aminopentyllinker_names:
            tmp = tmp + [formulas['aminopentyllinker'][i]]
        if custom_label_name != 'none':
            tmp = tmp + [custom_label_formdiff[i]]
        masses[atom] = list(tmp)
        #print("added to formula " + atom)
    masses['formula'] = "C" + pd.Series(masses["C"]).astype(str) + "H" + pd.Series(masses["H"]).astype(str) + "N" + pd.Series(masses["N"]).astype(str) + "O" + pd.Series(masses["O"]).astype(str) + "S" + pd.Series(masses["S"]).astype(str) + "P" + pd.Series(masses["P"]).astype(str)
    masses['formula'] = masses['formula'].astype(str).replace(r"\D0", "", regex=True)
    masses['formula'] = masses['formula'].astype(str).replace("N1O", "NO")
    masses['formula'] = masses['formula'].astype(str).replace("P1", "P")
    masses['formula'] = masses['formula'].astype(str).replace("S1$", "S", regex=True)
    masses = masses.drop(atom_names, axis = 1)
    del tmp, water_loss, molecule_numbers, molecules
    if 'anhydrobridge' in modifications and pent_option == True:
        indexDelete = masses[masses.hex < masses.anhydrobridge].index
        masses.drop(indexDelete, inplace=True)
        masses = masses.reset_index()
    #print("\nstep #5: configuring names to convention")
    #print("----------------------------------------------------------------\n")
    #reorder modifications
    if "sialicacid" in modifications: modifications.insert(0, modifications.pop(modifications.index('sialicacid')))
    if "deoxy" in modifications: modifications.insert(0, modifications.pop(modifications.index('deoxy')))
    if pent_option==False:
        if "none" not in modifications: molecules_names = ['hex'] + modifications
        else: molecules_names = ['hex']
    if pent_option==True:
        if "none" not in modifications: molecules_names = ['hex', 'pent'] + modifications
        else: molecules_names = ['hex', 'pent']
    #get numbers
    if unsaturated_option == 'y': molecules_names = ['unsaturated'] + molecules_names
    if alditol_option == 'y': molecules_names = ['alditol'] + molecules_names
    if dehydrated_option == 'y': molecules_names = ['dehydrated'] + molecules_names
    molecule_numbers = masses[molecules_names]
    if "IUPAC" in naming:
        masses['IUPAC name'] = molecule_numbers[molecules_names].apply(lambda row: ' '.join(
            f'{names_iupac[name]}' + str(int(row[name])) for name in molecules_names if row[name] != 0), axis=1)
        if unsaturated_option == 'y': masses['IUPAC name'] = masses['IUPAC name'].str.replace("Unsaturated1", "Unsaturated", regex=True)
        if alditol_option == 'y': masses['IUPAC name'] = masses['IUPAC name'].str.replace("Alditol1", "Alditol", regex=True)
        if dehydrated_option == 'y': masses['IUPAC name'] = masses['IUPAC name'].str.replace("Dehydrated1", "Dehydrated", regex=True)
        if label in proa_names: masses['IUPAC name'] = masses['IUPAC name'] + ' procA'
        if label in pa_names: masses['IUPAC name'] = masses['IUPAC name'] + ' 2-PA'
        if label in aba_names: masses['IUPAC name'] = masses['IUPAC name'] + ' 2-AA'
        if label in ab_names: masses['IUPAC name'] = masses['IUPAC name'] + ' 2-AB'
        if label in pmp_names: masses['IUPAC name'] = masses['IUPAC name'] + ' bis-PMP'
        if label in aba_nonreductive_names: masses['IUPAC name'] = masses['IUPAC name'] + ' 2-AA'
        if label in aq_nonreductive_names: masses['IUPAC name'] = masses['IUPAC name'] + ' 3-AQ'
        if label in aminopentyllinker_names: masses['IUPAC name'] = masses['IUPAC name'] + ' NH2 Pent1'
        if custom_label_name != 'none': masses['IUPAC name'] = masses['IUPAC name'] + ' ' + custom_label_name
    if "GlycoCT" in naming:
        masses['GlycoCT name'] = molecule_numbers[molecules_names].apply(lambda row: ''.join(
            f'{bracket_mapping[name][0]}{names_glycoct[name]}{bracket_mapping[name][1]}' + str(int(row[name])) for name in molecules_names if row[name] != 0), axis=1)
        if alditol_option == 'y': masses['GlycoCT name'] = masses['GlycoCT name'].str.replace("[ALD]1", "[ALD]", regex=False)
        if dehydrated_option == 'y': masses['GlycoCT name'] = masses['GlycoCT name'].str.replace("[Y]1", "[Y]", regex=False)
        if label in proa_names: masses['GlycoCT name'] = masses['GlycoCT name'] + ' procA'
        if label in pa_names: masses['GlycoCT name'] = masses['GlycoCT name'] + ' 2-PA'
        if label in aba_names: masses['GlycoCT name'] = masses['GlycoCT name'] + ' 2-AA'
        if label in ab_names: masses['GlycoCT name'] = masses['GlycoCT name'] + ' 2-AB'
        if label in pmp_names: masses['GlycoCT name'] = masses['GlycoCT name'] + ' bis-PMP'
        if label in aba_nonreductive_names: masses['GlycoCT name'] = masses['GlycoCT name'] + ' 2-AA'
        if label in aq_nonreductive_names: masses['GlycoCT name'] = masses['GlycoCT name'] + ' 3-AQ'
        if label in aminopentyllinker_names: masses['GlycoCT name'] = masses['GlycoCT name'] + ' NH2 Pent1'
        if custom_label_name != 'none': masses['GlycoCT name'] = masses['GlycoCT name'] + ' ' + custom_label_name
    if "Oxford" in naming:
        masses['Oxford name'] = molecule_numbers[molecules_names].apply(lambda row: ''.join(
            f'{bracket_mapping[name][0]}{names_oxford[name]}{bracket_mapping[name][1]}' + str(int(row[name])) for name in
            molecules_names if row[name] != 0), axis=1)
        if alditol_option == 'y': masses['Oxford name'] = masses['Oxford name'].str.replace("[o]1", "[o]", regex=False)
        if dehydrated_option == 'y': masses['Oxford name'] = masses['Oxford name'].str.replace("[Y]1", "[Y]", regex=False)
        if label in proa_names: masses['Oxford name'] = masses['Oxford name'] + ' procA'
        if label in pa_names: masses['Oxford name'] = masses['Oxford name'] + ' 2-PA'
        if label in aba_names: masses['Oxford name'] = masses['Oxford name'] + ' 2-AA'
        if label in ab_names: masses['Oxford name'] = masses['Oxford name'] + ' 2-AB'
        if label in pmp_names: masses['Oxford name'] = masses['Oxford name'] + ' bis-PMP'
        if label in aba_nonreductive_names: masses['Oxford name'] = masses['Oxford name'] + ' 2-AA'
        if label in aq_nonreductive_names: masses['Oxford name'] = masses['Oxford name'] + ' 3-AQ'
        if label in aminopentyllinker_names: masses['Oxford name'] = masses['Oxford name'] + ' NH2 Pent1'
        if custom_label_name != 'none': masses['Oxford name'] = masses['Oxford name'] + ' ' + custom_label_name
    #print("\nstep #7: calculating m/z values of ions")
    #print("----------------------------------------------------------------\n")
    if len(list(set(modifications).intersection(modifications_anionic))) >= 1:
        if "pos" in polarity:
            #create separate tables of sugars with (any) anionic modifications, and with (only) neutral modifications
            anionic_mod_used = list(set(modifications).intersection(modifications_anionic))
            masses_anionic = masses[masses[anionic_mod_used].gt(0).any(axis=1)]
            masses_neutral = masses[masses[anionic_mod_used].eq(0).all(axis=1)]
            # calculate m/z values for NEUTRAL molecules
            for a in adducts:
                if a == 'H':
                    masses_neutral['[M+H]+'] = masses_neutral.mass + ion_mdiff['H'] - e_mdiff
                if a == 'Na':
                    masses_neutral['[M+Na]+'] = masses_neutral.mass + ion_mdiff['Na'] - e_mdiff
                if a == 'NH4':
                    masses_neutral['[M+NH4]+'] = masses_neutral.mass + ion_mdiff['NH4'] - e_mdiff
                if a == 'K':
                    masses_neutral['[M+K]+'] = masses_neutral.mass + ion_mdiff['K'] - e_mdiff
            # calculate m/z values for ANIONIC molecules
            # calculate number of anionic groups
            if len(anionic_mod_used) > 1:
                masses_anionic['nmod_anionic'] = masses_anionic[anionic_mod_used].sum(axis=1)
                masses_anionic['nmod_anionic'] = masses_anionic.nmod_anionic.astype(int)
            elif len(anionic_mod_used) == 1:
                masses_anionic['nmod_anionic'] = masses_anionic[anionic_mod_used].astype(int)
            for a in adducts:
                if a == "H":
                    masses_anionic['[M+H]+'] = masses_anionic.mass + ion_mdiff['H'] - e_mdiff
                if a == "Na":
                    H_ions = list(range(0, masses_anionic.nmod_anionic.max() + 1))
                    Me_ions = [x + 1 for x in H_ions]
                    ions = list("[M-" + pd.Series(H_ions).astype(str) + "H+" + pd.Series(Me_ions).astype(str) + "Na]+")
                    for i in range(len(ions)):
                        masses_anionic[ions[i]] = masses_anionic.mass - (ion_mdiff['H'] * i) + (ion_mdiff['Na'] * (i + 1)) - e_mdiff
                        masses_anionic[ions[i]] = masses_anionic[ions[i]].where(masses_anionic['nmod_anionic'] >= (i))
                        masses_anionic = masses_anionic.rename({'[M-0H+1Na]+': '[M+Na]+'}, axis=1)
                        masses_anionic = masses_anionic.rename({'[M-1H+2Na]+': '[M-H+2Na]+'}, axis=1)
                if a == "K":
                    H_ions = list(range(0, masses_anionic.nmod_anionic.max() + 1))
                    Me_ions = [x + 1 for x in H_ions]
                    ions = list("[M-" + pd.Series(H_ions).astype(str) + "H+" + pd.Series(Me_ions).astype(str) + "K]+")
                    for i in range(len(ions)):
                        masses_anionic[ions[i]] = masses_anionic.mass - (ion_mdiff['H'] * i) + (
                                    ion_mdiff['K'] * (i + 1)) - e_mdiff
                        masses_anionic[ions[i]] = masses_anionic[ions[i]].where(masses_anionic['nmod_anionic'] >= (i))
                        masses_anionic = masses_anionic.rename({'[M-0H+1K]+': '[M+K]+'}, axis=1)
                        masses_anionic = masses_anionic.rename({'[M-1H+2K]+': '[M-H+2K]+'}, axis=1)
                if a == "NH4" :
                    H_ions = list(range(0, masses_anionic.nmod_anionic.max() + 1))
                    Me_ions = [x + 1 for x in H_ions]
                    ions = list("[M-" + pd.Series(H_ions).astype(str) + "H+" + pd.Series(Me_ions).astype(str) + "NH4]+")
                    for i in range(len(ions)):
                        masses_anionic[ions[i]] = masses_anionic.mass - (ion_mdiff['H'] * i) + (
                                    ion_mdiff['NH4'] * (i + 1)) - e_mdiff
                        masses_anionic[ions[i]] = masses_anionic[ions[i]].where(masses_anionic['nmod_anionic'] >= (i))
                        masses_anionic = masses_anionic.rename({'[M-0H+1NH4]+': '[M+NH4]+'}, axis=1)
                        masses_anionic = masses_anionic.rename({'[M-1H+2NH4]+': '[M-H+2NH4]+'}, axis=1)
            masses = pd.concat([masses_anionic, masses_neutral])
            del masses_anionic, masses_neutral
        if "neg" in polarity:
            #create separate tables of sugars with (any) anionic modifications, and with (only) neutral modifications
            anionic_mod_used = list(set(modifications).intersection(modifications_anionic))
            masses_anionic = masses[masses[anionic_mod_used].gt(0).any(axis=1)]
            masses_neutral = masses[masses[anionic_mod_used].eq(0).all(axis=1)]
            # calculate m/z values for NEUTRAL molecules
            for a in adducts:
                if a == 'H':
                    masses_neutral['[M-H]-'] = masses_neutral.mass - ion_mdiff['H'] + e_mdiff
                if a == 'Cl':
                    masses_neutral['[M+Cl]-'] = masses_neutral.mass + ion_mdiff['Cl'] + e_mdiff
                if a == 'CHOO':
                    masses_neutral['[M+CHOO]-'] = masses_neutral.mass + ion_mdiff['CHOO'] + e_mdiff
            # calculate m/z values for ANIONIC molecules
            # calculate number of anionic groups
            if len(anionic_mod_used) > 1:
                masses_anionic['nmod_anionic'] = masses_anionic[anionic_mod_used].sum(axis=1)
                masses_anionic['nmod_anionic'] = masses_anionic.nmod_anionic.astype(int)
            elif len(anionic_mod_used) == 1:
                masses_anionic['nmod_anionic'] = masses_anionic[anionic_mod_used].astype(int)
            if "MALDI" in ion_type:
                for a in adducts:
                    if a == "H":
                        masses_anionic['[M-H]-'] = masses_anionic.mass - ion_mdiff['H'] + e_mdiff
                    if a == "Na":
                        H_ions = list(range(2, masses_anionic.nmod_anionic.max() + 1))
                        Me_ions = [x - 1 for x in H_ions]
                        ions = list("[M-" + pd.Series(H_ions).astype(str) + "H+" + pd.Series(Me_ions).astype(str) + "Na]-")
                        for i in range(len(ions)):
                            masses_anionic[ions[i]] = masses_anionic.mass - (ion_mdiff['H'] * (i + 2)) + (ion_mdiff['Na'] * (i + 1)) + e_mdiff
                            masses_anionic[ions[i]] = masses_anionic[ions[i]].where(masses_anionic['nmod_anionic'] >= (i + 2))
                            masses_anionic = masses_anionic.rename({'[M-2H+1Na]-': '[M-2H+Na]-'}, axis=1)
                    if a == "K":
                        H_ions = list(range(2, masses_anionic.nmod_anionic.max() + 1))
                        Me_ions = [x - 1 for x in H_ions]
                        ions = list("[M-" + pd.Series(H_ions).astype(str) + "H+" + pd.Series(Me_ions).astype(str) + "K]-")
                        for i in range(len(ions)):
                            masses_anionic[ions[i]] = masses_anionic.mass - (ion_mdiff['H'] * (i+2)) + (ion_mdiff['K'] * (i+1))  + e_mdiff
                            masses_anionic[ions[i]] = masses_anionic[ions[i]].where(masses_anionic['nmod_anionic'] >= (i + 2))
                            masses_anionic = masses_anionic.rename({'[M-2H+1K]-': '[M-2H+K]-'}, axis=1)
                    if a == "NH4" :
                        H_ions = list(range(2, masses_anionic.nmod_anionic.max() + 1))
                        Me_ions = [x - 1 for x in H_ions]
                        ions = list("[M-" + pd.Series(H_ions).astype(str) + "H+" + pd.Series(Me_ions).astype(str) + "NH4]-")
                        for i in range(len(ions)):
                            masses_anionic[ions[i]] = masses_anionic.mass - (ion_mdiff['H'] * (i+2)) + (ion_mdiff['NH4'] * (i+1))  + e_mdiff
                            masses_anionic[ions[i]] = masses_anionic[ions[i]].where(masses_anionic['nmod_anionic'] >= (i + 2))
                            masses_anionic = masses_anionic.rename({'[M-2H+1NH4]-': '[M-2H+NH4]-'}, axis=1)
                    if a == 'Cl':
                        masses_anionic['[M+Cl]-'] = masses_anionic.mass + ion_mdiff['Cl'] + e_mdiff
                    if a == 'CHOO':
                        masses_anionic['[M+CHOO]-'] = masses_anionic.mass + ion_mdiff['CHOO'] + e_mdiff
            if "ESI" in ion_type:
                if "H" in adducts:
                    ions = list(range(1, masses_anionic.nmod_anionic.max() + 1))
                    ions = list("[M-" + pd.Series(ions).astype(str) + "H]-" + pd.Series(ions).astype(str))
                    for i in range(len(ions)):
                        masses_anionic[ions[i]] = (masses_anionic.mass - ion_mdiff['H'] * (i + 1) + e_mdiff * (i + 1)) / (i + 1)
                        masses_anionic[ions[i]] = masses_anionic[ions[i]].where(masses_anionic['nmod_anionic'] >= (i + 1))
                        masses_anionic = masses_anionic.rename({'[M-1H]-1': '[M-H]-'}, axis=1)
                for a in adducts:
                    if a == 'Cl':
                        masses_anionic['[M+Cl]-'] = masses_anionic.mass + ion_mdiff['Cl'] + e_mdiff
                    if a == 'CHOO':
                        masses_anionic['[M+CHOO]-'] = masses_anionic.mass + ion_mdiff['CHOO'] + e_mdiff
            masses = pd.concat([masses_anionic, masses_neutral])
            del masses_anionic, masses_neutral
        my_cols = list(masses.filter(like='[M', axis=1).columns)
        masses[my_cols] = masses[my_cols].where(masses[my_cols] >= scan_range[0])
        masses[my_cols] = masses[my_cols].where(masses[my_cols] <= scan_range[1])
        masses = masses.dropna(subset=my_cols, how='all')
        bad_cols = {'level_0','index','nmod','nmod_avg','nmod_anionic','_merge', 'k', 'x', 'name'}
        #bad_cols.update(modifications_anionic)
        #bad_cols.update(modifications_neutral)
        cols_del = list(set(masses.columns).intersection(bad_cols))
        masses = masses.drop(columns=cols_del)
    if len(list(set(modifications).intersection(modifications_anionic))) == 0:
        # calculate m/z values for neutral molecules
        if "neg" in polarity:
            for a in adducts:
                if a == 'H':
                    masses['[M-H]-'] = masses.mass - ion_mdiff['H'] + e_mdiff
                if a == 'CHOO':
                    masses['[M+CHOO]-'] = masses.mass + ion_mdiff['CHOO'] + e_mdiff
                if a == 'Cl':
                    masses['[M+Cl]-'] = masses.mass + ion_mdiff['Cl'] + e_mdiff
        if "pos" in polarity:
            for a in adducts:
                if a == 'H':
                    masses['[M+H]+'] = masses.mass + ion_mdiff['H'] - e_mdiff
                if a == 'Na':
                    masses['[M+Na]+'] = masses.mass + ion_mdiff['Na'] - e_mdiff
                if a == 'NH4':
                    masses['[M+NH4]+'] = masses.mass + ion_mdiff['NH4'] - e_mdiff
                if a == 'K':
                    masses['[M+K]+'] = masses.mass + ion_mdiff['K'] - e_mdiff
        # filter molecules based on scan range
        # set values outside range to NaN
        # remove rows where all ions are outside range
        my_cols = list(masses.filter(like='[M', axis=1).columns)
        masses[my_cols] = masses[my_cols].where(masses[my_cols] >= scan_range[0])
        masses[my_cols] = masses[my_cols].where(masses[my_cols] <= scan_range[1])
        masses = masses.dropna(subset=my_cols, how='all')
        # format nicely to only have useful columns
        bad_cols = {'level_0', 'index', 'nmod', 'nmod_avg', 'nmod_anionic', '_merge', 'k', 'x', 'name'}
        # bad_cols.update(modifications_anionic)
        # bad_cols.update(modifications_neutral)
        cols_del = list(set(masses.columns).intersection(bad_cols))
        masses = masses.drop(columns=cols_del)
    #print("\nstep #8: returning output")
    #print("----------------------------------------------------------------\n")
    masses = masses.reset_index(drop=True)
    masses = masses.sort_values(by = ["dp", "mass"])
    return(masses)
