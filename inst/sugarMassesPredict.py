#!/usr/bin/env python
# import modules
import pandas as pd
import numpy as np
import itertools
import math
import re
from operator import itemgetter

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
                          'sulphate']
# hexose and water masses to build molecule base
hex_mass = 180.06339
water_mass = 18.010565
# mass differences for modifications
pent_mdiff = -30.010566
modifications_mdiff = {
    "sulphate": 79.956817,
    "anhydrobridge": -water_mass,
    "omethyl": 14.01565,
    "carboxylicacid": 13.979265,
    "sialicacid" : 129.042594,
    "nacetyl": 41.026549,
    "oacetyl": 42.010565,
    "phosphate": 79.966333,
    "deoxy": -15.994915,
    "unsaturated": -water_mass,
    "alditol": 2.015650,
    "amino": -0.984016,
    "dehydrated": -water_mass
}

# mass differences for labels
proca_mdiff = 219.173546
pa_mdiff = 78.058183
aba_mdiff = 121.052764
pmp_mdiff = 330.148061
ab_mdiff = 120.068748

#possible names for labels
proa_names = {"procainamide", "proca", "procA", "ProA"}
pa_names = {"2-ap", "2-AP", "pa", "PA", "2-aminopyridine"}
aba_names = {"2-aa", "2-AA", "aba", "ABA", "2-aminobenzoic acid"}
ab_names = {"2-ab", "2-AB", "ab", "AB", "2-aminobenzamide"}
pmp_names = {"pmp", "PMP", "1-phenyl-3-methyl-5-pyrazolone"}


# mass differences for ions
ion_mdiff = {
    "H": 1.00782500000003,
    "Na": 22.98977,
    "Cl": 34.968853,
    "CHOO": 44.997655,
    "NH4": 18.034374,
    "K": 38.963708,
    "Ca": 39.962591
}
e_mdiff = 0.000548579909

# formulas

#chnosp
formulas = {
    "hex": [6, 12, 0, 6, 0, 0],
    "pent": [5, 10, 0, 5, 0, 0],
    "water": [0, -2, 0, -1, 0, 0],
    "sulphate": [0, 0, 0, 3, 1, 0],
    "anhydrobridge": [0, -2, 0, -1, 0, 0],
    "omethyl": [1, 2, 0, 0, 0, 0],
    "carboxylicacid": [0, -2, 0, 1, 0, 0],
    "sialicacid": [5, 7, 1, 3, 0, 0],
    "nacetyl": [2, 3, 1, 0, 0, 0],
    "oacetyl": [2, 2, 0, 1, 0, 0],
    "phosphate": [0, 1, 0, 3, 0, 1],
    "deoxy": [0, 0, 0, -1, 0, 0],
    "proca": [13, 21, 3, 0, 0, 0],
    "unsaturated":  [0, -2, 0, -1, 0, 0],
    "alditol": [0, +2, 0, 0, 0, 0],
    "amino": [0, +1, +1, -1, 0, 0],
    "dehydrated": [0, -2, 0, -1, 0, 0],
    "pa": [5, 6, 2, -1, 0, 0] ,
    "aba": [7, 7, 1, 1, 0, 0] ,
    "pmp": [20, 18, 4, 1, 0, 0],
    "ab": [7, 8, 2, 0, 0, 0]
}
# modification types
modifications_anionic = {"sulphate",
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
    "sulphate": 'Sulfate',
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
    "dehydrated": 'Dehydrated'
}
names_glycoct = {
    "hex": 'HEX',
    'pent': 'PEN',
    "sulphate": 'SO4',
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
    "dehydrated": 'Y'
}
names_oxford = {
    "hex": 'H',
    'pent': 'P',
    "sulphate": 'S',
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
    "dehydrated": 'Y'
}

#brackets
bracket_mapping = {
    "hex": '()',
    'pent': '()',
    "sulphate": '[]',
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
    "dehydrated": '[]'
}

def predict_sugars(dp= [1, 6], polarity='neg', scan_range=[175, 1400], pent_option=False, modifications='none', nmod_max=1, double_sulphate=False, label='none', ion_type = "ESI", format="long", adducts = ["all"], naming = "IUPAC"):
    if 'all' in adducts:
        adducts=['H', 'Cl', 'CHOO', 'Na', 'NH4', 'K']
    if type(adducts)==str:
        adducts = [adducts]
    dp_range_list = list(range(dp[0], dp[1] + 1))
    #print("step #1: getting arguments")
    #print("----------------------------------------")
    if "all" in modifications:
        modifications = possible_modifications
    if "sulphate" in modifications and len(modifications) > 1:
        modifications.append(modifications.pop(modifications.index('sulphate')))
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
    #print("\nstep #2: calculating all possible masses")
    #print("----------------------------------------\n")
    # build hexose molecules
    #print("--> getting hexose masses")
    def getHexMasses(dp_range_list):
        dp = pd.Series(dp_range_list)
        name = "hex-" + dp.astype(str)
        hex = dp
        mass = dp * hex_mass - (dp - 1) * water_mass
        masses = pd.DataFrame({'dp': dp,
                               'name': name,
                               'hex': hex.astype(int),
                               'mass': mass})
        return masses
    masses = getHexMasses(dp_range_list)
    #define additional functions
    def dpRepeats(dp_range_list):
        repeats_list = []
        for i in dp_range_list:
            repeats_list = repeats_list + list(range(0, i + 1))
        return repeats_list
    def getPentMasses(masses):
        dp = masses.dp.repeat(masses.dp.array + 1).reset_index(drop=True)
        pent = pd.Series(dpRepeats(dp_range_list))
        hex = dp - pent
        name = "hex-" + hex.astype(str) + "-pent-" + pent.astype(str)
        name = name.str.replace("hex-0-", "")
        name = name.str.replace("-pent-0", "")
        mass = masses.mass.repeat(masses.dp.array + 1).reset_index(drop=True)
        mass = mass + pent * pent_mdiff
        masses = pd.DataFrame({'dp': dp,
                               'name': name,
                               'hex': hex,
                               'pent': pent,
                               'mass': mass})
        return masses
    def getModificationNumbers(dp_range_list, m, pent_option, modifications):
        modification_numbers = []
        for i in dp_range_list:
            a = list(range(0, i + 1))
            if pent_option == True:
                modification_numbers = modification_numbers + \
                                       list(itertools.product(a, repeat=m)) * (i + 1)
            elif pent_option == False:
                modification_numbers = modification_numbers + \
                                       list(itertools.product(a, repeat=m))
        modification_numbers = pd.DataFrame(modification_numbers)
        modification_numbers.columns = modifications
        return modification_numbers
    if pent_option == 1:
        #print("--> getting pentose masses")
        masses = getPentMasses(masses)
    #add modifications
    if "none" not in modifications and pent_option == True:
        #print("--> adding modifications")
        m = len(modifications)
        dp = masses.dp.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        hex = masses.hex.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        pent = masses.pent.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        modification_numbers = getModificationNumbers(dp_range_list, m, pent_option, modifications)
        name = "hex-" + hex.astype(str) + "-pent-" + pent.astype(str)
        for i in range(m):
            name = name + "-" + modifications[i] + "-" + modification_numbers[modifications[i]].astype(str)
        name = name.str.replace("-\D+-0", "", regex=True)
        name = name.str.replace("hex-0-", "")
        mass = masses.mass.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        for i in range(m):
            mass = mass + modifications_mdiff[modifications[i]] * modification_numbers[modifications[i]]
        masses = pd.DataFrame({'dp': dp,
                               'name': name,
                               'hex': hex,
                               'pent': pent})
        masses = pd.concat([masses, modification_numbers], axis=1)
        masses['mass'] = mass
        del dp, hex, pent, modification_numbers, name, mass
    if "none" not in modifications and pent_option == False:
        #print("--> adding modifications")
        m = len(modifications)
        dp = masses.dp.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        hex = masses.hex.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        modification_numbers = getModificationNumbers(dp_range_list, m, pent_option, modifications)
        name = "hex-" + hex.astype(str)
        for i in range(m):
            name = name + "-" + modifications[i] + "-" + modification_numbers[modifications[i]].astype(str)
        name = name.str.replace("-\D+-0", "", regex=True)
        name = name.str.replace("hex-0-", "")
        mass = masses.mass.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        for i in range(m):
            mass = mass + modifications_mdiff[modifications[i]] * modification_numbers[modifications[i]]
        masses = pd.DataFrame({'dp': dp,
                               'name': name,
                               'hex': hex})
        masses = pd.concat([masses, modification_numbers], axis=1)
        masses['mass'] = mass
        del dp, hex, modification_numbers, name, mass
    if "sulphate" in modifications and double_sulphate == True:
        #print("--> adding extra sulphate groups")
        masses_s1 = masses.loc[masses['sulphate'] >= 1]
        masses_s2 = masses_s1
        masses_s2.sulphate = masses_s1.sulphate + masses_s1.dp
        masses_s2.name = masses_s2.name.str.replace("-sulphate-\d{1,2}", "", regex=True)
        masses_s2.name = masses_s2.name + '-sulphate-' + masses_s2.sulphate.astype(str)
        masses_s2.mass = masses_s2.mass + modifications_mdiff['sulphate'] * masses_s2.dp
        #masses = masses.append(masses_s2).reset_index()
        masses = pd.concat([masses, masses_s2]).reset_index()
        del masses_s1
        del masses_s2
    if label in proa_names:
        #print("--> adding procainamide label")
        masses['name'] = masses.name + ' procA'
        masses['mass'] = masses.mass + proca_mdiff
    if label in pa_names:
        #print("--> adding 2-aminopyridine")
        masses['name'] = masses.name + ' 2-PA'
        masses['mass'] = masses.mass + pa_mdiff
    if label in aba_names:
        #print("--> adding 2-aminobenzoic acid label")
        masses['name'] = masses.name + ' 2-AA'
        masses['mass'] = masses.mass + aba_mdiff
    if label in ab_names:
        #print("--> adding 2-aminobenzamide")
        masses['name'] = masses.name + ' 2-AB'
        masses['mass'] = masses.mass + ab_mdiff
    if label in pmp_names:
        #print("--> adding bis-PMP label")
        masses['name'] = masses.name + ' bis-PMP'
        masses['mass'] = masses.mass + pmp_mdiff
    if unsaturated_option == 'y':
        #print("--> adding unsaturated sugars")
        masses_a = masses.copy()
        masses_a.name = "unsaturated-" + masses.name
        masses_a['unsaturated'] = 1
        masses['unsaturated'] = 0
        masses_a.mass = masses.mass + modifications_mdiff['unsaturated']
        masses = pd.concat([masses, masses_a]).reset_index(drop=True)
        del masses_a
    if alditol_option == 'y':
        #print("--> adding alditol sugars")
        masses_a = masses.copy()
        masses_a.name = "alditol-" + masses_a.name
        masses_a['alditol'] = 1
        masses['alditol'] = 0
        masses_a.mass = masses_a.mass + modifications_mdiff['alditol']
        masses = pd.concat([masses, masses_a]).reset_index(drop=True)
        del masses_a
    if dehydrated_option == 'y':
        #print("--> adding dehydration to sugars")
        masses_a = masses.copy()
        masses_a.name = "dehydrated-" + masses_a.name
        masses_a['dehydrated'] = 1
        masses['dehydrated'] = 0
        masses_a.mass = masses_a.mass + modifications_mdiff['dehydrated']
        masses = pd.concat([masses, masses_a]).reset_index(drop=True)
        del masses_a
    #print("\nstep #3: building formulas")
    #print("----------------------------------------\n")
    molecules = list(masses.drop(['dp', "name", "mass"], axis=1).columns)
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
            masses[atom] = list(tmp)
    masses['formula'] = "C" + pd.Series(masses["C"]).astype(str) + "H" + pd.Series(masses["H"]).astype(str) + "N" + pd.Series(masses["N"]).astype(str) + "O" + pd.Series(masses["O"]).astype(str) + "S" + pd.Series(masses["S"]).astype(str) + "P" + pd.Series(masses["P"]).astype(str)
    masses['formula'] = masses['formula'].str.replace("\D0", "", regex=True)
    masses['formula'] = masses['formula'].str.replace("N1O", "NO")
    masses['formula'] = masses['formula'].str.replace("P1", "P")
    masses['formula'] = masses['formula'].str.replace("S1$", "S", regex=True)
    masses = masses.drop(atom_names, axis = 1)
    del tmp, water_loss, molecule_numbers, molecules
    #print("\nstep #4: filtering based on number of modifications per monomer")
    #print("----------------------------------------------------------------\n")
    if "none" not in modifications:
        masses['nmod'] = masses[modifications].sum(axis=1)
        masses['nmod_avg'] = masses.nmod / masses.dp
        masses = masses.drop(masses[masses.nmod_avg > nmod_max].index)
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
    #remove rows with deoxypentose
    if "deoxy" in modifications and pent_option==True:
        masses = masses[masses['deoxy'] <= masses['hex']]
    # remove rows with sialicacid pentose
    if "sialicacid" in modifications and pent_option == True:
        masses = masses[masses['sialicacid'] <= masses['hex']]
    #get numbers
    if unsaturated_option == 'y': molecules_names = ['unsaturated'] + molecules_names
    if alditol_option == 'y': molecules_names = ['alditol'] + molecules_names
    if dehydrated_option == 'y': molecules_names = ['dehydrated'] + molecules_names
    molecule_numbers = masses[molecules_names]
    #subtract deoxy from hex
    if "deoxy" in modifications:
        molecule_numbers["hex"] = molecule_numbers.hex - molecule_numbers.deoxy
    #subtract sialicacid from hex
    if "sialicacid" in modifications:
        molecule_numbers["hex"] = molecule_numbers.hex - molecule_numbers.sialicacid
    if "IUPAC" in naming:
        masses['IUPAC name'] = molecule_numbers[molecules_names].apply(lambda row: ' '.join(
            f'{names_iupac[name]}' + str(row[name]) for name in molecules_names if row[name] != 0), axis=1)
        if unsaturated_option == 'y': masses['IUPAC name'] = masses['IUPAC name'].str.replace("Unsaturated1", "Unsaturated", regex=True)
        if alditol_option == 'y': masses['IUPAC name'] = masses['IUPAC name'].str.replace("Alditol1", "Alditol", regex=True)
        if dehydrated_option == 'y': masses['IUPAC name'] = masses['IUPAC name'].str.replace("Dehydrated1", "Dehydrated", regex=True)
    if "GlycoCT" in naming:
        masses['GlycoCT name'] = molecule_numbers[molecules_names].apply(lambda row: ''.join(
            f'{bracket_mapping[name][0]}{names_glycoct[name]}{bracket_mapping[name][1]}' + str(row[name]) for name in molecules_names if row[name] != 0), axis=1)
        if unsaturated_option == 'y': masses['GlycoCT name'] = masses['GlycoCT name'].str.replace("UNS1", "UNS", regex=True)
        if alditol_option == 'y': masses['GlycoCT name'] = masses['GlycoCT name'].str.replace("ALD1", "ALD", regex=True)
        if dehydrated_option == 'y': masses['GlycoCT name'] = masses['GlycoCT name'].str.replace("Y1", "Y", regex=True)
    if "Oxford" in naming:
        masses['Oxford name'] = molecule_numbers[molecules_names].apply(lambda row: ''.join(
            f'{bracket_mapping[name][0]}{names_oxford[name]}{bracket_mapping[name][1]}' + str(row[name]) for name in
            molecules_names if row[name] != 0), axis=1)
        if unsaturated_option == 'y': masses['Oxford name'] = masses['Oxford name'].str.replace("U", "U", regex=True)
        if alditol_option == 'y': masses['Oxford name'] = masses['Oxford name'].str.replace("o1", "o", regex=True)
        if dehydrated_option == 'y': masses['Oxford name'] = masses['Oxford name'].str.replace("Y1", "Y", regex=True)
    #print("\nstep #6: calculating m/z values of ions")
    #print("----------------------------------------------------------------\n")
    if len(list(set(modifications).intersection(modifications_anionic))) >= 1:
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
        if "neg" in polarity:
            #create separate tables of sugars with (any) anionic modifications, and with (only) neutral modifications
            anionic_mod_used = list(set(modifications).intersection(modifications_anionic))
            masses_anionic = masses[masses['name'].str.contains('|'.join(anionic_mod_used))]
            masses_all = masses.merge(masses_anionic.drop_duplicates(), how='left', indicator=True)
            masses_neutral = masses_all[masses_all._merge == 'left_only']
            del masses_all
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
        bad_cols = {'level_0','index','hex','pent','alditol','nmod','nmod_avg','nmod_anionic','_merge', 'dehydrated', 'k', 'x', 'name'}
        bad_cols.update(modifications_anionic)
        bad_cols.update(modifications_neutral)
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
        bad_cols = {'level_0','index','alditol','hex','pent','nmod','nmod_avg','nmod_anionic','_merge', 'dehydrated','name'}
        bad_cols.update(modifications_neutral)
        cols_del = list(set(masses.columns).intersection(bad_cols))
        masses = masses.drop(columns=cols_del)
    #print("\nstep #7: returning ouput")
    #print("----------------------------------------------------------------\n")
    masses = masses.reset_index(drop=True)
    masses = masses.sort_values(by = ["dp", "mass"])
    return(masses)