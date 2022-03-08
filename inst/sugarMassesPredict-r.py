#!/usr/bin/env python
# import modules
import pandas as pd
import numpy as np
import itertools

# suppress warnings
pd.options.mode.chained_assignment = None  # default='warn'

#possible modifications
possible_modifications = ['carboxyl',
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
    "carboxyl": 13.979265,
    "nacetyl": 41.026549,
    "oacetyl": 42.010565,
    "phosphate": 79.966333,
    "deoxy": -15.994915,
    "unsaturated": -2.015650,
    "alditol": 2.015650,
    "amino": -0.984016,
    "dehydrated": -water_mass
}

# mass differences for labels
procainamide_mdiff = 219.173546
benzoic_acid_mdiff = 104.026215

# mass differences for ions
ion_mdiff = {
    "H": 1.00782500000003,
    "Na": 22.98977,
    "Cl": 34.968853,
    "CHOO": 44.997655,
    "NH4": 18.034374,
    "K": 38.963708
}

e_mdiff = 0.000548579909

# formulas

formulas = {
    "hex": [6, 12, 0, 6, 0, 0],
    "pent": [5, 10, 0, 5, 0, 0],
    "water": [0, -2, 0, -1, 0, 0],
    "sulphate": [0, 0, 0, 3, 1, 0],
    "anhydrobridge": [0, -2, 0, -1, 0, 0],
    "omethyl": [1, 2, 0, 0, 0, 0],
    "carboxyl": [0, -2, 0, 1, 0, 0],
    "nacetyl": [2, 3, 1, 0, 0, 0],
    "oacetyl": [2, 2, 0, 1, 0, 0],
    "phosphate": [0, 1, 0, 3, 0, 1],
    "deoxy": [0, 0, 0, -1, 0, 0],
    "procainamide": [13, 21, 3, 0, 0, 0],
    "benzoic_acid": [7, 4, 0, 1, 0, 0],
    "unsaturated": [0, -2, 0, 0, 0, 0],
    "alditol": [0, +2, 0, 0, 0, 0],
    "amino": [0, +1, +1, -1, 0, 0],
    "dehydrated": [0, -2, 0, -1, 0, 0]
}
# modification types
modifications_anionic = {"sulphate",
                         "phosphate",
                         "carboxyl"}
modifications_neutral = {"anhydrobridge",
                         "omethyl",
                         "nacetyl",
                         "oacetyl",
                         "deoxy",
                         "unsaturated",
                         "amino",
                         "dehydrated"}

def predict_sugars(dp1=1, dp2=2, ESI_mode='pos', scan_range1=100,
                   scan_range2=1000, pent_option=0, modifications='none', nmod_max=1, double_sulphate=0, label='none'):
    print("\nstep #1: getting and checking arguments")
    print("----------------------------------------\n")
    #get args and check that they are in the correct format
    dp_range = [dp1, dp2]
    if len(dp_range) != 2:
        raise ValueError('dp range does not have two values!')
    if dp_range[1] < dp_range[0]:
        raise ValueError('dp range needs to start with lower value!')
    dp_range_list = list(range(dp_range[0], dp_range[1] + 1))
    if not pent_option in [0,1]:
        raise ValueError('pent option needs to be 0 or 1!')
    if modifications != 'none':
        if not all(elem in possible_modifications for elem in modifications):
            raise ValueError('you have given a modification that is not allowed!')
    if type(ESI_mode) is str:
        if not ESI_mode in ['neg', 'pos']:
            raise ValueError('you have given an ESI mode that is not allowed! did you make a typo?')
    if type(ESI_mode) is list:
        if not all(elem in ['neg', 'pos'] for elem in ESI_mode):
            raise ValueError('you have given an ESI mode that is not allowed! did you make a typo?')
    scan_range = [scan_range1, scan_range2]
    if len(scan_range) != 2:
        raise ValueError('scan range does not have two values!')
    if type(scan_range[0]) is not int or type(scan_range[1]) is not int:
        raise ValueError('scan range includes non-integer values!')
    if scan_range[1] < scan_range[0]:
        raise ValueError('scan range needs to start with lower value!')
    if not label in ['none', 'procainamide']:
        raise ValueError('you have given a label that is not allowed!')
    if "all" in modifications:
        modifications = possible_modifications
    if "sulphate" in modifications:
        modifications.append(modifications.pop(modifications.index('sulphate')))
    if "alditol" in modifications:
        alditol_option = 'y'
        modifications.remove('alditol')
    elif "alditol" not in modifications:
        alditol_option = 'n'
    if "unsaturated" in modifications:
        unsaturated_option = 'y'
        modifications.remove('unsaturated')
    elif "unsaturated" not in modifications:
        unsaturated_option = 'n'
    if "dehydrated" in modifications:
        dehydrated_option = 'y'
        modifications.remove('dehydrated')
    elif "dehydrated" not in modifications:
        dehydrated_option = 'n'
    #calculate possible masses
    print("\nstep #2: calculating all possible masses")
    print("----------------------------------------\n")
    # build hexose molecules
    print("--> getting hexose masses")
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
    # calculate masses for pentose molecules if selected
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
            if pent_option == 1:
                modification_numbers = modification_numbers + \
                                       list(itertools.product(a, repeat=m)) * (i + 1)
            elif pent_option == 0:
                modification_numbers = modification_numbers + \
                                       list(itertools.product(a, repeat=m))
        modification_numbers = pd.DataFrame(modification_numbers)
        modification_numbers.columns = modifications
        return modification_numbers
    if pent_option == 1:
        print("--> getting pentose masses")
        masses = getPentMasses(masses)
    #add modifications
    if "none" in modifications and pent_option == 1:
        masses.name = masses.name.str.replace("hex-0-", "")
        masses.name = masses.name.str.replace("-pent-0", "")
    if "none" not in modifications and pent_option == 1:
        print("--> adding modifications")
        m = len(modifications)
        dp = masses.dp.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        hex = masses.hex.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        pent = masses.pent.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        modification_numbers = getModificationNumbers(dp_range_list, m, pent_option, modifications)
        name = "hex-" + hex.astype(str) + "-pent-" + pent.astype(str)
        for i in range(m):
            name = name + "-" + modifications[i] + "-" + modification_numbers[modifications[i]].astype(str)
        name = name.str.replace("-\D+-0", "")
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
    if "none" not in modifications and pent_option == 0:
        print("--> adding modifications")
        m = len(modifications)
        dp = masses.dp.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        hex = masses.hex.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        modification_numbers = getModificationNumbers(dp_range_list, m, pent_option, modifications)
        name = "hex-" + hex.astype(str)
        for i in range(m):
            name = name + "-" + modifications[i] + "-" + modification_numbers[modifications[i]].astype(str)
        name = name.str.replace("-\D+-0", "")
        name = name.str.replace("hex-0-", "")
        mass = masses.mass.repeat((masses.dp.array + 1) ** m).reset_index(drop=True)
        for i in range(m):
            mass = mass + modifications_mdiff[modifications[i]] * modification_numbers[modifications[i]]
        masses = pd.DataFrame({'dp': dp,
                               'name': name,
                               'hex': hex})
        masses = pd.concat([masses, modification_numbers], axis=1)
        masses['mass'] = mass
    if "sulphate" in modifications and double_sulphate == 1:
        print("--> adding extra sulphate groups")
        masses_s1 = masses.loc[masses['sulphate'] >= 1]
        masses_s2 = masses_s1
        masses_s2.sulphate = masses_s1.sulphate + masses_s1.dp
        masses_s2.name = masses_s2.name.str.replace("-sulphate-\d{1,2}", "")
        masses_s2.name = masses_s2.name + '-sulphate-' + masses_s2.sulphate.astype(str)
        masses_s2.mass = masses_s2.mass + modifications_mdiff['sulphate'] * masses_s2.dp
        masses = masses.append(masses_s2).reset_index()
        del masses_s1
        del masses_s2
    if "procainamide" in label:
        print("--> adding procainamide label")
        masses['name'] = masses.name + '-procA'
        masses['mass'] = masses.mass + procainamide_mdiff
    if unsaturated_option == 'y':
        print("--> adding unsaturated sugars")
        masses_a = masses.copy()
        masses_a.name = "unsaturated-" + masses.name
        masses_a['unsaturated'] = 1
        masses['unsaturated'] = 0
        masses_a.mass = masses.mass + modifications_mdiff['unsaturated']
        masses = masses.append(masses_a).reset_index()
        del masses_a
    if alditol_option == 'y':
        print("--> adding alditol sugars")
        masses_a = masses.copy()
        masses_a.name = "alditol-" + masses_a.name
        masses_a['alditol'] = 1
        masses['alditol'] = 0
        masses_a.mass = masses_a.mass + modifications_mdiff['alditol']
        masses = masses.append(masses_a).reset_index(drop=True)
        del masses_a
    if dehydrated_option == 'y':
        print("--> adding dehydration to sugars")
        masses_a = masses.copy()
        masses_a.name = "dehydrated-" + masses_a.name
        masses_a['dehydrated'] = 1
        masses['dehydrated'] = 0
        masses_a.mass = masses_a.mass + modifications_mdiff['dehydrated']
        masses = masses.append(masses_a).reset_index(drop=True)
        del masses_a
    print("\nstep #3: building formulas")
    print("----------------------------------------\n")
    if "none" in modifications and pent_option == 1:
        dp = masses.dp
        hex = masses.hex
        pent = masses.pent
        molecule_numbers = pd.DataFrame({'dp': dp,'hex': hex,'pent': pent})
        molecules = list(molecule_numbers.drop('dp', axis=1).columns)
        atom_names = ["C", "H", "N", "O", "S", "P"]
        atom_list = []
        for i in range(len(atom_names)):
            n = np.array([0] * len(masses.index))
            for j in range(len(molecules)):
                form_n = np.array([formulas[molecules[j]][i]] * len(masses.index))
                mol_n = np.array(molecule_numbers[molecules[j]])
                form_mol_n = form_n * mol_n
                n = n + form_mol_n
            if "procainamide" in label:
                p = np.array([formulas['procainamide'][i]] * len(masses.index))
                n = n + p
            atom_list.append(list(n))
        # remove molecules from formula for glycosidic bonds
        atom_list_2 = []
        for i in range(len(atom_names)):
            n = np.array(atom_list[i])
            form_n = np.array([formulas['water'][i]] * len(masses.index))
            mol_n = np.array(molecule_numbers['dp'] - 1)
            form_mol_n = form_n * mol_n
            n = n + form_mol_n
            atom_list_2.append(list(n))
        # concatenate to build formulas
        for i in range(len(atom_names)):
            if i == 0:
                formulas_final = atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
            else:
                formulas_final = formulas_final.astype(str) + atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
        # fix to remove atoms with zero
        formulas_final = formulas_final.str.replace("\D0", "")
        masses['formula'] = formulas_final
    if "none" in modifications and pent_option == 0:
        dp = masses.dp
        hex = masses.hex
        molecule_numbers = pd.DataFrame({'dp': dp,'hex': hex})
        molecules = list(molecule_numbers.drop('dp', axis=1).columns)
        atom_names = ["C", "H", "N", "O", "S", "P"]
        atom_list = []
        for i in range(len(atom_names)):
            n = np.array([0] * len(masses.index))
            for j in range(len(molecules)):
                form_n = np.array([formulas[molecules[j]][i]] * len(masses.index))
                mol_n = np.array(molecule_numbers[molecules[j]])
                form_mol_n = form_n * mol_n
                n = n + form_mol_n
            if "procainamide" in label:
                p = np.array([formulas['procainamide'][i]] * len(masses.index))
                n = n + p
            atom_list.append(list(n))
        # remove molecules from formula for glycosidic bonds
        atom_list_2 = []
        for i in range(len(atom_names)):
            n = np.array(atom_list[i])
            form_n = np.array([formulas['water'][i]] * len(masses.index))
            mol_n = np.array(molecule_numbers['dp'] - 1)
            form_mol_n = form_n * mol_n
            n = n + form_mol_n
            atom_list_2.append(list(n))
        # concatenate to build formulas
        for i in range(len(atom_names)):
            if i == 0:
                formulas_final = atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
            else:
                formulas_final = formulas_final.astype(str) + atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
        # fix to remove atoms with zero
        formulas_final = formulas_final.str.replace("\D0", "")
        masses['formula'] = formulas_final
    if "none" not in modifications and pent_option == 1:
        if unsaturated_option == 'y':
            modifications.append('unsaturated')
        if alditol_option == 'y':
            modifications.append('alditol')
        if dehydrated_option == 'y':
            modifications.append('dehydrated')
        dp = masses.dp
        hex = masses.hex
        pent = masses.pent
        molecule_numbers = pd.DataFrame({'dp': dp,'hex': hex,'pent': pent})
        modification_numbers = masses[modifications]
        molecule_numbers = pd.concat([molecule_numbers, modification_numbers], axis=1)
        molecules = list(molecule_numbers.drop('dp', axis=1).columns)
        atom_names = ["C", "H", "N", "O", "S", "P"]
        atom_list = []
        for i in range(len(atom_names)):
            n = np.array([0] * len(masses.index))
            for j in range(len(molecules)):
                form_n = np.array([formulas[molecules[j]][i]] * len(masses.index))
                mol_n = np.array(molecule_numbers[molecules[j]])
                form_mol_n = form_n * mol_n
                n = n + form_mol_n
            if "procainamide" in label:
                p = np.array([formulas['procainamide'][i]] * len(masses.index))
                n = n + p
            atom_list.append(list(n))
        # remove molecules from formula for glycosidic bonds
        atom_list_2 = []
        for i in range(len(atom_names)):
            n = np.array(atom_list[i])
            form_n = np.array([formulas['water'][i]] * len(masses.index))
            mol_n = np.array(molecule_numbers['dp'] - 1)
            form_mol_n = form_n * mol_n
            n = n + form_mol_n
            atom_list_2.append(list(n))
        # concatenate to build formulas
        for i in range(len(atom_names)):
            if i == 0:
                formulas_final = atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
            else:
                formulas_final = formulas_final.astype(str) + atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
        # fix to remove atoms with zero
        formulas_final = formulas_final.str.replace("\D0", "")
        masses['formula'] = formulas_final
    if "none" not in modifications and pent_option == 0:
        if unsaturated_option == 'y':
            modifications.append('unsaturated')
        if alditol_option == 'y':
            modifications.append('alditol')
        if dehydrated_option == 'y':
            modifications.append('dehydrated')
        dp = masses.dp
        hex = masses.hex
        molecule_numbers = pd.DataFrame({'dp': dp,'hex': hex})
        modification_numbers = masses[modifications]
        molecule_numbers = pd.concat([molecule_numbers, modification_numbers], axis=1)
        molecules = list(molecule_numbers.drop('dp', axis=1).columns)
        atom_names = ["C", "H", "N", "O", "S", "P"]
        atom_list = []
        for i in range(len(atom_names)):
            n = np.array([0] * len(masses.index))
            for j in range(len(molecules)):
                form_n = np.array([formulas[molecules[j]][i]] * len(masses.index))
                mol_n = np.array(molecule_numbers[molecules[j]])
                form_mol_n = form_n * mol_n
                n = n + form_mol_n
            if "procainamide" in label:
                p = np.array([formulas['procainamide'][i]] * len(masses.index))
                n = n + p
            atom_list.append(list(n))
        # remove molecules from formula for glycosidic bonds
        atom_list_2 = []
        for i in range(len(atom_names)):
            n = np.array(atom_list[i])
            form_n = np.array([formulas['water'][i]] * len(masses.index))
            mol_n = np.array(molecule_numbers['dp'] - 1)
            form_mol_n = form_n * mol_n
            n = n + form_mol_n
            atom_list_2.append(list(n))
        # concatenate to build formulas
        for i in range(len(atom_names)):
            if i == 0:
                formulas_final = atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
            else:
                formulas_final = formulas_final.astype(str) + atom_names[i] + pd.Series(atom_list_2[i]).astype(str)
        # fix to remove atoms with zero
        formulas_final = formulas_final.str.replace("\D0", "")
        masses['formula'] = formulas_final
    print("\nstep #4: filtering based on number of modifications per monomer")
    print("----------------------------------------------------------------\n")
    if "none" not in modifications:
        if unsaturated_option == 'y':
            modifications.remove('unsaturated')
        if alditol_option == 'y':
            modifications.remove('alditol')
        if dehydrated_option == 'y':
            modifications.remove('dehydrated')
        masses['nmod'] = masses[modifications].sum(axis=1)
        masses['nmod_avg'] = masses.nmod / masses.dp
        masses = masses.drop(masses[masses.nmod_avg > nmod_max].index)
    if 'anhydrobridge' in modifications and pent_option == 1:
        indexDelete = masses[masses.hex < masses.anhydrobridge].index
        masses.drop(indexDelete, inplace=True)
        masses = masses.reset_index()
    print("\nstep #5: calculating m/z values of ions")
    print("----------------------------------------------------------------\n")
    if len(list(set(modifications).intersection(modifications_anionic))) >= 1:
        # create separate tables of sugars with (any) anionic modifications, and with (only) neutral modifications
        anionic_mod_used = list(set(modifications).intersection(modifications_anionic))
        masses_anionic = masses[masses['name'].str.contains('|'.join(anionic_mod_used))]
        masses_all = masses.merge(masses_anionic.drop_duplicates(), how='left', indicator=True)
        masses_neutral = masses_all[masses_all._merge == 'left_only']
        # calculate m/z values for neutral molecules
        if "neg" in ESI_mode:
            masses_neutral['[M-H]-'] = masses_neutral.mass - ion_mdiff['H'] + e_mdiff
            masses_neutral['[M+Cl]-'] = masses_neutral.mass + ion_mdiff['Cl'] + e_mdiff
            masses_neutral['[M+CHOO]-'] = masses_neutral.mass + ion_mdiff['CHOO'] + e_mdiff
            masses_neutral['[M-2H]-2'] = (masses_neutral.mass - 2 * ion_mdiff['H'] + 2 * e_mdiff) / 2
            masses_neutral['[M+2Cl]-2'] = (masses_neutral.mass + 2 * ion_mdiff['Cl'] + 2 * e_mdiff) / 2
            masses_neutral['[M+2CHOO]-2'] = (masses_neutral.mass + 2 * ion_mdiff['CHOO'] + 2 * e_mdiff) / 2
            masses_neutral['[M+Cl-H]-2'] = (masses_neutral.mass + ion_mdiff['Cl'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
            masses_neutral['[M+CHOO-H]-2'] = (masses_neutral.mass + ion_mdiff['CHOO'] - ion_mdiff[
                'H'] + 2 * e_mdiff) / 2
            masses_neutral['[M+CHOO+Cl]-2'] = (masses_neutral.mass + ion_mdiff['CHOO'] + ion_mdiff[
                'Cl'] + 2 * e_mdiff) / 2
        if "pos" in ESI_mode:
            masses_neutral['[M+H]+'] = masses_neutral.mass + ion_mdiff['H'] - e_mdiff
            masses_neutral['[M+Na]+'] = masses_neutral.mass + ion_mdiff['Na'] - e_mdiff
            masses_neutral['[M+NH4]+'] = masses_neutral.mass + ion_mdiff['NH4'] - e_mdiff
            masses_neutral['[M+K]+'] = masses_neutral.mass + ion_mdiff['K'] - e_mdiff
        # filter neutral molecules based on scan range
        # set values outside range to NaN
        # remove rows where all ions are outside range
        my_cols = list(masses_neutral.filter(like='[M', axis=1).columns)
        masses_neutral[my_cols] = masses_neutral[my_cols].where(masses_neutral[my_cols] >= scan_range[0])
        masses_neutral[my_cols] = masses_neutral[my_cols].where(masses_neutral[my_cols] <= scan_range[1])
        masses_neutral = masses_neutral.dropna(subset=my_cols, how='all')
        # calculate m/z values for anionic molecules
        if len(anionic_mod_used) > 1:
            masses_anionic['nmod_anionic'] = masses_anionic[anionic_mod_used].sum(axis=1)
            masses_anionic['nmod_anionic'] = masses_anionic.nmod_anionic.astype(int)
        elif len(anionic_mod_used) == 1:
            masses_anionic['nmod_anionic'] = masses_anionic[anionic_mod_used].astype(int)
        if "neg" in ESI_mode:
            ions = list(range(1, masses_anionic.nmod_anionic.max() + 1))
            ions = list("[M-" + pd.Series(ions).astype(str) + "H]-" + pd.Series(ions).astype(str))
            for i in range(len(ions)):
                masses_anionic[ions[i]] = (masses_anionic.mass - ion_mdiff['H'] * (i + 1) + e_mdiff * (i + 1)) / (i + 1)
                masses_anionic[ions[i]] = masses_anionic[ions[i]].where(masses_anionic['nmod_anionic'] >= (i + 1))
            masses_anionic = masses_anionic.rename({'[M-1H]-1': '[M-H]-'}, axis=1)
            masses_anionic['[M+Cl]-'] = masses_anionic.mass + ion_mdiff['Cl'] + e_mdiff
            masses_anionic['[M+CHOO]-'] = masses_anionic.mass + ion_mdiff['CHOO'] + e_mdiff
            masses_anionic['[M+2Cl]-2'] = (masses_anionic.mass + 2 * ion_mdiff['Cl'] + 2 * e_mdiff) / 2
            masses_anionic['[M+2CHOO]-2'] = (masses_anionic.mass + 2 * ion_mdiff['CHOO'] + 2 * e_mdiff) / 2
            masses_anionic['[M+Cl-H]-2'] = (masses_anionic.mass + ion_mdiff['Cl'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
            masses_anionic['[M+CHOO-H]-2'] = (masses_anionic.mass + ion_mdiff['CHOO'] - ion_mdiff[
                'H'] + 2 * e_mdiff) / 2
            masses_anionic['[M+CHOO+Cl]-2'] = (masses_anionic.mass + ion_mdiff['CHOO'] + ion_mdiff[
                'Cl'] + 2 * e_mdiff) / 2
        if "pos" in ESI_mode:
            masses_anionic['[M+H]+'] = masses_anionic.mass + ion_mdiff['H'] - e_mdiff
            masses_anionic['[M+Na]+'] = masses_anionic.mass + ion_mdiff['Na'] - e_mdiff
            masses_anionic['[M+NH4]+'] = masses_anionic.mass + ion_mdiff['NH4'] - e_mdiff
            masses_anionic['[M+K]+'] = masses_anionic.mass + ion_mdiff['K'] - e_mdiff
        # filter anionic molecules based on scan range
        # set values outside range to NaN
        # remove rows where all ions are outside range
        my_cols = list(masses_anionic.filter(like='[M', axis=1).columns)
        masses_anionic[my_cols] = masses_anionic[my_cols].where(masses_anionic[my_cols] >= scan_range[0])
        masses_anionic[my_cols] = masses_anionic[my_cols].where(masses_anionic[my_cols] <= scan_range[1])
        masses_anionic = masses_anionic.dropna(subset=my_cols, how='all')
        # concatenate dataframes and format nicely to only have useful columns
        masses_final = pd.concat([masses_anionic, masses_neutral])
        bad_cols = {'level_0','index','hex','pent','alditol','nmod','nmod_avg','nmod_anionic','_merge', 'dehydrated'}
        bad_cols.update(modifications_anionic)
        bad_cols.update(modifications_neutral)
        cols_del = list(set(masses_final.columns).intersection(bad_cols))
        masses_final = masses_final.drop(columns=cols_del)
    if len(list(set(modifications).intersection(modifications_anionic))) == 0:
        # calculate m/z values for neutral molecules
        if "neg" in ESI_mode:
            masses['[M-H]-'] = masses.mass - ion_mdiff['H'] + e_mdiff
            masses['[M+Cl]-'] = masses.mass + ion_mdiff['Cl'] + e_mdiff
            masses['[M+CHOO]-'] = masses.mass + ion_mdiff['CHOO'] + e_mdiff
            masses['[M-2H]-2'] = (masses.mass - 2 * ion_mdiff['H'] + 2 * e_mdiff) / 2
            masses['[M+2Cl]-2'] = (masses.mass + 2 * ion_mdiff['Cl'] + 2 * e_mdiff) / 2
            masses['[M+2CHOO]-2'] = (masses.mass + 2 * ion_mdiff['CHOO'] + 2 * e_mdiff) / 2
            masses['[M+Cl-H]-2'] = (masses.mass + ion_mdiff['Cl'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
            masses['[M+CHOO-H]-2'] = (masses.mass + ion_mdiff['CHOO'] - ion_mdiff['H'] + 2 * e_mdiff) / 2
            masses['[M+CHOO+Cl]-2'] = (masses.mass + ion_mdiff['CHOO'] + ion_mdiff['Cl'] + 2 * e_mdiff) / 2
        if "pos" in ESI_mode:
            masses['[M+H]+'] = masses.mass + ion_mdiff['H']
            masses['[M+Na]+'] = masses.mass + ion_mdiff['Na']
            masses['[M+NH4]+'] = masses.mass + ion_mdiff['NH4'] - e_mdiff
            masses['[M+K]+'] = masses.mass + ion_mdiff['K'] - e_mdiff
        # filter neutral molecules based on scan range
        # set values outside range to NaN
        # remove rows where all ions are outside range
        my_cols = list(masses.filter(like='[M', axis=1).columns)
        masses[my_cols] = masses[my_cols].where(masses[my_cols] >= scan_range[0])
        masses[my_cols] = masses[my_cols].where(masses[my_cols] <= scan_range[1])
        masses = masses.dropna(subset=my_cols, how='all')
        # format nicely to only have useful columns
        masses_final = masses
        bad_cols = {'level_0','index','alditol','hex','pent','nmod','nmod_avg','nmod_anionic','_merge', 'dehydrated'}
        bad_cols.update(modifications_neutral)
        cols_del = list(set(masses_final.columns).intersection(bad_cols))
        masses_final = masses_final.drop(columns=cols_del)
    print("\nstep #6: returning ouput")
    print("----------------------------------------------------------------\n")
    return(masses_final)