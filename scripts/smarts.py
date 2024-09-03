from rdkit import Chem

structure=input('structure type: ')


## defining groups for ONO

if structure == "ONO" :

    # para groups
    para_Br = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#35]):[#6]:[#6]:[#6]:1-[#8]')
    para_N = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#7](-[#1])-[#1]):[#6]:[#6]:[#6]:1-[#8]')
    para_OH = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#8]-[#1]):[#6]:[#6]:[#6]:1-[#8]')
    para_Cl = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#17]):[#6]:[#6]:[#6]:1-[#8]')
    para_F = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#9]):[#6]:[#6]:[#6]:1-[#8]')
    para_Me = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#6](-[#1])(-[#1])-[#1]):[#6]:[#6]:[#6]:1-[#8]')
    para_OMe = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6]([#8]-[#6](-[#1])(-[#1])-[#1]):[#6]:[#6]:[#6]:1-[#8]')
    para_NMe2 = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#7](-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:[#6]:[#6]:1-[#8]')
    para_CN = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#6]#[#7]):[#6]:[#6]:[#6]:1-[#8]')
    para_tBu = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#6](-[#6](-[H])(-[H])-[H])(-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:[#6]:[#6]:1-[#8]')
    para_CF3 = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#6](-[#9])(-[#9])-[#9]):[#6]:[#6]:[#6]:1-[#8]')
    para_NO2 = Chem.MolFromSmarts('[#7]-[#6]1:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:[#6]:1-[#8]')

    #meta_groups

    meta_Br = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#35]):[#6]:[#6]:[#6]:1-[#7]')
    meta_N = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#7](-[#1])-[#1]):[#6]:[#6]:[#6]:1-[#7]')
    meta_OH = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#8]-[#1]):[#6]:[#6]:[#6]:1-[#7]')
    meta_Cl = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#17]):[#6]:[#6]:[#6]:1-[#7]')
    meta_F = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#9]):[#6]:[#6]:[#6]:1-[#7]')
    meta_Me = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#6](-[#1])(-[#1])-[#1]):[#6]:[#6]:[#6]:1-[#7]')
    meta_OMe = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6]([#8]-[#6](-[#1])(-[#1])-[#1]):[#6]:[#6]:[#6]:1-[#7]')
    meta_NMe2 = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#7](-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:[#6]:[#6]:1-[#7]')
    meta_CN = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#6]#[#7]):[#6]:[#6]:[#6]:1-[#7]')
    meta_tBu = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#6](-[#6](-[H])(-[H])-[H])(-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:[#6]:[#6]:1-[#7]')
    meta_CF3 = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#6](-[#9])(-[#9])-[#9]):[#6]:[#6]:[#6]:1-[#7]')
    meta_NO2 = Chem.MolFromSmarts('[#8]-[#6]1:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:[#6]:1-[#7]')

    #ortho groups

    ortho_Br = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#35]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_N = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#7](-[#1])-[#1]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_OH = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#8]-[#1]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_Cl = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#17]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_F = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#9]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_Me = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#6](-[#1])(-[#1])-[#1]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_OMe = Chem.MolFromSmarts('[#8]-[#6]1:[#6]([#8]-[#6](-[#1])(-[#1])-[#1]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_NMe2 = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#7](-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_CN = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#6]#[#7]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_tBu = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#6](-[#6](-[H])(-[H])-[H])(-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_CF3 = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#6](-[#9])(-[#9])-[#9]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    ortho_NO2 = Chem.MolFromSmarts('[#8]-[#6]1:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:[#6]:[#6]:1-[#7]')
    
elif structure == 'triarylboranes':
    ## defining groups for triarylboranes

    # para groups
    para_Br = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#35]):[#6]:[#6]:1')
    para_N = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#7](-[#1])-[#1]):[#6]:[#6]:1')
    para_OH = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#8]-[#1]):[#6]:[#6]:1')
    para_Cl = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#17]):[#6]:[#6]:1')
    para_F = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#9]):[#6]:[#6]:1')
    para_Me = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#6](-[#1])(-[#1])-[#1]):[#6]:[#6]:1')
    para_OMe = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]([#8]-[#6](-[#1])(-[#1])-[#1]):[#6]:[#6]:1')
    para_NMe2 = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#7](-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:[#6]:1')
    para_CN = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#6]#[#7]):[#6]:[#6]:1')
    para_tBu = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#6](-[#6](-[H])(-[H])-[H])(-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:[#6]:1')
    para_CF3 = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#6](-[#9])(-[#9])-[#9]):[#6]:[#6]:1')
    para_NO2 = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:[#6]:1')

    #meta_groups

    meta_Br = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#35]):[#6]:1')
    meta_N = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#7](-[#1])-[#1]):[#6]:1')
    meta_OH = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#8]-[#1]):[#6]:1')
    meta_Cl = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#17]):[#6]:1')
    meta_F = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#9]):[#6]:1')
    meta_Me = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#6](-[#1])(-[#1])-[#1]):[#6]:1')
    meta_OMe = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]([#8]-[#6](-[#1])(-[#1])-[#1]):[#6]:1')
    meta_NMe2 = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#7](-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:1')
    meta_CN = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#6]#[#7]):[#6]:1')
    meta_tBu = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#6](-[#6](-[H])(-[H])-[H])(-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]):[#6]:1')
    meta_CF3 = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#6](-[#9])(-[#9])-[#9]):[#6]:1')
    meta_NO2 = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6](-[#7+](-[#8-])=[#8]):[#6]:1')

    #ortho groups

    ortho_Br = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#35]')
    ortho_N = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7](-[#1])-[#1]')
    ortho_OH = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]-[#1]')
    ortho_Cl = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#17]')
    ortho_F = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#9]')
    ortho_Me = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](-[#1])(-[#1])-[#1]')
    ortho_OMe = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]-[#6](-[#1])(-[#1])-[#1]')
    ortho_NMe2 = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7](-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]')
    ortho_CN = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]#[#7]')
    ortho_tBu = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](-[#6](-[H])(-[H])-[H])(-[#6](-[H])(-[H])-[H])-[#6](-[H])(-[H])-[H]')
    ortho_CF3 = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6](-[#9])(-[#9])-[#9]')
    ortho_NO2 = Chem.MolFromSmarts('[#5]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7+](-[#8-])=[#8]')
    
else :
    pass

def id_substructure_para(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    if mol.HasSubstructMatch(para_Br):
        return('Br')
    elif mol.HasSubstructMatch(para_N):
        return('NH2')
    elif mol.HasSubstructMatch(para_OH):
        return('OH')
    elif mol.HasSubstructMatch(para_Cl):
        return('Cl')
    elif mol.HasSubstructMatch(para_F):
        return('F')
    elif mol.HasSubstructMatch(para_Me):
        return('Me')
    elif mol.HasSubstructMatch(para_OMe):
        return('OCH3')
    elif mol.HasSubstructMatch(para_NMe2):
        return('NMe2')
    elif mol.HasSubstructMatch(para_CN):
        return('CN')
    elif mol.HasSubstructMatch(para_tBu):
        return('tBu')
    elif mol.HasSubstructMatch(para_CF3):
        return('CF3')
    elif mol.HasSubstructMatch(para_NO2):
        return('NO2')
    else :
        return('H')

def id_substructure_meta(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    if mol.HasSubstructMatch(meta_Br):
        return('Br')
    elif mol.HasSubstructMatch(meta_N):
        return('NH2')
    elif mol.HasSubstructMatch(meta_OH):
        return('OH')
    elif mol.HasSubstructMatch(meta_Cl):
        return('Cl')
    elif mol.HasSubstructMatch(meta_F):
        return('F')
    elif mol.HasSubstructMatch(meta_Me):
        return('Me')
    elif mol.HasSubstructMatch(meta_OMe):
        return('OCH3')
    elif mol.HasSubstructMatch(meta_NMe2):
        return('NMe2')
    elif mol.HasSubstructMatch(meta_CN):
        return('CN')
    elif mol.HasSubstructMatch(meta_tBu):
        return('tBu')
    elif mol.HasSubstructMatch(meta_CF3):
        return('CF3')
    elif mol.HasSubstructMatch(meta_NO2):
        return('NO2')
    else :
        return('H')

def id_substructure_ortho(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    if mol.HasSubstructMatch(ortho_Br):
        return('Br')
    elif mol.HasSubstructMatch(ortho_N):
        return('NH2')
    elif mol.HasSubstructMatch(ortho_OH):
        return('OH')
    elif mol.HasSubstructMatch(ortho_Cl):
        return('Cl')
    elif mol.HasSubstructMatch(ortho_F):
        return('F')
    elif mol.HasSubstructMatch(ortho_Me):
        return('Me')
    elif mol.HasSubstructMatch(ortho_OMe):
        return('OCH3')
    elif mol.HasSubstructMatch(ortho_NMe2):
        return('NMe2')
    elif mol.HasSubstructMatch(ortho_CN):
        return('CN')
    elif mol.HasSubstructMatch(ortho_tBu):
        return('tBu')
    elif mol.HasSubstructMatch(ortho_CF3):
        return('CF3')
    elif mol.HasSubstructMatch(ortho_NO2):
        return('NO2')
    else :
        return('H')  