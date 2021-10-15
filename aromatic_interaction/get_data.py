import pynmrstar
from mmcif.io.PdbxReader import PdbxReader
import os
import numpy
import gzip

aromatic_atoms = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
    'W7F': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
    #'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2', 'xx', 'yy']  # if needed un comment
}


def get_coordinates(cif_file, use_auth_tag=False,nmrbox=False):
    """
    Extract coordinate information from cif file as a dictionary
    {model_id : {(seq_id,chain_id,res_id,atom_id) : array[x,y,x],...},...}
    :param cif_file: Input coordinate file
    :return: dictionary
    """
    cif_data = []
    if nmrbox:
        ifh = gzip.open(cif_file,'rt')
    else:
        ifh = open(cif_file, 'r')
    pRd = PdbxReader(ifh)
    pRd.read(cif_data)
    ifh.close()
    c0 = cif_data[0]
    seq_info = c0.getObj('entity_poly')
    c_names=seq_info.getAttributeList()
    ch_id = c_names.index('pdbx_strand_id')
    sq_id=c_names.index('pdbx_seq_one_letter_code')
    sq_info={}
    for d in seq_info.getRowList():
        sq_info[d[ch_id]]=(len(d[sq_id]),d[sq_id].count('W'))
    print (sq_info)
    atom_site = c0.getObj('atom_site')
    max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
    col_names = atom_site.getAttributeList()
    model_id = col_names.index('pdbx_PDB_model_num')
    x_id = col_names.index('Cartn_x')
    y_id = col_names.index('Cartn_y')
    z_id = col_names.index('Cartn_z')
    atom_id = col_names.index('label_atom_id')
    comp_id = col_names.index('label_comp_id')
    asym_id = col_names.index('label_asym_id')
    entity_id = col_names.index('label_entity_id')
    seq_id = col_names.index('label_seq_id')
    icode_id = col_names.index('pdbx_PDB_ins_code')
    alt_id = col_names.index('label_alt_id')
    aut_seq_id = col_names.index('auth_seq_id')
    aut_asym_id = col_names.index('auth_asym_id')
    aut_atom_id = col_names.index('auth_atom_id')
    aut_comp_id = col_names.index('auth_comp_id')
    pdb_models = {}
    atom_ids = {}
    for model in range(1, max_models + 1):
        pdb = {}
        aid = []
        for dat in atom_site.getRowList():
            if dat[comp_id] in aromatic_atoms and dat[atom_id] in aromatic_atoms[dat[comp_id]]:  # Only necessary coordinates for this
                # calculation
                if int(dat[model_id]) == model:
                    if use_auth_tag:
                        residue=(dat[aut_seq_id], dat[aut_asym_id], dat[aut_comp_id])
                        if residue not in aid:
                            aid.append(residue)
                        pdb[(dat[aut_seq_id], dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                            numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
                    else:
                       residue = (dat[seq_id], dat[asym_id], dat[comp_id])
                       if residue not in aid:
                           aid.append(residue)
                       pdb[(dat[seq_id], dat[asym_id], dat[comp_id], dat[atom_id])] = numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
        pdb_models[model] = pdb
        atom_ids[model] = aid
    return pdb_models, atom_ids, sq_info


def get_pdb_data(pdb_id, auth_tag=False, nmrbox=False, local_file=False):
    if local_file:
        pdb_data, ar_residues, sq_info = get_coordinates(pdb_id, auth_tag)
    else:
        if nmrbox:
            try:
                cif_file_path = '/reboxitory/2021/07/PDB/data/structures/all/mmCIF/{}.cif.gz'.format(pdb_id.lower())
                pdb_data, ar_residues,sq_info = get_coordinates(cif_file_path,auth_tag,nmrbox=nmrbox)
            except FileNotFoundError:
                pdb_data, ar_residues = None, None
        else:
            if not os.path.isdir('./data'):
                os.system('mkdir ./data')
            if not os.path.isdir('./data/PDB'):
                os.system('mkdir ./data/PDB')
            cif_file = './data/PDB/{}.cif'.format(pdb_id)
            if not os.path.isfile(cif_file):
                cmd = 'wget https://files.rcsb.org/download/{}.cif -O ./data/PDB/{}.cif'.format(pdb_id, pdb_id)
                os.system(cmd)
            pdb_data,ar_residues,sq_info = get_coordinates('./data/PDB/{}.cif'.format(pdb_id),auth_tag)
    return pdb_data, ar_residues, sq_info


if __name__ == "__main__":
    x,y,z=get_pdb_data('3K0N')
    print (x)
    print (y)
    print(z)