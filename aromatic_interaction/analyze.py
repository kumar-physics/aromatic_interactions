from mmcif.io.PdbxReader import PdbxReader
import os
import numpy
import gzip
import sys

aromatic_atoms = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
}


def get_coordinates(cif_file, use_auth_tag=False, nmrbox=False):
    """
    Extract coordinate information from cif file as a dictionary
    {model_id : {(seq_id,chain_id,res_id,atom_id) : array[x,y,x],...},...}
    :param cif_file: Input coordinate file
    :return: dictionary
    """
    cif_data = []
    if nmrbox:
        ifh = gzip.open(cif_file, 'rt')
    else:
        ifh = open(cif_file, 'r')
    pRd = PdbxReader(ifh)
    pRd.read(cif_data)
    ifh.close()
    c0 = cif_data[0]
    seq_info = c0.getObj('entity_poly')
    c_names = seq_info.getAttributeList()
    ch_id = c_names.index('pdbx_strand_id')
    sq_id = c_names.index('pdbx_seq_one_letter_code')
    sq_info = {}
    for d in seq_info.getRowList():
        sq_info[d[ch_id]] = (len(d[sq_id]), d[sq_id].count('W'))
    print(sq_info)
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
            if dat[comp_id] in aromatic_atoms and dat[atom_id] in aromatic_atoms[
                dat[comp_id]]:  # Only necessary coordinates for this
                # calculation
                if int(dat[model_id]) == model:
                    if use_auth_tag:
                        residue = (dat[aut_seq_id], dat[aut_asym_id], dat[aut_comp_id])
                        if residue not in aid:
                            aid.append(residue)
                        pdb[(dat[aut_seq_id], dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                            numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
                    else:
                        residue = (dat[seq_id], dat[asym_id], dat[comp_id])
                        if residue not in aid:
                            aid.append(residue)
                        pdb[(dat[seq_id], dat[asym_id], dat[comp_id], dat[atom_id])] = numpy.array(
                            [float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
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
                pdb_data, ar_residues, sq_info = get_coordinates(cif_file_path, auth_tag, nmrbox=nmrbox)
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
            pdb_data, ar_residues, sq_info = get_coordinates('./data/PDB/{}.cif'.format(pdb_id), auth_tag)
    return pdb_data, ar_residues, sq_info


ring_atoms = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'W7F': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', ]  # if needed un comment
}


def get_distance(c1, c2):
    """
    Calculates the distance between two coordinate points
    :param c1: array of x,y,z
    :param c2: array of x,y,z
    :return: distance between two ponts
    """
    return round(numpy.linalg.norm(c1 - c2), 4)


def get_centroid(p):
    # print (len(p),p)
    x = [i[0] for i in p]
    y = [i[1] for i in p]
    z = [i[2] for i in p]
    c = [sum(x) / len(p), sum(y) / len(p), sum(z) / len(p)]
    return numpy.array(c)


def find_angle(p, c, d):
    pc, p = p
    v1 = p[1] - pc
    v2 = p[2] - pc
    nv = numpy.cross(v1, v2)
    # nv2 = numpy.cross(v2, v1)
    cv = c - pc
    # cv2 = c - cn
    nnv = nv / numpy.linalg.norm(nv)
    ncv = cv / numpy.linalg.norm(cv)
    # ncv2 = cv2 / numpy.linalg.norm(cv2)
    dp = abs(numpy.dot(nnv, ncv))
    # dp2 = abs(numpy.dot(nnv, ncv2))
    ang = numpy.arccos(dp)
    # ang2 = numpy.arccos(dp2)
    ang_deg = (180 / numpy.pi) * ang
    # ang_deg2 = (180 / numpy.pi) * ang2
    # print(ang_deg)
    s_ang = solid_angle(ang_deg, d)
    return round(ang_deg, 4), s_ang


def find_angle2(p, c, d):
    pc, p = p
    v1 = p[1] - pc
    v2 = p[2] - pc
    # nv = numpy.cross(v1, v2)
    # nv2 = numpy.cross(v2, v1)
    cv = c - pc
    # cv2 = c - cn
    nv = p[3] - pc
    nnv = nv / numpy.linalg.norm(nv)
    ncv = cv / numpy.linalg.norm(cv)
    # ncv2 = cv2 / numpy.linalg.norm(cv2)
    dp = numpy.dot(nnv, ncv)
    # dp2 = abs(numpy.dot(nnv, ncv2))
    ang = numpy.arccos(dp)
    # ang2 = numpy.arccos(dp2)
    ang_deg = (180 / numpy.pi) * ang
    # ang_deg2 = (180 / numpy.pi) * ang2
    # print(ang_deg)
    s_ang = solid_angle(ang_deg, d)
    return round(ang_deg, 4), s_ang


def find_angle3(p, c, d):
    cz2, cent = p
    v1 = cz2 - cent
    v2 = c - cent
    nv1 = v1 / numpy.linalg.norm(v1)
    nv2 = v2 / numpy.linalg.norm(v2)
    # ncv2 = cv2 / numpy.linalg.norm(cv2)
    dp = numpy.dot(nv1, nv2)
    # dp2 = abs(numpy.dot(nnv, ncv2))
    ang = numpy.arccos(dp)
    # ang2 = numpy.arccos(dp2)
    ang_deg = (180 / numpy.pi) * ang
    # ang_deg2 = (180 / numpy.pi) * ang2
    # print(ang_deg)
    # s_ang = solid_angle(ang_deg, d)
    return round(ang_deg, 4)


def find_angle_between_plane(p1, p2):
    pc1, p1 = p1
    pc2, p2 = p2
    v11 = p1[1] - pc1
    v12 = p1[2] - pc1
    nv1 = numpy.cross(v11, v12)
    v21 = p2[1] - pc2
    v22 = p2[2] - pc2
    nv2 = numpy.cross(v21, v22)
    # nv2 = numpy.cross(v2, v1)
    # cv = c - pc
    # cv2 = c - cn
    nnv = nv1 / numpy.linalg.norm(nv1)
    ncv = nv2 / numpy.linalg.norm(nv2)
    # ncv2 = cv2 / numpy.linalg.norm(cv2)
    dp = abs(numpy.dot(nnv, ncv))
    # dp2 = abs(numpy.dot(nnv, ncv2))
    ang = numpy.arccos(dp)
    # ang2 = numpy.arccos(dp2)
    ang_deg = (180 / numpy.pi) * ang
    # ang_deg2 = (180 / numpy.pi) * ang2
    # print(ang_deg)
    return round(ang_deg, 4)


def solid_angle(a_deg, r):
    s = 1.4  # C-C bond length
    A = ((3.0 * numpy.sqrt(3)) / 2.0) * s * s  # area of the hexagonal plane
    a = (numpy.pi / 180) * a_deg
    # sa2=2*numpy.pi*(1.0-1.0/(numpy.sqrt(1+(A*numpy.cos(a)/(numpy.pi*r1**r1)))))
    sa = 2 * numpy.pi * (1.0 - 1.0 / (numpy.sqrt(1 + (A * numpy.cos(a) / (numpy.pi * r * r)))))
    # print (a_deg)
    sa_deg = (180.0 / numpy.pi) * sa
    # sa_deg2 = (180.0 / numpy.pi) * sa2
    # print (a_deg,sa_deg,sa_deg2)
    return round(sa_deg, 4)


def get_aromatic_info(pdb_data):
    aromtic_residues = sorted(
        list(set([(int(i[0]), i[1], i[2]) for i in pdb_data[1].keys() if i[2] in aromatic_atoms.keys()])))
    ar_info = {}
    for m in pdb_data.keys():
        ar_info[m] = {}
        for ar_res in aromtic_residues:
            p = []
            for atm in ring_atoms[ar_res[2]]:
                try:
                    p.append(pdb_data[m][(str(ar_res[0]), ar_res[1], ar_res[2], atm)])
                except KeyError:
                    pass
            if len(p) > 1:
                c = get_centroid(p)
                ar_info[m][(str(ar_res[0]), ar_res[1], ar_res[2])] = (c, p)
    return ar_info


def calculate_interaction(pdb_id, out_dir, nmrbox=False):
    fo = open('{}/{}.csv'.format(out_dir, pdb_id), 'w')
    pdb_data, aromatic_residues, sq_ifo = get_pdb_data(pdb_id, nmrbox=nmrbox, )
    for m in pdb_data:
        for r1 in range(len(aromatic_residues[m])):
            aro_res1 = aromatic_residues[m][r1]
            if aro_res1[2] == 'TRP':
                for r2 in range(len(aromatic_residues[m])):
                    aro_res2 = aromatic_residues[m][r2]
                    if aro_res2 != aro_res1:
                        try:
                            ring1 = [pdb_data[m][(aro_res1[0], aro_res1[1], aro_res1[2], i)] for i in
                                     ring_atoms[aro_res1[2]]]
                            ring2 = [pdb_data[m][(aro_res2[0], aro_res2[1], aro_res2[2], i)] for i in
                                     ring_atoms[aro_res2[2]]]
                            c1 = get_centroid(ring1)
                            c2 = get_centroid(ring2)
                            d = get_distance(c1, c2)
                            azi_ange1, s_ang1 = find_angle((c1, ring1), c2, d)
                            azi_ange2, s_ang2 = find_angle((c2, ring2), c1, d)
                            ang = find_angle_between_plane((c1, ring1), (c2, ring2))
                            p_cz2 = pdb_data[m][(aro_res1[0], aro_res1[1], aro_res1[2], 'CZ2')]
                            ang_cz2 = find_angle3((p_cz2, c1), c2, d)
                            # print (aro_res1,aro_res2,d,ang,azi_ange1,s_ang1,azi_ange2,s_ang2)
                            fo.write('{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.
                                     format(pdb_id, m,aro_res1[0],aro_res1[1],aro_res1[2],aro_res2[0], aro_res2[1],
                                            aro_res2[2],sq_ifo[aro_res1[1]][0],sq_ifo[aro_res1[1]][1],
                                            sq_ifo[aro_res2[1]][0],sq_ifo[aro_res2[1]][1],d, ang,azi_ange1,s_ang1,
                                            azi_ange2,s_ang2,ang_cz2))
                        except KeyError:
                            pass
    fo.close()


if __name__ == "__main__":
    pdb_id = sys.argv[1]
    out_dir = sys.argv[2]
    calculate_interaction(pdb_id, out_dir=out_dir)
