import logging
import requests
import get_data
import numpy
import os
import sys

aromatic_atoms = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2', 'xx', 'yy']  # if needed un comment
    }
ring_atoms = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', ]  # if needed un comment
    }

def get_distance(c1, c2):
    """
    Calculates the distance between two coordinate points
    :param c1: array of x,y,z
    :param c2: array of x,y,z
    :return: distance between two ponts
    """
    return numpy.linalg.norm(c1 - c2)


def get_centroid(p):
    # print (len(p),p)
    x = [i[0] for i in p]
    y = [i[1] for i in p]
    z = [i[2] for i in p]
    c = [sum(x) / len(p), sum(y) / len(p), sum(z) / len(p)]
    return numpy.array(c)


def find_angle( p, c, d):
    pc,p = p
    v1 = p[1] - pc
    v2 = p[2] - pc
    nv = numpy.cross(v1, v2)
    #nv2 = numpy.cross(v2, v1)
    cv = c - pc
    #cv2 = c - cn
    nnv = nv / numpy.linalg.norm(nv)
    ncv = cv / numpy.linalg.norm(cv)
    #ncv2 = cv2 / numpy.linalg.norm(cv2)
    dp = abs(numpy.dot(nnv, ncv))
    #dp2 = abs(numpy.dot(nnv, ncv2))
    ang = numpy.arccos(dp)
    #ang2 = numpy.arccos(dp2)
    ang_deg = (180 / numpy.pi) * ang
    #ang_deg2 = (180 / numpy.pi) * ang2
    # print(ang_deg)
    s_ang = solid_angle(ang_deg, d)
    return ang_deg, s_ang


def find_angle_between_plane( p1, p2):
    pc1,p1 = p1
    pc2,p2=p2
    v11 = p1[1] - pc1
    v12 = p1[2] - pc1
    nv1 = numpy.cross(v11, v12)
    v21 = p2[1] - pc2
    v22 = p2[2] - pc2
    nv2 = numpy.cross(v21, v22)
    #nv2 = numpy.cross(v2, v1)
    #cv = c - pc
    #cv2 = c - cn
    nnv = nv1 / numpy.linalg.norm(nv1)
    ncv = nv2 / numpy.linalg.norm(nv2)
    #ncv2 = cv2 / numpy.linalg.norm(cv2)
    dp = abs(numpy.dot(nnv, ncv))
    #dp2 = abs(numpy.dot(nnv, ncv2))
    ang = numpy.arccos(dp)
    #ang2 = numpy.arccos(dp2)
    ang_deg = (180 / numpy.pi) * ang
    #ang_deg2 = (180 / numpy.pi) * ang2
    # print(ang_deg)
    return ang_deg


def solid_angle(a_deg, r):
    s = 1.4 #C-C bond length
    A=((3.0*numpy.sqrt(3))/2.0)*s*s #area of the hexagonal plane
    a = (numpy.pi / 180) * a_deg
    # sa2=2*numpy.pi*(1.0-1.0/(numpy.sqrt(1+(A*numpy.cos(a)/(numpy.pi*r1**r1)))))
    sa = 2 * numpy.pi * (1.0 - 1.0 / (numpy.sqrt(1 + (A*numpy.cos(a) / (numpy.pi * r * r)))))
    # print (a_deg)
    sa_deg = (180.0 / numpy.pi) * sa
    # sa_deg2 = (180.0 / numpy.pi) * sa2
    # print (a_deg,sa_deg,sa_deg2)
    return sa_deg

def get_aromatic_info(pdb_data):
    aromatic_atoms = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2', 'xx', 'yy']  # if needed un comment
    }
    ring_atoms = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', ]  # if needed un comment
    }
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
            if len(p)>1:
                c = get_centroid(p)
                ar_info[m][(str(ar_res[0]), ar_res[1], ar_res[2])] = (c, p)
    return ar_info


def calculate_interaction(pdb_id):
    pdb_data,aromatic_residues = get_data.get_pdb_data(pdb_id)
    for m in pdb_data:
        for r1 in range(len(aromatic_residues[m])):
            for r2 in range(r1,len(aromatic_residues[m])):
                aro_res1 = aromatic_residues[m][r1]
                aro_res2 = aromatic_residues[m][r2]
                if aro_res1 != aro_res2:
                    ring1 = [pdb_data[m][(aro_res1[0],aro_res1[1],aro_res1[2],i)] for i in ring_atoms[aro_res1[2]]]
                    ring2 = [pdb_data[m][(aro_res2[0], aro_res2[1], aro_res2[2], i)] for i in ring_atoms[aro_res2[2]]]
                    c1 = get_centroid(ring1)
                    c2 = get_centroid(ring2)
                    d = get_distance(c1,c2)
                    azi_ange1, s_ang1 = find_angle((c1,ring1),c2,d)
                    azi_ange2, s_ang2 = find_angle((c2,ring2), c1, d)
                    ang = find_angle_between_plane((c1,ring1),(c2,ring2))
                    print (aro_res1,aro_res2,d,ang,azi_ange1,s_ang1,azi_ange2,s_ang2)


def run_on_nmrbox():
        id_pair=[]
        url = "http://api.bmrb.io/v2/mappings/bmrb/pdb?match_type=exact"
        r = requests.get(url).json()
        for ids_dict in r:
            bmrb_id = ids_dict['bmrb_id']
            pdb_ids = ids_dict['pdb_ids']
            for k in pdb_ids:
                id_pair.append((bmrb_id,k))
        for x in id_pair:
            print ('Calculating {} {}'.format(x[1],x[0]))
            calculate_interaction(x[1],x[0],nmrbox=True)


if __name__=="__main__":
    calculate_interaction('1AWV')
