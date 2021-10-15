import sys,os

def get_pml(pdb,r1,c1,r2,c2):
    fo = open('/Users/kumaran/aro_int/img/{}_{}_{}.pml'.format(pdb,r1,r2),'w')
    fo.write('reinitialize\n')
    fo.write('set cif_use_auth, off\n')
    fo.write('fetch {}\n'.format(pdb))
    fo.write('color yellow,resi {} and chain {}\n'.format(r1,c1))
    fo.write('color red,resi {} and chain {}\n'.format(r2,c2))
    #fo.write('show spheres, resi {} and name H\n'.format(r1))
    fo.write('show sticks, resi {} and chain {}\n'.format(r1,c1))
    fo.write('show sticks, resi {} and chain {}\n'.format(r2, c2))
    fo.write('zoom resi {} and chain {}+ resi {} and chain {}\n'.format(r1,c1,r2,c2))
    fo.write('png /Users/kumaran/aro_int/img/{}_{}_{}.png, width=10cm, dpi=300, ray=0\n'.format(pdb,r1,r2))
    fo.write('set all_states, on\n')
    fo.write('png /Users/kumaran/aro_int/img/{}_{}_{}_all.png, width=10cm, dpi=300, ray=0'.format(pdb, r1, r2))
    fo.close()
    return '/Users/kumaran/aro_int/img/{}_{}_{}.pml'.format(pdb,r1,r2)

if __name__ == "__main__":
    pdb = sys.argv[1]
    r1 = sys.argv[2]
    r2 = sys.argv[4]
    c1 = sys.argv[3]
    c2 = sys.argv[5]
    fname = get_pml(pdb,r1,c1,r2,c2)
    os.system('/Applications/PyMOL\ 2.app/Contents/bin/pymol -c {}'.format(fname))
