
import copy
import numpy as np
import isambard
import hub_mp_python3 as hub_mp
import sys, os


def rigid_transform_3D(A, B):
    # see http://nghiaho.com/?page_id=671
    assert len(A) == len(B)
    N = A.shape[0]; # total points
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    # centre the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))
    # dot is matrix multiplication for array
    H = np.transpose(AA) * BB
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T * U.T
    # special reflection case
    if np.linalg.det(R) < 0:
       #print "Reflection detected"
       Vt[2,:] *= -1
       R = Vt.T * U.T
    t = -R*centroid_A.T + centroid_B.T
    #print t
    return R, t


def transform_mon(mon, CG_com, CG_a):
    mon = copy.deepcopy(mon)
    coord0 = np.matrix( [mon[0].centre_of_mass, mon[0].primitive[-1]['CA'].array])
    h = np.linalg.norm ( mon[0].centre_of_mass - mon[0].primitive[-1]['CA'].array)
    coord1 = np.matrix( [CG_com, CG_com + CG_a*h])
    R, T = rigid_transform_3D( coord0, coord1)
    theta=np.arccos(0.5*(R[0,0]+R[1,1]+R[2,2] -1.))
    if theta>1e-10 :
        ax = np.array([(R[2,1]-R[1,2])/2/np.sin(theta), (R[0,2]-R[2,0])/2/np.sin(theta), (R[1,0]-R[0,1])/2/np.sin(theta)])
        T = np.array(T.T)[0]
        #print (ax,theta)
    else:
        print ('Error theta = %s' % theta )
        
    mon.rotate(theta, ax, radians=True)   
    mon.translate(T)
    return mon




if len(sys.argv) == 2 :
   traj_file = sys.argv[1]
   max_timestep = 1e10
   min_timestep = 0
elif len(sys.argv) == 4 :
   traj_file = sys.argv[1]
   min_timestep = float(sys.argv[2])
   max_timestep = float(sys.argv[3])
else:
    print ('Usage: %s dump_file [min_timestep max_timestep]' % (sys.argv[0]))
    sys.exit(1)



hub = isambard.ampal.Assembly()
hubA = isambard.ampal.convert_pdb_to_ampal(os.environ['HOME']+'/GoogleDrive/scratch/SAGE/scripts/SAGE/hubA.pdb')

dimerCC = hubA[0][:25]
trimerCC = hubA[0][25:50]
dimerCC.relabel_all()
trimerCC.relabel_all()
dimer_mon = isambard.ampal.Assembly()
dimer_mon.append(dimerCC)
trimer_mon = isambard.ampal.Assembly()
trimer_mon.append(trimerCC)


traj_data = hub_mp.read_dump(traj_file, min_timestep, max_timestep)
traj_data = traj_data[-1:] 


for step, snap in enumerate(traj_data):
    CGsnap = hub_mp.snap2CG(snap)
    sage_system = isambard.ampal.Assembly()
    for mol_id, mol in enumerate(CGsnap[:100]):   
        dimer = copy.deepcopy(dimer_mon)
        trimer = copy.deepcopy(trimer_mon)
        dimer = transform_mon(dimer, mol['di_x']*10, mol['di_a'] )
        trimer = transform_mon(trimer, mol['tri_x']*10, mol['tri_a'] )                                 
        sage_system.extend(dimer) 
        sage_system.extend(trimer)  
        print (mol_id)
    sage_system.relabel_all()
    break
with open('output.pdb', 'w') as outf:
    outf.write(sage_system.pdb)

