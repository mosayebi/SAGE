
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
   mols = -1
elif len(sys.argv) == 5 :
   traj_file = sys.argv[1]
   min_timestep = float(sys.argv[2])
   max_timestep = float(sys.argv[3])
   mols = float(sys.argv[4])
else:
    print ('Usage: %s dump_file [min_timestep max_timestep #molecules]' % (sys.argv[0]))
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
    for i in range(10): 
        sage_system = isambard.ampal.Assembly() 
        tcl_out = hub_mp.mols_tcl_header(snap['box'], False)
        cnt = 0
        print('%s mol_ids:' % i)
        for mol_id, mol in enumerate(CGsnap[0:int(mols)]): 
            if mol['visiblity'] != i : continue
            dimer = copy.deepcopy(dimer_mon)
            trimer = copy.deepcopy(trimer_mon)
            dimer = transform_mon(dimer, mol['di_x']*10, mol['di_a'] )
            trimer = transform_mon(trimer, mol['tri_x']*10, mol['tri_a'] )                                 
            sage_system.extend(dimer) 
            sage_system.extend(trimer)
            tcl_out += mol['di_tcl']  
            tcl_out += mol['tri_tcl']  
            print (mol_id, end=" ")
            cnt += 1
            #if (cnt>100): break
        if (cnt==0): break
        com_to_zero = - sage_system.centre_of_mass    
        sage_system.translate(com_to_zero)
        tcl_out += hub_mp.mols_tcl_footer()
        #sage_system.relabel_all() !can mess up with PDB for very large confs
        output = traj_file+'_'+str(i)+'_'+str(cnt)+'_'+str(snap['step'])+'.pdb'
        path='/'.join(traj_file.split('/')[:-1])
        if path=='': path='.'
        hub_mp.make_sure_path_exists(path+'/pdb')
        output = path+'/pdb/'+''.join(traj_file.split('/')[-1:])+'_'+str(i)+'_'+str(cnt)+'_'+str(snap['step'])+'.pdb'
        tcloutput = path+'/pdb/'+''.join(traj_file.split('/')[-1:])+'_'+str(i)+'_'+str(cnt)+'_'+str(snap['step'])+'.tcl'
        with open(output, 'w') as outf:
            outf.write(sage_system.pdb)
            print('\n(%s) %s molecules are converted to PDB' % (i, int(cnt)))
            print('The PDB file saved to %s' % output) 
        with open(tcloutput, 'w') as outf:
            outf.write(tcl_out)
            print('The tcl file saved to %s \n\n' % tcloutput) 
            hub_mp.render_tcl_file_LINUX (1024, 1024, png_file=output+".png", tcl_out_file=tcloutput)
        del sage_system         
    break



