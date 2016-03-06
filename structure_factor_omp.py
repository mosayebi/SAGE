from multiprocessing import Pool,  cpu_count, current_process
import contextlib
import hub_mp as hub
import time
import numpy as np


def random_unit_vector():
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return (x,y,z)


def get_structure_factor_q(snap, qmod):
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N_mol = snap['N']/16 * 2
    #print("Process %s going to get s(q) for q=%s and snapshot=%s" % (current_process().pid, qmod, step))
    Ndir = 40
    sum_sq = 0. + 0.j
    for idir in range(Ndir):
        q = -1j * np.array(random_unit_vector()) * qmod
        sq = 1.0
        for mol1_id in range(N_mol):
          for mol2_id in range(mol1_id+1, N_mol):
            if (mol1_id == mol2_id): continue 
            r1 = x [ hub.get_helix_COM_atom_id(mol1_id), :]
            r2 = x [ hub.get_helix_COM_atom_id(mol2_id), :]
            dr = hub.PBC (r2-r1, box)
            sq += np.exp( np.dot(dr,q) )
        sum_sq += sq
    #print sum_sq/Ndir    
    return sum_sq/Ndir







min_timestep = 0
max_timestep = 1e10
traj_file = '/Users/mm15804/scratch/SAGE/psi3_test/dump_0.05.lammpstrj'
traj_data = hub.read_dump(traj_file, min_timestep, max_timestep)



print("\nNumber of cores available equals %d" % cpu_count())

       



if __name__ == "__main__":

    qmin = 0.01
    qmax = 1.0
    Nq   = 10
    dq   = (qmax-qmin)/(Nq-1) 

    futures=[]
    q = []
    with contextlib.closing( Pool() ) as pool:
        for i in range(len(traj_data)):
            for iq in range(Nq):
                qmod = qmin + iq*dq
                snap = traj_data[i]
                futures.append( pool.apply_async( get_structure_factor_q, [snap, qmod] ) )
                q.append(qmod)

    futures[-1].wait()
    while True:
        all_finished = True
        print("\nHave the workers finished?")
        for i in range(0,len(futures)):
            #if futures[i].ready():
                #print("Worker %d has finished" % i)
            #else:
            if not futures[i].ready():
                all_finished = False
                print("Worker %d is running..." % i)
        if all_finished:
            print("All workers are done...")
            break
        time.sleep(2)

    sum_sq = {}
    N_sq = {}
    for i in range(0,len(futures)):
        if futures[i].successful():
            sq = futures[i].get()
            qmod = q[i]
            sum_sq [qmod] = (sum_sq.get(qmod, 0.0)) + sq
            N_sq [qmod] = (N_sq.get(qmod, 0)) + 1
        else:
            print("Worker %d failed!" % i)
            try:
                futures[i].get()
            except Exception as e:
                    print("Error = %s : %s" % (type(e), e))

   
for key in sorted(sum_sq, key=float) :
    if (N_sq[key] > 0): print (key, abs(sum_sq[key]/N_sq[key]))