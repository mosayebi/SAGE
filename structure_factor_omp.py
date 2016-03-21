from multiprocessing import Pool,  cpu_count, current_process
import contextlib
import hub_mp as hub
import time
import numpy as np
import sys
import os


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




def get_structure_factor_qv(xm, N_mol, box, qv):
    #q = np.array(random_unit_vector()) * qmod
    sq = 1.
    n_pairs = N_mol * (N_mol-1) / 2
    drs = np.zeros((n_pairs, 3))
    i = 0
    for mol1_id in range(N_mol):
      for mol2_id in range(mol1_id+1, N_mol):
        #if (mol1_id >= mol2_id): continue 
        dr = xm [ mol2_id, :] - xm [ mol1_id, :]
        drs[i,:] = dr
        i += 1
        # dr = hub.PBC (dr, box)
        # sq += 2*np.cos( np.dot(dr,qv) )
    
    drs = map(lambda x: hub.PBC(x, box), drs)
    sqs = 2*np.cos( np.sum(np.multiply(drs,qv), axis=1) )
    #sqs = map(lambda x: 2*np.cos(np.dot(x,qv)), drs)
    sq  += np.sum(sqs) / N_mol
    return sq  


def get_structure_factor_q(snap, qmod):
    np.random.seed()
    start = time.time()
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N_mol = snap['N']/16 * 2

    # xm =  np.zeros((N_mol,3))
    # for i in range(N_mol):
    #     xm [i,:] = x [ hub.get_helix_COM_atom_id(i), :]
    snap = snap2hubconf(snap)
    xm = snap['coords']
    N_mol = snap['N']

    Ndir = 30
    sum_sq = 0.0
    for idir in range(Ndir):
        q = np.array(random_unit_vector()) * qmod
        sum_sq += get_structure_factor_qv(xm, N_mol, box, q)        
    end = time.time()
    print("[trajectory timestep %s]: averaging s(q) for q = %s over %d directions for %d molec. took %s (s). {process %s}" \
        % (step, qmod, Ndir, N_mol, end-start, current_process().pid))
    return sum_sq/Ndir


def snap2hubconf(snap):
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N_mol = snap['N']/16*2
    N_hub = N_mol/6
    print "Nmol=%s , N_hub=%s"%(N_mol, N_hub)
    xm = np.zeros((N_hub,3))
    cnt = 0 
    for i in range(N_hub):  
        xx = [] 
        for p in range(3):
            while not p_type[cnt] == 14:
                cnt += 1
            #print cnt, i, p, p_type[cnt]
            xx.append(x[cnt,:]) 
            cnt += 1
        COM = [sum(p)/len(p) for p in zip(*xx)]   
        xm [i,:] = COM 
    snap['coords'] = xm
    snap['N'] = N_hub
    return snap



if len(sys.argv) == 2 :
   traj_file = sys.argv[1]
   max_timestep = 1e10
   min_timestep = 0
elif len(sys.argv) == 4 :
   traj_file = sys.argv[1]
   min_timestep = float(sys.argv[2])
   max_timestep = float(sys.argv[3])
else:
    print 'Usage: %s dump_file [min_timestep max_timestep]' % (sys.argv[0])
    sys.exit(1)




# min_timestep = 0
# max_timestep = 1e10
# traj_file = '/Users/mm15804/scratch/SAGE/psi3_test/dump_0.05.lammpstrj'
traj_data = hub.read_dump(traj_file, min_timestep, max_timestep)
traj_data = traj_data[-10:]
sq_file = traj_file+'.hubhub.sq.03'


print("\nNumber of cores available equals %d\n" % cpu_count())

if __name__ == "__main__":
    start = time.time()

    qmin = 0.01
    qmax = 7
    Nq   = 500
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



    #futures[-1].wait()
    cnt = 0
    while True:
        all_finished = True
        running = 0
        cnt += 1
        #print("\nHave the workers finished?")
        for i in range(0,len(futures)):
            #if futures[i].ready():
                #print("Worker %d has finished" % i)
            #else:
            if not futures[i].ready():
                all_finished = False
                running += 1
                #print("Worker %d is running..." % i)
        if all_finished:
            #print("All %d workers are done..."%len(futures))
            break
        time.sleep(10)
        print("(%d) running(queued) workers =  %d,     done = %d" % (cnt, running, len(futures)-running)) 

    end = time.time()
    print ("time: %s (s) for %d workers [%f]" %((end-start), len(futures), (end-start)/len(futures)))

    sum_sq = {}
    sum_sq2= {}
    N_sq = {}
    for i in range(0,len(futures)):
        if futures[i].successful():
            sq = futures[i].get()
            qmod = q[i]
            sum_sq [qmod] = (sum_sq.get(qmod, 0.0)) + sq
            sum_sq2[qmod] = (sum_sq2.get(qmod, 0.0)) + sq*sq
            N_sq [qmod] = (N_sq.get(qmod, 0)) + 1
        else:
            print("Worker %d failed!" % i)
            try:
                futures[i].get()
            except Exception as e:
                    print("Error = %s : %s" % (type(e), e))



N_mol = snap['N']/16 * 2
out = "#q s(q) err N\n"   
for key in sorted(sum_sq, key=float) :
    if N_sq[key]>0 : 
       avg = sum_sq[key]/N_sq[key]
       avg2 = sum_sq2[key]/N_sq[key]
       out += "%s %s %s %s\n" % (key, avg, np.sqrt((avg2 - avg*avg)/N_sq[key]), N_sq[key] ) 

print("\n%s\n"% out)
f = open(sq_file, 'w')
f.write(out)
f.close()
print("\n s(q) is written to %s") % sq_file

