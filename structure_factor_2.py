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
    sq = 1
    drs = []
    for mol1_id in range(N_mol):
      for mol2_id in range(mol1_id+1, N_mol):
        if (mol1_id >= mol2_id): continue 
        dr = xm [ mol2_id, :] - xm [ mol1_id, :]
        drs.append(dr)
        # dr = hub.PBC (dr, box)
        # sq += 2*np.cos( np.dot(dr,qv) )
    
    drs = map(lambda x: hub.PBC(x, box), drs)
    sqs = map(lambda x: 2*np.cos(np.dot(x,qv)), drs)
    sq  += np.sum(sqs) / N_mol
    return sq  


def get_structure_factor_q(snap, qmod):
    np.random.seed()

    start = time.time()
    xm = snap['coords']
    box = snap['box']
    step = snap['step']
    N_mol = snap['N']

    Ndir = 20
    sum_sq = 0.0
    for idir in range(Ndir):
        q = np.array(random_unit_vector()) * qmod
        sum_sq += get_structure_factor_qv(xm, N_mol, box, q)        
    end = time.time()
    print("[trajectory timestep %s]: averaging s(q) for q = %s over %d directions for %d molec. took %s (s). {process %s}" \
        % (step, qmod, Ndir, N_mol, end-start, current_process().pid))
    return sum_sq/Ndir






if len(sys.argv) == 3 :
   traj_file = sys.argv[1]
   mode = int(sys.argv[2])
elif len(sys.argv) == 2 :
   traj_file = sys.argv[1]
   mode = 0
else:
    print 'Usage: %s dump_file mode' % sys.argv[0]
    sys.exit(1)


print "mode %s"%mode
if mode==0:
   snap = hub.read_sesh_SAGE(traj_file)
   sq_file = traj_file+'.sq'
elif mode==1:
   url = 'http://neilsloane.com/ICOSP/ipack.3.932.txt'
   file = 'ipack.3.932.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'ipack.3.932.txt.sq'
elif mode==2:
   url = 'http://neilsloane.com/ICOSP/icover.3.932.7.4.txt'
   file = 'ipack.3.932.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'ipack.3.932.txt.sq'  
elif mode==3:
   url = 'http://neilsloane.com/ICOSP/ivol.3.932.7.4.txt'
   file = 'ivol.3.932.7.4.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'iivol.3.932.7.4.txt.sq'
elif mode==4:
   url = 'http://neilsloane.com/ICOSP/icover.3.482.4.4.txt'
   file = 'icover.3.482.4.4.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'icover.3.482.4.4.txt.sq'  
elif mode==5:
   snap = hub.make_sc_cube (4, 37)
   sq_file ='sc_cube.sq'
elif mode==6:
   snap = hub.make_sc_sheet (4, 100)
   sq_file ='sc_sheet.sq' 
elif mode==7:
   snap = hub.make_random_spherical_shell (32.1451, 620)
   sq_file ='random_spherical_shell_32.14_620.sq'
elif mode==8:
   url = 'http://neilsloane.com/ICOSP/icover.3.3002.10.10.txt'
   file = 'icover.3.3002.10.10.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'icover.3.3002.10.10.txt.sq'  
elif mode==9:
   url = 'http://neilsloane.com/ICOSP/icover.3.3242.18.0.txt'
   file = 'icover.3.3242.18.0.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'icover.3.3242.18.0.txt.sq'  
elif mode==10:
   url = 'http://neilsloane.com/ICOSP/icover.3.312.5.1.txt'
   file = 'icover.3.312.5.1.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'icover.3.312.5.1.txt.sq'
elif mode==11:
   url = 'http://neilsloane.com/ICOSP/icover.3.192.3.2.txt'
   file = 'iicover.3.192.3.2.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'icover.3.192.3.2.txt.sq' 
elif mode==12:
   url = 'http://neilsloane.com/ICOSP/icover.3.192.3.2.txt'
   file = 'icover.3.192.3.2.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'icover.3.192.3.2.txt.sq'  
elif mode==13:
   url = 'http://neilsloane.com/ICOSP/icover.3.912.6.5.txt'
   file = 'icover.3.912.6.5.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'icover.3.912.6.5.txt.sq' 
elif mode==14:
   snap = hub.make_random_spherical_shell (32.1451, 1000)
   sq_file ='random_spherical_shell_32.14_1000.sq'           
elif mode==15:
   snap = hub.read_sesh_SAGE(traj_file)
   x = snap ['coords']
   xt1 = np.zeros(x.shape)
   Translate = np.array([1000, 0, 0])
   for i in range(len(x)):
       xt1[i,:] = x[i,:] + Translate 
   xt2 = np.zeros(x.shape)
   Translate = np.array([0, 1500, 0])
   for i in range(len(x)):
       xt2[i,:] = x[i,:] + Translate 
   xt3 = np.zeros(x.shape)
   Translate = np.array([0, 0, 2800])
   for i in range(len(x)):
       xt3[i,:] = x[i,:] + Translate             
   xm = np.concatenate((x, xt1, xt2, xt3), axis=0)
   snap = {}
   snap ['coords'] = xm
   snap ['N'] = len(xm)
   snap ['step'] = '4_sesh_SAGE'
   snap ['traj'] =  traj_file
   snap['box'] = np.array ( [10*np.max(xm), 10*np.max(xm), 10*np.max(xm)] )
   sq_file = '4_'+traj_file+'.sq'

elif mode==16:
   snap = hub.make_random_spherical_shell (32.1451, 1860)
   sq_file ='random_spherical_shell_32.14_1860.sq' 

elif mode==17:
   url = 'http://neilsloane.com/ICOSP/ivol.3.1212.11.0.txt'
   file = 'ivol.3.1212.11.0.txt'
   snap = hub.read_spherical_packings_www_from_file(file)
   sq_file = 'ivol.3.1212.11.0.txt.sq'

elif mode == 18 :
    file='/projects/t3/mm15804/SAGE/model_5.3/hub_assembly/dump_1.50.lammpstrj'
    traj_data = hub.read_dump(file, 1, 1e10)
    snap = traj_data[-1]
    sq_file = 'model51hub_dump_1.50.lammpstrj.sq.03'

else:
   print "Error"       


# snap = hub.make_sc_sheet (4, 100)
# sq_file ='sc_sheet.saxs'    
# snap = hub.make_sc_cube (4, 40)
# sq_file ='sc_cube.saxs'
# snap['coords']=np.array([[1., 0, 0]])
# snap['N']=1
# sq_file ='single_sphere.saxs'
# snap['coords']=np.array([[0., 0, 0]])
# snap['N']=1
# sq_file ='single_cylindar.saxs'



print("\nNumber of cores available equals %d\n" % cpu_count())

if __name__ == "__main__":
    start = time.time()

    qmin = 0.01
    qmax = 7
    Nq   = 1000
    dq   = (qmax-qmin)/(Nq-1) 

    futures=[]
    q = []
    
    with contextlib.closing( Pool() ) as pool:
            for iq in range(Nq):
                qmod = qmin + iq*dq
                for replica in range(5):   
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



N_mol = snap['N']
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

