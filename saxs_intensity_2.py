from multiprocessing import Pool,  cpu_count, current_process
import contextlib
import hub_mp as hub
import sys
import time
import numpy as np


if len(sys.argv) == 2 :
   traj_file = sys.argv[1]
else:
    print 'Usage: %s sesh_SAGE' % (sys.argv[0])
    sys.exit(1)


print("\nNumber of cores available equals %d\n" % cpu_count())

if __name__ == "__main__":
    start = time.time()
    
    smin = 0.01
    smax = 8.0
    Ns   = 100
    ds   = (smax-smin)/(Ns-1) 

    futures=[]
    q = []


    mesh_flag = True 
    snap = hub.read_sesh_SAGE(traj_file)
    sq_file = traj_file+'.saxs'
    snap = hub.make_sc_sheet (4, 100)
    sq_file ='sc_sheet.saxs'    
    snap = hub.make_sc_cube (4, 40)
    sq_file ='sc_cube.saxs'
    snap['coords']=np.array([[1., 0, 0],[2., 0, 0]])
    snap['N']=1
    sq_file ='single_sphere.saxs'

    
    # print hub.get_saxs_intensity_mesh(0.5, snap, False)
    # sys.exit()

    with contextlib.closing( Pool() ) as pool:
            for ii in range(Ns):
                s_mag = smin + ii*ds
                if (mesh_flag):
                    futures.append( pool.apply_async( hub.get_saxs_intensity_mesh, [s_mag, snap, False] ) )
                else:
                    futures.append( pool.apply_async( hub.get_saxs_intensity, [s_mag, snap, False] ) )    
                q.append(s_mag)



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

N_mol = snap['N']
out = "#q I(q)\n"   
for key in sorted(sum_sq, key=float) :
    if (N_sq[key] > 0): out += "%s %s\n" % (key, abs(sum_sq[key])/N_sq[key]/N_mol) 

print("\n%s\n"% out)
f = open(sq_file, 'w')
f.write(out)
f.close()
print("\n s(q) is written to %s") % sq_file

