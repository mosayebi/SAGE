from multiprocessing import Pool,  cpu_count, current_process
import contextlib
import hub_mp as hub
import sys
import time


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


mesh_flag = True 


#print hub.cylinder_form_factor(0.5, 0.5, 4.)
#sys.exit()


# min_timestep = 0
# max_timestep = 1e10
# traj_file = '/Users/mm15804/scratch/SAGE/psi3_test/dump_0.05.lammpstrj'
traj_data = hub.read_dump(traj_file, min_timestep, max_timestep)
traj_data = traj_data[-20:]
sq_file = traj_file+'.saxs.02'


print("\nNumber of cores available equals %d\n" % cpu_count())

if __name__ == "__main__":
    start = time.time()
    
    smin = 0.01
    smax = 7.0
    Ns   = 200
    ds   = (smax-smin)/(Ns-1) 

    futures=[]
    q = []

    #print get_saxs_intensity_mesh(float(0.1), traj_data[-1])
    
    with contextlib.closing( Pool() ) as pool:
        for i in range(len(traj_data)):
            for ii in range(Ns):
                s_mag = smin + ii*ds
                snap = traj_data[i]
                if (mesh_flag):
                    futures.append( pool.apply_async( hub.get_saxs_intensity_mesh, [s_mag, snap, True, 'sphere'] ) )
                else:
                    futures.append( pool.apply_async( hub.get_saxs_intensity, [s_mag, snap, True, 'sphere'] ) ) 
                    #get_saxs_intensity_mesh(s_mag, snap, my_model_flag=True, form='sphere') 
                    #  my_model_flag=True : will use my_model conf
                    #  form='sphere' will use a unit solid sphere form_factor for each particle in conf.
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
        print("(%d) running(queued)/done workers =  %d / %d" % (cnt, running, len(futures)-running)) 

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

N_mol = snap['N']/16 * 2
out = "#q I(q)\n"   
for key in sorted(sum_sq, key=float) :
    if (N_sq[key] > 0): out += "%s %s\n" % (key, abs(sum_sq[key])/N_sq[key]/N_mol) 

print("\n%s\n"% out)
f = open(sq_file, 'w')
f.write(out)
f.close()
print("\n s(q) is written to %s") % sq_file

