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


# traj_file = '/Users/mm15804/scratch/SAGE/psi3_test/dump_0.05.lammpstrj'
traj_data = hub.read_dump(traj_file, min_timestep, max_timestep)
#traj_data = traj_data[-10:]


print("\nNumber of cores available equals %d\n" % cpu_count())

if __name__ == "__main__":
    start = time.time()

    futures = []
    with contextlib.closing( Pool() ) as pool:
        for s, snap in enumerate(traj_data):
             futures.append( pool.apply_async( hub.conf2tcl, [snap, False, True, False, []] ) ) 

    

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
        print("(%d) running(queued)/total workers =  %d / %d" % (cnt, running, len(futures))) 

    end = time.time()
    print ("time: %s (s) for %d workers [%f]" %((end-start), len(futures), (end-start)/len(futures)))





    print "\ntcl output generation took %f (s)\n" % (end-start)
    hub.make_sure_path_exists('tcl')
    for i in range(0,len(futures)):
        if futures[i].successful():
            tcl_out = futures[i].get()
            snap = traj_data[i]
            tcl_file = 'tcl/'+traj_file+'.'+str(snap['step'])+".tcl"
            hub.write_tcl_out (tcl_out, filename=tcl_file)
            hub.render_tcl_file_LINUX (1024, 1024, png_file=tcl_file+".png", tcl_out_file=tcl_file)
        else:
            print("Worker %d failed!" % i)
            try:
                futures[i].get()
            except Exception as e:
                print("Error = %s : %s" % (type(e), e))
    print "%s tcl files are written\n" % (len(futures))
            
