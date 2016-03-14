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


tcl_write_flag=True
psi3_angle_hist_flag=True
skip_snap = False
psi3_file = traj_file+'.psi3'


last_timestep = 0  
# get the last snapshot already written to psi3_file
if (not skip_snap): 
    try:
        f = open (psi3_file, 'w')
    except:
        print ("traj2psi3.py: cannot open output file. Exit")
        sys.exit(1)
else:        
    try:
        f = open (psi3_file, 'a+')
    except:
        print ("traj2psi3.py: cannot open output file. Exit")
        sys.exit(1)
    f.seek(0)
    lines=f.readlines()
    if len(lines)> 1 : 
        line = lines[-1].strip('\n')
        print '\n#last line in %s file is:'
        print line
        if line.split()[0].isdigit() : 
            last_timestep = float(line.split()[0]) 
            print '#skip_snap is True. skipping timesteps smaller or equal than %d\n\n' % last_timestep  
        else:
            skip_snap=False 
            print '\n#Warning: last line of %s is not as expected. Setting skip_snap to False.\n\n'% psi3_file






# traj_file = '/Users/mm15804/scratch/SAGE/psi3_test/dump_0.05.lammpstrj'
traj_data = hub.read_dump(traj_file, min_timestep, max_timestep)
#traj_data = traj_data[-60:]


print("\nNumber of cores available equals %d\n" % cpu_count())

if __name__ == "__main__":
    start = time.time()

    futures = []
    with contextlib.closing( Pool() ) as pool:
        for s, snap in enumerate(traj_data): 
             if (skip_snap and snap['step'] <= last_timestep ): 
                print 'skipping %s' % snap['step']
             else :   
                futures.append( pool.apply_async( hub.snap2psi3, \
                    [snap, skip_snap, tcl_write_flag, psi3_angle_hist_flag] ) ) 
   

    #futures[-1].wait()
    cnt = 0
    while True:
        all_finished = True
        running = 0
        cnt += 1
        #print("\nHave the workers finished?")
        for i in range(0,len(futures)):
            if futures[i].ready():
                out, psi3_vec, angles = futures[i].get()
                print("Worker %d has finished with following output %s" % (i, out))
            else:
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



    print "\ndata collection took %f (s)\n" % (end-start)


    all_angles=[]
    output= ''
    for i in range(0,len(futures)):
        if futures[i].successful():
            out, psi3_vec, angles = futures[i].get()
            all_angles += angles
            output += out
        else:
            print("Worker %d failed!" % i)
            try:
                futures[i].get()
            except Exception as e:
                print("Error = %s : %s" % (type(e), e))


    if (len(all_angles)>0) :
        make_sure_path_exists('psi3_angle_hist')
        pdf_file = 'psi3_angle_hist/'+traj_file+'.total_hist.'+str(traj_data[0]['step'])+'_'+str(traj_data[-1]['step'])+".pdf"
        hub.plot_1D_hist (all_angles, pdf_file) 


    if not output=='' :
        output = "#"+time.strftime("%c")+"\n#step psi3 N_angles count(arms>=2)\n"+output
        f.write(output)
        print "psi3 data is written to %s" % psi3_file

        


        

 


