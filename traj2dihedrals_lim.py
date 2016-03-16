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
traj_data = traj_data[-48:]


print("\nNumber of cores available equals %d\n" % cpu_count())

if __name__ == "__main__":
    start = time.time()

    futures = []
    with contextlib.closing( Pool() ) as pool:
        for s, snap in enumerate(traj_data): 
             futures.append( pool.apply_async( hub.snap2dihedrals_all, [snap] ) ) 

    

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





    print "\ndata collection took %f (s)\n" % (end-start)
    phi_mm=[]
    theta1_mm=[]
    theta2_mm=[]
    phi_hh=[]
    theta1_hh=[]
    theta2_hh=[]
    for i in range(0,len(futures)):
        if futures[i].successful():
            snap_phi_mm, snap_theta1_mm, snap_theta2_mm, snap_phi_hh, snap_theta1_hh, snap_theta2_hh = futures[i].get()
            phi_mm += snap_phi_mm
            theta1_mm += snap_theta1_mm
            theta2_mm += snap_theta2_mm        
            phi_hh += snap_phi_hh
            theta1_hh += snap_theta1_hh
            theta2_hh += snap_theta2_hh
        else:
            print("Worker %d failed!" % i)
            try:
                futures[i].get()
            except Exception as e:
                print("Error = %s : %s" % (type(e), e))

    if len(phi_hh)> 0 :	
        filename  =  traj_file+'.hubhub_dihedrals_hist_lim.pdf'         
        hub.plot_hist(phi_hh, theta1_hh, theta2_hh, filename, x_lim=[-30,30] ) 
        #print phi, theta1, theta2
        #plot_scatter_hist(phi, theta1)
        #plot_scatter_hist_sns(phi, theta1)
    else:
        print "no hub-hub angle to plot histograms"   
    
    if len(phi_mm)> 0 :	
        filename  =  traj_file+'.molmol_dihedrals_hist_lim.pdf'         
        hub.plot_hist(phi_mm, theta1_mm, theta2_mm, filename, x_lim=[-30,30] ) 
        #print phi, theta1, theta2
        #plot_scatter_hist(phi, theta1)
        #plot_scatter_hist_sns(phi, theta1)
    else:
        print "no mol_mol angle to plot histograms"




    # plot atomistic dihedrals
    atom_phi, atom_theta1, atom_theta2 = hub.read_atomistic_angles('/home/mm15804/SAGE/data/hub-hub_angles_lim.his')
    hub.plot_hist(atom_phi, atom_theta1, atom_theta2, filename='atomistic_hubhub_dihedrals_hist.pdf', x_lim=[-30,30])

    atom_phi, atom_theta1, atom_theta2 = hub.read_atomistic_angles('/home/mm15804/SAGE/data/hub-dimer_angles.his')
    hub.plot_hist(atom_phi, atom_theta1, atom_theta2, filename='atomistic_molmol_dihedrals_hist_lim.pdf', x_lim=[-30, 30])

