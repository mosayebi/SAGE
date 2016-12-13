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
hub.traj2movie(traj_data, movie_file=traj_file+'.movie', cluster_flag=False, rotation_flag=True)

            
