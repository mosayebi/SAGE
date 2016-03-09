import hub_mp as hub
import sys


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


traj_data = hub.read_dump(traj_file, min_timestep, max_timestep)
if __name__ == "__main__":
    for i in range(len(traj_data)):
        snap = traj_data[i]
        step = snap ['step']
        PDB_file = traj_file+'.'+str(step)+'.pdb'
        hub.conf2pdb_saxs(snap, filename)


