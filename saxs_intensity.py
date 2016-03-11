from multiprocessing import Pool,  cpu_count, current_process
import contextlib
import hub_mp as hub
import time
import numpy as np
import scipy
from scipy import special
import sys
import os


# calculates saxs intensity based on J Appl Cryst. (1995) 28, 768-773

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
    return np.array([x, y, z])


def Cartesian2Spherical(xyz):
    sph = np.zeros(xyz.shape)
    xy = xyz[:,0]**2 + xyz[:,1]**2
    sph[:,0] = np.sqrt(xy + xyz[:,2]**2)
    sph[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    #sph[:,1] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    sph[:,2] = np.arctan2(xyz[:,1], xyz[:,0])  # [-pi, pi] ??
    return sph

def sphere_form_factor(s, R):
    f = 9 * (  (np.sin(s*R) - s*R*np.cos(s*R)) / (s*R)**3 )**2
    return f

def form_factor(s):
    f =  sphere_form_factor(s, 1.0)
    return f

def A_lm(s_mag, spherical_coords, l, m):
    Alm_vector = form_factor(s_mag) * scipy.special.jn(l, s_mag * spherical_coords[:,0]) * \
                 np.conj ( scipy.special.sph_harm(m, l, spherical_coords[:,2], spherical_coords[:,1]) )
    Alm = 4.0 * np.pi * (complex(0,1)**l) * np.sum (Alm_vector) 
    return Alm

def get_Aa2(s_vec, spherical_coords, lmax):
    A = 0
    for l in range(lmax+1):
        for m in range(-l,l+1):
            A += A_lm( s_vec[0], spherical_coords, l, m) * scipy.special.sph_harm(m, l, s_vec[2], s_vec[1])
    return A * np.conj(A)

def get_saxs_intensity(s_mag, snap):
    start = time.time()
    lmax = 10
    Ndir = 10
    x = snap['coords']
    N_mol = snap['N']/16 * 2
    step = snap['step']

    spherical_coords = Cartesian2Spherical(x)
    
    sum_I = 0.0
    for i in range(Ndir):
        s = random_unit_vector() * s_mag
        s_vec = Cartesian2Spherical(np.array([s]))[0]
        sum_I += get_Aa2(s_vec, spherical_coords, lmax)
    end = time.time()
    print("[trajectory timestep %s]: averaging I(s) for s = %s over %d directions for %d molec. took %s (s). {process %s}" \
        % (step, s_mag, Ndir, N_mol, end-start, current_process().pid))
    return sum_I/Ndir




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
sq_file = traj_file+'.sq'


print("\nNumber of cores available equals %d\n" % cpu_count())

if __name__ == "__main__":
    start = time.time()

    smin = 0.01
    smax = 3.0
    Ns   = 80
    ds   = (smax-smin)/(Ns-1) 

    futures=[]
    q = []

    print get_saxs_intensity(float(0.1), traj_data[-1])
    
    with contextlib.closing( Pool() ) as pool:
        for i in range(len(traj_data)):
            for ii in range(Ns):
                s_mag = smin + ii*ds
                snap = traj_data[i]
                futures.append( pool.apply_async( get_saxs_intensity, [s_mag, snap] ) )
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

N_mol = snap['N']/16 * 2
out = "#q s(q)\n"   
for key in sorted(sum_sq, key=float) :
    if (N_sq[key] > 0): out += "%s %s\n" % (key, abs(sum_sq[key])/N_sq[key]/N_mol) 

print("\n%s\n"% out)
f = open(sq_file, 'w')
f.write(out)
f.close()
print("\n s(q) is written to %s") % sq_file

