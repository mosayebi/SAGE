import os
import sys
import fnmatch
import numpy as np
#import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
#import pandas as pd
#from scipy import stats, integrate
#from scipy.stats import kendalltau
#mlab.use('Agg')  # Force mpl to use Tk backend
#import seaborn as sns


def read_atomistic_angles(filer):
  phi=[]
  theta=[]
  with open(filer,'r') as f:
      while True:
        line=f.readline().strip('\n')
        line=f.readline().strip('\n').split()
        while line:           
          if line[0] == 'Time' :
            print line[2], len(phi)
            if float(line[2])>=500: return phi, theta
            line=f.readline().strip('\n').split()
            continue
          phi.append(float(line[3]))
          theta.append(float( line[4]) + float(line[5] ) - 180.0 )
          line=f.readline().strip('\n').split()
      print len(phi)    
      return phi, theta        



def read_dump(filer, max_step):
   data=[]  
   snapshot = {}
   snapshot['traj'] = filer
   with open(filer,'r') as f:
      while True:
        line=f.readline().strip('\n')
        if not line: break
        items = line.split(' ')
        if items[0] == 'ITEM:':
            if items[1] == 'TIMESTEP':
                step = int(f.readline().split(' ')[0])
                if step > float(max_step) :
                    print 'max_step reached (%s)' % max_step 
                    break
                #print 'reading TIMESTEP %d' % step  
                snapshot['step'] =  step 
            if items[1] == 'NUMBER':
                N = int(f.readline().split(' ')[0]) 
                snapshot['N'] =  N
            if items[1] == 'BOX':
                line = f.readline().split(' ')
                box_x = float(line[1]) - float(line[0])
                line = f.readline().split(' ')
                box_y = float(line[1]) - float(line[0])
                line = f.readline().split(' ')
                box_z = float(line[1]) - float(line[0])
                snapshot['box'] = np.array([box_x, box_y, box_z])
            if items[1] == 'ATOMS': 
                p_type = np.zeros(N, dtype=np.int)
                x = np.zeros((N,3))
                for i in range(N):
                      line = f.readline().split(' ')
                      p_type[ int(line[0])-1 ]     = int  (line[1])
                      x     [ int(line[0])-1 ] [0] = float(line[2])
                      x     [ int(line[0])-1 ] [1] = float(line[3])
                      x     [ int(line[0])-1 ] [2] = float(line[4])
                      snapshot['p_type'] = p_type  
                      snapshot['coords'] = x
                data.append(snapshot.copy())
                #print snapshot
                snapshot = {}
   print 'From %s, last TIMESTEP %d' % (filer, step) 
   return data                   


def read_topology(filer, N=None):
   if N is None: 
      print 'N is not set.'
      sys.exit(2)
   with open(filer,'r') as f:
      while True:
        line=f.readline().strip('\n')
        if not line : break
        if line.split(' ')[0] == 'HUB':
            f.readline()
            print 'reading nb list data ...'
            nb = []  
            for i in range(N/8):
                nb.append(  [int(x)-1 for x in f.readline().strip('\n').split()]  )
            return nb    
        if line.split(' ')[0] == 'End': break

        

def PBC(d, box):
   d[0] = d[0] - np.rint(d[0]/box[0])*box[0]
   d[1] = d[1] - np.rint(d[1]/box[1])*box[1]
   d[2] = d[2] - np.rint(d[2]/box[2])*box[2]
   return d



def get_points(hub1_id, hub2_id, x, p_type, box):
    i1 = (hub1_id)*8 + 0
    i2 = (hub1_id)*8 + 4
    i3 = (hub2_id)*8 + 0
    i4 = (hub2_id)*8 + 4
    p1 =  box * x[i1, :]
    p2 =  box * x[i2, :]
    pp1 = box * x[i3, :]
    pp2 = box * x[i4, :]
    
    #print p_type[i1], p_type[i2], p_type[i3],p_type[i4] 

    x1 = p1
    x2 = p1 + PBC(p2 - p1, box)/2
    x3 = pp1 + PBC(pp2 - pp1, box)/2
    x4 = pp1
    return x1, x2, x3, x4

def get_angles(x1,x2,x3,x4,box):
    b1 = (PBC(x2 - x1, box)) / np.linalg.norm (PBC(x2 - x1, box))
    b2 = (PBC(x3 - x2, box)) / np.linalg.norm (PBC(x3 - x2, box))
    b3 = (PBC(x4 - x3, box)) / np.linalg.norm (PBC(x4 - x3, box))

    b12 = np.cross(b1,b2)
    b23 = np.cross(b2,b3)
             
    dihedral = np.degrees ( np.arctan2 ( np.dot( np.cross(b12, b23), b2) ,    np.dot(b12,b23) )  ) 
    bend1 = np.degrees( np.arccos (np.dot(b1,b2))) + np.degrees( np.arccos (np.dot(b2,b3))) - 180
    bend2 = np.degrees( np.arccos (np.dot(b1,-b3)))

    return dihedral, bend1, bend2



def if_bonded(mol1_id, mol2_id, x, rc, box):
    bonded = False
    for i in range(3):
            j = i 
            ii = (mol1_id)*6 + i + 4 - 1
            jj = (mol2_id)*6 + j + 4 - 1
            pi =  box * x[ii, :]
            pj =  box * x[jj, :]
            #print pi, pj
            dij = np.linalg.norm (PBC (pj-pi, box))
            #print ii,jj,dij
            if dij < float(rc) :
                #print mol1_id,mol2_id,dij
                bonded = True 
                return bonded         
    return bonded       


def plot_hist(x,y,z):
    #the histogram of the data
    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75, label='$\phi$')
    n, bins, patches = plt.hist(y, 70, normed=1, facecolor='red', alpha=0.5, label='$\\theta$')
    #n, bins, patches = plt.hist(z, 30, normed=1, facecolor='blue', alpha=0.5, label='$\\theta_2$')
    plt.xlabel('angle')
    plt.ylabel('probability')
    plt.legend(loc='upper right')
    #plt.axis([40, 160, 0, 0.03])
    #plt.grid(True)
    #plt.show() 
    with PdfPages('plot1.pdf') as pdf:
         pdf.savefig()
    plt.close()        

def plot_scatter_hist(x,y):
    from matplotlib.ticker import NullFormatter

    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))


    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels for hists
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    axScatter.set_xlabel('$\phi$')
    axScatter.set_ylabel('$\\theta$')

    # the scatter plot:
    axScatter.scatter(x, y)



    # now determine nice limits by hand:
    binwidth = 1
    xmax = np.max(x)
    xmin = np.min(x)
    ymin = np.min(y)
    ymax = np.max(y)

    #lim = (int(xymax/binwidth) + 1) * binwidth

    axScatter.set_xlim( (int(xmin/binwidth) + 1) * binwidth, (int(xmax/binwidth) + 1) * binwidth)
    axScatter.set_ylim( (int(ymin/binwidth) + 1) * binwidth, (int(ymax/binwidth) + 1) * binwidth)

    bins_x = np.arange(xmin, xmax + binwidth, binwidth)
    bins_y = np.arange(ymin, ymax + binwidth, binwidth)
    axHistx.hist(x, bins=bins_x,normed=1)
    axHisty.hist(y, bins=bins_y,normed=1, orientation='horizontal')

    
    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())

    for label in axHistx.get_yticklabels()[::2]:
        label.set_visible(False)
    for label in axHisty.get_xticklabels()[::2]:
        label.set_visible(False)

    with PdfPages('plot2.pdf') as pdf:
         pdf.savefig()
    plt.close() 


def plot_op_vs_time (time, op):
    plt.plot(time,op)
    plt.xlabel('time [MD step]')
    plt.ylabel('fraction of bonds')
    #plt.legend(loc='upper right')
    #plt.axis([40, 160, 0, 0.03])
    #plt.grid(True)
    #plt.show() 
    with PdfPages('plot3.pdf') as pdf:
         pdf.savefig()
    plt.close()        

def plot_scatter_hist_sns(x, y):
    #sns.set(color_codes=True)
    #sns.set(style="darkgrid")
    sns.set(style="ticks")
    sns.jointplot(np.array(x), np.array(y), kind="hex", size=4, stat_func=None).set_axis_labels("$\phi$", "$\\theta$")
    with PdfPages('plot4.pdf') as pdf:
         pdf.savefig()
    sns.plt.close() 
              








if len(sys.argv) == 3 :
   traj_file = sys.argv[1]
   top = sys.argv[2]
   max_timestep = 1e10
elif len(sys.argv) == 4 :
   traj_file = sys.argv[1]
   top = sys.argv[2]
   max_timestep = sys.argv[3]
else:
    print 'Usage: %s dump_file topology [max_timestep]' % (sys.argv[0])
    sys.exit(1)




# atom_phi, atom_theta = read_atomistic_angles('/Users/mm15804/scratch/SAGE/atomistic_trajectory/hub-hub_angles.his')
# plot_scatter_hist_sns(atom_phi, atom_theta)
# sys.exit()


d = read_dump(traj_file, max_timestep)
N_mol = d[-1]['N']/6
#print d
for s, snap in enumerate(d): 
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']

nb  = [0]*N_mol
for mol1_id in range(N_mol):
    for mol2_id in range(N_mol):
        if (mol1_id >= mol2_id): continue
        if if_bonded( mol1_id, mol2_id,x, 0.5, box ):
             #print mol1_id, mol2_id
             nb[mol1_id] += 1
             nb[mol2_id] += 1 
#print
#print nb
#
#print

my_count_dimer = {}
my_count_trimer = {}
for i,val in enumerate(nb):
  if (i%2 == 1):
    my_count_dimer[val] = (my_count_dimer.get(val, 0)) + 1
  else:
    my_count_trimer[val] = (my_count_trimer.get(val, 0)) + 1
print 
print '#trimer_cluster_size, frequency (step=%s)' % (step)
#print '-------------------------------'
for key in sorted(my_count_trimer, key=int):
    print key+1,int(float(my_count_trimer[key])/int(key+1))
print 
print '#dimer_cluster_size, frequency(step=%s)' % (step)
#print '-------------------------------'
for key in sorted(my_count_dimer, key=int):
    print key+1, int(float(my_count_dimer[key])/int(key+1))

print
print

sys.exit()


nb = read_topology(top, d[0]['N'])

# print d[0]['step']
#print nb

print 'analysing...'

op = []
time = []
phi=[]
theta1=[]
theta2=[]
for s, snap in enumerate(d):
 
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    
    # first hub, atoms {1..8}
    cnt = 0
    cnt_bonds = 0
    for line in nb:
        hub1_id = line[0];
        #print line
        #print hub1_id, line[1:]
        for hub2_id in line[1:]:
             if hub1_id > hub2_id: continue
             cnt += 1
             if not if_bonded( hub1_id, hub2_id, 0.75, box ):
                #print '*',step, hub1_id, hub2_id
                continue
             cnt_bonds += 1   
             x1,x2,x3,x4 = get_points(hub1_id, hub2_id, x, p_type, box) 
             dihedral, bend1, bend2 = get_angles(x1,x2,x3,x4,box)
             phi.append(dihedral)
             theta1.append(bend1)
             theta2.append(bend2)

             #print bend1, bend2

    #print step, cnt_bonds, cnt
    op.append(cnt_bonds/(cnt+0.0))
    time.append(step)

plot_hist(phi,theta1, theta2) 
plot_scatter_hist(phi, theta1)
plot_scatter_hist_sns(phi, theta1)
plot_op_vs_time (time, op)







