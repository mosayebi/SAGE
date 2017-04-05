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


def PBC_wrap(x, box):
    while (x <  -box * 0.5):
     x = x + box
    while (x >=  box * 0.5):
     x = x - box
    return x 


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
                      x     [ int(line[0])-1 ] [0] = PBC_wrap(float(line[2]), box_x) #*float(box_x) - float(box_x)/2
                      x     [ int(line[0])-1 ] [1] = PBC_wrap(float(line[3]), box_y)#*float(box_y) - float(box_y)/2
                      x     [ int(line[0])-1 ] [2] = PBC_wrap(float(line[4]), box_z)#*float(box_z) - float(box_z)/2
                snapshot['p_type'] = p_type  
                snapshot['coords'] = x
                data.append(snapshot.copy())
                #print snapshot
                snapshot = {}
   print 'From %s, last TIMESTEP %d' % (filer, step) 
   return data                   

def read_dump_cluster(filer, max_step):
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
                p_type  = np.zeros(N, dtype=np.int)
                cluster = np.zeros(N, dtype=np.int)
                cluster2 = np.zeros(N, dtype=np.int)
                x       = np.zeros((N,3))
                for i in range(N):
                      line = f.readline().split(' ')
                      p_type [ int(line[0])-1 ]     = int  (line[1])
                      x      [ int(line[0])-1 ] [0] = PBC_wrap(float(line[2]), box_x)#*float(box_x) - float(box_x)/2
                      x      [ int(line[0])-1 ] [1] = PBC_wrap(float(line[3]), box_y)#*float(box_y) - float(box_y)/2
                      x      [ int(line[0])-1 ] [2] = PBC_wrap(float(line[4]), box_z)#*float(box_z) - float(box_z)/2
                      cluster[ int(line[0])-1 ]     = int (line[5])

                c = 1
                cluster2[0] = c
                for i in range(1, N):
                    if (cluster[i] == cluster[i-1]):
                	cluster2[i] = c
                    else:
                        c +=1

                print cluster2 
                print max(cluster[:])
                print max(cluster2[:])
                snapshot['p_type'] = p_type  
                snapshot['coords'] = x
                snapshot['cluster']=cluster
                data.append(snapshot.copy())
                #print snapshot
                snapshot = {}
   print 'From %s, last TIMESTEP %d' % (filer, step) 
   return data                

def PBC(d, box):
   d[0] = d[0] - np.rint(d[0]/box[0])*box[0]
   d[1] = d[1] - np.rint(d[1]/box[1])*box[1]
   d[2] = d[2] - np.rint(d[2]/box[2])*box[2]
   return d





def if_bonded(mol1_id, mol2_id, x, rc, box):
    bonded = False
    for p in range(2):
        for i in range(3):
                if (mol1_id%2 + mol2_id%2 == 0):
                    if (p==0):
                      j = i+3
                      ii = int(mol1_id/2)*16 + i + 3
                      jj = int(mol2_id/2)*16 + j + 3
                    else:
                      j = i-3
                      ii = int(mol1_id/2)*16 + i + 6
                      jj = int(mol2_id/2)*16 + j + 6                       
                elif (mol1_id%2 + mol2_id%2 == 2):
                    if (p>0): continue
                    j = i
                    ii = (int(mol1_id/2))*16 + 10 +i + 3
                    jj = (int(mol2_id/2))*16 + 10 +j + 3                
                else:
                    return bonded
                
                print mol1_id,mol2_id,ii,jj,p
                pi =   x[ii, :]
                pj =   x[jj, :]
                #print pi, pj
                dij = np.linalg.norm (PBC (pj-pi, box))
                print mol1_id,mol2_id,ii,jj,p,dij
                if dij < float(rc) :
                    #print mol1_id,mol2_id,dij
                    bonded = True 
                    return bonded         
    return bonded       


def vmd_bond(x1, x2, box, c=0):    
	n = PBC(x2 - x1,box)
	#mag_n = np.linalg.norm(n)
	r1 = x1 
	r2 = r1 + n
        ret=''
	if (not c==0):
            ret += "graphics 0 color %s\n" % c
	ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 30 filled yes\n" %(r1[0], r1[1], r1[2],  r2[0], r2[1], r2[2])
	return ret

def vmd_hub(x1, x2, box, c=0):
	n = PBC(x2 - x1,box)
	mag_n = np.linalg.norm(n)
	r1 = x1 - 0.5*n/mag_n
	r2 = r1 + 3*n/mag_n
	ret=''
	if (not c==0):
            ret += "graphics 0 color %s\n" % c  
	ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.4 resolution 30 filled yes\n" %(r1[0], r1[1], r1[2],  r2[0], r2[1], r2[2])
	return ret

def conf2tcl(snap):
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N = snap['N']

    ret = ""
    # ret = "color Display Background white\n"
    # ret += "mol new\n"
    ret += "graphics 0 delete all\n"

    box_radius = 0.1
    ret += "graphics 0 color 0\n"
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., -box[1]/2., -box[2]/2., box[0]/2., -box[1]/2., -box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., -box[1]/2., box[2]/2., box[0]/2., -box[1]/2., box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., +box[1]/2., -box[2]/2., -box[0]/2., -box[1]/2., -box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., +box[1]/2., box[2]/2., -box[0]/2., -box[1]/2., box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (box[0]/2., +box[1]/2., box[2]/2., box[0]/2., -box[1]/2., box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (box[0]/2., +box[1]/2., -box[2]/2., box[0]/2., -box[1]/2., -box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., +box[1]/2., -box[2]/2., box[0]/2., +box[1]/2., -box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., +box[1]/2., box[2]/2., box[0]/2., +box[1]/2., box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., -box[1]/2.,-box[2]/2.,-box[0]/2., -box[1]/2.,+box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (box[0]/2., -box[1]/2.,-box[2]/2.,box[0]/2., -box[1]/2.,+box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (box[0]/2., box[1]/2.,-box[2]/2.,box[0]/2., box[1]/2.,+box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., box[1]/2.,-box[2]/2.,-box[0]/2., box[1]/2.,+box[2]/2.)


    #N=96
    ret += "graphics 0 color 7\n"   # green
    for i in range(N):
    	if (i%16 ==0 and p_type[i]==1):	
    	   ret += vmd_hub (x[i,:], x[i+2,:], box)
    	else:
    	   continue



    ret += "graphics 0 color 0\n"   # blue
    for i in range(N):
    	#ret += '#%s  %s %s %s %s\n'% (i,p_type[i], x[i,0],x[i,1],x[i,2])  
    	if (i%16 == 10 and p_type[i+3]==8):	
    	   ret += vmd_hub (x[i,:], x[i+2,:], box)
    	else:
    	   continue
    	 
     	   
    ret += "graphics 0 color 1\n"   # red
    for i in range(N):
    	if (i%16 == 10 and p_type[i+3]==11):	
    	   ret += vmd_hub (x[i,:], x[i+2,:], box)
    	else:
    	   continue
    
    ret += "graphics 0 color 6\n"   # silver
    for i in range(N):
    	if (i%16==10 and p_type[i]==1):	
    	   ret += vmd_bond (x[i,:], x[i-10,:], box)
    	else:
    	   continue

    # print x[93,:]
    # print p_type[93]
    return ret

def conf2tcl_cluster(snap, cluster_flag=True):
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N = snap['N']
    cluster=snap['cluster']
    print cluster 

    ret = ""
    # ret = "color Display Background white\n"
    # ret += "mol new\n"
    ret += "graphics 0 delete all\n"

    box_radius = 0.1
    ret += "graphics 0 color 0\n"
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., -box[1]/2., -box[2]/2., box[0]/2., -box[1]/2., -box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., -box[1]/2., box[2]/2., box[0]/2., -box[1]/2., box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., +box[1]/2., -box[2]/2., -box[0]/2., -box[1]/2., -box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., +box[1]/2., box[2]/2., -box[0]/2., -box[1]/2., box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (box[0]/2., +box[1]/2., box[2]/2., box[0]/2., -box[1]/2., box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (box[0]/2., +box[1]/2., -box[2]/2., box[0]/2., -box[1]/2., -box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., +box[1]/2., -box[2]/2., box[0]/2., +box[1]/2., -box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., +box[1]/2., box[2]/2., box[0]/2., +box[1]/2., box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., -box[1]/2.,-box[2]/2.,-box[0]/2., -box[1]/2.,+box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (box[0]/2., -box[1]/2.,-box[2]/2.,box[0]/2., -box[1]/2.,+box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (box[0]/2., box[1]/2.,-box[2]/2.,box[0]/2., box[1]/2.,+box[2]/2.)
    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-box[0]/2., box[1]/2.,-box[2]/2.,-box[0]/2., box[1]/2.,+box[2]/2.)


    #N=96
    ret += "graphics 0 color 7\n"   # blue
    for i in range(N):
    	if (i%16 ==0 and p_type[i]==1):	
    	   if (cluster_flag): 	
    	     ret += vmd_hub (x[i,:], x[i+2,:], box, int(float(cluster[i])/max(cluster[:]) * 1023 +1) )
    	   else:
    	     ret += vmd_hub (x[i,:], x[i+2,:], box, 0 )    
    	else:
    	   continue



    ret += "graphics 0 color 0\n"   # blue
    for i in range(N):
    	#ret += '#%s  %s %s %s %s\n'% (i,p_type[i], x[i,0],x[i,1],x[i,2])  
    	if (i%16 == 10 and p_type[i+3]==8):	
    	   if (cluster_flag): 	
    	     ret += vmd_hub (x[i,:], x[i+2,:], box, int(float(cluster[i])/max(cluster[:]) * 1023 +1) )
    	   else:
    	     ret += vmd_hub (x[i,:], x[i+2,:], box, 0 ) 
    	else:
    	   continue
    	 
     	   
    ret += "graphics 0 color 1\n"   # red
    for i in range(N):
    	if (i%16 == 10 and p_type[i+3]==11):	
    	   if (cluster_flag): 	
    	     ret += vmd_hub (x[i,:], x[i+2,:], box, int(float(cluster[i])/max(cluster[:]) * 1023 +1) )
    	   else:
    	     ret += vmd_hub (x[i,:], x[i+2,:], box, 0 ) 
    	else:
    	   continue
    
    ret += "graphics 0 color 6\n"   # silver
    for i in range(N):
    	if (i%16==10 and p_type[i]==1):	
            if (cluster_flag):
      	     ret += vmd_bond (x[i,:], x[i-10,:], box, int(float(cluster[i])/max(cluster[:])*1023 + 1) )
      	    else:
      	     ret += vmd_bond (x[i,:], x[i-10,:], box, 0 )     	        	
    	else:
    	   continue

    # print x[93,:]
    # print p_type[93]
    return ret









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



d = read_dump_cluster(traj_file, max_timestep)
N_mol = d[-1]['N']/16
tcl_out  = "color Display Background white\n"
tcl_out  +="display projection orthographic\n"
tcl_out  +="axes location off\n"
tcl_out  += "mol new\n"
tcl_out  += "color scale method BGR\n"


# tcl_out  += "color scale max 0.9\n"
# tcl_out  += "color scale midpoint 0.3\n"
# for s, snap in enumerate(d): 
#     tcl_out += conf2tcl(snap)
#     tcl_out += "animate dup 0\n"
#     tcl_out += "animate next\n"

tcl_out += conf2tcl_cluster(d[-1],False)
tcl_out  += "display resetview\n"
tcl_out  += "rotate x by 45\n"
tcl_out  += "rotate z by 5\n"

outname = 'out.tcl'
try:
    f = open (outname, 'w')
except:
    print ("tcl_output: cannot open output file. Dying")
f.write (tcl_out)

sys.exit()





nb  = [0]*N_mol
for mol1_id in range(N_mol):
    for mol2_id in range(N_mol):
        if (mol1_id >= mol2_id): continue
        if if_bonded( mol1_id, mol2_id,x, 0.5, box ):
             #print mol1_id, mol2_id
             nb[mol1_id] += 1
             nb[mol2_id] += 1 
