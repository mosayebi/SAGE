import os
import sys
#import fnmatch
import numpy as np
import subprocess
import time
import errno 


def read_atomistic_angles(filer):
  phi=[]
  theta1=[]
  theta2=[]
  with open(filer,'r') as f:
      while True:
        line=f.readline().strip('\n')
        line=f.readline().strip('\n').split()
        while line:           
          if line[0] == 'Time' :
            #print line[2], len(phi)
            if float(line[2])>=500: return np.array(phi), np.array(theta1), np.array(theta2)
            line=f.readline().strip('\n').split()
            #continue
          phi.append(float(line[3]))
          theta1.append(float( line[4]) - 90. )
          theta2.append(float( line[5]) - 90. )
          line=f.readline().strip('\n').split()
      #print len(phi)    
      return np.array(phi), np.array(theta1), np.array(theta2)   



def read_dump(filer, min_step=0, max_step=1e10):
   data=[]  
   snapshot = {}
   snapshot['traj'] = filer
   cluster_flag = False
   minstep_flag = False
   read_flag = False
   if (max_step<min_step) : min_step, max_step = max_step, min_step
   cnt = 0
   with open(filer,'r') as f:
      while True:
        line=f.readline().strip('\n')
        if not line: break
        items = line.split(' ')
        if items[0] == 'ITEM:':
            if items[1] == 'TIMESTEP':
                step = int(f.readline().split(' ')[0])
                if step > max_step :
                    print 'max_step reached (%d)' % max_step 
                    print 'From %s, last TIMESTEP %d' % (filer, data[-1]['step'])
                    return data
                if ( step >= min_step and minstep_flag == False) :
                    print 'From %s, first TIMESTEP reached (%d)' % (filer, min_step)
                    minstep_flag = True
                if (step >= min_step and step <= max_step): 
                    read_flag = True
                    cnt += 1                      
                    print '%d : reading TIMESTEP %d' % (cnt,step)        
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
                if (len(items)>7):
                   if ( items[7] == "c_c1"): 
                      cluster_flag=True
                      cluster = np.zeros(N, dtype=np.int)
                   else:
                      cluster = []   
                p_type = np.zeros(N, dtype=np.int)
                x = np.zeros((N,3))
                for i in range(N):
                      line = f.readline().split(' ')
                      p_type[ int(line[0])-1 ]     = int  (line[1])
                      x     [ int(line[0])-1 ] [0] = PBC_wrap(float(line[2]), box_x) #*float(box_x) - float(box_x)/2
                      x     [ int(line[0])-1 ] [1] = PBC_wrap(float(line[3]), box_y)#*float(box_y) - float(box_y)/2
                      x     [ int(line[0])-1 ] [2] = PBC_wrap(float(line[4]), box_z)#*float(box_z) - float(box_z)/2
                      if (cluster_flag): cluster[ int(line[0])-1 ]     = int (line[5])

                snapshot['p_type'] = p_type  
                snapshot['coords'] = x
                snapshot['cluster']= cluster
                
                if (read_flag): data.append(snapshot.copy())

                #print snapshot
                snapshot = {}
                snapshot['traj'] = filer
   print 'From %s, last TIMESTEP %d' % (filer, data[-1]['step'])
   return data                   


def conf2tcl(snap, cluster_flag=True, box_flag=True, psi3_flag=False, psi3=[]):
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N = snap['N']
    cluster=snap['cluster']
    if (cluster_flag and len(cluster)==0): 
        print "Warning: no cluster data found. Cluster coloring is off."
        cluster_flag=False
        psi3_flag = False
    if (len(psi3)!=0):
        psi3_flag = True 
        cluster_flag = False 
        psi3 = np.nan_to_num(psi3)

    #print cluster 

    

    ret  = "color Display Background white\n"
    ret  +="display projection orthographic\n"
    ret  +="axes location off\n"
    ret  += "mol new\n"
    ret  += "color scale method RGB\n"
    ret  += "graphics 0 delete all\n"

    if (box_flag):
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



    if (psi3_flag): 
        cmax =  1.
        cmin =  0.
        # if len(psi3[np.nonzero(np.nan_to_num(psi3))]) > 0 :
        #    cmin = np.min(psi3[np.nonzero(psi3)])
        #    if cmin == 1 : 
        #     cmin = 0
        # print cmin, cmax   

    #N=96
    ret += "graphics 0 color 7\n"   # blue
    for i in range(N):
        if (i%16 ==0 and p_type[i]==1): 
           if (cluster_flag):   
              ret += vmd_hub (x[i,:], x[i+2,:], box, int(float(cluster[i])/max(cluster[:]) * 1023 +1) )
           elif (psi3_flag):
              c = np.nan_to_num(psi3[get_mol_id(i)]) ** 2
              if c == 0:
                c=1
              else:
                c=int((c-cmin)/(cmax-cmin) * 950 + 70)
              ret += vmd_hub (x[i,:], x[i+2,:], box, c )             
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
           elif (psi3_flag):
              c = np.nan_to_num(psi3[get_mol_id(i)]) ** 2
              if c == 0:
                c=1
              else:
                c=int((c-cmin)/(cmax-cmin) * 950 + 70)
              ret += vmd_hub (x[i,:], x[i+2,:], box, c ) 
           else:
              ret += vmd_hub (x[i,:], x[i+2,:], box, 0 ) 
        else:
           continue
                 
    ret += "graphics 0 color 1\n"   # red
    for i in range(N):
        if (i%16 == 10 and p_type[i+3]==11):    
           if (cluster_flag):   
              ret += vmd_hub (x[i,:], x[i+2,:], box, int(float(cluster[i])/max(cluster[:]) * 1023 +1) )
           elif (psi3_flag):
              c = np.nan_to_num(psi3[get_mol_id(i)]) ** 2
              if c == 0:
                c=1
              else:
                c=int((c-cmin)/(cmax-cmin) * 950 + 70)
              ret += vmd_hub (x[i,:], x[i+2,:], box, c )  
           else:
              ret += vmd_hub (x[i,:], x[i+2,:], box, 0 ) 
        else:
           continue
    
    ret += "graphics 0 color 6\n"   # silver
    for i in range(N):
        if (i%16==10 and p_type[i]==1): 
            if (cluster_flag):
              ret += vmd_bond (x[i,:], x[i-10,:], box, int(float(cluster[i])/max(cluster[:])*1023 + 1) )
            elif (psi3_flag):
              c = np.nan_to_num(psi3[get_mol_id(i-10)]) ** 2
              if c == 0:
                c=1
              else:
                c=int((c-cmin)/(cmax-cmin) * 950 + 70)               
              ret += vmd_bond (x[i,:], x[i-10,:], box, c )              
            else:
              ret += vmd_bond (x[i,:], x[i-10,:], box, 0 )                   
        else:
           continue

    ret += "display resetview\n"
    ret += "rotate x by 10\n"
    ret += "rotate y by -10\n"
    return ret


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

def write_tcl_out (tcl_out, filename="out.tcl"):
    try:
        f = open (filename, 'w')
    except:
        print ("tcl_output: cannot open output file. Dying")
    f.write (tcl_out)

def render_tcl_file (res_x, res_y, png_file="out.png", tcl_out_file="out.tcl"):
    tcl_out = "source %s\nrender Tachyon out.render\n" % (tcl_out_file)
    vmdin = os.popen("/Applications/VMD\ 1.9.2.app/Contents/Resources/VMD.app/Contents/MacOS/VMD -dispdev none","w")
    vmdin.write("%s"%tcl_out)
    vmdin.flush()

    args = "-aasamples 12 out.render -format TGA -res %s %s -o out.tga"%(res_x,res_y)
    command=[]
    command.append("/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86")
    command = command + args.split()
    exe = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = exe.communicate()  
    print stdout, stderr

    command = "convert out.tga -transparent black %s"% png_file
    exe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = exe.communicate()  
    print stdout, stderr
   
    # command = "rm %s out.tga out.render"%tcl_out_file
    # exe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout, stderr = exe.communicate()  
    # print stdout, stderr

def traj2movie(traj_data, movie_file='movie', cluster_flag=True, rotation_flag=True):
    image_dir = 'images'
    make_sure_path_exists(image_dir)
    for s, snap in enumerate(traj_data):
      tcl_out = conf2tcl (snap, cluster_flag)
      write_tcl_out (tcl_out, filename="out.tcl")
      png_file = image_dir+movie_file+".%05d.png"%int(s)
      print png_file
      render_tcl_file (1024, 1024, png_file, "out.tcl")
    #s = 415
    #tcl_out = conf2tcl (traj_data[-1], cluster_flag)
    if (rotation_flag):
       N_rot = 60
       for r in range(N_rot+1):
          rot_x = 10 + r * 360./N_rot
          rot_y = 0 #-10 + r*360/50
          write_tcl_out (tcl_out + "rotate x by %s\n rotate y by %s\n"%(rot_x,rot_y), filename="out.tcl")
          png_file = image_dir+movie_file+".%05d.png"%int(s)
          print png_file
          render_tcl_file (1024, 1024, png_file, "out.tcl")      
          s += 1
       for r in range(N_rot+1):
          rot_x = 0 #10 + r*360/50
          rot_y = -10 - r*360./N_rot
          write_tcl_out (tcl_out + "rotate x by %s\n rotate y by %s\n"%(rot_x,rot_y), filename="out.tcl")
          png_file = image_dir+movie_file+".%05d.png"%int(s)
          print png_file
          render_tcl_file (1024, 1024, png_file, "out.tcl")      
          s += 1
    command = 'ffmpeg -y -r 10 -i '+ image_dir+movie_file +r'.%05d.png'+' -c:v libx264 -r 30 -pix_fmt yuv420p %s.mp4'%movie_file         
    #command = ['ffmpeg', '-t', '60', '-i', r'out_%05d.png',  '-y', movie_file]
    exe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = exe.communicate()  
    print stdout, stderr
    print 
    print "trajectory movie saved to %s.mp4"%movie_file
    print 
    print 'ffmpeg command used: (-r 10, means 1/10 s per frame for input images. so you can set the movie duration.).' 
    print command
    print 'you may as well want to consider morphing'
    print r'convert *.png -delay 2 -morph 2 morph_%05d.png'



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

def PBC_wrap(x, box):
    while (x <  -box * 0.5):
     x = x + box
    while (x >=  box * 0.5):
     x = x - box
    return x 

def get_angles(x1,x2,x3,x4,box):
    #    
    #  2 3
    #  1 4
    #
    b1 = (PBC(x2 - x1, box)) / np.linalg.norm (PBC(x2 - x1, box))
    b2 = (PBC(x3 - x2, box)) / np.linalg.norm (PBC(x3 - x2, box))
    b3 = (PBC(x4 - x3, box)) / np.linalg.norm (PBC(x4 - x3, box)) #

    b12 = np.cross(b1,b2)
    b23 = np.cross(b2,b3)
             
    dihedral = np.degrees ( np.arctan2 ( np.dot( np.cross(b12, b23), b2) ,    np.dot(b12,b23) )  )
    bend1 = np.degrees( np.arccos ( round(np.dot(b1,b2),12))) - 90 
    bend2 = np.degrees( np.arccos ( round(np.dot(b2,b3),12))) - 90
    #bend2 = np.degrees( np.arccos ( round(np.dot(b1,-b3),12)))

    return dihedral, bend1, bend2



def get_helix_atom_ids(mol_id):
  atom_id = []
  atom1 = int(mol_id/2)*16
  if (mol_id%2 == 0):
      atom_id = range(atom1, atom1+3)
  else:
      atom_id = range(atom1+10, atom1+13) 
  return atom_id

def get_patch_atom_ids(mol_id):
  atom_id = []
  atom1 = int(mol_id/2)*16
  if (mol_id%2 == 0):
      atom_id = [range(atom1+3, atom1+6), range(atom1+6, atom1+9)]
      #atom_id = [atom1+3, atom1+6, atom1+4, atom1+7, atom1+5, atom1+8]
  else:
      atom_id = [range(atom1+13, atom1+16)] 
  return atom_id

def get_mol_id(atom_id):
    mol_id = int(atom_id / 16) * 2
    if (atom_id % 16) >= 10 :
        mol_id += 1
    return mol_id    

def get_points(coords1, coords2):
    x1 = coords1[0,:]
    x2 = coords1[1,:]
    x3 = coords2[1,:]
    x4 = coords2[0,:]
    return x1,x2,x3,x4





def build_nb_list(snap):
    rc = 0.6
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N_mol = snap['N']/16 * 2
    print '#building neighbour list.'
    #print 'N_mol= %s'%N_mol

    nb_no = np.zeros(N_mol, dtype=np.int)
    nb_list = - np.ones((N_mol,5), dtype=np.int)
    for mol1_id in range(N_mol):
      for mol2_id in range(mol1_id+1, N_mol):
        if ( (mol1_id%2 + mol2_id%2)==1 ): continue
        #print mol1_id, mol2_id
        ids1 = get_patch_atom_ids(mol1_id)
        ids2 = get_patch_atom_ids(mol2_id)
        if if_bonded(x[ids1,:], x[ids2,:], box, rc ):
             #print mol1_id, mol2_id
             nb_no[mol1_id] += 1
             nb_no[mol2_id] += 1
             for i, val in enumerate(nb_list[mol1_id,:]):
                 if (val == -1): break
             nb_list[mol1_id, i] = mol2_id 
             for i, val in enumerate(nb_list[mol2_id,:]):
                 if (val == -1): break
             nb_list[mol2_id, i] = mol1_id
    return nb_no, nb_list 

def if_bonded(xis, xjs, box, rc):
    bonded = False
    no_patch_groups = len(xis)
    for i in range(no_patch_groups):
      pi_set = xis[i]
      pj_set = xjs[(i+1)%no_patch_groups]
      #print '*'
      for ip in range(3):
       #print '+'
          dij = np.linalg.norm (PBC (pj_set[ip]-pi_set[ip] , box)) 
          #print dij, pi_set[ip], pj_set[ip]
          if dij < float(rc) :
             #print mol1_id,mol2_id,dij
             bonded = True 
             return bonded         
    return bonded       


def generate_HK_data(no_nb, nb_list):
    N_mol=len(no_nb)
    maxnb = 5
    linkid = 0
    NodeNext = np.zeros(N_mol, dtype=np.int)
    LinksOfNode = np.zeros((N_mol, maxnb), dtype=np.int)
    Nodes = np.ones(N_mol, dtype=np.int)   # change this if you want to define occupancy condition for Nodes

    maxind = 0

    for i in range(N_mol):
        for j in range(i+1, N_mol):
            if (j==i+1 or i in nb_list[0:no_nb[j]] and  NodeS[i]*NodeS[j]==1):
                linkid += 1
                ind = 0
                while (LinksOfNode[i, ind] != 0):
                    ind +=1
                if (maxind < ind ): maxind = ind
                LinksOfNode[i,ind] = linkid    
                NodeNext[i,ind] = j
                ind = 0
                while (LinksOfNode[j, ind] != 0):
                    ind += 1
                LinksOfNode[j,ind] = linkid    
                NodeNext[j,ind] = i
    return NodeS, NodeNext, LinksOfNode
    

    


def build_hub_hub_pairs(nb_list, nb_no):
    N_mol = len(nb_no)
    max_N_hub = N_mol/2
  
    hub_hub_pairs=[]

    flag = np.zeros(N_mol, dtype=np.int)
    for mol_i in range(0,N_mol,2):
        hub_arms = []
        if ( nb_no[ mol_i+1 ] > 1):
            print '#Warning mol %s has more than one arm (#nb = %s)' % (mol_i,  nb_no[ mol_i+1])
            #sys.exit(1)
        if ( flag[mol_i] == 0):
           for k in range(nb_no[ mol_i+1 ]):  
               paired_to_mol_i = nb_list[ mol_i+1, k] - 1  # find a trimer mol. id that is attached to the trimer mol_i
               hub_arms.append([mol_i, paired_to_mol_i])
               flag[mol_i] = 1    
        # loop over other members of the trimer       
        for j in range(nb_no[mol_i]):
            mol_j = nb_list[mol_i, j]      
            if ( nb_no[ mol_j+1 ] > 1):
                print '#Warning mol %s has more than one arm (#nb = %s)' % (mol_j,  nb_no[ mol_j+1])
                #sys.exit(1)
            if ( flag[mol_j] == 0):
               for k in range(nb_no[ mol_j+1 ]): 
                  paired_to_mol_j = nb_list[ mol_j+1, k] - 1
                  hub_arms.append([mol_j, paired_to_mol_j])
                  flag[mol_j] = 1
        if len(hub_arms)>0 : hub_hub_pairs.append(hub_arms)
    return hub_hub_pairs
    

def orientational_order(hub_hub_pairs, snap):
    x = snap['coords']
    N_mol = snap['N']/16*2
    box = snap['box']
    psi3 = 0.0
    n = 0
    my_count = {}
    #print N_mol
    psi3_vec = np.zeros(N_mol)
    psi3_vec_re = np.zeros(N_mol)
    psi3_vec_im = np.zeros(N_mol)
    cnt = np.zeros(N_mol, dtype=np.int)
    angles = []
    psi3_re = 0.
    psi3_im = 0.

    # print hub_hub_pairs
    # print

    for i, pairs in enumerate(hub_hub_pairs):
        for j in range(len(pairs)):
            for k in range(j+1,len(pairs)):
                pairj = pairs[j]
                pairk = pairs[k]
                vec_j = PBC( x[ get_helix_atom_ids(pairj[1])[1],:] - x[ get_helix_atom_ids(pairj[0])[1],:], box)
                vec_k = PBC( x[ get_helix_atom_ids(pairk[1])[1],:] - x[ get_helix_atom_ids(pairk[0])[1],:], box)
                vec_j = vec_j / np.linalg.norm(vec_j)
                vec_k = vec_k / np.linalg.norm(vec_k)
                cos =  np.dot (vec_j, vec_k)
                sin =  np.linalg.norm (np.cross (vec_j, vec_k))
                angles.append(np.degrees(np.arctan2(sin,cos)))
                psi3_rei = 4*(cos**3) - 3*cos
                psi3_imi = 3*sin - 4*(sin**3)
                 
                psi3_re +=  psi3_rei
                psi3_im +=  psi3_imi
                # print j, k, cos, 4*(cos**3) - 3*cos
                # print '****'
                # #sys.exit()
                n += 1
                val = len(pairs)
                my_count[val] = (my_count.get(val, 0)) + 1 

                mol_list = [pairk[0], pairk[1], pairj[0], pairj[1], pairk[0]+1, pairk[1]+1, pairj[0]+1, pairj[1]+1]
                psi3_vec_re [mol_list] += psi3_rei
                psi3_vec_im [mol_list] += psi3_imi

                cnt [mol_list] += 1
                # if (psi3_i)<0: 
                #     print 'Error. psi3 cannot be negative. orientational_order()'
                #     sys.exit()
                # val = len(pairs[j])
                # my_count[val] = (my_count.get(val, 0)) + 1

    if (n>0): 
        psi3_im /= n
        psi3_re /= n
        psi3 = np.sqrt(psi3_im**2 + psi3_re**2) 
    # n is the number of angles and psi3 the hony_comb order parameter
    for i in range(len(cnt)):
        if (cnt[i]>0): 
            psi3_vec_re[i] /= cnt[i]
            psi3_vec_im[i] /= cnt[i]
            psi3_vec[i] = np.sqrt(psi3_vec_re[i]**2 + psi3_vec_im[i]**2)
        else:
            psi3_vec[i] = 0.0
    return psi3, n, my_count, psi3_vec, angles



def traj2psi3(traj_data, filename='psi3_file', skip_snap=True, tcl_write_flag=True, psi3_angle_hist_flag=True):
    last_timestep = 0  
    # get the last snapshot already written to filename
    if (not skip_snap): 
        try:
            f = open (filename, 'w')
        except:
            print ("get_psi3(): cannot open output file. Dying")
            sys.exit(1)
    else:        
        try:
            f = open (filename, 'a+')
        except:
            print ("get_psi3(): cannot open output file. Dying")
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
                print '\n#Warning: last line of %s is not as expected. Setting skip_snap to False.\n\n'% filename

    

    psi3_traj=[]
    psi3_snap={}
    out = "#"+time.strftime("%c")+"\n" 
    out+= '#step psi3 N_angles count(arms>=2)\n'
    f.write(out)
    print '#step psi3 N_angles count(arms>=2)'
    #print traj_data
    for snap in traj_data:
        #print snap
        #print
        psi3_snap = {}
        step = snap['step']
        psi3_snap['step'] = step
        if (skip_snap and step <= last_timestep): 
            print 'skipping %s' %step
            continue
        nb_no, nb_list = build_nb_list(snap)
        hub_hub_pairs = build_hub_hub_pairs(nb_list, nb_no)
        psi3, N_angles, my_count, psi3_vec, angles = orientational_order(hub_hub_pairs, snap)
        out='%s %s %s %s\n' % (step, psi3, N_angles, my_count)
        print '%s %s %s %s' % (step, psi3, N_angles, my_count)
        # for i, val in enumerate(psi3_vec):
        #     print i, val

        if (tcl_write_flag):
            make_sure_path_exists('tcl')
            tcl_file = 'tcl/'+traj_file+'.psi3.'+str(snap['step'])+".tcl"
            tcl_out = conf2tcl (snap, cluster_flag=False, box_flag=True, psi3_flag = True, psi3 = psi3_vec) 
            write_tcl_out (tcl_out, filename=tcl_file)
            print '# wrote tcl output with psi3 coloring to %s'%tcl_file

        if (psi3_angle_hist_flag):
            make_sure_path_exists('psi3_angle_hist')
            pdf_file = 'psi3_angle_hist/'+traj_file+'.hist.'+str(snap['step'])+".pdf"
            plot_1D_hist (angles, pdf_file) 


        psi3_snap['psi3']=psi3_vec
        psi3_traj.append(psi3_snap.copy())
        f.write (out)
    print '#psi3 output is appended to %s' % filename     
    return psi3_traj



def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


# def if_bonded(mol1_id, mol2_id, x, rc, box):
#     bonded = False
#     for p in range(2):
#         for i in range(3):
#                 if (mol1_id%2 + mol2_id%2 == 0):
#                     if (p==0):
#                       j = i+3
#                       ii = int(mol1_id/2)*16 + i + 3
#                       jj = int(mol2_id/2)*16 + j + 3
#                     else:
#                       j = i-3
#                       ii = int(mol1_id/2)*16 + i + 6
#                       jj = int(mol2_id/2)*16 + j + 6                       
#                 elif (mol1_id%2 + mol2_id%2 == 2):
#                     if (p>0): continue
#                     j = i
#                     ii = (int(mol1_id/2))*16 + 10 +i + 3
#                     jj = (int(mol2_id/2))*16 + 10 +j + 3                
#                 else:
#                     return bonded
                
#                 print mol1_id,mol2_id,ii,jj,p
#                 pi =  x[ii, :]
#                 pj =  x[jj, :]
#                 #print pi, pj
#                 dij = np.linalg.norm (PBC (pj-pi, box))
#                 print mol1_id,mol2_id,ii,jj,p,dij
#                 if dij < float(rc) :
#                     #print mol1_id,mol2_id,dij
#                     bonded = True 
#                     return bonded         
#     return bonded       


def plot_hist(x,y,z, filename='hist.pdf'):
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.optimize import curve_fit
    from scipy import stats, integrate
    #print x,y,z
    (y,z) = (z,y)
    #the histogram of the data
    nx, binsx, patchesx = plt.hist(x, 72, normed=1, facecolor='green', alpha=0.75, label='$\phi$')
    ny, binsy, patchesy = plt.hist(y, 71, normed=1, facecolor='red', alpha=0.6, label=r'$\theta_1$')
    nz, binsz, patchesz = plt.hist(z, 70, normed=1, facecolor='cyan', alpha=0.3, label=r'$\theta_2$')

    (mu_z, sigma_z) = stats.norm.fit(z)
    (mu_y, sigma_y) = stats.norm.fit(y)
    (mu_x, sigma_x) = stats.norm.fit(x)

    fitx = mlab.normpdf( binsx, mu_x, sigma_x)
    fity = mlab.normpdf( binsy, mu_y, sigma_y)
    fitz = mlab.normpdf( binsz, mu_z, sigma_z)
    lx = plt.plot(binsx, fitx, 'g-', linewidth=1, alpha=0.85)
    ly = plt.plot(binsy, fity, 'r-', linewidth=1, alpha=0.85)
    lz = plt.plot(binsz, fitz, 'c-', linewidth=1, alpha=0.85)

    #n, bins, patches = plt.hist(z, 30, normed=1, facecolor='blue', alpha=0.5, label='$\\theta_2$')
    plt.xlabel(r'$\mathrm{angle}\ ^\circ$')
    plt.ylabel(r'$\mathrm {probability}$')
    plt.title(r'$\mu_{\phi}=%.2f,\ \sigma_{\phi}=%.2f\  \mu_{\theta_1}=%.2f,\ \sigma_{\theta_1}=%.2f\  \mu_{\theta_2}=%.2f,\ \sigma_{\theta_2}=%.2f$' %(mu_x, sigma_x, mu_y, sigma_y, mu_z, sigma_z))
    plt.legend(loc='upper right')
    plt.xlim([-15,15])
    #plt.grid(True)
    #plt.show() 
    with PdfPages(filename) as pdf:
         pdf.savefig()
    plt.close()   
    print "angel histograms are saved in %s" % filename     


def plot_1D_hist(x, filename='hist.pdf'):
    #import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # from scipy.optimize import curve_fit
    # from scipy import stats, integrate

    #the histogram of the data
    nx, binsx, patchesx = plt.hist(x, 50, normed=len(x), facecolor='green', alpha=0.75, label='$\psi$')

    #(mu_x, sigma_x) = stats.norm.fit(x)

    #fitx = mlab.normpdf( binsx, mu_x, sigma_x)
    #lx = plt.plot(binsx, fitx, 'g-', linewidth=1, alpha=0.85)


    #n, bins, patches = plt.hist(z, 30, normed=1, facecolor='blue', alpha=0.5, label='$\\theta_2$')
    plt.xlabel(r'$\mathrm{angle}\ ^\circ$')
    plt.ylabel(r'$\mathrm {probability}$')
    plt.title(r'$N=%d,\ \psi_{{\mathrm{avg}}}=%.2f$' % (len(x), sum(x)/len(x)))
    #plt.title(r'$\mu_{\psi}=%.2f,\ \sigma_{\phi}=%.2f\  \mu_{\theta_1}=%.2f,\ \sigma_{\theta_1}=%.2f\  \mu_{\theta_2}=%.2f,\ \sigma_{\theta_2}=%.2f$' %(mu_x, sigma_x, mu_y, sigma_y, mu_z, sigma_z))
    plt.legend(loc='upper right')
    #plt.xlim([-15,15])
    #plt.grid(True)
    #plt.show() 
    with PdfPages(filename) as pdf:
         pdf.savefig()
    plt.close()   
    print "histogram saved in %s" % filename   


def plot_1D_hist_noX(x, filename='hist.pdf'):
    import matplotlib as mpl
    mpl.use('pdf')
    import matplotlib.pyplot as plt

    # from scipy.optimize import curve_fit
    # from scipy import stats, integrate
    #the histogram of the data
    nx, binsx, patchesx = plt.hist(x, 60, normed=len(x), facecolor='green', alpha=0.75, label='$\psi$')
    #(mu_x, sigma_x) = stats.norm.fit(x)
    #fitx = mlab.normpdf( binsx, mu_x, sigma_x)
    #lx = plt.plot(binsx, fitx, 'g-', linewidth=1, alpha=0.85)
    #n, bins, patches = plt.hist(z, 30, normed=1, facecolor='blue', alpha=0.5, label='$\\theta_2$')
    plt.xlabel(r'$\mathrm{angle}\ ^\circ$')
    plt.ylabel(r'$\mathrm {probability}$')
    plt.title(r'$N=%d,\ \psi_{{\mathrm{avg}}}=%.2f$' % (len(x), sum(x)/len(x)))
    #plt.title(r'$\mu_{\psi}=%.2f,\ \sigma_{\phi}=%.2f\  \mu_{\theta_1}=%.2f,\ \sigma_{\theta_1}=%.2f\  \mu_{\theta_2}=%.2f,\ \sigma_{\theta_2}=%.2f$' %(mu_x, sigma_x, mu_y, sigma_y, mu_z, sigma_z))
    plt.legend(loc='upper right')
    #plt.xlim([-15,15])
    #plt.grid(True)
    #plt.show() 
    # with PdfPages(filename) as pdf:
    #      pdf.savefig()
    # plt.close()   
    # print "histogram saved in %s" % filename  
    fig = plt.figure()
    fig.savefig(filename)


def plot_scatter_hist(x,y):
    #import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    #from scipy.optimize import curve_fit
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
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
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
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import pandas as pd
    from scipy import stats, integrate
    from scipy.stats import kendalltau
    #mlab.use('Agg')  # Force mpl to use Tk backend
    import seaborn as sns
    #sns.set(color_codes=True)
    #sns.set(style="darkgrid")
    sns.set(style="ticks")
    sns.jointplot(np.array(x), np.array(y), kind="hex", size=4, stat_func=None).set_axis_labels("$\phi$", "$\\theta$")
    with PdfPages('plot4.pdf') as pdf:
         pdf.savefig()
    sns.plt.close() 
              
def traj2angles_plot(traj_data, plot_file = 'angle_hist.pdf' ):
    phi=[]
    theta1=[]
    theta2=[]
    for s, snap in enumerate(traj_data): 
        x = snap['coords']
        p_type = snap['p_type']
        box = snap['box']
        step = snap['step']
        N_mol = snap['N']/16*2
    
        for mol in range(0, N_mol, 2):
            x1,x2,x3,x4 = get_points( x[get_helix_atom_ids(mol),:] ,  x[get_helix_atom_ids(mol+1),:])
            #print x1,x2,x3,x4 
            dihedral, bend1, bend2 = get_angles(x1,x2,x3,x4,box)
            phi.append(dihedral)
            #print bend1,bend2
            theta1.append(bend1)
            theta2.append(bend2)
    filename =    plot_file+'.angles_hist.pdf'         
    plot_hist(phi, theta1, theta2, filename ) 
    #print phi, theta1, theta2
    #plot_scatter_hist(phi, theta1)
    #plot_scatter_hist_sns(phi, theta1) 


def read_lammps_cluster_log(logfile):   
    timestep=[]
    mean_s=[]
    N_cluster=[]
    max_cluster=[]
    min_cluster=[]
    with open(logfile,'r') as f:
      while True:
        line=f.readline().strip('\n').split() 
        while (line[0] != 'Step' and line[1]!='maxxluid'):
            line=f.readline().strip('\n').split()
            while len(line)==0 : line=f.readline().strip('\n').split()
        line=f.readline().strip('\n') 
        while line:
           line = line.split()
           if line[0]=='Loop' : 
              print 'from %s, last timestep read is %s' %(logfile, timestep[-1])
              return timestep, mean_s, N_cluster, min_cluster, max_cluster
           if line[0] != 'WARNING:' :
              timestep.append(float(line[0]))
              mean_s.append(float(line[2]))
              N_cluster.append(float(line[3]))
              min_cluster.append(float(line[4]))
              max_cluster.append(float(line[5]))
           line=f.readline().strip('\n')
      print 'from %s, last timestep read is %d' %(logfile, timestep[-1])
      return timestep, mean_s, N_cluster, min_cluster, max_cluster


def analyse_cluster_log(plotfile, logfile='log_rerun_0.75'):    
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    timestep, mean_s, N_cluster, min_cluster, max_cluster = read_lammps_cluster_log(logfile)


    plt.plot(timestep,mean_s, linewidth=2, alpha=1, label=r'mean cluster size, $s$')
    plt.plot(timestep,max_cluster, linewidth=2, alpha=1, label=r'max cluster size, $s_{max}$')
    plt.plot(timestep,min_cluster, linewidth=2, alpha=1, label=r'min cluster size, $s_{min}$')

    plt.xlabel('time [MD step]')
    plt.ylabel('size [#hub]')  
    plt.legend(loc='center right')
    #plt.axis([40, 160, 0, 0.03])
    #plt.grid(True)
    #plt.show() 
    with PdfPages(plotfile+'_cluster_size.pdf') as pdf:
         pdf.savefig()
    plt.close()
    print 'cluster info plotted in %s' % plotfile+'_cluster_size.pdf'  

           
  
def read_lammps_log(logfile):   
    timestep=[]
    temp=[]
    Epair=[]
    Emol=[]
    with open(logfile,'r') as f:
      while True : 
        line=f.readline()
        if not line:
            print 'from %s, last timestep read is %d' %(logfile, timestep[-1])
            return timestep, temp, Epair, Emol
        line = line.strip('\n').split()
        if (len(line) == 6):
            if( line[0].isdigit() and line[1].lstrip('-').replace('.','',1).isdigit() and line[2].lstrip('-').replace('.','',1).isdigit() and line[5].lstrip('-').replace('.','',1).isdigit() ): 
               print float(line[0]), float(line[2])
               if (float(line[0])>1.0):
                  timestep.append(float(line[0]))
                  temp.append(float(line[1]))
                  Epair.append(float(line[2]))
                  Emol.append(float(line[3]))
       

def analyse_log(plotfile, logfile='log_0.75'):    
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    timestep, temp, Epair, Emol = read_lammps_log(logfile)


    plt.plot(timestep,temp, linewidth=1, alpha=1, label=r'Temperature, $T$')
    plt.xlabel('time [MD step]')
    plt.ylabel('temperature [s.u.]')  
    with PdfPages(plotfile+'_temp.pdf') as pdf:
         pdf.savefig()
    plt.close()
    print 'cluster info plotted in %s' % plotfile+'_temp.pdf' 

    plt.plot(timestep,Epair, linewidth=1, alpha=1, label=r'Pair energy, $E_{p}$')
    plt.xlabel('time [MD step]')
    plt.ylabel('pair energy [s.u.]')  
    with PdfPages(plotfile+'_Epair.pdf') as pdf:
         pdf.savefig()
    plt.close()
    print 'cluster info plotted in %s' % plotfile+'_Epair.pdf' 




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


# logfile ='log_rerun_0.75'
# plotfile=logfile
# analyse_cluster_log(plotfile, logfile)


# logfile ='log_0.75'
# plotfile=logfile
# analyse_log(plotfile, logfile)

# sys.exit()


# atom_phi, atom_theta1, atom_theta2 = read_atomistic_angles('/Users/mm15804/scratch/SAGE/old/atomistic_trajectory/hub-hub_angles.his')
# plot_hist(atom_phi, atom_theta1, atom_theta2, filename='atomistic_angles.pdf')
#sys.exit()

traj_data = read_dump(traj_file, min_timestep, max_timestep)
psi3_file = traj_file+'.psi3'
psi3_traj = traj2psi3(traj_data, filename=psi3_file, skip_snap=False)
#psi3_traj=[]



#snap = traj_data[-1]
# plot_file=traj_file
# traj2angles_plot(traj_data, plot_file)

# nb_no, nb_list = build_nb_list(snap)
# # print
# # print 'mol_id nb_no'
# # for i in range(len(nb_no)):
# #     print i, nb_no[i]

# # print 
# # print 'mol_id nb_list(mol_id)'
# # for i in range(len(nb_list)):  
# #     print i, nb_list[i][0:nb_no[i]]


# #print
# hub_hub_pairs = build_hub_hub_pairs(nb_list, nb_no)
# #print hub_hub_pairs 
# #print
# #print len(hub_hub_pairs)
# psi3, N_angles = orientational_order(hub_hub_pairs, snap)



for i in  range(len(traj_data)):
    if len(psi3_traj)>0: 
    	psi3 = psi3_traj[i]
    	snap = traj_data[i]
        print i, snap['step']
        tcl_file = traj_file+'.psi3.'+str(snap['step'])+".tcl"
        tcl_out = conf2tcl (snap, cluster_flag=False, box_flag=True, psi3_flag = True, psi3 = psi3['psi3']) 
        write_tcl_out (tcl_out, filename=tcl_file)
    # tcl_file = traj_file+'.cluster.'+str(snap['step'])+".tcl"
    # tcl_out = conf2tcl (snap, cluster_flag=True)
    # write_tcl_out (tcl_out, filename=tcl_file)
    tcl_file = traj_file+str(snap['step'])+".tcl"
    png_file = traj_file+str(snap['step'])+".png"
    tcl_out = conf2tcl (snap, cluster_flag=False, box_flag=False)
    write_tcl_out (tcl_out, filename=tcl_file)
    render_tcl_file (4000, 4000, png_file=png_file, tcl_out_file=tcl_file)

#movie_file=traj_file
#traj2movie(traj_data, movie_file, cluster_flag=True, rotation_flag=True)




#launchargs = 'vmd -e out.tcl'
#print "%s", launchargs
#myinput = subprocess.Popen(launchargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#mystdout,mystderr = myinput.communicate()



sys.exit()

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
        print mol1_id, mol2_id

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







