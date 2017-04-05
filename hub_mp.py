import os
import sys
#import fnmatch
import numpy as np
import subprocess
import time
import errno 
import scipy
from   scipy import special
from   multiprocessing import Pool,  cpu_count, current_process


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


 
def read_atomistic_hub_vectors(filer):
  data = []
  N_hub = 620
  nb = np.zeros((N_hub,3), dtype=np.int)
  with open(filer,'r') as f:
      while True:
        line=f.readline().strip('\n')
        line=f.readline().strip('\n')
        for i in xrange(N_hub):
          line=f.readline().strip('\n').split()
          nb[i,:] = [int(line[1])-1, int(line[2])-1, int(line[3])-1] 
        line=f.readline().strip('\n')
        line=f.readline().strip('\n')  
        line=f.readline().strip('\n').split()  
        while line:
          snap = {}
          Nvec=[]
          Cvec=[]           
          if line[0] == 'Time' :
            snap['timestep'] = line[2]
            print "reading timestep %s ns"% line[2]
            for i in xrange(N_hub):
              line=f.readline().strip('\n').split()
              Nvec.append([ float(line[1]), float(line[2]), float(line[3]) ])
              Cvec.append([ float(line[4]), float(line[5]), float(line[6]) ])
            x = ( np.array(Nvec) + np.array(Cvec) ) / 2  
            snap['coords'] = x
            snap['N'] = N_hub
          data.append(snap.copy())  
          line=f.readline().strip('\n').split() 
        return nb, data   

def plot_atomistic_psi_hist(filer='/home/mm15804/SAGE/data/hub_vectors.his'):
    nb, traj_data = read_atomistic_hub_vectors(filer)
    angles = []
    for snap in traj_data:
        N_hub = snap['N']
        x = snap['coords']
        for i in xrange(N_hub):
            j, k = nb[i,0], nb[i,1]
            vec_j = x[j,:] - x[i,:] 
            vec_k = x[k,:] - x[i,:] 
            vec_j = vec_j / np.linalg.norm(vec_j)
            vec_k = vec_k / np.linalg.norm(vec_k)
            cos =  np.dot (vec_j, vec_k)
            sin =  np.linalg.norm (np.cross (vec_j, vec_k))
            angles.append(np.degrees(np.arctan2(sin,cos)))
            j, k = nb[i,0], nb[i,2]
            vec_j = x[j,:] - x[i,:] 
            vec_k = x[k,:] - x[i,:] 
            vec_j = vec_j / np.linalg.norm(vec_j)
            vec_k = vec_k / np.linalg.norm(vec_k)
            cos =  np.dot (vec_j, vec_k)
            sin =  np.linalg.norm (np.cross (vec_j, vec_k))
            angles.append(np.degrees(np.arctan2(sin,cos)))
            j, k = nb[i,1], nb[i,2]
            vec_j = x[j,:] - x[i,:] 
            vec_k = x[k,:] - x[i,:] 
            vec_j = vec_j / np.linalg.norm(vec_j)
            vec_k = vec_k / np.linalg.norm(vec_k)
            cos =  np.dot (vec_j, vec_k)
            sin =  np.linalg.norm (np.cross (vec_j, vec_k))
            angles.append(np.degrees(np.arctan2(sin,cos)))
    pdf_file = 'atomistic_psi_angle_hist.pdf'
    if len(angles)>1:
        plot_1D_hist (angles, pdf_file)         



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

    args = "-aasamples 24 out.render -format TGA -res %s %s -o out.tga"%(res_x,res_y)
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

def render_tcl_file_LINUX (res_x, res_y, png_file="out.png", tcl_out_file="out.tcl"):
    tcl_out = "source %s\nrender Tachyon out.render\n" % (tcl_out_file)
    vmdin = os.popen("vmd -dispdev none","w")
    vmdin.write("%s"%tcl_out)
    vmdin.flush()
    
    tga_file = png_file+'.tga'
    args = "-aasamples 24 out.render -format TGA -res %s %s -o %s"%(res_x,res_y,tga_file)
    command=[]
    command.append("/home/mm15804/Downloads/vmd-1.9.3beta1/lib/tachyon/tachyon_LINUXAMD64")
    command = command + args.split()
    exe = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = exe.communicate()  
    print stdout, stderr

    # command = "convert %s -transparent black %s"% (tga_file, png_file)
    # exe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout, stderr = exe.communicate()  
    # print stdout, stderr
    
    #print tga_file, png_file, tcl_out_file
    #command = "rm %s %s %s out.render" % (
    # exe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout, stderr = exe.communicate()  
    # print stdout, stderr


def traj2movie(traj_data, movie_file='movie', cluster_flag=True, rotation_flag=True):
    image_dir = 'images/'
    make_sure_path_exists(image_dir)
    for s, snap in enumerate(traj_data):
      tcl_out = conf2tcl (snap, cluster_flag)
      write_tcl_out (tcl_out, filename="out.tcl")
      png_file = image_dir+movie_file+".%05d.png"%int(s)
      print png_file
      render_tcl_file (2048, 2048, png_file, "out.tcl")
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
          render_tcl_file (2048, 2048, png_file, "out.tcl")      
          s += 1
       for r in range(N_rot+1):
          rot_x = 0 #10 + r*360/50
          rot_y = -10 - r*360./N_rot
          write_tcl_out (tcl_out + "rotate x by %s\n rotate y by %s\n"%(rot_x,rot_y), filename="out.tcl")
          png_file = image_dir+movie_file+".%05d.png"%int(s)
          print png_file
          render_tcl_file (2048, 2048, png_file, "out.tcl")      
          s += 1
    command = 'ffmpeg -y -r 15 -i '+ image_dir+movie_file +r'.%05d.png'+' -c:v libx264 -r 45 -pix_fmt yuv420p %s.mp4'%movie_file         
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


def conf2pdb_saxs(snap, filename):
    x = snap['coords']
    N = snap['N']
    N_mol = snap['N']/16*2
    box = snap['box']
    cluster = snap['cluster']

    if len(cluster)>0:
      count = {}
      for cid in cluster:
          count[cid] = (count.get(cid, 0)) + 1 
      max_key = max(count, key=lambda k: count[k]) 

    print 'largest cluster with index %s has %s atoms'%(max_key, count[max_key])
    xm =  np.zeros((N_mol,3))
    for i in range(N_mol):
      xm [i,:] = x [ get_helix_COM_atom_id(i), :]
    
    res = ''
    cnt = 0
    for i in range(N_mol):
        if len(cluster)>0:
           if not cluster[i]==max_key: continue 
        #xm [i,:] = xm [ i, :] - np.min(xm)
        res += "ATOM  %5d%4s  %3s %c%4d%c   %8.1f%8.1f%8.1f%6.2f%6.2f\n" % (cnt, "CA", "ASP", 'A', cnt % 10000,' ',xm[cnt,0]*10,xm[cnt,1]*10,xm[cnt,2]*10,1,20.0)
        cnt += 1
    f = open(filename, "w")
    f.write(res)
    print 'PDB file is written to %sfor SAXS analysis with CRYSOL'%filename
    print 'usage: crysol %s -sm 1.0 -ns 1000 -fb 18 -lm 25'%filename
    f.close()


def read_spherical_packings_www(url='http://neilsloane.com/ICOSP/ipack.3.932.txt'):
    import urllib
    f = urllib.urlopen(url)
    cnt = 0
    x = []
    lines=f.readlines()
    N = len(lines)/3
    x = np.zeros((N, 3))
    for i, line in enumerate(lines):
      x[i/3, i%3] = line.strip('\n')
    
    x = np.array(x)

    #move COM to origin   
    x = np.array(x) 
    COM = [sum(p)/len(p) for p in zip(*x)]
    print "     COM moved to origin (was %s)"%COM
    for i in range(len(x)):
        x[i,:] = x[i,:] - COM 

    x = x * 32.1451   # this is to get a same size sphere as Sesh's    

    snap = {}
    snap ['coords'] = x
    snap ['N'] = N
    snap ['step'] = 'sesh_SAGE with midpoints at every polygon face'
    snap['box'] = np.array ( [10*np.max(x), 10*np.max(x), 10*np.max(x)] )
    snap['traj'] = url
    return snap

def read_spherical_packings_www_from_file(filename):
    f = open(filename, 'r')
    cnt = 0
    x = []
    lines=f.readlines()
    N = len(lines)/3
    x = np.zeros((N, 3))
    for i, line in enumerate(lines):
      x[i/3, i%3] = line.strip('\n')
    
    x = np.array(x)

    #move COM to origin   
    x = np.array(x) 
    COM = [sum(p)/len(p) for p in zip(*x)]
    print "     COM moved to origin (was %s)"%COM
    for i in range(len(x)):
        x[i,:] = x[i,:] - COM 

    x = x * 32.1451   # this is to get a same size sphere as Sesh's    

    snap = {}
    snap ['coords'] = x
    snap ['N'] = N
    snap ['step'] = 'sesh_SAGE with midpoints at every polygon face'
    snap['box'] = np.array ( [10*np.max(x), 10*np.max(x), 10*np.max(x)] )
    snap['traj'] = filename
    return snap


def make_sc_sheet (a, L):
  Nrow = int(L/a)
  x = []
  for i in range(-Nrow/2, Nrow/2):
    for j in range(-Nrow/2, Nrow/2):
      x.append([i*a, j*a, 0])
  snap = {}
  snap ['coords'] = np.array(x)
  snap ['N'] = len(x)  
  snap ['step'] = 'sc_sheet on xy plane at z=0 with a=%s'%a

  x = np.array(x)
  #move COM to origin   
  x = np.array(x) 
  COM = [sum(p)/len(p) for p in zip(*x)]
  print "     COM moved to origin (was %s)"%COM
  for i in range(len(x)):
      x[i,:] = x[i,:] - COM 
  snap = {}
  snap ['coords'] = x
  snap ['N'] = len(x)
  snap ['step'] = 'sc_sheet on xy plane at z=0 with a=%s, N=%s'% (a, len(x))
  snap['box'] = np.array ( [10*np.max(x), 10*np.max(x), 10*np.max(x)] )
  snap['traj'] = 'sc_sheet'
  print 'sc_sheet with a=%s and N=%s is generated at z=0' % (a, snap['N'])
  return snap

def make_sc_cube (a, L):
  Nrow = int(L/a)
  x = []
  for i in range(-Nrow/2, Nrow/2):
    for j in range(-Nrow/2, Nrow/2):
       for k in range(-Nrow/2, Nrow/2):
          x.append([i*a, j*a, k*a])
  
  x = np.array(x)
  #move COM to origin   
  x = np.array(x) 
  COM = [sum(p)/len(p) for p in zip(*x)]
  print "     COM moved to origin (was %s)"%COM
  for i in range(len(x)):
      x[i,:] = x[i,:] - COM 
  snap = {}
  snap ['coords'] = x
  snap ['N'] = len(x)
  snap ['step'] = 'sc_cube with a=%s, N=%s'% (a, len(x))
  snap['box'] = np.array ( [10*np.max(x), 10*np.max(x), 10*np.max(x)] )
  snap['traj'] = 'sc_cube'

  print 'sc_cube with a=%s and N=%s is generated' % (a, snap['N'])
  return snap


def make_random_spherical_shell(R, N):
  np.random.seed()
  x = []
  for i in range(N):
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x1 = R* np.sin( theta) * np.cos( phi )
    x2 = R* np.sin( theta) * np.sin( phi )
    x3 = R* np.cos( theta )
    x.append([x1, x2, x3]) 
  x = np.array(x)
  snap = {}
  snap ['coords'] = x
  snap ['N'] = len(x)
  snap ['step'] = 'random spherical shell with R=%s, N=%s'% (R, len(x))
  snap['box'] = np.array ( [10*np.max(x), 10*np.max(x), 10*np.max(x)] )
  snap['traj'] = 'make_random_spherical_shell'
  print 'generated a random spherical shell with R=%s, N=%s'% (R, len(x))
  return snap    


def read_sesh_SAGE(filer):
    with open(filer,'r') as f:
      x = []
      while True:
        line=f.readline().strip('\n').split()
        if not line: break
        x.append([float(line[0]), float(line[1]), float(line[2])])

    #move COM to origin   
    x = np.array(x)
    COM = [sum(p)/len(p) for p in zip(*x)]
    print "     COM moved to origin (was %s)"%COM
    for i in range(len(x)):
        x[i,:] = x[i,:] - COM 

    snap = {}
    snap ['coords'] = x/10
    snap ['N'] = len(x)
    snap ['step'] = 'sesh_SAGE'
    snap ['traj'] =  filer
    snap['box'] = np.array ( [10*np.max(x), 10*np.max(x), 10*np.max(x)] )

    return snap


def sesh2pdb_saxs(sesh_file='sesh_SAGE.txt', pdb_file='sesh_SAGE.pdb'):
    snap = read_sesh_SAGE(sesh_file)
    xm = snap['coords']
    print np.min(xm)
    N_hub = snap['N']

    res = ''
    for i in range(N_hub):      
        res += "ATOM  %5d%4s  %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i, "CA", "ASP", 'A', i % 10000,' ',xm[i,0],xm[i,1],xm[i,2],1,20.0)
    f = open(pdb_file, "w")
    f.write(res)
    print 'PDB file is written to %sfor SAXS analysis with CRYSOL'%pdb_file
    print 'usage: crysol %s -sm 1.0 -ns 1000 -fb 18 -lm 25'%pdb_file
    f.close()     



        

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


def get_helix_COM_atom_id(mol_id):
  atom_id = []
  atom1 = int(mol_id/2)*16
  if (mol_id%2 == 0):
      atom_id = atom1+1 
  else:
      atom_id = atom1+11
  return atom_id   


def snap2dihedrals_all(snap):  
    phi_molmol, theta1_molmol, theta2_molmol = snap2dihedrals_molmol(snap)
    phi_hubhub, theta1_hubhub, theta2_hubhub = snap2dihedrals_hubhub(snap)
    return phi_molmol, theta1_molmol, theta2_molmol, phi_hubhub, theta1_hubhub, theta2_hubhub 

def snap2dihedrals_molmol(snap):
    phi=[]
    theta1=[]
    theta2=[]
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N_mol = snap['N']/16*2
    
    for mol in range(0, N_mol, 2):
         x1,x2,x3,x4 = get_points( x[get_helix_atom_ids(mol),:] ,  x[get_helix_atom_ids(mol+1),:])
         dihedral, bend1, bend2 = get_angles(x1,x2,x3,x4,box)
         phi.append(dihedral)
         theta1.append(bend1)
         theta2.append(bend2)
    return phi, theta1, theta2  


def snap2dihedrals_hubhub(snap):
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N_mol = snap['N']/16*2

    phi=[]
    theta1=[]
    theta2=[]

    nb_no, nb_list = build_nb_list(snap)
    hub_hub_pairs =  build_hub_hub_pairs(nb_list, nb_no)

    for i, pairs in enumerate(hub_hub_pairs):
         for j, pair in enumerate(pairs):
                  mol1 = pair[0]
                  mol2 = pair[1]
                  if (mol1 >= mol2): continue
                  x1,x2,x3,x4 = get_points( x[get_helix_atom_ids(mol1),:] ,  x[get_helix_atom_ids(mol2),:])
                  dihedral, bend1, bend2 = get_angles(x1,x2,x3,x4,box)
                  phi.append(dihedral)
                  theta1.append(bend1)
                  theta2.append(bend2)
    return phi, theta1, theta2              
   

def snap2CG(snap):
    x = snap['coords']
    p_type = snap['p_type']
    box = snap['box']
    step = snap['step']
    N_mol = int(snap['N']/16*2)
    cluster = snap['cluster']

    if len(cluster)>0:
      count = {}
      for cid in cluster:
          count[cid] = (count.get(cid, 0)) + 1 
      max_key = []    
      for key, value in sorted(count.iteritems(), key=lambda (k,v): (v,k)):
           max_key.append(key)
      #max_key = max(count, key=lambda k: count[k])

      print 'the largest cluster with index %s has %s atoms'%(max_key[0], count[max_key[0]])
      print 'the second largest cluster with index %s has %s atoms'%(max_key[1], count[max_key[1]])
   
      break
   

    CGsnap=[]
    # simple version
    for mol in range(0, N_mol, 2):
        mol_pair = {}
        # trimer CC coordinate



        CCt_a = x[get_helix_atom_ids(mol)[2],:] - x[get_helix_atom_ids(mol)[0],:]
        CCt_a = CCt_a /np.linalg.norm(CCt_a)
        CCt_x = x[get_helix_atom_ids(mol)[0],:]
        CCt_type = p_type[get_patch_atom_ids(mol)[0][0]] #is always 2 (trimer)
        # dimer CC coordinate
        CCd_a = x[get_helix_atom_ids(mol+1)[2],:] - x[get_helix_atom_ids(mol+1)[0],:]
        CCd_a = CCd_a /np.linalg.norm(CCd_a)
        CCd_x = x[get_helix_atom_ids(mol+1)[0],:]
        CCd_type = p_type[get_patch_atom_ids(mol+1)[0][0]] #is x or y for A and B CCdimers

        mol_pair['tri_a'] = CCt_a
        mol_pair['tri_x'] = CCt_x
        mol_pair['tri_type'] = CCt_type
        mol_pair['di_a'] = CCd_a
        mol_pair['di_x'] = CCd_x
        mol_pair['di_type'] = CCd_type
        mol_pair['visiblity'] = -1
        if len(cluster) > 0 and not cluster[get_helix_atom_ids(mol)[0]] == max_key:
             mol_pair['visiblity'] = 1


        CGsnap.append(mol_pair.copy())

    return CGsnap

    #nb_no, nb_list = build_nb_list(snap)
    #hub_hub_pairs =  build_hub_hub_pairs(nb_list, nb_no)  
    # for i, pairs in enumerate(hub_hub_pairs):
    #      for j, pair in enumerate(pairs):
    #               mol1 = pair[0]
    #               mol2 = pair[1]
    #               if (mol1 >= mol2): continue 
    #               print pair
    #               print x[get_helix_atom_ids(mol1),:] ,  x[get_helix_atom_ids(mol2),:]  
    #               sys.exit()




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



def snap2psi3(snap, traj_file, skip_snap=True, tcl_write_flag=True, psi3_angle_hist_flag=True):
        psi3_snap = {}
        step = snap['step']
        psi3_snap['step'] = step
        # if (skip_snap and step <= last_timestep): 
        #     print 'skipping %s' %step
        #     continue
        nb_no, nb_list = build_nb_list(snap)
        hub_hub_pairs = build_hub_hub_pairs(nb_list, nb_no)
        psi3, N_angles, my_count, psi3_vec, angles = orientational_order(hub_hub_pairs, snap)
        out='%s %s %s %s\n' % (step, psi3, N_angles, my_count)
        #print '%s %s %s %s' % (step, psi3, N_angles, my_count)
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
            if len(angles)>1:
               plot_1D_hist (angles, pdf_file) 

        return out, psi3_vec, angles




def traj2psi3_old(traj_data, filename='psi3_file', skip_snap=True, tcl_write_flag=True, psi3_angle_hist_flag=True):
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
    all_angles=[]
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
        all_angles += angles
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
            if len(angles)>1:
               angle_file = 'psi3_angle_hist/'+traj_file+'.angles'
               plot_1D_hist (angles, pdf_file) 
               g = open (angle_file, 'a+')
               for j, val in enumerate(angles):
                 g.write("%s %s\n"%(str(snap['step']), val)) 
               g.flush()
               g.close()  


        psi3_snap['psi3']=psi3_vec
        psi3_traj.append(psi3_snap.copy())
        f.write (out)
        f.flush ()



    if (psi3_angle_hist_flag and len(all_angles)>1):
        make_sure_path_exists('psi3_angle_hist')
        pdf_file = 'psi3_angle_hist/'+traj_file+'.total_hist.'+str(traj_data[0]['step'])+'_'+str(traj_data[-1]['step'])+".pdf"
        plot_1D_hist (all_angles, pdf_file)        
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


def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def plot_hist(x,y,z, filename, x_lim):
    import matplotlib as mpl
    mpl.use('Agg')  
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.optimize import curve_fit
    from scipy import stats, integrate


    import pylab 

    fig = pylab.figure(0, figsize = (6,4))
    fig_width_pt = 16*246.0  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]
    params = {'backend': 'ps',
              'axes.labelsize': 14,
              #'text.fontsize': 14,
              'legend.fontsize': 10,
              'xtick.labelsize': 11,
              'ytick.labelsize': 11,
              #'text.usetex': True,
              'figure.figsize': fig_size,
              'text.latex.preamble': [r"\usepackage{amstext}", r"\usepackage{mathpazo}"]}
    pylab.rcParams.update(params)

    (x, y, z) = (np.array(x), np.array(y), np.array(z))
    (x, y, z) = (x[~is_outlier(x)], y[~is_outlier(y)], z[~is_outlier(z)])

    #print x,y,z
    #(y,z) = (z,y)
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
    plt.title(r'$\mu_{\phi}=%.2f,\ \sigma_{\phi}=%.2f,\ \ \ \mu_{\theta_1}=%.2f,\ \sigma_{\theta_1}=%.2f,\ \ \ \mu_{\theta_2}=%.2f,\ \sigma_{\theta_2}=%.2f,\ \ \ N=%d$' \
              %(mu_x, sigma_x, mu_y, sigma_y, mu_z, sigma_z, len(x)),  fontsize=10)
    plt.legend(loc='upper right')
    if len(x_lim)>0 : plt.xlim(x_lim)
    #plt.grid(True)
    #plt.show() 


    fig.tight_layout()
    with PdfPages(filename) as pdf:
          pdf.savefig(transparent=True)
          plt.close()
          print "angel histograms are saved in %s" % filename     



def plot_1D_hist(x, filename='hist.pdf'):
    import matplotlib as mpl
    mpl.use('Agg')  
    #import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
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
    with PdfPages(filename) as pdf:
         pdf.savefig(transparent=True)
    plt.close()   
    print "histogram saved in %s" % filename    


def plot_1D_hist_noX(x, filename='hist.png'):
    import matplotlib as mpl
    mpl.use('Agg')
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
    #      pdf.savefig(transparent=True)
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
         pdf.savefig(transparent=True)
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
         pdf.savefig(transparent=True)
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
         pdf.savefig(transparent=True)
    sns.plt.close() 
              


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
         pdf.savefig(transparent=True)
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
         pdf.savefig(transparent=True)
    plt.close()
    print 'cluster info plotted in %s' % plotfile+'_temp.pdf' 

    plt.plot(timestep,Epair, linewidth=1, alpha=1, label=r'Pair energy, $E_{p}$')
    plt.xlabel('time [MD step]')
    plt.ylabel('pair energy [s.u.]')  
    with PdfPages(plotfile+'_Epair.pdf') as pdf:
         pdf.savefig(transparent=True)
    plt.close()
    print 'cluster info plotted in %s' % plotfile+'_Epair.pdf' 





#************************************************************************************************


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


def cylinder_form_factor_integrand (s, R, L, a):
    f =  2*scipy.special.jn(0, s*R*np.sin(a)) * np.sin(s*L*np.cos(a/2.0))/( s*R*np.sin(a) *  s*L*np.cos(a/2.0)   )
    f = (f ** 2) * np.sin(a)
    return f 

def cylinder_form_factor(s, R, L):
    from scipy import integrate
    # Na = 10
    # da = (np.pi/2-1e-8)/ (Na - 1)
    # a = [1e-8 + i*da for i in range(Na)]
    # integrand = map(lambda x: cylinder_form_factor_integrand(s, R, L, x), a)
    f = scipy.integrate.quad(lambda x: cylinder_form_factor_integrand(s, R, L, x), 1e-10 ,  np.pi/2, epsrel = 1e-10 )
    return f

def core_shell_form_factor(s, Rc, Rs):
    # see http://www.soft-matter.uni-tuebingen.de/teaching/SASTutorial.pdf
    rho_c = 0.
    rho_s = 1.
    rho_solv = 0.
    V_s = (4*np.pi/3)* Rs**3
    V_c = (4*np.pi/3)* Rc**3

    f = 3 * V_c * (rho_c - rho_s   ) * scipy.special.j1( s*Rc ) / (s*Rc) + \
        3 * V_s * (rho_s - rho_solv) * scipy.special.j1( s*Rs ) / (s*Rs)
    return f    

def plot_core_shell_form_factor():
    Rc = 319.4
    Rs = 323.4
    
    Ns = 10000
    ds = (2 - 0.001)/(Ns-1)
    s = [0.001+i*ds for i in range(Ns)]
    f = map(lambda x: core_shell_form_factor(x, Rc, Rs), s)

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    plt.plot(s,f)
    plt.xlabel('q [1/A]')
    plt.ylabel('core shell form factor [A.U.]')

    plt.title(r'$R_c=%.2f,\ R_s=%.2f$' %(Rc, Rs),  fontsize=10)

    #plt.legend(loc='upper right')
    #plt.axis([40, 160, 0, 0.03])
    #plt.grid(True)
    #plt.show() 
    with PdfPages('core_shell_form_factor.pdf') as pdf:
         pdf.savefig(transparent=True)
    plt.close()  

    file =open('core_shell_form_factor.txt', 'w')
    for i in range(len(s)):
      file.write("%s %s\n"% (s[i], f[i]))
    file.close()  



def form_factor(s, form='sphere'):
    if form=='sphere' :
       f =  sphere_form_factor(s, 1.0)
    elif form=='cylinder':
       f =  cylinder_form_factor(s, 1.0, 3.0)
    else:
       print "    unrecognized form_factor typr (%s)"%form   
    return f

def A_lm(s_mag, spherical_coords, l, m, form):
    Alm_vector = form_factor(s_mag, form) * scipy.special.jn(l, s_mag * spherical_coords[:,0]) * \
                 np.conj ( scipy.special.sph_harm(m, l, spherical_coords[:,2], spherical_coords[:,1]) )
    Alm = 4.0 * np.pi * (complex(0,1)**l) * np.sum (Alm_vector) 
    return Alm  

def get_Aa2(s_vec, spherical_coords, lmax, form):
    A = 0
    for l in range(lmax+1):
        for m in range(-l,l+1):
            A += A_lm( s_vec[0], spherical_coords, l, m, form) * scipy.special.sph_harm(m, l, s_vec[2], s_vec[1])
    return A * np.conj(A)

def get_saxs_intensity(s_mag, snap, my_model_flag=True, form='sphere'):
    np.random.seed()
    start = time.time()
    lmax = 15
    Ndir = 100

    x = snap['coords']
    if my_model_flag:
        N_mol = snap['N']/16 * 2
        step = snap['step']    
        # TODO : take all 3 particles in the helix into account not just the middle one
        xm =  np.zeros((N_mol,3))
        for i in range(N_mol):
            xm [i,:] = x [ get_helix_COM_atom_id(i), :]
    else: 
       N_mol = snap['N']
       step = snap['step']     
       xm = x

    #move COM to origin   
    COM = [sum(p)/len(p) for p in zip(*xm)]
    print "     COM moved to origin (was %s)"%COM
    for i in range(len(xm)):
        xm[i,:] = xm[i,:] - COM 

    spherical_coords = Cartesian2Spherical(xm)
    
    sum_I = 0.0
    for i in range(Ndir):
        print i
        s = random_unit_vector() * s_mag
        s_vec = Cartesian2Spherical(np.array([s]))[0]
        sum_I += get_Aa2(s_vec, spherical_coords, lmax, form)
    end = time.time()
    print("    [trajectory timestep %s]: averaging I(s) for s = %s over %d directions for %d molec. took %s (s). {process %s}" \
        % (step, s_mag, Ndir, N_mol, end-start, current_process().pid))
    return sum_I/Ndir





def creat_mesh(s_mag, spherical_coords, lmax, Ntheta=50, Nphi=100):
    start = time.time()
    dtheta = (np.pi - 0)/(Ntheta-1) 
    dphi = (np.pi + np.pi)/(Nphi-1)
    N_mol = len(spherical_coords)

    phi_vec = [- np.pi + j*dphi  for j in range(Nphi)]
    #theta_vec = [ j*dtheta  for j in range(Ntheta)]

    # f_mesh = np.zeros(Ns)
    # for i in range(Ns):
    #     s = smin + i*ds
    #     f_mesh[i] = form_factor(s)


    ylm_mesh = np.zeros((lmax+1, 2*lmax+1, Ntheta, Nphi), dtype=complex)
    for l in range(lmax+1):
        for m in range (-l, l+1):
                mi = m + lmax
                for i in range(Ntheta):
                        theta = i*dtheta
                        ylm_mesh[l,mi,i,:] = scipy.special.sph_harm(m, l, phi_vec, theta)

    ylm_conj_coords_mesh = np.zeros((lmax+1, 2*lmax+1, N_mol), dtype=complex)
    for l in range(lmax+1):
        for m in range (-l, l+1):
            mi = m + lmax
            ylm_conj_coords_mesh[l,mi,:] = np.conj ( scipy.special.sph_harm(m, l, spherical_coords[:,2], spherical_coords[:,1])  )


    Jl_mesh =  np.zeros((lmax+1, N_mol ), dtype=complex)  
    for l in range(lmax+1):
                #dummy = scipy.special.jn(l, s * spherical_coords[:,0])
                Jl_mesh [l, : ] = scipy.special.jn(l, s_mag * spherical_coords[:,0])

    end = time.time()
    print "    creating mesh took %s (s)" % (end-start)
    return (ylm_mesh, ylm_conj_coords_mesh, Jl_mesh)

def A_lm_mesh(s_mag, l, mi, mesh, s_inds, N_mol, form):
    (ylm_mesh, ylm_conj_coords_mesh, Jl_mesh) = mesh
    Alm_vector = form_factor(s_mag, form) * Jl_mesh[l, list(range(0, N_mol)) ] * \
                 ylm_conj_coords_mesh[l, mi, list(range(0, N_mol))  ] 
    Alm = 4.0 * np.pi * (complex(0,1)**l) * np.sum (Alm_vector) 
    return Alm      

def get_Aa2_mesh(s_mag, mesh, lmax, s_inds, N_mol, form):
    (ylm_mesh, ylm_conj_coords_mesh, Jl_mesh) = mesh
    A = 0
    # TODO : use map() => elementwise calculation of A
    for l in range(lmax+1):
        for m in range(-l,l+1):
            mi = m + lmax
            A += A_lm_mesh(s_mag, l, mi, mesh, s_inds, N_mol, form) * \
                 ylm_mesh[l, mi, s_inds[1], s_inds[2]]
    return A * np.conj(A)


def get_saxs_intensity_mesh(s_mag, snap, my_model_flag=True, form='sphere'):
    np.random.seed()
    start = time.time()
    lmax = 15
    Ndir = 200
    (Ntheta, Nphi) = (60, 120)   # mesh size
    dtheta = (np.pi - 0)/(Ntheta-1) 
    dphi = (np.pi + np.pi)/(Nphi-1)

    x = snap['coords']
    if my_model_flag:
        N_mol = snap['N']/16 * 2 * 3
        step = snap['step']    
        # take all 3 particles in the helix into account not just the middle one
        xm =  np.zeros((N_mol,3))
        cnt = 0
        for i in range(N_mol/3):   
            xm [cnt:cnt+3, :] = x [ get_helix_atom_ids(i), :]
            cnt += 3
    else: 
       N_mol = snap['N']
       step = snap['step']     
       xm = x

    #move COM to origin   
    COM = [sum(p)/len(p) for p in zip(*xm)]
    print "    COM moved to origin (was %s)"%COM
    for i in range(len(xm)):
        xm[i,:] = xm[i,:] - COM 
    

    spherical_coords = Cartesian2Spherical(xm)
    mesh = creat_mesh(s_mag, spherical_coords, lmax, Ntheta, Nphi)
    
    sum_I = 0.0
    for i in range(Ndir):
        s = random_unit_vector() * s_mag
        s_vec = Cartesian2Spherical(np.array([s]))[0]
        s_inds = [0, s_vec[1]/dtheta,  s_vec[2]/dphi]  # index of theta and phi on the grid. the first element will not be used.
        sum_I += get_Aa2_mesh(s_mag, mesh, lmax, s_inds, N_mol, form)
    end = time.time()
    print("    [trajectory timestep %s]: averaging I(s) for s = %s over %d directions for %d molec. took %s (s). {process %s}" \
        % (step, s_mag, Ndir, N_mol, end-start, current_process().pid))
    return sum_I/Ndir


#************************************************************************************************







