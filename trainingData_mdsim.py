""" Author: Hiqmet Kamberaj
    Copyright (2018)
    Contact: kamberaj@yahoo.co.uk   
"""
import numpy as np
import sys, random
###################################################################################
class Residue:
    def __init__(self, t, natoms, resname, atname, resid):
        self.T       = t
        self.Natoms  = natoms
        self.Ndof    = 3*natoms
        self.Resname = resname
        self.Atname = atname
        self.Z = np.zeros((self.T, self.Ndof), float)
#----------------------------------------------------------------------------------        
    def setX(self, x):
        self.X = x
    def getX(self):
        return self.X
#----------------------------------------------------------------------------------        
    def setY(self, y):
        self.Y = y
    def getY(self):
        return self.Y
       
###################################################################################
class Data:
    def __init__(self, T, fname):
        self.T = T
        self.Nres, self.Natoms, self.Residues = self.readHeaderPDB(T, fname)
        n = 0
        for residue in (self.Residues):
            n += residue.Natoms
        print n, self.Nres, self.Natoms
#----------------------------------------------------------------------------------
    def readInfoPDB(self,fname):
        try:
            finput = open(fname, 'r')
            t = -1
            prev_ires = 1
            i = 0
            for line in finput.xreadlines():
                if line == '': break
                fields= line.split()
                if fields[0] == 'TITLE' or fields[0] == 'REMARK' or fields[0] == 'CRYST1' or fields[0] == 'MODEL' or fields[0] == 'ENDMDL': continue
                if fields[0] == 'END' or fields[0] == 'TER':
                   t += 1
                   if t > self.T:
                     finput.close()
                     return 1

                   if t % 100 == 0:
                      print "# of frames so far=   ", t

                if fields[0] == 'ATOM':
                   resid = int(fields[5]) 
                   if prev_ires != resid:
                      prev_ires = resid
                      i = 0
  
                   if fields[2] == 'CA' or fields[2] == 'O' or fields[2] == 'C' or fields[2] == 'N':
                      self.Residues[resid-1].Z[t][i]   = float(fields[6])
                      self.Residues[resid-1].Z[t][i+1] = float(fields[7])
                      self.Residues[resid-1].Z[t][i+2] = float(fields[8])
                      i += 3
   
            finput.close()        
            return 1
        except IOError, (errno, strerror):
            print "(ReadInfoPdb) I/O error(%s):%s"%(errno,strerror)
        except ValueError:
            print "(ReadInfoPdb) Could not convert data to an integer."
        except:
            print "(ReadInfoPdb) Unexpected error:", sys.exc_info()[0]
            raise 
#----------------------------------------------------------------------------------        
    def readHeaderPDB(self, nframes, fname):
        N = 0
        nres=1
        atname=[]
        natoms = 0
        prev_resid = 1
        Residues=[]
        try:
            finput = open(fname, 'r')
            for line in finput.xreadlines():
                if line == '': break
                fields= line.split()
                if fields[0] == 'TITLE' or fields[0] == 'REMARK' or fields[0] == 'CRYST1' or fields[0] == 'MODEL' or fields[0] == 'ENDMDL': continue
                if fields[0] == 'END' or fields[0] == 'TER':
                   Residues.append( Residue(nframes, natoms, resname, atname, prev_resid) )
                   return nres, N, Residues
                if fields[0] == 'ATOM':
                   resid = int( fields[5] )
                   if resid != prev_resid:
                      Residues.append( Residue(nframes, natoms, resname, atname, prev_resid) )
                      prev_resid = resid
                      nres+=1 
                      atname=[]
                      natoms = 0

                   if fields[2] == 'CA' or fields[2] == 'O' or fields[2] == 'C' or fields[2] == 'N':
                      N += 1
                      natoms += 1
                      resname = fields[3]
                      atname.append( fields[2] )

            finput.close()
         
            return nres, N, Residues
        except IOError, (errno, strerror):
            print "(ReadHeaderPdb) I/O error(%s):%s"%(errno,strerror)
        except ValueError:
            print "(ReadHeaderPdb) Could not convert data to an integer."
        except:
            print "(ReadHeaderPdb) Unexpected error:", sys.exc_info()[0]
            raise 
#----------------------------------------------------------------------------------        
    def readHeaderXYZ(self, T, fname):
        try:
           finput = open(fname, 'r')
           t = 0
           nres = 1
           natoms = 0
           prev_ires = 1
           Residues = []
           atname = []
           for line in finput.xreadlines():
               fields= line.split()
               if len(fields) <= 0: 
                  break
               elif len(fields) == 1:
                  N = int(fields[0])
                  t += 1
                  if t > 1:
                     finput.close()
                     Residues.append( Residue(T, natoms, t_resname, atname, prev_ires) )
                     return nres, N, Residues
               else:
                  atIndex = int(fields[0])
                  resid   = int(fields[1]) 
                  resname = fields[2] 
                  if prev_ires != resid:
                     Residues.append( Residue(T, natoms, t_resname, atname, prev_ires) )
                     nres += 1
                     natoms = 1
                     prev_ires = resid
                     atname = []
                     atname.append(fields[3])
                  else:
                     t_resname = resname
                     atname.append(fields[3])
                     natoms += 1
                 
           finput.close()
           Residues.append( Residue(T, natoms, t_resname, atname, prev_ires) )    
           return nres, N, Residues
    
        except IOError, (errno, strerror):
               print "(readHeaderXYZ - function) I/O error(%s):%s"%(errno,strerror)
        except ValueError:
               print "(readHeaderXYZ - function) Could not convert data to an integer."
        except:
               print "(readHeaderXYZ - function) Unexpected error:", sys.exc_info()[0]
               raise
#-------------------------------------------------------------------------------
    def readInfoXYZ(self, fname):
        try:
           finput = open(fname, 'r')
           t = -1
           prev_ires = 1
           i = 0
           for line in finput.xreadlines():
               fields= line.split()
               if len(fields) <= 0: 
                  break
               elif len(fields) == 1:
                  N = int(fields[0])
                  t += 1
                  if t > self.T:
                     finput.close()
                     return 1
               else:
                  resid   = int(fields[1]) 
                  if prev_ires != resid:
                     prev_ires = resid
                     i = 0   
                  self.Residues[resid-1].Z[t][i]   = float(fields[4])
                  self.Residues[resid-1].Z[t][i+1] = float(fields[5])
                  self.Residues[resid-1].Z[t][i+2] = float(fields[6])
                  i += 3
                 
           finput.close()
           return 1
    
        except IOError, (errno, strerror):
               print "(readInfoXYZ - function) I/O error(%s):%s"%(errno,strerror)
        except ValueError:
               print "(readInfoXYZ - function) Could not convert data to an integer."
        except:
               print "(readInfoXYZ - function) Unexpected error:", sys.exc_info()[0]
               raise
#-------------------------------------------------------------------------------
    def calcTimeLag(self):
        for residue in self.Residues:
            Nk = residue.T
            mu_A = np.mean( residue.Z, 0 )
            sigma2 = np.sum( (residue.Z - mu_A) * (residue.Z - mu_A), 0 ) / float(Nk - 1)
            Tau = np.zeros(residue.Ndof, int)
            for idim in range(residue.Ndof):
                incr = 1
                t = 1
                for tt in range(Nk):
                    if t >= Nk-1: break
                    C = 0.0
                    for k in range(Nk-t):
                        C += residue.Z[k][idim]*residue.Z[k+t][idim]  
                    C /= float(Nk-t)
                    
                    C = (C - mu_A[idim]*mu_A[idim]) / sigma2[idim]
                    if C <= 0.0: break
                    t += incr
                    incr += 1
                    print idim,",",tt,",",C
                Tau[idim] = tt
                
            tau = np.max(Tau)
            print residue.Resname, " Time lag=   ", tau
            residue.tau = tau
            residue.X = np.zeros((residue.T-residue.tau, residue.Ndof), float)
            residue.Y = np.zeros((residue.T-residue.tau, residue.Ndof), float)
      
#-------------------------------------------------------------------------------
    def removeMean(self):
       for residue in self.Residues:
           xa = np.zeros( residue.Ndof, float)
           ya = np.zeros( residue.Ndof, float)
           for t in range(residue.T-residue.tau):
               xa += residue.Z[t]
               ya += residue.Z[t+residue.tau]
           xa /= (residue.T-residue.tau)
           ya /= (residue.T-residue.tau)

           for t in range(residue.T-residue.tau):
               residue.X[t] = residue.Z[t] - xa
               residue.Y[t] = residue.Z[t+residue.tau] - ya
#-------------------------------------------------------------------------------
    def scaleData(self):
        for residue in self.Residues:
            mn =  1000
            mx = -1000
            for i in range(residue.T-residue.tau):
                for j in range(residue.Ndof):
                    mn = min(mn, residue.Y[i][j])
                    mx = max(mx, residue.Y[i][j])

            for i in range(residue.T-residue.tau):
                for j in range(residue.Ndof):
                    residue.Y[i][j] = self.scale(residue.Y[i][j], mx, mn)
            residue.Ymax = mx
            residue.Ymin = mn 
            print "(Ymax, Ymin):   ", residue.Ymax, residue.Ymin
#-----------------------------------------------------------------------------------
    def reScale(self, y, ymax, ymin):
        return y*(ymax - ymin) + ymin
#-----------------------------------------------------------------------------------
    def scale(self, y, ymax, ymin):
        return  (y - ymin) / (ymax - ymin)
        
#----------------------------------------------------------------------------------        
    def packData(self,x,y):
        return zip(x, y)
#----------------------------------------------------------------------------------  
    def unpackData(self, data, T, N):  
        X = np.zeros((T, N), float)    
        Y = np.zeros((T, N), float)    
        i=0
        for x, y in data:
            X[i] = x
            Y[i] = y
            i += 1
        return X, Y
#----------------------------------------------------------------------------------        
    def prepareData(self, data, T, N, t=10, nshuffles=100, p=1.0):
        """ Prepare the data """
        for i in range(nshuffles):
            np.random.shuffle(data) 
  
        trainData_size = t
        validData_size = T - t

        X1=np.zeros((trainData_size,N), float)
        Y1=np.zeros((trainData_size,N), float)
        X2=np.zeros((T-trainData_size,N), float)
        Y2=np.zeros((T-trainData_size,N), float)

        i1=0
        i2=0
        for x, y in data:
            rn = np.random.rand(1,1)
            if rn <= p and i1 < trainData_size:
               X1[i1] = x
               Y1[i1] = y
               i1+=1
            elif i2 < validData_size:
               X2[i2] = x
               Y2[i2] = y
               i2+=1
   
        if i1 < trainData_size or i2 < validData_size:
           print "(train, validation, data - sizes):  ", i1, i2, trainData_size, validData_size, T
           print "Problem with data preparation! "
           sys.exit(1)
        else:
           training_data   = zip(X1,Y1)
           validation_data = zip(X2,Y2)
           return training_data, trainData_size, validation_data, validData_size
#-------------------------------------------------------------------------------
