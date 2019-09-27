from trainingData_mdsim import *
from nnclass import *
import numpy as np
####################################################################
#-------------------------------------------------------------------
def embeddedDimension(X, n, T, Mopt_Min, Mopt_Max):
    nsize = Mopt_Max - Mopt_Min + 1
    Rtol = 10.0
    count_fnn = np.zeros(nsize, int)
    im=-1
    for m in xrange(Mopt_Min, Mopt_Max+1):
        im = im + 1
        Neff = n - m * T
        for I in range(Neff-2):
            min_dist1 = 0
            for k in range(m):
                min_dist1 += square( X[I+k*T] - X[I+1+k*T] )
            ipos=I+1
            for J in xrange(I+2, Neff):
                dist1 = 0
                for k in range(m): 
                    dist1 += square( X[I+k*T] - X[J+k*T] )
                if min_dist1 > dist1:
                    ipos = J
                    min_dist1 = dist1
   
            dist2 = square( X[I+m*T] - X[ipos+m*T] )
            if (min_dist1 <= 0): min_dist1=1.0e-20
    
            R = np.sqrt( dist2 / min_dist1 )
            if (R >= Rtol): count_fnn[im] += 1

        print m, count_fnn[im]  

    Mopt_arg=Mopt_Max
    im=-1
    for m in range(Mopt_Min, Mopt_Max+1):
        im = im+1
        if (count_fnn[im] <= 0):
            Mopt_arg = m
            break
   
    return Mopt_arg
####################################################################
def timeLags(Z,N,T):
    mu_A = np.mean( Z, 0 )
    sigma2 = np.sum( (Z - mu_A) * (Z - mu_A), 0 ) / float(T - 1)
    Tau = np.zeros(N, int)
    xc = np.zeros(T, float)
    for idim in range(N):
        incr = 1
        t = 1
        for tt in range(T):
            if t >= T-1: break
            C = 0.0
            for k in range(T-t):
                C += Z[k][idim]*Z[k+t][idim]  
            C /= float(T-t)          
            C = (C - mu_A[idim]*mu_A[idim]) / sigma2[idim]
            if C <= 0.0: 
               print idim,",",tt,",",C
               break
            t += incr
            incr += 1
            print idim,",",tt,",",C
        Tau[idim] = tt

    return Tau
########################################################################
def removeMean(Z,N,T,tau):
    xa = np.zeros( N, float)
    ya = np.zeros( N, float)
    for t in range(T-tau):
        xa += Z[t]
        ya += Z[t+tau]
    xa /= float(T-tau)
    ya /= float(T-tau)

    X = np.zeros((T-tau, N), float)
    Y = np.zeros((T-tau, N), float)
    for t in range(T-tau):
        X[t] = Z[t]     - xa
        Y[t] = Z[t+tau] - ya

    return X, Y
######################################################################
def scaleOutputs(Y,N,T,a,b):
    mn =  1000
    mx = -1000
    for i in range(T):
        for j in range(N):
            mn = min(mn, Y[i][j])
            mx = max(mx, Y[i][j])

    print "(Ymin, Ymax)   ", mn, mx
    for i in range(T):
        for j in range(N):
            Y[i][j] = (b-a)*((Y[i][j] - mn) / (mx - mn)) + a

    return Y, mn, mx
######################################################################
def rescaleOutputs(Y,T,f1,f2,a,b):
    f = f1/(b-a)
    for i in range(T):
        Y[i] = f*(Y[i] - a) + f2
    return Y
#########################################################################
def prepareData(T,t,tau,N,d,X,Y,p):
    iframe=0
    iframe1=0
    x = np.zeros((t, N), float)
    y = np.zeros((t, d), float)
    x1 = np.zeros((T-tau-t, N), float)
    y1 = np.zeros((T-tau-t, d), float)
    for i in range(T-tau):
        r = np.random.rand(1,1)
        if r < p and iframe < t:
           x[iframe] = X[i]
           y[iframe] = Y[i]
           iframe += 1
        elif iframe1 < (T-tau-t):
           x1[iframe1] = X[i]
           y1[iframe1] = Y[i]
           iframe1 += 1
           

    if iframe < t or iframe1 < (T-tau-t):
       print iframe, iframe1, t, T-tau-t
       sys.exit("Decrease NSKIP!!!")
    return x,y,x1,y1
####################################################################
def myTest1(fname, TFmode):
    T = 10000
    N = 3
    d = 1
    Z = np.zeros( (T, N), float)
    try:
           finput = open(fname, 'r')
           t = 0
           for line in finput.xreadlines():
               fields= line.split()
               if len(fields) <= 0: break
               for k in range(N):
                   Z[t][k] = float(fields[k])
               t += 1
                 
           finput.close()
     
    except IOError, (errno, strerror):
           print "(readHeaderXYZ - function) I/O error(%s):%s"%(errno,strerror)
    except ValueError:
           print "(readHeaderXYZ - function) Could not convert data to an integer."
    except:
           print "(readHeaderXYZ - function) Unexpected error:", sys.exc_info()[0]
           raise

    mu_A = np.mean( Z, 0 )
    sigma2 = np.sum( (Z - mu_A) * (Z - mu_A), 0 ) / float(T - 1)
    Tau = np.zeros(N, int)
    xc = np.zeros(T, float)
    for idim in range(N):
        incr = 1
        t = 1
        for tt in range(T):
            if t >= T-1: break
            C = 0.0
            for k in range(T-t):
                C += Z[k][idim]*Z[k+t][idim]  
            C /= float(T-t)          
            C = (C - mu_A[idim]*mu_A[idim]) / sigma2[idim]
            if C <= 0.0: break
            t += incr
            incr += 1
            print idim,",",tt,",",C
        Tau[idim] = tt

    tau = np.max(Tau)
    print " Time lag=   ", tau
    xa = np.zeros( N, float)
    ya = np.zeros( N, float)
    for t in range(T-tau):
        xa += Z[t]
        ya += Z[t+tau]
    xa /= float(T-tau)
    ya /= float(T-tau)

    X = np.zeros((T-tau, N), float)
    Y = np.zeros((T-tau, N), float)
    for t in range(T-tau):
        X[t] = Z[t]     - xa
        Y[t] = Z[t+tau] - ya

    mn =  1000
    mx = -1000
    for i in range(T-tau):
        for j in range(N):
            mn = min(mn, Y[i][j])
            mx = max(mx, Y[i][j])
            mn = min(mn, X[i][j])
            mx = max(mx, X[i][j])

    print "(Ymin, Ymax)   ", mn, mx
    for i in range(T-tau):
        for j in range(N):
            Y[i][j] = (Y[i][j] - mn) / (mx - mn)
            X[i][j] = (X[i][j] - mn) / (mx - mn)

    t=1000
    nskip = int( (T-tau)/t )
    iframe=0
    x = np.zeros((t, N), float)
    y = np.zeros((t, N), float)
    for i in range(T-tau):
        if i % nskip == 0 and iframe < t:
           x[iframe] = X[i]
           y[iframe] = Y[i]
           iframe += 1

    if iframe < t:
       print iframe, t
       sys.exit("Decrease NSKIP!!!")

    Nnets = 1
    epochs= 10000
    Nets=[]
    #               0   -1   -2  -3   -4   -5  -6
    myLayer=-3
    LayerNeurons = [N,  10,   5,  d,   5,  10,  N]
    dimensions   = [N,   d,   d,  d,   N,   N,  N] 
    for inet in range(Nnets):
        Nets.append( NeuralNetwork(x, y, dimensions, LayerNeurons, TFmode) )

############# Run Independently the Ensemble of Neural Networks
    train_output_opt   = np.zeros((Nnets, t, N), float)
    encoded_output_opt = np.zeros((Nnets, t, d), float)
    inet = 0
    for net in Nets:
        net.NNrun(epochs)
        encoded_output = net.getActivation(myLayer)
        for i in range(t): 
            train_output_opt[inet][i]   = net.output_opt[i]
            encoded_output_opt[inet][i] = encoded_output[i]
        inet += 1
        print "Optimized so far  ", inet, "  Neural Networks out of ", Nnets

    train_output_mean   = np.mean(train_output_opt,0)
    train_output_std    = np.std(train_output_opt,0)
    encoded_output_mean = np.mean(train_output_opt,0)
    encoded_output_std  = np.std(train_output_opt,0)

    fid = open("data/benchmark1/io/training.csv", 'w')
    for i in range(t):
        print >> fid, \
                 y[i][0] * (mx - mn) + mn, ",", \
                 train_output_mean[i][0] * (mx - mn) + mn, ",",  \
                 train_output_std[i][0]  * (mx - mn), ",", \
                 y[i][1] * (mx - mn) + mn, ",", \
                 train_output_mean[i][1] * (mx - mn) + mn, ",",  \
                 train_output_std[i][1]  * (mx - mn), ",", \
                 y[i][2] * (mx - mn) + mn, ",", \
                 train_output_mean[i][2] * (mx - mn) + mn, ",",  \
                 train_output_std[i][2]  * (mx - mn), ",", \
                 encoded_output_mean[i][0] * (mx - mn) + mn, ",",  \
                 encoded_output_std[i][0]  * (mx - mn)
    fid.close()

    fid = open("data/benchmark1/io/encoded.csv", 'w')
    for i in range(t-2):
        print >> fid, \
                 y[i][0] * (mx - mn) + mn,  \
                 y[i][1] * (mx - mn) + mn,  \
                 y[i][2] * (mx - mn) + mn,  \
                 encoded_output_mean[i  ][0] * (mx - mn) + mn, \
                 encoded_output_mean[i+1][0] * (mx - mn) + mn, \
                 encoded_output_mean[i+2][0] * (mx - mn) + mn
    fid.close()
####################################################################
def myTest2(name, TFmode):
    T = 2000
    dt = Data(T, name)
    print dt.Nres, dt.Natoms
    if dt.readInfoPDB(name) == 1:
       dt.calcTimeLag()
    else:
       print "Error reading data information in (readInfoData)"
       sys.exit(1)
    dt.removeMean()
    dt.scaleData()
    d = 2
    ires = 0
    for residue in dt.Residues:
        ires += 1
        X = residue.getX()
        Y = residue.getY()
        data = dt.packData(X,Y)
        T = residue.T - residue.tau
        t=100
        nshuffles=1000
        p = 1.6*float(t)/float(T)

        training_data,trainData_size,validation_data,validData_size=dt.prepareData(data, T, residue.Ndof, t, nshuffles, p)
        x, y   = dt.unpackData(training_data, trainData_size, residue.Ndof)
        x1, y1 = dt.unpackData(validation_data, validData_size, residue.Ndof)
 
############# Prepare the Ensemble of Neural Networks
        Nnets = 1
        epochs= 20000
        Nets=[]
        #                    0       -1  -2 -3 -4 -5  -6  -7
        LayerNeurons = [residue.Ndof,  5, 10, d, d,  5, 10, residue.Ndof]
        dimensions   = [residue.Ndof,  d,  d, d, d,  d,  d, residue.Ndof] 
        for inet in range(Nnets):
            Nets.append( NeuralNetwork(x, y, dimensions, LayerNeurons, TFmode) )

############# Run Independently the Ensemble of Neural Networks
        train_output_opt = np.zeros((Nnets, t, residue.Ndof), float)
        encoded_output_opt = np.zeros((Nnets, t, d), float)
        inet = 0
        for net in Nets:
            net.NNrun(epochs)
            encoded_output = net.getActivation(-3)
            for i in range(t): 
                train_output_opt[inet][i] = net.output_opt[i]
                encoded_output_opt[inet][i] = encoded_output[i]
            inet += 1
            print "Resid:   ", ires, "  Optimized so far  ", inet, "  Neural Networks out of ", Nnets

        train_output_mean = np.mean(train_output_opt,0)
        train_output_std = np.std(train_output_opt,0)
        encoded_output_mean = np.mean(train_output_opt,0)
        encoded_output_std = np.std(train_output_opt,0)

        for idim in range(residue.Ndof):
            fid = open("data/benchmark2/io/training_"+str(ires)+"_"+str(idim)+".csv", 'w')
            for i in range(t):
                print >> fid, \
                         dt.reScale(y[i][idim],                 residue.Ymax, residue.Ymin), ",", \
                         dt.reScale(train_output_mean[i][idim], residue.Ymax, residue.Ymin), ",",  \
                         dt.reScale(train_output_std[i][idim],  residue.Ymax - residue.Ymin, 0.0)
            fid.close()

        fid = open("data/benchmark2/io/encoded_"+str(ires)+".csv", 'w')
        for i in range(t):
            print >> fid, \
                     dt.reScale(encoded_output_mean[i][0], residue.Ymax, residue.Ymin), ",",  \
                     dt.reScale(encoded_output_std[i][0],  residue.Ymax - residue.Ymin, 0.0), ",",  \
                     dt.reScale(encoded_output_mean[i][1], residue.Ymax, residue.Ymin), ",",  \
                     dt.reScale(encoded_output_std[i][1],  residue.Ymax - residue.Ymin, 0.0)
        fid.close()

####################################################################
if __name__ == '__main__':
    myTest1("data/benchmark1/lorenzmapode45.dat",3)
    myTest2("data/benchmark2/coords.pdb",3)
    sys.exit(1)
