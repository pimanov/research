import math
import numpy as np

class Rrt:

    def __init__(self, y0, y01):

        yp = lambda a,b: 0.01 * math.sqrt( 0.3*(b+1) ) * (a*(58.*b-21.-21.*b*b)+100.)
        y1p = lambda a,b: 1. + a * ( b - 0.45 * (b+1.)**2 )
        dyda = lambda a,b: 0.01 * math.sqrt( 0.3*(b+1) ) * ( 58.*b - 21. - 21.*b*b )
        dydb = lambda a,b: 0.01 * math.sqrt( 0.3*(b+1) ) * ( 0.5 * ( a * (58.*b - 21. - 21.*b*b) + 100. ) / (b+1) + a*(58.-42*b) )
        dy1da = lambda a,b: b - 0.45 * (b+1.)**2
        dy1db = lambda a,b: a * ( 1. - 0.9*(b+1.) )

        a = 0.0
        b = 10. / 3.*y0*y0 - 1.

        while True:
            c1 = dyda(a,b)
            c2 = dydb(a,b)
            c3 = y0 - yp(a,b)
            c4 = dy1da(a,b)
            c5 = dy1db(a,b)
            c6 = y01 - y1p(a,b)
            d  = c1*c5 - c2*c4
            da = (c3*c5 - c2*c6) / d
            db = (c1*c6 - c3*c4) / d
            a = a + da
            b = b + db
            if abs(da) < 1.e-5 and abs(db) < 1.e-5 : break


        self.aset = a
        self.bset = b

        
    i0 = lambda self, x: x * (1. + self.aset * (1.-x*x) * (self.bset-x*x) )
    
    i1 = lambda self, x: 1. + self.aset * (self.bset - 3.*x*x*(1.+self.bset) + 5.*x*x*x*x)
    

    
    
class xgrid:
    
    def __init__(self,N,L):
        self.m = N
        self.h = L / self.m
        self.n = np.array( [ i*self.h for i in range(0, self.m+2) ] )
        self.f = np.array( [ (i-0.5)*self.h for i in range(0, self.m+2) ] )

        
        
class rgrid:
    
    def __init__(self,N,epsr):
        self.m = N
        self.epsr = epsr
        h = 1. / self.m
        rrt = Rrt( 1.0, epsr )

        ros = np.array( [ j*h for j in range( 0, self.m+2 ) ] )
        self.n = np.array( [ rrt.i0(ro) for ro in ros ] )
        self.n1 = np.array( [ rrt.i1(ro)*h for ro in ros ] )

        ros = np.array( [ (j-0.5)*h for j in range( 0, self.m+2 ) ] )
        self.f = np.array( [ rrt.i0(ro) for ro in ros ] )
        self.f1 = np.array( [ rrt.i1(ro)*h for ro in ros ] )

        
        
class symthgrid:
    
    def __init__(self,N,nsym):
        self.m = N
        self.h = math.pi / (nsym * self.m)
        self.n = np.array( [ k*self.h for k in range(0, self.m+2) ] )
        self.f = np.array( [ (k-0.5)*self.h for k in range(0, self.m+2) ] )

        
        
class thgrid:
    
    def __init__(self,N):
        self.m = N
        self.h = 2.0 * math.pi / self.m
        self.n = np.array( [ k*self.h for k in range(0, self.m+2) ] )
        self.f = np.array( [ (k-0.5)*self.h for k in range(0, self.m+2) ] )

