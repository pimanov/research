import struct
import os
import numpy as np


_get_i = lambda(f): struct.unpack('i', f.read(4))[0]
_get_d = lambda(f): struct.unpack('d', f.read(8))[0]

def get_line(f,line):
    n = len(line)
    _get_i(f)
    for i in range(n): 
        line[i] = _get_d(f)
    _get_i(f)


def get_vfield(f,Im,Jm,Km):
    vel = np.zeros((3, Km+2, Jm+2, Im+2), order='C')
    
    for k in range(1,Km+1):
        for j in range(1,Jm+1):
            get_line(f, vel[0,k,j,1:-1])
            get_line(f, vel[1,k,j,1:-1])
            get_line(f, vel[2,k,j,1:-1])

    return vel

def get_dcp(fn):
    i = _get_i
    d = _get_d

    with open(fn, "rb") as f:
        _,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,_ = i(f),d(f),d(f),d(f),d(f),d(f),d(f),i(f),i(f),i(f),d(f),i(f)
        vel = get_vfield(f, 2**lx, Jm, 2**lt)

    return t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,vel


def get_scp(fn):
    i = _get_i
    d = _get_d

    with open(fn, "rb") as f:
        _,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,_ = i(f),d(f),d(f),d(f),d(f),d(f),d(f),i(f),i(f),i(f),i(f),i(f)
        vel = get_vfield(f, 2**lx, Jm, 2**lt)

    return t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,vel


def get_u_Pois(lx, Xmax, Jm, epsr, lt, nsym):
    X = com.xgrid( 2**lx, Xmax)
    R = com.rgrid( Jm, epsr)
    Th = com.symthgrid( 2**lt, nsym)
    
    vel = np.zeros((3, Th.m+2, R.m+2, X.m+2), order='C')
    for k in xrange(1, Th.m+1):
        for j in xrange(1, R.m+1):
            for i in xrange(1, X.m+1):
                vel[0,k,j,i] = 1.0 - R.f[j]**2
    return vel
    

def put_line(f,line):
    n = line.shape[0]
    f.write(struct.pack("i",n*8))
    for l in line: 
        f.write(struct.pack("d",l))
    f.write(struct.pack("i",n*8))
    return 

def put_vfield(f,vel):
    Km = vel.shape[1] - 2
    Jm = vel.shape[2] - 2
    
    for k in range(1,Km+1):
        for j in range(1,Jm+1):
            put_line(f,vel[0,k,j,1:-1])
            put_line(f,vel[1,k,j,1:-1])
            put_line(f,vel[2,k,j,1:-1])
            
    return 

def put_header(f,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym):
    m = 3*4+7*8
    f.write(struct.pack("i",m))
    f.write(struct.pack("d",t))
    f.write(struct.pack("d",dt))
    f.write(struct.pack("d",Dp))
    f.write(struct.pack("d",Re))
    f.write(struct.pack("d",Xmax))
    f.write(struct.pack("d",epsr))
    f.write(struct.pack("i",lx))
    f.write(struct.pack("i",Jm))
    f.write(struct.pack("i",lt))
    f.write(struct.pack("d",nsym))
    f.write(struct.pack("i",m))
    return

def put_dcp(fn,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,vel):
    
    with open(fn,"wb") as f:
        put_header(f,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym)
        put_vfield(f,vel)
    
    return


def put_car(tmax,dt,cf,cpfn):
    f = open("pipe.car","w")
    f.write("1.e-2         -tol\n")
    f.write("10         -kprt\n")
    f.write("10000         -kwrt\n")
    f.write("%22.15f         -tmax\n" % tmax)
    f.write("%18.15f         -dtmax\n" % dt)
    f.write("%18.15f         -cf\n" % cf)
    f.write("%s\n" % cpfn)
    f.write("a0.dat")
    f.close()
    return



import com


def cs_mean_(u,(R,Th)):
    res = 0.0
    ss = 0.0
    for k in xrange(1, Th.m+1):
        for j in xrange(1, R.m+1):
            res += R.f1[j] * R.f[j] * Th.h * u[k,j]
            ss  += R.f1[j] * R.f[j] * Th.h
    return res / ss


def pipe_mean_(u, (X,R,Th)):
    res = 0.0
    ll = 0.0
    for i in xrange(1, X.m+1):
        res += cs_mean_(u[:,:,i], (R,Th)) * X.h
        ll  += X.h
    return res / ll


def cs_mean(u, Jm, epsr, lt, nsym):
    R = com.rgrid( Jm, epsr)
    Th = com.symthgrid( 2**lt, nsym)
    return cs_mean_(u, (R,Th))


def pipe_mean(u, lx, Xmax, Jm, epsr, lt, nsym):
    X = com.xgrid( 2**lx, Xmax)
    R = com.rgrid( Jm, epsr)
    Th = com.symthgrid( 2**lt, nsym)
    return pipe_mean_(u, (X,R,Th))
    


