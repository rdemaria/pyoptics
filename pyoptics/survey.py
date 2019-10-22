import numpy as np
from . import tfsdata

def rot_z(psi):
    """ Rotation around z axis. Positive angle move x into y.
    """
    c=np.cos(psi); s=np.sin(psi)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]])

def rot_x(phi):
    """ Rotation around x axis. Positive angle move y into z.
    """
    c=np.cos(phi); s=np.sin(phi)
    return np.array([[1,0,0],[0,c,-s],[0,s,c]])

def rot_y(theta):
    """ Rotation around y axis. Positive angle move z into x.
    """
    c=np.cos(theta); s=np.sin(theta)
    return np.array([[c,0,s],[0,1,0],[-s,0,c]])

def rot_mad(theta=0,phi=0,psi=0):
    """Intrinsic rotaion. Order psi(roll), phi(pitch), theta(jaw)
    """
    return rot_y(theta)@rot_x(-phi)@rot_z(psi)

def get_w_from_angles(theta,phi,psi,w):
    costhe = np.cos(theta)
    cosphi = np.cos(phi)
    cospsi = np.cos(psi)
    sinthe = np.sin(theta)
    sinphi = np.sin(phi)
    sinpsi = np.sin(psi)
    w[0, 0] = + costhe * cospsi - sinthe * sinphi * sinpsi
    w[0, 1] = - costhe * sinpsi - sinthe * sinphi * cospsi
    w[0, 2] = sinthe * cosphi
    w[1, 0] = cosphi * sinpsi
    w[1, 1] = cosphi * cospsi
    w[1, 2] = sinphi
    w[2, 0] = - sinthe * cospsi - costhe * sinphi * sinpsi
    w[2, 1] = + sinthe * sinpsi - costhe * sinphi * cospsi
    w[2, 2] = costhe * cosphi
    return w

def get_global_orbit(xs,ys,zs,theta,phi,psi,x,y,s):
    w=w_from_angles(theta,phi,psi)
    return np.array([xs,ys,ss])+w@np.array([x,y,s])

class Survey:
    mad_columns=['name', 's', 'l', 'angle',
                 'x', 'y', 'z', 'theta', 'phi', 'psi',
                 'globaltilt', 'mech_sep', 'assembly_id', 'slot_id']
    @classmethod
    def from_data(cls,data):
        self=cls()
        for cc in cls.mad_columns:
            setattr(self,cc,data[cc])
        return self

    @classmethod
    def from_mad_table(cls,mad):
        return cls.from_data(mad.table.survey)

    def get_global_orbit(self,idx,x,y,s):
        vv=p.array([self.x[idx],self.y[idx],self.z[idx]])
        ww=np.zeros((3, 3))
        get_w_from_angles(self.theta[idx],self.phi[idx],self.psi[idx],ww)
        ll=np.zeros([x,y,s])
        return vv+ww@ll

    @classmethod
    def open(cls,filename):
        return cls.from_data(tfsdata.open(filename))

