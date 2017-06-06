import numpy as np
from mpmath import mp
from .AngularMomentum import Jrange, Jplus, Jminus, Jx, Jy, Jz



def _Omn(J, m, jjz):
    if m > 0:
        lp = Jplus(J)**m + Jminus(J)**m
        return (jjz*lp + lp*jjz)/4
    elif m < 0:
        m = -m
        lp = Jplus(J)**m - Jminus(J)**m
        return (jjz*lp + lp*jjz)/4j
    else:
        return jjz


def O20(J, no_constant_term=False):
    if no_constant_term:
        return mp.diag([3*Jz**2 for Jz in Jrange(J)])
    else:
        JJ = J*(J+1)
        return mp.diag([3*Jz**2 - JJ for Jz in Jrange(J)])


def O21(J):
    jx = Jx(J)
    jz = Jz(J)
    return 0.5*(jx*jz + jz*jx)


def O21s(J):
    jy = Jy(J)
    jz = Jz(J)
    return 0.5*(jy*jz + jz*jy)


def O22(J):
    return 0.5*(Jplus(J)**2 + Jminus(J)**2)


def O22s(J):
    return -0.5j*(Jplus(J)**2 - Jminus(J)**2)


def O40(J, no_constant_term=False):
    JJ = J*(J+1)
    if no_constant_term:
        return mp.diag([35*Jz**4 + 25*Jz**2 - 30*Jz**2*JJ for Jz in Jrange(J)])
    else:
        return mp.diag([35*Jz**4 + 25*Jz**2 - 30*Jz**2*JJ
                     - 6*JJ + 3*JJ**2 for Jz in Jrange(J)])


def O41(J):
    JJ = J*(J+1)
    jp = mp.diag([7*Jz**3 - (3*JJ+1)*Jz for Jz in Jrange(J)])
    jx = Jx(J)

    return 0.5*(jx*jp + jp*jx)


def O41s(J):
    JJ = J*(J+1)
    jp = mp.diag([7*Jz**3 - (3*JJ+1)*Jz for Jz in Jrange(J)])
    jy = Jy(J)

    return 0.5*(jy*jp + jp*jy)


def O42(J):
    JJ = J*(J+1)
    jp = mp.diag([7*Jz**2 - JJ - 5 for Jz in Jrange(J)])
    lp = Jplus(J)**2 + Jminus(J)**2

    return (jp*lp + lp*jp)/4


def O42s(J):
    JJ = J*(J+1)
    jp = mp.diag([7*Jz**2 - JJ - 5 for Jz in Jrange(J)])
    lp = Jplus(J)**2 - Jminus(J)**2

    return -0.25j*(jp*lp + lp*jp)


def O43(J):
    jp = Jz(J)
    lp = Jplus(J)**3 + Jminus(J)**3
    return (jp*lp + lp*jp)/4


def O43s(J):
    jpm = Jplus(J)**3 - Jminus(J)**3
    jz = Jz(J)
    return -0.25j*(jpm*jz+jz*jpm)


def O44(J):
    return 0.5*(Jplus(J)**4 + Jminus(J)**4)


def O44s(J):
    return -0.5j*(Jplus(J)**4 - Jminus(J)**4)


def O60(J, no_constant_term=False):
    JJ = J*(J+1)
    if no_constant_term:
        return mp.diag([231*Jz**6 - 315*JJ*Jz**4 + 735*Jz**4
                        + 105*JJ**2*Jz**2 - 525*JJ*Jz**2 + 294*Jz**2
                        for Jz in Jrange(J)])
    else:
        JJJ = - 5*JJ**3 + 40*JJ**2 - 60*JJ
        return mp.diag([231*Jz**6 - 315*JJ*Jz**4 + 735*Jz**4
                        + 105*JJ**2*Jz**2 - 525*JJ*Jz**2 + 294*Jz**2 + JJJ
                        for Jz in Jrange(J)])


def O61(J):
    JJ = J*(J+1)
    jp = mp.diag([33*Jz**5 - (30*JJ-15)*Jz**3
                  + (5*JJ**2 - 10*JJ + 12)*Jz for Jz in Jrange(J)])
    jx = Jx(J)

    # Note to self:
    # O61 = 1/2 * [jp, j++j-]+
    #     = (jp*(j++j-) + (j++j-)*jp)/4
    #     = (jx*jp + jp*jx)/2

    return (jx*jp + jp*jx)/2


def O61s(J):
    JJ = J*(J+1)
    jp = mp.diag([33*Jz**5 - (30*JJ-15)*Jz**3
                           + (5*JJ**2 - 10*JJ + 12)*Jz for Jz in Jrange(J)])
    jy = Jy(J)

    return (jy*jp + jp*jy)/2


def O62(J):
    JJ = J*(J+1)
    jp = mp.diag([33*Jz**4 - (18*JJ+123)*Jz**2
                           + (JJ**2 + 10*JJ + 102) for Jz in Jrange(J)])
    lp = Jplus(J)**2 + Jminus(J)**2

    return (jp*lp + lp*jp)/4


def O62s(J):
    JJ = J*(J+1)
    jp = mp.diag([33*Jz**4 - (18*JJ+123)*Jz**2
                           + (JJ**2 + 10*JJ + 102) for Jz in Jrange(J)])
    lp = Jplus(J)**2 - Jminus(J)**2

    return (jp*lp + lp*jp)/4j


def O63(J):
    JJ = J*(J+1)
    jp = mp.diag([11*Jz**3 - 3*JJ*Jz - 59*Jz for Jz in Jrange(J)])
    lp = Jplus(J)**3 + Jminus(J)**3

    return (jp*lp + lp*jp)/4


def O63s(J):
    JJ = J*(J+1)
    jp = mp.diag([11*Jz**3 - 3*JJ*Jz - 59*Jz for Jz in Jrange(J)])
    lp = Jplus(J)**3 - Jminus(J)**3

    return (jp*lp + lp*jp)/4j


def O64(J):
    JJ = J*(J+1)
    jp = mp.diag([11*Jz**2 - JJ - 38 for Jz in Jrange(J)])
    lp = Jplus(J)**4 + Jminus(J)**4

    return (jp*lp + lp*jp)/4


def O64s(J):
    JJ = J*(J+1)
    jp = mp.diag([11*Jz**2 - JJ - 38 for Jz in Jrange(J)])
    lp = Jplus(J)**4 - Jminus(J)**4

    return (jp*lp + lp*jp)/4j


def O65(J):
    return _Omn(J, 5, Jz(J))


def O65s(J):
    return _Omn(J, -5, Jz(J))


def O66(J):
    return 0.5*(Jplus(J)**6 + Jminus(J)**6)


def O66s(J):
    return -0.5j*(Jplus(J)**6 - Jminus(J)**6)

ofuncs = {
    #2
    (2,0):O20,
    
    (2,1):O21,
    (2,2):O22,

    (2,-1):O21s,
    (2,-2):O22s,

    #4
    (4,0):O40,

    (4,1):O41,
    (4,2):O42,
    (4,3):O43,
    (4,4):O44,

    (4,-1):O41s,
    (4,-2):O42s,
    (4,-3):O43s,
    (4,-4):O44s,
    
    #6
    (6,0):O60,
    
    (6,1):O61,
    (6,2):O62,
    (6,3):O63,
    (6,4):O64,
    (6,5):O65,
    (6,6):O66,

    (6,-1):O61s,
    (6,-2):O62s,
    (6,-3):O63s,
    (6,-4):O64s,
    (6,-5):O65s,
    (6,-6):O66s,
    }

def O(J, n, q, no_constant_term=False):
    if (n, q) in ofuncs:
        if q == 0:
            return ofuncs[(n,q)](J, no_constant_term)
        else:
            return ofuncs[(n,q)](J)
    else:
        raise NotImplementedError("O{}{}".format(n, q))

if __name__ == "__main__":

    def print_np_matrix(m):
        for xs in np.array(m):
            print(" ".join("{:12.9g}".format(xxs) for xxs in xs))

    print_np_matrix(O20(4))
    print_np_matrix(O43(4))
    print_np_matrix(O66(4))
