from numpy.polynomial.polynomial import Polynomial, polyval
import numpy as np

def integral_radial(poly1, poly2):
    p = poly1 * poly2 * Polynomial([0,0,1])
    ip = p.integ()
    I = ip(1) - ip(0)
    return I

def zernike_dict(nmax):
    zerDic = {}
    for n in range(nmax+1):
        for l in range(n,-1,-2):
            # print(n,l)
            key = (n,l)
            if l == n:
                # coef = [0 for _ in range(n+1)]
                # coef[n] = 1
                poly = xn(n)
                # zerDic[key] = poly
            elif l == n-2:
                if n >=2:
                    poly = (n+0.5)*zerDic[(n,n)] - (n-0.5)*zerDic[(n-2,n-2)]
                else:
                    poly = (n+0.5)*zerDic[(n,n)]
                # zerDic[key] = poly
            elif l == 0:
                n2 = 2*n
                M1 = (n2+1)*(n2-1)/(n+l+1)/(n-l)
                M2 = -0.5*((2*l+1)**2*(n2-1) + (n2+1)*(n2-1)*(n2-3))/(n+l+1)/(n-l)/(n2-3)
                M3 = -1*(n2+1)*(n+l-1)*(n-l-2)/(n+l+1)/(n-l)/(n2-3)
                poly = (M1*Polynomial([0,0,1]) + M2)*zerDic[(n-2,l)] + M3*zerDic[(n-4,l)]
                # zerDic[key] = poly
            else:
                L1 = (2*n+1)/(n+l+1)
                L2 = -1*(n-l)/(n+l+1)
                poly = L1*Polynomial([0,1])*zerDic[(n-1,l-1)] + L2*zerDic[(n-2,l)]

            zerDic[key] = poly

    return zerDic

def save_zernike_dict(fname, zerDic):
    pass

def div_xn(poly, n):
    """
    Divide the polynomial by x^n
    n <= m, m is the minimal index
    """
    c = poly.coef
    assert np.all(c[:n] == 0)
    return Polynomial(c[n:])

def xn(n):
    coef = [0 for _ in range(n+1)]
    coef[n] = 1
    return Polynomial(coef)

def G1nl_dict(zerDic):
    G1Dic = {}
    for key, val in zerDic.items():
        _, l = key
        iPoly = val * xn(l+2)
        iPoly = iPoly.integ()
        gPoly = div_xn(iPoly, l+1)
        G1Dic[key] = gPoly

    return G1Dic

def G2nl_dict(zerDic):
    G2Dic = {}
    for key, val in zerDic.items():
        _, l = key
        if l <= 1:
            iPoly = val * xn(1-l)
        else:
            iPoly = div_xn(val, l-1)
        iPoly = iPoly.integ()
        iPoly = -1*iPoly + iPoly(1)
        gPoly = iPoly * xn(l)
        G2Dic[key] = gPoly

    return G2Dic



if __name__ == '__main__':
    nmax = 8
    zerDic = zernike_dict(nmax)
    G1Dic = G1nl_dict(zerDic)
    G2Dic = G2nl_dict(zerDic)

    for (n,l), poly in zerDic.items():
        print((n,l),integral_radial(poly,poly), 1/(2*n+3))

    for (n,l), poly in zerDic.items():
        for i in range(n+2, nmax+1, 2):
            print((n,l),(i,l),integral_radial(poly, zerDic[(i,l)]))