from numpy.polynomial.polynomial import Polynomial

def zernike_radial(poly1, poly2):
    p = poly1 * poly2 * Polynomial([0,0,1])
    ip = p.integ()
    I = ip.coef.sum()
    return I

nmax = 4
zerDic = {}
for n in range(nmax+1):
    for l in range(n,-1,-2):
        # print(n,l)
        key = (n,l)
        if l == n:
            coef = [0 for _ in range(n+1)]
            coef[n] = 1
            poly = Polynomial(coef)
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
            M2 = -0.5*(2*l+1)**2*(n2-1) + (n2+1)*(n2-1)*(n2-3)/(n+l+1)/(n-l)/(n2-3)
            M3 = -1*(n2+1)*(n+l+1)*(n-l-2)/(n+l+1)/(n-l)/(n2-3)
            poly = (M1*Polynomial([0,0,1]) + M2)*zerDic[(n-2,l)] + M3*zerDic[(n-4,l)]
            # zerDic[key] = poly
        else:
            L1 = (2*n+1)/(n+l+1)
            L2 = -1*(n-l)/(n+l+1)
            poly = L1*Polynomial([0,1])*zerDic[(n-1,l-1)] + L2*zerDic[(n-2,l)]

        zerDic[key] = poly