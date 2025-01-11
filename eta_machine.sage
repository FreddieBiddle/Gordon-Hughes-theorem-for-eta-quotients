NN = 1000
R.<q> = LaurentSeriesRing(QQ, NN)
from sage.modular.etaproducts import qexp_eta
eta = qexp_eta(ZZ[['q']],NN)

#Returns eta-product eta(r1*z)^s1*eta(r2*z)^s2 ... eta(rn*z)^sn. 
#Input L = [(r1, s1),..,(rn, sn)] and MM = desired order of q-expansion.

def EtaProd(L, MM):
    
    #Level of eta-product:
    N = lcm([L[i][0] for i in range(len(L))])
    
    #Weight of eta-product:
    k = (1/2)*sum(L[i][1] for i in range(len(L)))
    
    #Relevant sums involved in calculations:
    sum1 = sum(L[i][0]*L[i][1] for i in range(len(L)))
    sum2 = sum((N/L[i][0])*L[i][1] for i in range(len(L)))
    
    #To compute the character of the eta-product:
    s = prod(L[i][0]^(L[i][1]) for i in range(len(L)))
    
    #Increase level N if sum2 not divisible by 24:
    if sum2%24 != 0:
        N = N*lcm(24, sum2)
    
    #Making sure powers of q are integral:
    if sum1%24 != 0:
        print(f"Sum of d*r_d not divisible by 24.")
    
    else:
        
        #Getting denominators of cusps of congruence subgroup Gamma0(N):
        cusps = Gamma0(N).cusps()
        cusp_denoms = [cusps[i].denominator() for i in range(len(cusps)-1)]
        
        #The eta-product's q-expansion:
        Eta = q^(sum1/24)*prod(eta.subs(q=q^(L[i][0]))^(L[i][1]) for i in range(len(L)))
        
        print(f"Level {N}, weight {k}, character mod {(-1)^(k)*s}, q-expansion:")
        print('')
        print(Eta.O(MM))
        print('')
        
        #The orders of eta-product at the cusps of Gamma0(N):
        for d in cusp_denoms:
            cusp_order = (N/24)*sum((gcd(d, L[i][0])^2*L[i][1])/(gcd(d, N/d)*d*L[i][0]) for i in range(len(L)))
            print(f"Order at cusp denominator {d}: {cusp_order}")
            print("")
        
            
