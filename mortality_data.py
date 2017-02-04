import numpy as np
import scipy.interpolate as interp
import scipy.integrate as integ
import scipy as sp

G_vec = np.array([ 100000,
                   99450 ,
                   99409 ,
                   99380 ,
                   99360 ,
                   99341 ,
                   99326 ,
                   99313 ,
                   99302 ,
                   99292 ,
                   99284 ,
                   99277 ,
                   99268 ,
                   99256 ,
                   99238 ,
                   99209 ,
                   99171 ,
                   99123 ,
                   99064 ,
                   98992 ,
                   98907 ,
                   98809 ,
                   98699 ,
                   98581 ,
                   98457 ,
                   98331 ,
                   98204 ,
                   98076 ,
                   97946 ,
                   97814 ,
                   97679 ,
                   97543 ,
                   97404 ,
                   97262 ,
                   97118 ,
                   96970 ,
                   96817 ,
                   96658 ,
                   96492 ,
                   96317 ,
                   96133 ,
                   95938 ,
                   95730 ,
                   95507 ,
                   95264 ,
                   94998 ,
                   94706 ,
                   94387 ,
                   94037 ,
                   93650 ,
                   93224 ,
                   92759 ,
                   92255 ,
                   91710 ,
                   91123 ,
                   90491 ,
                   89809 ,
                   89078 ,
                   88297 ,
                   87467 ,
                   86587 ,
                   85658 ,
                   84675 ,
                   83635 ,
                   82537 ,
                   81377 ,
                   80150 ,
                   78852 ,
                   77473 ,
                   76007 ,
                   74446 ,
                   72776 ,
                   70982 ,
                   69056 ,
                   67002 ,
                   64831 ,
                   62541 ,
                   60115 ,
                   57539 ,
                   54813 ,
                   51921 ,
                   48892 ,
                   45743 ,
                   42473 ,
                   39088 ,
                   35629 ,
                   32148 ,
                   28633 ,
                   25138 ,
                   21721 ,
                   18443 ,
                   15362 ,
                   12530 ,
                   9992  ,
                   7775  ,
                   5893  ,
                   4344  ,
                   3108  ,
                   2156  ,
                   1448  ,
                   940 ]) / 100000

Δt = 1
π_vec = np.hstack([np.diff(1-G_vec)/Δt, np.zeros((1))])
tf = len(π_vec)-1
tf_pp = len(π_vec)
π_interp = interp.interp1d(range(tf_pp), π_vec, kind='quadratic')
t_vec = np.linspace(0,tf,10000)
Δt = t_vec[1] - t_vec[0]
π_vec = [π_interp(t) for t in t_vec]
π_vec = [np.max(π_vec[k],0) for k in range(len(t_vec))]
π_vec /= np.trapz(π_vec)*Δt
G_vec = 1-integ.cumtrapz(π_vec)*Δt
G_vec = np.hstack([np.ones((1)), G_vec])
G_vec = [np.max(G_vec[k],0) for k in range(len(t_vec))]

##tf = 110
##t_vec = list(range(0, tf))
##G_vec = G_vec[:tf]
#t_vec = list(range(0, 120))
#tf = 120
#
#π_vec = -np.diff(G_vec)
#π_vec = abs(π_vec)
#π_vec = np.hstack([-np.diff(G_vec), 1/sum(π_vec)*np.ones((1))])

#G = interp.interp1d(t_vec, G_vec, kind='quadratic')

def π(t):
    if t > t_vec[-1]:
        return 0
    if t < t_vec[0]:
        return 0
    else:
        return np.interp(t, t_vec, π_vec)

def G(t):
    if t > t_vec[-1]:
        return 0
    if t < t_vec[0]:
        return 1
    else:
        return np.interp(t, t_vec, G_vec)
