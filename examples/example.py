import investor_lifespan_model as mdl
import numpy as np
import matplotlib.pyplot as plt

# Define problem parameters:
np.random.seed(0)
t0 = 30    # Current age.
W0 = 100   # Initial wealth.
δ  = 10    # Relative risk aversion parameter.
K  = 10    # Number of lifetime simulations to run.
Δt = 1/12  # Duration of each time step.
tf = 100   # Maximum age considered.

# Define after-tax income over time.
def Y(t):
    t_vec = [  29,  30,  40,  50,  60,  61,  70, 100 ] # age
    Y_vec = [  70,  70,  80,  90, 100,   0,   0,   0 ] # income
    return np.interp(t, t_vec, Y_vec)

# Define time-weighting functions for utility.
def h(t):
    t_vec = [  30,   40,   50,  60,  70,  80,  90, 100 ] # age
    h_vec = [   1,    1,    1,   1,   1,   1,   1,   1 ] # no children
    h_vec = [   1,    2,    2,   1,   1,   1,   1,   1 ] # with children
    return np.interp(t, t_vec, h_vec)**(δ-1)

# Define time-weighting functions for bequeathment.
def m(t):
    t_vec = [  30,  40,  50,  60,  70,  80,  90, 100 ] # age
    m_vec = [   0,   0,   0,    0,  0,   0,   0,   0 ] # no children
    m_vec = [   0,  20,  10,   8,   5,   3,   2,   1 ] # with children
    return np.interp(t, t_vec, m_vec)**(-(1-δ))

# Set up problem:
inv = mdl.Investor(mdl.π, mdl.G, δ, h, m, Y)
ins = mdl.Insurer(inv)
mkt = mdl.Market()
mdl = mdl.LifespanModel(inv, ins, mkt, t0=t0, tf=tf, Δt=Δt)

# Simulate several lifetimes:
res = []
for k in range(K):
    res += [mdl.simulate(W0)]

print('\nRECOMMENDATIONS:')
print('Consumption / year  :  $', int(1000*mdl.C(W0, 0)))
print('Stock/bond ratio    :   ', int(100*mdl.w(W0, 0)), '%')
print('Insurance premium   :  $', int(1000*np.max([mdl.P(W0, 0),0])))
print('Annuity income      :  $', int(1000*np.max([-mdl.P(W0, 0),0])))
print('NPV of future wages :  $', int(1000*(mdl.b_vec[0])))
print('Rel. risk aversion  :  ',  δ)
print('Abs. risk aversion  :  ',  res[0]['ARA'][0])
print('Discount factor     :  ',  int(res[0]['discount']*1000)/10, '%'  )


#PLOTTING
if True:
    plt.close('all')

    # Plot mortality statistics:
    if True:

        plt.figure()
        plt.title('MORTALITY STATISTICS')
        plt.subplot(311)
        plt.step(res[0]['t'], res[0]['π'] / res[0]['G'][0], color='k')
        plt.ylim(ymin=0)
        plt.xlim(xmin=t0)
        plt.ylabel('Death PDF')

        plt.subplot(312)
        plt.step(res[0]['t'], res[0]['G'], color='k')
        plt.ylim(ymin=0)
        plt.xlim(xmin=t0)
        plt.ylabel('Survival function')

        plt.subplot(313)
        plt.step(res[0]['t'], res[0]['λ'], color='k')
        plt.ylim(ymin=0, ymax=1)
        plt.xlim(xmin=t0)
        plt.ylabel('Force of mortality')


    # Plot risk aversion, marginal utility:
    if True:
        plt.figure()
        plt.title('PERSONAL ECONOMIC DATA')
        plt.subplot(211)
        for k in range(K):
            plt.step(res[k]['t'], res[k]['ARA'], color='.75', zorder=1)
            plt.scatter(res[k]['t'][res[k]['k_death']],
                        res[k]['ARA'][res[k]['k_death']],
                        marker='x', color='r', zorder=2)
        plt.xlim(xmin=t0)
        plt.ylim(ymin=0)
        plt.ylabel('Absolute risk aversion')

        plt.subplot(212)
        for k in range(K):
            plt.step(res[k]['t'], res[k]['dJdW'], color='.75', zorder=1)
            plt.scatter(res[k]['t'][res[k]['k_death']],
                        res[k]['dJdW'][res[k]['k_death']],
                        marker='x', color='r', zorder=2)
        plt.step(res[k]['t'], res[k]['dJdW_bar'], color='k', zorder=1)
        plt.xlim(xmin=t0)
        plt.yscale('log')
        plt.ylabel('Marginal utility of wealth')


    # Plot investment fraction, life insurance:
    if True:
        plt.figure()
        plt.title('DECISIONS')
        plt.subplot(311)
        for k in range(K):
            plt.step(res[k]['t'], res[k]['w'], color='.75', zorder=1)
            plt.scatter(res[k]['t'][res[k]['k_death']],
                        res[k]['w'][res[k]['k_death']],
                        marker='x', color='r', zorder=2)
        plt.xlim(xmin=t0)
        plt.ylabel('Stock/bond ratio')

        plt.subplot(312)
        for k in range(K):
            plt.step(res[k]['t'], res[k]['P'], color='.75', zorder=1)
            plt.scatter(res[k]['t'][res[k]['k_death']],
                        res[k]['P'][res[k]['k_death']],
                        marker='x', color='r', zorder=2)
        plt.xlim(xmin=t0)
        plt.ylabel('Insurance premium (k$)')

        plt.subplot(313)
        for k in range(K):
            plt.step(res[k]['t'], res[k]['P']/res[k]['μ'], color='.75', zorder=1)
            plt.scatter(res[k]['t'][res[k]['k_death']],
                        res[k]['P'][res[k]['k_death']] / res[k]['μ'][res[k]['k_death']],
                        marker='x', color='r', zorder=2)
        plt.xlim(xmin=t0)
        plt.ylabel('Insurance Payout (k$)')
        plt.xlabel('Age')


    # Plot wealth, consumption, and bequeathment:
    if True:
        plt.figure()
        plt.title('WEALTH USAGE')
        plt.subplot(311)
        for k in range(K):
            plt.step(res[k]['t'], res[k]['W'], color='.75', zorder=1)
            plt.scatter(res[k]['t'][res[k]['k_death']],
                        res[k]['W'][res[k]['k_death']],
                        marker='x', color='r', zorder=2)
        plt.step(res[k]['t'], res[k]['Wbar'], color='k', zorder=1)
        plt.ylim(ymin=0)
        plt.xlim(xmin=t0)
        plt.ylabel('Wealth (k$)')

        plt.subplot(312)
        for k in range(K):
            plt.step(res[k]['t'], res[k]['C'], color='.75', zorder=1)
            plt.scatter(res[k]['t'][res[k]['k_death']],
                        res[k]['C'][res[k]['k_death']],
                        marker='x', color='r', zorder=2)
        plt.ylim(ymin=0)
        plt.xlim(xmin=t0)
        plt.ylabel('Consumption (k$)')

        plt.subplot(313)
        for k in range(K):
            plt.step(res[k]['t'], res[k]['Z'], color='.75', zorder=1)
            plt.scatter(res[k]['t'][res[k]['k_death']],
                        res[k]['Z'][res[k]['k_death']],
                        marker='x', color='r', zorder=2)
        plt.xlim(xmin=t0)
        plt.ylim(ymin=0)
        plt.ylabel('Bequeathment (k$)')
        plt.xlabel('Age')



    plt.show()
