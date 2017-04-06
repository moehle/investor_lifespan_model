import numpy as np

import matplotlib.pyplot as plt

class LifespanModel():

    def __init__(self, inv, ins, mkt, t0=30, tf=120, Δt=1/3):
        self.inv = inv
        self.ins = ins
        self.mkt = mkt
        self.t0 = t0
        self.tf = tf
        self.Δt = Δt
        self.t_vec = np.arange(t0, tf, Δt)
        self.K = len(self.t_vec)
        
        self._solve()


    def _solve(self):

        r = self.mkt.r
        α = self.mkt.α
        σ = self.mkt.σ
        γ = self.inv.γ
        δ = self.inv.δ
        K = self.K

        t_vec, π_vec, h_vec, m_vec, Y_vec, G_vec, λ_vec = \
            self.inv.get_data_vectors(self.t0, self.tf, self.Δt)

        t_vec, μ_vec, η_vec, H_vec = \
            self.ins.get_data_vectors(self.t0, self.tf, self.Δt)
       
        # AS IT APPEARS IN THE PAPER:
        # k_vec = (λ_vec / μ_vec)**(γ/δ) * (λ_vec * m_vec)**(1/δ) + h_vec**(1/δ)

        # MY VERSION
        k_vec_me = (λ_vec / μ_vec)**(γ/δ) * λ_vec * m_vec**(1/δ) + h_vec**(1/δ)

        # SCOTT VERSION
        #k_vec = μ_vec**(-γ/δ) * (λ_vec * m_vec)**(1/δ) + h_vec**(1/δ)
        k_vec = (1/μ_vec)**(γ/δ) * (λ_vec * m_vec)**(1/δ) + h_vec**(1/δ)

        #t1 = λ_vec**(γ/δ) * λ_vec
        t1 = λ_vec**(γ/δ+1)
        t2 = λ_vec**(1/δ)
        print(γ/δ + 1)
        print(1/δ)

        #print(k_vec_me - k_vec)
        print(np.linalg.norm(t1 - t2))

        A = np.zeros((K,K))
        v = (α - r)**2 / (2*δ*σ**2)
        for k,t in enumerate(t_vec):
            for j,s in enumerate(t_vec):
                if s >= t:
                    A[k,j] = (k_vec[j]*G_vec[j]/G_vec[k] 
                              * np.exp((γ/δ)*(v+r)*(s-t) + (γ/δ)*(H_vec[k] - H_vec[j])))
                    #A[k,j] = (k_vec[j]*G_vec[j]/G_vec[k])
        a_vec = (self.Δt * np.sum(A, 1))**δ

        B = np.zeros((K,K))
        for k,t in enumerate(t_vec):
            for j,s in enumerate(t_vec):
                if s >= t:
                    B[k,j] = (Y_vec[j]*G_vec[j]/G_vec[k] 
                              * np.exp(-r*(s-t) - (H_vec[k] - H_vec[j])) )
        b_vec = self.Δt * np.sum(B, 1)

        self.J = lambda W, k:  a_vec[k]/γ * (W + b_vec[k])**γ
        self.C = lambda W, k:  (h_vec[k]/a_vec[k])**(1/δ) * (W + b_vec[k])
        self.Z = lambda W, k:  (m_vec[k]*λ_vec[k] / (a_vec[k]*μ_vec[k]))**(1/δ) \
                                * (W + b_vec[k])
        self.P = lambda W, k:  μ_vec[k] * (self.Z(W,k) - W)
        self.w = lambda W, k:  ((α - r)/(δ*σ**2)*(W + b_vec[k]))/W

        self.t_vec = t_vec
        self.t_vec = t_vec
        self.π_vec = π_vec
        self.h_vec = h_vec
        self.m_vec = m_vec
        self.Y_vec = Y_vec
        self.G_vec = G_vec
        self.λ_vec = λ_vec
        self.t_vec = t_vec
        self.μ_vec = μ_vec
        self.η_vec = η_vec
        self.H_vec = H_vec
        self.a_vec = a_vec
        self.b_vec = b_vec
        self.k_vec = k_vec
    

    def plot_functions(self, W_vec):
        J_mat = np.zeros(len(W_vec), self.K)
        C_mat = np.zeros(len(W_vec), self.K)
        P_mat = np.zeros(len(W_vec), self.K)
        w_mat = np.zeros(len(W_vec), self.K)
        for k,t in enumerate(np.arange(self.t0, self.tf, self.Δt)):
            for j,W in enumerate(W_vec):
                J_mat[j,k] = self.J[W,k]
                C_mat[j,k] = self.C[W,k]
                P_mat[j,k] = self.P[W,k]
                w_mat[j,k] = self.w[W,k]
        return J_mat, C_mat, P_mat, w_mat


    def simulate(self, W0, t_death=None):
        if t_death==None:
            k_death, t_death = self.inv.sample_mortality()

        K = len(self.t_vec)
        W_vec = np.nan * np.ones(K)
        r_vec = np.nan * np.ones(K)
        α_vec = np.nan * np.ones(K)
        w_vec = np.nan * np.ones(K)
        C_vec = np.nan * np.ones(K)
        P_vec = np.nan * np.ones(K)
        U_vec = np.nan * np.ones(K)
        Z_vec = np.nan * np.ones(K)
        J_vec = np.nan * np.ones(K)
        Xbar_vec = np.nan * np.ones(K)
        Wbar_vec = np.nan * np.ones(K)
        Xbar_vec[0] = W0 + self.b_vec[0]
        W_vec[0] = W0

        for k,t in enumerate(self.t_vec[:k_death+1]):
            # Update controls:
            C_vec[k] = self.C(W_vec[k], k)
            P_vec[k] = self.P(W_vec[k], k)
            w_vec[k] = self.w(W_vec[k], k)

            # Update state:
            r_vec[k] = (1/self.Δt) * ((1+self.mkt.r)**self.Δt - 1)
            α_vec[k] = (1/self.Δt) * ((1+self.mkt.α)**self.Δt - 1) \
                       + np.sqrt(1/self.Δt)*self.mkt.σ*np.random.randn()
            W_vec[k+1] = W_vec[k] + self.Δt * (r_vec[k] * (1-w_vec[k])*W_vec[k] 
                          + α_vec[k]*w_vec[k]*W_vec[k]
                          + self.Y_vec[k] - P_vec[k] - C_vec[k])
            Z_vec[k] = P_vec[k] / self.μ_vec[k] + W_vec[k]

            # Update utility:
            U_vec[k] = self.inv.U(C_vec[k], t)
            J_vec[k] = self.J(W_vec[k], k)

        for k,t in enumerate(self.t_vec[:-1]):
            # Mean process
            r = (1/self.Δt) * ((1+self.mkt.r)**self.Δt - 1)
            α = (1/self.Δt) * ((1+self.mkt.α)**self.Δt - 1)
            σ = self.mkt.σ
            δ = self.inv.δ
            Xbar_vec[k+1] = Xbar_vec[k] + self.Δt * Xbar_vec[k] \
                            * ((α - r)**2/(δ*σ**2) + r + self.μ_vec[k]
                               - self.k_vec[k] / self.a_vec[k]**(1/δ))
            Wbar_vec[k] = Xbar_vec[k] - self.b_vec[k]
            
        total_utility = sum(self.Δt*U_vec[:k_death]) + self.inv.B(Z_vec[k_death],t)
        X_vec = W_vec + self.b_vec
        ARA_vec = self.inv.δ / X_vec
        dJdW_vec = self.a_vec * X_vec**(-self.inv.δ)
        dJdW_bar_vec = self.a_vec * Xbar_vec**(-self.inv.δ)
        discount = (dJdW_bar_vec[0] / dJdW_bar_vec[1])**(1/self.Δt) - 1

        results = {'t' : self.t_vec,
                   'W' : W_vec,
                   'X' : X_vec,
                   'w' : w_vec,
                   'r' : r_vec,
                   'α' : α_vec,
                   'Y' : self.Y_vec,
                   'P' : P_vec,
                   'C' : C_vec,
                   'J' : J_vec,
                   'π' : self.π_vec,
                   'λ' : self.λ_vec,
                   'G' : self.G_vec,
                   'U' : U_vec,
                   'Z' : Z_vec,
                   'b' : self.b_vec,
                   'μ' : self.μ_vec,
                   'X' : X_vec,
                   'Wbar' : Wbar_vec,
                   'Xbar' : Xbar_vec,
                   'ARA' : ARA_vec,
                   'dJdW': dJdW_vec,
                   'dJdW_bar': dJdW_bar_vec,
                   'discount': discount,
                   'total_utility' : total_utility,
                   't_death' : t_death,
                   'k_death' : k_death }
        return results
