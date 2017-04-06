import numpy as np

EPS = 0

class Insurer():
    
    def __init__(self, investor, μ=None):
        self.μ = μ
        self.inv = investor


    def get_data_vectors(self, t0, tf, Δt):
        t_vec, π_vec, h_vec, m_vec, Y_vec, G_vec, λ_vec = \
            self.inv.get_data_vectors(t0, tf, Δt)
        if self.μ:
            μ_vec = np.array([self.μ(t) for t in t_vec])
        else:
            μ_vec = λ_vec + EPS
            #μ_vec = λ_vec/10 # TODO
            #μ_vec = λ_vec/2 # TODO
            #μ_vec = λ_vec*2 # TODO
        η_vec = μ_vec - λ_vec
        H_vec = np.cumsum(η_vec[::-1])[::-1] * Δt
        return t_vec, μ_vec, η_vec, H_vec
