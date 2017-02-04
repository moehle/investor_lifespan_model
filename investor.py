import numpy as np

class Investor():

    def __init__(self, π, G, δ, h, m, Y):
        self.π = π
        self.G = G
        self.δ = δ
        self.γ = 1-δ
        self.h = h
        self.m = m
        self.Y = Y


    # Utility function:
    def U(self, C, t):
        if self.h(t) <= 0:
            return 0
        elif self.γ<1:
            return self.h(t) / self.γ * C**self.γ
        elif self.γ==0:
            return self.h(t)*np.log(C)
        else:
            raise ValueError("δ must be positive.")


    # Bequeathment utility function:
    def B(self, Z, t):
        if self.m(t) <= 0:
            return 0
        elif self.γ<1:
            return self.m(t) / self.γ * Z**self.γ
        elif self.γ==0:
            return self.m(t)*np.log(Z)
        else:
            raise ValueError("δ must be positive.")


    def sample_mortality(self):
        k_death = (sum(np.random.rand() < self.G_vec)) + 1
        return k_death, self.t_vec[k_death]



    def get_data_vectors(self, t0, tf, Δt):
        t_vec = np.arange(t0, tf, Δt)
        K = len(t_vec)
        π_vec = np.array([self.π(t) for t in t_vec])
        G_vec = np.array([self.G(t) for t in t_vec])
        h_vec = np.array([self.h(t) for t in t_vec])
        m_vec = np.array([self.m(t) for t in t_vec])
        Y_vec = np.array([self.Y(t) for t in t_vec])
        λ_vec = np.abs(np.array([π_vec[k] / G_vec[k] if G_vec[k]>0 else 1.0
                          for k in range(K)]))
        self.G_vec = G_vec
        self.t_vec = t_vec

        #m_vec /= π_vec
        return t_vec, π_vec, h_vec, m_vec, Y_vec, G_vec, λ_vec
