

class Market():
    
    def __init__(self, r=.02, α=.07, σ=.17):
        if σ <= 0:
            raise ValueError('σ must be positive.')
        self.r = r
        self.α = α
        self.σ = σ

    def get_realization(self):
        return self.r, self.α + self.σ*np.random.randn()
