import mmh3
import math

class HyperLogLog:
    def __init__(self, error_rate=0.01):
        self.b = self.determine_precision(error_rate)
        self.m = 1 << self.b
        self.M = [0] * self.m
        self.alpha_mm = self.alpha_m() * self.m * self.m

    def determine_precision(self, error_rate):
        m = (1.04 / error_rate) ** 2
        return math.ceil(math.log2(m))

    def alpha_m(self):
        if self.m == 16:
            return 0.673
        elif self.m == 32:
            return 0.697
        elif self.m == 64:
            return 0.709
        else:
            return 0.7213 / (1 + 1.079 / self.m)

    def hash(self, value):
        # Use MMH3 to hash the value to a 128-bit hash, then take the lower 64 bits
        value_bytes = str(value).encode('utf8')
        hash_value = mmh3.hash128(value_bytes)
        return hash_value & ((1 << 64) - 1)

    def add(self, value):
        x = self.hash(value)
        j = x >> (64 - self.b)
        w = x & ((1 << (64 - self.b)) - 1)
        self.M[j] = max(self.M[j], self.rho(w, 64 - self.b))

    def rho(self, w, max_width):
        return max_width - w.bit_length() + 1

    def estimate(self):
        Z = 1.0 / sum([2.0 ** -m for m in self.M])
        E = self.alpha_mm * Z

        # Small range correction
        if E <= 5.0 / 2.0 * self.m:
            V = self.M.count(0)
            if V > 0:
                E = self.m * math.log(self.m / V)
        
        # Large range correction
        elif E > 1.0 / 30.0 * (1 << 32):
            E = -(1 << 32) * math.log(1.0 - E / (1 << 32))

        return E

    def merge(self, other):
        if self.m != other.m:
            raise ValueError("Cannot merge HLLs with different precisions")
        for i in range(self.m):
            self.M[i] = max(self.M[i], other.M[i])

    def __str__(self):
        return f'HyperLogLog(b={self.b}, m={self.m})'