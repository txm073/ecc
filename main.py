import numpy as np
import random


class EllipticCurve:

    params = {
        "name": "secp256k1",
        # `p` is a large prime number used as a modulo
        # To cause the points on the curve to cycle back around to prevent huge numbers
        "p": 2 ** 256 - 2 ** 32 - 2 ** 9 - 2 ** 8 - 2 ** 7 - 2 ** 6 - 2 ** 4 - 1, 
        # `a` and `b` are the constants in the elliptic curve general formula: `y²=x³+ax+b`
        "a": 0, 
        "b": 7, 
        # `x` and `y` are the co-ordinates of the generator point
        "x": 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798, 
        "y": 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8, 
        # `n` is the number of operations required 
        # For n * the generator point to yield infinity, therefore defining
        # The maximum number of points that are on the curve and can be used
        "n": 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141,
        # The percentage of the curve that can be used for points a.k.a (p / n)
        "h": 1
    }

    def add(self, p, q):
        """Calculate the addition of two points on an elliptic curve"""
        if p is None:
            return q
        elif q is None:
            return p
        # Check for vertical line
        if p[0] == q[0] and p[1] != q[1]:
            return (np.inf, np.inf)
        # If the points are the same use implicit differentiation
        if p == q:
            m = ((3 * p[0] ** 2 + self.params["a"]) * self.modinv(2 * p[1], self.params["p"])) % self.params["p"]        
            x = (m ** 2 - 2 * p[0]) % self.params["p"]
            y = (m * (p[0] - x) - p[1]) % self.params["p"]
        # Calculate the line between two points (rise / run)
        else:        
            m = (q[1] - p[1]) * self.modinv(q[0] - p[0], self.params["p"]) % self.params["p"]
            x = (m ** 2 - p[0] - q[0]) % self.params["p"]
            y = (m * (p[0] - x) - p[1]) % self.params["p"]
        return x, y 

    def mult(self, p, n):
        """Multiply point `p` by `n` using double and add algorithm"""
        if n == 0:
            return (np.inf, np.inf)
        # Result of addition
        result = None
        # List to store the result of doubling previous result
        stored = [p]
        # Binomial expansion of n
        powers = [i for i, bit in enumerate(bin(n)[2:][::-1]) if bit == "1"]
        for i in range(1, powers[-1] + 1):
            stored.append(self.add(stored[-1], stored[-1]))
        # Add powers of 2 that make up n
        for power in powers:
            result = self.add(result, stored[power])
        return result

    def modinv(self, a, m):
        """Calculate the multiplicative inverse of `a` under modulo `m`"""
        mod = m
        x, y = 1, 0
        if (m == 1):
            return 0
        while (a > 1):
            q = a // m 
            t = m
            m = a % m
            a = t
            t = y
            y = x - q * y
            x = t
        return x % mod

    def derive_keys(self):
        private = random.randrange(1, self.params["n"])
        public = self.mult((self.params["x"], self.params["y"]), private)
        return private, public

# Curve: y^2 = x^3 - 7x + 10
curve = EllipticCurve()
a0, a1 = curve.derive_keys()
b0, b1 = curve.derive_keys()
print(a0)

exit()
alice_private_key, alice_public_key = curve.derive_keys()
print("Alice's private key:", hex(alice_private_key))
print("Alice's public key: (0x{:x}, 0x{:x})".format(*alice_public_key))

# Bob generates his own key pair.
bob_private_key, bob_public_key = curve.derive_keys()
print("Bob's private key:", hex(bob_private_key))
print("Bob's public key: (0x{:x}, 0x{:x})".format(*bob_public_key))
# Alice and Bob exchange their public keys and calculate the shared secret.
s1 = curve.mult(bob_public_key, alice_private_key)
s2 = curve.mult(alice_public_key, bob_private_key)
print(s1)
print(s2)