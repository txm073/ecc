#!/usr/bin/env python3
import collections
import random

EllipticCurve = collections.namedtuple('EllipticCurve', 'name p a b g n h')

curve = EllipticCurve(
    'secp256k1',
    # Field characteristic.
    p=0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f,
    # Curve coefficients.
    a=0,
    b=7,
    # Base point.
    g=(0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
       0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8),
    # Subgroup order.
    n=0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141,
    # Subgroup cofactor.
    h=1,
)


# Modular arithmetic ##########################################################

def inverse_mod(k, p):
    """Returns the inverse of k modulo p.
    This function returns the only integer x such that (x * k) % p == 1.
    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError('division by zero')

    if k < 0:
        # k ** -1 = p - (-k) ** -1  (mod p)
        return p - inverse_mod(-k, p)

    # Extended Euclidean algorithm.
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, k

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    gcd, x, y = old_r, old_s, old_t

    assert gcd == 1
    assert (k * x) % p == 1

    return x % p


# Functions that work on curve points #########################################

def is_on_curve(point):
    """Returns True if the given point lies on the elliptic curve."""
    if point is None:
        # None represents the point at infinity.
        return True

    x, y = point

    return (y * y - x * x * x - curve.a * x - curve.b) % curve.p == 0


def point_neg(point):
    """Returns -point."""
    assert is_on_curve(point)

    if point is None:
        # -0 = 0
        return None

    x, y = point
    result = (x, -y % curve.p)

    #assert is_on_curve(result)

    return result


def point_add(point1, point2):
    """Returns the result of point1 + point2 according to the group law."""
    #assert is_on_curve(point1)
    #assert is_on_curve(point2)

    if point1 is None:
        # 0 + point2 = point2
        return point2
    if point2 is None:
        # point1 + 0 = point1
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        # point1 + (-point1) = 0
        return None

    if x1 == x2:
        # This is the case point1 == point2.
        m = (3 * x1 * x1 + curve.a) * inverse_mod(2 * y1, curve.p)
    else:
        # This is the case point1 != point2.
        m = (y1 - y2) * inverse_mod(x1 - x2, curve.p)

    x3 = m * m - x1 - x2
    y3 = y1 + m * (x3 - x1)
    result = (x3 % curve.p,
              -y3 % curve.p)

    #assert is_on_curve(result)

    return result


def scalar_mult(k, point):
    """Returns k * point computed using the double and point_add algorithm."""
    #assert is_on_curve(point)

    if k % curve.n == 0 or point is None:
        return None

    if k < 0:
        # k * point = -k * (-point)
        return scalar_mult(-k, point_neg(point))

    result = None
    addend = point

    while k:
        if k & 1:
            # Add.
            result = point_add(result, addend)

        # Double.
        addend = point_add(addend, addend)

        k >>= 1

    #assert is_on_curve(result)

    return result


# Keypair generation and ECDHE ################################################

def make_keypair():
    """Generates a random private-public key pair."""
    private_key = random.randrange(1, curve.n)
    public_key = scalar_mult(private_key, curve.g)

    return private_key, public_key


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
ec = EllipticCurve()

print('Curve:', curve.name)

#g = ec.params["x"], ec.params["y"]
a0, a1 = ec.derive_keys()
print(ec.mult(a1, 1))
print(scalar_mult(1, a1))
exit(0)
#print("Alice's private key:", hex(alice_private_key))
#print("Alice's public key: (0x{:x}, 0x{:x})".format(*alice_public_key))

# Bob generates his own key pair.
bob_private_key, bob_public_key = make_keypair()
print("Bob's private key:", hex(bob_private_key))
print("Bob's public key: (0x{:x}, 0x{:x})".format(*bob_public_key))

# Alice and Bob exchange their public keys and calculate the shared secret.
s1 = scalar_mult(alice_private_key, bob_public_key)
s2 = scalar_mult(bob_private_key, alice_public_key)
print(s1)
print(s2)
