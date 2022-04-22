"""Elliptic Curve Diffie-Hellman Encryption algorithm implemented in Python"""

import random


class EllipticCurve:

    params = {
        "name": "secp256k1",
        # `p` is a large prime number used as a modulo
        # To cause the points on the curve to cycle back around to prevent huge numbers
        "p": 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f, 
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
            return (None, None)
        # If the points are the same use implicit differentiation
        if p[0] == q[0]:
            m = (3 * p[0] * p[0] + self.params["a"]) * self.modinv(2 * p[1], self.params["p"])        
        # Calculate the line between two points (rise / run)
        else:        
            m = (p[1] - q[1]) * self.modinv(p[0] - q[0], self.params["p"]) 
        x = m ** 2 - p[0] - q[0]
        y = p[1] + m * (x - p[0])
        return x % self.params["p"], -y % self.params["p"] 

    def mult(self, p, n):
        """Multiply point `p` by `n` using double and add algorithm"""
        if n == 0:
            return (None, None)
        # Result of addition
        result = None
        point = p
        while n:
            if n & 1:
                # Add on the last iteratiion
                result = self.add(result, point)
            # Double
            point = self.add(point, point)
            # Bitwise right shift of n by 1 bit
            n >>= 1
        return result

    def modinv(self, a, m):
        """Calculate the multiplicative inverse of `a` under modulo `m`"""
        s, old_s = 0, 1
        t, old_t = 1, 0
        r, old_r = m, a
        while r != 0:
            quotient = old_r // r
            old_r, r = r, old_r - quotient * r
            old_s, s = s, old_s - quotient * s
            old_t, t = t, old_t - quotient * t
        gcd, x, y = old_r, old_s, old_t
        return x % m

    def derive_keys(self):
        """Derive a private/public keypair"""
        private = random.randrange(1, self.params["n"])
        public = self.mult((self.params["x"], self.params["y"]), private)
        return private, public

    def hcf(self, i, j):
        """Find the highest common factor (gcd) of two integers"""
        rem = i % j
        while rem:
            i = j
            j = rem
            rem = i % j
        return j


def main():
    import hashlib

    """
    An overview of the algorithm:

    - Both the client and server create a public/private key pair
      using the same elliptic curve and the same parameters: (p, a, b, G, n, h)
    - They then exchange their public keys (an x-y co-ordinate) and multiply the public key recieved
      together with their own private key (a random integer between 1 and 2^256)
    - This allows both the client and the server to establish a shared large number, which they can 
      combine with a hashing function (i.e. SHA-256) to create a fixed-length shared key
    - This shared key can be used to send information through a faster symmetric encryption 
      method, such as AES (This is how it is used in the TLS/SSL handshake)
    - The advantage to using Elliptic Curve Diffie-Hellman over traditional DH or RSA is that ECDHE
      is less computationally intensive, so you can obtain the same level of security as RSA with much
      shorter keys sizes, making encryption much faster
    """
    ec = EllipticCurve()
    a0, a1 = ec.derive_keys()
    b0, b1 = ec.derive_keys()
    s1 = ec.mult(a1, b0)
    s2 = ec.mult(b1, a0)
    assert s1 == s2
    secret = hashlib.sha256(str(s1[0]).encode()).hexdigest()
    print("Shared secret:", secret)

if __name__ == "__main__":
    main()