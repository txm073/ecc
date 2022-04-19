import numpy as np

def add(p, q, round=True, dp=5):
    if p is None:
        return q
    elif q is None:
        return p
    # Coefficients of elliptic curve
    curve = {"a": -7, "b": 10}
    # Check for vertical line
    if p[0] == q[0] and p[1] != q[1]:
        return (np.inf, np.inf)
    # If the points are the same use implicit differentiation
    if p == q:
        m = (3 * p[0] ** 2 + curve["a"]) / (2 * p[1])        
    # Calculate the line between two points (rise / run)
    else:        
        m = (p[1] - q[1]) / (p[0] - q[0])
    # Find co-ords of 3rd intersection point and reflect in the x-axis
    x = m ** 2 - p[0] - q[0]
    y = p[1] + m * (x - p[0])
    if round:
        x, y = np.round([x], dp)[0], np.round([y], dp)[0]
    return x, y * -1

def mult(p, n, round=True, dp=5):
    if n == 0:
        return (np.inf, np.inf)
    # Result of addition
    result = None
    # List to store the result of doubling previous result
    stored = [p]
    # Binomial expansion of n
    powers = [i for i, bit in enumerate(bin(n)[2:][::-1]) if bit == "1"]
    for i in range(1, powers[-1] + 1):
        stored.append(add(stored[-1], stored[-1], round=round, dp=dp))
    # Add powers of 2 that make up n
    for power in powers:
        result = add(result, stored[power], round=round, dp=dp)
    return result

# Curve: y^2 = x^3 - 7x + 10
print(mult((3, 4), 1000000, round=False))