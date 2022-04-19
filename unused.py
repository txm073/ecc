from decimal import Decimal


def _round(n, threshold=5):
    if isinstance(n, int) or str(n).endswith(".0"):
        return n
    is_negative = 1 if n > 0 else -1
    string = str(n)
    integer = int(string[:string.index(".")])
    decimal = string[string.index("."):]
    if decimal.count("0") > threshold:
        point = decimal[:decimal.index("0" * threshold)]
        if point == ".":
            return (integer) * is_negative
        return (integer + float(decimal[:decimal.index("0" * threshold)])) * is_negative
    elif decimal.count("9") > threshold:
        point = decimal[:decimal.index("9" * threshold)]
        if point == ".":
            return (integer + 1) * is_negative
        return (integer + float(point[:-1] + str(int(point[-1]) + 1))) * is_negative
    return n

def add_points(p, q):
    # Coefficients of elliptic curve
    curve = {"a": -7, "b": 10}
    # Check for vertical line
    if p[0] == q[0] and p[1] != q[1]:
        return (np.inf, np.inf)
    # Starting point requires tangent
    if p == q:
        grad = (3 * p[0] ** 2 + curve["a"]) / (2 * p[1]) # Implicit differentiation
        line = {"grad": grad, "y_int": p[1] - (grad * p[0])}
    # Calculate the line between two points (rise / run)
    else:        
        line = {"grad": (p[1] - q[1]) / (p[0] - q[0])}
        line["y_int"] = p[1] - (line["grad"] * p[0])    
    # Substitute for y^2 and rearrange
    squared = {"a": line["grad"] ** 2, "b": line["grad"] * line["y_int"] * 2, "c": line["y_int"] ** 2}
    cubic = {
        "a": 1, 
        "b": squared["a"] * -1, 
        "c": curve["a"] - squared["b"],
        "d": curve["b"] - squared["c"]
    }        
    # Solve cubic and find the 3 intersection points
    roots = [_round(np.real(root)) for root in np.roots(list(cubic.values())) if np.isreal(root)]
    # If duplicates are present it means there are only 2 roots which will throw an error
    try:
        x = _round(list(filter(lambda x: x not in (p[0], q[0]), roots))[0])
    except IndexError:
        duplicate = max(set(roots), key=roots.count)
        x = _round(duplicate)
    y = abs((x ** 3 + curve["a"] * x + curve["b"]) ** (1 / 2)) * -1
    return x, y