import numpy as np
import matplotlib.pyplot as plt
import random

_min, _max = -10, 10
#a, b = random.randint(_min, _max), random.randint(_min, _max)
a, b = 27, 2    
validate_coefficients = lambda a, b: (4 * (a ** 3) + 27 * (b ** 2) != 0)
while not validate_coefficients(a, b):
    a, b = random.randint(_min, _max), random.randint(_min, _max)
f = lambda x: ((x ** 3) + (a * x) + b) 

print(f"Coefficients: a={a}, b={b}")
points = []
for x in np.linspace(-10, 10, 100):  
    y_squared = pow(x, 3) + a * x + b
    if y_squared <= 0:
        continue
    y = pow(y_squared, 0.5)
    points.append((x, y))
    points.append((x, y * -1))
x, y = zip(*points)

plt.figure()
plt.plot(x, y)
plt.show()      