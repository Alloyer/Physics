from scipy import integrate
from math import pi, pow, sqrt
 
# func = lambda x: 1/cos(x)
func = lambda x: x / (sqrt(pow(x, 4) + 16))
 
answer = integrate.quad(func, 0, sqrt(3))
print(answer)