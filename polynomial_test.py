import matplotlib.pyplot as plt
import numpy as np

x = np.arange(-1, 1, 0.02)

p0 = 1
p1 = x
p2 = 0.5*(3*x**2-1)
p3 = 0.5*(5*x**3 -3*x)
p4 = (35*x**4 -30*x**2 +3)/8
p5 = (63*x**5 -70*x**3 +15*x)/8

p3_dot  = 0.5*(15*x**2-3)
p4_dot  = 0.5*x*(35*x**2-15)
p5_dot  = (315*x**4 -210*x**2 +15)/8
Gdg4    = (-1)**4*0.5*(p4 - p3)
Gdg5    = (-1)**5*0.5*(p5 - p4)
Gdg4_dot= (-1)**4*0.5*(p4_dot - p3_dot)
Gdg5_dot= (-1)**5*0.5*(p5_dot - p4_dot)

#plt.plot(x,p1)
#plt.plot(x,p2)
#plt.plot(x,p3)
plt.plot(x,p4)
#plt.plot(x,p5)
plt.plot(x,Gdg5)
#plt.plot(x,Gdg5_dot)
plt.show()

