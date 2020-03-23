from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt

xl = np.array([-0.875, 0, 0.875])
yl = np.array([-0.875, 0, 0.875])
zl = np.array([0   , 0, 0])

fig = plt.figure()
ax = plt.axes(projection="3d")

def z_function(x, y, k, l):
    faik = 1.0
    for i in xl:
        if i==xl[k]:
            continue
        faik = faik*(x-i)/(xl[k]-i) 

    fail = 1.0
    for i in yl:
        if i==yl[l]:
            continue
        fail = fail*(y-i)/(yl[l]-i) 

    return faik*fail

x = np.linspace(-1, 1, 30)
y = np.linspace(-1, 1, 30)

X, Y = np.meshgrid(x, y)
Z = z_function(X, Y, 1, 2)

#fig = plt.figure()
#ax = plt.axes(projection="3d")
#ax.plot_wireframe(X, Y, Z, color='green')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')

ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap='winter', edgecolor='none', alpha=0.3)

ax.set_title('surface')
X2, Y2 = np.meshgrid(xl,yl)
ax.scatter(X2, Y2, 0)

plt.show()
