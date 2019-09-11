#!/usr/bin/python
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

gamma = 1; C1 = 10; zero = 1e-50

def solvr(X, t):
    return [X[1], -(gamma - X[3]*X[3])*X[1] - C1, X[3], -(2.0*X[1]/(X[0] + zero) - gamma)*X[3]]

t = np.arange(0, 0.4, 1e-2)
asol = integrate.odeint(solvr, [3, 1, 100, 1], t)


ax1 = plt.axes((0, 0, 0.5, 0.5))
ax1.plot(t, asol[:,0], color = 'k', marker = 'o', markersize = 3, linestyle = '-', label = r'$r$', markerfacecolor = 'None', linewidth = 1)
ax1.set_xlabel(r'$t$', fontsize = 12)
ax1.set_ylabel(r'$r(t)$', fontsize = 12)
ax1.tick_params(labelsize = 12)
#ax1.legend(loc = 'best', fontsize = 12)

ax2 = plt.axes((0, -0.6, 0.5, 0.5))
ax2.plot(t, asol[:,2], color = 'r', marker = 'o', markersize = 3, linestyle = '-', label = r'$\theta$', markerfacecolor = 'None', linewidth = 1)
ax2.set_xlabel(r'$t$', fontsize = 12)
ax2.set_ylabel(r'$\theta(t)$', fontsize = 12)
ax2.tick_params(labelsize = 12)
#ax2.legend(loc = 'best', fontsize = 12)

x = asol[:,0]*np.cos(asol[:,2]); y = asol[:,0]*np.sin(asol[:,2])
ax3 = plt.axes((0.65, -0.6, 1.0, 1.0 + 0.1))
ax3.plot(x, y, marker = 'o', markersize = 3, color = 'g', linestyle = '-', markerfacecolor = 'None', linewidth = 1)
ax3.set_xlabel(r'$x$', fontsize = 12)
ax3.set_ylabel(r'$y$', fontsize = 12)
ax3.tick_params(labelsize = 12)

plt.savefig('figure.pdf', dpi = 300, bbox_inches = 'tight')
