#%%
import matplotlib.pyplot as plt

from matplotlib.patches import Ellipse

o_x   = 94
o_y   = 132
alpha = 53.758720
beta = 35
theta = -1.402590

fig, ax = plt.subplots()
ax.set(xlim=(0, 300), ylim=(0, 300), aspect="equal")
ellipse = Ellipse((o_x, o_y), alpha, beta, angle=theta, alpha=1)
ax.add_artist(ellipse)
ellipse.set_edgecolor('k')
ellipse.set_facecolor('w')

plt.show()
# %%
