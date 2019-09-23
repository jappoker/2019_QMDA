import numpy as np
import matplotlib.pyplot as plt

C = (0,0)
B = (6,0)
D = (np.cos(25/180*np.pi) * 10 , np.sin(25/180*np.pi) * 10 )


plt.scatter(*C)
plt.scatter(*B)
plt.scatter(*D)
plt.show()
