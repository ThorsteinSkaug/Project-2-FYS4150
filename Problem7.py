import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fig,ax = plt.subplots(1,2,figsize=(10,5))
plt.suptitle(r'Vector elements as a function of position $\hat{x}$')
j = 0
for i in [9,99]:
    df = pd.read_csv(f'vec_val_{i}.txt', sep=' ', header=None)
    x = df[0].to_numpy()
    vec1 = df[1].to_numpy()
    vec2 = df[2].to_numpy()
    vec3 = -df[3].to_numpy()
    if i == 9:
        vec3 = -vec3
    ax[j].plot(x,vec1, label=r'$\vec{v_1}$ for $\lambda_1$')
    ax[j].plot(x,vec2, label=r'$\vec{v_2}$ for $\lambda_2$')
    ax[j].plot(x,vec3, label=r'$\vec{v_3}$ for $\lambda_3$')
    ax[j].set_xlabel(r'Position $\hat{x}$')
    ax[j].set_ylabel('Vector elements')
    ax[j].set_title(f'n={i+1}')
    ax[j].legend()
    ax[j].grid()
    j += 1
plt.tight_layout()
plt.show()
