import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Importing and plotting for both n=10 and n=100
fig,ax = plt.subplots(1,2,figsize=(10,5))
plt.suptitle(r'Vector elements as a function of position $\hat{x}$')
j = 0
for i in [9,99]:
    df = pd.read_csv(f'vec_val_{i}.txt', sep=' ', header=None)
    x = df[0].to_numpy()
    vec1 = df[1].to_numpy()
    vec2 = df[2].to_numpy()
    vec3 = df[3].to_numpy()

    # Since a scalar times an eigenvector is another eigenvector, we multiply the third eigenvector for
    # n=100 with -1, so that it has the same sign as for n=10
    if i == 99:
        vec3 = -vec3

    ax[j].plot(x,vec1, label=r'$\vec{v_1}$ for $\lambda_1$')
    ax[j].plot(x,vec2, label=r'$\vec{v_2}$ for $\lambda_2$')
    ax[j].plot(x,vec3, label=r'$\vec{v_3}$ for $\lambda_3$')
    ax[j].set_xlabel(r'Position $\hat{x}$'
    ax[j].set_ylabel('Vector elements')
    ax[j].set_title(f'n={i+1}')
    ax[j].legend()
    ax[j].grid()
    j += 1
plt.tight_layout()
plt.savefig('vecel_func.pdf', dpi=900)
plt.show()



# Importing and plotting both analytical and numerical solution for n=10 and n=100
fig,ax = plt.subplots(2,3,figsize=(10,5))
plt.suptitle(r'Vector elements as a function of position $\hat{x}$')
j = 0
for i in [9,99]:
    df = pd.read_csv(f'vec_val_{i}.txt', sep=' ', header=None)
    x = df[0].to_numpy()
    vec1 = df[1].to_numpy()
    vec2 = df[2].to_numpy()
    vec3 = df[3].to_numpy()

    # Since a scalar times an eigenvector is another eigenvector, we multiply the third eigenvector for
    # n=10 with -1, so that it has the same sign as the analytical solution
    if i == 9:
        vec3 = -vec3
    vec_list = [vec1,vec2,vec3]

    vec1_an = df[4].to_numpy()
    vec2_an = df[5].to_numpy()
    vec3_an = df[6].to_numpy()
    vec_list_an = [vec1_an,vec2_an,vec3_an]


    for k in range(len(vec_list)):
        ax[j][k].plot(x,vec_list[k], label=r'$\vec{v}_%s$ for $\lambda$' %(k+1))
        ax[j][k].plot(x,vec_list_an[k], label=r'Analytical $\vec{v}_%s$ for $\lambda$' %(k+1))
        ax[j][k].set_xlabel(r'Position $\hat{x}$')
        ax[j][k].set_ylabel('Vector elements')
        ax[j][k].set_title(f'n={i+1}')
        ax[j][k].legend()
        ax[j][k].grid()

    j += 1

plt.tight_layout()
plt.savefig('vec_vec.pdf', dpi=900)
plt.show()
