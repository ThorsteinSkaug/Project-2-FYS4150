import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Importing from txt file made in the cpp code
df = pd.read_csv('iterations_per_N.txt', sep=' ', header=None)
itr = df[0].to_numpy()
N = df[1].to_numpy()

#Plotting
plt.plot(N,itr)
plt.grid()
plt.ylabel('Number of iterations')
plt.xlabel('Matrix size')
plt.title('Number of iterations needed as a function of matrix size')
plt.savefig('N_iterations.pdf', dpi=900)
plt.show()
