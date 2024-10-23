#%%
# Import Dependencies
from pyvlm.tools import airfoil_from_dat

#%%
# Create Airfoil
datfilepath = '../files/rhofw_root.dat'
airfoil = airfoil_from_dat(datfilepath)

#%%
# Plot Profile
ax = airfoil.plot_normalised_aifoil()
