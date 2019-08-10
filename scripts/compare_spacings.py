#%% Import Dependencies
from pyvlm.tools import equal_spacing
from pyvlm.tools import semi_cosine_spacing, semi_elliptical_spacing
from pyvlm.tools import full_cosine_spacing, full_elliptical_spacing

#%% Create Spacings

num = 10

s1 = equal_spacing(num)
s2 = semi_cosine_spacing(num)
s3 = semi_elliptical_spacing(num)

num = 20

s4 = equal_spacing(num)
s5 = full_cosine_spacing(num)
s6 = full_elliptical_spacing(num)