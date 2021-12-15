#!/usr/bin/env python
# coding: utf-8

# In[29]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from astropy.io import fits
from pynufft import NUFFT
import finufft


# In[30]:


# Defiying routes
antenna_route = '/home/seba/Desktop/alma.C34-2.cfg'
sky_image_route = '/home/seba/Downloads/cameraman(1).fits'
# Open cfg
antenna_pos = np.loadtxt(antenna_route, skiprows=3, usecols=(0, 1, 2))
# Open ffits
with fits.open(sky_image_route) as image:
    image_data = image[0].data
sky_image = image_data


# In[31]:


# Draw the entrys
fig,ax =  plt.subplots(1, 2, figsize=(10,10))
# Sky image
ax[0].set_title('Sky Image')
ax[0].set_xlabel('x(pixel)')
ax[0].set_ylabel('y(pixel)')
im_si = ax[0].imshow(sky_image, vmin=0, vmax=np.max(sky_image), aspect='equal')
# Antenna config.
ax[1].set_title('Antenna Config')
ax[1].set_xlabel('x(m)')
ax[1].set_ylabel('y(m)')
x1 = np.min(antenna_pos[:, 0])
x2 = np.max(antenna_pos[:, 0])
y1 = np.min(antenna_pos[:, 1])
y2 = np.max(antenna_pos[:, 1])
im_ac = ax[1].scatter(antenna_pos[:,0], antenna_pos[:,1], s=10, vmin=0, vmax=np.max(antenna_pos))
ax[1].set_xlim([x1 * 1.1, x2 * 1.1])
ax[1].set_ylim([y1 * 1.1, y2 * 1.1])
ax[1].set_aspect('equal', adjustable='box')
ax[1].axis([x1,x2,y1,y2])

fig.colorbar(im_si, ax=ax[0], orientation='horizontal')
fig.colorbar(im_ac, ax=ax[1], orientation='horizontal')
fig.tight_layout()
plt.show()


# In[32]:


# Define input variables of interferometer
obs_latitude = -23.0234 # latitude of observatory
source_declination = -60 # declination of the source
ha_start = -6 # in ha from -12 to 12
ha_end = 6 # in ha from -12 to 12
sample_interval = 600 # in seconds
frequency = 90e9 # in Hz equal to 90 GHz


# In[33]:


# antpos[] are XYZ coordinates in meters. XY aligned with celestial equator
# Z pointing to NCP
# antpos[1,:] -> East,  antpos[2, :] -> North, antpos[3,:] -> North Celestial Pole
# B : baselines vectors [:, 1] = x, [:, 2] = y, [:, 3] = z
# antpos[] antenna positions in meters [:, 1] = x, [:, 2] = y, [:, 3] = z

#
ant_set1 = np.repeat(antenna_pos, repeats=len(antenna_pos), axis=0)
ant_set2 = np.concatenate([antenna_pos] * len(antenna_pos), axis=0)
baseline = ant_set1 - ant_set2
# take off the [0,0,0] baseline
baseline = baseline[baseline[:, 0].argsort(kind='mergesort')]
valid_baselines = int(len(baseline) / 2) - int(len(antenna_pos) / 2)
baseline = np.concatenate((baseline[:valid_baselines, :], baseline[valid_baselines + len(antenna_pos):, :]), axis=0)
baseline = baseline.transpose()
print('the amount of baselines is ', len(baseline[0]))
print('the max baseline is ', np.max(abs(baseline[0])))


# In[34]:


#
# ha_rad[] vector with hour angle samplings in radians!!!
# ha_star, ha_end in Hours (for instance -1 hora, 1 hora)
# dt sampling interval in seconds (for instance 60 seconds)
# check hour angle limits
if ha_start < -12.0:
    ha_start = -12.0
if ha_end > 12.0:
    ha_end = 12.0

# convert from degrees to radian
ha_start_rad = 15 * ha_start * np.pi / 180.0
ha_end_rad = 15 * ha_end * np.pi / 180.0
# convert from seconds to radian
dt_rad = 15 * sample_interval * np.pi / (3600 * 180)
# generate array of hour angle in rads
ha_rad = np.arange(ha_start_rad, ha_end_rad, dt_rad)
print('the amount of ha obtained is ', len(ha_rad))


# In[35]:



# obs_freq is observing frequency in Hz (continuum ??)
# telescope_lat in degrees. Example -23
# source_decl Source declination in degrees. For instance 18

# transform from degree to radian
decl_rad = source_declination * np.pi / 180
lat_rad = (-90 + obs_latitude) * np.pi / 180

# Rotate around x to rise w -Dec degrees
R2 = np.array([[1, 0, 0],
               [0, np.cos((np.pi / 2) - decl_rad), -np.sin((np.pi / 2) - decl_rad)],
               [0, np.sin((np.pi / 2) - decl_rad), np.cos((np.pi / 2) - decl_rad)]])

# Rotate around x to correct for telescope latitude
R4 = np.array([[1, 0, 0],
               [0, np.cos(lat_rad), -np.sin(lat_rad)],
               [0, np.sin(lat_rad), np.cos(lat_rad)]])


# In[36]:


UVW = np.empty([3,3])
for k in range(0, len(ha_rad)):
    # Rotate around z to get w pointing to source
    # This matrix, R3, depends on HA, so it is computed for each
    # sampling
    # create a R3 for every HA
    R3 = np.array([[np.cos(ha_rad[k]), -np.sin(ha_rad[k]), 0],
                   [np.sin(ha_rad[k]), np.cos(ha_rad[k]), 0],
                   [0, 0, 1]])

    uvw = R4 @ R3 @ R2
    uvw = uvw.transpose()
    uvw = uvw @ baseline
    UVW = np.concatenate((UVW, uvw), axis=1)

# define lambda, being c the light speed constant
lambda_num = c / frequency

UVW = UVW/ lambda_num
max_uv = np.max(np.abs(UVW[:2]))
print('The number of visibilities get is ', len(UVW[0]))
print('The max uv coordinate is ', max_uv)


# In[37]:


# transform with fast fourier transform in 2 dimensions
fft_image = np.fft.fft2(sky_image)
fft_image = np.fft.fftshift(fft_image)

# transform with not uniform fast fourier transform in 2 dimensions
#       1st method with pynufft
NufftObj = NUFFT()                  #creation of the object
uv = UVW.transpose()                # trasposing to use in pynufft
uv = uv[:, :2]                      # only the u and v value excluding w value
uv = uv * np.pi / max_uv            # convert to radians
m, n = np.shape(sky_image)
Nd = (m, n)                         # image size
Kd = (1024, 1024)                   # k-space size
Jd = (10, 10)
NufftObj.plan(uv, Nd, Kd, Jd)
nufft_valuesv1 = NufftObj.forward(sky_image)

#       2nd method with finufft (the currently used)
u = UVW[0] * np.pi / max_uv
v = UVW[1] * np.pi / max_uv
complex_sky_image = sky_image + 1j * np.zeros_like(sky_image)
nufft_valuesv2 = finufft.nufft2d2(u, v, complex_sky_image, eps=1e-12)


# In[38]:


# inputs of the function gridder
imagesize = max(sky_image.shape)
scheme = 1
robust = 0
# the max uv used
epsilon = 1e-5
maxuv = max_uv + epsilon
# defines deltas
deltau = -2 * maxuv / (imagesize - 1)
deltav = - deltau
deltax = 1 / (imagesize * deltau)
deltay = 1 / (imagesize * deltav)

print('delta u is ', deltau, 'delta v is ', deltav)
print('delta x is ', deltax, 'delta y is ', deltay)


# In[39]:


# Finding the weight values
z = len(nufft_valuesv2)
uvw_weights = np.ones(z)
pos_u_index = np.zeros(z).astype(int)
pos_v_index = np.zeros(z).astype(int)
gridded_weights = np.zeros((imagesize, imagesize))
uv = uv * max_uv / np.pi

# Gridding weights if ROBUST OR UNIFORM
for k in range(0, z):
    j = int(np.round(uv[k][0] / deltau) + imagesize / 2)
    i = int(np.round(uv[k][1] / deltav) + imagesize / 2)
    pos_v_index[k] = i
    pos_u_index[k] = j
    gridded_weights[i][j] += uvw_weights[k]

# Selecting a scheme for weights
for k in range(0, z):
    if scheme == 0:
        uvw_weights[k] = uvw_weights[k]
    elif scheme == 1:
        uvw_weights[k] /= gridded_weights[pos_v_index[k]][pos_u_index[k]]
    elif scheme == 2:
        uvw_weights[k] *= np.sqrt(uv[k][0] ** 2 + uv[k][1] ** 2)
    else:
        average_weights = np.sum(gridded_weights ** 2) / np.sum(uvw_weights)
        f2 = (5.0 * np.power(10.0, -robust)) ** 2 / average_weights
        uvw_weights[k] /= (1.0 + gridded_weights[pos_v_index[k]][pos_u_index[k]] * f2)


# In[40]:


# Gridding visibilities and weights with a scheme
gridded_weights = np.zeros((imagesize, imagesize))
gridded_Vo = np.zeros((imagesize, imagesize)) + 1.0j * np.zeros((imagesize, imagesize))

for k in range(0, z):
    j = int(np.round(uv[k][0] / deltau) + imagesize / 2)
    i = int(np.round(uv[k][1] / deltav) + imagesize / 2)
    gridded_weights[i][j] += uvw_weights[k]
    gridded_Vo[i][j] += uvw_weights[k] * nufft_valuesv2[k]

# Check rows and columns where gridded weights are greater than 1
rows, columns = np.where(gridded_weights > 0)
# Divide the gridded visibilities by the gridded weights
gridded_Vo[rows, columns] /= gridded_weights[rows, columns]


# In[41]:


# get dirty image with the inverse fft in 2d
dirty_image = np.fft.ifftshift(gridded_Vo)
dirty_image = np.fft.ifft2(dirty_image)
dirty_image = np.fft.fftshift(dirty_image)


# In[45]:


# create the subplots
fig,ax =  plt.subplots(1, 4, figsize=(30,10))
# transform the complex values
fft_image = np.log(abs(fft_image) + 1)
nufft_valuesv1 = np.log(abs(nufft_valuesv1) + 1)
nufft_valuesv2 = np.log(abs(nufft_valuesv2) + 1)
gridded_Vo = np.log(abs(gridded_Vo) + 1)
dirty_image = abs(dirty_image)

# fast fourier transform image
ax[0].set_title('Transform Image')
ax[0].set_xlabel('u (λ)')
ax[0].set_ylabel('v (λ)')
limit_hor = deltax * len(fft_image[0]) / 2
limit_vert = deltay * len(fft_image) / 2
im_fft = ax[0].imshow(fft_image, aspect='equal', vmin=0, vmax=np.max(fft_image), extent=[-limit_hor, limit_hor, limit_vert, -limit_vert], origin='lower')
# not uniform fast fourier transform image
ax[1].set_title('Visibilities finufft')
ax[1].set_xlabel('u (λ)')
ax[1].set_ylabel('v (λ)')
im_nufft2 = ax[1].scatter(UVW[0] * deltau, UVW[1] * deltav, s=0.01, c=nufft_valuesv2, vmin=0, vmax=np.max(nufft_valuesv2))
ax[1].set_aspect('equal')
# gridded image
ax[2].set_title('Gridded Image')
ax[2].set_xlabel('u (λ)')
ax[2].set_ylabel('v (λ)')
limit_hor = deltau * max_uv
limit_vert = deltau * max_uv
im_gv = ax[2].imshow(gridded_Vo, aspect='equal', vmin=0, vmax=np.max(gridded_Vo), extent=[-limit_hor, limit_hor, -limit_vert, limit_vert], origin='lower')
# dirty image
ax[3].set_title('Dirty Image')
ax[3].set_xlabel('x (v)')
ax[3].set_ylabel('y (v)')
limit_hor = deltax * len(dirty_image[0]) / 2
limit_vert = deltay * len(dirty_image) / 2
im_di = ax[3].imshow(dirty_image, aspect='equal', vmin=0, vmax=np.max(dirty_image), extent=[-limit_hor, limit_hor, -limit_vert, limit_vert], origin='lower')
# colorbars
fig.colorbar(im_fft, ax=ax[0], orientation='horizontal')
fig.colorbar(im_nufft2, ax=ax[1], orientation='horizontal')
fig.colorbar(im_gv, ax=ax[2], orientation='horizontal')
fig.colorbar(im_di, ax=ax[3], orientation='horizontal')
fig.tight_layout()
plt.show()


# In[ ]:




