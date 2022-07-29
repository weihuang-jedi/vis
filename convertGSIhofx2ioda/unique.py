#=========================================================================
import sys
import netCDF4 as nc4

#-----------------------------------------------------------------------------------------
debug = 1
dirname = 'sfc_ps_out'
#filename = 'sfc_ps_obs_2020011006_0000.nc4'
filename = 'sfcship_ps_obs_2020011006_0000.nc4'
#filename = 'sondes_ps_obs_2020011006_0000.nc4'

jedipath = '%s/%s' %(dirname, filename)

print('jedipath = ', jedipath)

nco = nc4.Dataset(jedipath, 'r')

blat = nco.groups['MetaData'].variables['latitude'][:]
blon = nco.groups['MetaData'].variables['longitude'][:]

print('len(blat) = ', len(blat))
#print('blat = ', blat)
#print('blon = ', blon)

nco.close()

dirname = 'ensmean'
filename = filename.replace('_0000.nc', '.nc')
gsipath = '%s/%s' %(dirname, filename)

print('gsipath = ', gsipath)

nco = nc4.Dataset(gsipath, 'r')

slat = nco.groups['MetaData'].variables['latitude'][:]
slon = nco.groups['MetaData'].variables['longitude'][:]

print('len(slat) = ', len(slat))
#print('slat = ', slat)
#print('slon = ', slon)

nco.close()

idx = []
blist = list(range(len(blat)))

for n in range(len(slat)):
  found = 0
  for i in range(len(blist)):
    k = blist[i]
    if((abs(slat[n] - blat[k]) < 0.01) and
       (abs(slon[n] - blon[k]) < 0.01)):
      blist.pop(i)
      idx.append(k)
      found = 1
      break
  if(not found):
    print('could not find slat[%d] = %f, slon[%d] = %f' %(n, slat[n], n, slon[n]))

print('len(idx) = ', len(idx))
#print('idx = ', idx)

#print('blat[idx] = ', blat[idx])
#print('blon[idx] = ', blon[idx])

