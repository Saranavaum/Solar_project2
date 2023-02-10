# Packages
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

# Plots parameters
plt.rcParams["figure.figsize"] = (10, 6)    # Size of the figures
plt.rc('xtick', labelsize = 10)             # Number size in X
plt.rc('ytick', labelsize = 10)             # Number size in Y
plt.rcParams['axes.titlesize'] = 15         # Title size
plt.rcParams['axes.labelsize'] = 10         # Label size on axes


#----Functions----#
# Least squares polynomial fit
def ajuste(I,landa,nlambda_l,nlambda_h):
    
    nl=nlambda_l
    nh=nlambda_h
        
    coef=np.polyfit(landa[nl:nh],I[nl:nh],deg=2)
    xadj=np.linspace(landa[nl],landa[nh],(nh-nl)*20) # The last number multiplied indicates how much we increase the resolution
    pol=np.poly1d(coef)
    Iadj= pol(xadj)
    
    return xadj, Iadj

# Search for the position of the minimum value
def min_arg(vec,in1,in2):
    vmin = min(vec[in1:in2])
    indmin=np.argmin(vec[in1:in2])+in1
    
    return vmin, indmin

# Search for the nearest point in an array
def find_closest(arr, val):  
    diff = abs(arr - val)
    closest = np.min(diff)
    index = np.where(diff == closest)[0][0]
    return index

def der(y, z):  # Derivada centered
    y1 = y[2:]
    y2 = y[:-2]
    z1 = z[2:]
    z2 = z[:-2]
    der = (y1-y2)/((z1-z2))
    return der


#---Main Program---#

#plt.close('all')

# Getting the data from the .fits
data=fits.getdata('Stokes_sunspot_HINODE.fits')
cal=fits.getdata('atlas_6301_6302.fits')

# Defining the Stokes parameters
I=data[0,:,:,:]
Q=data[1,:,:,:]
U=data[2,:,:,:]
V=data[3,:,:,:]

# Defining the calibration parameters
landa_cal=cal[:,0]
int_cal=cal[:,1]

# Enumeration of the frequencies to calibrate
landa=np.arange(len(I[0,0,:]))

# Cube dimensions: different positions of the slit in time in its sweep

# Plotting the sunspot
plt.figure(1)
mancha = plt.imshow(I[:, :, 0],origin='lower', cmap = 'inferno')
cbar = plt.colorbar(mancha)
cbar.set_label('I [cuentas]')
plt.title(' I (Mancha solar)')
plt.xlabel('Desplazamiento de la rendija en el tiempo [px]')
plt.ylabel('Puntos a lo largo de la rendija [px]')



# Selecting the calm section
xc=2
yc=2
dltc=50

# Selecting the umbra section
xu=160
yu=210
dltu=35


plt.gca().add_patch(Rectangle((xc,yc),dltc,dltc,
                    edgecolor='black',
                    facecolor='none',
                    lw=2))
plt.text(xc,yc+dltc+20, 'Calma', c = 'k', fontsize=15)

plt.gca().add_patch(Rectangle((160,210),35,35,
                    edgecolor='yellow',
                    facecolor='none',
                    lw=2))
plt.text(xu,yu+dltu+20, 'Umbra', c = 'y', fontsize=15)

#plt.savefig('1.eps',format='eps')

# Medium spectrum of I in calm section
I_calm=np.mean(I[xc:xc+dltc,yc:yc+dltc,:],axis=(0,1)) 
Q_calm=np.mean(Q[xc:xc+dltc,yc:yc+dltc,:],axis=(0,1)) 
U_calm=np.mean(U[xc:xc+dltc,yc:yc+dltc,:],axis=(0,1)) 
V_calm=np.mean(V[xc:xc+dltc,yc:yc+dltc,:],axis=(0,1)) 

plt.figure(2)
espect = plt.plot(landa,I_calm)
plt.plot(landa[22:27],I_calm[22:27],'r',markersize=1)
plt.plot(landa[68:74],I_calm[68:74],'r',markersize=1)
plt.title('Espectro I zona en calma (sin calibrar)')
plt.ylabel('Intensidad [cuentas]')
plt.xlabel(r'$\lambda$')
plt.grid()
#plt.savefig('2.eps',format='eps')


# Medium spectrum of I in umbra section
I_um=np.mean(I[xu:xu+dltu,yu:yu+dltu,:],axis=(0,1)) 
Q_um=np.mean(Q[xu:xu+dltu,yu:yu+dltu,:],axis=(0,1)) 
U_um=np.mean(U[xu:xu+dltu,yu:yu+dltu,:],axis=(0,1)) 
V_um=np.mean(V[xu:xu+dltu,yu:yu+dltu,:],axis=(0,1)) 

#plt.figure(3)
#plt.plot(landa,I_um)




#---Calibrating the Spectrum---#
plt.figure(4)
espect_cal = plt.plot(landa_cal,int_cal)
# the data in this atlas are multiplied by 1.e4 and have a sampling of 2 mA.
espect_cal = plt.plot(landa_cal[745:765],int_cal[745:765],'r',markersize=1)
espect_cal = plt.plot(landa_cal[1245:1260],int_cal[1245:1260],'r',markersize=1)
plt.title('Espectro de calibración')
plt.ylabel('Intensidad [cuentas]')
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.grid()
#plt.savefig('4.eps',format='eps')

# Search for the sunspot minimum

# Line 1 of the calibration
# search range
in1= 745
in2= 765

v_m_1 = min(int_cal[in1:in2])
ind_min1=np.argmin(int_cal[in1:in2])+in1
long_min=landa_cal[756]
plt.figure(5)
espect_cal = plt.plot(landa_cal[in1:in2],int_cal[in1:in2],'r.',label='datos')
espect_cal1 = plt.plot(landa_cal[756],int_cal[756],'k.',markersize=9,label='Minimo de los datos')


# Fitting the minimum
aju_atlas1=ajuste(int_cal,landa_cal,in1,in2)
aju_atlas1=np.array(aju_atlas1)
vmim1=min(aju_atlas1[1,:])
indmi1=np.where(aju_atlas1[1,:]==vmim1)
plt.plot(aju_atlas1[0,:],aju_atlas1[1,:],color='blue',linewidth=2,label='Ajuste de parábola')
plt.plot(aju_atlas1[0,indmi1],aju_atlas1[1,indmi1],'g.',markersize= 11,label='Minimo ajustado')
plt.title('Ajuste Linea 1 FTS')
plt.grid()
plt.legend()
plt.ylabel('Intensidad [cuentas]')
plt.xlabel(r'$\lambda$ [$\AA$]')
#plt.savefig('5.eps',format='eps')


# Line 2 of the calibration
in1= 1245
in2= 1260

plt.figure(6)
v_m_2 = min(int_cal[in1:in2])
ind_min2=np.argmin(int_cal[in1:in2])+in1
long_min=landa_cal[ind_min2]
plt.plot(landa_cal[in1:in2],int_cal[in1:in2],'r.',label='datos')
plt.plot(landa_cal[ind_min2],int_cal[ind_min2],'k.',markersize=9,label='Minimo de los datos')


#hacemos el ajuste
aju_atlas2c=ajuste(int_cal,landa_cal,in1,in2)
aju_atlas2c=np.array(aju_atlas2c)
vmim2=min(aju_atlas2c[1,:])
indmi2=np.where(aju_atlas2c[1,:]==vmim2)
plt.plot(aju_atlas2c[0,:],aju_atlas2c[1,:],color='blue',linewidth=2,label='Ajuste de parábola')
plt.plot(aju_atlas2c[0,indmi2],aju_atlas2c[1,indmi2],'g.',markersize= 11,label='Minimo ajustado')
plt.title('Ajuste Linea 2 FTS')
plt.grid()
plt.legend()
plt.ylabel('Intensidad [cuentas]')
plt.xlabel(r'$\lambda$ [$\AA$]')
#plt.savefig('6.eps',format='eps')

#linea 1 Modelo

in1= 22
in2= 27



plt.figure(7)
v_m_1m = min(I_calm[in1:in2])
ind_min1m=np.argmin(I_calm[in1:in2])+in1
long_min1m=landa[ind_min1m]
plt.plot(landa[in1:in2],I_calm[in1:in2],'r',label='datos')
plt.plot(landa[ind_min1m],I_calm[ind_min1m],'k.',markersize=9,label='Minimo de los datos')


#hacemos el ajuste
aju_atlas1m=ajuste(I_calm,landa,in1,in2)
aju_atlas1m=np.array(aju_atlas1m)
vmim1m=min(aju_atlas1m[1,:])
indmi1m=np.where(aju_atlas1m[1,:]==vmim1m)
plt.plot(aju_atlas1m[0,:],aju_atlas1m[1,:],'b',linewidth=2,label='Ajuste de parábola')
plt.plot(aju_atlas1m[0,indmi1m],aju_atlas1m[1,indmi1m],'g.',markersize= 11,label='Minimo ajustado')
plt.title('Ajuste Linea 1 HINODE')
plt.grid()
plt.legend()
plt.ylabel('Intensidad [cuentas]')
plt.xlabel(r'$\lambda$ [$\AA$]')
#plt.savefig('7.eps',format='eps')


#linea 2 Modelo

in1= 68
in2= 74



plt.figure(8)
v_m_2m = min(I_calm[in1:in2])
ind_min2m=np.argmin(I_calm[in1:in2])+in1
long_min2m=landa[ind_min2m]
plt.plot(landa[in1:in2],I_calm[in1:in2],'r',label='datos')
plt.plot(landa[ind_min2m],I_calm[ind_min2m],'k.',markersize=9,label='Minimo de los datos')


#hacemos el ajuste
aju_atlas2m=ajuste(I_calm,landa,in1,in2)
aju_atlas2m=np.array(aju_atlas2m)
vmim2m=min(aju_atlas2m[1,:])
indmi2m=np.where(aju_atlas2m[1,:]==vmim2m)
plt.plot(aju_atlas2m[0,:],aju_atlas2m[1,:],'b',linewidth=2,label='Ajuste de parábola')
plt.plot(aju_atlas2m[0,indmi2m],aju_atlas2m[1,indmi2m],'g.',markersize= 11,label='Minimo ajustado')
plt.title('Ajuste Linea 2 HINODE')
plt.grid()
plt.legend()
plt.ylabel('Intensidad [cuentas]')
plt.xlabel(r'$\lambda$ [$\AA$]')
#plt.savefig('8.eps',format='eps')


#Resultados del ajuste de los mínimos
m=[aju_atlas1m[0,indmi1m][0][0],aju_atlas1m[1,indmi1m][0][0]],[aju_atlas2m[0,indmi2m][0][0],aju_atlas2m[1,indmi2m][0][0]]
c=[aju_atlas1[0,indmi1][0][0],aju_atlas1[1,indmi1][0][0]],[aju_atlas2c[0,indmi2][0][0],aju_atlas2c[1,indmi2][0][0]]

#ecuacion de la recta
pend=(c[1][0]-c[0][0])/(m[1][0]-m[0][0])
orde=c[0][0]-pend*m[0][0]

long_cal=pend*landa+orde


plt.figure(9)
plt.plot(long_cal,I_calm,label='HINODE')
plt.plot(landa_cal,int_cal*2.95,label='FTS x2.95')
plt.grid()
plt.ylabel('Intensidad [cuentas]')
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.legend()
#plt.savefig('9.eps',format='eps')

plt.figure(10)
plt.plot(long_cal,I_um)
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel('I [cuentas]')
plt.grid()
plt.title('Espectro I zona umbra')
#plt.savefig('10.eps',format='eps')


########################Ensanchamiento doppler
IpV_um=I_um+V_um
ImV_um=I_um-V_um

in1= 61
in2= 69

aju_IpV=ajuste(IpV_um,long_cal,in1,in2)
aju_IpV=np.array(aju_IpV)
vmin_IV=min(aju_IpV[1,:])
vmax_IV=max(aju_IpV[1,:])
vmed=(vmax_IV-vmin_IV)/2+vmin_IV
vmed1=find_closest(aju_IpV[1,0:80],vmed)
vmed2=find_closest(aju_IpV[1,80:-1],vmed)+80

ind_mIpV=np.where(aju_IpV[1,:]==vmin_IV)

anc_dop=(aju_IpV[0,vmed2]-aju_IpV[0,vmed1])/2


plt.figure(11)
plt.plot(long_cal,IpV_um,'k.',label='I+V ')
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel('I [cuentas]')
plt.title('I+V zona umbra')
plt.grid()
#plt.plot(long_cal[in1:in2],IpV_um[in1:in2],'r',markersize=1)
plt.plot(aju_IpV[0,:],aju_IpV[1,:],label='ajuste')
plt.legend()
#plt.savefig('11.eps',format='eps')


##################Landa_b
plt.figure(12)
plt.plot(long_cal,IpV_um,'r-.',label='I + V')
plt.plot(long_cal,ImV_um,'b-.',label='I - V')
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel('I [cuentas]')

plt.grid()

plt.plot(aju_IpV[0,:],aju_IpV[1,:],'k',label='ajuste')


in1= 72
in2= 79

aju_ImV=ajuste(ImV_um,long_cal,in1,in2)
aju_ImV=np.array(aju_ImV)
vmin_ImV=min(aju_ImV[1,:])
ind_mImV=np.where(aju_ImV[1,:]==vmin_ImV)

plt.plot(aju_ImV[0,:],aju_ImV[1,:],'k')

anc_zeeman=(aju_ImV[0,ind_mImV]-aju_IpV[0,ind_mIpV])/2

plt.title('Ajuste Zeeman')
plt.legend()
#plt.savefig('12.eps',format='eps')

####################Campo magnético

C = 4.67e-13 
geff = np.array([1.667, 2.5])
lines = np.array([6301.52, 6302.51])
line1 = find_closest(aju_atlas2c[0,:], lines[0])
line2 = find_closest(aju_atlas2c[0,:], lines[1])

B_umb = anc_zeeman/(C*geff[1]*aju_atlas2c[0,line2]**2) 
################2377.47007631
########################Calculo de los ángulos


def angulos(Q,U,V):
    phi=0.5*np.arctan2(U,Q)
    C=1/(np.sqrt((Q**2+U**2))/V)
    gamma=np.arccos(np.divide(-1+np.sqrt(4*(C**2)+1),2*C))
    phi[np.isnan(phi)]=0
    
    return(phi*180/np.pi,gamma*180/np.pi)
    
long=64

phi,gamma=angulos(Q[:,:,long],U[:,:,long],V[:,:,long])


plt.figure(13)
#azimuth = plt.imshow(np.swapaxes(phi[:, :],0,1), cmap = 'rainbow')
azimuth = plt.imshow(phi[:, :],origin='lower' ,cmap = 'rainbow')
cbar = plt.colorbar(azimuth)
cbar.set_label(r'$\phi$ [º]')
plt.xlabel('Desplazamiento de la rendija en el tiempo [px]')
plt.ylabel('Puntos a lo largo de la rendija [px]')

#plt.savefig('13.eps',format='eps')



plt.figure(14)
gamma = plt.imshow(gamma[:, :],origin='lower', cmap = 'rainbow')
cbar = plt.colorbar(gamma)
cbar.set_label(r'$\gamma$ [º]')
plt.xlabel('Desplazamiento de la rendija en el tiempo [px]')
plt.ylabel('Puntos a lo largo de la rendija [px]')
#plt.savefig('14.eps',format='eps')


####################################Campo débil


########################Ensanchamiento doppler
IpV_calm=I_calm+V_calm
ImV_calm=I_calm-V_calm

in1= 20
in2= 28

aju_IpVc=ajuste(IpV_calm,long_cal,in1,in2)
aju_IpVc=np.array(aju_IpVc)
vmin_IVc=min(aju_IpVc[1,:])
vmax_IVc=max(aju_IpVc[1,:])
vmed=(vmax_IVc-vmin_IVc)/2+vmin_IVc
vmed1c=find_closest(aju_IpVc[1,0:80],vmed)
vmed2c=find_closest(aju_IpVc[1,80:-1],vmed)+80

ind_mIpVc=np.where(aju_IpVc[1,:]==vmin_IVc)

anc_dopc=(aju_IpVc[0,vmed2c]-aju_IpVc[0,vmed1c])/2



plt.figure(15)
plt.plot(long_cal,IpV_calm,'k.',label='I+V')
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel('I [cuentas]')
plt.title('I+V zona en Calma')
plt.grid()
#plt.plot(long_cal[in1:in2],IpV_calm[in1:in2],'r',markersize=1)
plt.plot(aju_IpVc[0,:],aju_IpVc[1,:],label='ajuste')
plt.legend()
#plt.savefig('15.eps',format='eps')

##################Landa_b
plt.figure(16)
plt.plot(long_cal,IpV_calm,'r-',label='I + V')
plt.plot(long_cal,ImV_calm,'b-.',label='I - V')
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel('I [cuentas]')
plt.title('Espectro I+-V Zona en calma')



plt.plot(aju_IpVc[0,:],aju_IpVc[1,:],'k')

in1= 20
in2= 28

aju_ImVc=ajuste(ImV_calm,long_cal,in1,in2)
aju_ImVc=np.array(aju_ImVc)
vmin_ImVc=min(aju_ImVc[1,:])
ind_mImVc=np.where(aju_ImVc[1,:]==vmin_ImVc)

plt.plot(aju_ImVc[0,:],aju_ImVc[1,:],'k',label='ajuste')

anc_zeemanc=(aju_ImVc[0,ind_mImVc]-aju_IpVc[0,ind_mIpVc])/2


plt.legend()
plt.grid()
#plt.savefig('16.eps',format='eps')

##################Campo  magnético longitudinal



paso=np.abs(long_cal[1]-long_cal[2])
dIdlan=np.gradient(I_calm,paso)
#dIdlan=der(I_calm, long_cal)
B_calm = -(1/(C*aju_atlas2c[0,line1]**2*geff[0]))*(np.sum(V_calm*dIdlan)/np.sum(dIdlan**2))
sigma = np.std(V_calm[40:60])
B_calm_err = sigma/(C*aju_atlas2c[0,line1]**2*geff[0]*np.sqrt(np.sum(dIdlan**2)))

print('Campo longitudinal, zona en calma',B_calm,'+-',B_calm_err)


##################################Temperatura

h = 6.63e-34    # J·s
k = 1.38e-23    # J/K
c = 3e8         # m/s
T_media = 6300  # K
ind_land=64

exp_min = np.exp(h*c/k*T_media*long_cal[ind_land]*1e-10)
Ic=28920

T=(1/T_media-(long_cal[ind_land]*1e-10*k)/(h*c)*np.log(I[:,:,ind_land]/Ic))**(-1)

plt.figure(17)
Temperatura = plt.imshow(T[:, :],origin='lower', cmap = 'inferno')
cbar = plt.colorbar(Temperatura)
cbar.set_label('T [K]')
plt.xlabel('Desplazamiento de la rendija en el tiempo [px]')
plt.ylabel('Puntos a lo largo de la rendija [px]')
#plt.savefig('17.eps',format='eps')

##################################Tmapa de velocidades


'''
bol_mat=np.zeros((np.shape(I)[0],np.shape(I)[1]))
for i in range(np.shape(I)[0]):
    for j in range(np.shape(I)[1]):
        if (max(V[i,j,:])<900):
            #bol_mat[i,j]=True
            I_min=min(I[i,j,:])
            ind_Imin=np.where(I[i,j,:]==I_min)
            l1=long_cal[ind_Imin]
            I_min_c=min(I_calm[40:60])
            ind_Imin2=np.where(I_calm[40:60]==I_min_c)
            l2=long_cal[ind_Imin2]
            bol_mat[i,j]=float(l1[0])-float(l2)

            
        else:
            V_min=min(V[i,j,:])
            V_max=max(V[i,j,:])
            V_mean=(V_min+V_max)/2
            ind_V=find_closest(V[i,j,:],V_mean)
            l1=long_cal[ind_V]
            I_min_c=min(I_calm[40:60])
            ind_Imin2=np.where(I_calm[40:60]==I_min_c)
            l2=long_cal[ind_Imin2]
            bol_mat[i,j]=float(l1)-float(l2)
'''
    
'''
fig=plt.figure(18)
ax=fig.add_subplot(111)
ax.imshow(bol_mat,aspect='auto')
'''

c = 3e8         # m/s
#v = (bol_mat*c/long_cal[line1])*1e-3  # km/s
v = np.loadtxt('velocity.txt')

v2=np.zeros(np.shape(v))

v2[100:300][30:150]=v[100:300][30:150]*1.00034

v3=v2+v

v=v-(v**(-1)*0.001)

plt.figure(18)
velmap = plt.imshow(v, cmap = 'seismic',origin='lower', vmin = -2, vmax = 2)
cbar = plt.colorbar(velmap)
cbar.set_label('v [km/s]')
plt.title('Mapa de velocidad ')
plt.xlabel('Desplazamiento de la rendija en el tiempo [px]')
plt.ylabel('Puntos a lo largo de la rendija [px]')
plt.savefig('18.eps',format='eps')
