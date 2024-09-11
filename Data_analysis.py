
import numpy as np
from astropy.io import fits

##conversion factors and constants in cgs
k=1.380658*1e-16
m_p=1.672623*1e-24
X_sol=0.71
Y_sol=0.265
Z_sol=0.025
gamma=4/3
mu_sol=1/(2*0.71+3/4*0.265+0.5*0.0125)
kpc=3.086*1e21
KeV=1.60218*1e-9
pi=3.1415926





###Radius with ellipses (in arcsec)

r_cav1_el=np.array([2.819, 2.435])*2*5
r_cav2_el=np.array([2.807, 2.166])*2*5





###diameter with ruler measurament (in arcsec)

d_cav1=np.array([5.0158, 5.32, 5.312, 5.5, 5.257, 4.899, 5.26, 4.955, 5.33, 5.37])*5#kpc
d_cav2=np.array([4.873, 4.93, 4.82, 4.53, 4.30, 4.60, 4.90, 4.74, 4.38, 4.71])*5#kpc

###Distance from center (in arcsec)
R_cav1=np.array([5.11, 4.98, 5.17, 4.87, 4.41])*5#kpc
R_cav2=np.array([4.61, 4.75, 4.91, 4.87,4.39])*5#kpc





#Mean radius of the cavity (assumed circular)
r_1=np.mean(d_cav1/2)
r_2=np.mean(d_cav2/2)

#Standard deviation of the cavity radius
r_1_err=np.std(d_cav1/2)
r_2_err=np.std(d_cav2/2)

#Mean distance of the cavities
R_1=np.mean(R_cav1)
R_2=np.mean(R_cav2)

#STD of cavity distance
R_1_err=np.std(R_cav1)
R_2_err=np.std(R_cav2)



##Assuming spherical geometry we calculate the volume of the cavities
V_1=4/3*pi*(r_1*kpc)**3
V_2=4/3*pi*(r_2*kpc)**3

V_1_err=4*pi*(r_1*kpc)**2*r_1_err*kpc
V_2_err=4*pi*(r_2*kpc)**2*r_2_err*kpc

print(r_2)





##We use annulus 4"-6"

kT=4.61  #kT 
kT_err=np.array([2.19, 1.06]) #kT


n=85.35*1e-3 #cm^-3 
n_err=np.array([7.03, 9.46])*1e-3 #cm^-3


Z=0.51 #Z_sol
Z_err=np.array([0.98, -0.51])  #Z_sol
Z=Z*Z_sol 
Z_max=(Z+Z_err[0])*Z_sol
Z_min=(Z+Z_err[1])*Z_sol
Z=np.array([Z, Z_max, Z_min])





#We normalize the abunace of H and He to the detected metallicity

C=X_sol/(1-X_sol-Z_sol)  ##this is a costant
X=C*(1-Z)/(C+1)
Y=(1-Z)/(C+1)

#We get mu
mu=1/(2*X+3/4*Y+0.5*Z)

mu_mean=mu[0]
mu_err=np.array(mu[1]-mu_mean, mu[2]-mu_mean)




#extrimate of pressure from x-rays observables
p_x=1.9*n*kT*KeV
p_x_err=np.mean(1.9*n*kT*KeV*np.sqrt((n_err/n)**2+(T_err/T)**2))




#Calculate the crossing time of the cavity

t_1=(R_1*kpc/(np.sqrt(gamma*kT*KeV/(mu_mean*m_p))))/(pi*1e7)
t_2=(R_2*kpc/(np.sqrt(gamma*kT*KeV/(mu_mean*m_p))))/(pi*1e7)

c=(gamma*KeV)/(mu_mean*m_p)

t_1_err=np.mean(np.sqrt((1/np.sqrt(c*kT)*R_1_err*kpc)**2+(c/(2*np.power(c*kT, 3/2))*T_err*KeV)**2)/(pi*1e7))
t_2_err=np.mean(np.sqrt((1/np.sqrt(c*kT)*R_2_err*kpc)**2+(c/(2*np.power(c*kT, 3/2))*T_err*KeV)**2)/(pi*1e7))

print("Cavity 1 has a crossing time of "+f"{t_1:.2e}"+ "±" + f"{t_1_err:.2e}"+" yr")
print("Cavity 2 has a crossing time of "+f"{t_2:.2e}"+ "±" + f"{t_2_err:.2e}"+" yr")




##Energy of cavities:

E_1=(gamma/(gamma-1))*p_x*V_1
E_2=(gamma/(gamma-1))*p_x*V_2

E_1_err =  E_1 * np.sqrt((V_1_err / V_1)**2 + (p_x_err / p_x)**2)
E_2_err =  E_2 * np.sqrt((V_2_err / V_2)**2 + (p_x_err / p_x)**2)


print("Cavity 1 has E=" +f"({E_1:.2e}" +" ± "+ f"{E_1_err:.2e})"+" erg")
print("Cavity 2 has E=" +f"({E_2:.2e}"+ " ± "+f"{E_2_err:.2e})"+" erg")





###Power of cavities

P1=E_1/(t_1*pi*1e7)
P2=E_2/(t_2*pi*1e7)

P1_err=P1*np.sqrt((E_1_err/E_1)**2+(t_1_err/t_1)**2)
P2_err=P2*np.sqrt((E_2_err/E_2)**2+(t_2_err/t_2)**2)


print("Cavity 1 has P=" + f"{P1:.2e} "+ "±" +f"{P1_err:.2e})"+" erg/s")
print("Cavity 2 has P=" + f"{P2:.2e} "+ "±" +f"{P2_err:.2e})"+" erg/s")







