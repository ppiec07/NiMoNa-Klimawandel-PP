import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import cm


def Saltz(x_old, *function_parameters):
  x = x_old[0]
  y = x_old[1]
  z = x_old[2]
  xp = (-x) - x
  yp = p*z + r*y - (y*z*z)
  zp = q*(x - z)
  
  
  xnew = np.array([xp, yp, zp])
  return xnew

def RKF(func, x_old, stepsize, *function_parameters):
    alpha = np.array([0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0])
    beta = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0/4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0, 0.0],
        [1932.0/2197.0, (-7200.0)/2197.0, 7296.0/2197.0, 0.0, 0.0, 0.0],
        [439.0/216.0, -8.0, 3680.0/513.0, (-845.0)/4104.0, 0.0, 0.0],
        [(-8.0)/27.0, 2.0, (-3544.0)/2565.0, 1859.0/4104.0, (-11.0)/40.0, 0.0]])
    c = np.array([25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, (-1.0)/5.0, 0.0]) # coefficients for 4th order method
    c_star = np.array([16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, (-9.0)/50.0, 2.0/55.0]) # coefficients for 5th order method
    cerr=c-c_star # build the difference of both c and c_star for the error estimation
    """
    Berechnung der Konstanten nach dem Butcher-Tableau, um x_new und x_newstar zu schätzen
    anschließend Berechnung von x_new und x_newstar
    """
    
    k1 = func(x_old, function_parameters)
    k2 = func(x_old + beta[1][0]*k1, function_parameters)
    k3 = func(x_old + beta[2][0]*k1 + beta[2][1]*k2, function_parameters)
    k4 = func(x_old + beta[3][0]*k1 + beta[3][1]*k2 + beta[3][2]*k3, function_parameters)
    k5 = func(x_old + beta[4][0]*k1 + beta[4][1]*k2 + beta[4][2]*k3 + beta[4][3]*k4, function_parameters)
    k6 = func(x_old + beta[5][0]*k1 + beta[5][1]*k2 + beta[5][2]*k3 + beta[5][3]*k4 + beta[5][4]*k5, function_parameters)
    
    x_newstar = x_old + c_star[0]*k1 + c_star[1]*k2 + c_star[2]*k3 + c_star[3]*k4 + c_star[4]*k5 + c_star[5]*k6
    x_new = x_old + c[0]*k1 + c[1]*k2 + c[2]*k3 + c[3]*k4 + c[4]*k5 + c[5]*k6
    """
    fehler epsilon ist die Differenz aus x_new und x_newstar
    anschließend werden die beträge der werte sortiert und der höchste wert epsilon mit epsilon_tol verglichen
    Die Schrittweite wird entsprechend dem Ergebnis angepasst
    """
    
    epsilon = abs(x_new - x_newstar)
    epsilon = sorted(epsilon)
    
    if epsilon[2] < epsilon_tol :
        stepsize = 1*stepsize*((epsilon_tol/epsilon[2])**(0.25))
    else: 
        stepsize = 1*stepsize*((epsilon_tol/epsilon[2])**(0.2))
    """
    anschließend wird die neue Position x_newstar und die neue Schrittweite weitergegeben
    """
    
    return x_newstar, stepsize
stepsize = 0.025
a = 0.
b = 200.
epsilon_tol = 0.01 #maximal erlaubte Abweichung

p = 0.1
q = 1.5
r = 1.5
function_parameters = [p,r,q]

x0 = -np.sqrt(r-p)+0.0001
y0 = np.sqrt(r-p)
z0 = -np.sqrt(r-p)

tp = np.array([a]) #liste mit zeitpunkten
t = a*1 #startzeit bzw aktuelle zeit
result = np.array([x0, y0, z0]) # liste mit Lösungen
x_old = 1*result # Startpunkt

while t < b :
    """
    der Loop läuft bis zum erreichen des ende des Intervalls (b)
    als erstes wird die zeit aktualisiert und in die liste mit Zeitpunkten eingetragen
    danach wird ein neuer punkt sowie eine neue Schrittweite berechnet
    der neue punkt wird in die Lösungsliste eingetragen
    """
    t = t + stepsize
    tp += t
    x_old, stepsize = RKF(Saltz, x_old, stepsize, *function_parameters)
    result += x_old
    
    
    
x = [result[i][0] for i in range(len(tp))]
y = [result[i][1] for i in range(len(tp))]
z = [result[i][2] for i in range(len(tp))]

plt.plot(tp,z,label="Ocean Temp" , color = 'tab:brown')
plt.plot(tp,y,label="CO2" , color = 'tab:purple')
plt.plot(tp,x,label="Ice Mass", color = 'tab:green' )
plt.xlabel('$t$')
#plt.ylabel('$relation$')
plt.legend()
plt.title("Saltzmann and Sutera")
plt.grid()
plt.show()
