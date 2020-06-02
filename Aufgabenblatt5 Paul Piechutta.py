import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def rungekutta(func1, func2, func3, var, parameter, schrittweite):
    """
    Integriert 3 gekoppelte Dgls nach einem Runge-Kutta-Verfahren
    func1/func2/func3 : gekoppelte Dgls
    var : liste mit Startwerten
    parameter: Listen mit funktionsparametern
    schrittweite: Länge eines Integrationsschrittes


    """
    varhelp = [0,0,0]
    
    
    
    k1 = func1(var, parameter)[0]
    m1 = func2(var, parameter)[1]
    n1 = func3(var, parameter)[2]
    varhelp[0] = var[0] + 0.5*schrittweite*k1
    varhelp[1] = var[1] + 0.5*schrittweite*m1
    varhelp[2] = var[2] + 0.5*schrittweite*n1
    k2 = func1(varhelp, parameter)[0]
    m2 = func2(varhelp, parameter)[1]
    n2 = func3(varhelp, parameter)[2]
    varhelp[0] = var[0] + 0.5*schrittweite*k2
    varhelp[1] = var[1] + 0.5*schrittweite*m2
    varhelp[2] = var[2] + 0.5*schrittweite*n2
    k3 = func1(varhelp, parameter)[0]
    m3 = func2(varhelp, parameter)[1]
    n3 = func3(varhelp, parameter)[2]
    varhelp[0] = var[0] + schrittweite*k3
    varhelp[1] = var[1] + schrittweite*m3
    varhelp[2] = var[2] + schrittweite*n3
    k4 = func1(varhelp, parameter)[0]
    m4 = func2(varhelp, parameter)[1]
    n4 = func3(varhelp, parameter)[2]


    
    var[0] = var[0] + (schrittweite/6)*(k1 +2*k2 + 2*k3 + k4)
    var[1] = var[1] + (schrittweite/6)*(m1 +2*m2 + 2*m3 + m4)
    var[2] = var[2] + (schrittweite/6)*(n1 +2*n2 + 2*n3 + n4)
    return var

def xfunc(var,parameter):
   

    xneu = -var[0]-var[1]
    var[0] = xneu
    return var

def yfunc(var,parameter):
 
    yneu = parameter[0]*var[2] + parameter[2]*var[1] - var[1]*(var[2]**2)
    var[1] = yneu
    return var

def zfunc(var,parameter):

    zneu = parameter[1]*(var[0]-var[2])
    var[2] = zneu
    return var

p = 0.1 #Eingabe der Funktionsparameter und Startwerte, werden als liste gespeichert
q = 1.5
r = 1.5
parameter = [p,q,r]
x = -1*(math.sqrt(r-p)) + 0.0001 
y = math.sqrt(r-p)
z = -1*math.sqrt(r-p)  
var = [x,y,z]
h = 0.1 #schrittweite
ts = np.arange(0.0, 40, h) #liste mit Zeitpunkten
xs = [x] #listen der Funktionswerte
ys = [y]
zs = [z]
for t in ts[1:]:
    """
    berechnet neue Funktionswerte für jeden Schritt und hängt diese den listen an
    """
    var = rungekutta(xfunc, yfunc, zfunc, var, parameter, h)
    xs += [var[0]]
    ys += [var[1]]   
    zs += [var[2]]
    
  
plt.plot(ts, xs, label = "X") #plottet die ergebnisse als graphen
plt.plot(ts, ys, label = "Y")
plt.plot(ts, zs, label = "Z")
plt.xlabel("Zeit")
plt.ylabel("Funktionswerte")
plt.legend(loc = "upper right")
plt.title("Zeitreihe")
plt.show()

plt.axes(projection="3d") #plottet den phasenraum
plt.plot(xs, ys, zs)

plt.title("Phasenraum")

plt.show()