
import numpy as np
import math 
import matplotlib.pyplot as plt


def euler(func, x_old, function_parameters, stepsize):
    """
    Euler-Integration
    func : Funktion, die integriert wird
    x_old : aktueller Funktionswert(bei t_n)
    function_parameters : Liste mit Parametern der Funktion
    stepsize : Schrittweite der Euler-Integration

    Returns: z_new: Funktionswert zum Zeitpunkt t = t_n + h 
    
    """
    z_new = x_old + stepsize*func(x_old, *function_parameters)
    return z_new

def heun(func, x_old, z_new, function_parameters, newfunction_parameters, stepsize):
    x_new = x_old + 0.5*stepsize*(func(x_old, *function_parameters) + func(z_new, *newfunction_parameters))
    """
    Heun Integration
    da die euler integration der zweiten Ableitung (pendel_dgl1) den funktionswert x beinhält, muss dieser als parameter gespeichert werden, also nicht als x_old
    also wird das unnötig kompliziert(man braucht zwei parameterlisten), es wäre einfacher gewesen einfach von vornherein 2 eulerfunktionen zu haben
    """
    return x_new

def pendel_dgl1(x_old, omega, x ):
    """
    Erster Schritt der Integration der Schwingungsgleichung, Integration liefert erste zeitlich Ableitung der Schwingungsgeichung
    x_old : aktueller Funktionswert(bei t_n)
    omega: Kreisfrequenz
    """
    
    return (-1*(omega**2)*x)

def pendel_dgl2(x_old, x_punkt):
    """
    
    Zweiter Schritt der Integration der Schwingungsgleichung
    Integration der ersten zeitlichen Ableitung liefert die Schwingungsgleichung
    x_old: aktueller Wert der Schwingungsgleichung, fließt nicht in die Änderung mit rein
    x_punkt: aktueller Wert der ersten Ableitung

    """
    
    return 1*x_punkt

#Aufrufen der Funktion:
x = 0.0 #Startwert
x_punkt = 1.0 #Startwert der ersten Ableitung
omega = 1.0 #Kreisfrequenz
h = 0.1 #Schrittweite der Integration
ts = np.arange(0., 2*math.pi, h) 
x_punkts = [x_punkt] # Liste mit den x-Werten der ersten Ableitung, wird gefüllt
xs = [x] # Liste mit den x-Werten der Schwingungsgleichung , wird gefüllt

#plt.plot(ts[0], xs[-1], 'ro') #plottet den ersten punkt, wäre zu unübersichtlich bei kleinem h

for t in ts[1:]:
    z_punkt = euler(pendel_dgl1, x_punkt, [omega, x], h) #Eulerintegration zum Abschätzen des Wertes der Ableitung bei t_n = t_n +h
    z = euler(pendel_dgl2, x, [x_punkt], h) #Eulerintegration zum Abschätzen des Wertes der Schwingungsfunktion bei t_n = t_n +h
    x_punkt = heun(pendel_dgl1, x_punkt, z_punkt, [omega, x], [omega, z], h) #Berechnung des Wertes der Ableitung bei t = t_n + h
    x_punkts += [x_punkt]   #listet die x-Werte der ersten Ableitung auf 
    x = heun(pendel_dgl2, x, z, [x_punkt], [z_punkt], h) #Berechnung des Wertes der Schwingungsgleichung bei t = t_n + h
    xs += [x] #listet die x-Werte der Schwingungsgleichung auf 
    
    #plt.plot(t, xs[-1], 'ro') #plottet den neu errechneten punkt, also bei t_n+1, wäre zu unübersichtlich bei kleinem h
    
    
plt.plot(ts, xs) #plottet die ergebnisse als graph
plt.show()
