import numpy as np
import matplotlib.pyplot as plt

def rungekutta(func1, func2, p1, p2, eps, gam, schrittweite):
    """
    Integriert 2 gekoppelte Dgls nach einem Runge-Kutta-Verfahren
    func1/func2 : gekoppelte Dgls
    p1/p2 : Startwerte von func1/func2
    eps, gam : Listen mit funktionsparametern für func1/func2
    schrittweite: Länge eines Integrationsschrittes
    
    die Funktion rechnet die Steigungen am Anfang(k1/m1), in der Mitte (k2,k3/m2/m3) und am Ende(k4/m4) des Integrationsintervalls aus, um daraus den neuen funktionswert zu errechnen
    die dafür benötigten funktionswerte(zB p1_fürk2) werden jeweils per Euler-verfahren geschätzt
    Dabei gehört k1, p1_fürk2, ... zu func1 und m1, p2_fürm2... zu func2


    """
    
    k1 = func1(p1, p2, eps[0], gam[0])
    m1 = func2(p1, p2, eps[1], gam[1])
    p1_fürk2 = p1 + 0.5*schrittweite*k1
    p2_fürm2 = p2 + 0.5*schrittweite*m1
    k2 = func1(p1_fürk2, p2_fürm2, eps[0], gam[0])
    m2 = func2(p1_fürk2, p2_fürm2, eps[1], gam[1])
    p1_fürk3 = p1 + 0.5*schrittweite*k2
    p2_fürm3 = p2 + 0.5*schrittweite*m2
    k3 = func1(p1_fürk3, p2_fürm3, eps[0], gam[0])
    m3 = func2(p1_fürk3, p2_fürm3, eps[1], gam[1])
    p1_fürk4 = p1 + schrittweite*k3
    p2_fürm4 = p2 + schrittweite*m3
    k4 = func1(p1_fürk4, p2_fürm4, eps[0], gam[0])
    m4 = func2(p1_fürk4, p2_fürm4, eps[1], gam[1])
    
    p1neu = p1 + (schrittweite/6)*(k1 +2*k2 + 2*k3 + k4)
    p2neu = p2 + (schrittweite/6)*(m1 +2*m2 + 2*m3 + m4)
    return p1neu, p2neu


def beute(p1, p2, eps1, gam1):
    """
    Differentialgleichung der Beute-Population
    p1: Beute-Population
    p2: Räuber-Population
    eps1: Konstante, beschreibt die Reproduktionsrate der Beute
    gam1: Konstante, beschreibt die Fressrate pro Beute
    """
    p1neu = p1*(eps1-gam1*p2)
    return p1neu

def räuber(p1, p2, eps2, gam2):
    """
    Differentialgleichung der Räuber-Population
    p1: Beute-Population
    p2: Räuber-Population
    eps2: Konstante, beschreibt die Sterberate der Räuber
    gam2: Konstante, beschreibt die Reproduktionsrate der Räuber pro Beutelebewesen
    """
    p2neu = -1*p2*(eps2-gam2*p1)
    return p2neu

p1 = 10000.0 #Eingabe der Funktionsparameter und Startwerte
p2 = 100.0
eps =[2.0, 0.8] # eps und gam als liste, da je zwei werte existieren
gam = [0.02, 0.0002]
h = 0.025 #schrittweite
ts = np.arange(0.0, 40, h) #liste mit Zeitpunkten
p1s = [p1] #listen der Funktionswerte
p2s = [p2]

for t in ts[1:]:
    """
    berechnet neue Populationswerte für jeden Schritt und hängt diese den listen an
    """
    p1, p2 = rungekutta(beute, räuber, p1, p2, eps, gam, h)
    p1s += [p1]
    p2s += [p2]
    
  
plt.plot(ts, p1s, label = "Beute") #plottet die ergebnisse als graphen
plt.plot(ts, p2s, label = "Räuber")
plt.xlabel("Zeit")
plt.ylabel("Anzahl Tiere")
plt.legend(loc = "upper right")
plt.title("Zeitreihe")
plt.show()

plt.plot(p2s, p1s) #plottet den phasenraum
plt.xlabel("Anzahl Räuber")
plt.ylabel("Anzahl Beute")
plt.title("Phasenraum")

plt.show()