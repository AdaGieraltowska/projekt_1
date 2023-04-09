# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 22:10:05 2023

@author: Ada
"""
import numpy as np
class Transformation:
    
    def __init__(self,X,Y,Z,a,e2):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.a = a
        self.e2 = e2
        
    def Npu(self,fi):
        """
        I wertykal. W kierunku poziomym. 
        JEDNOSTKA -- metr
        ----------
        fi : float
            Szerokosc geodezyjna punktu. 
            JEDNOSTKA -- RAD
        a : dlugosc wielkiej polosi jednostka - metr optional
            DESCRIPTION. The default is 6378137.000.
        e2 : 1 mimosrod jednostka -- metr, optional
            DESCRIPTION. The default is 0.00669438002290.
    
        Returns
        -------
        N : float
            dlugosc w 1 wertykale. 
            JEDNOSTKA -- metr
    
        """
        N = self.a/np.sqrt(1-self.e2*(np.sin(fi))**2)
        return N 
    
    def hirvonen(self,X,Y,Z):

        p=np.sqrt(X**2+Y**2)
        fi = np.arctan(Z/(p*(1-self.e2)))
        
        # Krok 2 pętla 
        while True:
            N = self.Npu(fi)
            h = (p/np.cos(fi))-N
            fi_poprzednia = fi
            fi = np.arctan((Z/p)/((1-(N*self.e2/(N+h)))))
            if abs(fi_poprzednia-fi)<(0.000001/206265):
                break
            
        N = self.Npu(fi)
        h = p/np.cos(fi) - N
        lam= np.arctan(Y/X)
        
        #Kontrola
        Xk = (N+h)*np.cos(fi)*np.cos(lam)
        Yk = (N+h)*np.cos(fi)*np.sin(lam)
        Zk = (N*(1-self.e2)+h)*np.sin(fi)
        roznicax=X-Xk
        roznicay = Y - Yk
        roznicaz = Z - Zk 
        return fi,lam ,h
    
    
    
## SPRAWDZENIE CZY DZIAŁA

obiekt = Transformation(1,1,1,6378137, 0.00669438002290078)
result = obiekt.Npu(0.5)
hirek = obiekt.hirvonen(1,1,1)
print(result)