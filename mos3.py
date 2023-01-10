#    MIT License
#  Copyright (c) 2022 Mikhail Tegin
#  michail3110@gmail.com
#  
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.

import math
from scipy.optimize import fsolve, minimize, least_squares
import numpy as np

# Common parameters:
Esi = 1.035943139907e-10 # Permittivity of Silicon F/m
Eox = 3.45314379969e-11 # Permittivity of SiO2 F/m
q = 1.6E-19 # Electron charge
k = 1.3806226E-23 # Boltzmann constant
Vsb = 0 # Voltage substrate-source. It is fixed for optimization purpose.
Vtm = k*(273.15+20)/q # Temperature potential

# Class realizes NMOS LEVEL=3 model
class mos3():
    def __init__(self):
        # Parameters of model:
        self.VT0 = 0.0 # zero-bias threshold voltage
        self.KP = 0.0 # Transconductance parameter, A/V^2
        self.U0 = 600.0 # Surface mobility cm^2/(V*s)
        self.Weff = 100E-06 # in m
        self.Leff = 100E-06 # in m
        self.TOX = 1E-07 # Gate oxide thickness, m
        self.VMAX = 0.0 # Maximum carrier drift velocity, m/s
        self.THETA = 0.0 # mobility degradation factor, 1/V
        self.NSUB = 0.0 # Bulk surface doping cm-3
        self.XJ = 0.0 # Metallurgical junction depth, m
        self.PHI = 0.576 # Surface inversion potential, V
        self.DELTA = 0.0 # Narrow width factor for adjusting threshold
        self.GAMMA = 0.5276 # Body effect factor, sqrt(V)
        self.ETA = 0.0 # static feedback factor for adjusting threshold
        self.KAPPA = 0.2 # saturation field factor, 1/V
        self.NFS = 0.0 # fast surface state density cm^-2 V^-1
        self.update()
    # Update calculated parameters
    def update(self):
        self.ciNSUB = self.NSUB*1E+6
        self.ciNFS = self.NFS*1E+04
        self.ciU0 = self.U0*1E-04
        if self.ciNSUB != 0.0:
            self.Xd = math.sqrt(2.0*Esi/(q*self.ciNSUB)) # Depletion layer width, m
        else:
            self.Xd = 0.0
        self.COX = Eox/self.TOX # oxide capacitance per unit gate area, F/m^2
        self.fn = self.DELTA/self.Weff/4.0*2.0*math.pi*Esi/self.COX
        Wp = self.Xd*math.sqrt(self.PHI + Vsb)
        if self.XJ != 0:
            Wc = self.XJ*(0.0631353 + 0.8013292*Wp/self.XJ - 0.01110777*(Wp/self.XJ)**2)
            self.fs = 1.0 - self.XJ/self.Leff*(Wc/self.XJ*math.sqrt(1.0 - (Wp/(self.XJ + Wp)**2)))
        else:
            self.fs = 1.0
        self.fb = self.fn + self.GAMMA*self.fs/4.0/math.sqrt(self.PHI + Vsb) #  Начальный изгиб зон вблизи затвора, В
        self.fast = Vtm*(1.0 + q*self.ciNFS/self.COX + (self.GAMMA*self.fs*math.sqrt(self.PHI + Vsb) + self.fn*(self.PHI + Vsb))/2.0/(self.PHI + Vsb))
    # Returns channel length modulation, meters
    def __deltaL(self, Vds, Vdsat, Vc):
        if Vds <= Vdsat:
            return 0.0
        if self.VMAX <= 0.0:
            return self.Xd*math.sqrt(self.KAPPA*(Vds - Vdsat))
        else:
            Ep = self.KAPPA*Vc*(Vc + Vdsat)/self.Leff/Vdsat   # Multiplication on KAPPA is taken from ngspice (mos3load.c: 761)
            var = Ep*self.Xd**2/2.0
            return -var + math.sqrt(var**2 + self.KAPPA*self.Xd**2*(Vds - Vdsat))
    # Returns On Region current
    def __Id(self, Vgs, Vds, Vth):
        # Calculation of Vc
        us = self.ciU0/(1.0 + self.THETA*(Vgs - Vth))
        Vc = self.VMAX*self.Leff/us
        # Calculation of saturation voltage
        Vsat = (Vgs - Vth)/(1.0 + self.fb)
        if self.VMAX > 0.0:
            Vdsat = Vsat + Vc - math.sqrt(Vsat**2 + Vc**2)
        else:
            Vdsat = Vsat
        vde = min(Vds, Vdsat)
        # Calculation of effective mobility
        if self.VMAX > 0:
            Ueff = us/(1.0 + vde/Vc)
        else:
            Ueff = us
        if self.KP != 0.0:
            Ueff = Ueff/self.ciU0/self.COX*self.KP
        # Calculation of channel length modulation
        deltaL = self.__deltaL(Vds, Vdsat, Vc)
        if deltaL > self.Leff/2.0:
            deltaL = self.Leff - (self.Leff/2.0)**2/deltaL # HSPICE limit
        # Calculation of drain current
        return Ueff*self.COX*self.Weff/self.Leff*(Vgs - Vth - (1.0 + self.fb)/2.0*vde)*vde/(1.0 - deltaL/self.Leff)
    # Returns drain current
    def Id(self, Vgs, Vds):
        # Calculation of threshold voltage
        if self.VT0 == 0.0:
            Vbi = self.fb + self.PHI
        else:
            Vbi = self.VT0 - self.GAMMA*math.sqrt(self.PHI)
        Vth = Vbi - 8.14E-22*self.ETA/self.COX/self.Leff**3*Vds + self.GAMMA*self.fs*math.sqrt(self.PHI + Vsb) + self.fn*(self.PHI + Vsb)
        Von = Vth + self.fast
        if Vgs < Von:
            if self.ciNFS <= 0.0:
                return 0.0
            return self.__Id(Von, Vds, Vth)*math.exp((Vgs - Von)/self.fast)
        else:
            return self.__Id(max(Vgs, Von), Vds, Vth)
        
# Class realizes transistor that uses NMOS LEVEL=3
class Transistor(mos3):
    def __init__(self):
        super().__init__()
        self.Rd = 0.0 # Drain resistance
        self.Rs = 0.0 # Source resistance
        self.accuracy = 0.000001 # Maximum relative tolerance Vds
    def __error(self, Vds, Vgs, Vcc):
        I = (Vcc - Vds)/(self.Rd + self.Rs)
        return abs(super().Id(Vgs - I*self.Rs, Vds[0]) - I)
    def Id(self, Vgs, Vds):
        if self.Rs + self.Rd == 0.0:
                return super().Id(Vgs, Vds)
        Vcc = Vds
        Vds = fsolve(self.__error, [Vds], args = (Vgs, Vcc), xtol = self.accuracy)
        I = (Vcc - Vds[0])/(self.Rd + self.Rs)
        return I

# Class realizes optimization of parameters
# pointsSubTh ans pointsTh - tuple of points (Vgs, Vds, Id) that lies on experimental curves
# pointsSubTh - points on Cutoff Region
# pointsTh - points on On Region
class Optimizer():
    def __init__(self, model, pointsSubTh, pointsTh):
        self.model = model
        self.pointsSubTh = pointsSubTh
        self.pointsTh = pointsTh
        self.Id_vect = np.vectorize(model.Id)
        self.VT0min = 0.5
        self.VT0max = 10.0
        self.NFSmin = 1.0E+6
        self.NFSmax = 1.0E+16
        self.KPmin = 0.0
        self.KPmax = 10000.0
        self.THETAmin = 0.0
        self.THETAmax = 1000.0
        self.KAPPAmin = 0.0
        self.KAPPAmax = 200.0
    def __funResiduals(self, x):
        self.model.VT0 = x[0]
        self.model.NFS = x[1]
        self.model.KP = x[2]
        self.model.THETA = x[3]
        self.model.KAPPA = x[4]
        self.model.update()
        return [math.log(1.0 + abs(self.Id_vect(point[0], point[1]) - point[2])) for point in self.pointsSubTh] + [self.Id_vect(point[0], point[1]) - point[2] for point in self.pointsTh]
    def run(self):
        ret = True
        result = least_squares(self.__funResiduals, x0 = [self.model.VT0, self.model.NFS, self.model.KP, self.model.THETA, self.model.KAPPA],
                               bounds = [(self.VT0min, self.NFSmin, self.KPmin, self.THETAmin, self.KAPPAmin), (self.VT0max, self.NFSmax, self.KPmax, self.THETAmax, self.KAPPAmax)],
                               x_scale = [0.1, 100.0, 1.0, 0.1, 1.0])
        if not result.success:
            ret = False
        self.residuals = 0.0
        for r in self.__funResiduals(result.x):
            self.residuals = self.residuals + r*r
        return ret
