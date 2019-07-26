#!/usr/bin/python3
# 2-CLAUSE BSD LICENCE
#Copyright 2010-2018 Hugo Tabernero, Jonay Gonzalez Hernandez, David Montes, and Emilio Gomez Marfil 
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE

from scipy.interpolate import griddata as intrp
import numpy as np
import _pickle as pic
#import matplotlib.pyplot as plt
from time import time
from os import system
from scipy import optimize

def interpol_mod(x, metal,fesun):
    Teff = np.round(x[0],0)
    logg = np.round(x[1],2)
    micro = np.round(x[2],2)
    metal = np.round(metal,2)
    sel=np.where((np.abs(tmod-Teff) <= 251.) & (np.abs(mmod-metal) <= 0.251) & (np.abs(gmod-logg) <= 0.5))
    lT=np.shape(tmod[sel])[0]
    lm=np.shape(mmod[sel])[0]
    lg=np.shape(gmod[sel])[0]
    if lT > 1 and lm > 1 and lg > 1 and micro <= 3. and micro >= 0.:
        dT=max(tmod[sel])-min(tmod[sel])
        dm=max(mmod[sel])-min(mmod[sel])
        dg=max(gmod[sel])-min(gmod[sel])
        testt= (min(tmod[sel]) <= Teff) and (max(tmod[sel]) >= Teff) and dT > 0.
        testm= (min(mmod[sel]) <= metal) and (max(mmod[sel]) >= metal) and dm > 0.
        testg= (min(gmod[sel]) <= logg) and (max(gmod[sel]) >= logg) and dg > 0.
        if testt and testm and testg:
            testint=True
        else:
            testint=False
    else:
        testint=False
    if testint:
        teint=intrp((tmod[sel],gmod[sel],mmod[sel]),Temod[sel],(Teff,logg,metal),method='linear',fill_value=np.nan, rescale=True)
        if np.isfinite(teint[0]):
            lpgint=intrp((tmod[sel],gmod[sel],mmod[sel]),lpgmod[sel],(Teff,logg,metal),method='linear',fill_value=np.nan, rescale=True)
            lpeint=intrp((tmod[sel],gmod[sel],mmod[sel]),lpemod[sel],(Teff,logg,metal),method='linear',fill_value=np.nan, rescale=True)
            rhoxint=intrp((tmod[sel],gmod[sel],mmod[sel]),rhoxmod[sel],(Teff,logg,metal),method='linear',fill_value=np.nan, rescale=True)
            kint=intrp((tmod[sel],gmod[sel],mmod[sel]),kmod[sel],(Teff,logg,metal),method='linear',fill_value=np.nan,rescale=True)    
            nlayers=len(kint)
            model=open('out.mod','w')
            model.write('KURUCZ\n')
            model.write('MODEL MOOG Teff = {0:5.0f} log g = {1:4.2f}\n'.format(Teff,logg))
            model.write('NTAU        {0:3.0f}\n'.format(nlayers))
            for i in range(nlayers):
                model.write(' {0:15.8E}  {1:7.1f} {2:10.3E}  {3:10.3E} {4:10.3E}\n'.format(rhoxint[i],teint[i],10.**lpgint[i],10.**lpeint[i],kint[i]))
            model.write('    {0:5.3f}e+05\n'.format(micro))
            model.write('NATOMS     1  {0:5.2f}\n'.format(metal))
            model.write('      26.0   {0:4.2f}\n'.format(metal+fesun))
            model.write('NMOL      19\n')
            model.write('      606.0    106.0    607.0    608.0    107.0    108.0    112.0    707.0\n')
            model.write('      708.0    808.0     12.1  60808.0  10108.0    101.0      6.1      7.1\n')
            model.write('        8.1    822.0     22.1\n')
            model.close()
            return True
        else:
            return False
    else:
        return False

def moog_run(starname):
    modelname = 'out.mod'
    moog_par_file  = 'PAR/'  + starname + '_abfind.par'
    moog_txt_file  = 'TXT/'  + starname + '_abfind.rsp'
    moog_out2_file = 'ABUN/' + starname + '_Feout2'
    # Running MOOG.
    system('./MOOGnointro < ' + 'TXT/' + starname + '_abfind.rsp  > MOOG_last_call')
    # These lines delete MOOG messages (OH NO...) using sed, so MOOG output can be properly read.
    system('sed \'/OH/d\' ' + moog_out2_file + ' > MOOG_cleanFeout2')
    system('mv MOOG_cleanFeout2 ' + moog_out2_file)

def abun_read(fich):
    abfile=open(fich,'r')
    abuns=abfile.readlines()
    abfile.close()
    abunds=[]
    for i in range(len(abuns)):
        temp=abuns[i].split()
        if len(temp) == 8:
                if temp[0] != 'MINE' and temp[0] != 'wavelength':
                    abunds.append([float(j) for j in temp])
    return abunds

def moog_read(starname):
    moog_out2_file = 'ABUN/' + starname + '_Feout2'
    out2_data=abun_read(moog_out2_file)
    nlines=len(out2_data)
    ab  = np.array([out2_data[i][6] for i in range(nlines)]) 
    ele = np.array([out2_data[i][1] for i in range(nlines)])
    EP  = np.array([out2_data[i][2] for i in range(nlines)])
    RW  = np.array([out2_data[i][5] for i in range(nlines)])
    sel1 = ele < 26.05
    sel2 = ele > 26.05
    return ab[sel1],EP[sel1],RW[sel1],ab[sel2]

def moog_eval(starname,fesun):
    moog_run(starname)
    ab1,EP1,RW1,ab2 = moog_read(starname) 
    fit1  = np.polyfit(EP1, ab1, 1)
    fit2  = np.polyfit(RW1, ab1, 1)
    slope1=np.round(fit1[0],4)
    slope2=np.round(fit2[0],4)
    mab1=np.round(np.mean(ab1)-fesun,3)
    mab2=np.round(np.mean(ab2)-fesun,3)
    return slope1,slope2,mab1,mab2

def objective_function_vec(x,met,starname,fesun):
    evalmod = interpol_mod(x, met,fesun)
    if evalmod:
        return moog_eval(starname,fesun)
    else:
        return 10.**20.,10.**20.,10.**20.,10.**20. 

def objective_function(x,met,starname,fesun):
    evalmod = interpol_mod(x, met,fesun)
    if evalmod:
        slp1, slp2, AFe1, AFe2 = moog_eval(starname,fesun)
        #function given by Nuno Santos (priv comm)
        return 5*((3.5* slp1)**2.+(1.3*slp2)**2.)+2*(AFe1 - AFe2)**2.

    else: 
        return 10.**20.
# THE NELDER-MEAD -> PRESS ET AL. -> NUMERICAL RECIPES 
def simplex(S, met, starname, fesun):
    Xm = [0, 0, 0]
    Xr = [0, 0, 0]
    Xe = [0, 0, 0]
    Xc = [0, 0, 0]

    for i in range(0, 3):
        for j in range(0, 3):
            Xm[i] += S[j][1][i]
        Xm[i] /= 3
        Xr[i] = 2 * Xm[i] - S[3][1][i]

    fr = objective_function(Xr, met, starname, fesun)

    if S[0][0] <= fr and S[2][0] > fr:
        S[3][1] = Xr
        S[3][0] = fr

    elif fr < S[0][0]:

        for i in range(0, 3):
            Xe[i] = (3 * Xm[i]) - (2 * S[3][1][i])
        fe = objective_function(Xe, met, starname, fesun)

        if fe < fr:
            S[3][1] = Xe
            S[3][0] = fe
        else:
            S[3][1] = Xr
            S[3][0] = fr

    else:

        for i in range(0, 3):
            Xc[i] = 0.5 * (Xm[i] + S[3][1][i])
        fc = objective_function(Xc, met, starname, fesun)

        if fc <= S[3][0]:
            S[3][1] = Xc
            S[3][0] = fc
        else:
            for i in range(1, 4):
                for j in range(0, 3):
                    S[i][1][j] = 0.5 * (S[0][1][j] + S[i][1][j])
                S[i][0] = objective_function(S[i][1], met, starname, fesun)
    S.sort(key=lambda x: x[0])
    return S
#outlier razorblade
def outlier_razor(starname,x,met,fesun):
        evalmod = interpol_mod(x,met,fesun)
        if evalmod:
            moog_run(starname)
            ab1,EP1,RW1,ab2 = moog_read(starname)
            n1=len(ab1)
            n2=len(ab2)
            mab1 = np.median(ab1)
            mab2 = np.median(ab2)
            sab1=np.median(np.abs(ab1-mab1))
            sab2=np.median(np.abs(ab2-mab2))
            wFe1=np.ones(n1)
            wFe2=np.ones(n2)
            wFe1[np.abs(ab1 - mab1) > 3.*sab1] = 0.
            wFe2[np.abs(ab2 - mab2) > 3.*sab2] = 0.
            ewfile=open('EW/' + starname + 'Fe.l','r')
            ewlist = ewfile.readlines()
            ewfile.close()
            ewfilenew=open('EW/'+ starname + 'newFe.l', 'w')
            ewfilenew.write(ewlist[0])
            n1=len(ab1)
            n2=len(ab2)
            for i in range(n1):
                if wFe1[i] > 0.:
                    ewfilenew.write(ewlist[i + 1])
            for i in range(n2):
                if wFe2[i] > 0.:
                    ewfilenew.write(ewlist[i + 1 + n1])
            ewfilenew.close()
            return True
        else:    
            return False
# Slope errors, etc, etc.

def get_err_slopes(x,met,starname,fesun):
    evalmod = interpol_mod(x,met,fesun)
    if evalmod:
        moog_run(starname)
        ab1,EP1,RW1,ab2 = moog_read(starname)
        n1=len(ab1)
        n2=len(ab2)
        slope1, orden1 = np.polyfit(EP1, ab1, 1)
        slope2, orden2 = np.polyfit(RW1, ab1, 1)
        fact=(n1-1.)*(n1-2.)
        serr1   = np.sum((ab1-(slope1*EP1)-orden1)**2.)
        serr2   = np.sum((ab1-(slope2*RW1)-orden2)**2.)
        sslope1 = np.sqrt(serr1/(fact*np.var(EP1)))
        sslope2 = np.sqrt(serr2/(fact*np.var(RW1)))
        
        AFe1  = np.round(np.mean(ab1)-fesun,3)
        AFe2  = np.round(np.mean(ab2)-fesun,3)
        sAFe1 = np.round(np.std(ab1)/np.sqrt(n1),3)
        sAFe2 = np.round(np.std(ab2)/np.sqrt(n2),3)
        
        return slope1,slope2,AFe1,AFe2,sslope1,sslope2,sAFe1,sAFe2
    else:
        return -9999.9,0.,0.,0.,0.,0.,0.

# THE ERROR FUNCTIONS (aka, THE KRAKEN)

def error_v(vmic, T, logg, met, starname,  fesun, dslope2):
    x = [T, logg, vmic]
    sl1, sl2, a1, a2 = objective_function_vec(x, met, starname,  fesun)
    return sl2 - dslope2


def error_tm(T, logg, vmic, met, starname,  fesun, dslope1):
    x = [T, logg, vmic]
    sl1, sl2, a1, a2 = objective_function_vec(x, met, starname,  fesun)
    return sl1 - dslope1


def error_tv(T, logg, vmic, met, starname,  fesun):
    x = [T, logg, vmic]
    sl1, sl2, a1, a2 = objective_function_vec(x, met, starname,  fesun)
    return sl1


def error_gd(logg, T, vmic, met, starname,  fesun, dAFe2):
    x = [T, logg, vmic]
    sl1, sl2, a1, a2 = objective_function_vec(x, met, starname,  fesun)
    return a2 - dAFe2


def error_gp(logg, T, vmic, met, starname,  fesun):
    x = [T, logg, vmic]
    sl1, sl2, a1 ,a2 = objective_function_vec(x, met, starname,  fesun)
    return a1 - a2

# geterrors
def check_error_v(x,met,starname,fesun,dslope2):
      T, logg, vmic =  x
      delta_vmic    =  0.1
      fa, fb = [1.,1.]
      n = 1
      root = False
      while fa*fb > 0. and n < 10.:
          ndelta_vmic=n*delta_vmic
          vmica = max([0., vmic-ndelta_vmic])
          vmicb = max([3., vmic+ndelta_vmic])
          xa = [T, logg, vmica]
          xb = [T, logg, vmicb]
          evalmoda = interpol_mod(xa, met, fesun)
          evalmodb = interpol_mod(xb, met, fesun)
          if evalmoda and evalmodb:
              fa = error_v(vmica, T, logg, met, starname,  fesun, dslope2)
              fb = error_v(vmicb, T, logg, met, starname,  fesun, dslope2)
          n += 1
      if fa*fb <= 0.:
          return True,vmica,vmicb
      else:
          return False,0.,0.

def check_error_tm(x,met,starname,fesun,dslope1):
      T, logg, vmic =  x
      delta_T    =  100.
      fa, fb = [1.,1.]
      n = 1
      root = False
      while fa*fb > 0. and n < 10.:
          ndelta_T=n*delta_T
          Ta = T-ndelta_T
          Tb = T+ndelta_T
          xa = [Ta, logg, vmic]
          xb = [Tb, logg, vmic]
          evalmoda = interpol_mod(xa, met, fesun)
          evalmodb = interpol_mod(xb, met, fesun)
          if evalmoda and evalmodb:
              fa = error_tm(Ta, logg, vmic, met, starname, fesun, dslope1)
              fb = error_tm(Tb, logg, vmic, met, starname, fesun, dslope1)
          n += 1
      if fa*fb <= 0.:
          return True,Ta,Tb
      else:
          return False,0.,0.

def check_error_tv(x,met,starname,fesun):
      T, logg, vmic =  x
      delta_T    =  100.
      fa, fb = [1.,1.]
      n = 1
      root = False
      while fa*fb > 0. and n < 10.:
          ndelta_T=n*delta_T
          Ta = T-ndelta_T
          Tb = T+ndelta_T
          xa = [Ta, logg, vmic]
          xb = [Tb, logg, vmic]
          evalmoda = interpol_mod(xa, met, fesun)
          evalmodb = interpol_mod(xb, met, fesun)
          if evalmoda and evalmodb:
              fa = error_tv(Ta, logg, vmic, met, starname, fesun)
              fb = error_tv(Tb, logg, vmic, met, starname, fesun)
          n += 1
      if fa*fb <= 0.:
          return True,Ta,Tb
      else:
          return False,0.,0

def check_error_gd(x,met,starname,fesun,dAFe2):
      T, logg, vmic =  x
      delta_g    =  0.1
      fa, fb = [1.,1.]
      n = 1
      root = False
      while fa*fb > 0. and n < 10.:
          ndelta_g=n*delta_g
          logga = logg - ndelta_g
          loggb = logg + ndelta_g
          xa = [T, logga, vmic]
          xb = [T, loggb, vmic]
          evalmoda = interpol_mod(xa, met, fesun)
          evalmodb = interpol_mod(xb, met, fesun)
          if evalmoda and evalmodb:
              fa = error_gd(logga, T, vmic, met, starname, fesun, dAFe2)
              fb = error_gd(loggb, T, vmic, met, starname, fesun, dAFe2)
          n += 1
      if fa*fb <= 0.:
          return True,logga,loggb
      else:
          return False,0.,0.


def check_error_gp(x,met,starname,fesun):
      T, logg, vmic =  x
      delta_g    =  0.1
      fa, fb = [1.,1.]
      n = 1
      root = False
      while fa*fb > 0. and n < 10.:
          ndelta_g=n*delta_g
          logga = logg - ndelta_g
          loggb = logg + ndelta_g
          xa = [T, logga, vmic]
          xb = [T, loggb, vmic]
          evalmoda = interpol_mod(xa, met, fesun)
          evalmodb = interpol_mod(xb, met, fesun)
          if evalmoda and evalmodb:
              fa = error_gp(logga, T, vmic, met, starname, fesun)
              fb = error_gp(loggb, T, vmic, met, starname, fesun)
          n += 1
      if fa*fb <= 0.:
          return True,logga,loggb
      else:
          return False,0.,0.

def get_error_v(x, met, starname,fesun,slope2,sslope2):
    T, logg, vmic = x
    rootvp,vmica,vmicb = check_error_v(x,met,starname,fesun,slope2+sslope2)
    if rootvp:
        vmicn = optimize.brentq(error_v, vmica, vmicb, xtol=0.001, args=(T, logg, met, starname, fesun, slope2+sslope2))
        root  = True
        log_string="Root found +delta slope2\n"
    else:
        rootvn,vmica,vmicb = check_error_v(x,met,starname,fesun,slope2-sslope2)
        if rootvn:
            vmicn = optimize.brentq(error_v, vmica, vmicb, xtol=0.001, args=(T, logg, met, starname, fesun, slope2-sslope2))
            root = True
            log_string = "Root found -delta slope2\n"
        else:
            root = False
            log_string = "Could not find a root\n"

    if root: 
        evmic = np.abs(vmicn-vmic)
    else:
        evmic = -9.99

    return evmic,log_string


def get_error_tm(x, met, starname,fesun,slope1,sslope1):
    T, logg, vmic = x
    roottp,Ta,Tb = check_error_tm(x,met,starname,fesun,slope1+sslope1)
    if roottp:
        Tn = optimize.brentq(error_tm, Ta, Tb, xtol=0.001, args=(logg, vmic, met, starname, fesun, slope1+sslope1))
        root  = True
        log_string="Root found +delta slope1\n"
    else:
        roottn,Ta,Tb = check_error_tm(x,met,starname,fesun,slope1-sslope1)
        if roottn:
            vmicn = optimize.brentq(error_tm, Ta, Tb, xtol=0.001, args=(logg, vmic, met, starname, fesun, slope1-sslope1))
            root = True
            log_string = "Root found -delta slope1\n"
        else:
            root = False
            log_string = "Could not find a root\n"

    if root:
        eT = np.abs(Tn-T)
    else:
        eT = -9999
    return eT,log_string

def get_error_tv(x, met, starname,fesun,dvmic):
    T, logg, vmic = x
    xp=[T, logg, vmic+dvmic]
    roottp,Ta,Tb = check_error_tv(xp,met,starname,fesun)
    if roottp:
        Tn = optimize.brentq(error_tv, Ta, Tb, xtol=0.001, args=(logg, vmic+dvmic, met, starname, fesun))
        root  = True
        log_string="Root found +delta vmicro\n"
    else:
        xn=[T, logg, vmic-dvmic]
        roottn,Ta,Tb = check_error_tv(xn,met,starname,fesun)
        if roottn:
            Tn = optimize.brentq(error_tm, Ta, Tb, xtol=0.001, args=(logg, vmic-dvmic, met, starname, fesun))
            root = True
            log_string = "Root found -delta vmicro\n"
        else:
            root = False
            log_string = "Could not find a root\n"

    if root:
        eT = np.abs(Tn-T)
    else:
        eT = -9999
    return eT,log_string

def get_error_gd(x, met, starname,fesun, AFe2, sAFe2):
    T, logg, vmic = x
    rootgp,logga,loggb = check_error_gd(x,met,starname,fesun,AFe2+sAFe2)
    if rootgp:
        loggn = optimize.brentq(error_gd, logga, loggb, xtol=0.001, args=(T, vmic, met, starname, fesun,AFe2+sAFe2))
        root  = True
        log_string="Root found +delta sAFe2\n"
    else:
        rootgn,logga,loggb = check_error_gd(x,met,starname,fesun,AFe2-sAFe2)
        if rootgn:
            loggn = optimize.brentq(error_gd, logga, loggb, xtol=0.001, args=(T, vmic, met, starname, fesun, AFe-sAFe2))
            root = True
            log_string = "Root found -delta sAFe2\n"
        else:
            root = False
            log_string = "Could not find a root\n"

    if root:
        elg = np.abs(loggn-logg)
    else:
        elg = -9.99
    return elg,log_string

def get_error_gp(x, met, starname,fesun,spar):

    xp = np.array(x) + np.array(spar)
    xn = np.array(x) - np.array(spar)
    rootgp,logga,loggb = check_error_gp(xp,met,starname,fesun)
    if rootgp:
        T, logg, vmic = xp
        loggn = optimize.brentq(error_gp, logga, loggb, xtol=0.001, args=(T, vmic, met, starname, fesun))
        root  = True
        log_string="Root found +delta PAR\n"
    else:
        T, logg, vmic = xn
        rootgn,logga,loggb = check_error_gp(xn,met,starname,fesun)
        if rootgn:
            loggn = optimize.brentq(error_gp, logga, loggb, xtol=0.001, args=(T, vmic, met, starname, fesun))
            root = True
            log_string = "Root found -delta PAR\n"
        else:
            root = False
            log_string = "Could not find a root\n"

    if root:
        elg = np.abs(loggn-logg)
    else:
        elg = -9.99
    return elg,log_string

def get_error_metal(x, met, starname, fesun, sx):
    Teff, logg , vmic = x
    eTeff, elogg , evmic = sx
    xT = [Teff+eTeff, logg, vmic]
    xg = [Teff, logg+elogg, vmic]
    xv = [Teff, logg, vmic+evmic]
    evalmodel = interpol_mod(x,met,fesun)
    if evalmodel:
        slp1, slp2, af1, af2 = objective_function_vec(x, met, starname,fesun)
        evalmodelpT  = interpol_mod(xT,met,fesun)
        evalmodelpg = interpol_mod(xg,met,fesun)
        evalmodelpv  = interpol_mod(xv,met,fesun)
 
        if evalmodelpT:
            evalmodelT = True
        else:
            xT = [Teff-eTeff, logg, vmic]
            evalmodelT = interpol_mod(xT,met,fesun)

        if evalmodelpg:
            evalmodelg = True
        else:
            xg = [Teff, logg-elogg, vmic]
            evalmodelg = interpol_mod(xg,met,fesun)

        if evalmodelpv:
            evalmodelv = True
        else:
            xv = [Tef, logg, vmic-evmic]
            evalmodelv = interpol_mod(xv,met,fesun)

        if evalmodelT and evalmodelg and evalmodelv:
            slp1, slp2, af1t, af2t = objective_function_vec(xT, met, starname,fesun)
            slp1, slp2, af1g, af2g = objective_function_vec(xg, met, starname,fesun)
            slp1, slp2, af1v, af2v = objective_function_vec(xv, met, starname,fesun)
            eAFe1t = np.abs(af1t-af1)
            eAFe2t = np.abs(af2t-af2)
            eAFe1g = np.abs(af1g-af1)
            eAFe2g = np.abs(af2g-af2)
            eAFe1v = np.abs(af1v-af1)
            eAFe2v = np.abs(af2v-af2)
            eAFe1  = np.sqrt(eAFe1t**2. + eAFe1g**2. + eAFe1v**2.)
            eAFe2  = np.sqrt(eAFe2t**2. + eAFe2g**2. + eAFe2v**2.)
            #last evaluation is the final solution (you do not have to rerun MOOG/and recreate the atmospheric model)
            slp1, slp2, af1t, af2t = objective_function_vec(x, met, starname,fesun)
            return eAFe1,eAFe2
        else:
            return -9.99,-9.99
    else:
        return -9.99,-9.99
# Master error function:
def get_errors(starname, x, met, fesun):
    Teff, logg, vmic =x 
    slope1,slope2,AFe1,AFe2,sslope1,sslope2,sAFe1,sAFe2=get_err_slopes(x, met, starname, fesun)
    file_error_log=open('LOGS/'+starname+'_error.log','w')
    if slope1 > -9000.:
       evmic, log_string_vmic =get_error_v(x, met, starname, fesun, slope2, sslope2)
       eTm,log_string_Tm = get_error_tm(x, met, starname, fesun, slope1,sslope1)
       file_error_log.write(log_string_vmic)
       if evmic >= 0.:
           file_error_log.write("error on vmicro: {0:4.2f}\n".format(evmic))
           eTv,log_string_Tv = get_error_tv(x, met ,starname, fesun, evmic)
       else:
           eTv = -9999
       if eTv >= 0. and eTm >= 0.:
           eTeff = np.sqrt(eTv**2.+eTm**2.)
       else :
           eTeff = -9999

       file_error_log.write(log_string_Tm)
       file_error_log.write("error on Teff (slope1): {0:4.0f}\n".format(eTm))
       file_error_log.write(log_string_Tv)
       file_error_log.write("error on Teff (vmicro): {0:4.0f}\n".format(eTv))
       file_error_log.write("Total error on Teff: {0:4.0f}\n".format(eTeff))

       if evmic >= 0. and eTeff >= 0.:
           elgd,log_string_gd = get_error_gd(x, met, starname, fesun, AFe2, sAFe2)   
           dvecT=[eTeff,0.,0.]
           elgt,log_string_gt = get_error_gp(x, met, starname, fesun, dvecT)
           dvecV=[0.,0.,evmic]
           elgv,log_string_gv = get_error_gp(x, met, starname, fesun, dvecV)          
           if elgd >= 0. and  elgt >= 0. and elgv >= 0.:
               elogg = np.sqrt(elgd**2.+elgt**2.+elgv**2.)
           else:
               elogg = -9.99
       else:
           elogg = -9.99

       if elogg >= 0. and eTeff >= 0. and evmic >= 0.:
            sx=[eTeff,elogg,evmic]
            eAFe1,eAFe2 = get_error_metal(x, met, starname, fesun, sx)
            eAFe1 = np.sqrt(eAFe1**2.+sAFe1**2.)
            eAFe2 = np.sqrt(eAFe1**2.+sAFe2**2.)
       else:
            eAFe1,eAFe2 = [-9.99, -9.99]
       
       errs=[eTeff,elogg,evmic,eAFe1,eAFe2]
       string = '{0:5.0f} +- {1:5.0f} {2:5.2f} +- {3:5.2f} {4:5.2f} +- {5:5.2f} {6:5.2f} +- {7:5.2f} '.format(Teff,eTeff,logg,elogg,vmic,evmic,AFe1,eAFe1,AFe2,eAFe2)
       eAFe12=np.sqrt(eAFe1**2.+eAFe2**2.)
       string += ' {0:6.4f} +- {1:6.4f} {2:6.4f} +- {3:6.4f} {4:6.4f} +- {5:6.4f}'.format(slope1,sslope1,slope2,sslope2,AFe1-AFe2,eAFe12)
       file_error_log.close()
    else:
       errs=[-9999.9,-9.9,-9.9,-9.9,-9.9]
       string='NOPE'
       file_error_log.write("No errors could be found")
    file_error_log.close()
    return errs,string

#functions to write MOOG input

def write_par(starname):
    moog_par_file  = 'PAR/'  + starname + '_abfind.par'
    parfile=open(moog_par_file,'w')
    parfile.write('abfind\n')
    parfile.write('standard_out "STD_OUT/' + starname + '_Feout1"\n')
    parfile.write('summary_out "ABUN/' + starname + '_Feout2"\n')
    parfile.write('model_in    "out.mod"\n')
    parfile.write('lines_in    "EW/' + starname + 'Fe.l"\n')
    parfile.write('terminal     null\n')
    parfile.write('lines          1\n')
    parfile.write('freeform       0\n')
    parfile.write('flux/int       0\n')
    parfile.write('plot           0\n')
    parfile.write('damping        1\n')
    parfile.write('units          0\n')
    parfile.write('molecules      1\n')
    parfile.close()

def write_txt(starname):
    moog_txt_file  = 'TXT/'  + starname + '_abfind.rsp'
    moog_par_file  = 'PAR/'  + starname + '_abfind.par'
    txtfile=open(moog_txt_file,"w")
    txtfile.write(moog_par_file+'\n')
    txtfile.write('y\ny\ny\ny\n')
    txtfile.close()

#Open the model grid before anything else

def nelder_optimizer(starname,xin,metin,fesun,it_simp,it_res_simp):
    file_log=open('LOGS/'+starname+'.log','w')
    # We make sure to write down computation times.
    to = time()
    # We write the initial MOOG files to run StePar
    write_par(starname)
    write_txt(starname)
    # We print these to check everything's fine.
    print('Proceeding to analyse the following star: ' + str(starname))
    log_string = starname + '\n' + 'Initial parameters:\n' + \
        ' - Effective temperature (T eff): ' + str(xin[0]) + '\n' + \
        ' - Surface gravity (log g): ' + str(xin[1]) + '\n' + \
        ' - Micro-turbulence velocity (v micro): ' + str(xin[2]) + '\n'
    file_log.write(log_string)
    print(log_string)
    counter = 0
    T, logg, vmic = xin
    met  =  metin
    for i in range(0, it_res_simp):
        xin=[T,logg,vmic]
        metin = met
        slp1, slp2, af1, af2 = objective_function_vec(xin, metin, starname,fesun)
        # Here, the simplex will be oriented towards the mininum (acording to slope1, slope2, and Fe1-Fe2)
        # Steps on each parameter: Teff,logg, and vmic -> 200 K, 0.2 dex, and 0.2 km/s -> l1,l2,l3
        # 
        # Feel free to do as you please with the simplex initialization

        #Teff:
        
        if slp1 <= 0:
            l1 = -200
        else:
            l1 = +200
        
        #logg:
        
        if af1-af2 <= 0:
            l2 = -0.20
        else:
            l2 = +0.20

        # vmicro:
        
        if slp2 <= 0:
            l3 = -0.20
        else:
            l3 = +0.20

        # Generate the simplex. X0 to X3 are the points where we want to evaluate our objective function.
        
        X0 = [T, logg, vmic]
        X1 = [T + l1, logg, vmic]
        X2 = [T, logg + l2, vmic]
        X3 = [T, logg, vmic + l3]

        # We evaluate these simplex points.
        
        f0 = objective_function(X0, met, starname, fesun)
        f1 = objective_function(X1, met, starname, fesun)
        f2 = objective_function(X2, met, starname, fesun)
        f3 = objective_function(X3, met, starname, fesun)

        # 4-point matrix is then ordered.
        
        S = [[f0, X0], [f1, X1], [f2, X2], [f3, X3]]

        # New way of sorting S. S is sorted according to f0 to f4 values (minimum f 'number' value goes first).
       
        S.sort(key=lambda x: x[0])
        count_simp = 0
        [slp1, slp2, af1, af2] = objective_function_vec(S[0][1], met, starname,fesun)
        
        # While s1, s2, s31, s32 values are too big and number of iterations through this loop is less than 'it', then do this:
        while (np.abs(slp1) >= 0.001 or np.abs(slp1) >= 0.002 or np.abs(af1-af2) >= 0.005) and count_simp < it_simp:
            S = simplex(S, met, starname,fesun)
            [slp1, slp2, af1, af2] = objective_function_vec(S[0][1], met, starname,fesun)
            log_string = 'Iteration #' + str(count_simp + 1) + '\n'
            log_string += str(S[0][0]) + ' ' + format(S[0][1][0], '4.0f') + ' ' + \
                format(S[0][1][1], '4.2f') + ' ' + format(S[0][1][2], '5.3f') + ' ' \
                + str(slp1) + ' ' + str(slp2) + ' ' + str(af1-af2) + '\n'
            file_log.write(log_string)
            print(log_string)
            count_simp += 1
            # This is to know the values of the objective function on the other three evaluation points (to see where
            # the simplex is drifting to!).
            for j in range(1, 4):
                log_string = str(S[j][0]) + ' ' + format(S[j][1][0], '4.0f') + ' ' + \
                    format(S[j][1][1], '4.2f') + ' ' + format(S[j][1][2], '5.3f') + '\n'
                file_log.write(log_string)
                print(log_string)

        # Check if values fall into tolerance range already during this iteration.
        # Tolerance value check below:
        # - 0.001 for Fe 1 slope on the EP - ab plane for Fe 1 lines.
        # - 0.002 for Fe 2 slope on the EW - ab plane for Fe 1 lines.
        # - 0.005 for Fe 1 - Fe 2 abundance difference.
        if np.abs(slp1) <= 0.001 and np.abs(slp2) <= 0.002 and np.abs(af1-af2) <= 0.005 and np.abs(met-af1) < 0.01:
            counter += 1
        # If values have already fallen into the tolerance range twice, then the method is over
        # and we get out of the loop.
        if counter > 1:
            log_string = 'I\'m getting out of the loop at last!'
            log_string += 'The minimum of the objective function lies at: ' + \
                str(slp1) + ' ' + str(slp2) + ' ' + str(af1-af2) + ' ' + str(af1 - met) + '\n'
            log_string += 'Total computation time: ' + format((time() - to) / 60, '.1f') + ' minutes.'
            T=S[0][1][0]
            logg=S[0][1][1]
            vmic=S[0][1][2]
            met=af1
            log_string = '\n' + ' '.ljust(20) + '###############################################\n'
            log_string += ' '.ljust(20) + '#   Converged #' + str(i + 1) + ': '
            log_string += format(T, '4.0f') + ' ' + format(logg, '4.2f') + ' ' + format(vmic, '5.3f')+' '+ format(af1, '5.2f')
            log_string += '       #\n' + ' '.ljust(20) + '###############################################\n\n'
            file_log.write(log_string)
            print(log_string)
            return S[0][1],met
            break
        T=S[0][1][0]
        logg=S[0][1][1]
        vmic=S[0][1][2]
        met=af1
        log_string = '\n' + ' '.ljust(20) + '###############################################\n'
        log_string += ' '.ljust(20) + '#       Reset #' + str(i + 1) + ': '
        log_string += format(T, '4.0f') + ' ' + format(logg, '4.2f') + ' ' + format(vmic, '5.3f')+' '+ format(af1, '5.2f')
        log_string += '       #\n' + ' '.ljust(20) + '###############################################\n\n'
        file_log.write(log_string)
        print(log_string)
    file_log.close()
    return S[0][1],met

star = input()
gridMODS=open("MARCS1M.bin","rb")
tmod=pic.load(gridMODS)
gmod=pic.load(gridMODS)
mmod=pic.load(gridMODS)
ltaumod=pic.load(gridMODS)
Temod=pic.load(gridMODS)
lpgmod=pic.load(gridMODS)
lpemod=pic.load(gridMODS)
rhoxmod=pic.load(gridMODS)
kmod=pic.load(gridMODS)
gridMODS.close()
#test interpolator
xin=[5777.,4.438,1.00]
metin=0.0

#this is the solar abundance for Fe (Grevesse et al., 2007) compatible with MARCS. 
#You may want to alter this value if needed.
fesun = 7.45

xbad=[-9999.,-99.99,-99.99]
sxbad=[-9999.,-99.99,-99.99]
mbad=[-99.99]
xnew, metnew=nelder_optimizer(star, xin, metin, fesun, 60, 6)
tout  = outlier_razor(star, xnew, metnew, fesun)
sxnew, strnew = get_errors(star, xnew, metnew, fesun)
system('echo '+star+' '+strnew+' >> StePar_results' )
if tout:
    xres, metres = nelder_optimizer(star + 'new', xnew, metnew, fesun, 60, 6)
    sxres, strres = get_errors(star +'new',xres,metres, fesun)
    system('echo '+star+' '+strres+' >> StePar_results' )
else:
    system('echo '+star+'new  NOPE >> StePar_results')
