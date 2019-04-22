# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:21:12 2018

@author: Breno
"""
import glob
import sys
import numpy as np
import math
import os
import Options
import matplotlib.pyplot as plt
import time


def read_txt(Arquivo):
    f = open(Arquivo, 'r')
    
    txt = f.readlines()
    i = 0
    for line in txt:
        if line == '\n': ##### if testes do cobem ou anteriores: '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n'; Os novos usar penas '\n'
            begin_data = i + 1
        i+=1
    
    parametros = []
    data = {}
    
    lines = txt[begin_data].split('\t')
    
    for parametro in lines:
        if parametro == lines[-1]:
            parametro =  (parametro).replace('\n','')
            parametros.append(parametro)
        else:
            parametros.append(parametro)
        data[parametro] = []
    
    
    for i in range((begin_data+1),len(txt)):
        lines = txt[i].split('\t')
        j = 0 
        for line in lines:
            if j == (len(lines)-1):
                try:
                    data[parametros[j]].append(float(line.replace('\n','')))
                except ValueError:
                    data[parametros[j]].append(line.replace('\n',''))
            else:
                try:
                    data[parametros[j]].append(float(line))
                except ValueError:
                    data[parametros[j]].append(line)
            j+=1
    return data, parametros

def find_CB_HB(data,parametros,t_total):
    ###############################
    ### encontrar início do CB e HB
    ###############################
    
    # diferença entre PC in e PH out
    dp = data[' dP - CB (bar)']
    
    # São tomados n ciclos, entao utilizaremos somente o primeiro
    #t_total = len(data['P_Cin (bar)'])
    t_ciclo = int(t_total/Options.no_cycles)
    
    #determinando dp e tempo para t = t+ (delta t)
    dp2 = np.delete(dp,0) 
    time1 = data[' Run Time (s)']
    time2 = np.delete(data[' Run Time (s)'],0)
    
         
    
    # derivada da pressao com o tempo
    variacao_dp = (dp2[:t_ciclo]-dp[:t_ciclo])/(time2[:t_ciclo]-time1[:t_ciclo])
    
    sorted_dp = sorted(variacao_dp,reverse=True)
    max_dps = []
    min_dps = []
    indice_max = []
    indice_min = []
    i = 0
       
    while (len(min_dps) < 2):
        min_dp = sorted_dp[i*-1 -1]
        indice = np.where(variacao_dp == min_dp)[0][0]
    
        if ((variacao_dp[indice] - variacao_dp[indice -1]) < 0 and (variacao_dp[indice+1]- variacao_dp[indice]) > 0):
            min_dps.append(min_dp)
            indice_min.append(indice)
        i=i+1
    
    i=0
    while (len(max_dps) < 2):
        max_dp = sorted_dp[i]
        indice = np.where(variacao_dp == max_dp)[0][0]
    
        if ((variacao_dp[indice] - variacao_dp[indice -1])  > 0 and (variacao_dp[indice+1]- variacao_dp[indice]) < 0):
            max_dps.append(max_dp)
            indice_max.append(indice)
        i=i+1
    
    indice_min = np.hstack(indice_min)
    indice_max = np.hstack(indice_max)
    indices = np.sort(np.append(indice_max,indice_min))
    
    CB = np.zeros(2)
    HB = np.zeros(2)

    for i in range(len(indices)-1):
        result = variacao_dp[indices]
        result[i] = 0
        if (variacao_dp[indices[i]] > 0 and result[i+1] < 0):
            CB = np.array((indices[i],indices[i+1]))
        elif (variacao_dp[indices[i]] < 0 and result[i+1] > 0):
            HB = np.array((indices[i],indices[i+1]))
    
    if (CB[0] == 0 and CB[1] == 0):
        CB = np.array((indices[-1],indices[0]+(len(variacao_dp)-1)))
    elif (HB[0] == 0 and HB[1] == 0):
        HB = np.array((indices[-1],indices[0]+(len(variacao_dp)-1)))

    for i in range(len(parametros)):
        reorder = data[parametros[i]][0:(CB[0])]
        data[parametros[i]] = np.delete(data[parametros[i]],range((CB[0])))
        data[parametros[i]] = np.append(data[parametros[i]],reorder)

    HB = np.array(((HB[0]+Options.points_per_cycle-CB[0]),(HB[1]+Options.points_per_cycle-CB[0])))
    if (HB[0] >= 300):
        HB = HB - Options.points_per_cycle
    CB = np.array((0,(CB[1]-CB[0])))
    return data, CB, HB


def Simpson_integral(parameter,interval,cycles,samples,samp_freq):
    parameter_array = np.zeros(cycles)

    for cycle in range(cycles):
        SUM = 0
        
        boundary_1 = interval[0] + cycle*samples
        boundary_2 = interval[1] + cycle*samples
        
        SUM = np.sum(parameter[(boundary_1+3):(boundary_2-3)])
        
        for i in range(3):
            value1 = (parameter[(boundary_1+i)])
            value2 = (parameter[(boundary_2-i-1)])
            if i == 0:
                SUM += (3/8)*value1
                SUM += (3/8)*value2
            elif i == 1:
                SUM += (7/6)*value1
                SUM += (7/6)*value2     
            else:
                SUM += (23/24)*value1
                SUM += (23/24)*value2
                
        parameter_array[cycle] = (1/samp_freq)*SUM
        
    return parameter_array

###############################################################################################################################

######################
### Working Folder ###
######################

#path ="C:\\Users\\Alan\\Documents\\Alan Nakashima - Documentos\\AMR\\Eletrovalvula\\Análise_Teste_Cobem\\Dados_pydpp\\025"
#os.chdir(path)

####################
### Main Routine ###
####################


print('O programa está rodando.\n')

start_time = time.time()

Arquivos =  glob.glob("*.txt")
aviso = 'O(s) arquivo(s) :'


for Arquivo in Arquivos:
    if'_results' in Arquivo:
        Arquivos.remove(Arquivo)
        Arquivos.remove(Arquivo.replace("_results",""))
    else:
        aviso = aviso + ' ' + Arquivo + ';'
        
        
if Arquivos == []:
    print('Nenhum .txt foi encontrado ou todos os .txts já foram processados.')
    sys.exit()
aviso = aviso + ' foram encontrados. Iniciando pós-processamento! \n\n'

print(aviso)

for Arquivo in Arquivos:
    data, parametros = read_txt(Arquivo)

    data, CB, HB = find_CB_HB(data, parametros,Options.no_cycles*Options.points_per_cycle)

    no_samples = Options.no_cycles * Options.points_per_cycle
    if Options.sampling_freq == 0:
        Options.sampling_freq = (data[' Run Time (s)'][-1] - data[' Run Time (s)'][0])/no_samples
    
    
    
    Tamb_avg = (np.average(data['Tamb1 (K)']) + np.average(data['Tamb2 (K)']) + np.average(data['Tamb3 (K)']))/3.0
   
    
    #################################
    # Uncertainties
    #################################
    
    # Thermocouple
    
    In_Tpar = 0.15 #K
    In_Tpar_rect = 1
    if In_Tpar_rect == 1:
        In_Tpar = In_Tpar/math.sqrt(3)
    else:
        In_Tpar = In_Tpar/2
    
    u_dT = 2*In_Tpar
    
    # Pressure transducer
    
    In_PT = 0.5 #kPa
    In_PT_rect = 1
    if In_PT_rect == 1:
        In_PT = In_PT/math.sqrt(3)
    else:
        In_PT = In_PT/2
    
    u_dP = 2*In_PT
    
    # System Pressure
    
    In_SysP = 1 # %FS
    In_SysP_rect = 1
    
    if In_SysP_rect == 1:
        In_SysP = (In_SysP*20.68/100)/math.sqrt(3)
    else:
        In_SysP = (In_SysP*20.68/100)/2
        
    In_SysP = 2*In_SysP
    
    # Balança
    
    In_bal = 1 # % of the measure
    In_bal_rect = 1
    if In_bal_rect == 1:
        In_bal = (In_bal/1000)/math.sqrt(3)
    else:
        In_bal = (In_bal/1000)/2
    
    # Other uncertainties
        
    In_pq = 0.05/1000 #mm
    In_MoPa = 0.18 #deg
    In_Dpst = 0.025/1000 #mm
    In_Shft = 0.011/1000 #mm
    In_DReg = 0.05/1000  #mm
    In_LReg = 0.05/1000 #mm
    In_Stk = 0.05/1000 #mm
    In_Ang = 0.1 #deg
    In_EE = 0.005/1000 #mm
    
    u_pq = In_pq/math.sqrt(3)
    u_MoPa = (In_MoPa/360)/(math.sqrt(3)*Options.sampling_freq)
    u_Dpst = In_Dpst/math.sqrt(3)
    u_Shft = In_Shft/math.sqrt(3)
    u_DReg = In_DReg/math.sqrt(3)
    u_LReg = In_LReg/math.sqrt(3)
    u_Stk = In_Stk/math.sqrt(3)
    u_Ang = In_Ang/math.sqrt(3)
    u_geo = In_EE/math.sqrt(3)
    
    u_DReg = math.sqrt(u_DReg**2 + u_pq**2)
    u_LReg = math.sqrt(u_LReg**2 + u_pq**2)
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Period and Frequency Uncertainties
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    frequency = Options.frequency
    
    Period = 1/frequency
    
    TimeStep = Period/Options.points_per_cycle
    
    u_period = u_MoPa
    
    u_TimeStep = u_MoPa/Options.points_per_cycle
    
    u_freq = (u_period/(Period**2))
    
    #################################################
    # Define Hot, Cold and no blow periods
    #################################################
    
    No_flow_period_blow = (1 - Options.CB_frac - Options.HB_frac)/2
    
    begin_HB = HB[0]
    end_HB = HB[1]
    HB_period = end_HB - begin_HB
    
    begin_CB = CB[0]
    end_CB = CB[1]
    CB_period = end_CB - begin_CB
    
    CB_interval = CB_period*(1/Options.sampling_freq)
    HB_interval = HB_period*(1/Options.sampling_freq)
    
    ##################################################
    # Fluid's temperature leaving the HHEX (Hot Blow)
    ##################################################
    
    T_HHEX_out = Simpson_integral(data['T_Hhex-out  (K)'],HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_HHEX_out_avg = np.mean(T_HHEX_out/(HB_interval-1/Options.sampling_freq))
    
    SUM_DP = np.sum((T_HHEX_out - T_HHEX_out_avg)**2)
    
    DesvPad_T_HHEX = math.sqrt(SUM_DP/(Options.no_cycles-1))
    
    u_T_HHEX = math.sqrt(In_Tpar**2 + (DesvPad_T_HHEX**2)/(Options.no_cycles))
    
    
    ##################################################
    # Fluid's temperature leaving the CHEX (Cold Blow)
    ##################################################
    
    T_CHEX_out = Simpson_integral(data['T_Chex-out (K)'],CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_CHEX_out_avg = np.mean(T_CHEX_out/(CB_interval-1/Options.sampling_freq))
    
    SUM_DP = np.sum((T_CHEX_out - T_CHEX_out_avg)**2)
    
    DesvPad_T_CHEX = math.sqrt(SUM_DP/(Options.no_cycles-1))
    
    u_T_CHEX = math.sqrt(In_Tpar**2 + (DesvPad_T_CHEX**2)/(Options.no_cycles))
    
    #################################################
    # Fluid's temperature on the cold end (Hot Blow)
    #################################################
    
    parameter = data['Tint Cold End (K)']
    
    T_CE_HB = Simpson_integral(parameter,HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_CE_HB_avg = np.mean(T_CE_HB/(HB_interval-1/Options.sampling_freq))
    
    SUM_DP = np.sum((T_CE_HB - T_CE_HB_avg)**2)
    
    DesvPad_T_CE_HB = math.sqrt(SUM_DP/(Options.no_cycles-1))
    
    u_T_CE_HB = math.sqrt(In_Tpar**2 + (DesvPad_T_CE_HB**2)/(Options.no_cycles))
    
    #################################################
    # Fluid's temperature on the hot end (Cold Blow)
    #################################################
    
    parameter = data['Tint Hot End (K)']
    
    T_HE_CB = Simpson_integral(parameter,CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_HE_CB_avg = np.mean(T_HE_CB/(CB_interval-1/Options.sampling_freq))
    
    SUM_DP = np.sum((T_HE_CB - T_HE_CB_avg)**2)
    
    DesvPad_T_HE_CB = math.sqrt(SUM_DP/(Options.no_cycles-1))
    
    u_T_HE_CB = math.sqrt(In_Tpar**2 + (DesvPad_T_HE_CB**2)/(Options.no_cycles))

    #################################################
    # Fluid Properties Calculations (Etilenoglicol)
    #################################################
    
    glycol_wt = (0.0633642 + 1.10771*Options.EG - 0.00108914*(Options.EG**2))/100
    
    # Average temperature between hot and cold sides
    T_f = (T_HHEX_out_avg + T_CHEX_out_avg)/2
    
    # Density calculation
    C1 = (-517.378 + 14.3157*(T_f) - 0.0440498*(T_f**2) + 0.0000438395*(T_f**3))
    C2 = (1.089040385 - 0.01022673*(T_f) + 0.0000312579*(T_f**2) - 0.000000032614*(T_f**3) - 0.02518839*glycol_wt - 0.00060141*(glycol_wt**2))
    C3 = (1 - 0.009426*(T_f) + 0.0000289266*(T_f**2) - 0.000000030449*(T_f**3) - 0.01825238*glycol_wt)
    
    rho_f = C1*C2/C3
    
    # Viscosity calculation
    C1 = (0.79913 - 0.00976317*(T_f) + 0.0000449597*(T_f**2) - 0.000000092347*(T_f**3) + 0.0000000000713175*(T_f**4))
    C2 = (- 0.48083137 + 0.003165299*(T_f) - 0.0000055754*(T_f**2) - 0.04388738*glycol_wt - 0.0831036*(glycol_wt**2))
    C3 = (1 - 0.01052098*(T_f) + 0.0000362597*(T_f**2) - 0.000000042298*(T_f**3) + 0.026624841*glycol_wt)
    
    mi_f = C1*C2/C3
    
    # Conductivity calculation
    C1 = (- 0.680944 + 0.00687609*(T_f) - 0.00000871265*(T_f**2))
    C2 = (0.769354184 - 0.00648959*(T_f) - 0.0000040541*(T_f**2) - 2.86147298*glycol_wt + 0.060202185*(glycol_wt**2) + 0.010604554*(T_f)*glycol_wt)
    C3 = (1 - 0.00780393*(T_f) - 0.0000019933*(T_f**2) + 0.109891962*glycol_wt - 0.20530181*(glycol_wt**2) - 0.00296141*(T_f)*glycol_wt)
    
    k_f = C1*C2/C3
    
    # Specific heat calculation
    C1 = (15.60442824 - 13270.6349/(T_f) + 4017470/(T_f**2) - 405000000/(T_f**3) + 1.791510444*glycol_wt - 0.43299511*(glycol_wt**2) - 0.16418117*(glycol_wt**3) - 1061.16827*glycol_wt/(T_f) + 125.020179*(glycol_wt**2)/(T_f) + 130141.358*glycol_wt/(T_f**2))
    if T_f < 293:
        C2 = (245103 - 2503.25*(T_f) + 8.67151*(T_f**2) - 0.0100147*(T_f**3))
        cp_f = C2*C1
    else:
        cp_f = 4184*C1
    
    # kinematic density, thermal diffusivity and Prandtl number calculation
    ni_f = mi_f/rho_f
    dif_f = k_f/(rho_f*cp_f)
    Pr_f = ni_f/dif_f
    
    ##################
    # Solid Properties
    ##################
    # Gd
    rho_s = 7900
    k_s = 9
    c_s = 350
    
    ###################################
    # Incertezas das áreas e volumes
    ###################################
    
    A_reg = (math.pi/4)*((Options.D_casing/1000)**2)
    
    u_Areg = (math.pi/2)*math.sqrt((Options.D_casing*u_DReg)**2)
    
    Vol_reg = A_reg*Options.L_casing/1000
    
    u_VolReg = math.sqrt((Options.L_casing*u_Areg + A_reg*u_LReg)**2)
        
    ############
    # Porosidade
    ############

    mass = (Options.m_esf/1000)

    Porosity = 1 - mass/(rho_s*Vol_reg)
    
    u_por = math.sqrt((In_bal/(rho_s*Vol_reg))**2 + ((Options.m_esf*u_VolReg)/(rho_s*(Vol_reg)**2))**2)
    
    ###########################
    # Cold side flow (coriolis)
    ###########################
    
    MF_CS = Simpson_integral(data[' MFR - CB (kg/h)'],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    MF_CS_avg = np.mean(MF_CS/(CB_interval-1/Options.sampling_freq))
    
    SUM_DP = np.sum((MF_CS - MF_CS_avg)**2)
    
    DesvPad_MF_CS = math.sqrt(SUM_DP/(Options.no_cycles-1))
    
    In_Cor = 2 #%
    In_Cor_rect = 1
    if In_Cor_rect == 1:
        In_Cor = In_Cor/(100*math.sqrt(3))
    else:
        In_Cor = In_Cor/(100*2)
    
    u_MF_CS = math.sqrt((1/(Period**2))*((-MF_CS_avg*u_period)**2) + (1/(TimeStep**2))*((MF_CS_avg*u_TimeStep)**2) + (2/(Period*TimeStep))*(MF_CS_avg**2)*u_period*u_TimeStep + ((TimeStep/Period)*(Options.points_per_cycle-1)*In_Cor)**2 + (DesvPad_MF_CS**2)/(Options.no_cycles))
        
    ###########################
    # Hot side flow (coriolis)
    ###########################
    MF_HS  = Simpson_integral(data[' MFR - HB (kg/h)'],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    MF_HS_avg = np.mean(MF_HS/(HB_interval-1/Options.sampling_freq))
    
    SUM_DP = np.sum((MF_HS - MF_HS_avg)**2)
    
    DesvPad_MF_HS = math.sqrt(SUM_DP/(Options.no_cycles-1))
    
    u_MF_HS = math.sqrt((1/(Period**2))*((-MF_HS_avg*u_period)**2) + (1/(TimeStep**2))*((MF_HS_avg*u_TimeStep)**2) + (2/(Period*TimeStep))*(MF_HS_avg**2)*u_period*u_TimeStep + ((TimeStep/Period)*(Options.points_per_cycle-1)*In_Cor)**2 + (DesvPad_MF_HS**2)/(Options.no_cycles))

    ###################################
    # desbalanceamento da massa, utilizacao, vazao media, velocidade no regenerador...
    ###################################
    
    #MFR_avg = (MF_HS_avg + MF_CS_avg)/2
    
    #Mass_unbalance = 100*(MF_HS_avg - MF_CS_avg)/MFR_avg
    
    MFR_avg = (MF_HS + MF_CS)/2
    
    Mass_unbalance = 100*(MF_HS - MF_CS)/MFR_avg
    
    Mass_unbalance_avg = np.mean(Mass_unbalance)
    
    u_Mass_unbalance = Mass_unbalance_avg*math.sqrt(u_MF_HS**2 + u_MF_CS**2)*((MF_HS_avg - MF_CS_avg)**2+(1/(MF_HS_avg + MF_CS_avg)**2))
    
    #MFR_calc = np.mean(MFR_avg*(1/(3600)))
    
    #VFR_avg = MFR_avg*(1/(3600*1000))
    
    #u_MFRcalc = u_Ref_FR*rho_f*(1/(3600*1000))
    
    # fluid velocity inside the regenerator
    Vel_reg_HS = (MF_HS_avg/3600)/(A_reg*rho_f)
    
    u_Velreg_HS = math.sqrt(((u_MF_HS/3600)/(A_reg*rho_f))**2 + ((MF_HS_avg/3600*u_Areg)/((A_reg**2)*rho_f))**2)

    Vel_reg_CS = (MF_CS_avg/3600)/(A_reg*rho_f)
    
    u_Velreg_CS = math.sqrt(((u_MF_CS/3600)/(A_reg*rho_f))**2 + ((MF_CS_avg/3600*u_Areg)/((A_reg**2)*rho_f))**2)

    # Utilizacao ## REFAZER INCERTEZA DO PERIODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Utilization_HB = (HB_interval*(MF_HS_avg/3600)*cp_f)/(mass*c_s)

    u_Utlz_HB = 2*math.sqrt(((cp_f*(u_MF_HS/3600))/(Period*c_s*mass))**2 + ((cp_f*MF_HS_avg/3600*In_bal)/(Period*c_s*(mass**2)))**2 + ((cp_f*MF_HS_avg/3600*u_period)/((Period**2)*c_s*mass))**2)
    
    
    Utilization_CB = (CB_interval*(MF_CS_avg/3600)*cp_f)/(mass*c_s)

    u_Utlz_CB = 2*math.sqrt(((cp_f*(u_MF_CS/3600))/(Period*c_s*mass))**2 + ((cp_f*MF_CS_avg/3600*In_bal)/(Period*c_s*(mass**2)))**2 + ((cp_f*MF_CS_avg/3600*u_period)/((Period**2)*c_s*mass))**2)
   
        
    ###########################
    # Pressure loss (Cold Blow)
    ###########################
    
    Pres_Loss_CB = Simpson_integral(data[' dP - CB (bar)'],CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    Pres_Loss_CB_avg = np.mean(Pres_Loss_CB/(CB_interval-1/Options.sampling_freq))*100
    
    SUM_DP = np.sum((Pres_Loss_CB*100 - Pres_Loss_CB_avg)**2)
    
    DesvPad_Pres_Loss_CB = math.sqrt(SUM_DP/(Options.no_cycles-1))
    
    u_Pres_Loss_CB = math.sqrt((1/(Period**2))*((-Pres_Loss_CB_avg*u_period)**2) + (1/(TimeStep**2))*((Pres_Loss_CB_avg*u_TimeStep)**2) + (2/(Period*TimeStep))*(-DesvPad_Pres_Loss_CB**2)*u_period*u_TimeStep + (2*(TimeStep/Period)*(Options.points_per_cycle*Options.CB_frac-1)*u_dP)**2 + (DesvPad_Pres_Loss_CB**2)/(Options.no_cycles))
    
    ##########################
    # Pressure loss (Hot Blow)
    ##########################
    
    Pres_Loss_HB = Simpson_integral(data[' dP - HB (bar)'],HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    Pres_Loss_HB_avg = np.mean(Pres_Loss_HB/(HB_interval-1/Options.sampling_freq))*100
    
    SUM_DP = np.sum((Pres_Loss_HB*100 - Pres_Loss_HB_avg)**2)
    
    DesvPad_Pres_Loss_HB = math.sqrt(SUM_DP/(Options.no_cycles-1))
    
    u_Pres_Loss_HB = math.sqrt((1/(Period**2))*((-Pres_Loss_HB_avg*u_period)**2) + (1/(TimeStep**2))*((Pres_Loss_HB_avg*u_TimeStep)**2) + (2/(Period*TimeStep))*(-DesvPad_Pres_Loss_HB**2)*u_period*u_TimeStep + (2*(TimeStep/Period)*(Options.points_per_cycle*Options.HB_frac-1)*u_dP)**2 + (DesvPad_Pres_Loss_HB**2)/(Options.no_cycles))
    
    ######################
    # System Pressure drop
    ######################
    
    SysPres_Loss = Simpson_integral(data['DP_sys (bar)'],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    SysPres_Loss_avg = np.mean(SysPres_Loss/(Period-1/Options.sampling_freq))*100
    
    SUM_DP = np.sum((SysPres_Loss*100 - SysPres_Loss_avg)**2)
    
    DesvPad_SysPres_Loss = math.sqrt(SUM_DP/(Options.no_cycles-1))

    u_SysPres_Loss = math.sqrt((1/(Period**2))*((-SysPres_Loss_avg*u_period)**2) + (1/(TimeStep**2))*((SysPres_Loss_avg*u_TimeStep)**2) + (2/(Period*TimeStep))*(-DesvPad_SysPres_Loss**2)*u_period*u_TimeStep + (2*(TimeStep/Period)*(Options.points_per_cycle*(Options.HB_frac+Options.CB_frac)-1)*In_SysP)**2 + (DesvPad_SysPres_Loss**2)/(Options.no_cycles))
    
    
    Supply_Pressure = Simpson_integral(data['P_sysin (bar)'],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    Supply_Pressure_avg = np.mean(Supply_Pressure/(Period-1/Options.sampling_freq))*100
    
    
    #################################
    # Reference Flow Rate Calculation
    #################################
    
       
    Ref_flow = Simpson_integral(data['HighFlow (L/h)'],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    Ref_flow_avg = np.mean(Ref_flow/(Period-1/Options.sampling_freq))
    
    if Ref_flow_avg < 13:
        In_Ref_flow = 0.092/2
    else:
        In_Ref_flow = 0.463/2
    
    SUM_DP = np.sum((Ref_flow - Ref_flow_avg)**2)
    
    DesvPad_Ref_flow = math.sqrt(SUM_DP/(Options.no_cycles-1))

    u_Ref_flow = math.sqrt((1/(Period**2))*((-Ref_flow_avg*u_period)**2) + (1/(TimeStep**2))*((Ref_flow_avg*u_TimeStep)**2) + (2/(Period*TimeStep))*(-DesvPad_Ref_flow**2)*u_period*u_TimeStep + (2*(TimeStep/Period)*(Options.points_per_cycle*(Options.HB_frac+Options.CB_frac)-1)*In_Ref_flow)**2 + (Ref_flow_avg**2)/(Options.no_cycles))
   
    ########
    # Torque
    ########
    
    Torque = Simpson_integral(data['Torque (Nm)'],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    Torque_avg = np.mean(Torque/(Period-1/Options.sampling_freq))
    
    SUM_DP = np.sum((Torque - Torque_avg)**2)
    
    DesvPad_Torque = math.sqrt(SUM_DP/(Options.no_cycles-1))
    
    In_Torque_scale = 0.2 # xScale
    In_Torque_meas = 0.4 #%
    In_Torque_rect = 1
    if In_Torque_rect == 1:
        In_Torque = In_Torque_scale*Torque_avg/(5*math.sqrt(3)) + In_Torque_meas*Torque_avg/(100*math.sqrt(3))
    else:
        In_Torque = In_Torque_scale*Torque_avg/(5*2) + In_Torque_meas*Torque_avg/(100*2)
    
    u_Torque = math.sqrt((1/(Period**2))*((-Torque_avg*u_period)**2) + (1/(TimeStep**2))*((Torque_avg*u_TimeStep)**2) + (2/(Period*TimeStep))*((Torque_avg**2))*u_period*u_TimeStep + ((TimeStep/Period)*(Options.points_per_cycle-1)*In_Torque)**2 + (DesvPad_Torque**2)/(Options.no_cycles))
    
    Torque_total = Torque_avg
    
    ############################
    ### Valve Power and heat ###
    ############################
        
            
    Pot_Valve_CB_high = Simpson_integral(data["Power VCB High (W)"],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    Pot_Valve_CB_high_avg = np.mean(Pot_Valve_CB_high/(CB_interval-1/Options.sampling_freq))
    
    Pot_Valve_CB_low = Simpson_integral(data["Power VCB Low (W)"],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    Pot_Valve_CB_low_avg = np.mean(Pot_Valve_CB_low/(CB_interval-1/Options.sampling_freq))

    Pot_Valve_HB_high = Simpson_integral(data["Power VHB High (W)"],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    Pot_Valve_HB_high_avg = np.mean(Pot_Valve_HB_high/(HB_interval-1/Options.sampling_freq))    

    Pot_Valve_HB_low = Simpson_integral(data["Power VHB Low (W)"],[0,Options.points_per_cycle],Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    Pot_Valve_HB_low_avg = np.mean(Pot_Valve_HB_low/(HB_interval-1/Options.sampling_freq))
    
    
    ########################
    # Porous Media Variables
    ########################
    
    Dp = Options.D_esf/1000
    u_Dparticle = u_geo
    Sv_PM = 6/Dp
    u_Sv = 6*math.sqrt((u_Dparticle/(Dp**2))**2)
    Beta_PM = Sv_PM*(1-Porosity)
    u_beta = math.sqrt((Sv_PM*u_por)**2 + u_Sv**2 - 2*Porosity*(u_Sv**2) + (Porosity*u_Sv)**2)
    Dh_PM = 4*Porosity/Beta_PM
    u_Dh = math.sqrt((u_por/Beta_PM)**2 + (Porosity*u_beta/(Beta_PM**2))**2)
    A_int = Sv_PM*mass/rho_s
    u_Aint = (Sv_PM*In_bal)/rho_s
    
    # Dimensionless groups
    Re_w_PM = 2*math.pi*frequency*(Dh_PM**2)/ni_f
    
    Re_Dp_HS = Vel_reg_HS*Dp/ni_f
    Re_Dh_HS = Vel_reg_HS*Dh_PM/ni_f
    
    u_Rew = 2*math.pi*math.sqrt((((Dh_PM**2)*u_period)/(ni_f*(Period**2)))**2 + 4*((Dh_PM*u_Dh/(Period*ni_f))**2))
    u_ReDp = math.sqrt((Dp*u_Velreg_HS/ni_f)**2 + (Vel_reg_HS*u_Dparticle/ni_f)**2)
    u_ReDh = math.sqrt((Dh_PM*u_Velreg_HS/ni_f)**2 + (Vel_reg_HS*u_Dh/ni_f)**2)
    
    Re_Dp_CS = Vel_reg_CS*Dp/ni_f
    Re_Dh_CS = Vel_reg_CS*Dh_PM/ni_f
    
    u_Rew = 2*math.pi*math.sqrt((((Dh_PM**2)*u_period)/(ni_f*(Period**2)))**2 + 4*((Dh_PM*u_Dh/(Period*ni_f))**2))
    u_ReDp = math.sqrt((Dp*u_Velreg_CS/ni_f)**2 + (Vel_reg_CS*u_Dparticle/ni_f)**2)
    u_ReDh = math.sqrt((Dh_PM*u_Velreg_CS/ni_f)**2 + (Vel_reg_CS*u_Dh/ni_f)**2)
   
    # NTU (Correlação de Pallares)
    
    Nu_HS = 2*(1 + 4*(1-Porosity)/Porosity) + (1-Porosity)**0.5*Re_Dp_HS**0.6*Pr_f**(1/3)
    heat_coeff_HS = Nu_HS*k_f/Dp
    NTU_HS = heat_coeff_HS*A_int/(MF_CS_avg/3600*cp_f)
       
    Nu_CS = 2*(1 + 4*(1-Porosity)/Porosity) + (1-Porosity)**0.5*Re_Dp_CS**0.6*Pr_f**(1/3)
    heat_coeff_CS = Nu_CS*k_f/Dp
    NTU_CS = heat_coeff_CS*A_int/(MF_CS_avg/3600*cp_f)
    
    #############################
    ### HOT BLOW MEASUREMENTS ###
    #############################
    
    # Fluid's temperature entering the HHEX
    T_HHEX_in_HB = Simpson_integral(data['T_Hhex-in (K)'],HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_HHEX_in_HB_avg = np.mean(T_HHEX_in_HB/(HB_interval-1/Options.sampling_freq))
    
    # Fluid's temperature leaving the HHEX
    T_HHEX_out_HB = Simpson_integral(data['T_Hhex-out  (K)'],HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_HHEX_out_HB_avg = np.mean(T_HHEX_out_HB/(HB_interval-1/Options.sampling_freq))
    
    # Fluid's temperature entering hot end
    T_Hend_in_HB = Simpson_integral(data['T_Hot Side_in (K)'],HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_Hend_in_HB_avg = np.mean(T_Hend_in_HB/(HB_interval-1/Options.sampling_freq))
    
    # Fluid's temperature leaving cold end
    T_Cend_out_HB = Simpson_integral(data['T_Cold Side_out (K)'],HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_Cend_out_HB_avg = np.mean(T_Cend_out_HB/(HB_interval-1/Options.sampling_freq))
    
    # Manifold temperature
    
    T_Man_HB = Simpson_integral(data['T_MAN (K)'],HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_Man_HB_avg = np.mean(T_Man_HB/(HB_interval-1/Options.sampling_freq))
    
    # Valve exit temperature 
    
    T_Valve_HB = Simpson_integral(data['T_VHB (K)'],HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_Valve_HB_avg = np.mean(T_Valve_HB/(HB_interval-1/Options.sampling_freq))
    
    # Valve dT
    
    dT_valve_HB = Simpson_integral(data['T_VHB (K)']-data['T_MAN (K)'],HB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    dT_valve_HB_avg = np.mean(dT_valve_HB/(HB_interval-1/Options.sampling_freq))
    
    ##############################
    ### COLD BLOW MEASUREMENTS ###
    ##############################
    
    # Fluid's temperature entering the CHEX
    T_CHEX_in_CB = Simpson_integral(data['T_Chex-in (K)'],CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_CHEX_in_CB_avg = np.mean(T_CHEX_in_CB/(CB_interval-1/Options.sampling_freq))
    
    # Fluid's temperature leaving the CHEX
    T_CHEX_out_CB = Simpson_integral(data['T_Chex-out (K)'],CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_CHEX_out_CB_avg = np.mean(T_CHEX_out_CB/(CB_interval-1/Options.sampling_freq))
    
    # Fluid's temperature entering cold end
    T_Cend_in_CB = Simpson_integral(data['T_Cold Side_in (K)'],CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_Cend_in_CB_avg = np.mean(T_Cend_in_CB/(CB_interval-1/Options.sampling_freq))
    
    # Fluid's temperature leaving got end
    T_Hend_out_CB = Simpson_integral(data['T_Hot Side_out (K)'],CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_Hend_out_CB_avg = np.mean(T_Hend_out_CB/(CB_interval-1/Options.sampling_freq))
    
    # Manifold temperature
    
    T_Man_CB = Simpson_integral(data['T_MAN (K)'],CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_Man_CB_avg = np.mean(T_Man_CB/(CB_interval-1/Options.sampling_freq))
    
    # Valve exit temperature 
    
    T_Valve_CB = Simpson_integral(data['T_VCB (K)'],CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    T_Valve_CB_avg = np.mean(T_Valve_CB/(CB_interval-1/Options.sampling_freq))
    
    # Valve dT
    
    dT_valve_CB = Simpson_integral(data['T_VCB (K)']-data['T_MAN (K)'],CB,Options.no_cycles,Options.points_per_cycle,Options.sampling_freq)
    
    dT_valve_CB_avg = np.mean(dT_valve_CB/(CB_interval-1/Options.sampling_freq))
            
    ##############
    # Calculations
    ##############
    
    # Mechanical Power due to Torque
    Total_Pot_tq = Torque_total*2*math.pi*frequency/Options.no_polos
    u_Pot_tq = (2/Options.no_polos)*math.pi*math.sqrt((u_Torque/Period)**2 + (Torque_total*u_period/(Period**2))**2)
    
    # Pumping Power
    MP_Pot_HB = Pres_Loss_HB_avg*MF_HS_avg*(1000/3600)/rho_f
    MP_Pot_CB = Pres_Loss_CB_avg*MF_CS_avg*(1000/3600)/rho_f
    Total_Pump_Pot = MP_Pot_HB*HB_interval/Period + MP_Pot_CB*CB_interval/Period
    
    u_MP_Pot_HB = math.sqrt((Pres_Loss_HB_avg*(u_MF_HS/3600)/rho_f)**2+(MF_HS_avg/3600*u_Pres_Loss_HB/rho_f))
    u_MP_Pot_CB = math.sqrt((Pres_Loss_CB_avg*(u_MF_CS/3600)/rho_f)**2+(MF_CS_avg/3600*u_Pres_Loss_CB/rho_f))
    
    Sys_Pot = SysPres_Loss_avg*Ref_flow_avg*(1000/(3600*1000)) 
    u_Sys_Pot = math.sqrt((SysPres_Loss_avg*(u_Ref_flow/(3600*1000))/rho_f)**2 + (u_SysPres_Loss*(Ref_flow_avg/(3600*1000))/rho_f)**2) 
    
    # Valve Power
    
    Total_Valve_Pot = (Pot_Valve_CB_high_avg + Pot_Valve_CB_low_avg)*CB_interval/Period + (Pot_Valve_HB_high_avg + Pot_Valve_HB_low_avg)*HB_interval/Period
    
    ### COP and efficiencies
    # Temperature inside the regenerator
    T_cold_reg = T_CE_HB_avg
    T_hot_reg = T_HE_CB_avg
    
    # Temperature difference between sources
    dT_HEX = (T_HHEX_out_avg - T_CHEX_out_avg)
    u_dTHEX = math.sqrt(u_T_HHEX**2 + u_T_CHEX**2)
    
    # Temperature difference on the regenerator
    dT_Reg = (T_hot_reg - T_cold_reg)
    u_dTReg = math.sqrt(u_T_CE_HB**2 + u_T_HE_CB**2)
    
    # Temperature spans at the regenerator's entrances
    dT_Hend = T_hot_reg - T_HHEX_out_avg
    u_dTHend = math.sqrt(dT_Hend**2 + u_T_HHEX**2)
    
    dT_Cend = T_cold_reg - T_CHEX_out_avg
    u_dTCend = math.sqrt(dT_Cend**2 + u_T_CHEX**2)
    
    # Heat exchanged calculated at the regenerator's entrances
    Q_Hend = (MF_CS_avg/3600)*cp_f*(dT_Hend)*(CB_interval/Period)
    Q_Cend = -(MF_HS_avg/3600)*cp_f*(dT_Cend)*(HB_interval/Period)
    
    u_QHend = (1/2)*(1/3600)*math.sqrt((MF_HS_avg*cp_f*u_dTHend)**2 + (u_MF_HS*cp_f*dT_Hend)**2)
    u_QCend = (1/2)*(1/3600)*math.sqrt((MF_CS_avg*cp_f*u_dTCend)**2 + (u_MF_CS*cp_f*dT_Cend)**2)
    
    # Temperature spans at the regenerator's ports
    dT_Hport = T_Hend_out_CB_avg - T_Hend_in_HB_avg
    u_dTHend = math.sqrt(dT_Hend**2 + u_T_HHEX**2)
    
    dT_Cport = T_Cend_out_HB_avg - T_Cend_in_CB_avg
    u_dTCend = math.sqrt(dT_Cend**2 + u_T_CHEX**2)
    
    # Heat exchanged calculated at the regenerator's ports
    Q_Hport = (MF_CS_avg/3600)*cp_f*(dT_Hport)*(CB_interval/Period)
    Q_Cport = -(MF_HS_avg/3600)*cp_f*(dT_Cport)*(HB_interval/Period)
    
    u_QHend = (1/2)*(1/3600)*math.sqrt((MF_HS_avg*cp_f*u_dTHend)**2 + (u_MF_HS*cp_f*dT_Hend)**2)
    u_QCend = (1/2)*(1/3600)*math.sqrt((MF_CS_avg*cp_f*u_dTCend)**2 + (u_MF_CS*cp_f*dT_Cend)**2)
    
    
    
    # Tdm_work
    Tdm_work = Q_Hend-Q_Cend
    u_Tdmwork = math.sqrt(u_QHend**2 + u_QCend**2)
    
    # Total work performed
    W_total = Total_Pot_tq +Total_Pump_Pot + Total_Valve_Pot
    u_Wtotal = math.sqrt(u_Pot_tq**2 + u_Sys_Pot**2)
    
    # COP calculation
    COP_total = Q_Cend/W_total
    u_COPtotal = math.sqrt((Q_Cend*u_Wtotal/(W_total**2))**2 + (u_QCend)**2)
    
    # COP Carnot
    COP_Carnot = T_CHEX_out_avg/dT_HEX
    COP_Carnot_reg = T_cold_reg/dT_Reg
    
    u_COPCarnot = math.sqrt((T_CHEX_out_avg*u_dTHEX/(dT_HEX**2))**2 + (u_T_CHEX/dT_HEX)**2)
    u_COPCarnotreg = math.sqrt((T_cold_reg*u_dTReg/(dT_Reg**2))**2 + (u_T_CE_HB/dT_Reg)**2)
    
    rend_Tdm = COP_Carnot_reg/COP_Carnot
    rend_2nd = COP_total/COP_Carnot
    
    u_rendTdm = math.sqrt((COP_Carnot_reg*u_COPCarnot/(COP_Carnot**2))**2 + (u_COPCarnotreg/COP_Carnot)**2)
    u_rend2nd = math.sqrt((COP_total*u_COPCarnot/(COP_Carnot**2))**2 + (u_COPtotal/COP_Carnot)**2)
    
    # Valve Heat Gain
    
    Q_HB_valve = (MF_HS_avg/3600)*cp_f*(dT_valve_HB_avg)*(HB_interval/Period)
    Q_CB_valve = (MF_CS_avg/3600)*cp_f*(dT_valve_CB_avg)*(CB_interval/Period)
    
    Q_valve_total = Q_HB_valve*HB_interval/Period + Q_CB_valve*CB_interval/Period
    
    ###########################
    ### Main Results Writing###
    ###########################
    
    # Results Selection
    
    Parameters_names = ["Hex Span [K] =","Frequency [Hz] =","CB Blow Fraction [%] =","CB Utilization [-] =",
                        "HB Blow Fraction [%] =","HB Utilization [-] =", "T Amb [K] ="] 
    
    Parameters_results = [dT_HEX, frequency, CB_interval/Period, Utilization_CB,
                          HB_interval/Period, Utilization_HB,Tamb_avg]
    
    Results_names = ["Cooling Capacity [W] =","Rejected Heat [W] =","Blows Pumping Power [W] =","Torque Power [W] =","Valve Power [W] =","Total Power [W] =",
                     "COP [-] =","Carnot COP [-] =", "eta_2nd [%] ="]
    
    Results_results = [Q_Cend, Q_Hend, Total_Pump_Pot, Total_Pot_tq, Total_Valve_Pot, W_total, COP_total, COP_Carnot, rend_2nd]
    
    Flow_names = ["CB Pressure Drop [kPa] =","CB MFR [kg/h] =","CB Pumping Power [W] =","HB Pressure Drop [kPa] =","HB MFR [kg/h] =","HB Pumping power [W] =",
                  "System Pressure Drop [kPa] =","System Flow Rate [L/h] =","System Pumping Power [W] =","Supply Pressure [kPa] =","Mass Imbalance [%] ="]
    
    Flow_results = [Pres_Loss_CB_avg, MF_CS_avg, MP_Pot_CB, Pres_Loss_HB_avg, MF_HS_avg, MP_Pot_HB,
                    SysPres_Loss_avg, Ref_flow_avg, Sys_Pot, Supply_Pressure_avg, Mass_unbalance_avg]
    
    Actuation_names = ["Torque [Nm] =","CB High Valve Power [W] =","CB Low Valve Power [W] =","HB High Valve Power [W] =","HB Low Valve Power [W] =",
                       "Valve Heat gain [W] ="]
    
    Actuation_results = [Torque_total, Pot_Valve_CB_high_avg , Pot_Valve_CB_low_avg , Pot_Valve_HB_high_avg , Pot_Valve_HB_low_avg , Q_valve_total]
    
    Temperature_names = ["DT Cold end [K] =","T Cold Hex [K] =","T Cold end [K] =","T Cold in [K] =","T Cold out [K] =",
                         "DT Hot end [K] =","T Hot Hex [K] =","T Hot end [K] =","T Hot in [K] =","T Hot out [K] =",
                         "T_manifold_CB [K] =","T_CB_high_valve [K] =","T_manifold_CB [K] =","T_HB_high_valve [K] ="] 
    
    Temperature_results = [dT_Cend, T_CHEX_out_avg, T_CE_HB_avg, T_Cend_in_CB_avg, T_Cend_out_HB_avg,
                         dT_Hend, T_HHEX_out_avg, T_HE_CB_avg, T_Hend_in_HB_avg, T_Hend_out_CB_avg,
                         T_Man_CB_avg,T_Valve_CB_avg,+T_Man_HB_avg,T_Valve_HB_avg]\
    
    Regenerator_names = ["Casing Diameter [mm] =","Casing Length [mm] =","Particle Diameter [mm] =","Porosity [-] =","Compactness Factor [m^2/m^3] =", 
                         "CB Re_dp [-] =", "CB NTU [-] =", "HB Re_dp [-] =", "HB NTU [-] ="]    
    
    Regenerator_results = [Options.D_casing,Options.L_casing,Options.D_esf,Porosity,Beta_PM,
                           Re_Dp_CS,NTU_CS,Re_Dp_HS,NTU_HS]
        
    filename = Arquivo.replace(".txt","") + "_results.txt"
    file_object  = open(filename, "w")
    
    file_object.write("Test Results\n")
    
    file_object.write("\nParameters:\n\n")
    for i in range(len(Parameters_names)):
        file_object.write(Parameters_names[i]+" "+str(round(Parameters_results[i],2))+"\n")
    
    file_object.write("\nMain Results:\n\n")
    for i in range(len(Results_names)):
        file_object.write(Results_names[i]+" "+str(round(Results_results[i],2))+"\n")
           
    file_object.write("\nFlow results:\n\n")
    for i in range(len(Flow_names)):
        file_object.write(Flow_names[i]+" "+str(round(Flow_results[i],2))+"\n")
    
    file_object.write("\nTransmission and valves results:\n\n")
    for i in range(len(Actuation_names)):
        file_object.write(Actuation_names[i]+" "+str(round(Actuation_results[i],2))+"\n")
    
    file_object.write("\nSystem temperatures:\n\n")
    for i in range(len(Temperature_names)):
        file_object.write(Temperature_names[i]+" "+str(round(Temperature_results[i],2))+"\n")
    
    file_object.write("\nRegenerator Parameters:\n\n")
    for i in range(len(Regenerator_names)):
        file_object.write(Regenerator_names[i]+" "+str(round(Regenerator_results[i],2))+"\n")
    
    
    file_object.close()
    
  
    
print("--- %s seconds ---" % (time.time() - start_time))
#import pdb; pdb.set_trace()
