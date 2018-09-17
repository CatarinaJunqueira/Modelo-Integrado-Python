# -*- coding: utf-8 -*-

from __future__ import print_function
from scipy.io import loadmat
from connect_ifttt import email_alert
import numpy as np
import cplex
#from compare_models import compare_models
#def b[i][j][n][t]

def ModeloIntegrado(casename):
    
    data = loadmat(casename)
    C = data['C'][0][0]
    R = data['R'][0][0]
    Patios = data['Patios'].tolist()
    TT =data['TT']
    phi = data['phi'].tolist()
    Npatios = len(Patios)
    P=Npatios+1 # numero de portos
    
    omega=[ [] for i in range(Npatios) ] # omega = conjunto dos indices dos conteineres em cada patio
    for i in range(Npatios):
        Patios[i]=Patios[i][0]
        omega[i]=np.extract(Patios[i]!= 0 , Patios[i]).tolist()
        
    for o in range(Npatios):
        for d in range(P):
            if phi[o][d].shape[1] != 0 :
                phi[o][d] = phi[o][d].tolist()[0]         
            else:
                phi[o][d] = []
    
    
    N=[ 0 for i in range(Npatios) ] # N = quantidade de conteineres em cada patio
    for i in range(Npatios):
        N[i]=np.count_nonzero(Patios[i])
    
    T=N
    
    H=[] # H = numero de linhas de cada patio
    for i in range(Npatios):
        H.append(Patios[i].shape[0])
    
    W= [] # W = numero de colunas de cada patio
    for i in range(Npatios):
        W.append(Patios[i].shape[1])
    
    
    model = cplex.Cplex()
    start_time = model.get_time()
    model.objective.set_sense(model.objective.sense.minimize)
    
    #------------------------------------------------------------#
    #--------------------  Variaveis  ---------------------------#
    nvar = 0 
    
    model,xnames,nvar = variavel_x(model,omega,N,H,W,T,nvar)
    model,bnames,nvar = variavel_b(model,omega,N,H,W,T,nvar)
    model,vnames,nvar = variavel_v(model,omega,N,T,nvar)
    model,ynames,nvar = variavel_y(model,omega,N,H,W,T,nvar)
    model,znames,nvar = variavel_z(model,omega,N,T,R,C,nvar)
    model,qnames,nvar = variavel_q(model,N,R,C,nvar)
    model,wnames,nvar = variavel_w(model,N,R,C,nvar)
    model,unames,nvar = variavel_u(model,N,R,C,nvar)
    
    
    #------------------------------------------------------------#
    #-------------------  Restricoes  ---------------------------#
    #Patio:
    model = restricao_P0(model,omega,N,H,W,Patios)
    model = restricao_P1(model,omega,N,H,W,T)
    model = restricao_P2(model,omega,N,H,W,T)
    model = restricao_P3(model,omega,N,H,W,T)
    model = restricao_P6(model,omega,N,H,W,T)
    model = restricao_P7(model,omega,N,H,W,T)
    model = restricao_P8(model,omega,N,H,W,T)
    model = restricao_P9(model,omega,N,H,W,T)
    model = restricao_P10(model,omega,N,H,W,T)
    model = restricao_PA(model,omega,N,H,W,T)
    model = restricao_I1(model,omega,N,T)
    model = restricao_I2(model,omega,N,T)
    model = restricao_I3(model,omega,N,R,C,T)
    model = restricao_I4(model,omega,N,R,C,T,P)
    model = restricao_I5(model,omega,N,R,C,T)
    model = restricao_I6(model,phi,R,C,T,P,N)
    model = restricao_I7(model,omega,N,R,C,T)
    model = restricao_I8(model,omega,N,R,C,T,P)
    model = restricao_I9(model,omega,N,R,C,P)
    model = restricao_I10(model,omega,N,R,C,TT)
    model = restricao_N1(model,P,R,C,TT)
    model = restricao_N2(model,R,C,P)
    model = restricao_N3(model,N,R,C)    
    model = restricao_N4(model,P,R,C)
    
    #model.write("model_python.lp")
    
   # histogram = model.variables.get_histogram()
   # histogram = np.sum(CHist,axis=1)
  #  print(histogram)
    
    variaveis = model.variables.get_num()
    print('Numero de Variaveis = ',variaveis)
    
    restricoes = model.linear_constraints.get_num()
    print('Numero de Restricoes = ',restricoes)
    
    z = 'No modelo ha %s variaveis e %s restricoes' %(variaveis,restricoes)
        
    out_file = open("Resultados.txt",'w+')
   
    model.set_results_stream(out_file)
   
    #model.set_results_stream("Resultados_InstanciaModeloIntegrado_II.txt")
    
    model.solve()
    
    out_file.seek(0)
    out_string =out_file.read()
    out_file.close()
    
    print(out_string)
   
    end_time = model.get_time()
    solvetime = end_time-start_time
   # print('Duracao = ',solvetime)
    
    status = model.solution.get_status_string()
    fobj = model.solution.get_objective_value()
    
    print('\nSolution status = ',status)
    print('Valor da Funcao Objetivo: ',fobj )
    
    Y = ' Solution status: %s  <br> Valor da Funcao Objetivo: %s  <br> Duracao: %s <br>' %(status,fobj,solvetime)
   
    email_alert(Y, z, out_string)
    return model

#------------------------------------------------------------------------------------------------------------------#
#--------------------  Variaveis  ---------------------------#
#------------------------------------------------------------------------------------------------------------------#  
def variavel_x(model,omega,N,H,W,T,nvar):
    xnames = []
    
    for o in range(len(N)):  
        for i in range(1,W[o]+1): 
            for j in range(1,H[o]+1):
                for k in range(1,W[o]+1):
                    for l in range(1,H[o]+1):
                        for n in omega[o]:
                            for t in range(1,T[o]+1):
                                xnames.append('x_'+str(i)+'_'+str(j)+'_'+str(k)+'_'+str(l)+'_'+str(n)+'_'+str(t))
    
    nx = len(xnames)
    
    nvar += nx
    
    lb = [0.0]*nx
    ub = [1.0]*nx
    ctypes =[model.variables.type.binary]*nx
    obj = [1.0]*nx
    model.variables.add(obj=obj, lb=lb, ub=ub, types=ctypes,names=xnames)
    
    return model,xnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_b(model,omega,N,H,W,T,nvar):
    bnames = []
    
    for o in range(len(N)):  
        for i in range(1,W[o]+1): 
            for j in range(1,H[o]+1):
                for n in omega[o]:
                    for t in range(1,T[o]+1):
                        bnames.append('b_'+str(i)+'_'+str(j)+'_'+str(n)+'_'+str(t))
    
    nb = len(bnames)
    
    nvar += nb
    
    lb = [0.0]*nb
    ub = [1.0]*nb
    ctypes =[model.variables.type.binary]*nb
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=bnames)
    
    return model,bnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#
def variavel_v(model,omega,N,T,nvar):
    vnames = []
    
    for o in range(len(N)):  
        for n in omega[o]:
            for t in range(1,T[o]+1):
                vnames.append('v_'+str(n)+'_'+str(t))
    
    nv = len(vnames)
    
    nvar += nv
    
    lb = [0.0]*nv
    ub = [1.0]*nv
    ctypes =[model.variables.type.binary]*nv
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=vnames)
    
    return model,vnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_y(model,omega,N,H,W,T,nvar):
    ynames = []
    
    for o in range(len(N)):  
        for i in range(1,W[o]+1): 
            for j in range(1,H[o]+1):
                for n in omega[o]:
                    for t in range(1,T[o]+1):
                        ynames.append('y_'+str(i)+'_'+str(j)+'_'+str(n)+'_'+str(t))
    
    ny = len(ynames)
    
    nvar += ny
    
    lb = [0.0]*ny
    ub = [1.0]*ny
    ctypes =[model.variables.type.binary]*ny
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=ynames)
    
    return model,ynames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_z(model,omega,N,T,R,C,nvar):
    znames = []
    
    for o in range(len(N)):  
        for n in omega[o]:
            for t in range(1,T[o]+1):
                for r in range(1,R+1):
                    for c in range(1,C+1):
                        znames.append('z_'+str(n)+'_'+str(t)+'_'+str(r)+'_'+str(c))
    
    nz = len(znames)
    
    nvar += nz
    
    lb = [0.0]*nz
    ub = [1.0]*nz
    ctypes =[model.variables.type.binary]*nz
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=znames)
    
    return model,znames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_q(model,N,R,C,nvar):
    qnames = []
    
    for o in range(1,len(N)+1):  
        for d in range(o+1,len(N)+2):
              for r in range(1,R+1):
                  for c in range(1,C+1):
                      qnames.append('q_'+str(o)+'_'+str(d)+'_'+str(r)+'_'+str(c))
    
    nq = len(qnames)
    
    nvar += nq
    
    lb = [0.0]*nq
    ub = [1.0]*nq
    ctypes =[model.variables.type.binary]*nq
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=qnames)
    
    return model,qnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_w(model,N,R,C,nvar):
    wnames = []                        
    obj = []

    for o in range(1,len(N)+1):  
        for d in range(o+1,len(N)+2):
            for a in range(o+1,d+1):
                for r in range(1,R+1):
                    for c in range(1,C+1):
                        
                        wnames.append('w_'+str(o)+'_'+str(d)+'_'+str(a)+'_'+str(r)+'_'+str(c))
                        if a == d:
                            obj.append(0.0)
                        else:
                            obj.append(1.0)
    
    nw = len(wnames)
    
    nvar += nw
    
    lb = [0.0]*nw
    ub = [1.0]*nw
    ctypes =[model.variables.type.binary]*nw
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=wnames)
    
    return model,wnames,nvar

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_u(model,N,R,C,nvar):
    unames = []
    
    for o in range(1,len(N)+1):  
        for r in range(1,R+1):
            for c in range(1,C+1):
                unames.append('u_'+str(o)+'_'+str(r)+'_'+str(c))
    
    nu = len(unames)
    
    nvar += nu
    
    lb = [0.0]*nu
    ub = [1.0]*nu
    ctypes =[model.variables.type.binary]*nu
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=unames)
    
    return model,unames,nvar


#------------------------------------------------------------------------------------------------------------------#
#-------------------  Restricoes  ---------------------------#
#------------------------------------------------------------------------------------------------------------------#
def restricao_P0(model,omega,N,H,W,Patios):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for j in range(H[o]):
                for n in omega[o]:
                    rest_names.append('restP0_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(o+1))

                    if Patios[o][H[o]-j-1,i] == n :
                        rhs.append(1.0)
                    else :
                        rhs.append(0.0)

                    #add b coeficientes
                    cols = ['b_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(1)]
                    coefs = [1.0]
                    expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 1: In each time period, each block must either be within the stack or in the outside region:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P1(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for n in omega[o]:
            for t in range(T[o]):
                rest_names.append('restP1_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(1.0)
                cols = ['v_'+str(n)+'_'+str(t+1)]
                coefs = [1.0]
                
                for i in range(W[o]):  
                    for j in range(H[o]):
                        #add b coeficientes
                        cols.append('b_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t+1))
                        coefs.append(1.0)
                expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 2: In each time period, each slot (i,j) must be occupied by at most one block:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P2(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]):  
            for j in range(H[o]):
                for t in range(T[o]):
                    rest_names.append('restP2_'+str(i+1)+'_'+str(j+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(1.0)
                    cols = []
                    coefs = []
                    for n in omega[o]:          
                        #add b coeficientes
                        cols.append('b_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t+1))
                        coefs.append(1.0)
                        
                    expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model


#------------------------------------------------------------------------------------------------------------------#
#Constraint 3: garante que não hajam ‘buracos’ no pátio ao restringir que se há um contêiner posição $(i,j+1)$, 
#então a posição $(i,j)$ abaixo também deve estar ocupada:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P3(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]):  
            for j in range(H[o]-1):
                for t in range(T[o]):
                    rest_names.append('restP3_'+str(i+1)+'_'+str(j+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(0.0)
                    cols = []
                    coefs = []
                    for n in omega[o]:          
                        #add b coeficientes
                        cols.append('b_'+str(i+1)+'_'+str(j+2)+'_'+str(n)+'_'+str(t+1))
                        coefs.append(1.0)
                        cols.append('b_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t+1))
                        coefs.append(-1.0)
                        
                    expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 4: restrição de equilíbrio de fluxo entre as variáveis de configuração e de movimento no pátio.  
#Vincula o layout no período t com o layout no período t + 1 através das retiradas e realocações executadas:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P6(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]):  
            for j in range(H[o]):
                for n in omega[o]:
                    for t in range(1,T[o]):
                        rest_names.append('restP6_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                        rhs.append(0.0)
                        cols = ['b_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t+1)] #b_{i,j,n,t}
                        coefs = [1.0]
                        cols.append('b_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t)) #-b_{i,j,n,t-1}
                        coefs.append(-1.0)
                        cols.append('y_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t)) #y_{i,j,n,t-1}
                        coefs.append(1.0)                        
                        for k in range(W[o]):          
                            for l in range(H[o]):
                                if k != i or l !=j :
                                    cols.append('x_'+str(k+1)+'_'+str(l+1)+'_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t)) #-\sum_{k=1}^{W}\sum_{l=1}^{H}x_{k,l,i,j,n,t-1}
                                    coefs.append(-1.0)
                                    cols.append('x_'+str(i+1)+'_'+str(j+1)+'_'+str(k+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t)) #\sum_{k=1}^{W}\sum_{l=1}^{H}x_{i,j,k,l,n,t-1}
                                    coefs.append(1.0)
                        
                        expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 5:  define a variável $v_{nt}$ e assegura que todos os contêineres sejam retirados do pátio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P7(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for n in omega[o]:
            for t in range(1,T[o]):
                rest_names.append('restP7_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(0.0)
                cols = ['v_'+str(n)+'_'+str(t+1)] #v_{nt}
                coefs = [1.0]                       
                for i in range(W[o]):          
                    for j in range(H[o]):
                        for tt in range(t): #\sum_{k=1}^{W}\sum_{l=1}^{H}x_{k,l,i,j,n,t-1}
                            cols.append('y_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(tt+1)) #-\sum_{i=1}^{W}\sum_{j=1}^{H}\sum_{t'=1}^{t-1}y_{ijnt'} 
                            coefs.append(-1.0)
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 6:  garante a política LIFO, ou seja, se no período t, o contêiner $n$ está abaixo do contêiner $q$  
# e o contêiner $n$ é remanejado, então no período $t + 1$ o contêiner $n$ não pode estar alocado em uma posição 
# abaixo do contêiner $q$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P8(model,omega,N,H,W,T):
    
    M=0
    for o in range(len(N)):
        M=M+(N[o]*((H[o]-1)**2))   
        
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for k in range(W[o]):
                for j in range(H[o]-1):
                    for l in range(H[o]-1):
                        for t in range(T[o]):
                            rest_names.append('restP8_'+str(i+1)+'_'+str(j+1)+'_'+str(k+1)+'_'+str(l+1)+'_'+str(t+1)+'_'+str(o+1))
                            rhs.append(M)
                            cols = [] 
                            coefs = []
                            
                            for n in omega[o]: 
                                cols.append('x_'+str(i+1)+'_'+str(j+1)+'_'+str(k+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t+1)) #\sum_{n=1}^{N}x_{i,j,k,l,n,t}
                                coefs.append(M)
                            
                                for jj in range(j+1,H[o]):
                                    for ll in range(l+1,H[o]):
                                        cols.append('x_'+str(i+1)+'_'+str(jj+1)+'_'+str(k+1)+'_'+str(ll+1)+'_'+str(n)+'_'+str(t+1))
                                        coefs.append(1.0)               
                        
                            expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 7:  garante que sejam remanejados apenas os contêineres que estão acima, ou seja, na mesma coluna, 
# de um contêiner a ser retirado:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P9(model,omega,N,H,W,T):
    
    M=0
    for o in range(len(N)):
        M=M+((H[o]**2)*(W[o]**2)*N[o]*N[o]+1)   
        
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]):
            for t in range(T[o]):
                for n in omega[o]: 
                    rest_names.append('restP9_'+str(n)+'_'+str(i+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(M)
                    cols = [] 
                    coefs = []
                    
                    for j in range(H[o]):
                        cols.append('b_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t+1)) 
                        coefs.append(M)
                        for k in range(W[o]):
                            for l in range(H[o]):
                                for ii in range(i):
                                    cols.append('x_'+str(ii+1)+'_'+str(j+1)+'_'+str(k+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t+1)) 
                                    coefs.append(1.0)
                                for iii in range(i+1,W[o]):
                                    cols.append('x_'+str(iii+1)+'_'+str(j+1)+'_'+str(k+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t+1)) 
                                    coefs.append(1.0)          
                        
                    expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 8:  garante que nenhum contêiner pode ser remanejado para outra posição que esteja da mesma coluna na
# qual ele se encontra:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P10(model,omega,N,H,W,T):
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for j in range(H[o]):
                for l in range(H[o]):
                    for n in omega[o]: 
                        for t in range(T[o]):
                            rest_names.append('restP10_'+str(i+1)+'_'+str(j+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                            rhs.append(0.0)
                            cols= ['x_'+str(i+1)+'_'+str(j+1)+'_'+str(i+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t+1)]
                            coefs = [1.0]        
                        
                            expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 9:  garante que um contêiner na posição $(i,j)$ só pode ser movido depois que o contêiner na posição 
# $(i,j+1)$ é movido. Se o contêiner na posição $(i,j+1)$ não é movido então temos que $b_{i(j+1)nt} = 1$ e 
# $x_{i(j+1)klnt} = 0$, e o lado esquerdo da equação se torna 0. Consequentemente o lado direito da equação também 
# deve ser 0. Dessa forma, nenhuma realocação ou remanejamento é permitido para o contêiner na posição $(i,j)$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_PA(model,omega,N,H,W,T):
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for j in range(H[o]-1):
                for t in range(T[o]):
                    rest_names.append('restA_'+str(i+1)+'_'+str(j+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(1.0)
                    cols = [] 
                    coefs = []
                    for n in omega[o]:
                        cols.append('b_'+str(i+1)+'_'+str(j+2)+'_'+str(n)+'_'+str(t+1)) #\sum_{n=1}^{N}b_{i,j+1,n,t}
                        coefs.append(1.0)
                        cols.append('y_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t+1)) #\sum_{n=1}^{N}y_{i,j,n,t}
                        coefs.append(1.0)
                        for k in range(W[o]):
                            for l in range(H[o]):
                                cols.append('x_'+str(i+1)+'_'+str(j+2)+'_'+str(k+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t+1)) #\sum_{k=1}^{W}\sum_{l=1}^{H}\sum_{n=1}^{N}x_{i,j+1,k,l,n,t}
                                coefs.append(-1.0)  
                                cols.append('x_'+str(i+1)+'_'+str(j+1)+'_'+str(k+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t+1)) #\sum_{k=1}^{W}\sum_{l=1}^{H}\sum_{n=1}^{N}x_{i,j,k,l,n,t}
                                coefs.append(1.0)

                    expr.append(cplex.SparsePair(cols,coefs))  

    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model


#------------------------------------------------------------------------------------------------------------------#
# Constraint 10: garante que em cada período de tempo um contêiner seja retirado do pátio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I1(model,omega,N,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for t in range(T[o]):
            rest_names.append('restI1_'+str(t+1)+'_'+str(o+1))
            rhs.append(t)
            cols = []
            coefs = []
            for n in omega[o]:          
                #add v coeficientes
                cols.append('v_'+str(n)+'_'+str(t+1))
                coefs.append(1.0)
                        
            expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model



#------------------------------------------------------------------------------------------------------------------#
# Constraint 11: define a variável $v_{nt}$. Quando um contêiner $n$ é retirado do pátio, a variável $v_{nt}$ se 
# torna 1 e se mantém igual a 1 nos períodos de tempo seguintes:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I2(model,omega,N,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for n in omega[o]: 
            for t in range(T[o]-1):
                rest_names.append('restI2_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(0.0)     
                    #add v coeficientes
                cols = ['v_'+str(n)+'_'+str(t+1)]
                coefs = [1.0]
                cols.append('v_'+str(n)+'_'+str(t+2))
                coefs.append(-1.0)
                        
                expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model
#------------------------------------------------------------------------------------------------------------------#
# Constraint 12: garante que o contêiner $n$ seja carregado no navio no período de tempo $t$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I3(model,omega,N,R,C,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for n in omega[o]: 
            for t in range(T[o]-1):
                rest_names.append('restI3_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(0.0)     
                    #add v coeficientes
                cols = ['v_'+str(n)+'_'+str(t+2)]
                coefs = [-1.0]
                for r in range(R):
                    for c in range(C):
                        cols.append('z_'+str(n)+'_'+str(t+1)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(1.0)
                    
                expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 13: assegura que uma posição $(r,c)$ no navio só pode ser ocupada por um contêiner, seja ele um 
# contêiner que foi carregado no porto atual (porto $o$), em algum porto anterior (porto $o-1$) ou um contêiner
# que já estava no navio e está sendo remanejado em $o$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I4(model,omega,N,R,C,T,P):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for r in range(R):
            for c in range(C):
                for t in range(T[o]):
                    rest_names.append('restI4_'+str(r+1)+'_'+str(c+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(1.0)
                    cols = [] 
                    coefs = []
                    
                    for oo in range(o):
                        for d in range(o+1,P):
                            for a in range(o+1,d+1):
                                cols.append('w_'+str(oo+1)+'_'+str(d+1)+'_'+str(a+1)+'_'+str(r+1)+'_'+str(c+1))
                                coefs.append(1.0)
                    
                    for d in range(o+1,P):
                        cols.append('q_'+str(o+1)+'_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(1.0)                            
                    
                    for n in omega[o]: 
                        cols.append('z_'+str(n)+'_'+str(t+1)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(1.0)
                    
                    expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 14: certifica que o contêiner $n$, depois de carregado, não mude de posição enquanto o navio estiver
# parado no mesmo porto:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I5(model,omega,N,R,C,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for n in omega[o]: 
            for t in range(T[o]-1):
                for r in range(R):
                    for c in range(C):
                        rest_names.append('restI5_'+str(n)+'_'+str(t+1)+'_'+str(r+1)+'_'+str(c+1)+'_'+str(o+1))
                        rhs.append(0.0)     
                       #add z coeficientes
                        cols=['z_'+str(n)+'_'+str(t+1)+'_'+str(r+1)+'_'+str(c+1)]
                        coefs = [1.0]
                        cols.append('z_'+str(n)+'_'+str(t+2)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(-1.0)
                        
                        expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 15:  garante que se há um contêiner na posição $(r,c)$ do navio, ele deve ser um contêiner que acabou 
# de ser embarcado, ou um contêiner de remanejamento:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I6(model,phi,R,C,T,P,N):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):    
          for r in range(R):
              for c in range(C):
                  for d in range(o+1,P):
                      rest_names.append('restI6_'+str(r+1)+'_'+str(c+1)+'_'+str(o+1)+'_'+str(d+1))
                      rhs.append(0.0)     
                      #add q coeficientes
                      cols = ['q_'+str(o+1)+'_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1)]
                      coefs = [1.0]
                      for n in phi[o][d]:
                          cols.append('z_'+str(n)+'_'+str(T[o])+'_'+str(r+1)+'_'+str(c+1))
                          coefs.append(1.0)                       
                      for a in range(o+1,d+1):
                          cols.append('w_'+str(o+1)+'_'+str(d+1)+'_'+str(a+1)+'_'+str(r+1)+'_'+str(c+1))
                          coefs.append(-1.0)                   
                        
                      expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 16: assegura que todos os $N_{o}$ contêineres do pátio $o$ já foram embarcados no navio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I7(model,omega,N,R,C,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):    
          for n in omega[o]:
              for t in range(T[o]-1):
                      rest_names.append('restI7_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                      rhs.append(1.0)
                      cols = [] 
                      coefs = []
                      for r in range(R):
                          for c in range(C):
                              cols.append('z_'+str(n)+'_'+str(T[o])+'_'+str(r+1)+'_'+str(c+1))
                              coefs.append(1.0)                                     
                        
                      expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 17:garante que, durante o processo de carregamento do navio, nenhum contêiner seja alocado em uma 
# posição flutuante ou que ocupe a posição de um contêiner que já estava no navio ou foi remanejado:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I8(model,omega,N,R,C,T,P):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for r in range(R-1):
            for c in range(C):
                for t in range(T[o]):
                    rest_names.append('restI8_'+str(r+1)+'_'+str(c+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(0.0)
                    cols = [] 
                    coefs = []
                    for n in omega[o]: 
                        cols.append('z_'+str(n)+'_'+str(t+1)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(-1.0)
                        cols.append('z_'+str(n)+'_'+str(t+1)+'_'+str(r+2)+'_'+str(c+1))
                        coefs.append(1.0)                      
                    
                    for oo in range(o):
                        for d in range(o+1,P):
                            for a in range(o+1,d+1):
                                cols.append('w_'+str(oo+1)+'_'+str(d+1)+'_'+str(a+1)+'_'+str(r+1)+'_'+str(c+1))
                                coefs.append(-1.0)
                    
                    for d in range(o+1,P):
                        cols.append('q_'+str(o+1)+'_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(-1.0)                            
                    
                    expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 18: contabiliza o número total de contêineres que foram remanejados no porto $o$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I9(model,omega,N,R,C,P):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for d in range(o+1,P):
            rest_names.append('restI9_'+str(o+1)+'_'+str(d+1))
            rhs.append(0.0)
            cols = [] 
            coefs = []
            for oo in range(o):
                for r in range(R):
                    for c in range(C):
                        cols.append('w_'+str(oo+1)+'_'+str(d+1)+'_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(1.0)
            for r in range(R):
                for c in range(C):
                    cols.append('q_'+str(o+1)+'_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1))
                    coefs.append(-1.0)
                            
                    
            expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 19: mantém a estabilidade do navio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I10(model,omega,N,R,C,TT):
    
    rhs = []
    rest_names = []
    expr = []
    
    v=np.sum(TT,axis=0)
    theta=[0.0]*len(N)
    theta[0]= N[0]-v[0]
    
    for i in range(1,len(N)):
        theta[i]= theta[i-1]+ N[i]-v[i]
        
    
    for o in range(len(N)):
        temp = int(np.ceil(float(theta[o])/float(C)))
        
        if temp < R :
            rest_names.append('restI10_'+str(o+1))
            rhs.append(0.0)
            cols = [] 
            coefs = []
            
            for c in range(C):
                for r in range(temp,R):
                    cols.append('u_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1))
                    coefs.append(1.0)
            
            expr.append(cplex.SparsePair(cols,coefs)) 
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 20: restrição de conservação de fluxo, e indica que o número total de contêineres no porto $o$ deve ser
# igual ao número de contêineres que foram embarcados nos portos $p=1,2,...,o$ menos os contêineres que foram 
# desembarcados nos portos $p=1,2,...,o$ :
#------------------------------------------------------------------------------------------------------------------#
def restricao_N1(model,P,R,C,TT):   
    rhs = []
    rest_names = []
    expr = []
    for o in range(P-1):
        for d in range(o+1,P):
            rest_names.append('restN1_'+str(o+1)+'_'+str(d+1))
            rhs.append(float(TT[o,d]))
            cols = [] 
            coefs = []
            for a in range(o+1,d+1):        
                for r in range(R):
                    for c in range(C):
                        cols.append('w_'+str(o+1)+'_'+str(d+1)+'_'+str(a+1)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(1.0)
            
            for m in range(o):
                for r in range(R):
                    for c in range(C):            
                        cols.append('w_'+str(m+1)+'_'+str(d+1)+'_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(-1.0)                                
                    
            expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model  
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 21: garante que cada posição $(r, c)$ tenha no máximo um único contêiner:
#------------------------------------------------------------------------------------------------------------------#
def restricao_N2(model,R,C,P):   
    rhs = []
    rest_names = []
    expr = []
    for o in range(P-1):
        for r in range(R):
            for c in range(C):
                rest_names.append('restN2_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1))
                rhs.append(0.0)
                cols = ['u_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1)]
                coefs = [-1.0]
                for m in range(o+1):
                    for d in range(o+1,P):
                        for a in range(o+1,d+1):
                            cols.append('w_'+str(m+1)+'_'+str(d+1)+'_'+str(a+1)+'_'+str(r+1)+'_'+str(c+1))
                            coefs.append(1.0)                                                
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model      
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 22: garante que existem contêineres embaixo do contêiner que ocupa a célula $(r, c)$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_N3(model,N,R,C):   
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for r in range(R-1):
            for c in range(C):
                rest_names.append('restN3_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1))
                rhs.append(0.0)
                cols = ['u_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1)]
                coefs = [-1.0]
                cols.append('u_'+str(o+1)+'_'+str(r+2)+'_'+str(c+1))
                coefs.append(1.0)
                                              
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model      
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 23: é responsável por definir como um contêiner pode ser desembarcado no porto $d$ ao impor que se
# um contêiner que ocupa a posição $(r, c)$, então ele será desembarcado no porto $d$, se não houver um contêiner
# na posição $(r+1, c)$ acima dele:
#------------------------------------------------------------------------------------------------------------------#
def restricao_N4(model,P,R,C):   
    rhs = []
    rest_names = []
    expr = []
    
    for d in range(1,P):
        for r in range(R-1):
            for c in range(C):
                rest_names.append('restN4_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1))
                rhs.append(1.0)
                cols = [] 
                coefs = []
                for o in range(d):
                    for e in range(d,P):
                        cols.append('w_'+str(o+1)+'_'+str(e+1)+'_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1))
                        coefs.append(1.0)
                    for u in range(d+1,P):
                        for a in range(d+1,u+1):
                            cols.append('w_'+str(o+1)+'_'+str(u+1)+'_'+str(a+1)+'_'+str(r+2)+'_'+str(c+1))
                            coefs.append(1.0)           
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model    

if __name__ == "__main__":
     model1 = ModeloIntegrado('InstanciaModeloIntegrado_I.mat')
     #model2 = cplex.Cplex("model_Matlab_Inst_II.lp")
     #compare_models(model1,model2)