# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 17:29:50 2018

@author: catar
"""

def compare_models(model1,model2):
    
    
    if model1.variables.get_num() != model2.variables.get_num() :
        print("numero de variaveis nao sao iguais: model1[ ",model1.variables.get_num()," ] vs model 2[ ",model2.variables.get_num()," ]")
    
    namesv1 = set(model1.variables.get_names())
    namesv2 = set(model2.variables.get_names())
    
    if namesv1 != namesv2 :
        print("nomes de variaveis nao iguais")
        
        print("variaveis no modelo 1 que nao esta no modelo 2:", namesv1 -namesv2 )
        print("variaveis no modelo 2 que nao esta no modelo 1:", namesv2 -namesv1 )
        
    
    if model1.linear_constraints.get_num() != model2.linear_constraints.get_num() :
        print("numero de restricoes nao sao iguais: model1[ ",model1.linear_constraints.get_num()," ] vs model 2[ ",model2.linear_constraints.get_num()," ]")
    
    namesc1 = set(model1.linear_constraints.get_names())
    namesc2 = set(model2.linear_constraints.get_names())
    
    if namesc1 != namesc2 :
        print("nomes de retricoes nao iguais")
        
        print("restricoes no modelo 1 que nao esta no modelo 2:", namesc1 -namesc2 )
        print("restricoes no modelo 2 que nao esta no modelo 1:", namesc2 -namesc1 )
        
    
    for namec in namesc1:
        
        #rhs
        rhs1 = model1.linear_constraints.get_rhs(namec)
        rhs2 = model2.linear_constraints.get_rhs(namec)
        
        if rhs1 != rhs2 :
            print("Diferente rhs em restricao: ",namec,", rhs no modelo 1 = ",rhs1,", rhs no modelo2 = ",rhs2)
        
        #senses
        sense1 = model1.linear_constraints.get_senses(namec)
        sense2 = model2.linear_constraints.get_senses(namec)
        
        if sense1 != sense2 :
            print("Diferente sentido em restricao: ",namec,", sentido no modelo 1 = ",sense1,", sentido no modelo2 = ",sense2)
        
        #expr
        expr1 = model1.linear_constraints.get_rows(namec)
        expr2 = model2.linear_constraints.get_rows(namec)
        
        inds1, vals1 = expr1.unpack()
        inds2, vals2 = expr2.unpack()
        
        names1 = model1.variables.get_names(inds1)
        names1s = set(names1)
        
        names2 = model2.variables.get_names(inds2)
        names2s = set(names2)
        
        if names1s != names2s :
            print(" restricao: ",namec, " diferentes variaveis")
            print("modelo 1 tem e modelo 2 nao tem :", names1s -names2s )
            print("modelo 2 tem e modelo 1 nao tem :" ,names2s -names1s )
            
        for c in range(len(inds1)) :
            vname1 = names1[c]
            coeff1 = vals1[c]
            
            c2 = names2.index(vname1)
            coeff2 = vals2[c2]
            
            if coeff1 != coeff2 :
                print("coeficente de variavel ",vname1, " na restricao ", namec, " eh differente")
                print("no modelo1 eh: ",coeff1," , no modelo 2 eh: ",coeff2)
            
        
        
    return