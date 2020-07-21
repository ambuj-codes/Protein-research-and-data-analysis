# # distace between methyl methyl
# import glob
# with open('input.txt','a') as input:
#     inputfile = glob.glob('*.pdb')
#     for i in inputfile:
#         input.write(i)
#         input.write('\n')
from math import*
import math
import re
class atom:
    
    aid=0    
    atype='' 
    x=0.0  
    y=0.0    
    z=0.0    
    rid=0    
    rtype='' 
    model=[]
    chainid=''

A =[]
B =[]
X=[]
modelno=[]
m= 0.0
N = 3.4
def getlen(atm1,atm2):
    dist=sqrt(pow(atm1.x-atm2.x,2)+pow(atm1.y-atm2.y,2)+pow(atm1.z-atm2.z,2)) 
    return dist
with open('input.txt','r') as file_txt:

    files=file_txt.read().split('\n')

    with open('out_Me...Me_Interaction.txt','w') as output:
        with open('error_Me...Me_Interaction.txt','w') as f1:


            for i in range(len(files)):
    
                filename=files[i]
                if filename=='':
                    continue
                
                print('%.2f'%((i+1)*100.0/(len(files)-1))+'% ('+str(i+1)+'/'+str(len(files)-1)+')  Executing for:'+filename)
                with open(filename,'r') as newfile:
                    lines=newfile.read().split('\n')


                try:
            
                    for line in lines:
            
                        if len(line) >= 6 and ( line[0:4]=='ATOM' or line[0:6]=='HETATM'):
                            atm=atom()
                            atm.rid=int(line[22:26])
                            atm.atype = line[12:16].strip()
                            atm.rtype = line[17:20].strip()
                            atm.aid = int(line[6:11])
                            atm.x = float(line[30:38])
                            atm.y =float(line[38:46])
                            atm.z = float(line[46:54])
                            atm.chainid = line[21]
                            atm.model = modelno                                              
                            if atm.rtype =='ALA':
                                if atm.atype=='CB':          
                                    A.append(atm)
                                    B.append(atm)
                            if atm.rtype=='VAL':
                                if atm.atype=='CG1'or atm.atype=='CG2':
                                    A.append(atm)
                                    B.append(atm)
                            if atm.rtype=='LEU':
                                if atm.atype=='CD1' or atm.atype=='CD2':
                                    A.append(atm)
                                    B.append(atm)
                            if atm.rtype=='ILE':
                                if atm.atype=='CG2' or atm.atype=='CD':    
                                    A.append(atm)
                                    B.append(atm)
                            if atm.rtype=='MET':
                                if atm.atype=='CE':
                                    A.append(atm)
                                    B.append(atm)
                            if atm.rtype=='THR':
                                if atm.atype=='CG2':
                                    A.append(atm)
                                    B.append(atm)
                
                        elif len(line) >=5 and line[0:5]=='MODEL':
                            modelno=int(line[12:]) 
                                                         
                except:
                    f1.write(filename+'\n')
                    print('error')
               
                for a in A:
                    for b in B:

                    
                        
                        if getlen(a,b)<=N and getlen(a,b)>m and a.rid<b.rid:
                    
                        
                

                        
                            a1=getlen(a,b)
                            
                            
                            X.append([])
                            X[len(X)-1].append(filename)
                            X[len(X)-1].append(a.rtype)
                            X[len(X)-1].append(b.rtype)
                            X[len(X)-1].append(a.rid)
                            X[len(X)-1].append(b.rid)
                            X[len(X)-1].append(a.atype)
                            X[len(X)-1].append(b.atype)
                            X[len(X)-1].append(a1)


                        
                                
with open('out_Me...Me_Interaction.txt','a') as output:

    for abc in X:
        for efg in abc:
            output.write(str(efg))
            output.write('\t')
        output.write('\n')
        X=[]  
#     output.close()
#     output=open('out_Me...Me_Interaction.txt','w')
# output.close()
# f1.close()
    


                                     
          




#     with open('output_files.txt','a') as O:
#         R =T.read().split('/n')
#         print(type(R))
                          
# with open('error_file.txt','w') as f1:
#     with open('output_files.txt','a') as O:    
#         for k in X:
#             O.write(str(k))
#             O.write('\n') 
            
            
            
         
                           
                        
                                 
                          
            
               
                    
            

     

