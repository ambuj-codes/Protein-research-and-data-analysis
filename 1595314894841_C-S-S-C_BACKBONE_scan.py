from math import *
import math
import re
import numpy as np
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

def getlen(atm1,atm2):
    dist=sqrt(pow(atm1.x-atm2.x,2)+pow(atm1.y-atm2.y,2)+pow(atm1.z-atm2.z,2)) 
    return dist

def getangle(atm1,atm2,atm3):
    dist1=sqrt(pow(atm1.x-atm2.x,2)+pow(atm1.y-atm2.y,2)+pow(atm1.z-atm2.z,2)) 
    dist2=sqrt(pow(atm3.x-atm2.x,2)+pow(atm3.y-atm2.y,2)+pow(atm3.z-atm2.z,2)) 
    dotp=(atm1.x-atm2.x)*(atm3.x-atm2.x)+(atm1.y-atm2.y)*(atm3.y-atm2.y)+(atm1.z-atm2.z)*(atm3.z-atm2.z) 
    angle=acos(dotp/(dist1*dist2))*180/pi 
    return angle

def getangledihedral(atm1,atm2,atm3,atm4):
    ab=np.zeros(3)
    bc=np.zeros(3)
    cd=np.zeros(3)
    p=[]
    q=[]
    ab[0]=atm2.x-atm1.x
    ab[1]=atm2.y-atm1.y
    ab[2]=atm2.z-atm1.z
    bc[0]=atm3.x-atm2.x
    bc[1]=atm3.y-atm2.y
    bc[2]=atm3.z-atm2.z
    cd[0]=atm4.x-atm3.x
    cd[1]=atm4.y-atm3.y
    cd[2]=atm4.z-atm3.z
    p.append(ab[1]*bc[2]-ab[2]*bc[1])
    p.append(ab[2]*bc[0]-ab[0]*bc[2])
    p.append(ab[0]*bc[1]-ab[1]*bc[0])
    q.append(bc[1]*cd[2]-bc[2]*cd[1])
    q.append(bc[2]*cd[0]-bc[0]*cd[2])
    q.append(bc[0]*cd[1]-bc[1]*cd[0])


    r1=0
    r2=0
    dp=0
    dpcd=0
    for i in range(0,3):
        r1 += math.pow(p[i],2)
        r2 += math.pow(q[i],2)
        dp += p[i]*q[i]
        dpcd += p[i]*cd[i]

    dih=(dpcd/abs(dpcd))*math.acos(dp/(math.sqrt(r1)*math.sqrt(r2)))*180/math.pi
    

    return dih

def getdihedralstrain(a1,a2,a3,a4,a5):
    dse=8.37*(1+math.cos(3*a1*math.pi/180))+8.37*(1+math.cos(3*a5*math.pi/180))+4.18*(1+math.cos(3*a2*math.pi/180))+4.18*(1+math.cos(3*a4*math.pi/180))+14.64*(1+math.cos(2*a3*math.pi/180))+2.51*(1+math.cos(3*a3*math.pi/180))
    return dse

s_s_l=1.6
s_s_u=2.5

filetxt=open('filelist.txt') 
txt_lines=filetxt.read().split('\n') 
filetxt.close()
fileout=open('out_C-S-S-C_BACKBONE_scan.txt','w')
f1=open('error_C-S-S-C_scan.txt','w')
intr=[]
lenlines=len(txt_lines)
for ppp in range(lenlines):
    filename=txt_lines[ppp]
    if filename=='':
        continue
    print('%.2f'%((ppp+1)*100.0/(lenlines-1))+'% ('+str(ppp+1)+'/'+str(lenlines-1)+')  Executing for:'+filename)
    file=open(filename,'r')
    lines=file.read().split('\n')
    file.close()
    T=[]
    D=[]
    S=[] 
    C=[]
    SX=[]
    TX=[]
    A=[]
    B=[]
    E=[]
    F=[]
    modelno=[]

 
    try:
        for ln in lines:
            if len(ln)>=6 and (ln[0:4]=='ATOM' or ln[0:6]=='HETATM'):
                atm=atom()
                atm.aid=int(ln[6:11]) 
                atm.atype=ln[12:16].strip() 
                atm.rtype=ln[17:20].strip() 
                atm.chainid=ln[21]
                atm.rid=int(ln[22:26]) 
                atm.x=float(ln[30:38]) 
                atm.y=float(ln[38:46]) 
                atm.z=float(ln[47:54]) 
                atm.model=modelno
                symb=ln[13].strip()
                if atm.atype=='CB' and (modelno==1 or modelno==A or modelno==[])  :
                    if atm.rtype=='CYS' : 
                        C.append(atm)
                        D.append(atm)
                if atm.atype=='SG'and (modelno==1 or modelno==A or modelno==[]) :
                    if atm.rtype=='CYS': 
                        SX.append(atm)
                        TX.append(atm)
                if atm.atype=='CA' and (modelno==1 or modelno==A or modelno==[]) :
                    if atm.rtype=='CYS':
                        B.append(atm)
                        E.append(atm)
                if atm.atype=='N' and (modelno==1 or modelno==A or modelno==[]) :
                    if atm.rtype=='CYS'  :
                        A.append(atm)
                        F.append(atm)
            elif len(ln)>=5 and ln[0:5]=='MODEL':
                modelno=int(ln[12:])

    except:
        f1.write(filename+'\n')


    for k in SX:
        for k1 in SX:
            if k1.chainid==k.chainid:                
                if k1.rid==k.rid and k1.aid!=k.aid :
                    break
        else:
            S.append(k)
    for m in TX:
        for m1 in TX:
            if m1.chainid==m.chainid:
                if m1.rid==m.rid and m1.aid!=m.aid :
                    break
        else:
            T.append(m)
    
    for a in range(len(A)):
        for b in range(len(B)):
            if A[a].rid==B[b].rid:
                for j in range(len(C)):
                    for k in range(len(S)):
                        if C[j].rid==S[k].rid and C[j].rid==B[b].rid and C[j].chainid==B[b].chainid==S[k].chainid==A[a].chainid :
                            for m in range(len(T)):
                                if getlen(S[k],T[m])>=s_s_l and getlen(S[k],T[m])<=s_s_u and S[k].rid<T[m].rid :
                                    for n in range(len(D)):
                                        for e in range(len(E)):
                                            if E[e].rid==D[n].rid:
                                                for f in range(len(F)):
                                                    if D[n].rid==T[m].rid and E[e].rid==F[f].rid and D[n].chainid==T[m].chainid==E[e].chainid==F[f].chainid :
                                                        a1=getangledihedral(A[a],B[b],C[j],S[k])
                                                        a2=getangledihedral(B[b],C[j],S[k],T[m])
                                                        a3=getangledihedral(C[j],S[k],T[m],D[n])
                                                        a4=getangledihedral(S[k],T[m],D[n],E[e])
                                                        a5=getangledihedral(T[m],D[n],E[e],F[f])
                                                        dse=getdihedralstrain(a1,a2,a3,a4,a5)
                                                        intr.append([])
                                                        intr[len(intr)-1].append(filename)                                                     
                                                        intr[len(intr)-1].append(C[j].chainid)
                                                        intr[len(intr)-1].append(C[j].rid)                
                                                        intr[len(intr)-1].append(T[m].rid)
                                                        intr[len(intr)-1].append(T[m].chainid)
                                                        intr[len(intr)-1].append(getlen(C[j],S[k]))        
                                                        intr[len(intr)-1].append(getlen(T[m],S[k]))         
                                                        intr[len(intr)-1].append(getlen(T[m],D[n]))     
                                                        intr[len(intr)-1].append(a1)
                                                        intr[len(intr)-1].append(a2)
                                                        intr[len(intr)-1].append(a3)
                                                        intr[len(intr)-1].append(a4)
                                                        intr[len(intr)-1].append(a5)
                                                        intr[len(intr)-1].append(dse)

    
    C=[]
    T=[]
    D=[]
    S=[]
    SX=[]
    TX=[]
    A=[]
    B=[]
    E=[]
    F=[]
    for line in intr:
        for xxd in line:
            fileout.write(str(xxd))
            fileout.write('\t')
        fileout.write('\n')
    intr=[]
    fileout.close()
    fileout=open('out_C-S-S-C_BACKBONE_scan.txt','a')
fileout.close()
f1.close()
