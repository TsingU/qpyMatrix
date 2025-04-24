import numpy as np
class Matrix:#矩阵对象(Matrix objects)
    def __init__(self,m,n,p=0):#初始化矩阵(Initialize the matrix)
        self.row=m
        self.col=n
        self.matrix=[]
        for i in range(m):
            tmp=[]
            for j in range(n):
                tmp.append(p)
            self.matrix.append(tmp)
    #运算符重载(2025-04-24更新)
    def __add__(self,other):#重载加法+为矩阵加法
        return plus(self,other)
    def __mul__(self,other):#重载乘法*为矩阵乘法
        return multiply(self,other)
    def __pow__(self,other):#重载**为张量积
        return kron(self,other)
    
    def show(self):#打印矩阵(Print the matrix)
        print()
        for i in range(self.row):
            print(self.matrix[i])
        print()
    def element(self,i,j,p,rw=False):#访问矩阵元素，包括修改元素值(Access matrix elements, including modifying element values)
        #rw为控制变量，如果为True,则修改元素值；否则直接打印元素值
        #(rw is the control variable, if it is True, the element value is modified; Otherwise, the element values are printed directly)
        if rw==True:
            self.matrix[i][j]=p
        else:
            print(self.matrix[i][j])

class qProg:#量子线路
    def __init__(self):#初始化量子线路
        self.prog=[]
    def add(self,qprog):#向线路中添加操作
        self.prog.append(qprog)
    def draw(self):#绘制线路图
        for i in range(len(self.prog[0])):
            s=""
            for j in range(len(self.prog)):
                s+="--\t"+self.prog[j][i]+"\t"
            print("q"+str(i)+'\t'+s+"--")
    def help(self):
        print("\t量子线路门符号")
        print("\t--------------------\t")
        print("I:\t代表不操作")
        print("X:\t代表X门")
        print("Y:\t代表Y门")
        print("Z:\t代表Z门")
        print("H:\t代表H门")
        print("S:\t代表S门")
        print("T:\t代表qT门")
        print("SWAP:\t代表SWAP门\t用法:add([\"SWAP\",\"SWAP\"])")
        print("iSWAP:\t代表iSWAP门\t用法:add([\"iSWAP\",\"iSWAP\"])")
        print("*,+:\t代表CNOT门\t用法:add([\"*\",\"+\"])")
        print("+,*:\t代表CNOT_H门\t用法:add([\"+\",\"*\"])")
        print("CZ:\t代表CZ门\t用法:add([\"CZ\",\"CZ\"])")
        print("\tRX,RY,RZ门开发中")
        print("\t--------------------\t")
    def subMatrix(self):#计算线路的合并矩阵
        result=diag(4,1)
        for i in range(len(self.prog))[::-1]:
            if self.prog[i]==["X","I"]:#X-I
                result=multiply(result,kron(I,X))
            elif self.prog[i]==["I","X"]:#I-X
                result=multiply(result,kron(X,I))
            elif self.prog[i]==["Y","I"]:#Y-I
                result=multiply(result,kron(I,Y))
            elif self.prog[i]==["I","Y"]:#I-Y
                result=multiply(result,kron(Y,I))
            elif self.prog[i]==["Z","I"]:#Z-I
                result=multiply(result,kron(I,Z))
            elif self.prog[i]==["I","Z"]:#I-Z
                result=multiply(result,kron(Z,I))
            elif self.prog[i]==["H","I"]:#H-I
                result=multiply(result,kron(I,H))
            elif self.prog[i]==["I","H"]:#I-H
                result=multiply(result,kron(H,I))
            elif self.prog[i]==["S","I"]:#S-I
                result=multiply(result,kron(I,S))
            elif self.prog[i]==["I","S"]:#I-S
                result=multiply(result,kron(S,I))
            elif self.prog[i]==["T","I"]:#T-I
                result=multiply(result,kron(I,qT))
            elif self.prog[i]==["I","T"]:#I-T
                result=multiply(result,kron(qT,I))
            elif self.prog[i]==["I","I"]:#I-I
                result=multiply(result,kron(I,I))
            elif self.prog[i]==["X","X"]:#X-X
                result=multiply(result,kron(X,X))
            elif self.prog[i]==["Y","Y"]:#Y-Y
                result=multiply(result,kron(Y,Y))
            elif self.prog[i]==["Z","Z"]:#Z-Z
                result=multiply(result,kron(Z,Z))
            elif self.prog[i]==["H","H"]:#H-H
                result=multiply(result,kron(H,H))
            elif self.prog[i]==["S","S"]:#S-S
                result=multiply(result,kron(S,S))
            elif self.prog[i]==["T","T"]:#T-T
                result=multiply(result,kron(qT,qT))
            elif self.prog[i]==["*","+"]:
                result=multiply(result,CNOT)
            elif self.prog[i]==["+","*"]:
                result=multiply(result,CNOT_H)
            elif self.prog[i]==["SWAP","SWAP"]:
                result=multiply(result,SWAP)
            elif self.prog[i]==["iSWAP","iSWAP"]:
                result=multiply(result,iSWAP)
            elif self.prog[i]==["CZ","CZ"]:
                result=multiply(result,CZ)
            ###
            elif self.prog[i]==["X","Y"]:#X-Y
                result=multiply(result,kron(Y,X))
            elif self.prog[i]==["X","Z"]:#X-Z
                result=multiply(result,kron(Z,X))
            elif self.prog[i]==["X","H"]:#X-H
                result=multiply(result,kron(H,X))
            elif self.prog[i]==["X","S"]:#X-S
                result=multiply(result,kron(S,X))
            elif self.prog[i]==["X","T"]:#X-T
                result=multiply(result,kron(qT,X))
            ###
            elif self.prog[i]==["Y","X"]:#Y-X
                result=multiply(result,kron(X,Y))
            elif self.prog[i]==["Y","Z"]:#Y-Z
                result=multiply(result,kron(Z,Y))
            elif self.prog[i]==["Y","H"]:#Y-H
                result=multiply(result,kron(H,Y))
            elif self.prog[i]==["Y","S"]:#Y-S
                result=multiply(result,kron(S,Y))
            elif self.prog[i]==["Y","T"]:#Y-T
                result=multiply(result,kron(qT,Y))
            ###
            elif self.prog[i]==["Z","X"]:#Z-X
                result=multiply(result,kron(X,Z))
            elif self.prog[i]==["Z","Y"]:#Z-Y
                result=multiply(result,kron(Y,Z))
            elif self.prog[i]==["Z","H"]:#Z-H
                result=multiply(result,kron(H,Z))
            elif self.prog[i]==["Z","S"]:#Z-S
                result=multiply(result,kron(S,Z))
            elif self.prog[i]==["Z","T"]:#Z-T
                result=multiply(result,kron(qT,Z))
            ###
            elif self.prog[i]==["H","X"]:#H-X
                result=multiply(result,kron(X,H))
            elif self.prog[i]==["H","Y"]:#H-Y
                result=multiply(result,kron(Y,H))
            elif self.prog[i]==["H","Z"]:#H-Z
                result=multiply(result,kron(Z,H))
            elif self.prog[i]==["H","S"]:#H-S
                result=multiply(result,kron(S,H))
            elif self.prog[i]==["H","T"]:#H-T
                result=multiply(result,kron(qT,H))
            ###
            elif self.prog[i]==["S","X"]:#S-X
                result=multiply(result,kron(X,S))
            elif self.prog[i]==["S","Y"]:#S-Y
                result=multiply(result,kron(Y,S))
            elif self.prog[i]==["S","Z"]:#S-Z
                result=multiply(result,kron(Z,S))
            elif self.prog[i]==["S","H"]:#S-H
                result=multiply(result,kron(H,S))
            elif self.prog[i]==["S","T"]:#S-T
                result=multiply(result,kron(qT,S))
            ###
            elif self.prog[i]==["T","X"]:#T-X
                result=multiply(result,kron(X,T))
            elif self.prog[i]==["T","Y"]:#T-Y
                result=multiply(result,kron(Y,T))
            elif self.prog[i]==["T","Z"]:#T-Z
                result=multiply(result,kron(Z,T))
            elif self.prog[i]==["T","H"]:#T-H
                result=multiply(result,kron(H,T))
            elif self.prog[i]==["T","S"]:#T-S
                result=multiply(result,kron(S,qT))
            ###
                #RX,RY,RZ门开发中···#
            ###
        return result




def plus(A,B):#矩阵加法(Matrix addition)
        if A.row!=B.row or A.col!=B.col:
            print("矩阵不符合矩阵加法法则")
        else:
            for i in range(A.row):
                for j in range(A.col):
                    A.matrix[i][j]=A.matrix[i][j]+B.matrix[i][j]
        return A

def cmul(c,A):#矩阵数乘(2025-04-24更新)
    result=Matrix(A.row,A.col)
    for i in range(A.row):
        for j in range(A.col):
            result.matrix[i][j]=A.matrix[i][j]*c
    return result

def multiply(A,B):#矩阵乘法(Matrix multiplication)
    if A.col!=B.row:
        print("矩阵不符合矩阵乘法法则")
    else:
        result=Matrix(A.row,B.col)
        for i in range(A.row):
            for j in range(B.col):
                s=0
                for k in range(B.row):
                    s=s+A.matrix[i][k]*B.matrix[k][j]
                result.element(i,j,s,True)
        return result

def det(A):#计算矩阵行列式(Calculate the matrix determinant)
    return np.linalg.det(A)

def inverse(A):#求逆矩阵(Find the inverse matrix)
    if A.row!=A.col:
        print("矩阵不是方阵")
    else:
        result=Matrix(A.row,A.col)
        r=np.linalg.inv(A.matrix)
        for i in range(A.row):
            for j in range(A.col):
                result.element(i,j,r[i][j],True)
        return result

def T(A):#计算转置矩阵(Calculate the transpose matrix)
    result=Matrix(A.col,A.row)
    for i in range(A.row):
        for j in range(A.col):
            result.element(j,i,A.matrix[i][j],True)
    return result

def diag(n,p):#构造对角矩阵(Construct diagonal matrices)
    result=Matrix(n,n)
    for i in range(n):
        result.matrix[i][i]=p
    return result

def kron(A,B):#计算矩阵的张量积(Calculate the tensor product of the matrix)
    result=Matrix(A.row*B.row,A.col*B.col)
    for i in range(A.row):
        for j in range(A.col):
            a_ij=A.matrix[i][j]
            s=cmul(a_ij,B)
            for m in range(B.row):
                for n in range(B.col):
                    result.matrix[m+(i)*B.row][n+(j)*B.col]=s.matrix[m][n]
    return result
########################################## begin 量子计算预设矩阵(Quantum computing preset matrices)

I=diag(2,1)#单位矩阵

X=Matrix(2,2)#泡利X门(Pauli X-Gate)
X.matrix=[[0,1],[1,0]]

Y=Matrix(2,2)#泡利Y门(Pauli Y-Gate)
Y.matrix[0][1]=-1j
Y.matrix[1][0]=1j

Z=diag(2,1)#泡利Z门(Pauli Z-Gate)
Z.matrix[1][1]=-1

H=Matrix(2,2)#阿达玛门(Hadamard Gate)
H.matrix=[[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]]

S=Matrix(2,2)#S门/相位门(S-Gate/Phase Gate)
S.matrix=[[1,0],[0,1j]]

qT=Matrix(2,2)#T门(T-Gate)
qT.matrix=[[1,0],[0,(1+1j)/np.sqrt(2)]]

CNOT=Matrix(4,4)#CNOT门,低位比特为控制比特(CNOT gate, the low bit is the control bit)
CNOT.matrix=[[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]]

CNOT_H=Matrix(4,4)#CNOT门，高位比特为控制比特(CNOT gate, the high bit is the control bit)
CNOT_H.matrix=[[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]

SWAP=Matrix(4,4)#SWAP门，交换两个比特的状态(SWAP gate, which swaps the state of two bits)
SWAP.matrix=[[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]

iSWAP=Matrix(4,4)#SWAP门，增加一个相对相位，如果两个比特的状态不同，则相位发生变化(SWAP gate, adding a relative phase, if the state of the two bits is different, the phase changes)
iSWAP.matrix=[[1,0,0,0],[0,0,1j,0],[0,1j,0,0],[0,0,0,1]]

CZ=Matrix(4,4)#控制Z门(Controlled Z-Gate)
CZ.matrix=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]]

CCNOT=Matrix(8,8)#Toffoli门当两个高位控制比特均为1时翻转低位比特(Toffoli Gate, flips the low bit when both high control bits are 1)
CCNOT.matrix=[[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0]]

########################################## end 量子计算预设矩阵(Quantum computing preset matrices)
#----------------------------------------#
########################################## begin 量子计算预设函数(Quantum computing preset functions)

def RX(p):#旋转X门(Revolving X-Gate)
    rx=Matrix(2,2)
    rx.matrix=[[np.cos(p/2),-1j*np.sin(p/2)],[-1j*np.sin(p/2),np.cos(p/2)]]
    return rx

def RY(p):#旋转Y门(Revolving Y-Gate)
    ry=Matrix(2,2)
    ry.matrix=[[np.cos(p/2),-np.sin(p/2)],[np.sin(p/2),np.cos(p/2)]]
    return ry

def RZ(p):#旋转Z门(Revolving Z-Gate)
    rz=Matrix(2,2)
    rz.matrix=[[np.cos(p/2)-1j*np.sin(p/2),0],[0,np.cos(p/2)+1j*np.sin(p/2)]]

def CT(A):#共轭转置(Conjugate transpose)
    result=Matrix(A.col,A.row)
    for i in range(A.row):
        for j in range(A.col):
            s=A.matrix[i][j]
            result.element(j,i,s.conjugate(),True)
    return result

def Measure(m):#测量函数，m为测量基(Measurement function, m is the measurement base)
    result=multiply(m,CT(m))
    return result

def measure_p(m,qubit):#测量量子比特qubit在测量基m上的概率(Measure the probability of qubit qubit on the measurement basis m)
    M=Measure(m)
    s=multiply(CT(qubit),M)
    result=multiply(s,qubit)
    return result.matrix[0][0]

def bloch(v):#返回一个量子比特在布洛赫球上的z坐标(2025-04-24更新)
    s=measure_p(v,v_zero)
    z=2*s-1
    return z


########################################## end 量子计算预设函数(Quantum computing preset functions)
#----------------------------------------#
########################################## begin 量子计算预设基向量与测量基(2025-04-24更新)

v_zero=Matrix(2,1)#基态|0>
v_zero.matrix=[[1],[0]]

v_one=Matrix(2,1)#基态|1>
v_one.matrix=[[0],[1]]

v_currect=Matrix(2,1)#基态|+>
v_currect.matrix=[[1/np.sqrt(2)],[1/np.sqrt(2)]]

v_minus=Matrix(2,1)#基态|->
v_minus.matrix=[[1/np.sqrt(2)],[-1/np.sqrt(2)]]

v_i=Matrix(2,1)#基态|i>
v_i.matrix=[[1/np.sqrt(2)],[1j/np.sqrt(2)]]

v_minus_i=Matrix(2,1)#基态|-i>
v_minus_i.matrix=[[1/np.sqrt(2)],[-1j/np.sqrt(2)]]

M_0=v_zero*CT(v_zero)#|0>的测量基矩阵
M_1=v_one*CT(v_one)#|1>的测量基矩阵
M_current=v_currect*CT(v_currect)#|+>的测量基矩阵
M_minus=v_minus*CT(v_minus)#|->的测量基矩阵
M_i=v_i*CT(v_i)#|i>的测量基矩阵
M_minus_i=v_minus_i*CT(v_minus_i)#|-i>的测量基矩阵

########################################## end 量子计算预设基向量与测量基(2025-04-24更新)