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

def plus(A,B):#矩阵加法(Matrix addition)
        if A.row!=B.row or A.col!=B.col:
            print("矩阵不符合矩阵加法法则")
        else:
            for i in range(A.row):
                for j in range(A.col):
                    A.matrix[i][j]=A.matrix[i][j]+B.matrix[i][j]
        return A

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
            a_ij=diag(B.row,A.matrix[i][j])
            s=multiply(a_ij,B)
            for m in range(B.row):
                for n in range(B.col):
                    result.matrix[m+(i-1)*B.row][n+(j-1)*B.col]=s.matrix[m][n]
    return result
########################################## begin 量子计算预设矩阵(Quantum computing preset matrices)

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
S.matrix=[[0,1],[0,1j]]

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
    result.show()

########################################## end 量子计算预设函数(Quantum computing preset functions)