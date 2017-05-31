import numpy as np
from mpmath import mp
from math import sqrt

def Jrange(J):
	return [x-J+0.0 for x in range(int(2*J+1))]

def J2range(J):
	return [J*(J+1)]*int(2*J+1)
	
def Jminusrange(J):
	return [x-J+0.0 for x in range(1,int(2*J+1))]

def Jplusrange(J):
	return [x-J+0.0 for x in range(0,int(2*J))]

def Jdiagrange(J,k):
	if k > 0: # minus
		return [x-J+0.0 for x in range(k,int(2*J+1))]
	if k < 0:
		return [x-J+0.0 for x in range(0,int(2*J+1+k))]
	return Jrange(J)

def Jz(J):
	return mp.diag(Jrange(J))
	
def Jplus(J):
	return mp.matrix(np.diag([sqrt((J-Jz)*(J+Jz+1)) for Jz in Jplusrange(J)],-1))

def Jminus(J):
	return mp.matrix(np.diag([sqrt((J+Jz)*(J-Jz+1)) for Jz in Jminusrange(J)],1))

def Jx(J):
	return 0.5*(Jplus(J) + Jminus(J))
	
def Jy(J):
	return -0.5j*(Jplus(J) - Jminus(J))

def J2(J):
	return mp.diag(J2range(J))
	
if __name__ == "__main__":
	import numpy as np
	def print_np_matrix(m):
		for xs in np.array(m):
			print(" ".join("{:5.2g}".format(xxs) for xxs in xs))

	print("J+ =")
	print_np_matrix(Jplus(4))
	print("J- =")
	print_np_matrix(Jminus(4))
	print("Jz =")
	print_np_matrix(Jz(4))
	print("J-^3 =")
	print_np_matrix(Jminus(4)**3)
	print("J-^3 * Jz = ")
	print_np_matrix(Jminus(4)**3*Jz(4))

