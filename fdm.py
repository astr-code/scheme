import numpy as np
from fractions import Fraction
import scheme_class as sc
import taylor as ty
import plot_scheme as ps

# def Factorial(n): # return factorial
#     result = 1
#     for i in range (1,n+1):
#         result = result * i
#     # print( "factorial is ",result)
#     return result


# the function is to derive compact scheme by giving the order of accuracy
def compact_scheme_analysis(order_of_accuracy):

	print(" The order of accuracy:",order_of_accuracy,type(order_of_accuracy))

	a = np.zeros([order_of_accuracy+1,order_of_accuracy+1], dtype = float) 
	b = np.zeros([order_of_accuracy+1], dtype = float) 
	c = np.zeros([order_of_accuracy+1], dtype = float) 

	lhs = np.zeros([3], dtype = float) 
	rhs = np.zeros([order_of_accuracy-1], dtype = float) 

	c[1] = 1.0

	# print(c)

	# print(' Shape of array a: ',np.shape(a))

	first_node_stencil = int(-(order_of_accuracy+1)/2+1)
	last_node_stencil  = int(first_node_stencil + order_of_accuracy - 2)

	print("{}{}{}{}{}".format(" The stencil is from: [i",first_node_stencil," ~ i+",last_node_stencil,"]"))

	j = -1

	for n in range(first_node_stencil,last_node_stencil+1):

		j = j + 1

		a[:,j] = ty.TaylorSeriesCoef(n,order_of_accuracy ,'f')

		# print(n,':',a[:,j])
	
	for n in range(-1, 2, 2):

		j = j+1

		a[:,j] = ty.TaylorSeriesCoef(n,order_of_accuracy ,'df')
		
		# print(n,':',a[:,j])

	b = np.matmul(np.linalg.inv(a), c)

	rhs = b[0:order_of_accuracy-1]
	lhs[0]=-b[order_of_accuracy-1]
	lhs[1]=1.0
	lhs[2]=-b[order_of_accuracy]
	# print(rhs,lrh)
	# print(Fraction(lrh[0]).limit_denominator())
	print(' scheme of',order_of_accuracy,'th order of accuracy: ')
	
	scheme=sc.fdm_compact_tridiagonal(first_node_stencil,last_node_stencil,lhs,rhs)
	scheme.display()
	scheme.spectra_property()

	ps.plot_wavenumber(scheme.wavenumber,scheme.modified_wavenumber_real,scheme.modified_wavenumber_imag)


	# np.savetxt('myfile.dat', np.c_[scheme.wavenumber,scheme.modified_wavenumber_real,scheme.modified_wavenumber_imag])


	# WaveMK=spectral_analyis(rhs,lrh)


	
