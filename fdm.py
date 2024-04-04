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
def scheme_analysis(order_of_accuracy,form):

	print(" The order of accuracy:",order_of_accuracy,type(order_of_accuracy))
	print(" The scheme is :",form)

	a = np.zeros([order_of_accuracy+1,order_of_accuracy+1], dtype = float) 
	b = np.zeros([order_of_accuracy+1], dtype = float) 
	c = np.zeros([order_of_accuracy+1], dtype = float) 

	if form == 'compact':
		lhs = np.zeros([3], dtype = float) 
		rhs = np.zeros([order_of_accuracy-1], dtype = float)
		first_node_stencil = int(-(order_of_accuracy+1)/2+1)
		last_node_stencil  = int(first_node_stencil + order_of_accuracy - 2)
	elif form == 'explicit':
		lhs = np.zeros([1], dtype = float) 
		rhs = np.zeros([order_of_accuracy+1], dtype = float)
		first_node_stencil = int(-(order_of_accuracy+1)/2)
		last_node_stencil  = int(first_node_stencil + order_of_accuracy)

	c[1] = 1.0

	print("{}{}{}{}{}".format(" The stencil is from: [i",first_node_stencil," ~ i+",last_node_stencil,"]"))

	j = -1

	for n in range(first_node_stencil,last_node_stencil+1):

		j = j + 1

		a[:,j] = ty.TaylorSeriesCoef(n,order_of_accuracy ,'f')

		# print(n,':',a[:,j])

	if form == 'compact':

		for n in range(-1, 2, 2):
			j = j+1
			a[:,j] = ty.TaylorSeriesCoef(n,order_of_accuracy ,'df')
			# print(n,':',a[:,j])

	b = np.matmul(np.linalg.inv(a), c)

	rhs = b[0:order_of_accuracy-1]

	if form == 'compact':
		rhs = b[0:order_of_accuracy-1]
		lhs[0]=-b[order_of_accuracy-1]
		lhs[1]=1.0
		lhs[2]=-b[order_of_accuracy]
	elif form == 'explicit':
		rhs = b[0:order_of_accuracy+1]
		lhs[0]=1.0
	# print(rhs,lrh)
	# print(Fraction(lrh[0]).limit_denominator())
	print(rhs)
	print(' scheme of',order_of_accuracy,'th order of accuracy: ')
	
	scheme=sc.fdm_scheme(first_node_stencil,last_node_stencil,lhs,rhs,form)
	scheme.display()
	scheme.spectra_property()

	filename='modified_wavenumber_stanadard_'+form+'_'+str(order_of_accuracy)+'.dat'
	ps.file_write_wavenumber(filename,scheme.wavenumber,scheme.modified_wavenumber_real,scheme.modified_wavenumber_imag)
	ps.plot_wavenumber(scheme.wavenumber,scheme.modified_wavenumber_real,scheme.modified_wavenumber_imag)


	# np.savetxt('myfile.dat', np.c_[scheme.wavenumber,scheme.modified_wavenumber_real,scheme.modified_wavenumber_imag])


	# WaveMK=spectral_analyis(rhs,lrh)


	
