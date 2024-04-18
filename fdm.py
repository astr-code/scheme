import numpy as np
from fractions import Fraction
import scheme_class as sc
import taylor as ty
import plot_scheme as ps

# giving the stencil and the required order of accuracy
def scheme_derive_parameter(first_node_left,last_node_left,first_node_right,last_node_right,order_of_accuracy):
	
	length_of_stencil = last_node_left-first_node_left+1
	lhs_stencil = np.zeros(length_of_stencil,dtype=np.int32)

	length_of_stencil = last_node_right-first_node_right+1
	rhs_stencil = np.zeros(length_of_stencil,dtype=np.int32)

	j=-1
	for i in range(first_node_left,first_node_left+1):
		j=j+1
		lhs_stencil[j] = i
    
	j=-1
	for i in range(first_node_right,last_node_right+1):
		j=j+1
		rhs_stencil[j] = i

	max_order_of_accuray = last_node_left-first_node_left+last_node_right-first_node_right

	print('           LHS stencil:',lhs_stencil)
	print('           RHS stencil:',rhs_stencil)
	print('     Order of accuracy:',order_of_accuracy)
	print(' max Order of accuracy:',max_order_of_accuray)

# the function is to derive the standard scheme by giving the order of accuracy
def scheme_analysis_standard(order_of_accuracy,form):

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


# the function is to derive the standard scheme by giving the stencil
def scheme_analysis_stencil(lhs_stencil,rhs_stencil):

	left_length  = rhs_stencil.size
	right_length = lhs_stencil.size

	order_of_accuracy = right_length + left_length -2

	a = np.zeros([order_of_accuracy+1,order_of_accuracy+1], dtype = float) 
	b = np.zeros([order_of_accuracy+1], dtype = float) 
	c = np.zeros([order_of_accuracy+1], dtype = float) 

	lhs = np.zeros([lhs_stencil.size], dtype = float)
	rhs = np.zeros([rhs_stencil.size], dtype = float)

	first_node_stencil = rhs_stencil[0]
	last_node_stencil  = rhs_stencil[-1]

	c[1] = 1.0


	if (lhs_stencil != 0).any():
		form = 'compact'
	else:
		form = 'explicit'

	# print("{}{}{}{}{}".format(" The stencil is from: [i",first_node_stencil," ~ i+",last_node_stencil,"]"))

	j = -1

	for n in range(first_node_stencil,last_node_stencil+1):

		j = j + 1

		a[:,j] = ty.TaylorSeriesCoef(n,order_of_accuracy ,'f')

		# print(n,':',a[:,j])

	if (lhs_stencil != 0).any(): #compact scheme

		for n in range(lhs_stencil[0], lhs_stencil[-1]+1):

			if n == 0 :
				continue

			j = j+1
			a[:,j] = ty.TaylorSeriesCoef(n,order_of_accuracy ,'df')
			# print(n,':',a[:,j])

	b = np.matmul(np.linalg.inv(a), c)

	rhs = b[0:left_length]

	i = -1
	k = left_length - 1
	for j in lhs_stencil:

		i = i + 1

		if j == 0:
			lhs[i] = 1.
		else:
			k = k + 1
			lhs[i] = -b[k]
	
	scheme=sc.fdm_scheme(lhs_stencil,lhs,rhs_stencil,rhs)
	scheme.display()
	scheme.spectra_property()

	filename='modified_wavenumber_stanadard_'+form+'_'+str(order_of_accuracy)+'.dat'
	ps.file_write_wavenumber(filename,scheme.wavenumber,scheme.modified_wavenumber_real,scheme.modified_wavenumber_imag)
	ps.plot_wavenumber(scheme.wavenumber,scheme.modified_wavenumber_real,scheme.modified_wavenumber_imag)