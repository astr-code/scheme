import numpy as np


# the function is to calculate the modified wavenumber from a finit-difference scheme
def spectral_analyis(lhs_stencil,lhs_coefficient,rhs_stencil,rhs_coefficient):
	
	kmax=128

	modified_wavenumber = np.zeros([kmax+1], dtype = np.csingle) 
	wavenumber = np.zeros([kmax+1], dtype = np.single) 

	for i in range(0,kmax+1):

		wavenumber[i] = np.pi/kmax*i

		# print(wavenumber[i])
		n = -1

		for j in range(rhs_stencil[0],rhs_stencil[-1]+1):
			n = n + 1
			c=rhs_coefficient[n]*np.exp(float(j)*(0.0+1.0j)*wavenumber[i])
			# print(c,type(c))
			modified_wavenumber[i] = modified_wavenumber[i] + c

		modified_wavenumber[i]=modified_wavenumber[i]*(0.0-1.0j)
		# print(modified_wavenumber[i])

	for i in range(0,kmax+1):

		c=(1.0+0.0j)
		n = -1
		for j in range(lhs_stencil[0],lhs_stencil[-1]+1):
			n = n + 1

			if j==0:
				continue
			c = c + lhs_coefficient[n]*np.exp(float(j)*(0.0+1.0j)*wavenumber[i])

		modified_wavenumber[i]=modified_wavenumber[i]/c
		# print(wavenumber[i],modified_wavenumber[i])

	return wavenumber,modified_wavenumber.real, modified_wavenumber.imag


# the function is to calculate the modified wavenumber from a finit-difference scheme
def spectral_analyis_filter(lhs_stencil,lhs_coefficient,rhs_stencil,rhs_coefficient,rhs_alpha_term,alpha):
	
	kmax=128

	modified_wavenumber = np.zeros([kmax+1], dtype = np.csingle) 
	wavenumber = np.zeros([kmax+1], dtype = np.single) 
	rhs_coefficient_full =  np.zeros([rhs_coefficient.size], dtype = float) 

	n=-1
	for j in range(rhs_stencil[0],rhs_stencil[-1]+1):
		n=n+1
		rhs_coefficient_full[n]=rhs_coefficient[n]+rhs_alpha_term[n]*alpha

	for i in range(0,kmax+1):

		wavenumber[i] = np.pi/kmax*i

		# print(wavenumber[i])
		n = -1

		for j in range(rhs_stencil[0],rhs_stencil[-1]+1):
			n = n + 1
			c=rhs_coefficient_full[n]*np.exp(float(j)*(0.0+1.0j)*wavenumber[i])
			# print(c,type(c))
			modified_wavenumber[i] = modified_wavenumber[i] + c
		
		var_alpha=0.0
		n = -1
		for j in range(lhs_stencil[0],lhs_stencil[-1]+1):
			n = n + 1
			if lhs_stencil[n] != 0:
				var_alpha = var_alpha + alpha*np.exp((0.0+1.0j)*j*wavenumber[i])

		modified_wavenumber[i] = modified_wavenumber[i]/(1.0+var_alpha)

	return wavenumber,modified_wavenumber.real, modified_wavenumber.imag