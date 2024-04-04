import numpy as np

def TaylorSeriesCoef(delta,n,char):

	c = np.zeros(n+1)

	if char == 'f':
		for i in range(0,n+1):
			c[i] = np.power(delta,i)/np.math.factorial(i)
		# for i in c:
		# 	print(Fraction(i)),
		# i=n+1
		# res = np.power(delta,i)/Factorial(i)
	elif char == 'df':
		c[0] = 0
		for i in range(1,n+1):
			c[i] = np.power(delta,i-1)/np.math.factorial(i-1)
	else:
		print(" error in TaylorSeriesCoef",char)
		quit()

	# print(char)
	return c

def truncation_error_analysis(lhs_stencil,lhs_coefficient,rhs_stencil,rhs_coefficient):
	
    max_order_of_accuracy = len(lhs_stencil) + len(rhs_stencil) + 2
	
    first_node = rhs_stencil[0]
    last_node  = rhs_stencil[-1]
	
    coef = np.zeros(max_order_of_accuracy+1)
	
    # right hand side term
    n = -1   
    for i in range(first_node,last_node+1):
        c = TaylorSeriesCoef(i,max_order_of_accuracy ,'f')
        n = n +1
        coef = coef + c*rhs_coefficient[n]
	# 
	# left hand side term
    n = -1   
    for i in range(lhs_stencil[0],lhs_stencil[-1]+1):
        c = TaylorSeriesCoef(i,max_order_of_accuracy,'df')
        n = n +1
        coef = coef - c*lhs_coefficient[n]

    # search for the max non-zero term
    for i in range(0,max_order_of_accuracy+1):
    	if np.abs(coef[i]) > 1.e-12:
    		order_of_accuracy = i - 1
    		truncation_error  = np.abs(coef[i])
    		break

    # print(' Order of accuracy: ',order_of_accuracy)
    # print('  Truncation_error: ',truncation_error,'*O(Delta**',order_of_accuracy,')')
		
    return order_of_accuracy