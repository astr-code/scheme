import numpy as np
from fractions import Fraction
import taylor as ty
import spectra as sp

class fdm_compact_tridiagonal(object):
    
    lhs_stencil   = np.array([-1, 0, 1])

    def __init__(self,first_node,last_node,coefficient_lhs,coefficient_rhs):
        self.first_node        = first_node
        self.last_node         = last_node
        self.length_of_stencil = last_node-first_node+1
        self.coefficient_lhs   = coefficient_lhs
        self.coefficient_rhs   = coefficient_rhs

        self.rhs_stencil = np.zeros(self.length_of_stencil,dtype=np.int32)
        for i in range(first_node,last_node+1):
            self.rhs_stencil[i-first_node] = i
        
        self.order_of_accuracy  = ty.truncation_error_analysis(self.lhs_stencil,self.coefficient_lhs,self.rhs_stencil,self.coefficient_rhs)

    def display(self):
        # print(' length of LHS stencil:',self.length_of_stencil,type(self.length_of_stencil))
        print('')
        print(' ---------------------------------- information of the scheme ----------------------------------')
        print('           LHS stencil:',self.lhs_stencil)
        print('       LHS coefficient:',self.coefficient_lhs)
        print('           RHS stencil:',self.rhs_stencil)
        print('       RHS coefficient:',self.coefficient_rhs)
        print('     Order of accuracy:',self.order_of_accuracy)
        print('      format of scheme:',"{}{}{}{}".format(Fraction(self.coefficient_lhs[0]).limit_denominator(),'*df(i-1) + df(i) + ',Fraction(self.coefficient_lhs[2]).limit_denominator(),'*df(i+1) = '))
        j = -1
        for i in range(self.first_node,self.last_node+1):
            j = j+1
            if self.coefficient_rhs[j] > 0.0:
                char1='+'
            else:
                char1=''
            
            if i > 0:
                char2='*f(i+'+str(i)+')'
            elif i == 0:
                char2='*f(i)'
            else:
                char2='*f(i'+str(i)+')'
            
            print(char1,"{}{}".format(Fraction(self.coefficient_rhs[j]).limit_denominator(),char2),end='')
        print('')
        print(' ----------------------------------------------------------------------------------------------')

    def spectra_property(self):
        self.wavenumber,self.modified_wavenumber_real,self.modified_wavenumber_imag = sp.spectral_analyis(self.lhs_stencil,self.coefficient_lhs,self.rhs_stencil,self.coefficient_rhs)
        # c = np.savetxt('spectra.txt', wavenumer, self.modified_wavenumber_real,delimiter =', ') 
        # print(wavenumer,self.modified_wavenumber_real,self.modified_wavenumber_imag)