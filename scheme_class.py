import numpy as np
from fractions import Fraction
import taylor as ty
import spectra as sp

# this is to convert a fraction number string to a latex string
def latex_fraction_number(fraction_number):
    str_break = fraction_number.split('/')

    if len(str_break) == 1:
        frac_numb = str_break[0]
        if str_break[0] == '1':
            frac_numb=''
        if str_break[0] == '-1':
            frac_numb='-'
        if str_break[0] == '0':
            frac_numb=''
    elif len(str_break) == 2:
        if int(str_break[0])<0:
            pvalue =str(-int(str_break[0]))
            frac_numb = "-\\frac{"+pvalue+'}'+'{'+str_break[1]+'}'
        else:
            frac_numb = "\\frac{"+str_break[0]+'}'+'{'+str_break[1]+'}'
    return frac_numb

class filter_scheme(object):
    
    def __init__(self,lhs_stencil,lhs_coefficient,rhs_stencil,rhs_coefficient,rhs_alpha_term):

        if (lhs_stencil != 0).any():
            self.form = 'compact'
        else:
            self.form = 'explicit'

        self.lhs_stencil       = lhs_stencil
        self.rhs_stencil       = rhs_stencil

        self.first_node        = self.rhs_stencil[0]
        self.last_node         = self.rhs_stencil[-1]
        self.length_of_stencil = self.rhs_stencil.size
        self.lhs_coefficient   = lhs_coefficient
        self.rhs_coefficient   = rhs_coefficient
        self.rhs_alpha_term    = rhs_alpha_term

    def display(self):
        # print(' length of LHS stencil:',self.length_of_stencil,type(self.length_of_stencil))
        print('')
        print(' ---------------------------------- information of the scheme ----------------------------------')
        print('           Scheme type:',self.form)
        print('           LHS stencil:',self.lhs_stencil)
        print('       LHS coefficient:',self.lhs_coefficient)
        print('           RHS stencil:',self.rhs_stencil)
        print('       RHS coefficient:',self.rhs_coefficient,self.rhs_alpha_term)

        j = -1
        for i in self.lhs_stencil:
            j = j+1
            if j==0:
                char1=''
            else:
                char1='+'
            
            if i > 0:
                char2='α*ff(i+'+str(i)+')'
            elif i == 0:
                char2='ff(i)'
            else:
                char2='α*ff(i'+str(i)+')'
            
            print(char1,char2,end='')

        print(' = ',end='')

        j = -1
        oper=''
        for i in range(self.first_node,self.last_node+1):
            j = j+1
            if self.rhs_coefficient[j] > 0.0:
                char1='+'
            else:
                char1=''
            
            if i > 0:
                char2='*f(i+'+str(i)+')'
            elif i == 0:
                char2='*f(i)'
            else:
                char2='*f(i'+str(i)+')'

            char3=Fraction(self.rhs_coefficient[j]).limit_denominator()
            char4=Fraction(self.rhs_alpha_term[j]).limit_denominator()

            if self.rhs_alpha_term[j]>0.0:
                oper='+'
            else:
                oper=''

            char5='('+str(char3)+str(oper)+str(char4)+'*α)'
            print(char1,char5,char2,end='')

        print('')
        print(' ----------------------------------------------------------------------------------------------')


        j = -1
        for i in range(self.first_node,self.last_node+1):
            j = j+1
            if self.rhs_alpha_term[j]>0.0:
                oper='+'
            else:
                oper=''
            char3=Fraction(self.rhs_coefficient[j]).limit_denominator()
            numer1=char3.numerator
            char4=Fraction(self.rhs_alpha_term[j]).limit_denominator()
            numer2=char4.numerator
            denom1 = char3.denominator
            denom2 = char4.denominator
            if denom1 > denom2:
                numer2=int(numer2/denom2*denom1)
                denom=denom1
            elif denom2 > denom1:
                numer1=int(numer1/denom1*denom2)
                denom=denom2
            else:
                denom=denom1

            char5='('+str(numer1)+'.d0'+str(oper)+str(numer2)+'.d0'+'*alfa)/'+str(denom)+'.d0'
            print(char5)

    def spectra_property(self):
        alpha=0.499
        self.wavenumber,self.modified_wavenumber_real,self.modified_wavenumber_imag = sp.spectral_analyis_filter(self.lhs_stencil,self.lhs_coefficient,  + \
                                                                                      self.rhs_stencil,self.rhs_coefficient,self.rhs_alpha_term,alpha)
        # c = np.savetxt('spectra.txt', wavenumer, self.modified_wavenumber_real,delimiter =', ') 
        # print(wavenumer,self.modified_wavenumber_real,self.modified_wavenumber_imag)

    def get_formula(self):

        j = -1
        for i in self.lhs_stencil:
            j = j+1
            if self.lhs_coefficient[j] >= 0.0:
                char1='+'
                if j==0:
                    char1=''
                    formula=''
            else:
                char1=''
            
            if i > 0:
                char2='{α\widetilde{f}_{i+'+str(i)+'}}'
            elif i == 0:
                char2='{\widetilde{f}_i}'
            else:
                char2='{α\widetilde{f}_{i'+str(i)+'}}'

            frac_numb=latex_fraction_number(str(Fraction(self.lhs_coefficient[j]).limit_denominator()))

            formula = formula + char1 + frac_numb+char2

        formula = formula+'='
        j = -1
        for i in range(self.first_node,self.last_node+1):
            j = j+1
            if self.rhs_coefficient[j] > 0.0:
                char1='+'
            else:
                char1=''
            
            if i > 0:
                char2='f_{i+'+str(i)+'}'
            elif i == 0:
                char2='f_{i}'
            else:
                char2='f_{i'+str(i)+'}'

            if abs(self.rhs_coefficient[j]) < 1.e-10:
                char1=''
                char2=''
            
            frac_numb=latex_fraction_number(str(Fraction(self.rhs_coefficient[j]).limit_denominator()))
            frac_alpha=latex_fraction_number(str(Fraction(self.rhs_alpha_term[j]).limit_denominator()))
            if(self.rhs_alpha_term[j]>0):
                frac_alpha='+'+frac_alpha+'α'

            formula = formula + char1+'('+frac_numb+frac_alpha+')'+char2
        
        print(formula)
        return formula

class fdm_scheme(object):
    
    def __init__(self,lhs_stencil,lhs_coefficient,rhs_stencil,rhs_coefficient):

        if (lhs_stencil != 0).any():
            self.form = 'compact'
        else:
            self.form = 'explicit'

        self.lhs_stencil       = lhs_stencil
        self.rhs_stencil       = rhs_stencil

        self.first_node        = self.rhs_stencil[0]
        self.last_node         = self.rhs_stencil[-1]
        self.length_of_stencil = self.rhs_stencil.size
        self.rhs_coefficient   = rhs_coefficient
        self.rhs_coefficient   = rhs_coefficient
        self.lhs_coefficient   = lhs_coefficient

        self.order_of_accuracy  = ty.truncation_error_analysis(self.lhs_stencil,self.lhs_coefficient,self.rhs_stencil,self.rhs_coefficient)

    def display(self):
        # print(' length of LHS stencil:',self.length_of_stencil,type(self.length_of_stencil))
        print('')
        print(' ---------------------------------- information of the scheme ----------------------------------')
        print('           LHS stencil:',self.lhs_stencil)
        print('       LHS coefficient:',self.lhs_coefficient)
        print('           RHS stencil:',self.rhs_stencil)
        print('       RHS coefficient:',self.rhs_coefficient)
        print('     Order of accuracy:',self.order_of_accuracy)

        j = -1
        for i in self.lhs_stencil:
            j = j+1
            if self.lhs_coefficient[j] > 0.0:
                char1='+'
                if j==0:
                    char1=''
            else:
                char1=''
            
            if i > 0:
                char2='*df(i+'+str(i)+')'
            elif i == 0:
                char2='*df(i)'
            else:
                char2='*df(i'+str(i)+')'
            
            print(char1,"{}{}".format(Fraction(self.lhs_coefficient[j]).limit_denominator(),char2),end='')

        print(' = ',end='')

        j = -1
        for i in range(self.first_node,self.last_node+1):
            j = j+1
            if self.rhs_coefficient[j] > 0.0:
                char1='+'
            else:
                char1=''
            
            if i > 0:
                char2='*f(i+'+str(i)+')'
            elif i == 0:
                char2='*f(i)'
            else:
                char2='*f(i'+str(i)+')'
            
            print(char1,"{}{}".format(Fraction(self.rhs_coefficient[j]).limit_denominator(),char2),end='')
        print('')
        print(' ----------------------------------------------------------------------------------------------')

        j = -1
        for i in self.lhs_stencil:
            j = j+1
            if self.lhs_coefficient[j] > 0.0:
                char1='+'
                if j==0:
                    char1=''
            else:
                char1=''
            
            if i > 0:
                char2='*f(i+'+str(i*2+1)+'/2)'
            elif i == 0:
                char2='*f(i+1/2)'
            else:
                char2='*f(i'+str(i*2+1)+'/2)'
            
            print(char1,"{}{}".format(Fraction(self.lhs_coefficient[j]).limit_denominator(),char2),end='')

        print(' = ',end='')

        # calculate flux form at i+1/2

        flux_coefficient = np.zeros([len(self.rhs_coefficient)-1], dtype = float)

        for j in range(0,len(flux_coefficient)):

            if j==0 :
                flux_coefficient[j] = -self.rhs_coefficient[j]
            else:
                flux_coefficient[j] = flux_coefficient[j-1] - self.rhs_coefficient[j]

        j = -1
        for i in range(self.first_node+1,self.last_node+1):
            j = j+1
            if flux_coefficient[j] > 0.0:
                char1='+'
            else:
                char1=''
            
            if i > 0:
                char2='*f(i+'+str(i)+')'
            elif i == 0:
                char2='*f(i)'
            else:
                char2='*f(i'+str(i)+')'
            
            print(char1,"{}{}".format(Fraction(flux_coefficient[j]).limit_denominator(),char2),end='')

        print('')
        print(' ----------------------------------------------------------------------------------------------')

        # print(len(self.rhs_coefficient),self.rhs_coefficient)

    def get_formula(self):

        j = -1
        formula=''
        for i in self.lhs_stencil:
            j = j+1
            if self.lhs_coefficient[j] > 0.0:
                char1='+'
                if j==0:
                    char1=''
                    formula=''
            else:
                char1=''
            
            if i > 0:
                char2='{f_{i+'+str(i)+'}}\''
            elif i == 0:
                char2='{f_i}\''
            else:
                char2='{f_{i'+str(i)+'}}\''

            frac_numb=latex_fraction_number(str(Fraction(self.lhs_coefficient[j]).limit_denominator()))

            formula = formula + char1 + frac_numb+char2

        formula = formula+'='
        j = -1
        for i in range(self.first_node,self.last_node+1):
            j = j+1
            if self.rhs_coefficient[j] > 0.0:
                char1='+'
            else:
                char1=''
            
            if i > 0:
                char2='f_{i+'+str(i)+'}'
            elif i == 0:
                char2='f_{i}'
            else:
                char2='f_{i'+str(i)+'}'

            if abs(self.rhs_coefficient[j]) < 1.e-10:
                char1=''
                char2=''
            
            frac_numb=latex_fraction_number(str(Fraction(self.rhs_coefficient[j]).limit_denominator()))

            formula = formula + char1+frac_numb+char2
        
        formula = formula + '+O(\\Delta^{'+str(self.order_of_accuracy)+'})'
        # print(formula)
        return formula
        
    def spectra_property(self):
        self.wavenumber,self.modified_wavenumber_real,self.modified_wavenumber_imag = sp.spectral_analyis(self.lhs_stencil,self.lhs_coefficient,self.rhs_stencil,self.rhs_coefficient)
        # c = np.savetxt('spectra.txt', wavenumer, self.modified_wavenumber_real,delimiter =', ') 
        # print(wavenumer,self.modified_wavenumber_real,self.modified_wavenumber_imag)
