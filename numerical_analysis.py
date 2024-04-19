import fdm as fdm
import numpy as np

def generate_consecutive_array(first, last):
    return list(range(first, last + 1))

def scheme_analysis():
	print( "   You are analysing a numerical scheme through giving a set of stencils **")
	print( "     A scheme should have a L.H.S stencil and a R.H.S. stencil **")
	print( "     First you need to input the range of the stencils, **")
	print( "       for example, for a stencil ranging from i-1 to i+2, you need to first input -1, and then input 1")
	
	print(" >> The first node of the L.H.S. stencil:  ", end='')
	number_1    = int(input())
	print(" >> The last node of the L.H.S. stencil:  ", end='')
	number_2    = int(input())
	rhs_stencil = generate_consecutive_array(number_1,number_2)
	rhs_stencil = np.array(rhs_stencil)


	print(" >> The first node of the R.H.S. stencil:  ", end='')
	number_1    = int(input())
	print(" >> The last node of the R.H.S. stencil:  ", end='')
	number_2    = int(input())
	lhs_stencil = generate_consecutive_array(number_1,number_2)
	lhs_stencil = np.array(lhs_stencil)

	fdm.scheme_analysis_stencil(rhs_stencil,lhs_stencil)

def scheme_derive():
	first_node_left,last_node_left = input(" >> input the first and last nodes on the LHS of the scheme:  ").split()
	first_node_left = int(first_node_left)
	last_node_left  = int(last_node_left)
	first_node_right,last_node_right = input(" >> input the first and last nodes on the RHS of the scheme:  ").split()
	first_node_right = int(first_node_right)
	last_node_right  = int(last_node_right)
	order_of_accuracy = input(" >> input the order of the scheme:  ")
	order_of_accuracy = int(order_of_accuracy)

	fdm.scheme_derive_parameter(first_node_left,last_node_left,first_node_right,last_node_right,order_of_accuracy)
# this is to show the scheme definition:

def main():
	scheme_analysis()

if __name__ == '__main__':
	main()
