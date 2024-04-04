import fdm as fdm

def scheme_analysis():
	order_of_accuracy = input(" >> input the order of the scheme:  ")
	order_of_accuracy = int(order_of_accuracy)
	form = input(" >> the scheme is compact or explicit:  ")
	fdm.scheme_analysis(order_of_accuracy,form)

def main():
	scheme_analysis()

if __name__ == '__main__':
	main()
