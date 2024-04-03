import numpy as np
import matplotlib.pyplot as plt

def plot_wavenumber(k,mkr,mki):
	plt.plot(k,k, color="grey", linewidth=0.5, linestyle="-", label="exact")
	plt.plot(k,mkr, color="red", linewidth=1.5, linestyle="-", label="real")
	plt.plot(k,mki, color="blue", linewidth=1.5, linestyle="-", label="image")
	plt.legend()
	plt.xlabel('wavenumber')
	plt.xlim(0, np.pi)
	plt.ylabel('modified wavenumber')
	plt.ylim(min(mki)*1.01, np.pi)
	plt.axhline(linewidth=0.5, color='black')
	plt.show()