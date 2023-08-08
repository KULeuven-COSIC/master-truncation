"""

Data saved to Player 0 input

"""

import numpy as np

#Writes the data to each player's input file in the order below
def main():
	#Clears all existing data in the Players' input files
	open("../../Player-Data/Input-P0-0", 'w').close() 
	open("../../Player-Data/Input-P1-0", 'w').close()
	open("../../Player-Data/Input-P2-0", 'w').close()
	WriteFile('X_train.txt', 0)
	WriteFile('Y_train.txt', 1)
	WriteFile('X_val.txt', 2)
	WriteFile('Y_val.txt', 2)
	
#Loads data from a file where the first 2 rows always correspond to the dimension of the (multi)array
def LoadFile(filename):
	dim = []
	with open(filename) as f:
		dim.append(f.readline())
		dim.append(f.readline())
	data = np.loadtxt(filename, skiprows = 2)
	return dim, data 

#Writes data into the chosen player's input file so that before each array the 2 values in separate lines stand for the dimension of the array
def WriteFile(filename, player):
	data = LoadFile(filename)
	with open("../../Player-Data/Input-P"+str(player)+"-0", 'a') as f:
		f.write(data[0][0]+ data[0][1])
		np.savetxt(f, data[1])

from mnist import get_data_files
get_data_files()
main()

