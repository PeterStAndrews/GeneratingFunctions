#  GeneratingFunctions

## GeneratingFunctions.cpp

`GeneratingFunctions.cpp` solves the generating functions associated to the percolation mapping of an SIR process over a network [1]. The infected fraction of the network versus the transmissibility is saved in a file; once plotted, the phase transition is obtained: 


![plot](https://user-images.githubusercontent.com/29250174/45455013-9c3bea00-b6dd-11e8-9a3a-c69a96fb00e9.png)

The file requires the degree distribution for your network to be stored in an input text file, an example for a Erdos & Renyi random network has been provided. The output will be a text file containing (x,y) scatter data which can be plotted using external software. 

[1] M. E. J. Newman. The spread of epidemic disease on networks (2002)

## Author & license 
Copyright (c) 2017-2018, Peter Mann 
Licensed under the GNU General Public Licence v.2.0 <https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html>.
