"""
   EE2703 Applied Programming Lab - 2021
   Assignment 1
"""
from sys import * #imports everything from the module sys

if len(argv)!=2 :
    print('\nUsage: %s <inputfile>' % argv[0])
    exit()
"""To check whether the user have entered only the required inputs ,if not then to show the expected usage"""
circuit='.circuit' #initialising so that we only need to change the constant if we decide to change the command line 
end='.end'
try:
    with(open(argv[1])) as f:
        lines=f.readlines()
        start=-1;stop=-2
        for line in lines:
            if circuit==line[:len(circuit)].lower():
                start=lines.index(line)
            elif end==line[:len(end)].lower():
                stop=lines.index(line)
                break
        if start >= stop:
            print('INVALID CIRCUIT DEFINITION')
            exit(0)
        revstr = []
        for line in lines[start+1:stop]:
            str = line.split('#')[0].split()
            revstr.append(' '.join(reversed(str)))

        for line in reversed(revstr):
            print(line)
except IOError:
    print('INVALID FILE')
    exit()
"""
try-catch is used because if we input a wrong file then open function will give an IOError.
"""
