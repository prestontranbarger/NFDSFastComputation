import matplotlib as mpl
mpl.use('Agg')
print(mpl.get_cachedir())

from NFDS import *
import matplotlib.pyplot as plt

#plt.rcParams['text.usetex'] = True

lines = []
readFile = open("temp.txt", 'r')
linesIn = readFile.readlines()

maxHeight = 1000000
for line in linesIn[1:]:
    (ac, nfds) = line.split(": ")
    nfds = float(QQ(nfds))
    if abs(nfds) < maxHeight:
        lines.append([float(QQ(ac)), nfds])

lines.sort()

xVals = [line[0] for line in lines]
yVals = [line[1] for line in lines]

plt.scatter(xVals, yVals, s = 1)
plt.xlabel(r'$a/c$')
plt.ylabel(r'$h_{\gamma,\chi_1,\chi_2,k}(a/c)$')

plt.savefig("h" + linesIn[0][:-1] + ".png", dpi = 500, bbox_inches='tight')