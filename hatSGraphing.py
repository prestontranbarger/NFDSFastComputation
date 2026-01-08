import math

import matplotlib as mpl
mpl.use('Agg')
print(mpl.get_cachedir())

from NFDS import *
import matplotlib.pyplot as plt

#plt.rcParams['text.usetex'] = True

k = 4
dChar1 = allDCharacters(5)[2]
dChar2 = allDCharacters(5)[2]
print(isQuadratic(dChar1), isPrimitive(dChar1), isQuadratic(dChar2), isPrimitive(dChar2), dChar1(-1) * dChar2(-1) == (-1) ** k)
q_1, q_2 = modulus(dChar1), modulus(dChar2)
N = q_1 * q_2

j = 50
maxHeight = 2
lines = []
pairs = [(am, cm) for cm in range(1, j) for am in range(0, cm)]
for am, cm in tqdm(pairs):
        a, c = 1 + am * N, cm * N
        ac = a / c
        if 1 == math.gcd(a, c) :
            d = inverse_mod(a, c)
            b = (a * d - 1) // c
            gamma = matrix(ZZ, [[a, b], [c, d]])
            nfds = higherWeightDedekindSum(k, dChar1, dChar2, gamma)
            if abs(nfds) < maxHeight:
                lines.append([float(QQ(ac)), nfds])
lines.sort()

xVals = [line[0] for line in lines]
yVals = [line[1] for line in lines]

plt.scatter(xVals, yVals, s = 1)
plt.xlabel(r'$a/c$')
plt.ylabel(r'$\widehat{S}_{\chi_1,\chi_2,k}(a/c)$')

plt.savefig("hatS-k" + str(k) + "," + dCharString(dChar1) + ";" + dCharString(dChar2) + ".png", dpi = 500, bbox_inches='tight')