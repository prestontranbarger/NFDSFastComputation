from NFDS import *

import time

dChar1 = allDCharacters(5)[3]
dChar2 = allDCharacters(7)[5]

print(dChar1)
print(dChar2)
print(isEven(dChar1) == isEven(dChar2), isPrimitive(dChar1), isPrimitive(dChar2))
print("preliminary checks complete!")

gamma = matrix(ZZ, [[46741638, 43234369],
                    [43234205, 39990117]])
print("matrix construction complete!")

chpr = readAllChpr(chprPathFinder(dChar1, dChar2))
print("chpr read complete!")

print("beginning fast computation...")
bT = time.time()
print(newFormDedekindSumFast(dChar1, dChar2, gamma, chpr))
bT = time.time() - bT
print("finished fast computation!")
print("final time: " + str(bT) + "s")

print("beginning slow computation...")
bT = time.time()
print(newFormDedekindSum(dChar1, dChar2, gamma, True))
bT = time.time() - bT
print("finished slow computation!")
print("final time: " + str(bT) + "s")