from NFDS import *

m = matrix(ZZ, [[17, 32],
                [ 9, 17]])

dChar1 = allDCharacters(3)[1]
dChar2 = allDCharacters(3)[1]
print(isEven(dChar1) == isEven(dChar2))
print(isPrimitive(dChar1))
print(isPrimitive(dChar2))

NFDS = newFormDedekindSumFast(dChar1, dChar2, m, chprPathFinder(dChar1, dChar2))
print(float(NFDS.real()) + float(NFDS.imag()))