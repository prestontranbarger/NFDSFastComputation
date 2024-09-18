from NFDS import *
lines = []

k = 4
dChar1 = allDCharacters(5)[2]
dChar2 = allDCharacters(5)[2]
n = modulus(dChar1) * modulus(dChar2)

print(isPrimitive(dChar1), isQuadratic(dChar1))
print(isPrimitive(dChar2), isQuadratic(dChar2))
print((dChar1(-1) * dChar2(-1)) == ((-1) ** k))

gammaFixed = matrix(ZZ, [[3926, 155], [1925, 76]])
print(str(gammaFixed))

lines.append(str(k) + "," + dCharString (dChar1) + ";" + dCharString(dChar2) + ","
             + str(gammaFixed[0][0]) + "_" + str(gammaFixed[0][1]) + "_"
             + str(gammaFixed[1][0]) + "_" + str(gammaFixed[1][1]) + "\n")

depth = 15
for c in tqdm(range(n, depth * n, n)):
    for a in tqdm(range(1, c, n)):
        d, u, v = xgcd(a, c)
        if d == 1:
            mtrx = matrix(ZZ, [[a, -1 * v], [c, u]])
            hFunction = higherWeightDedekindSum(k, dChar1, dChar2, mtrx) - (gammaFixed[1][0] * a / c + gammaFixed[1][1]) ** (k - 2) * higherWeightDedekindSum(k, dChar1, dChar2, gammaFixed * mtrx)
            lines.append(str(a) + "/" + str(c) + ": " + str(hFunction) + "\n")
print(dChar1, dChar2)

writeFile = open("hFunctionPrecomputation.txt", 'w')
writeFile.writelines(lines)
writeFile.close()