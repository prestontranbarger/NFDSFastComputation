from NFDS import *

def lagrangeInterpolate(pts):
    polyRing = PolynomialRing(QQ, names=('x',)); (x,) = polyRing._first_ngens(1)
    polyOut = 0 * x
    for i in range(len(pts)):
        term = 1 + 0 * x
        for j in range(len(pts)):
            if i != j:
                term *= (x - pts[j][0]) / (pts[i][0] - pts[j][0])
        polyOut += pts[i][1] * term
    return polyOut

def hFunctionInterpolate(gamma, dChar1, dChar2, k):
    q1, q2 = modulus(dChar1), modulus(dChar2)
    n = q1 * q2
    pts = []
    for c in tqdm(range(n, k * n, n)):
        mtrx = matrix(ZZ, [[1, 0], [c, 1]])
        hFunction = higherWeightDedekindSum(k, dChar1, dChar2, mtrx) - (gamma[1][0] * 1 / c + gamma[1][1]) ** (k - 2) * higherWeightDedekindSum(k, dChar1, dChar2, gamma * mtrx)
        pts.append((QQ(1 / c), hFunction))
    return lagrangeInterpolate(pts)

# gamma = matrix(ZZ, [[1351, 2755], [1300, 2651]])
# dChar1 = allDCharacters(5)[2]
# dChar2 = allDCharacters(5)[2]
# print(hFunctionInterpolate(gamma, dChar1, dChar2, 4))

dChar1 = allDCharacters(8)[3]
dChar2 = allDCharacters(4)[1]
print(dCharString(dChar1), dCharString(dChar2))

for k in range(4, 9, 2):
    lines = []
    for gen in tqdm(Gamma1(modulus(dChar1) * modulus(dChar2)).gens()):
        hF = hFunctionInterpolate(gen, dChar1, dChar2, k)
        print("\n" + str(gen), "\n", hF)
        print()
        lines.append(str(gen[0][0]) + "_" + str(gen[0][1]) + "_" + str(gen[1][0]) + "_" + str(gen[1][1]) + ":" + str(hF) + "\n")

    fileOut = open('hF-k' + str(k) + ',' + dCharString(dChar1) + ';' + dCharString(dChar2) + '.txt', 'w')
    fileOut.writelines(lines)
    fileOut.close()

print(dCharString(dChar1), dCharString(dChar2))