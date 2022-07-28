from precomputation import *

def sawtooth(x):
    #sawtooth function
    if x == floor(x):
        return 0
    return x - floor(x) - 1 / 2

def newFormDedekindSum(dChar1, dChar2, gamma):
    #computes the new form dedekind sum of gamma given two primative characters with similar parity
    #this is in accordance with SVY's definition of a finite double sum formula
    sum = 0
    q1, q2 = modulus(dChar1), modulus(dChar2)
    a, c = gamma[0][0] if gamma[1][0] > 0 else -1 * gamma[0][0],\
           gamma[1][0] if gamma[1][0] > 0 else -1 * gamma[1][0]
    for j in tqdm(range(c)):
        for n in range(q1):
            sum += dChar2(j).conjugate() * dChar1(n).conjugate() * (sawtooth(j / c) + 0 * I) * (sawtooth(n / q1 + a * j / c) + 0 * I)
    return sum

def newFormDedekindSumFast(dChar1, dChar2, gamma, chprPath):
    n = modulus(dChar1) * modulus(dChar2)
    repsDict, gRepsDict = cosetRepsSLTwoZOverGammaOneFast(n), cosetRepsGammaZeroOverGammaOneFast(n)
    g = gRepsDict[gamma[1][1] % n]
    m = gamma * g ** (-1)
    precomputedDedekindSums = readAllChpr(chprPath)
    rwt, sum = TSDecompToRewritingTape(TSDecomp(m)), precomputedDedekindSums[matrixString(g)]
    for i in tqdm(range(0, len(rwt), 2)):
        cosetRep1, cosetRep2, a = repsDict[(rwt[i][0][1][0] % n, rwt[i][0][1][1] % n)],\
                                  repsDict[(rwt[i + 1][0][1][0] % n, rwt[i + 1][0][1][1] % n)],\
                                  int(rwt[i][1][0][1])
        q, r = a // n, a % n
        cosetRep1n, cosetRep1r, cosetRep2t = cosetRep1 * T ** n, cosetRep1 * T ** r, cosetRep2 * rwt[i + 1][1]
        sum += precomputedDedekindSums[matrixString(cosetRep1n * repsDict[(cosetRep1n[1][0] % n, cosetRep1n[1][1] % n)] ** (-1))] * q +\
               precomputedDedekindSums[matrixString(cosetRep1r * repsDict[(cosetRep1r[1][0] % n, cosetRep1r[1][1] % n)] ** (-1))] +\
               precomputedDedekindSums[matrixString(cosetRep2t * repsDict[(cosetRep2t[1][0] % n, cosetRep2t[1][1] % n)] ** (-1))]
    return sum