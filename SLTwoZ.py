from base import *

S = matrix(ZZ, [[0, -1], [1, 0]])
T = matrix(ZZ, [[1, 1], [0, 1]])

S.set_immutable()
T.set_immutable()

def inSLTwoZ(m):
    #checks if m is in sl2z
    return True if m.determinant() == 1 else False

def inGammaZero(m, n):
    #checks if m is in Gamma0(n)
    return inSLTwoZ(m) and (m[1][0] % n == 0)

def inGammaOne(m, n):
    #checks if m is in Gamma1(n)
    return inGammaZero(m, n) and (m[0][0] % n == 1) and (m[1][1] % n == 1)

def TSDecomp(m):
    #returns the TS decomposition of a matrix
    TS = []
    while m[1][0] != 0:
        exp = -1 * floor(m[0][0] / m[1][0])
        m = S * T ** exp * m
        TS.append(-1 * exp)
    pm = int(m[0][0] == 1)
    TS.append((2 * pm - 1) * m[0][1])
    TS.append((len(TS) + pm) % 2)
    return TS

def projectiveLift(c, d, N):
    if c != 0 and d != 0:
        while gcd(c, d) != 1:
            d += N
        b = (-1 * inverse_mod(c, d)) % d
        a = (1 + b * c) // d
        return matrix(ZZ, [[a, b],
                           [c, d]])
    elif c == 0:
        if d == 1:
            return matrix.identity(2)
        else:
            di = inverse_mod(d, N)
            return matrix(ZZ, [[di, (d * di - 1) // N],
                               [N, d]])
    elif d == 0:
        ci = (-1 * inverse_mod(c, N)) % N
        return matrix(ZZ, [[(c * ci + 1) // N, ci],
                           [c, N]])

def PRing(n):
    elements = []
    for c in range(n):
        for d in range(n):
            if gcd(gcd(c, d), n) == 1:
                elements.append((c, d))
    return elements

def TSDecompToRewritingTape(tsd):
    tape = [[matrix.identity(2), T ** tsd[0]]]
    for exp in tsd[1:-1]:
        tape.append([tape[-1][0] * tape[-1][1], S])
        tape.append([tape[-1][0] * tape[-1][1], T ** exp])
    tape.append([tape[-1][0] * tape[-1][1], S ** (2 * tsd[-1])])
    return tape

def cosetRepsSLTwoZOverGammaOneFast(n):
    # returns a set of coset representatives of SL2Z/Gamma1(n), uses the projective lift for all n and creates a dictionary for fast coset finding
    return {(pair[0], pair[1]): projectiveLift(pair[0], pair[1], n) for pair in PRing(n)}

def cosetRepsGammaZeroOverGammaOneFast(n):
    # returns a set of coset representatives of Gamma0(n)/Gamma1(n)
    reps = {1: matrix.identity(2)}
    for i in range(2, n):
        if gcd(i, n) == 1:
            ii = inverse_mod(i, n)
            reps[i] = matrix(ZZ, [[ii, (i * ii - 1) // n],
                                  [ n,                 i]])
    return reps