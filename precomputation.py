from NFDS import *
from dirichletCharacters import *
from SLTwoZ import *

import NFDS

import os
from os.path import isdir

gensWritePath = "/precomputation"
charWritePath = "/precomputation"
reltWritePath = "/precomputation"

def modPairs(N):
    pairs = []
    for j in range(3, N // 3 + 1):
        for k in range(3, N // j + 1):
            pairs.append((j, k))
    return pairs

def precomputeMatrices(n):
    subPath = os.getcwd() + gensWritePath + str("/mtrxs/")
    filePath = subPath + str(n) + ".mtrxs"
    if not isdir(subPath):
        os.mkdir(subPath)
    if os.path.exists(filePath):
        f = open(filePath, "r")
        mtrxs = [stringMatrix(line[:-1]) for line in f.readlines()]
        f.close()
    else:
        repsDict, gRepsDict = cosetRepsSLTwoZOverGammaOneFast(n), cosetRepsGammaZeroOverGammaOneFast(n)
        mtrxs = []
        for key in gRepsDict:
            entry = gRepsDict[key]
            entry.set_immutable()
            mtrxs.append(entry)
        for key in repsDict:
            for exp in range(1, n + 1):
                tempMtrx = repsDict[key] * T ** exp
                entry = tempMtrx * repsDict[(tempMtrx[1][0] % n, tempMtrx[1][1] % n)] ** (-1)
                entry.set_immutable()
                mtrxs.append(entry)
            for exp in range(0, 3):
                tempMtrx = repsDict[key] * S ** exp
                entry = tempMtrx * repsDict[(tempMtrx[1][0] % n, tempMtrx[1][1] % n)] ** (-1)
                entry.set_immutable()
                mtrxs.append(entry)
        mtrxs = list(set(mtrxs))
        f = open(filePath, "w")
        f.writelines([matrixString(mtrx) + "\n" for mtrx in mtrxs])
        f.close()
    return mtrxs

def filterGammaOne(mtrxs, n):
    outMtrxsIn, outMtrxsOut = [], []
    for mtrx in mtrxs:
        if inGammaOne(mtrx, n):
            outMtrxsIn.append(mtrx)
        else:
            outMtrxsOut.append(mtrx)
    return outMtrxsIn, outMtrxsOut

def precomputeRelations(n):
    subPath = os.getcwd() + gensWritePath + str("/relations/")
    # TODO: CHANGE BACK FILE DIRECTORY
    filePath = subPath + str(n) + ".relts"
    if not isdir(subPath):
        os.mkdir(subPath)
    if os.path.exists(filePath):
        f = open(filePath, "r")
        rref = [[float(element) for element in line.split(",")[:-1]] for line in f.readlines()]
        f.close()
    else:
        U, out = filterGammaOne(precomputeMatrices(n), n)
        m = []
        for URow in tqdm(U):
            lC = NFDS.linearCombo(URow, U, n)
            row = [lC[matrixString(UCol)] for UCol in U]
            m.append(row)
        sys = matrix(ZZ, m) - matrix.identity(len(m))
        rrefSys = sys.rref()
        lines = []
        for row in rrefSys:
            line = ""
            for element in row:
                line += str(float(element)) + ","
            lines.append(line + "\n")
        f = open(filePath, "w")
        f.writelines(lines)
        f.close
        rref = [[float(element) for element in row] for row in rrefSys]
    return rref

def precomputeRelationsLowMem(n):
    subPath = os.getcwd() + gensWritePath + str("/relations/")
    # TODO: CHANGE BACK FILE DIRECTORY
    filePath = subPath + str(n) + "lm.relts"
    if not isdir(subPath):
        os.mkdir(subPath)
    if os.path.exists(filePath):
        f = open(filePath, "r")
        lines = f.readlines()
        f.close()
        numPivots = int(lines[0])
        pivots, free = [list(map(int, line[:-1].split(", "))) for line in lines[1:numPivots + 1]], list(map(int, lines[numPivots + 1].split(", ")))
        relts = [list(map(float, line[:-2].split(",")))for line in lines[numPivots + 2:]]
    else:
        U, out = filterGammaOne(precomputeMatrices(n), n)
        m = []
        for URow in tqdm(U):
            lC = NFDS.linearCombo(URow, U, n)
            row = [lC[matrixString(UCol)] for UCol in U]
            m.append(row)
        sys = matrix(ZZ, m) - matrix.identity(len(m))
        rrefSys = sys.rref()
        pivots, free = pivotsCompSage(rrefSys)
        lines = [str(len(pivots)) + "\n"] + [str(pivot)[1:-1] + "\n" for pivot in pivots] + [str(free)[1:-1] + "\n"]
        relts = []
        for row in rrefSys:
            line = ""
            reltsRow = []
            for i in free:
                line += str(float(row[i])) + ","
                reltsRow.append(float(row[i]))
            lines.append(line + "\n")
            relts.append(reltsRow)
        f = open(filePath, "w")
        f.writelines(lines)
        f.close
    return pivots, free, relts

def precomputeCharacterPairs(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    subPath, filePath = createCharacterPairFile(dChar1, dChar2)
    if filePath != True:
        f = open(filePath, "w")
        f.writelines([matrixString(mtrx) + complexString(newFormDedekindSum(dChar1, dChar2, mtrx)) + "\n" for mtrx in tqdm(precomputeMatrices(n))])
        f.close()

def pivotsComp(m):
    rc, pivots, free = [0, 0], [], []
    while rc[0] < len(m) and rc[1] < len(m[0]):
        if m[rc[0]][rc[1]] == 0:
            free.append(rc[1])
            rc[1] += 1
        else:
            pivots.append([rc[0], rc[1]])
            rc[0] += 1
            rc[1] += 1
    free.extend(range(rc[1], len(m[0])))
    return pivots, free

def pivotsCompSage(m):
    rc, pivots, free = [0, 0], [], []
    while rc[0] < m.nrows() and rc[1] < m.ncols():
        if m[rc[0]][rc[1]] == 0:
            free.append(rc[1])
            rc[1] += 1
        else:
            pivots.append([rc[0], rc[1]])
            rc[0] += 1
            rc[1] += 1
    free.extend(range(rc[1], m.ncols()))
    return pivots, free

def precomputeCharacterPairsFast(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    subPath, filePath = createCharacterPairFile(dChar1, dChar2)
    if filePath != True:
        mtrxsIn, mtrxsOut = filterGammaOne(precomputeMatrices(n), n)
        relts = precomputeRelations(n)
        pivots, free = pivotsComp(relts)
        lines, lookup = [], {}
        for mtrx in tqdm(mtrxsOut):
            lines.append(matrixString(mtrx) + complexString(newFormDedekindSum(dChar1, dChar2, mtrx)) + "\n")
        for i in tqdm(free):
            nfds = newFormDedekindSum(dChar1, dChar2, mtrxsIn[i])
            lines.append(matrixString(mtrxsIn[i]) + complexString(nfds) + "\n")
            lookup[i] = nfds
        for pivot in pivots:
            sum = 0
            for j in free:
                sum += -1 * relts[pivot[0]][j] * lookup[j]
            lines.append(matrixString(mtrxsIn[pivot[1]]) + complexString(sum) + "\n")
        f = open(filePath, "w")
        f.writelines(lines)
        f.close()

def precomputeCharacterPairsFastLowMem(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    subPath, filePath = createCharacterPairFile(dChar1, dChar2)
    if filePath != True:
        mtrxsIn, mtrxsOut = filterGammaOne(precomputeMatrices(n), n)
        pivots, free, relts = precomputeRelationsLowMem(n)
        lines, lookup = [], {}
        for mtrx in tqdm(mtrxsOut):
            lines.append(matrixString(mtrx) + complexString(newFormDedekindSum(dChar1, dChar2, mtrx)) + "\n")
        for i in tqdm(free):
            nfds = newFormDedekindSum(dChar1, dChar2, mtrxsIn[i])
            lines.append(matrixString(mtrxsIn[i]) + complexString(nfds) + "\n")
            lookup[i] = nfds
        for pivot in pivots:
            sum = 0
            for i in range(len(free)):
                sum += -1 * relts[pivot[0]][i] * lookup[free[i]]
            lines.append(matrixString(mtrxsIn[pivot[1]]) + complexString(sum) + "\n")
        f = open(filePath, "w")
        f.writelines(lines)
        f.close()

def createCharacterPairFile(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    str1, str2 = dCharString(dChar1), dCharString(dChar2)
    subPath = os.getcwd() + gensWritePath + str("/characterPairsLM/")
    #TODO: CHANGE BACK FILE DIRECTORY
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str1.split("c")[0] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str1.split(";")[0].split("c")[1] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str1.split(";")[1] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str2.split("c")[0] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str2.split(";")[0].split("c")[1] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str2.split(";")[1] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    filePath = subPath + str(n) + ".chpr"
    if not os.path.exists(filePath):
        f = open(filePath, 'w')
        f.close()
    else:
        filePath = True
    return subPath, filePath

def chprPathFinder(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    str1, str2 = dCharString(dChar1), dCharString(dChar2)
    subPath = os.getcwd() + gensWritePath + \
              str("/characterPairs/") + \
              str1.split("c")[0] + "/" + \
              str1.split(";")[0].split("c")[1] + "/" + \
              str1.split(";")[1] + "/" + \
              str2.split("c")[0] + "/" + \
              str2.split(";")[0].split("c")[1] + "/" + \
              str2.split(";")[1] + "/"
    filePath = subPath + str(n) + ".chpr"
    return filePath

def readAllChpr(path):
    dict = {}
    f = open(path, 'r')
    for line in f.readlines():
        splitted = line[:-1].split(":")
        dict[splitted[0] + ":"] = stringComplex(splitted[1])
    f.close()
    return dict