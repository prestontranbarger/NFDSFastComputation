from NFDS import *
from dirichletCharacters import *
from SLTwoZ import *

import os
from os.path import isdir

gensWritePath = "/precomputation"
charWritePath = "/precomputation"

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

def precomputeCharacterPairs(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    subPath, filePath = createCharacterPairFile(dChar1, dChar2)
    if filePath != True:
        f = open(filePath, "w")
        f.writelines([matrixString(mtrx) + complexString(newFormDedekindSum(dChar1, dChar2, mtrx)) + "\n" for mtrx in tqdm(precomputeMatrices(n))])
        f.close()

def createCharacterPairFile(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    str1, str2 = dCharString(dChar1), dCharString(dChar2)
    subPath = os.getcwd() + gensWritePath + str("/characterPairs/")
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