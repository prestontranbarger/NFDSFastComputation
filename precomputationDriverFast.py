from precomputation import *

N = 40
for pair in modPairs(N):
    for dChar1 in allDCharacters(pair[0]):
        for dChar2 in allDCharacters(pair[1]):
            if isEven(dChar1) == isEven(dChar2):
                if isPrimitive(dChar1) and isPrimitive(dChar2):
                    print(dChar1, r'&', dChar2)
                    precomputeCharacterPairsFastLowMem(dChar1, dChar2)