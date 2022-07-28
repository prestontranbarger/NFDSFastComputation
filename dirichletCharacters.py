from base import *

def allDCharacters(n):
    return DirichletGroup(n)

def modulus(dChar):
    return dChar.modulus()

def isPrimitive(dChar):
    return dChar.is_primitive()

def isEven(dChar):
    return dChar.is_even()

def dCharString(dChar):
    l = str(dChar).split()
    out = l[3] + "c" + l[6] + ";"
    Tn = totient(modulus(dChar))
    zTn = CyclotomicField(Tn).gen()
    ns = [ss.split(" |--> ")[0].split()[-1] for ss in str(dChar).split(",")]
    for n in ns:
        p = 0
        for j in range(Tn):
            if zTn ** j == dChar(int(n)):
                p = j
        out += n + "-" + str(p) + "-"
    return out[:-1]