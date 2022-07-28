from base import *

def allDCharacters(n):
    return DirichletGroup(n)

def modulus(dChar):
    return dChar.modulus()

def isPrimitive(dChar):
    return dChar.is_primitive()

def isEven(dChar):
    return dChar.is_even()