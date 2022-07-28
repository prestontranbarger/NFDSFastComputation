from sage.all import *
from tqdm import tqdm

def matrixString(m):
  return str(m[0][0]) + "," + str(m[0][1]) + "," + str(m[1][0]) + "," + str(m[1][1]) + ":"

def stringMatrix(s):
  l = s[:-1].split(",")
  return matrix(ZZ, [[int(l[0]), int(l[1])],
                     [int(l[2]), int(l[3])]])

def complexString(z):
  z += 0 * I
  return "(" + str(float(z.real())) + "," + str(float(z.imag())) + ")"

def stringComplex(s):
  l = s[1:-1].split(",")
  return float(l[0]) + float(l[1]) * I

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