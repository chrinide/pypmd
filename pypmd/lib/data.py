#!/usr/bin/env python

import numpy

EPSOCC = 1e-6
CUTTZ = 1e-14
nlm = numpy.zeros((56,3), dtype=numpy.int32)

# p's
nlm[1,0] =  1 #px
nlm[2,1] =  1 #py
nlm[3,2] =  1 #pz

# d's
nlm[4,0] =  2  #xx
nlm[5,1] =  2  #yy
nlm[6,2] =  2  #zz
nlm[7,0] =  1  #xy
nlm[7,1] =  1
nlm[8,0] =  1  #xz
nlm[8,2] =  1
nlm[9,1] =  1  #yz
nlm[9,2] =  1

# f's
nlm[10,0] = 3 #xxx
nlm[11,1] = 3 #yyy
nlm[12,2] = 3 #zzz    
nlm[13,0] = 2 #xxy
nlm[13,1] = 1
nlm[14,0] = 2 #xxz
nlm[14,2] = 1
nlm[15,1] = 2 #yyz
nlm[15,2] = 1 
nlm[16,0] = 1 #xyy
nlm[16,1] = 2 
nlm[17,0] = 1 #xzz
nlm[17,2] = 2 
nlm[18,1] = 1 #yzz
nlm[18,2] = 2
nlm[19,0] = 1 #xyz
nlm[19,1] = 1 
nlm[19,2] = 1 

# g's
nlm[20,0] = 4 #xxxx
nlm[21,1] = 4 #yyyy
nlm[22,2] = 4 #zzzz
nlm[23,0] = 3 #xxxy
nlm[23,1] = 1
nlm[24,0] = 3 #xxxz
nlm[24,2] = 1 
nlm[25,0] = 1 #xyyy
nlm[25,1] = 3 
nlm[26,1] = 3 #yyyz
nlm[26,2] = 1 
nlm[27,0] = 1 #xzzz 
nlm[27,2] = 3 
nlm[28,1] = 1 #yzzz
nlm[28,2] = 3
nlm[29,0] = 2 #xxyy
nlm[29,1] = 2 
nlm[30,0] = 2 #xxzz
nlm[30,2] = 2
nlm[31,1] = 2 #yyzz 
nlm[31,2] = 2 
nlm[32,0] = 2 #xxyz 
nlm[32,1] = 1
nlm[32,2] = 1
nlm[33,0] = 1 #xyyz
nlm[33,1] = 2 
nlm[33,2] = 1
nlm[34,0] = 1 #xyzz
nlm[34,1] = 1
nlm[34,2] = 2

# h's
nlm[35,0] = 0
nlm[35,1] = 0
nlm[35,2] = 5
nlm[36,0] = 0
nlm[36,1] = 1
nlm[36,2] = 4
nlm[37,0] = 0
nlm[37,1] = 2
nlm[37,2] = 3
nlm[38,0] = 0
nlm[38,1] = 3
nlm[38,2] = 2
nlm[39,0] = 0
nlm[39,1] = 4
nlm[39,2] = 1
nlm[40,0] = 0
nlm[40,1] = 5
nlm[40,2] = 0
nlm[41,0] = 1
nlm[41,1] = 0
nlm[41,2] = 4
nlm[42,0] = 1 
nlm[42,1] = 1
nlm[42,2] = 3
nlm[43,0] = 1
nlm[43,1] = 2
nlm[43,2] = 2
nlm[44,0] = 1
nlm[44,1] = 3
nlm[44,2] = 1
nlm[45,0] = 1
nlm[45,1] = 4
nlm[45,2] = 0
nlm[46,0] = 2
nlm[46,1] = 0
nlm[46,2] = 3
nlm[47,0] = 2
nlm[47,1] = 1
nlm[47,2] = 2
nlm[48,0] = 2
nlm[48,1] = 2
nlm[48,2] = 1
nlm[49,0] = 2
nlm[49,1] = 3
nlm[49,2] = 0
nlm[50,0] = 3
nlm[50,1] = 0
nlm[50,2] = 2
nlm[51,0] = 3
nlm[51,1] = 1
nlm[51,2] = 1
nlm[52,0] = 3
nlm[52,1] = 2
nlm[52,2] = 0
nlm[53,0] = 4
nlm[53,1] = 0
nlm[53,2] = 1
nlm[54,0] = 4
nlm[54,1] = 1
nlm[54,2] = 0
nlm[55,0] = 5
nlm[55,1] = 0
nlm[55,2] = 0

