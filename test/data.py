#!/usr/bin/env python

import numpy

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

