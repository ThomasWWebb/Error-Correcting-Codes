import random
import numpy
import math
        
def randomMessage(r):
    length = 2**r - r -1 
    m = [None] * length
    for stepper in range(0, length):
        m[stepper] = round(random.uniform(0, 1))
    return m

#function HammingG
#input: a number r
#output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code
def hammingGeneratorMatrix(r):
    n = 2**r-1    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    
    return G


#function decimalToVector
#input: numbers n and r (0 <= n<2**r)
#output: a string v of r bits representing n
def decimalToVector(n,r): 
    v = []
    for s in range(r):
        v.insert(0,n%2)
        n //= 2
    return v

def encoder(m):
    m = numpy.matrix(m)
    r = math.floor(math.log2(m.size)+ 1)    
    G = numpy.matrix(hammingGeneratorMatrix(r))
    I = numpy.identity(2**r-1)
    ones = numpy.ones(((2**r -1), 1), dtype = int)
    M = numpy.append(I, ones, axis = 1)
    g = G*M
    C = (m*g) % 2
    C = C.astype(int)
    c = []
    for stepper in range(0,C.size):
        c.append(C[0,stepper])    
    return c


def BSC(c,p):
    c = numpy.matrix(c)
    for stepper in range(0,c.size): 
        chance = random.random()
        if (chance <= p):
            c[0,stepper] = (c[0,stepper] + 1) % 2
    V = c
    v = []
    for stepper in range(0,V.size):
        v.append(V[0,stepper ])
    return v

def syndrome(v):
    hasFailed = False
    v = numpy.matrix(v)
    r = math.floor(math.log2(v.size))
    colOnes = numpy.ones(((2**r), 1), dtype = int)
    rowZeros = numpy.zeros((1, r), dtype = int)
    
    h = []
    for i in range(1,2**r):
        h.append(decimalToVector(i,r))
    H = numpy.matrix(h)
    
    H = numpy.append(H, rowZeros, axis = 0)
    H = numpy.append(H, colOnes, axis = 1)
    
    synd = (v*H)% 2
    print("syndrome = ")
    s = []
    for stepper in range(0,synd.size):
        s.append(synd[0,stepper ])
    print(s)

    fullSyndDecimal = 0 
    for stepper in range(synd.size -1, -1,-1):
        if (synd[0, stepper] == 1):
            fullSyndDecimal = fullSyndDecimal + 2**(r - (stepper))
    print("syndrome = %d"  % (fullSyndDecimal))
    
    
    syndDecimal = 0 
    for stepper in range(synd.size -2, -1,-1):
        if (synd[0, stepper] == 1):
            syndDecimal = syndDecimal + 2**(r - (stepper +1))
    print("i = %d" % (syndDecimal))
    
    if not (syndDecimal == 0):
        v[0,syndDecimal -1] = (v[0,syndDecimal-1] + 1) % 2

    V = []
    for stepper in range(0,v.size):
        V.append(v[0,stepper ])
        
    fail = (V*H)% 2
    for stepper in range(0,fail.size):
            if not (fail[0, stepper] == 0):
                hasFailed = True
    if hasFailed:
        f = []
        for stepper in range(0,v.size - 1):
            f.append(0)
        f.append(1)
        #if a failure has occurred it will always return a non-codeword of [0,0,...,1] which will signify failure
        return f
    else:
        return V
                
def retrieveMessage(c):
    c = numpy.matrix(c)
    m = []
    for stepper in range(0,c.size):
        position = stepper + 1
        if not(math.log2(position) == math.floor(math.log2(position))):
            m.append(c[0,stepper ])
    return m
            
          
          

def simulation(r,N,p):
    success = 0
    failure = 0
    error = 0
    notFailed = False
    print(" r = %d; N = %d; p = %s" % (r, N, p))
    for i in range(1,N+1):
        notFailed = False
        print("Experiment %d of %d" % (i,N))
        print("")
        
        print("Source:")
        print("Message")
        m = randomMessage(r)
        print(m)
        print("")
        print("Codeword")
        c = encoder(m)
        print(c)
        print("")
        
        print("Channel:")
        print("Received vector")
        v =BSC(c,p)
        print(v)
        print("")

        print("Destination:")
        print("Decoding by syndrome")
        c = syndrome(v)

        C = numpy.matrix(c)  
        for stepper in range(0,C.size - 1):
            if not (C[0, stepper] == 0):              
                notFailed = True
        if not (C[0, C.size -1 ] == 1):
            notFailed = True
        if notFailed:
            
            print("")
            print("Codeword estimate")
            print(c)
            print("")

            print("Message estimate")
            M = retrieveMessage(c)
            print(M)
            print("")
            m = numpy.matrix(m)
            n = []
            for i in range(0,m.size):
                n.append(m[0,i])
           
            if (n == M):
                success += 1
            else :
                error +=1
        else :
            if not notFailed:
                failure += 1
            print("A decoder failure has occurred")
            print("")
        
        
    print("successes: %d" % (success))
    print("failures: %d" % (failure))
    print("errors: %d" % (error))

    print("")
    dep = error / N
    print("Experimental DEP:")
    print(round(dep, 4)* 100)
       


    
                     
