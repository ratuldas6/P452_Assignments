import numpy as np

#----------------------------------------------------------------------------------------------------
#General function set
def integerRound(x):
    from math import floor, ceil
    a = floor(x)
    b = ceil(x)
    if x-a<=.1:
        result = a
    elif b-x<=.1:
        result = b
    else:
        print("Please apply the function to a smaller difference")
    return result

def dec2Round(x): #rounds a number up to2 decimal places
    y = 100*x
    remainder = y%1
    y = y - remainder
    if remainder > 0.1:
        modified_digit = 1
    else:
        modified_digit = 0
    y = y + modified_digit
    y = y/100
    return y

def sumAP(a,d,n): #finds sum of an arithmetic progression upto n terms
    """
    a : first term of AP
    d : difference between terms of AP
    n : number of terms in AP
    """
    n = float(n)
    if n >= 1 and float(n)==int(n):
        n = int(n)    # Alotting specific data types
        a = float(a)
        d = float(d)
        sum_terms = 0
        for i in range(1,n+1):
            sum_terms += a
            a += d
        return sum_terms
    else:
        print("Invalid input. Please enter an integer for n!")
        return None
    
def sumGP(a,r,n): #finds sum of a geometric progression upto n terms
    """
    a : first term of AP
    d : difference between terms of AP
    n : number of terms in AP
    """
    n = float(n)
    if n >= 1 and float(n)==int(n):
        n = int(n)
        a = float(a)
        r = float(r)
        sum_g = 0
        for i in range(1,n+1):
            sum_g += a
            a = a*r
        return sum_g
    else:
        print("Invalid input. Please enter an integer for n!")
        return None
    
def sumHP(a,d,n): #finds sum of a harmonic progression upto n terms
    n = float(n)
    if n >= 1 and float(n)==int(n):
        n = int(n)
        a = float(a)
        d = float(d)
        sum_h = 0
        for i in range(1,n+1):
            a = 1/a
            sum_h += a
            a = 1/a
            a = a + d     
        return sum_h
    else:
        print("Invalid Input!")
        return None
    
def derivative(f, n, x, h=1e-3): #derivative finder
    #n decides 1st or 2nd order derivative of f(x)
    if int(n)==1:   #at x; h is the tolerance value
        d1f = round((f(x+h)-f(x))/h,2)
        return d1f
    elif int(n)==2:
        d2f = round((f(x+h)+f(x-h)-2*f(x))/(h**2),2)
        return d2f

def dfdx_sym(P, n, x, dx = 10**(-6)):                                   # symmetric derivative of Legendre polynomials
    return (P(n, x + dx) - P(n, x - dx))/(2 * dx)

def P(n, x):
    # for defining Legendre polynomial
    if n < 0:
        print("provide positive order of Legendre Polynomial Pn(x)!")
        return None
    elif n == 0:                                                        
        return 1
    elif n == 1:                                                        
        return x
    
    func = ((2*n - 1)/n) * x * P(n-1, x) - ((n-1)/n) * P(n-2, x) 

    return func

def showConvergenceTable(steps, values):
    print("# of iterations (i)     ","Absolute Error (|b-a|)\n")
    for j in range(len(steps)):
        print("        ", steps[j], "             ", values[j])

#----------------------------------------------------------------------------------------------------
#Matrix function set
def displayMatrix(matrix): #display given matrix
    for row in matrix:
        print(row)
        
def storeInvMatrix(matrix): #display inverse GJ eliminated matrix without the identity in the first half
    r = len(matrix)
    c = len(matrix)
    new_mat = []
    for i in range(r):
        new_mat.append(matrix[i][c:2*c])
    return new_mat

def transpose(A):
    n = len(A)
    A_T = [[0 for i in range(n)] for j in range(n)]
    for i in range(len(A)):
        for j in range(len(A[0])):
            A_T[j][i] = A[i][j]
    return A_T               
        
def applyGJ(A,B): #function to apply GJ elim
    for k in range(len(B)):
        #rows are pivotted
        if abs(A[k][k]) < 1.0e-6:
            #defining upper limit for element = 0
            for i in range(k+1, len(B)): 
                if abs(A[i][k]) > abs(A[k][k]):
                    #swap check
                    for j in range(k, len(B)):
                        #obtaining swapped rows
                        A[k][j], A[i][j] = A[i][j], A[k][j] 
                    B[k], B[i] = B[i], B[k] 
                    break
        A_kk = A[k][k]
        if A_kk == 0:        
            print("No distinct solution was found for this system of equations")
            return
        for j in range(k, len(B)): #column index of row is pivotted
            A[k][j] /= A_kk         #pivot row division for row ech form
        B[k] /= A_kk
        for i in range(len(B)):    #changed rows assigned new indices
            if i == k or A[i][k] == 0: continue
            factor = A[i][k]
            for j in range(k, len(B)): 
                #columns for subtraction assigned indices
                A[i][j] -= factor*A[k][j]
            B[i] -= factor * B[k]
    return B

def productMatrix(B,A): #finds product of two matrices
    try:               
        if  len(A) != len(B[0]): 
            print("Multiplication is undefined for given matrices!") 
        else:
            C = [[0 for i in range(len(A[0]))] for j in range(len(B))]
            for i in range(len(B)):       #rows of the matrix product.
                for j in range(len(A[0])):#columns of the matrix product.
                    for k in range(len(A)):
                        C[i][j] += dec2Round(B[i][k]*A[k][j])
            return C
    except TypeError:
        print("Invalid entries found in matrix!")
        return None

def findDeterminant(A):     #determinant function for use in applyGJ()
    if len(A) != len(A[0]): #square matrix check
        print("Determinant is undefined for square matrices")
    else:
        count = 0
        for i in range(len(A) - 1): 
            if abs(A[i][i]) < 1.0e-6:
                for j in range(i+1 , len(A)): 
                    if  abs(A[i][i]) < abs(A[j][i]): 
                        for k in range(i , len(A)):
                            #swapping the rows of the matrix
                            A[i][k], A[j][k] = A[j][k], A[i][k]
                            count += 1
            for j in range(i+1 , len(A)):
                if A[j][i] == 0: continue 
                result = A[j][i]/A[i][i] #dividing the rows.
                for k in range(i , len(A)):
                    A[j][k] -= result*A[i][k]
        initial = 1
        for j in range(len(A)):
            initial *= A[j][j] #product of diagonal elements of matrix
        initial *= (-1)**count
        print(initial) 

def inverseMatrix(A):
    #function for inverse matrix (n x n)
    M = [[ 0.00 for i in range(len(A))] for j in range(len(A))] 
    for i in range(len(A)):     #run loop in rows
        for j in range(len(A)): #run loop in columns
            M[j][j] = 1.00
    for i in range(len(A)):
        A[i].extend(M[i])
    for k in range(len(A)):
        if abs(A[k][k]) < 1.0e-6:
            #GJ elimination segment of the function
            for i in range(k+1, len(A)):
                if abs(A[i][k]) > abs(A[k][k]):
                    for j in range(k,2*len(A)):
                        #swap rows
                        A[k][j], A[i][j] = A[i][j], A[k][j]
                    break
        count = A[k][k] #element is pivotted
        if count == 0:  #checking if pivot = 0
            print("The matrix does not have a defined inverse")
            return
        else:
            for j in range(k, 2*len(A)):
                #pivotted row columns
                A[k][j] /= count
            for i in range(len(A)):
                #substracted rows indiced
                if i == k or A[i][k] == 0: continue
                result = A[i][k] 
                for j in range(k, 2*len(A)): 
                    #columns for subtraction indiced
                    A[i][j] -= result*A[k][j]
                    
    A_inv = []
    
    for i in range(len(A)):
        blank_row = []
        for j in range(len(A),len(A[0])):
            blank_row.append(A[i][j])
        A_inv.append(blank_row)
    return A_inv

#----------------------------------------------------------------------------------------------------
#Advanced matrix function set

def swapRows(M, row_old, row_new, columns):
    #to swap rows, wherever needed
    tmp = []
    for i in range (0, int(columns)):
        tmp.append(M[int(row_old)][i])
        M[int(row_old)][i] = M[int(row_new)][i]
        M[int(row_new)][i] = tmp[i]
        
def unSwap(A,B,A_orig):
    #to undo swapping from inversion
    N = len(A)
    fix_list = []
    B_fixed = []
    
    for i in range(N):
        for j in range(N):
            if A[i]==A_orig[j]:
                fix_list.append(j)
    
    for elem in fix_list:
        B_fixed.append(B[elem])
        
    B = B_fixed
        
def crout(A):
    n = len(A)
    L = [[0.0 for i in range(n)] for j in range(n)]
    U = [[0.0 for i in range(n)] for j in range(n)]
    
    for i in range(n):
        U[i][i] = 1
        for j in range(i, n):
            tmp1 = A[j][i] 
            for k in range(i):
                tmp1 -= L[j][k]*U[k][i]
            L[j][i] = tmp1
        for j in range(i+1, n):
            tmp2 = A[i][j]
            for k in range(i):
                tmp2 -= L[i][k]*U[k][j]
            U[i][j] = tmp2/L[i][i]
    return (L, U)

def doolittle(A):
    n = len(A)
    L = [[0.0 for i in range(n)] for j in range(n)]
    U = [[0.0 for i in range(n)] for j in range(n)]
    
    for z in range(n):
        L[z][z] = 1
        f1 = 0
        f2 = 0
        f3 = 0
        for p in range(z):
            f1 = f1 + L[z][p]*U[p][z]
        U[z][z] = (A[z][z] - f1)
        for i in range(z+1, n):
            for p in range(z):
                f2 = f2 + L[z][p]*U[p][i]
            U[z][i] = (A[z][i] - f2)
        for k in range(z+1, n):
            for p in range(z):
                f3 = f3 + L[k][p]*U[p][z]
            L[k][z] = (A[k][z] - f3)/U[z][z]
    return (L, U)

def forwardSolve(L, b): #solves L.y = b
    y = [0 for i in range(len(b))]
    for i in range(len(b)):
        sumj = 0
        for j in range(i):
            sumj = sumj + L[i][j]*y[j]
        y[i] = (b[i] -sumj)/L[i][i]
    return y

def backwardSolve(U, y): #solves U.x = y
    n = len(y)
    x = [0 for i in range(len(y))]
    for i in range(n-1, -1, -1):
        sumj = 0
        for j in range(i+1, n):
            sumj = sumj + U[i][j] * x[j]
        x[i] = (y[i] - sumj)/U[i][i]
    return x 

def partPivot(M, m, r, c):
    #to get partially pivotted matrix for LU decomp
    global n, swaps
    n = 0
    swaps = 0
    pivot = M[int(m)][int(m)]
    for i in range (int(r)):         
        if pivot < M[int(i)][int(m)]:
            #same-column elements are checked
            pivot = M[int(i)][int(m)]
            n += 1
            swaps = i
    if swaps != 0: #swap if allowed
        swapRows(M, m, swaps, c)
    if int(pivot) == 0:
        print ("No unique solution")
        return None

def solveLU(A, b, method):
    if method==1:
        L, U = crout(A)
    elif method==2:
        L, U = doolittle(A)
    y = forwardSolve(L, b)
    x = backwardSolve(U, y)
    return x

def LU_inverse(M):
    #to find inverse of matrix using LU decomposition
    M_orig = [[M[j][i] for i in range(len(M))] for j in range(len(M))]
    
    I = [[0.00 for i in range(len(M))] for j in range(len(M))] 
    for i in range(len(M)):
        for j in range(len(M)):
            #creates an identity matrix
            I[j][j] = 1.00

    if M[1][1] == 0 and M[0][1] != 0:
        #swaps rows of first submatrices to prevent det=0
        swapRows(M, 0, 1, 4)

    partPivot(M, 0, len(M), len(M[0]))
    L_i, U_i = doolittle(M)
    
    y = [[0 for c in range(len(I))] for r in range(len(I))]
    for i in range(len(I)):
        for k in range (len(I[0])):
            y[i][k] = I[i][k]
            for j in range(i):
                y[i][k]=y[i][k]-(L_i[i][j]*y[j][k])
            y[i][k] = y[i][k]/L_i[i][i]
            
    n = len(y)
    x = [[0,0,0,0] for r in range(len(I))]
    if U_i[n-1][n-1] == 0: 
        #check for diagonal elements = 0
        raise ValueError
    
    for i in range(n-1, -1, -1):
        for k in range (len(I[0])): 
            x[i][k] = y[i][k]
            for j in range(i+1,n):
                x[i][k] = x[i][k] -(U_i[i][j]*x[j][k])
            x[i][k] = x[i][k]/U_i[i][i]
    
    x = transpose(x)
    unSwap(M, x, M_orig)
    x = transpose(x)
    
    return(x)

def choleskydecomp(A):
    from math import sqrt
    #finds L using Cholesky's algorithm
    n = len(A) #null matrix for L
    L = [[0.0] * n for i in range(n)]
    for i in range(n):
        for k in range(i+1):
            tmp_sum = sum(L[i][j] * L[k][j] for j in range(k))
            if (i == k): # Diagonal elements
                L[i][k] = sqrt(A[i][i] - tmp_sum)
            else:
                L[i][k] = (1.0 / L[k][k] * (A[i][k] - tmp_sum))
    return L
    
def solveCholesky(L, U, b):
    n = len(L)
    y = [0 for i in range(n)]
    x = [0 for i in range(n)]
    for i in range(n):
        sumj = 0
        for j in range(i):
            sumj += L[i][j]*y[j]
        y[i] = (b[i]-sumj)/L[i][i]
    for i in range(n-1, -1, -1):
        sumj = 0
        for j in range(i+1, n):
            sumj += U[i][j]*x[j]
        x[i] = dec2Round((y[i]-sumj)/U[i][i])
    return x

#----------------------------------------------------------------------------------------------------
#Root-finding function set
def bracket(func, a, b): #defines the bracket
    """
    func : function in question
    a : lower limit of interval
    b : upper limit of interval
    """
    beta = 0.05
    i = 0 #number of iterations
    if func(a)*func(b) > 0:
        while i < 12 and func(a)*func(b) > 0:
            #loop until sign AND iteration requirement is met
            if abs(func(a)) < abs(func(b)):
                #a is moved slightly to left if root is in (a,b)
                a = a - beta*(b-a)
                i += 1
            elif abs(func(b)) < abs(func(a)):
                #b is moved slightly to right if root is in (a,b)
                b = b + beta*(b-a)
                i += 1    
        if i == 12:
            return "Root not found. Please try a different interval."
        else:
            return a, b  
    else:
        return a, b

def divideSyn(coeff,guess): #synthetic division function
    ans = []   #coefficients of Q(x)
    temp = []  #creating the matrix/addend to be added
    temp.append(0)
    for i in range(len(coeff)-1):
        #loop through each coeff
        ans.append(coeff[i]+temp[i])        
        temp.append(guess*ans[i]) #loops until second last place
    return ans

def bisection(func, a, b, prec): #roots using bisection method
    """
    func : function in question
    a : lower limit of interval
    b : upper limit of interval
    prec : precision value for finding root
    """
    i = 0
    list_i = []
    list_f_i = []
    abs_err = []
    if func(a)*func(b)>0: #check for proper bracketing
        return("Missing","Proper","Bracketing",".")
    else: #go ahead with brackets containing root
        while abs(b-a) > prec and i < 15:
            #loop until precision AND iteration requirement is met
            abs_err.append(abs(b-a))
            c = (a + b)/2
            if func(c) == 0:
                break
            if func(c)*func(a)<0:
                b = c
            else:
                a = c
            i += 1
            list_i.append(i)
            list_f_i.append(func(c))
    return c, list_i, list_f_i, abs_err

def regulafalsi(func, a, b, prec): #roots by Regula-Falsi method
    """
    func : function in question
    a : lower limit of interval
    b : upper limit of interval
    prec : precision value for finding root
    """
    i = 0
    list_i = []
    list_f_i = []
    abs_err = []
    c = a
    d = b
    if func(a)*func(b)>0: #check for proper bracketing
        return("Missing","Proper","Bracketing",".")
    else:
        while abs(d-c) > prec and i < 15:
            #loop until 15 steps or until root is found, whichever is earlier
            d = c
            c = b - ((b-a)*func(b))/(func(b)-func(a))    
            if func(a)*func(c)<0:
                b = c
            else:
                a = c
            i += 1
            list_i.append(i)
            list_f_i.append(func(c))
            abs_err.append(abs(d-c))                
    return c, list_i, list_f_i, abs_err

def newtonraphson(func, x_0, prec): #roots using Newton-Raphson method
    """
    f : function in question
    x_0 : Initial guess
    prec : precision value for finding root
    """
    i = 0
    list_i = []
    abs_err = []
    if abs(func(x_0)) == 0:
        return x_0
    x = x_0-func(x_0)/dfdx_sym(func, 1, x_0)
    while abs(x_0-x) > prec and i < 15:
        #loop until precision AND iteration requirement is met
        x = x_0
        #applying the formula for N-R method
        x_0 = x_0 - func(x_0)/dfdx_sym(func, 1, x_0)
        if func(x_0) == 0:
            break
        i+=1
        list_i.append(i)
        abs_err.append(abs(x_0-x))
    return x_0, list_i, abs_err

def findRoot(coeff, degree, alpha, prec): #in Laguerre, find root by taking a guess as input
    from math import sqrt
    #coeff = list of coefficients for polynomial, deg = degree of poly., alpha = guess, err = tolerance
    f = lambda x: sum([coeff[i]*x**(degree-i) for i in range(degree+1)])
    roots = []
    while True:  #break if alpha is a root
        if abs(f(alpha)) < prec: #compare with prec to decide if alpha is the root
            alpha, l1, l2 = newtonraphson(f, alpha, prec)
            ans = divideSyn(coeff, alpha) #performing synthetic division for +ve response
            return alpha, ans
            break
        else: #adjusting the root
            d1f = derivative(f,1,alpha,prec) #value of 1st derivative of f at alpha
            d2f = derivative(f,2,alpha,prec) #value of 2nd derivative of f at alpha
            G = d1f/f(alpha)
            H = G**2 - d2f/f(alpha)
            g1 = G + sqrt((degree-1)*(degree*H-G**2))
            g2 = G - sqrt((degree-1)*(degree*H-G**2))
            sq = max(g1, g2) #pick the larger abs value for the denominator
            alpha = alpha-degree/sq        

def laguerre(coeff, alpha, prec): #roots by Laguerre's method
    """
    coeff : list of coefficients of polynomial in question,
            corresponding to greatest to least degree
    alpha : initial guess
    prec : precision value for finding root
    """
    degree = len(coeff)-1  #degree of given polynomial
    roots = []             #declaring empty list of roots
    while degree > 1:
        #keep looping until linear polynomial is reached
        alpha, coeff = findRoot(coeff, degree, alpha, prec)
        #storing new alpha and coeff
        roots.append(integerRound(alpha)) #appending roots
        degree = degree - 1
    roots.append(integerRound(-coeff[1]/coeff[0]))
    #extracting last root from remaining linear part of polynomial
    return roots


#----------------------------------------------------------------------------------------------------
#Numerical integration function set
def MonteCarlo(a,b,f,N):  #Monte-Carlo method for integration
    from random import uniform
    """
    a : lower bound for integration
    b : upper bound for integration
    f : integrand
    N : number of iterations for integrating
    result : value of the integral using Monte-Carlo method
    """
    count_int = 0       #number of points contributing to the integral
    
    for i in range(N):
        x = uniform(a,b) #value between a and b is chosen randomly
        if f(a)>f(b):    #range of y is obtained
            f_maximum = f(a)
            y = uniform(0,f(a))
        elif f(b)>f(a):
            f_maximum = f(b)
            y = uniform(0,f(b))
        area = (f_maximum-0)*(b-a)
        if y <= f(x):    #condition for contribution to integral satisfied
            count_int += 1
            
    result = (count_int/N)*area   #result of the integration
    return result, N

def Midpoint(a,b,f,N):  #Midpoint method for integration
    """
    a : lower bound for integration
    b : upper bound for integration
    f : integrand
    N : number of iterations for integrating
    integral : value of the integral using Midpoint method
    """
    h = (b-a)/N         #height of rectangular element
    integral = 0
    
    for index in range(N):
        integral += h*(f(a+index*h)+f(a+(index+1)*h))/2
    return integral

def Trapezoid(a,b,f,N): #Trapezoidal method for integration
    """
    a : lower bound for integration
    b : upper bound for integration
    f : integrand
    N : number of iterations for integrating
    integral : value of the integral using Trapezoidal method
    """
    h = (b-a)/N         #defining the trapezoidal width
    integral = 0
    
    for index in range(N+1):
        if index==0 or index==N:         #for the first and last term, weight function = 1 or h/2
            integral += h*f(a+index*h)/2
        else:
            integral += h*f(a+index*h)   #for other terms, weight function = 2 or h
    return integral

def Simpson(a,b,f,N):   #Simpson method for integration
    """
    a : lower bound for integration
    b : upper bound for integration
    f : integrand
    N : number of iterations for integrating
    integral : value of the integral using Simpson method
    """
    h = (b-a)/N
    integral = 0
    
    for index in range(N+1):
        if index == 0 or index == N:         #for first and last term, weight = 1
            integral += h*f(a + index*h)/3 
        elif index%2 == 1:                   #for odd indices, weight = 4
            integral += 4*h*f(a + index*h)/3     
        else:
            integral += 2*h*f(a + index*h)/3 #for even indices, weight = 2
    return integral

def gauss_quadrature(func, a, b, n, h):                              # Code for Gauss-Legendre Quadrature based definite integration
    guess_list = np.linspace(-1, 1, n)                                     # generating n partitions of interval [-1,1]
    root_list = [newtonraphson(P, x, 1e-6) for x in guess_list]           # finding roots of Legendre polynomial by providing the equal partitions of [-1,1] to be guess values 
    weight_list = [2/((1 - x**2) * (dfdx_sym(P, n, x, h))**2) for x in root_list]
                                                                        # obtaining weight functions for each root of the Legendre polynmial
    t_list = ((b - a)/2)*(np.array(root_list)) + (a + b)/2                 # scaling the interval [-1,1] to [a,b]
    sum = 0
    for i in range(len(weight_list)):
        sum += func(t_list[i]) * weight_list[i]                            # performing sum over w(t_i) * f(t_i) where t_i is in [a,b]

    return sum * ((b - a)/2)

#----------------------------------------------------------------------------------------------------
#Differential equation function set

def LagrangeInterpolation(zh, zl, yh, yl, y): #function for Lagrange Interpolation
    z = zl + (zh - zl)*(y - yl)/(yh - yl)
    return z

def Euler(h, x, x_lim, y, func): #function for Euler method
    """
    h : step size for RK4
    x : initial value for x
    x_lim : final value for x
    y : value of y(0)
    func : dy/dx (function of both x and y, but mandatorily x)
    
    X : array of values for x
    Y : solution curve values
    """
    X = []
    Y = [] 
    
    while x <= x_lim: #compute y for x < x_lim
        y = y + h*func(x, y)
        x = x + h
        X.append(x)
        Y.append(y) 
    return X, Y

def Predictor(x, x_lim, y, h, f): #function for Predictor method
    """
    x : lower bound
    x_lim : upper bound
    y : value of y at initial value of x
    h : step size
    f : function of x and f, also dy/dx
    
    X : array of values for x
    Y : solution curve values
    """
    X = [x]
    Y = [y]
    while x < x_lim:                        
        y_1 = y + h*f(x, y)
        y += h*0.5*(f(x, y) + f(x + h, y_1))
        x += h  
        Y.append(y)
        X.append(x)
    return X, Y

def RK4(h, x, x_lim, y, func): #code for Runge-Kutta method
    """
    h : step size for RK4
    x : initial value for x
    x_lim : final value for x
    y : value of y(0)
    func : 1st derivative of y wrt x
    
    X : array of values for x
    Y : solution curve values
    """
    X = []
    Y = []
    
    if x < x_lim: #start iff initial x is less than x_lim
        while x <= x_lim: #iterate until x reaches x_lim
            #increments for func (ki) calculated
            k1 = h*func(x, y)
            k2 = h*func(x + h/2, y + h*k1/2)
            k3 = h*func(x + h/2, y + h*k2/2)
            k4 = h*func(x + h, y + h*k3)
            
            #update y and x
            y = y + (h*(k1 + 2*k2 + 2*k3 + k4))/2
            x = x + h
            
            #produce x-y pairs from x to x_lim
            X.append(x)
            Y.append(y)
            
        return X, Y 
    
    else:
        e = "Please keep the upper bound for x higher than the lower bound"
        return e
    