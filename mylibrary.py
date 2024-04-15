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

def sumAP(a, d, n): #finds sum of an arithmetic progression upto n terms
    """
    Parameters:
    a : first term of AP
    d : difference between terms of AP
    n : number of terms in AP

    Returns:
    sum_terms : sum of all terms in AP
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
    Parameters:
    a : first term of GP
    d : difference between terms of GP
    n : number of terms in GP

    Returns:
    sum_g : sum of all terms in GP
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
    """
    Parameters:
    a : first term of HP
    d : difference between terms of HP
    n : number of terms in HP

    Returns:
    sum_h : sum of all terms in HP
    """
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
    """
    Parameters:
    f : function to find derivative for
    n : decides between 1st or 2nd order derivative ('1' for 1st, '2' for second)
    x : value of function to find derivative at
    h : tolerance

    Returns:
    d1f : returns value of 1st order derivative
    d2f : returns value of 2nd order derivative
    """

    #n decides 1st or 2nd order derivative of f(x)
    if int(n)==1:
        d1f = round((f(x+h)-f(x))/h,2)
        return d1f
    elif int(n)==2:
        d2f = round((f(x+h)+f(x-h)-2*f(x))/(h**2),2)
        return d2f

def dfdx_sym(P, n, x, dx = 10**(-6)):   # symmetric derivative of Legendre polynomials
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

def is_symmetric(M):
    for i in range(0,len(M)):
        for j in range(0,len(M)):
            if M[i][j] != M[j][i]:
                return False
    return True
        
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
        x[i] = (y[i]-sumj)/U[i][i]
    return x

def gauss_seidel(a,b,tol=0.000001,iter_count=100):
    L = []
    v = []

    for i in range(0,len(a)):
        L.append(0)                                       # make two zero row matrices that has a guess solution of 0
        v.append(0)

    k = 0
    e = 10

    while e > tol:                                      # setting epsilon from tolerance
        for i in range(0,len(a)):
            sum1 = 0
            sum2 = 0
            for j in range(0,i): 
                sum1 = sum1 + a[i-1][j-1]*L[j-1] 
            for j in range(i+1,len(a)):
                sum2 = sum2 + a[i-1][j-1]*L[j-1]
            v[i-1] = L[i-1]                                 # storing the values of previous iteration
            L[i-1] = (b[i-1] - sum1 - sum2)/a[i-1][i-1]
        e = (abs(L[3])-abs(v[3]))/abs(L[3])               # calculation of epsilon after each iteration
        k = k + 1
        #print(L)
        #print(k)
        if k > iter_count:                                # stopping the loop after 30 iterations
            return "The solution does not converge"
            break
    return L

def create_augment(z,x):
    for i in range(0,len(z)):
        z[i].append(x[i])
    return z

def max_swap(z,a):
  for i in range(a,len(z)): #for swapping(sorting) the rows with largest element and smallest element 
    for j in range(i+1,len(z)):
      if abs(z[i][a]) < abs(z[j][a]):
        m = z[i]
        z[i] = z[j]
        z[j] = m
      else:
        break

def norm_row(z,a):
  x = []  #for creating All the diagonal elements 1
  for i in range(0,len(z)+1):
    c = z[a][i]
    x.append(c/z[a][a])
  z[a] = x

def column(z,a):
  for i in range(a+1,len(z)):
    x = []
    if abs(z[i][a]) > 0:  #for making a elements 0 that are below the diagonal element
      for j in range(0,len(z[i])):
        c = z[i][j] - z[a][j]*z[i][a]
        x.append(c)
      z[i] = x
  
def reverse_column(z,a): #for making a elements 0 that are above the diagonal elemnt
  for i in range(a-1,-1,-1):
    if abs(z[i][a]) > 0:
      for j in range(len(z[i])-1,-1,-1):
        c = z[i][j]-z[a][j]*z[i][a]
        z[i][j] = c
  return z


def give_solutions(z):
    for i in range(0,len(z)): # the output matrix is an augumented matrix 
        print('for x',(i+1),end =" root is equal to ") # for printing the last element of each row.
        print(z[i][len(z)])

def gauss_jordan(A, b):
    n = len(b)
    # Augmenting the matrix A with vector b
    Ab = np.hstack((A, np.expand_dims(b, axis=1)))
    
    # Forward elimination
    for i in range(n):
        # Partial pivoting for numerical stability
        pivot_row = np.argmax(np.abs(Ab[i:, i])) + i
        Ab[[i, pivot_row]] = Ab[[pivot_row, i]]
        
        pivot = Ab[i, i]
        Ab[i] /= pivot
        for j in range(i + 1, n):
            factor = Ab[j, i]
            Ab[j] -= factor * Ab[i]
    
    # Backward elimination
    for i in range(n - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            factor = Ab[j, i]
            Ab[j] -= factor * Ab[i]
    
    return Ab[:, -1]

def lu_factorization(A):
    n = A.shape[0]
    U = np.copy(A)
    L = np.eye(n)
    P = np.eye(n)

    for i in range(n - 1):
        pivot_row = np.argmax(np.abs(U[i:, i])) + i
        if pivot_row != i:
            U[[i, pivot_row]] = U[[pivot_row, i]]
            P[[i, pivot_row]] = P[[pivot_row, i]]
            if i >= 1:
                L[[i, pivot_row]] = L[[pivot_row, i]]
        for j in range(i + 1, n):
            factor = U[j, i] / U[i, i]
            L[j, i] = factor
            U[j] -= factor * U[i]

    return L, U, P

def forward_substitution(L, b):
    n = len(b)
    y = np.zeros_like(b)
    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= L[i, j] * y[j]
    return y

def backward_substitution(U, y):
    n = len(y)
    x = np.zeros_like(y)
    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= U[i, j] * x[j]
        x[i] /= U[i, i]
    return x

def calculate_norm(r):
    """
    Parameters:
    r : an n-dimensional vector

    Returns:
    norm: the Euclidean norm of the vector r
    """
    norm_squared = sum(element ** 2 for element in r)
    norm = norm_squared ** 0.5
    return norm

def conjugate_gradient(A, b, x0, tol=10**(-4), max_iter=100):
    """
    Parameters:
    A : coefficient matrix with n rows and columns
    b : an n-dimensional vector
    x0 : initial guess vector
    tol : tolerance
    max_iter : max no. of iterations allowed

    Returns:
    x : solution vector 
    """
    x = np.array(x0, dtype=float)  # Initial guess
    r = b - np.dot(A, x)  # Initial residual
    p = r.copy()  # Initial search direction
    r_r = np.dot(r, r)
   
    
    for k in range(max_iter):
        Ap = np.dot(A, p)
        alpha = r_r / np.dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        
        r_r_next = np.dot(r, r)
        beta = r_r_next / r_r
        p = r + beta * p
        r_r = r_r_next
        
        if calculate_norm(r) < tol:
            break
        
    return x


def conjugate_on_fly(matrix, b, x0, tol=10**(-6), max_iter=100):
    """
    Parameters:
    matrix : an input matrix
    b : a vector
    x0 : initial guess vector
    tol : tolerance
    max_iter : max no. of iterations allowed

    Returns:
    x : solution vector
    residue_norms : final residue from calculation
    iter_count : maximum number of iterations allowed
    """
    iter_count = 0
    x = x0.copy()  # Initial guess
    r = b - matrix(x)  # Initial residual
    p = r.copy()  # Initial search direction
    residue_norms = [calculate_norm(r)]  # List to store residue norms

    for k in range(max_iter):
        Ap = matrix(p)
        alpha = np.dot(r, r) / np.dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        
        beta = np.dot(r, r) / np.dot(r - alpha * Ap, r - alpha * Ap)
        p = r - alpha * Ap + beta * p

        residue_norm = calculate_norm(r)
        residue_norms.append(residue_norm)
        iter_count = iter_count +1
        if residue_norm < tol:
            break

    return x, residue_norms, iter_count

def conjugate_inverse(matrix, b, x0, tol=10**(-6), max_iter=100):
    """
    Parameters:
    matrix : an input matrix
    b : a vector
    x0 : initial guess vector
    tol : tolerance
    max_iter : max no. of iterations allowed

    Returns:
    A_inv : inverse of the matrix
    """
    N = len(b)
    inverse_columns = []
    
    for i in range(N):
        # Create the right-hand side vector for solving Ax = e_i
        ei = np.zeros(N)
        ei[i] = 1
        
        # Solve the equation Ax = e_i using Conjugate Gradient method
        x, _, _ = conjugate_on_fly(matrix, ei, x0, tol, max_iter)
        
        # Append the solution (column of the inverse matrix) to the list
        inverse_columns.append(x)
    
    # Stack the columns of the inverse matrix horizontally to form the complete inverse matrix
    A_inv = np.column_stack(np.round(inverse_columns,4))
    return A_inv

#----------------------------------------------------------------------------------------------------
#Root-finding function set
def bracket(func, a, b): #defines the bracket
    """
    Parameters:
    func : function in question
    a : lower limit of interval
    b : upper limit of interval

    Returns:
    a : 'a' after adjustment
    b : 'b' after adjustment
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

def divideSyn(coeff, guess): #synthetic division function
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

def power_iter(A, num_iterations=1000, tol=1e-6):
    n = A.shape[0]
    x = np.random.rand(n)  # Random initial guess for the eigenvector

    for _ in range(num_iterations):
        x1 = np.dot(A, x)
        eigenvalue = np.linalg.norm(x1)
        x1 = x1 / eigenvalue  # Normalize the eigenvector estimate
        if np.linalg.norm(x - x1) < tol:
            break
        x = x1

    return eigenvalue, x

def qr(A):
    """
    Parameters:
    A : matrix for performing QR factorization using Gram-Scmidt orthogonalization
    
    Returns:
    Q : matrix Q for QR
    R : matrix R for QR
    """
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))

    for j in range(n):
        v = A[:, j]
        for i in range(j):
            R[i, j] = np.dot(Q[:, i], A[:, j])
            v = v - R[i, j] * Q[:, i]
        R[j, j] = np.linalg.norm(v)
        Q[:, j] = v / R[j, j]

    return Q, R

def eigenvalue(A, iterations = 50000):
    A_k = np.copy(A)
    n = A.shape[0]
    QQ = np.eye(n)
    for k in range(iterations):
        Q, R = qr(A_k)
        A_k = R @ Q
        ev = []
    for i in range (0, A_k.shape[0]):
        ev.append(A_k[i, i])
    ev = np.array(ev)
    return ev, A_k, QQ


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
    
#----------------------------------------------------------------------------------------------------
#Least-square fit function set

def poly_basis(x):
    return np.array([np.ones_like(x), x**1, x**2, x**3])

def cheby_basis(x):
    return np.array([np.ones_like(x), 2*x - 1, 8*(x**2)-8*x+1, 32*(x**3)-48*(x**2)+18*x-1])

def poly_fit(xlist, ylist, basis):
    avg = 0
    for x in (ylist):
        avg += x
    avg = avg/(xlist.shape)
    lhs = basis(xlist) @ basis(xlist).T
    rhs = basis(xlist) @ ylist.T
    par = np.linalg.inv(lhs)@rhs
    return par, np.linalg.cond(lhs)

#----------------------------------------------------------------------------------------------------
#Random number function set

def linear_congruential_gen(a=587982351, c=54312, m = 2**32, n=10):
    """
    Parameters:
    a : seed value
    c : constant to go along with seed value
    m : period parameter
    n : no. of random nos. required
    
    output : a 1-d array containing the lower and the upper bound of the random numbers generated.
    """
    x0 = 123
    arr = np.zeros((n))
    arr[0] = x0

    for i in range (1, n):
        arr[i] = ((a*arr[i-1] + c)%m)

    output = arr/m

    return output

def monte_carlo_lcg(func, a, b, step, lc_a = 1103515245, lc_c = 12345, lc_m = 2**32):
    x = np.array(linear_congruential_gen(a = lc_a, c = lc_c, m = lc_m, n = step))
    x = a + (b-a)*x
    my_func = np.vectorize(func)
    y = my_func(x)
    return (b-a)*(np.sum(y))/step

class MLCG:
    def __init__(self, seed = 30):
        self.state = seed

    def rand(self, coeff = 979, multiplier = 32311):
        self.state = (coeff * self.state) % multiplier
        return (self.state)/multiplier
    
    def rand_range(self, a= 0, b= 1): 
        return (b - a) * self.rand() + a
    
    def rand_array(self, array_size):
        return np.array([self.rand() for _ in range(array_size)])
    
    def rand_range_array(self, a, b, array_size):
        return np.array([self.rand_range(a, b) for _ in range(array_size)])

def mc_importance_sampling(imp_samples, int_func, pd_func):
    int_sum = 0
    ind_val_list = []
    for val in imp_samples:
        ind_val = int_func(val)/pd_func(val)
        int_sum += ind_val
        ind_val_list.append(ind_val)

    mean = int_sum/ len(imp_samples)
    variance = (sum([(val - mean)**2 for val in ind_val_list])/len(imp_samples)) ** 0.5
    return mean, variance

def inv_transform_sampling(a, b, N, inverse_trans_func):
    mlcg = MLCG()
    return inverse_trans_func(mlcg.rand_range_array(a, b, N))

def accept_reject_sampling(sample_size, low_lim, up_lim, target_dist, samp_dist):
    sample_array = []
    mlcg = MLCG()

    k = 0

    while len(sample_array) < sample_size:        
        x = mlcg.rand_range(low_lim, up_lim)
        u = mlcg.rand_range(low_lim, up_lim)
        if u <= target_dist(x)/((1/3) * samp_dist(x)):
            sample_array.append(x)
        
    return np.array(sample_array)