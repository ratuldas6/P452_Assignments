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

def newtonraphson(func,x_0,prec): #roots using Newton-Raphson method
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
    x = x_0-func(x_0)/derivative(func, 1, x_0, prec*10**(-3))
    while abs(x_0-x) > prec and i < 15:
        #loop until precision AND iteration requirement is met
        x = x_0
        #applying the formula for N-R method
        x_0 = x_0 - func(x_0)/derivative(func, 1, x_0, prec)
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

def derivative(f, n, x, h=1e-3): #derivative finder
    #n decides 1st or 2nd order derivative of f(x)
    if int(n)==1:   #at x; h is the tolerance value
        d1f = round((f(x+h)-f(x))/h,2)
        return d1f
    elif int(n)==2:
        d2f = round((f(x+h)+f(x-h)-2*f(x))/(h**2),2)
        return d2f

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

def gauss_quadrature(func, a, b, n, h):                              # Code for Gauss-Legendre Quadrature based definite integration
    dxdt = (b - a)/2                                                    # defining interval scaling from [-1,1] to [a,b]

    guess_list = np.linspace(-1, 1, n)                                     # generating n partitions of interval [-1,1]
    root_list = [newtonraphson(func, x, 1e-6) for x in guess_list]           # finding roots of Legendre polynomial by providing the equal partitions of [-1,1] to be guess values 
    weight_list = [2/((1 - x**2) * (derivative(func, n, x, h))**2) for x in root_list]
                                                                        # obtaining weight functions for each root of the Legendre polynmial
    t_list = ((b - a)/2)*(array(root_list)) + (a + b)/2                 # scaling the interval [-1,1] to [a,b]
    sum = 0
    for i in range(len(weight_list)):
        sum += f(t_list[i]) * weight_list[i]                            # performing sum over w(t_i) * f(t_i) where t_i is in [a,b]

    return sum * dxdt

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
    