import numpy as np

def thomas_solve(A, D):
    '''
    Solves a tridiagonal system of equation of n variables using the Thomas Algorithm
    '''
    if(A.shape[0] != D.shape[0] or A.shape[1] != D.shape[0]):
        print("Invalid matrix passed to Thomas algorithm!")
    
    cdashes = []
    ddashes = []
    n = A.shape[0]
    for i in range(n):
        a = 0
        b = 0
        c = 0
        d = 0
        if(i != 0):
            a = A[i][i - 1]
        
        b = A[i][i]

        if(i != n - 1):
            c = A[i][i + 1]
        
        d = D[i]

        cdash = (c / b) if (i == 0) else (c / (b - a*cdashes[-1]))
        ddash = (d / b) if (i == 0) else (d - a*ddashes[-1]) / (b - a*cdashes[-1])

        if(i != n-1):
            cdashes.append(cdash)
        ddashes.append(ddash)
    
    xvals = []
    for i in range(n):
        x = 0
        if(i == 0):
            x = ddashes[-1]
        else:
            x = ddashes[n - i - 1] - cdashes[n - i - 1] * xvals[-1]
        xvals.append(x)
    
    xvals.reverse()
    return xvals

if __name__ == "__main__":
    A = np.array([(2,1,0), (1,2,1), (0, 1, 2)])
    D = np.array([1,2,3])

    print(thomas_solve(A, D))