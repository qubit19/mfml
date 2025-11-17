def cholesky_decomposition(A):
    n = len(A)
    
    # Create empty L matrix
    L = [[0.0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(i + 1):  # Only compute lower triangle
            sum_val = sum(L[i][k] * L[j][k] for k in range(j))
            
            if i == j:  # Diagonal values
                val = A[i][i] - sum_val
                if val <= 0:
                    raise ValueError("Matrix is not positive definite.")
                L[i][j] = val ** 0.5
            else:  # Off-diagonal
                L[i][j] = (A[i][j] - sum_val) / L[j][j]
    
    return L


# Example usage
A = [
    [25, 15, -5],
    [15, 18,  0],
    [-5, 0, 11]
]

L = cholesky_decomposition(A)

print("L matrix:")
for row in L:
    print(row)
