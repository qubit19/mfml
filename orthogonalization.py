def dot(u, v):
    return sum(ui * vi for ui, vi in zip(u, v))

def scalar_mult(c, v):
    return [c * vi for vi in v]

def vector_sub(u, v):
    return [ui - vi for ui, vi in zip(u, v)]

def gram_schmidt(V):
    """
    V = list of input vectors (each vector is a list)
    Returns orthonormal basis vectors
    """
    U = []  # Orthonormal vectors

    for v in V:
        # Start with the original vector
        u = v[:]

        # Subtract projections onto previous U vectors
        for e in U:
            proj_coeff = dot(v, e)
            proj = scalar_mult(proj_coeff, e)
            u = vector_sub(u, proj)

        # Normalize
        norm = (dot(u, u)) ** 0.5
        if norm == 0:
            raise ValueError("Input vectors are linearly dependent.")

        u = scalar_mult(1 / norm, u)

        U.append(u)

    return U


# Example usage
V = [
    [1.0, 1.0, 0.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0]
]

U = gram_schmidt(V)

print("Orthonormal basis:")
for vec in U:
    print(vec)
