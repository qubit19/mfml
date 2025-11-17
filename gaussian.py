def gaussian_elimination_all_cases(A, b):
    n = len(A)
    m = len(A[0])   # number of variables
    
    # Form augmented matrix [A | b]
    aug = [A[i] + [b[i]] for i in range(n)]

    row = 0
    for col in range(m):
        # Find pivot row
        pivot = max(range(row, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < 1e-12:
            continue     # no pivot in this column → free variable

        # Swap
        aug[row], aug[pivot] = aug[pivot], aug[row]

        # Normalize pivot row
        pivot_val = aug[row][col]
        aug[row] = [x / pivot_val for x in aug[row]]

        # Eliminate all other rows
        for r in range(n):
            if r != row:
                factor = aug[r][col]
                aug[r] = [aug[r][c] - factor * aug[row][c] for c in range(m+1)]
        
        row += 1
        if row == n:
            break

    # After elimination, check special cases
    # 1. Check for inconsistency (0 ... 0 | nonzero)
    for r in range(n):
        if all(abs(aug[r][c]) < 1e-12 for c in range(m)) and abs(aug[r][m]) > 1e-12:
            return "No solution (inconsistent system)"

    # 2. Determine rank
    rank = sum(any(abs(aug[r][c]) > 1e-12 for c in range(m)) for r in range(n))

    if rank == m:
        # Unique solution
        x = [0]*m
        for r in range(m):
            x[r] = aug[r][m]
        return x
    
    # 3. Infinite solutions → parametric form
    free_vars = []
    pivot_cols = set()

    # Identify pivot columns
    for r in range(n):
        for c in range(m):
            if abs(aug[r][c]) > 1e-12:
                pivot_cols.add(c)
                break

    free_vars = [c for c in range(m) if c not in pivot_cols]

    # Build parametric solution
    solution = {f"x{c}": {} for c in range(m)}

    # Express pivot variables
    for r in range(n):
        leading = None
        for c in range(m):
            if abs(aug[r][c]) > 1e-12:
                leading = c
                break
        if leading is None:
            continue

        sol = aug[r][m]
        expr = {"constant": sol}

        # Terms for free variables
        for fv in free_vars:
            expr[f"t{fv}"] = -aug[r][fv]

        solution[f"x{leading}"] = expr

    # Free variables = parameters
    for fv in free_vars:
        solution[f"x{fv}"] = {f"t{fv}": 1, "constant": 0}

    return {"Infinite solutions": solution}
