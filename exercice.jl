using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
    x = similar(b)

    n = size(R, 1)
    for i in n:-1:1
        c_i = b[i]

        if i < n
            for j in i+1:n
                c_i -= R[i, j]*x[j]
            end
        end

        x[i] = c_i / R[i, i]
    end

    return x
end

# 2. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système Hessenberg supérieur Hx = b ou du problème aux moindres
#    carrés min ‖Hx - b‖ à l'aide de rotations ou
#    réflexions de Givens et d'une remontée triangulaire.
#    Votre fonction peut modifier H et b si nécessaire.
#    Il n'est pas nécessaire de garder les rotations en mémoire et la
#    fonction ne doit pas les renvoyer.
#    Seul le cas réel sera testé ; pas le cas complexe.
function hessenberg_solve(H::UpperHessenberg, b)
    # Transformation de la matrice de Hessenberg à une matrice triangulaire supérieure
    R = copy(H)
    m, n = size(H)

    # Rotations de Givens
    for i in 1:n-1
        Q_i_j = Matrix{Float64}(I, m, m)
        j = i+1

        x, y = R[i, i], R[j, i]
        rho = norm([x, y])

        c, s = x/rho, y/rho

        Q_i_j[i, i] = c
        Q_i_j[j, j] = c
        Q_i_j[i, j] = s
        Q_i_j[j, i] = -s

        R = Q_i_j * R
        b = Q_i_j * b
    end
    
    # Réflexions de Givens
    if m > n
        e1 = zeros(m-n+1)
        e1[1] = 1
        v = copy(R[n:m, n])
        
        u = v - sign(v[1])*norm(v)*e1
        u_ = u/(u'*u)
        for j in 1:n
            R[n:m, j] .-= 2*u_*(u'*R[n:m, j])
        end
        b[n:m] .-= 2*u_*(u'*b[n:m])

    end

    # Remontée triangulaire
    R_thin = UpperTriangular(R[1:n, :])
    b_thin = b[1:n]
    x = backsolve(R_thin, b_thin)

    return x
end

# vérification
using Test
for n ∈ (10, 20, 30)
    # square system
    A = rand(n, n)
    A[diagind(A)] .+= 1
    b = rand(n)
    R = UpperTriangular(A)
    x = backsolve(R, b)
    @test norm(R * x - b) ≤ sqrt(eps()) * norm(b)
    H = UpperHessenberg(A)
    x = hessenberg_solve(copy(H), copy(b))
    @test norm(H * x - b) ≤ sqrt(eps()) * norm(b)
    # slightly overdetermined least squares
    A = rand(n + 1, n)
    A[diagind(A)] .+= 1
    b = rand(n+1)
    H = UpperHessenberg(A)
    x_ls = hessenberg_solve(copy(H), copy(b))
    x_qr = H \ b
    @test norm(x_ls - x_qr) ≤ sqrt(eps()) * norm(x_qr)
end
