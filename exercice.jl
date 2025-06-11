using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
    x = similar(b)
    ### votre code ici ; ne rien modifier d'autre
    n = length(b)
    x[n] = b[n] / R[n,n]
    for i in n-1:-1:1
        x[i] = (b[i] - sum(R[i,j]*x[j] for j in i+1:n; init=zero(eltype(b)))) / R[i,i]
    end
    ###
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
    ### votre code ici ; ne rien modifier d'autre
    m, n = size(H)

    for i in 1:m-1
        h₁, h₂ = H[i,i], H[i+1,i]
        ρ = √(h₁^2+h₂^2)
        c, s = h₁/ρ, h₂/ρ

        for j in i:n
            Hij, Hip1j = H[i,j], H[i+1,j]
            H[i,j] = c*Hij + s*Hip1j
            H[i+1,j] = c*Hip1j - s*Hij
        end

        bi, bip1 = b[i], b[i+1]
        b[i] = c*bi + s*bip1
        b[i+1] = c*bip1 - s*bi
    end

    x = backsolve(UpperTriangular(H[1:n,:]),b[1:n])
    ###
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
    H = UpperHessenberg(A)
    b = rand(n + 1)
    x_ls = hessenberg_solve(copy(H), copy(b))
    x_qr = H \ b
    @test norm(x_ls - x_qr) ≤ sqrt(eps()) * norm(x_qr)
end
