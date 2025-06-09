using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
    x = similar(b)
    ### votre code ici ; ne rien modifier d'autre
    # ...
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
    # ...
    # x = ...
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
