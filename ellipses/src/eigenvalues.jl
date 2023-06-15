function calc_eigvals_2x2(A::Matrix; ϵ=1e-6)
    @assert size(A) == (2,2) "size(A)=$(size(A)), specialised function only works on 2x2 matrices"
    a = A[1, 1]
    b = A[1, 2]
    c = A[2, 1]
    d = A[2, 2]
    
    if abs(c) < ϵ && abs(b) < ϵ
        λs = [a, d]
        W = [[1.0 0.0]; [0.0 1.0]]
        idxs = sortperm(λs)
        return λs[idxs], W[:, idxs]
    end
        
    Δ = 0.5 * sqrt((a - d) * (a - d) + 4 * b * c) 
    λ1 = (a + d) / 2 - Δ
    λ2 = (a + d) / 2 + Δ
    
    if abs(c) < ϵ
        W1 = [b; (λ1 - a)]
        W2 = [b; (λ2 - a)]
    else
        W1 = [(λ1 - d) ; c]
        W2 = [(λ2 - d) ; c]   
    end
    W1 = W1 / sqrt(sum(W1 .^ 2))
    W2 = W2 / sqrt(sum(W2 .^ 2))
    W = hcat(W1, W2) 
    
    [λ1, λ2], W
end