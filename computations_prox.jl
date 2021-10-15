# Computation prox of the 1D functions

function proxCosh(x, gamma =1., tol = 1e-10)

    # We want to solve
    # y + gamma*sinh(y) = x

    if x > 0
        # Newton's iteration

        # Counter to check it does not go too far
        counter = 0

        # Do an initialization depending if x >> 1 or not
        y = min(x/(1+gamma), asinh(x/(gamma+tol)))


        f = y + gamma * sinh(y) - x

        while abs(f) > tol && counter <= 10
            df = 1 + gamma * cosh(y)
            y -= f/df
            f = y + gamma * sinh(y) - x

            counter += 1
        end

        if counter == 11
            println("Problem computation prox coshh, too many iterations")
        end

        return y

    elseif x < 0
        return -proxCosh(-x,gamma,tol)
    else
        return 0.
    end
end
