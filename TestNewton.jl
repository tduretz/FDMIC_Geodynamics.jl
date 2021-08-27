function f(x)
    return x^2 - 2
end

function f1(x)
    return 2x
end

method = "Newton"
method = "Steffensen"
niter  = 10
x      = 1.0
fv     = zeros(niter)
for iter = 1:niter
    println( f(x) ) 
    fv[iter] = abs( f(x) )
    if method == "Newton"
        x -= f(x)/f1(x)
    elseif method == "Steffensen"
        h  =  f(x)
        g  =  (f(x + h) - f(x)) / h
        x -= f(x)/g
    end
end

display( plot(1:niter, log10.(fv), markershape=:cross) )