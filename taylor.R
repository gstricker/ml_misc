## Taylor expansion
###################

f <- function(x) {
    1e-6*x^4 - 1e-3*x^3 + x
}

f_prime <- function(x) {
    4e-6*x^3 - 3e-3*x^2 + 1
}

f_2prime <- function(x) {
    12e-6*x^2 - 6e-3*x
}

f_3prime <- function(x) {
    24e-6*x - 6e-3
}

taylor <- function(x, a) {
    f(a) + f_prime(a)*(x - a)
}

taylor2 <- function(x, a) {
    taylor(x, a) + (f_2prime(a)*(x-a)^2)/2
}

taylor3 <- function(x,a) {
    taylor2(x, a) + (f_3prime(a)*(x-a)^3)/factorial(3)
}

x <- -100:100

plot(x, f(x), type = "l")
lines(x, taylor(x, 50), col = "red")
lines(x, taylor2(x, 50), col = "blue")
lines(x, taylor3(x, 50), col = "green")
