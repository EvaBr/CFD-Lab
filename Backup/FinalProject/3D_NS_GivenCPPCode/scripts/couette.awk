#!@gawk@ -f
# Get an analytical slution for Couette flow
# Input parameters:
# V0: velocity of the wall
# nu: kinimatic viscosity [L^2 T^-1]
# L: the length of the domain
# t: the time 
# provide set of y as an input stream 
# the wall with y=L is moving and has velocity V0
# 
# Output:
# the the two columns
# y v_y

function f(y, n) {
    return 2*V0 / (n*pi) * (-1)^n * sin(n*pi/L*y) * exp (-nu*n^2*pi^2/L^2*t)
}

function abs(x) {
    if (x>0) {
	return x
    } 
    else {
	return -x
    }
}

BEGIN {
    if (!V0) {V0=1.25e-5}
    if (!nu) {nu=1.0e-6}
    if (!L) {L=1e-3}
    if (!t) {t=0.1}

    pi=3.1416
    eps=1e-18
}


{
    y = $1
    $1 = "" 
    rest = $0

    if (y>L) {
	print "y cannot be bigger then L" 
	exit(-1)
    }
    fn = 1e10
    fsum = 0.0
    n=1
    if (t==0) {
	fsum=-V0/L*y
    }
    else {
	while (abs(fn) > eps) {
	    fn = f(y, n)
	    fsum+=fn
	    n++
	}
    }
    
    print y, rest, V0/L*y + fsum
}