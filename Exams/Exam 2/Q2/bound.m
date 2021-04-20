% bound.m
% Used to bound calls to sd. normal density/distribution function.
function y=bound(x);
    bval=6;
    bvec=abs(x)>bval;
    y=bvec.*(bval.*sign(x))+(1-bvec).*x;