% acvf()
%
% return autocovariance function at lag(s) x given a variance, a
% correlation length, and a function type ('E' for exponential, 'G' for
% Gaussian)

function acv=acvf(x,var,corrlength,acvftype)

switch acvftype
    case {'E','e'}  % exponential covariance
        exparg=-(1/corrlength)*abs(x);
    case {'G','g'}  % gaussian bell covariance
        exparg=-(1/(2*corrlength^2))*x.^2;
    otherwise
        error('acvf: priortype %s unknown',acvftype);
end
acv=var*exp(exparg);