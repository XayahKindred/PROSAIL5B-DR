%********************************************************************************
%*                         function out = Jfunc3(k, l, t)                           
%*     
%*    Computes the J-function for given parameters k, l, and t. 
%*    Specifically, it calculates the integral form used in radiative transfer models.                                                                                                                                                               
%********************************************************************************
function out=Jfunc3(k,l,t)
out=(1.-exp(-(k+l)*t))/(k+l);
