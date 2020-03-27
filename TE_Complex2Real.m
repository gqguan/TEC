function [xout] = TE_Complex2Real(xin, imag_err)
%% get the positive real part of "xin" if imag(xin)<imag_err
%  notes of I/O arguments
%  xin      - (i complex) input complex number
%  imag_err - (i real scalar) error of imaginary part
%  xout     - (o real) real part of "xin"
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-07
% 
%%
% get the array where each element's real part is greater than 0
x = xin(real(xin)>0);
% get the real part of the arry where each element's imaginary part is less
% than imag_err
xout = double(real(x(imag(x)<imag_err)));
%
end