function [IM] = InterpolationMatrix(N,r,V)
 Vr = Vandermonde1D(N, r);
 IM = Vr/V;
 return
end

