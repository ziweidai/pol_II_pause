function P = steady_P_before_TES(x,F,v0,k,b,TES)
d2 = (x-TES).^2;
P = (F/v0)*((d2+b)./(d2+k*b));
end