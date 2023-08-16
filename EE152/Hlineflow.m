%Compute the reactive and real power loss on each harmonic frequency 
SLT = 0;
HLL = zeros(1,length(linedata(:,1)));
for n = 1:nbus
busprt = 0;
   for L = 1:nbr;
       if busprt == 0
       busprt = 1;
       else, end
       if nl(L)==n      k = nr(L);
       In = (V(n) - a(L)*V(k))*y(L)/a(L)^2 + Bc(L)/a(L)^2*V(n);
       Ik = (V(k) - V(n)/a(L))*y(L) + Bc(L)*V(k);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       HLL(1,n)=SL; 
       elseif nr(L)==n  k = nl(L);
       In = (V(n) - V(k)/a(L))*y(L) + Bc(L)*V(n);
       Ik = (V(k) - a(L)*V(n))*y(L)/a(L)^2 + Bc(L)/a(L)^2*V(k);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;                                              %<----- contains the harmonic line loss
                                                            
       else, end
       
  end
end
SLT = SLT/2;

