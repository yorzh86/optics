function [abs, ems] = Spectral(wl, Abs)
nn = length(wl);
wl_0 = 0.4;
wl_1 = 2;
wl_2 = 2;
wl_3 = 20;

for ii = 1:nn-1
    abs_dum(ii) = (Abs(ii) + Abs(ii+1))/2*(wl(ii+1) - wl(ii));
end

abs = 0;
ems = 0;
for jj = 1:nn
    if (wl(jj)>wl_0)&&(wl(jj)<wl_1)
        abs = abs_dum(jj)+abs;
    elseif (wl(jj)>wl_2)&&(wl(jj)<wl_3)
        ems = abs_dum(jj)+ems;
    end
end

abs = abs/(wl_1-wl_0);
ems = ems/(wl_3-wl_2);
end