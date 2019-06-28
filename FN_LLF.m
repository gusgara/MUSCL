function [fr,fl] = FN_LLF(flux,df,swi,sr,sl,M)

beta = max(abs(df(sr,swi,M)),abs(df(sl,swi,M)));
f = 0.5*(flux(sr,swi,M) + flux(sl,swi,M) - beta.*(sr -sl));
fr=f(2:end); fl=f(1:end-1);

end