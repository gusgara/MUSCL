function [fr,fl] = MUSCL(s0,f,df,Swi,M,dx)
           
            N=length(s0);
            s0=[s0(1) s0 s0(end)];
            sl=zeros(1,N+1); sr=zeros(1,N+1);
            
            %extrapolação
            
            rl =(s0(3:N+2)-s0(2:N+1) +dx^2)./(s0(2:N+1)-s0(1:N)+dx^2);
            sl(2:end)=s0(2:N+1) + (1/6).*VanAlbada(rl).*(2.*s0(3:N+2)-s0(2:N+1)-s0(1:N));
            
            rr = (s0(2:N+1)-s0(1:N)+dx^2)./(s0(3:N+2)-s0(2:N+1)+dx^2);
            sr(1:N) = s0(2:N+1) - (1/6).*VanAlbada(rr).*(s0(3:N+2)+s0(2:N+1)-2*s0(1:N));

   
    %% Fluxo númerico
    Sw=s0(2:N+1);
    sl(1) =0.9; sr(end) = sr(end-1); % condição de contorno
    [fr,fl] = FN_LLF(f,df,Swi,sr,sl,M);
end
