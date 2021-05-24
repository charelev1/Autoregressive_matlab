function [a,sigma] = myar_methods(x, p, str )
    N = length(x);
    x =x(:);
    if  strcmp(str,'Burg')
        
        a = zeros(p+1,1);
        a(1) = 1;
        rho = x'*x./N;
        f = x(2:end);
        b = x(1:end-1);
        
        for i = 1: p 
            k = -(2.*b'*f)/(f'*f + b'*b);
            rho = (1-k*k')*rho;
            ftmp = f(2:end)    + k  .*b(2:end);
            b = b(1:end-1) + k' .* f(1:end-1);
            f = ftmp;	
            a(2:i+1) = a(2:i+1) + k .* conj(a(i:-1:1));	
        end 
        sigma = rho;
    else
        
        rxx = zeros(p,1);
        xx = zeros(1,2*N-1);
        xx(1:N) = x;
        yy = xx';

        for i = 1: p+1
            rxx(i) = xx(1:end-p-1)*yy(1:end-p-1);
            yy = circshift(yy,1);
        end
        rxx = rxx/N;
        R = toeplitz(rxx(1:end-1));
        a = R\-rxx(2:end);
        sigma = [1; a]' * rxx;
        a = [1;a];
    end
end

