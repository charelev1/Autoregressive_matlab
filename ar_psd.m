function Sxx = ar_psd(sigma, a, p, w)
    a = a(:);
	Sxx = zeros(length(w),1);
    for i = 1: length(w)
        ee = exp(-1i*(0:p)*w(i));
        dem = ee*a;
        Sxx(i) = 1/ (dem*dem');
    end
    Sxx = abs(Sxx);
    Sxx = Sxx*sigma;
end

