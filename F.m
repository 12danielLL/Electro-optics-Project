function F_x=F(x)
F_x=fftshift(fft2(ifftshift(x)));
end