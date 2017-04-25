theta = this.layers.fits.angles.para_range_1;
phi = this.layers.fits.angles.para_range_2;
lk = this.layers.fits.angles.likelihood;

[Phi, Theta] = meshgrid(phi, theta);
figure;
surf(Phi, Theta, lk);

figure;
surf(Phi, Theta, this.layers.fits.angles.chi2);

figure;
contourf(Phi, Theta, this.layers.fits.angles.chi2);