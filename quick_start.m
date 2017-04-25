
close all;
files = {'g20_30ps_1.ref', 'g20_30ps_mfeg8_1.ref', 'g25_30ps_mfeg8_1.ref', 'g30_30ps_mfeg8_1.ref', 'f30_30ps_mfge8_1.ref'};
this = XeRef(files);
% pdbfile = 'mfge8_trial1_protonly.pdb';
% pdbfile = 'hmmm-end-frame-protien.pdb';
pdbfile = 'mfge8_trial5bfearlyf199.pdb';
this.control('load-pdb', pdbfile);
a = this.data{1};
s = this.layers;
% this.layers.protein.contourHeight();