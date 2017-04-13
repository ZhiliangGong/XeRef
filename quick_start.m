
close all;
files = {'g20_30ps_1.ref', 'g20_30ps_mfeg8_1.ref', 'g30_30ps_mfeg8_1.ref', 'f30_30ps_mfge8_1.ref'};
this = XeRef(files);
this.control('load-pdb', 'MFGE8_NMR.pdb');
a = this.data{1};
s = this.layers;