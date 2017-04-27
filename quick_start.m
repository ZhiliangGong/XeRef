
close all;
files = {'g20_30ps_1.ref', 'g20_30ps_mfeg8_1.ref', 'g20_30ps_mfeg8_2.ref', 'g20_30ps_mfeg8_3.ref', 'g25_30ps_mfeg8_1.ref', 'g30_30ps_mfeg8_1.ref', 'f30_30ps_mfge8_1.ref'};
this = XeRef(files);
% pdbfile = 'mfge8_trial1_protonly.pdb';
% pdbfile = 'hmmm-end-frame-protien.pdb';
% pdbfile = 'mfge8_trial5bfearlyf199.pdb';
% pdbfile = 'mfge8_trial5bflatef199.pdb';
% pdbfile = 'mfge8-nmr-histag.pdb';
% pdbfile = 'protein-lowest-rmsd-frame-239.pdb';
% pdbfile = 'protein-frame-2000-aligned.pdb';
pdbfile = 'mfge8-nmr-histag-aligned.pdb';
this.control('load-pdb', pdbfile);
a = this.data{1};
s = this.layers;

% this.control('toggle-protein');
% this.gui.parametersTable.Data(:, 3) = parameters_0425(:, 3);

% this.layers.protein.contourHeight();