function pdb = processPdbFile(text)

raw = textscan(text, '%s %d %s %*s %*s %*d %f %f %f %*f %*f %*s', 'HeaderLines', 1);

labels = raw{1};
anisous = raw{2};

selects = (strcmpi('atom', labels) | strcmpi('hetatm', labels)) & anisous;
pdb.x = raw{4}(selects)';
pdb.y = raw{5}(selects)';
pdb.z = raw{6}(selects)';

atomNames = raw{3}(selects)';
atoms = repmat(' ', 1, length(atomNames));

for i = 1 : length(atoms)
    atoms(i) = atomNames{i}(1);
end

pdb.atoms = atoms;

end