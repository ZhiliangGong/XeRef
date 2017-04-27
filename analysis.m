traj = readdcd('protein.dcd');
rmsd = superimpose(traj(1, :), traj);

n_frame = size(traj, 1);

figure;
plot(rmsd, 'o');

t_seconds = 1 : n_frame;
f_indinces = 1 : n_frame;

[F, G] = meshgrid(f_indinces, f_indinces);

rmsd_mat = zeros(n_frame, n_frame);
parfor i = 1 : n_frame
    rmsd_mat(:, i) = superimpose(traj(i, :), traj);
end

rmsd_sum = sum(rmsd_mat, 2);
figure;
plot(t_seconds, rmsd_sum, 'o');

