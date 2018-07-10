function velocity_filtering(vmax, x, dx, tx)
kp = [2,2]';
kv = sqrt(kp);
lamb = kp ./ kv
x_tilde = x - tx
sat = vmax ./ (lamb .* abs(x_tilde))
scale = ones(2,1)
if any(sat < 1)
    [val, index] = min(sat)
    unclipped = kp .* x_tilde(index)
    clipped = kv .* vmax .* sign(x_tilde(index))
    scale = ones(2,1) .* clipped ./ unclipped
    scale(index) = 1
end
clipped = sat ./ scale
clipped(clipped > 1) = 1;
u_xyz = -kv .* (dx + clipped .* scale .* lamb .* x_tilde)
end