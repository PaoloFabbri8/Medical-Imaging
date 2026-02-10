% Sarvas Function
function B = sarvas_function(r, r0, Q, mu0)

a = r - r0;
a_n = norm(a);
r_n = norm(r);

F = a_n*(r_n*a_n + r_n^2 - dot(r0,r));

grad_F = r*(r_n^(-1)*a_n^2 + a_n^(-1)*dot(a,r) + 2*a_n + 2*r_n) - r0*(a_n + 2*r_n + a_n^(-1)*dot(a,r));

B = (mu0 / (4*pi*F^2) ) * (F*cross(Q,r0) - dot(cross(Q,r0), r)* grad_F);

end
