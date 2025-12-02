%build the mass and stiffness matrices that describe the 2nd order system.
%INPUTS
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
%OUTPUTS
%M_mat: the n x n mass (inertia) matrix
%K_mat: the n x n stiffness matrix
function [M_mat,K_mat, ec] = construct_2nd_order_matrices(string_params)
    n = string_params.n;
    M = string_params.M;
    dx = string_params.dx;
    Tf = string_params.Tf;
    % —— Persistent cached matrices ——————————————
    persistent cached_n tri_diag Mm Km ecv

    % Recompute only if not defined or size changed
    if isempty(cached_n) || cached_n ~= n
        cached_n = n;

        % Build tri-diagonal
        tri_diag = -2*eye(n);
        tri_diag(2:n,1:n-1) = tri_diag(2:n,1:n-1) + eye(n-1);
        tri_diag(1:n-1,2:n) = tri_diag(1:n-1,2:n) + eye(n-1);

        % Mass matrix
        Mm = (M/n) * eye(n);

        Km = -Tf/dx*tri_diag;

        % End condition vector
        cond_end = zeros(n,1);
        cond_end(end) = 1;
        ecv = Tf/dx*cond_end;
    end
    M_mat = Mm;
    K_mat = Km;
    ec = ecv;
    % ————————————————————————————————
end