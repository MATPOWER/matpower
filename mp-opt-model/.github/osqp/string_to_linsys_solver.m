function [linsys_solver] = string_to_linsys_solver(linsys_solver_string)
linsys_solver_string = lower(linsys_solver_string);
switch linsys_solver_string
    case 'qdldl'
        linsys_solver = osqp.constant('QDLDL_SOLVER');
    case 'mkl pardiso'
        linsys_solver = osqp.constant('MKL_PARDISO_SOLVER');
    % Default solver: QDLDL
    case ''
        linsys_solver = osqp.constant('QDLDL_SOLVER');
    otherwise
        warning('Linear system solver not recognized. Using default solver QDLDL.')
        linsys_solver = osqp.constant('QDLDL_SOLVER');
end
