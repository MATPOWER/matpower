function [linsys_solver_string] = linsys_solver_to_string(linsys_solver)
% Convert linear systme solver integer to stringh
switch linsys_solver
    case osqp.constant('QDLDL_SOLVER')
        linsys_solver_string = 'qdldl';
    case osqp.constant('MKL_PARDISO_SOLVER')
        linsys_solver_string = 'mkl pardiso';
    otherwise
        error('Unrecognized linear system solver.');
end
