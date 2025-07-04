function obj = struct2object(s)
% mp.struct2object - Convert a struct back to the object from which it was created.
% ::
%
%   obj = mp.struct2object(s)
%
% Input:
%   s (struct) : a struct of the form created by a :meth:`to_struct` method,
%       containing the data from the object plus a char array field naned
%       ``class_`` with the class name of the desired object, and an optional
%       cell array field named ``constructor_args_`` with arguments to pass
%       to the constructor.
%
% Output:
%   obj (classdef object) : an instance of the object identical to the one
%       used to create the input struct
%
% As of version 10.x, Octave is still not able to save and load classdef
% objects. To aid in creating workarounds, this function allows objects to
% implement the following pattern with appropriately coded :meth:`to_struct`
% and :meth:`from_struct` methods::
%
%   s = obj.to_struct();            % convert to normal struct
%   new_obj = mp.struct2object(s);  % convert back to classdef object
%   isequal(new_obj, obj)           % returns true
%
% The :meth:`to_struct` method of the object must create a struct containing
% all of the data in the object, plus a char array field naned ``class_``
% with the class name of the desired object, and an optional cell array field
% named ``constructor_args_`` with arguments to pass to the object constructor.
%
% The :meth:`from_struct` method takes a freshly constructed object and the
% struct above and copies the data from the struct back to the object.
%
% This function creates an instance of the specified class, by calling
% its constructor with any specified arguments, then calling its
% :meth:`from_struct` method.

%   MP-Opt-Model
%   Copyright (c) 2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if isfield(s, 'class_')
    try
        if isfield(s, 'constructor_args_')
            obj = feval(s.class_, s.constructor_args_{:});
            s = rmfield(s, {'class_', 'constructor_args_'});
        else
            obj = feval(s.class_);
            s = rmfield(s, 'class_');
        end
    catch me
        error('mp.struct2object : unable to create object of class "%s" : %s', s.class_, me.message);
    end
else
    error('mp.struct2object : unable to create object of class "%s"', s.class_);
end

obj.from_struct(s);
