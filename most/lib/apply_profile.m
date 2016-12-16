function argout = apply_profile(profile, argin, dim)
%APPLY_PROFILE  Applies changes defined in a profile to a data structure.
%
%   CHGTABS = APPLY_PROFILE( PROFILE, CHGTABSI )
%   XGD     = APPLY_PROFILE( PROFILE, XGDI, DIM )
%   SD      = APPLY_PROFILE( PROFILE, SDI, DIM )
%   CTSETS  = APPLY_PROFILE( PROFILE, CTSETS, DIM ) (not yet implemented)
%   
%   Applies a single profile of the given type to the given ARGIN. There are
%   4 different types of profiles, and each one affects differently the
%   input to produce the output. Profile input must contain the following
%   fields:
%
%   Inputs:
%       PROFILE: a single-dimensional Profile struct with the following fields:
%           .type       (string)
%           .table      (string or scalar)
%           .rows       (vector)
%           .col        (scalar)
%           .chgtype    (scalar)
%           .values     (array) with at most 3 dimensions assumed to be
%                         [ (1 or nt) by (1 or nj_max) or (1 or length(rows)) ]
%           See IDX_PROFILE for details on the Profile struct.
%
%       CHGTABI:    cell array of change tables to be modified
%       XGDI:       xGenData struct to be modified
%       STORAGEI:   StorageData struct to be modified
%       CTSETSI:    array with ContingencyData (not yet implemented)
%
%       DIM :   (scalar) indicates the total number of elements in the
%               table or field being modified. Elements here refers to the
%               3rd dimension, not time nor scenarios, but rather elements
%               such as generators, number of different contingencies in
%               master chgtab matrix (different labels), and storage units.
%               DIM required to be able to expand, to a full DIM dimension,
%               the data to be modified when it is summarized by a
%               singleton dimension representing all the elements in that
%               particular data set. It is ignored for type 'mpcData',
%               mandatory for all other types.
%
%   Outputs:
%       CHGTABS : cell array of modified change tables (nr x 7)
%       XGD:      modified xGenData struct
%       STORAGE:  modified StorageData struct
%       CTSETS:   (not yet implemented)
%
%   Additional notes:
%
%   In general, field 'values' does not need to match dimensions of
%   dim = [nt nj_max n], where 'n' represents the subset of elements being
%   affected by the profile, i.e., the elements indicated by 'rows', but it
%   does need to be smaller or equal. Each dimension of 'values' is allowed
%   to be either the indicated above, or a singleton dimension, in which
%   case a singleton meaning that the profile "applies to all" elements in
%   that dimension, with the exception that the third dimension may be a
%   singleton also in the case affecting a single element (as opposed to all
%   elements in the third dimension).
%
%   type == 'mpcData'
%       Generates/adds contingency-like tables to a cell array can be
%       used to apply a change 'chgtype' to values in column 'col' of
%       elements 'rows' on table 'table'. 'values' is a numeric array
%       with up to 3 dimensions organized necessarily as in [nt nj_max
%       n]. The third dimension indicates the subset of elements to
%       which the profile is to be applied. Output CHGTABS is a (nt by
%       nj) cell array of chgtab matrices (7 cols) with unspecified
%       labels nor probabilities. CHGTABI must always be provided, even
%       if it's a cell array with (nt x nj_max) empty entries. These
%       dimensions are required in order to be able to expand changes
%       correctly across time periods and scenarios. Dimensions of
%       'values' are expanded if required (i.e., if inconsistent with nt,
%       nj_max, or length of 'rows', resp.).
%
%   type == 'xGenData'
%       Profile modifies the field of XGD struct indicated by the
%       string 'table'. 'rows' indicates gens to modify, 'col' is
%       ignored, and 'chgtype' the type of change using 'values'.
%       Dimensions of 'values' are expanded if required (ie, if
%       inconsistent with nt, nj_max, or length of 'rows', resp.).
%
%   type == 'StorageData'
%       Profile modifies the field of STORAGE struct indicated by the
%       string 'table'. 'rows' indicates storage units to modify (using
%       storage-dedicated idx's as opposed to gen idx's), 'col' is ignored,
%       and 'chgtype' the type of change using 'values'. Dimensions of
%       'values' are expanded if required (ie, if inconsistent with nt,
%       nj_max, or length of 'rows', resp.).
%
%   type == 'ContingencyData' (not yet implemented)
%       Profile modifies the provided 'indicative' 3-dim array of
%       binary variables: 1st dim spans the labels of contingencies,
%       2nd dim spans time periods, and 3rd dim spans scenarios. Thus,
%       'rows' indicates which labeled contingencies are to be modified
%       by the profile, 'col' is ignored, and 'cghtype' the type of
%       change using 'values'. Dimensions of 'values' expanded if
%       required  (ie, if inconsistent with nt, nj_max, or length of
%       'rows', resp.).

% Created by Daniel Munoz-Alvarez (4/18/2013)

%   MOST
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Daniel Munoz-Alvarez, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 2
    error('apply_profile: insufficient arguments')
end

% (A) Preliminary checking

typ     = profile.type;
tbl     = profile.table;
rows    = profile.rows;
col     = profile.col;
chgtyp  = profile.chgtype;
val     = profile.values;

[PR_REP, PR_REL, PR_ADD, PR_TCONT, PR_TYPES, PR_TMPCD,...
    PR_TXGD, PR_TCTD, PR_TSTGD, PR_CHGTYPES] = idx_profile;

if length(profile)~=1
    error('apply_profile: must input a single profile')
end
if length(rows) > 1 && any(rows == 0)
    error('apply_profile: rows field of a profile must not contain zero unless it is the only entry')
end
    

switch typ
% (B) Type mpcData profile
    case 'mpcData'
        chgtabs = argin;

        nt = size(chgtabs, 1);
        nj_max = size(chgtabs, 2);

% (B.1) Check dimensions and fields of PROFILE
        if length(profile) > 1
            error('apply_profile: multiple profiles should be added separately')
        end

        % (B) Check consistency of IDX, VALUES and CHGTABI
        if size(val,3) ~= length(rows)
            error('apply_profile: third dimension of profile.values should match length of profile.rows')
        end

        if isempty(chgtabs)
            error('apply_profile: chgtabs cell array should have dimensions nt by nj_max')
        end


        % (C) Generate contingency-like rows to add to CHGTABI
        
%       At this point, val can only have dimensions (1 or nt) by
%       (1 or nj_max) by (1 or length(rows)), so before transforming into a
%       chgtab, val needs to be expanded to full dimensions [nt nj_max length(rows)]
        
        if size(val,1) == 1 && nt > 1
            val = repmat(val, [nt 1 1]);
        end
        if size(val,2) == 1 && nj_max > 1
            val = repmat(val, [1 nj_max 1]);
        end
        if size(val,3) == 1 && length(rows) > 1
            val = repmat(val, [1 1 length(rows)]);
        end
        
        if any(tbl == PR_TMPCD)
            for t = 1:nt
                for j = 1:nj_max
                    new_rows = [];
                    for i = 1:length(rows)
                            new_rows = [ new_rows ; 0 0 tbl rows(i) col chgtyp val(t,j,i) ];
                    end
                    chgtabs{t,j} = [ chgtabs{t,j} ; new_rows ];
                end
            end
        else
            error('apply_profile: indicated profile.table not supported for profile changes')
        end

        argout = chgtabs;
    
    case 'xGenData'
% (C) Type xGenData profile
        xgd = argin;
        ng = dim;
        
        nt_adhoc = max( size(val,1), size(xgd.(tbl),2) ); % nt_adhoc equals 1 or nt
        nj_adhoc = max( size(val,2), size(xgd.(tbl),3) ); % nj_adhoc equals 1 or nj_max

% (C.1) Check fields of PROFILE
        if ~ischar(tbl)
            error('apply_profile: table field of xGenData profile must be a string')
        end
        if ~any(strcmp(tbl, PR_TXGD))
            error('apply_profile: field %s of xgd struct cannot be modified through a profile',tbl)
        end

% (C.2) Check consistency of ROWS, VALUES and affected field of XGD

        if ~(size(xgd.(tbl),1) == 1 || size(xgd.(tbl),1) == ng)
             error('apply_profile: rows of xgd.%s must equal 1 or ng',tbl)
        end
        if ~(size(xgd.(tbl),2) == 1 || size(xgd.(tbl),2) == nt_adhoc)
            error('apply_profile: time dimension mismatch between xgd.%s and profile.values',tbl)
        end
        if ~(size(xgd.(tbl),3) == 1 || size(xgd.(tbl),3) == nj_adhoc)
            error('apply_profile: scenario dimension mismatch between xgd.%s and profile.values',tbl)
        end
        
        
        if ~(size(val,1) == 1 || size(val,1) == nt_adhoc)
            error('apply_profile: time dimension mismatch between xgd.%s and profile.values',tbl)
        end
        if ~(size(val,2) == 1 || size(val,2) == nj_adhoc)
            error('apply_profile: scenario dimension mismatch between xgd.%s and profile.values',tbl)
        end
        if length(rows) == 1
            if size(val,3) ~= 1
                error('apply_profile: 3rd dimension of values field must equal 1 when rows is scalar')
            end
        else
            if ~(size(val,3) == length(rows) || size(val,3) == 1)
                error('apply_profile: 3er dimension of values field must equal 1 or length of rows field when rows is not scalar')
            end
        end

        
% (C.3) Verify validity of changes and expand field in question to full dimensions if required
%       Important: Notice the permutation of values dimensions from 
%       [1 2 3] to [3 1 2].

          val = permute(val,[3 1 2]);
          
          % From here and on, val dimensions are
          % (1 or legnth(rows)) by (1 or nt) by (1 or nj)
          % Also, xgd.(tbl) dimensions are
          % (1 or ng) by (1 or nt) by (1 or nj)


          % Expand cols (time) of xgd.(tbl)
          % if val has a time dimension and xgd.(tbl) does not (only necessary if nt_adhoc > 1)
          if size(val,2) > size(xgd.(tbl),2)
              xgd.(tbl) = repmat( xgd.(tbl), [1 size(val,2) 1]);
          
          % Expands cols (time) of val to match xgd.(tbl) time dimension
          elseif size(xgd.(tbl),2) > size(val,2)
              val = repmat( val, [ 1 nt_adhoc 1]);
          end
          
          % Expand 3rd dim (scenarios) of xgd.(tbl)
          % if val has a scenario dimension and xgd.(tbl) does not (only necessary if nj_adhoc > 1)
          if size(val,3) > size(xgd.(tbl),3)
              xgd.(tbl) = repmat( xgd.(tbl), [ 1 1 size(val,3)] ) ;
          
          % Expands 3rd dim (scenarios) of val to match xgd.(tbl) scenarios dimension
          elseif size(xgd.(tbl),3) > size(val,3)
              val = repmat(val, [ 1 1 size(xgd.(tbl),3)]);
          end

          % Expand rows (gens) of xgd.(tbl)
          % if profile modifies subset of gens (i.e. if rows~=0)
          if size(xgd.(tbl),1) == 1 && any(rows ~= 0)
              xgd.(tbl) = repmat( xgd.(tbl), [ ng 1 1]);
          end

          % Expand rows of val to match gens to modify
          if size(val,1) == 1 && size(xgd.(tbl), 1) == ng && length(rows) > 1
              val = repmat(val, [length(rows) 1 1]);
          elseif size(val,1) == 1 && size(xgd.(tbl), 1) == ng && length(rows) == 1 && rows == 0
              val = repmat(val, [ng 1 1]);
          end

          % Error if chgtype is CT_REL or CT_ADD and field involved is empty
          if (chgtyp == PR_REL || chgtyp == PR_ADD) && isempty(xgd.(tbl))
              error('apply_profile: PR_REL or PR_ADD modification cannot be done to xgd.%s if it is empty',tbl)
          end

% (C.4) Apply change
          if length(rows) == 1 && rows == 0                 %% modify all rows
            if chgtyp == PR_REP                                    %% replace
              xgd.(tbl) = val;
            elseif chgtyp == PR_REL                                %% scale
              xgd.(tbl) = val .* xgd.(tbl);
            elseif chgtyp == PR_ADD                                %% shift
              xgd.(tbl) = val + xgd.(tbl);
            else
              error('apply_profile: modification type %d for xgd table not supported', chgtyp);
            end
          else                                              %% modify single row
            if chgtyp == PR_REP                                    %% replace
              xgd.(tbl)(rows,:,:) = val;
            elseif chgtyp == PR_REL                                %% scale
              xgd.(tbl)(rows,:,:) = val .* xgd.(tbl)(rows,:,:);
            elseif chgtyp == PR_ADD                                %% shift
              xgd.(tbl)(rows,:,:) = val + xgd.(tbl)(rows,:,:);
            else
              error('apply_profile: modification type %d for xgd table not supported', chgtyp);
            end
          end

        
        argout = xgd;

    case 'ContingencyData'
% (D) Type ContingencyData profile (not yet implemented)
        ct_subset = argin;
        error('apply_profile: type ContingencyData is not yet supported')
        argout = ct_subset;

    case 'StorageData'
% (E) Type StorageData profile
        storage = argin;
        ns = dim;   % total number of storage units in the system
        nt_adhoc = max( size(val,1), size(storage.(tbl),2) ); % nt_adhoc equals 1 or nt

% (E.1) Check fields of PROFILE
        if ~ischar(tbl)
            error('apply_profile: table field of storage profile must be a string')
        end
        if ~any(strcmp(tbl, PR_TSTGD))
            error('apply_profile: field %s of storage struct cannot be modified through a profile',tbl)
        end

% (E.2) Check consistency of ROWS, VALUES and affected field of STORAGE
        
        if ~(size(storage.(tbl),1) == 1 || size(storage.(tbl),1) == ns)
             error('apply_profile: first dimension of field %s of storage struct must equal 1 or ns',tbl)
        end
        if ~(size(storage.(tbl),2) == 1 || size(storage.(tbl),2) == nt_adhoc)
            error('apply_profile: time dimension mismatch between storage.%s and profile.values',tbl)
        end
        if size(storage.(tbl),3) ~= 1
            error('apply_profile: no scenario dimension (3rd) allowed for storage.%s field',tbl)
        end
        
        
        if ~(size(val,1) == 1 || size(val,1) == nt_adhoc)
            error('apply_profile: time dimension mismatch between storage.%s and profile.values',tbl)
        end
        if size(val,2) ~= 1
            error('apply_profile: 2nd dimension of values field must equal 1 since no scenario dependent changes allowed ')
        end
        if length(rows) == 1
            if size(val,3) ~= 1
                error('apply_profile: 3rd dimension of values field must equal 1 when rows is scalar')
            end
        else
            if ~(size(val,3) == length(rows) || size(val,3) == 1)
                error('apply_profile: 3er dimension of values field must equal 1 or length of rows field when rows is a vector')
            end
        end

        
% (E.3) Verify validity of changes and expand field in question to full dimensions if required
%       Important: Notice the permutation of values dimensions from 
%       [1 2 3] to [3 1 2]. No scenario dimension allowed.

          val = permute(val,[3 1 2]); % Squeeze not use to avoid problems when nt=1
          val = val(:,:,1); % From here and on, val dimensions are (1 or legnth(rows) by 1 or nt)


          % Expand cols (time) of field involved
          if size(val,2) > size(storage.(tbl),2)
              storage.(tbl) = storage.(tbl) * ones(1, size(val,2));
          end

          % Expand rows (ess units) if change modifies submatrix of the parameter that is being modified (not all rows)
          if size(storage.(tbl), 1) == 1 && any(rows ~= 0)
              storage.(tbl) = ones(ns, 1) * storage.(tbl);
          end

          % Expand rows of val to match ess units to change
          if size(val,1) == 1 && size(storage.(tbl), 1) == ns && length(rows) > 1
              val = ones(length(rows), 1) * val;
          end
          if size(val,1) == 1 && size(storage.(tbl), 1) == ns && length(rows) == 1 && rows == 0
              val = ones(ns, 1) * val;
          end

          % Expands cols of val to match field's time dimension
          if size(val, 2) == 1 && size(storage.(tbl), 2) == nt_adhoc
              val = val * ones(1, nt_adhoc);
          end

          % Error if chgtype is CT_REL or CT_ADD and field involved is empty
          if (chgtyp == PR_REL || chgtyp == PR_ADD) && isempty(storage.(tbl))
              error('apply_profile: PR_REL or PR_ADD modification cannot be done to storage.%s if it is empty',tbl)
          end

% (E.4) Apply change
          if length(rows) == 1 && rows == 0                 %% modify all rows
            if chgtyp == PR_REP                                    %% replace
              storage.(tbl) = val;
            elseif chgtyp == PR_REL                                %% scale
              storage.(tbl) = val .* storage.(tbl);
            elseif chgtyp == PR_ADD                                %% shift
              storage.(tbl) = val + storage.(tbl);
            else
              error('apply_profile: modification type %d for storage table not supported', chgtyp);
            end
          else                                              %% modify single row
            if chgtyp == PR_REP                                    %% replace
              storage.(tbl)(rows,:) = val;
            elseif chgtyp == PR_REL                                %% scale
              storage.(tbl)(rows,:) = val .* storage.(tbl)(rows,:);
            elseif chgtyp == PR_ADD                                %% shift
              storage.(tbl)(rows,:) = val + storage.(tbl)(rows,:);
            else
              error('apply_profile: modification type %d for storage table not supported', chgtyp);
            end
          end
        argout = storage;
end
