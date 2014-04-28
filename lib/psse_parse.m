function data = psse_parse(records, sections, verbose, rev)
%PSSE_PARSE  Parses the data from a PSS/E RAW data file.
%   DATA = PSSE_PARSE(RECORDS, SECTIONS, VERBOSE)
%
%   Parses the data from a PSS/E RAW data file (as read by PSSE_READ)
%   into a struct.
%
%   Input:
%       RECORDS : cell array of strings, corresponding to the lines
%                 in the RAW file
%       SECTIONS : struct array with indexes marking the beginning
%                  and end of each section, and the name of the
%                  section, fields are:
%           first   : index into RECORDS of first line of section
%           last    : index into RECORDS of last line of section
%           name    : name of the section, as extracted from the
%                     END OF ... DATA comments
%       VERBOSE      :  1 (default) to display progress info, 0 otherwise
%       REV          :  (optional) assume the input file is of this
%                       PSS/E revision number, attempts to determine
%                       REV from the file by default
%
%   Output:
%       DATA :  a struct with the following fields, each with two
%               sub-fields, 'num' and 'txt' containing the numeric and
%               text data read from the file for the corresponding section
%           id
%           bus
%           load
%           gen
%           shunt
%           branch
%           trans2
%           trans3
%           area
%           twodc
%           swshunt
%
%   See also PSSE2MPC, PSSE_READ, PSSE_PARSE_SECTION, PSSE_PARSE_LINE

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Based on mpreadraw.m, written by: Yujia Zhu, Jan 2014, yzhu54@asu.edu.
%   Copyright (c) 2014 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% default args
if nargin < 4
    rev = 0;
    if nargin < 3
        verbose = 0;
    end
end

%% inititialize section counter
s = 1;

%%-----  parse case identification data  -----
%% v29-30:
%%  IC, SBASE
%% v31-33:
%%  IC, SBASE, REV, XFRRAT, NXFRAT, BASFRQ
%% In some files, SBASE is followed by a comment with the REV number such as:
%%  / PSS/E-29.0
%%  / PSS(tm)E-30 RAW
if verbose
    if rev
        fprintf('Forcing interpretation as PSS/E revision %d\n', rev);
    else
        fprintf('Attempting to determine PSS/E revision from content.\n');
    end
    fprintf('Parsing case identification data ...');
end
[d, c] = psse_parse_line(records{1}, 'dfdfff');
nn = length(d);
data.id.IC = d{1};
data.id.SBASE = d{2};
data.id.REV = 0;
data.id.XFRRAT = 0;
data.id.NSFRAT = 0;
data.id.BASFRQ = 0;
if nn > 2
    data.id.REV = d{3};
    if nn > 3
        data.id.XFRRAT = str2double(d{4});
        if nn > 4
            data.id.NXFRAT = str2double(d{5});
            if nn > 5
                data.id.BASFRQ = str2double(d{6});
            end
        end
    end
else    %% attempt to extract revision from comment
    tmp = regexp(c, 'PSS(/|\(tm\))E-(?<rev>\d+)', 'tokens');
    if ~isempty(tmp) && size(tmp{1}, 2) == 2
        data.id.REV = str2num(tmp{1}{2});
    end
end
data.id.comment0 = c;
data.id.comment1 = records{2};
data.id.comment2 = records{3};
if verbose
    if data.id.REV
        fprintf('......... rev %2d format detected ... done.\n', data.id.REV);
    else
        fprintf('.............. rev not specified ... done.\n');
    end
end
if ~rev
    rev = data.id.REV;
end
if data.id.IC ~= 0
    fprintf('WARNING: IC = %d indicates that this may be a change case, rather than base case\n         PSSE2MPC is NOT designed to handle change cases.\n', data.id.IC);
end
if data.id.XFRRAT > 0
    fprintf('WARNING: PSSE2MPC does not correctly handle XFRRAT > 0 [%d].\n', data.id.XFRRAT);
end
if data.id.NSFRAT > 0
    fprintf('WARNING: PSSE2MPC does not correctly handle NSFRAT > 0 [%d].\n', data.id.NSFRAT);
end
s = s + 1;

%%-----  parse bus data  -----
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I, IDE, PL, QL, GL, BL, IA, VM, VA, 'NAME', BASKL, ZONE
%% v29-30:
%%  I, 'NAME', BASKV, IDE, GL, BL, AREA, ZONE, VM, VA, OWNER
%% v31:
%%  I, 'NAME', BASKV, IDE, AREA, ZONE, OWNER, VM, VA
%% v33:
%%  I, 'NAME', BASKV, IDE, AREA, ZONE, OWNER, VM, VA, NVHI, NVLO, EVHI, EVLO
%% Note: Some v33 files seem to follow v31 format
if rev == 1
    data.bus = psse_parse_section(records, sections, s, verbose, ...
        'bus', 'ddffffdffsfd');
elseif rev < 31
    data.bus = psse_parse_section(records, sections, s, verbose, ...
        'bus', 'dsfdffddffd');
else
    data.bus = psse_parse_section(records, sections, s, verbose, ...
        'bus', 'dsfddddffff..');
%       'bus', 'dsfddddffffff');
end
s = s + 1;

%%-----  parse load data  -----
%% v29-31:
%%  I, ID, STATUS, AREA, ZONE, PL, QL, IP, IQ, YP, YQ, OWNER
%% v33:
%%  I, ID, STATUS, AREA, ZONE, PL, QL, IP, IQ, YP, YQ, OWNER, SCALE, INTRPT
%% Note: Some v33 files seem to end with SCALE
if rev > 1
    data.load = psse_parse_section(records, sections, s, verbose, ...
        'load', 'd.d..ffffff...');
%       'load', 'dsdddffffffddd');
    s = s + 1;
end

%%-----  parse fixed shunt data  -----
%% v31-33:
%%  I, ID, STATUS, GL, BL
if rev > 30     %% fixed shunt data is included in bus data for rev <= 30
    data.shunt = psse_parse_section(records, sections, s, verbose, ...
        'fixed shunt', 'd.dff');
%       'fixed shunt', 'dsdff');
    s = s + 1;
end

%%-----  parse generator data  -----
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I,ID,PG,QG,QT,QB,VS,IREG,MBASE,ZR,ZX,RT,XT,GTAP,STAT,RMPCT,PT,PB
%% v29-30:
%%  I,ID,PG,QG,QT,QB,VS,IREG,MBASE,ZR,ZX,RT,XT,GTAP,STAT,RMPCT,PT,PB,O1,F1,...,O4,F4
%% v31-33:
%%  I,ID,PG,QG,QT,QB,VS,IREG,MBASE,ZR,ZX,RT,XT,GTAP,STAT,RMPCT,PT,PB,O1,F1,...,O4,F4,WMOD,WPF
data.gen = psse_parse_section(records, sections, s, verbose, ...
    'generator', 'd.fffff.f.....d.ff...........');
%   'generator', 'dsfffffdffffffdfffdfdfdfdfsdf');
s = s + 1;

%%-----  parse non-transformer branch data  -----
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I,J,CKT,R,X,B,RATEA,RATEB,RATEC,RATIO,ANGLE,GI,BI,GJ,BJ,ST
%% v29-30:
%%  I,J,CKT,R,X,B,RATEA,RATEB,RATEC,GI,BI,GJ,BJ,ST,LEN,O1,F1,...,O4,F4
%% v31-33:
%%  I,J,CKT,R,X,B,RATEA,RATEB,RATEC,GI,BI,GJ,BJ,ST,MET,LEN,O1,F1,...,O4,F4
if rev == 1
    data.branch = psse_parse_section(records, sections, s, verbose, ...
        'branch', 'dd.ffffffffffffd');
%       'branch', 'dddffffffffffffd');
elseif 1    %% skip everything after ST for speed
    data.branch = psse_parse_section(records, sections, s, verbose, ...
        'branch', 'dd.ffffffffffd');
%       'branch', 'ddsffffffffffd');
else    %% read all available data
    if rev < 31
        data.branch = psse_parse_section(records, sections, s, verbose, ...
            'branch', 'dd.ffffffffffd.........');
%           'branch', 'ddsffffffffffdfdfdfdfdf');
    else
        data.branch = psse_parse_section(records, sections, s, verbose, ...
            'branch', 'dd.ffffffffffd..........');
%           'branch', 'ddsffffffffffddfdfdfdfdf');
    end
end
s = s + 1;

%%-----  skip impedance correction data  -----
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I,J,CKT,ICONT,RMA,RMI,VMA,VMI,STEP,TABLE
%% v29-33:
%%  (precedes multi-terminal DC section)
if rev == 1
    s = psse_skip_section(sections, s, verbose, 'impedance correction');
end

%%-----  parse transformer data  -----
if rev > 1
    %% PSS/E stores two winding and three winding transformer data in the same
    %% section in RAW file. We read in 2 passes, first pass determines the type of
    %% each, second pass reads the data.
    label = 'transformer';
    if ~isempty(sections(s).name) && ~strcmp(upper(label), sections(s).name)
        fprintf('-----  WARNING:  Expected section labeled: ''%s''\n', upper(label));
        fprintf('-----            Found section labeled:    ''%s''\n', sections(s).name);
    end

    %% Step 1 : Count and collect transformer types
    if verbose
        fprintf('Analyzing transformer types ...');
    end

    %% estimate max number of transformers by number of lines in section
    nt2 = round((sections(s).last - sections(s).first + 1) / 4);
    nt3 = round((sections(s).last - sections(s).first + 1) / 5);

    %% initialize indexes for transformer types
    idx2 = zeros(nt2, 1);
    idx3 = zeros(nt3, 1);

    %% set up counters
    i = sections(s).first;      %% initialize record index
    i2 = 0;
    i3 = 0;

    while i <= sections(s).last
        %% determine transformer type
        pat = '[^''",\s/]+\s*(,|\s)\s*[^''",\s/]+\s*(,|\s)\s*([^''",\s/]+)';
        m = regexp(records{i}, pat, 'tokens', 'once');
        if length(m) ~= 3
            m
            error('m should be length 3');
        end
        if length(m{3}) == 1 && m{3}(1) == '0'  %% two-winding
            i2 = i2 + 1;
            idx2(i2) = i;
            i = i + 4;
        else                                    %% three-winding
            i3 = i3 + 1;
            idx3(i3) = i;
            i = i + 5;
        end
    end
    nt2 = i2;
    nt3 = i3;

    if verbose
%        str = sprintf(' %d two-winding, %d three-winding.', nt2, nt3);
        str = sprintf(' %d(%d) two(three)-winding.', nt2, nt3);
        spacers = repmat('.', 1, 36-length(str));
        fprintf('%s %s ... done.\n', spacers, str);
    end

    %% trim index vectors down to size
    idx2 = idx2(1:nt2);
    idx3 = idx3(1:nt3);

    %% parse record 1 (cols 1-20)
    %% v29-31:
    %%  I,J,K,CKT,CW,CZ,CM,MAG1,MAG2,NMETR,'NAME',STAT,O1,F1,...,O4,F4
    %% v33:
    %%  I,J,K,CKT,CW,CZ,CM,MAG1,MAG2,NMETR,'NAME',STAT,O1,F1,...,O4,F4,VECGRP
    t2_1 = psse_parse_section(records(idx2), verbose, ...
        '2-winding transformers (1)', 'dd..ddd....d........');
%       '2-winding transformers (1)', 'dddsdddffdsddfdfdfdf');
    t3_1 = psse_parse_section(records(idx3), verbose, ...
        '3-winding transformers (1)', 'ddd.ddd....d........');
%       '3-winding transformers (1)', 'dddsdddffdsddfdfdfdf');

    %% two-winding
    %% parse record 2 (cols 21-23)
    %% v29-33:
    %%  R1-2,X1-2,SBASE1-2
    t2_2 = psse_parse_section(records(idx2+1), verbose, ...
        '2-winding transformers (2)', 'fff');

    %% parse record 3 (cols 24-39)
    %% v29-31:
    %%  WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1
    %% v33:
    %%  WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1,CNXA1
    %% parse up to CX1
    %% warn if CNXA1 is present and non-zero
    t2_3 = psse_parse_section(records(idx2+2), verbose, ...
        '2-winding transformers (3)', 'ffffff..........');
%       '2-winding transformers (3)', 'ffffffddffffddff');

    %% parse record 4 (cols 40-41)
    %% v29-33:
    %%  WINDV2,NOMV2
    t2_4 = psse_parse_section(records(idx2+3), verbose, ...
        '2-winding transformers (4)', 'ff');

    %% three-winding
    %% parse record 2 (cols 21-31)
    %% v29-33:
    %%  R1-2,X1-2,SBASE1-2,R2-3,X2-3,SBASE2-3,R3-1,X3-1,SBASE3-1,VMSTAR,ANSTAR
    t3_2 = psse_parse_section(records(idx3+1), verbose, ...
        '3-winding transformers (2)', 'fffffffffff');
%       '3-winding transformers (2)', 'fffffffffff');

    %% parse record 3 (cols 32-47)
    %% v29-31:
    %%  WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1
    %% v33:
    %%  WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1,CNXA1
    %% parse up to CX1
    %% warn if CNXA1 is present and non-zero
    t3_3 = psse_parse_section(records(idx3+2), verbose, ...
        '3-winding transformers (3)', 'ffffff..........');
%       '3-winding transformers (3)', 'ffffffddffffddff');

    %% parse record 4 (cols 48-63)
    %% v29-31:
    %%  WINDV2,NOMV2,ANG2,RATA2,RATB2,RATC2,COD2,CONT2,RMA2,RMI2,VMA2,VMI2,NTP2,TAB2,CR2,CX2
    %% v33:
    %%  WINDV2,NOMV2,ANG2,RATA2,RATB2,RATC2,COD2,CONT2,RMA2,RMI2,VMA2,VMI2,NTP2,TAB2,CR2,CX2,CNXA2
    %% parse up to CX2
    t3_4 = psse_parse_section(records(idx3+3), verbose, ...
        '3-winding transformers (4)', 'ffffff..........');
%       '3-winding transformers (4)', 'ffffffddffffddff');

    %% parse record 5 (cols 64-79)
    %% v29-31:
    %%  WINDV3,NOMV3,ANG3,RATA3,RATB3,RATC3,COD3,CONT3,RMA3,RMI3,VMA3,VMI3,NTP3,TAB3,CR3,CX3
    %% v33:
    %%  WINDV3,NOMV3,ANG3,RATA3,RATB3,RATC3,COD3,CONT3,RMA3,RMI3,VMA3,VMI3,NTP3,TAB3,CR3,CX3,CNXA3
    %% parse up to CX3
    t3_5 = psse_parse_section(records(idx3+4), verbose, ...
        '3-winding transformers (5)', 'ffffff..........');
%       '3-winding transformers (5)', 'ffffffddffffddff');

    %% assemble two-winding transformer records
    data.trans2.num = [t2_1.num t2_2.num t2_3.num(:, 1:16) t2_4.num];
    data.trans2.txt = [t2_1.txt t2_2.txt t2_3.txt(:, 1:16) t2_4.txt];

    %% assemble three-winding transformer records
    data.trans3.num = [t3_1.num t3_2.num t3_3.num(:, 1:16) t3_4.num(:, 1:16) t3_5.num(:, 1:16)];
    data.trans3.txt = [t3_1.txt t3_2.txt t3_3.txt(:, 1:16) t3_4.txt(:, 1:16) t3_5.txt(:, 1:16)];

    % if verbose
    %     fprintf('%s\n', upper(label));
    %     fprintf('%s\n', sections(s).name);
    % end
    s = s + 1;
end

%%-----  parse area interchange data  -----
%% v1, 29-33:
%%  I, ISW, PDES, PTOL, 'ARNAME'
data.area = psse_parse_section(records, sections, s, verbose, 'area', 'ddffs');
s = s + 1;

%%-----  parse two-terminal DC transmission line data  -----
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I,MDC,RDC,SETVL,VSCHD,VCMOD,RCOMP,DELTI,METER (1-9)
%%  IPR,NBR,ALFMAX,ALFMN,RCR,XCR,EBASR,TRR,TAPR,TPMXR,TPMNR,TSTPR (13-24)
%%  IPI,NBI,GAMMX,GAMMN,RCI,XCI,EBASI,TRI,TAPI,TPMXI,TPMNI,TSTPI (30-41)
%% v29-30:
%%  I,MDC,RDC,SETVL,VSCHD,VCMOD,RCOMP,DELTI,METER,DCVMIN,CCCITMX,CCCACC (1-12)
%%  IPR,NBR,ALFMX,ALFMN,RCR,XCR,EBASR,TRR,TAPR,TMXR,TMNR,STPR,ICR,IFR,ITR,IDR,XCAPR (13-29)
%%  IPI,NBI,GAMMX,GAMMN,RCI,XCI,EBASI,TRI,TAPI,TMXI,TMNI,STPI,ICI,IFI,ITI,IDI,XCAPI (30-46)
%% v31-33:
%%  'NAME',MDC,RDC,SETVL,VSCHD,VCMOD,RCOMP,DELTI,METER,DCVMIN,CCCITMX,CCCACC (1-12)
%%  IPR,NBR,ANMXR,ANMNR,RCR,XCR,EBASR,TRR,TAPR,TMXR,TMNR,STPR,ICR,IFR,ITR,IDR,XCAPR (13-29)
%%  IPI,NBI,ANMXI,ANMNI,RCI,XCI,EBASI,TRI,TAPI,TMXI,TMNI,STPI,ICI,IFI,ITI,IDI,XCAPI (30-46)
label = 'two-terminal DC';
if ~isempty(sections(s).name) && ~strcmp(upper(label), sections(s).name)
    fprintf('-----  WARNING:  Expected section labeled: ''%s''\n', upper(label));
    fprintf('-----            Found section labeled:    ''%s''\n', sections(s).name);
end
idx = sections(s).first:3:sections(s).last;
if rev < 31
    dc1 = psse_parse_section(records(idx), verbose, ...
        'two-terminal DC (1)', '.d.ff.......');
%       'two-terminal DC (1)', 'ddffffffsfdf');
    dc2 = psse_parse_section(records(idx+1), verbose, ...
        'two-terminal DC (2)', 'd.ff.............');
%       'two-terminal DC (2)', 'ddffffffffffdddsf');
    dc3 = psse_parse_section(records(idx+2), verbose, ...
        'two-terminal DC (3)', 'd.ff.............');
%       'two-terminal DC (3)', 'ddffffffffffdddsf');
else
    dc1 = psse_parse_section(records(idx), verbose, ...
        'two-terminal DC (1)', '.d.ff.......');
%       'two-terminal DC (1)', 'sdffffffsfdf');
    dc2 = psse_parse_section(records(idx+1), verbose, ...
        'two-terminal DC (2)', 'd.ff.............');
%       'two-terminal DC (2)', 'ddffffffffffdddDf');
    dc3 = psse_parse_section(records(idx+2), verbose, ...
        'two-terminal DC (3)', 'd.ff.............');
%       'two-terminal DC (3)', 'ddffffffffffdddDf');
end

%% assemble two-terminal DC transmission line
data.twodc.num = [dc1.num dc2.num dc3.num];
data.twodc.txt = [dc1.txt dc2.txt dc3.txt];
% if verbose
%     fprintf('%s\n', upper(label));
%     fprintf('%s\n', sections(s).name);
% end
s = s + 1;

%%-----  skip voltage source converter data  -----
%% v29-33:
%%  'NAME', MDC, RDC, O1, F1, ... O4, F4<\n>IBUS,TYPE,MODE,DCSET,ACSET,ALOSS,BLOSS,MINLOSS,SMAX,IMAX,PWF,MAXQ,MINQ,REMOT,RMPCT<\n>IBUS,TYPE,MODE,DCSET,ACSET,ALOSS,BLOSS,MINLOSS,SMAX,IMAX,PWF,MAXQ,MINQ,REMOT,RMPCT
if rev > 1
    s = psse_skip_section(sections, s, verbose, 'voltage source converter');
end

%%-----  parse switched shunt data  -----
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I, MODSW, VSWHI, VSWLO, SWREM, BINIT, N1, B1, N2, B2, ... N8, B8
%% v29:
%%  I, MODSW, VSWHI, VSWLO, SWREM, RMIDNT, BINIT, N1, B1, N2, B2, ... N8, B8
%% v30:
%%  I, MODSW, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
%% v31-33:
%%  (switched shunt section follows FACTS device data)
if rev < 31
    %% parse up to B1
    if rev == 1
        data.swshunt = psse_parse_section(records, sections, s, verbose, ...
            'switched shunt', 'd....f');
%           'switched shunt', 'ddffdfdfdfdfdfdfdfdfdf');
    elseif rev <= 29
        data.swshunt = psse_parse_section(records, sections, s, verbose, ...
            'switched shunt', 'd.....f');
%           'switched shunt', 'ddffdsfdfdfdfdfdfdfdfdf');
    else    %%  rev == 30
        data.swshunt = psse_parse_section(records, sections, s, verbose, ...
            'switched shunt', 'd......f');
%           'switched shunt', 'ddffdfsfdfdfdfdfdfdfdfdf');
    end
    s = s + 1;
end

%%-----  skip impedance correction data  -----
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  (follows branch data)
%% v29-33:
%%  I, T1, F1, T2, F2, T3, F3, ... T11, F11
if rev > 1
    s = psse_skip_section(sections, s, verbose, 'impedance correction');
end

%%-----  skip multi-terminal DC data  -----
%% v29-30:
%%       I, NCONV, NDCBS, NDCLN, MDC, VCONV, VCMOD, VCONVN<\n>IB,N,ANGMX,ANGMN,RC,XC,EBAS,TR,TAP,TPMX,TPMN,TSTP,SETVL,DCPF,MARG,CNVCOD<\n>IDC, IB, IA,   ZONE, 'NAME',   IDC2, RGRND, OWNER<\n>IDC, JDC, DCCKT, RDC, LDC
%% v31-33:
%%  'NAME', NCONV, NDCBS, NDCLN, MDC, VCONV, VCMOD, VCONVN<\n>IB,N,ANGMX,ANGMN,RC,XC,EBAS,TR,TAP,TPMX,TPMN,TSTP,SETVL,DCPF,MARG,CNVCOD<\n>IDC, IB, AREA, ZONE, 'DCNAME', IDC2, RGRND, OWNER<\n>IDC, JDC, DCCKT, MET, RDC, LDC
if rev > 1
    s = psse_skip_section(sections, s, verbose, 'multi-terminal DC');
end

%%-----  skip multi-section line data  -----
%% v29-30:
%%  I, J, ID, DUM1, DUM2, ... DUM9
%% v31-33:
%%  I, J, ID, MET, DUM1, DUM2, ... DUM9
if rev > 1
    s = psse_skip_section(sections, s, verbose, 'multi-section line');
end

%%-----  skip zone data  -----
%% v29-33:
%%  I, 'ZONAME'
if rev > 1
    s = psse_skip_section(sections, s, verbose, 'zone');
end

%%-----  skip inter-area transfer data  -----
%% v29-33:
%%  ARFROM, ARTO, TRID, PTRAN
if rev > 1
    s = psse_skip_section(sections, s, verbose, 'inter-area transfer');
end

%%-----  skip owner data  -----
%% v29-33:
%%  I, 'OWNAME'
if rev > 1
    s = psse_skip_section(sections, s, verbose, 'owner');
end

%%-----  skip FACTS control device data  -----
%% v29:
%%  N,I,J,MODE,PDES,QDES,VSET,SHMX,TRMX,VTMN,VTMX,VSMX,IMX,LINX,OWNER,SET1,SET2,VSREF
%% v30:
%%  N,I,J,MODE,PDES,QDES,VSET,SHMX,TRMX,VTMN,VTMX,VSMX,IMX,LINX,RMPCT,OWNER,SET1,SET2,VSREF
%% v31-33:
%%  'NAME',I,J,MODE,PDES,QDES,VSET,SHMX,TRMX,VTMN,VTMX,VSMX,IMX,LINX,RMPCT,OWNER,SET1,SET2,VSREF,REMOT,'MNAME'
if rev > 1
    s = psse_skip_section(sections, s, verbose, 'FACTS control device');
end

%%-----  parse switched shunt data  -----
%% v29-30:
%%  switched shunt section follows Voltage Source Converter (VSC) DC Transmission Line data
%% v31:
%%  I, MODSW, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
%% v32:
%%  ? (we will assume v32 is the same as v33, until we find out otherwise)
%% v33:
%%  I, MODSW, ADJM, STAT, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
if rev > 30
    %% parse up to B1
    if rev < 32     %% assume 32 is same as 33
        data.swshunt = psse_parse_section(records, sections, s, verbose, ...
            'switched shunt', 'd......f');
%           'switched shunt', 'ddffdfsfdfdfdfdfdfdfdfdf');
    else
        data.swshunt = psse_parse_section(records, sections, s, verbose, ...
            'switched shunt', 'd........f');
%           'switched shunt', 'ddddffdfsfdfdfdfdfdfdfdfdf');
    end
    s = s + 1;
end

%%-----  skip GNE device data  -----
%% v32: (includes this section, probably same as v33)
%% v33:
%%  'NAME','MODEL',NTERM, BUS1, ...,BUSNTERM,NREAL,NINTG,NCHAR, STATUS, OWNER,NMETR
%%    REAL1, ..., REALmin(10,NREAL)
%%    INTG1, ..., INTGmin(10,NINTG)
%%    CHAR1, ..., CHARmin(10,NCHAR)
if rev > 31
    s = psse_skip_section(sections, s, verbose, 'GNE device');
end

%%-----  skip induction machine data  -----
%% v33:
%%  I,ID,STAT,SCODE,DCODE,AREA,ZONE,OWNER,TCODE,BCODE,MBASE,RATEKV, PCODE,PSET,H,A,B,D,E,RA,XA,XM,R1,X1,R2,X2,X3,E1,SE1,E2,SE2,IA1,IA2, XAMULT
if rev > 32
    s = psse_skip_section(sections, s, verbose, 'induction machine');
end

%%-----  check for extra sections  -----
if s <= length(sections)
    fprintf('-----  WARNING:  Found %d additional section(s):\n', length(sections)-s+1);
end
while s <= length(sections)
    n = sections(s).last - sections(s).first + 1;
    if n
        str = sprintf('with %d line(s)', n);
    else
        str = sprintf('(empty)');
    end
    if isempty(sections(s).name)
        fprintf('-----            unlabeled section %s\n', str);
    else
        fprintf('-----            ''%s DATA'' %s\n', sections(s).name, str);
    end
    s = s + 1;
end



%%---------------------------------------------------------------------
function s = psse_skip_section(sections, s, verbose, label)
if s > length(sections)
    if verbose
        spacers = repmat('.', 1, 58-length(label));
        fprintf('No %s data read %s done.\n', label, spacers);
    end
else
    if ~isempty(sections(s).name) && ~strcmp(upper(label), sections(s).name)
        fprintf('-----  WARNING:  Expected section labeled: ''%s''\n', upper(label));
        fprintf('-----            Found section labeled:    ''%s''\n', sections(s).name);
    end
    if verbose
        spacers = repmat('.', 1, 42-length(label));
        nr = sections(s).last - sections(s).first + 1;
        fprintf('Skipping%6d lines of %s data %s done.\n', nr, label, spacers);
    %     fprintf('%s\n', upper(label));
    %     fprintf('%s\n', sections(s).name);
    end
    s = s + 1;
end
