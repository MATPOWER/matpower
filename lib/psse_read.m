function data = psse_read(rawfile_name)
%PSSE_READ  Reads the data from a PSS/E RAW data file.
%   DATA = PSSE_READ(RAWFILE_NAME)
%
%   Reads the data from a PSS/E RAW data file into a struct.
%
%   Input:
%       RAWFILE_NAME : the name of the PSS/E RAW file to be converted
%           (opened directly with fopen)
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
%   See also PSSE2MPC.

%   MATPOWER
%   $Id$
%   by Yujia Zhu, PSERC ASU
%   and Ray Zimmerman, PSERC Cornell
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

%% universal textscan options
tsopt0 = {'Whitespace', ' \b\t,\r\n'};
tsopt  = {tsopt0{:}, 'CommentStyle', ' /'};
tsopt1 = {tsopt{:}, 'CollectOutput', 1};
verbose = 1;

%% open the file
fid = fopen(rawfile_name);
if fid == -1
    error('psse_read: The file ''%s'' could not be found.', rawfile_name);
end

%%-----  read case identification data  -----
%% v29-30:
%%  IC, SBASE
%% v31-33:
%%  IC, SBASE, REV, XFRRAT, NXFRAT, BASFRQ
%% In some files, SBASE is followed by a comment with the REV number such as:
%%  / PSS/E-29.0
%%  / PSS(tm)E-30 RAW
if verbose
    fprintf('Reading case identification data ...');
end
ln = fgets(fid);
id = textscan(ln,'%d %f %[^\n]', 1, tsopt0{:});
%id = textscan(ln,'%d %f %d %f %f %f', 1, tsopt{:});
data.id.IC = id{1}(1);
data.id.SBASE = id{2}(1);
data.id.REV = 0;
data.id.XFRRAT = 0;
data.id.NSFRAT = 0;
data.id.BASFRQ = 0;
data.id.comment0 = '';
if ~isempty(id{3})
    rev = textscan(id{3}{1}, '%d %[^\n]', 1, tsopt0{:});
    if ~isempty(rev{1})
        data.id.REV = rev{1}(1);
        r = textscan(rev{2}{1}, '%f %f %f %[^\n]', 1, tsopt0{:});
        if ~isempty(r{1})
            data.id.XFRRAT = r{1}(1);
            data.id.NSFRAT = r{2}(1);
            data.id.BASFRQ = r{3}(1);
            data.id.comment0 = r{4}(1);
        end
    else
        %% attempt to extract revision from comment
        data.id.comment0 = id{3}{1};
        rev = regexp(data.id.comment0, 'PSS(/|\(tm\))E-(?<rev>\d+)', 'tokens');
        if ~isempty(rev) && size(rev{1}, 2) == 2
            data.id.REV = str2num(rev{1}{2});
        end
    end
end
data.id.comment1 = fgetl(fid);
data.id.comment2 = fgetl(fid);
[data.id.num, data.id.txt] = psse_extract_data(id);
rev = data.id.REV;
if verbose
    if rev
        fprintf(' rev %d format ...', rev);
    else
        fprintf(' rev not specified ...');
    end
    fprintf(' done.\n');
end
if data.id.IC ~= 0
    fprintf('WARNING: IC = %d indicates that this may be a change case, rather than base case\n         PSSE2MPC does NOT handle change cases.\n');
end
if data.id.XFRRAT > 0
    fprintf('WARNING: PSSE2MPC does not correctly handle XFRRAT > 0.\n');
end
if data.id.NSFRAT > 0
    fprintf('WARNING: PSSE2MPC does not correctly handle NSFRAT > 0.\n');
end

%%-----  read bus data  -----
%% v29-30:
%%  I, 'NAME', BASKV, IDE, GL, BL, AREA, ZONE, VM, VA, OWNER
%% v31:
%%  I, 'NAME', BASKV, IDE, AREA, ZONE, OWNER, VM, VA
%% v33:
%%  I, 'NAME', BASKV, IDE, AREA, ZONE, OWNER, VM, VA, NVHI, NVLO, EVHI, EVLO
%% Note: Some v33 files seem to follow v31 format
if rev < 31
    [bus, cnt] = psse_read_section(fid, 'bus', '%d ''%[^'']'' %f %d %f %f %d %d %f %f %d', verbose);
else
    [bus, cnt] = psse_read_section(fid, 'bus', {'%d ''%[^'']'' %f %d %d %d %d %f %f %[^\n]', '%f %f %[^\n]', '%f %f'}, verbose);
%   [bus, cnt] = psse_read_section(fid, 'bus', '%d %*[''"]%[^''"]%*[''"] %f %d %d %d %d %f %f', verbose);
end
[data.bus.num, data.bus.txt] = psse_extract_data(bus);

%%-----  read load data  -----
%% v29-31:
%%  I, ID, STATUS, AREA, ZONE, PL, QL, IP, IQ, YP, YQ, OWNER
%% v33:
%%  I, ID, STATUS, AREA, ZONE, PL, QL, IP, IQ, YP, YQ, OWNER, SCALE, INTRPT
%% Note: Some v33 files seem to end with SCALE
[load, cnt] = psse_read_section(fid, 'load', {'%d ''%[^'']'' %d %d %d %f %f %f %f %f %f %d %[^\n]', '%d %[^\n]', '%d'}, verbose);
[data.load.num, data.load.txt] = psse_extract_data(load);

%%-----  read fixed shunt data  -----
%% v31-33:
%%  I, ID, STATUS, GL, BL
if rev > 30     %% fixed shunt data is included in bus data for rev <= 30
    [shunt, cnt] = psse_read_section(fid, 'fixed bus shunt', '%d ''%[^'']'' %d %f %f', verbose);
    [data.shunt.num, data.shunt.txt] = psse_extract_data(shunt);
end

%%-----  read generator data  -----
%% v29-30:
%%  I,ID,PG,QG,QT,QB,VS,IREG,MBASE,ZR,ZX,RT,XT,GTAP,STAT,RMPCT,PT,PB,O1,F1,...,O4,F4
%% v31-33:
%%  I,ID,PG,QG,QT,QB,VS,IREG,MBASE,ZR,ZX,RT,XT,GTAP,STAT,RMPCT,PT,PB,O1,F1,...,O4,F4,WMOD,WPF
%[gen, cnt] = psse_read_section(fid, 'generator', '%d ''%[^'']'' %f %f %f %f %f %d %f %f %f %f %f %f %d %f %f %f %d %f %d %f %d %f %d %f %d %f', verbose);
[gen, cnt] = psse_read_section(fid, 'generator', {'%d ''%[^'']'' %f %f %f %f %f %d %f %f %f %f %f %f %d %f %f %f %d %f %d %f %d %f %d %f %[^\n]', '%d %f'}, verbose);
[data.gen.num, data.gen.txt] = psse_extract_data(gen);

%%-----  read non-transformer branch data  -----
%% v29-30:
%%  I,J,CKT,R,X,B,RATEA,RATEB,RATEC,GI,BI,GJ,BJ,ST,LEN,O1,F1,...,O4,F4
%% v31-33:
%%  I,J,CKT,R,X,B,RATEA,RATEB,RATEC,GI,BI,GJ,BJ,ST,MET,LEN,O1,F1,...,O4,F4
if 1    %% skip everything after ST for speed
    [branch, cnt] = psse_read_section(fid, 'non-transformer branch', '%d %d ''%[^'']'' %f %f %f %f %f %f %f %f %f %f %d %[^\n]', verbose);
else    %% read all available data
    if rev < 31
        [branch, cnt] = psse_read_section(fid, 'non-transformer branch', {'%d %d ''%[^'']'' %f %f %f %f %f %f %f %f %f %f %d %[^\n]', '%f %d %f %d %f %d %f %d %f'}, verbose);
    else
        [branch, cnt] = psse_read_section(fid, 'non-transformer branch', {'%d %d ''%[^'']'' %f %f %f %f %f %f %f %f %f %f %d %[^\n]', '%d %f %d %f %d %f %d %f %d %f'}, verbose);
    end
end
[data.branch.num, data.branch.txt] = psse_extract_data(branch);

%%-----  read transformer data  -----
%% PSS/E stores two winding and three winding transformer data in the same
%% section in RAW file. We read in 2 passes, first pass determines the type of
%% each, second pass reads the data.
%% Step 1 : Count and collect transformer types
if verbose
    fprintf('Reading transformer data ...');
end
%position = ftell(fid); % save pointer to beginning of transformer data

i2 = 0;
i3 = 0;

%% initialize transformer data
data.trans2.num = NaN(1000, 41);
data.trans2.txt = cell(1000, 41);
data.trans3.num = NaN(1000, 79);
data.trans3.txt = cell(1000, 79);

%% read record 1
[ln, lt] = fgets(fid);

while isempty(regexp(ln(1:end-length(lt)), '^\s*0[\s]*(/.*)?$'))   %% look for 0 record (END)
    %% parse record 1
    %% v29-31:
    %%  I,J,K,CKT,CW,CZ,CM,MAG1,MAG2,NMETR,'NAME',STAT,O1,F1,...,O4,F4
    %% v33:
    %%  I,J,K,CKT,CW,CZ,CM,MAG1,MAG2,NMETR,'NAME',STAT,O1,F1,...,O4,F4,VECGRP
    %% parse everything up to STAT
    r1 = textscan(ln, '%d %d %d %4c %d %d %d %f %f %d %10c %d %[^\n]', 1, tsopt{:});
    %% parse owner data
    if ~isempty(r1{13})
        r1b = textscan(r1{13}{1}, '%d %f %d %f %d %f %d %f %[^\n]', 1, tsopt{:});
    else
        r1b = [];
    end
    if r1{3}(1) == 0    %% K == 0, two-winding
        %% double size if necessary
        if i2 == size(data.trans2.num, 1)
            data.trans2.num = [ data.trans2.num; NaN(size(data.trans2.num)) ];
            data.trans2.txt{2*i2, 1} = [];
        end
        i2 = i2 + 1;
        
        %% store record 1
        for k = 1:length(r1)-1
            if ischar(r1{k})
                data.trans2.txt{i2, k} = r1{k};
                data.trans2.num(i2, k) = NaN;
            elseif iscell(r1{k})
                data.trans2.txt{i2, k} = r1{k}{1};
                data.trans2.num(i2, k) = NaN;
            else
                data.trans2.txt{i2, k} = '';
                data.trans2.num(i2, k) = r1{k};
            end
        end
        if ~isempty(r1b)
            kk = find(~cellfun(@isempty, r1b));
            data.trans2.num(i2, 12+(kk)) = cellfun(@double, r1b(kk));
        end
        
        %% read, parse & store record 2
        %% v29-33:
        %%  R1-2,X1-2,SBASE1-2
        [ln, lt] = fgets(fid);
        r2 = textscan(ln, '%f %f %f', 1, tsopt1{:});
        data.trans2.num(i2, 20+(1:3)) = r2{1};
        
        %% read, parse & store record 3
        %% v29-31:
        %%  WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1
        %% v33:
        %%  WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1,CNXA1
        [ln, lt] = fgets(fid);
        %% parse up to CX1
        r3 = textscan(ln, '%f %f %f %f %f %f %d %d %f %f %f %f %d %d %f %f %[^\n]', 1, tsopt{:});
        kk = find(~cellfun(@isempty, r3));
        data.trans2.num(i2, 23+(kk)) = cellfun(@double, r3(kk));
%% warn if CNXA1 is present and non-zero
        
        %% read, parse & store record 4
        %% v29-33:
        %%  WINDV2,NOMV2
        [ln, lt] = fgets(fid);
        r4 = textscan(ln, '%f %f', 1, tsopt1{:});
        data.trans2.num(i2, 39+(1:2)) = r4{1};
    else                %% K ~= 0, three-winding
        %% double size if necessary
        if i3 == size(data.trans3.num, 1)
            data.trans3.num = [ data.trans3.num; NaN(size(data.trans3.num)) ];
            data.trans3.txt{2*i3, 1} = [];
        end
        i3 = i3 + 1;

        %% store record 1
        for k = 1:length(r1)-1
            if ischar(r1{k})
                data.trans3.txt{i3, k} = r1{k};
                data.trans3.num(i3, k) = NaN;
            elseif iscell(r1{k})
                data.trans3.txt{i3, k} = r1{k}{1};
                data.trans3.num(i3, k) = NaN;
            else
                data.trans3.txt{i3, k} = '';
                data.trans3.num(i3, k) = r1{k};
            end
        end
        if ~isempty(r1b)
            kk = find(~cellfun(@isempty, r1b));
            data.trans3.num(i3, 12+(kk)) = cellfun(@double, r1b(kk));
        end
        
        %% read, parse & store record 2
        %% v29-33:
        %%  R1-2,X1-2,SBASE1-2,R2-3,X2-3,SBASE2-3,R3-1,X3-1,SBASE3-1,VMSTAR,ANSTAR
        [ln, lt] = fgets(fid);
        r2 = textscan(ln, '%f %f %f %f %f %f %f %f %f %f %f', 1, tsopt1{:});
        data.trans3.num(i3, 20+(1:11)) = r2{1};
        
        %% read, parse & store record 3
        %% v29-31:
        %%  WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1
        %% v33:
        %%  WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1,CNXA1
        [ln, lt] = fgets(fid);
        %% parse up to CX1
        r3 = textscan(ln, '%f %f %f %f %f %f %d %d %f %f %f %f %d %d %f %f %[^\n]', 1, tsopt{:});
        kk = find(~cellfun(@isempty, r3));
        data.trans3.num(i3, 31+(kk)) = cellfun(@double, r3(kk));
        
        %% read, parse & store record 4
        %% v29-31:
        %%  WINDV2,NOMV2,ANG2,RATA2,RATB2,RATC2,COD2,CONT2,RMA2,RMI2,VMA2,VMI2,NTP2,TAB2,CR2,CX2
        %% v33:
        %%  WINDV2,NOMV2,ANG2,RATA2,RATB2,RATC2,COD2,CONT2,RMA2,RMI2,VMA2,VMI2,NTP2,TAB2,CR2,CX2,CNXA2
        [ln, lt] = fgets(fid);
        %% parse up to CX2
        r4 = textscan(ln, '%f %f %f %f %f %f %d %d %f %f %f %f %d %d %f %f %[^\n]', 1, tsopt{:});
        kk = find(~cellfun(@isempty, r4));
        data.trans3.num(i3, 47+(kk)) = cellfun(@double, r4(kk));
        
        %% read, parse & store record 5
        %% v29-31:
        %%  WINDV3,NOMV3,ANG3,RATA3,RATB3,RATC3,COD3,CONT3,RMA3,RMI3,VMA3,VMI3,NTP3,TAB3,CR3,CX3
        %% v33:
        %%  WINDV3,NOMV3,ANG3,RATA3,RATB3,RATC3,COD3,CONT3,RMA3,RMI3,VMA3,VMI3,NTP3,TAB3,CR3,CX3,CNXA3
        [ln, lt] = fgets(fid);
        %% parse up to CX3
        r5 = textscan(ln, '%f %f %f %f %f %f %d %d %f %f %f %f %d %d %f %f %[^\n]', 1, tsopt{:});
        kk = find(~cellfun(@isempty, r5));
        data.trans3.num(i3, 63+(kk)) = cellfun(@double, r5(kk));
    end

    %% read record 1
    [ln, lt] = fgets(fid);
end
%% trim down to size
data.trans2.num = data.trans2.num(1:i2,:);
data.trans2.txt = data.trans2.txt(1:i2,:);
data.trans3.num = data.trans3.num(1:i3,:);
data.trans3.txt = data.trans3.txt(1:i3,:);
if verbose
    fprintf(' %d records ... done.\n', i2+i3);
end

%%-----  read area interchange data  -----
%% v29-33:
%%  I, ISW, PDES, PTOL, 'ARNAME'
[area, cnt] = psse_read_section(fid, 'area interchange', '%d %d %f %f ''%[^'']''', verbose);
[data.area.num, data.area.txt] = psse_extract_data(area);

%%-----  read two-terminal DC transmission line data  -----
%% v29-30:
%%  I,MDC,RDC,SETVL,VSCHD,VCMOD,RCOMP,DELTI,METER,DCVMIN,CCCITMX,CCCACC<\n>IPR,NBR,ALFMX,ALFMN,RCR,XCR,EBASR,TRR,TAPR,TMXR,TMNR,STPR,ICR,IFR,ITR,IDR,XCAPR<\n>IPI,NBI,GAMMX,GAMMN,RCI,XCI,EBASI,TRI,TAPI,TMXI,TMNI,STPI,ICI,IFI,ITI,IDI,XCAPI
%% v31-33:
%%  'NAME',MDC,RDC,SETVL,VSCHD,VCMOD,RCOMP,DELTI,METER,DCVMIN,CCCITMX,CCCACC<\n>IPR,NBR,ANMXR,ANMNR,RCR,XCR,EBASR,TRR,TAPR,TMXR,TMNR,STPR,ICR,IFR,ITR,IDR,XCAPR<\n>IPI,NBI,ANMXI,ANMNI,RCI,XCI,EBASI,TRI,TAPI,TMXI,TMNI,STPI,ICI,IFI,ITI,IDI,XCAPI
if rev < 31
    [twodc, cnt] = psse_read_section(fid, 'two-terminal DC transmission line', '%d %d %f %f %f %f %f %f ''%[^'']'' %f %d %f %d %d %f %f %f %f %f %f %f %f %f %f %d %d %d ''%[^'']'' %f %d %d %f %f %f %f %f %f %f %f %f %f %d %d %d ''%[^'']'' %f', verbose, 3);
elseif rev < 32
    [twodc, cnt] = psse_read_section(fid, 'two-terminal DC transmission line', '''%[^'']'' %d %f %f %f %f %f %f ''%[^'']'' %f %d %f %d %d %f %f %f %f %f %f %f %f %f %f %d %d %d %d %f %d %d %f %f %f %f %f %f %f %f %f %f %d %d %d %d %f', verbose, 3);
else
    [twodc, cnt] = psse_read_section(fid, 'two-terminal DC transmission line', '''%[^'']'' %d %f %f %f %f %f %f %1c %f %d %f %d %d %f %f %f %f %f %f %f %f %f %f %d %d %d %d %f %d %d %f %f %f %f %f %f %f %f %f %f %d %d %d %d %f', verbose, 3);
end
[data.twodc.num, data.twodc.txt] = psse_extract_data(twodc);

%%-----  skip voltage source converter data  -----
%% v29-33:
%%  'NAME', MDC, RDC, O1, F1, ... O4, F4<\n>IBUS,TYPE,MODE,DCSET,ACSET,ALOSS,BLOSS,MINLOSS,SMAX,IMAX,PWF,MAXQ,MINQ,REMOT,RMPCT<\n>IBUS,TYPE,MODE,DCSET,ACSET,ALOSS,BLOSS,MINLOSS,SMAX,IMAX,PWF,MAXQ,MINQ,REMOT,RMPCT
if verbose
    fprintf('Skipping %s data ...', 'voltage source converter');
end
[cnt, endpos, llnn] = psse_count_lines(fid, 'dont_rewind');  %% voltage source converter data
if verbose
    fprintf(' %d lines.\n', cnt);
end

%%-----  read switched shunt data  -----
%% v29:
%%  I, MODSW, VSWHI, VSWLO, SWREM, RMIDNT, BINIT, N1, B1, N2, B2, ... N8, B8
%% v30:
%%  I, MODSW, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
%%  I, MODSW, VSWHI, VSWLO, SWREM,        'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
%% v31-33:
%%  switched shunt section follows FACTS device data
if rev < 31
    %% read up to B1
    if rev < 30
        [swshunt, cnt] = psse_read_section(fid, 'switched shunt', '%d %d %f %f %d %10c %f %d %f %[^\n]', verbose);
    else
        [swshunt, cnt] = psse_read_section(fid, 'switched shunt', '%d %d %f %f %d %f %10c %f %d %f %[^\n]', verbose);
    end
    [data.swshunt.num, data.swshunt.txt] = psse_extract_data(swshunt);
end

%%-----  skip impedance correction data  -----
%% v29-33:
%%  I, T1, F1, T2, F2, T3, F3, ... T11, F11
if verbose
    fprintf('Skipping %s data ...', 'impedance correction');
end
[cnt, endpos, llnn] = psse_count_lines(fid, 'dont_rewind');  %% impedance correction data

%%-----  skip multi-terminal DC data  -----
%% v29-30:
%%       I, NCONV, NDCBS, NDCLN, MDC, VCONV, VCMOD, VCONVN<\n>IB,N,ANGMX,ANGMN,RC,XC,EBAS,TR,TAP,TPMX,TPMN,TSTP,SETVL,DCPF,MARG,CNVCOD<\n>IDC, IB, IA,   ZONE, 'NAME',   IDC2, RGRND, OWNER<\n>IDC, JDC, DCCKT, RDC, LDC
%% v31-33:
%%  'NAME', NCONV, NDCBS, NDCLN, MDC, VCONV, VCMOD, VCONVN<\n>IB,N,ANGMX,ANGMN,RC,XC,EBAS,TR,TAP,TPMX,TPMN,TSTP,SETVL,DCPF,MARG,CNVCOD<\n>IDC, IB, AREA, ZONE, 'DCNAME', IDC2, RGRND, OWNER<\n>IDC, JDC, DCCKT, MET, RDC, LDC
if verbose
    fprintf(' %d lines.\n', cnt);
    fprintf('Skipping %s data ...', 'multi-terminal DC transmission line');
end
[cnt, endpos, llnn] = psse_count_lines(fid, 'dont_rewind');  %% multi-terminal DC data

%%-----  skip multi-section line data  -----
%% v29-30:
%%  I, J, ID, DUM1, DUM2, ... DUM9
%% v31-33:
%%  I, J, ID, MET, DUM1, DUM2, ... DUM9
if verbose
    fprintf(' %d lines.\n', cnt);
    fprintf('Skipping %s data ...', 'multi-section line grouping');
end
[cnt, endpos, llnn] = psse_count_lines(fid, 'dont_rewind');  %% multi-section line data

%%-----  skip zone data  -----
%% v29-33:
%%  I, 'ZONAME'
if verbose
    fprintf(' %d lines.\n', cnt);
    fprintf('Skipping %s data ...', 'zone');
end
[cnt, endpos, llnn] = psse_count_lines(fid, 'dont_rewind');  %% zone data

%%-----  skip inter-area transfer data  -----
%% v29-33:
%%  ARFROM, ARTO, TRID, PTRAN
if verbose
    fprintf(' %d lines.\n', cnt);
    fprintf('Skipping %s data ...', 'inter-area transfer');
end
[cnt, endpos, llnn] = psse_count_lines(fid, 'dont_rewind');  %% inter-area transfer data

%%-----  skip owner data  -----
%% v29-33:
%%  I, 'OWNAME'
if verbose
    fprintf(' %d lines.\n', cnt);
    fprintf('Skipping %s data ...', 'owner');
end
[cnt, endpos, llnn] = psse_count_lines(fid, 'dont_rewind');  %% owner data

%%-----  skip FACTS control device data  -----
%% v29:
%%  N,I,J,MODE,PDES,QDES,VSET,SHMX,TRMX,VTMN,VTMX,VSMX,IMX,LINX,OWNER,SET1,SET2,VSREF
%% v30:
%%  N,I,J,MODE,PDES,QDES,VSET,SHMX,TRMX,VTMN,VTMX,VSMX,IMX,LINX,RMPCT,OWNER,SET1,SET2,VSREF
%% v31-33:
%%  'NAME',I,J,MODE,PDES,QDES,VSET,SHMX,TRMX,VTMN,VTMX,VSMX,IMX,LINX,RMPCT,OWNER,SET1,SET2,VSREF,REMOT,'MNAME'
if verbose
    fprintf(' %d lines.\n', cnt);
    fprintf('Skipping %s data ...', 'FACTS control device');
end
[cnt, endpos, llnn] = psse_count_lines(fid, 'dont_rewind');  %% FACTS control device data
if verbose
    fprintf(' %d lines.\n', cnt);
end

%%-----  read switched shunt data  -----
%% v29-30:
%%  switched shunt section follows Voltage Source Converter (VSC) DC Transmission Line data
%% v31:
%%  I, MODSW, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
%% v32:
%%  ? (we will assume v32 is the same as v33, until we find out otherwise)
%% v33:
%%  I, MODSW, ADJM, STAT, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
if rev > 30
    %% read up to B1
    if rev < 32     %% assume 32 is same as 33
        [swshunt, cnt] = psse_read_section(fid, 'switched shunt', '%d %d %f %f %d %f %10c %f %d %f %[^\n]', verbose);
    else
        [swshunt, cnt] = psse_read_section(fid, 'switched shunt', '%d %d %d %d %f %f %d %f %10c %f %d %f %[^\n]', verbose);
    end
    [data.swshunt.num, data.swshunt.txt] = psse_extract_data(swshunt);
end

%% finish parsing switched shunt data
c = size(data.swshunt.txt, 2);
data.swshunt.num(:, c) = [];   %% get rid of NaNs in last col
for i = 1:size(data.swshunt.txt, 1)
    tmp = textscan(data.swshunt.txt{i, c}, '%d %f', tsopt{:});
    n = size(tmp{1}, 1);
    for j = 1:n
        data.swshunt.num(i, c+2*j-2) = tmp{1}(j);
        data.swshunt.num(i, c+2*j-1) = tmp{2}(j);
    end
end

fclose(fid);
