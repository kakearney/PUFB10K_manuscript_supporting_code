% popup_analysis_final
%
% This script covers the data analysis and plotting underlying the
% NPRB-funded (Rogers et al.) P.Cod popup Obj. 2 manuscript.  
% 
% While currently written with ROMS model output in mind, portions can also
% be adapted to rerun a similar analysis for future ROMS and MOM6
% simulations. I've noted those locations where necessary.
%
% The following non-native Matlab functions are used in this script, and
% should be downloaded and installed on the user's path in before running.
% The majority of these are called in the context of creating custom
% graphics for the paper's figures. 
%
%  From the MatlabCentral File Exchange:
%
%   3rd-party utilities:
%
%   - Climate Data Toolbox (borders, cmocean, ncstruct, shadem, rgb)
%     FEX #70338
%   - RunLength
%     FEX #41813
%   - Crameri perceptually uniform scientific colormaps (crameri)
%     FEX #68546
%   - hex2rgb & rgb2hex (hex2rgb)
%     FEX #46289 (not required for R2024a+)
%   - Hungarian Algorithm for Linear Assignment Problems (munkres)
%     FEX #20652
%   - export_fig
%     FEX #23629
%   - xy2sn, sn2xy
%     FEX #39796
%   - ConsoleProgressBar
%     FEX #30297
%   - arclength
%     FEX #34871
%   - axescoord2figurecoord
%     FEX #13634
%   - distance2curve
%     FEX #34869
%   - interparc
%     FEX #34874
%   - subaxis
%     FEX #3696
%
%   Author's utilities:
%
%   - labelaxes
%   - mergeaxes
%   - minmax
%   - plotboxpos
%     FEX #9615
%   - plotgrid
%   - shapeprjread
%   - skillstats
%
% From GitHub:
%
%   - roms (calcromsz, plotromsrho, stretching)
%   - boxmap (boxworldmap)
%   - Bering10KPostprocessing (coldpool_surveyreplicate)
%
%  Other
%   - SeaWater v3.2 (sw_dens0, sw_dpth, sw_smow)


%% Setup

%% ... helper functions

dfun = @(a,b) distance(a, b, referenceEllipsoid('earth', 'km'));

ff = @(x) reshape(fullfile({x.folder}, {x.name}), size(x));

mergestruct = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

%% ... locations of various datasets

% Working directory used as base for filenames

wd = pwd;

%----------------------------
% Data accessed via web
%----------------------------

% Post-processed ROMS simulation output

ppsimdir = 'https://data.pmel.noaa.gov/aclim/thredds/dodsC/B10K-K20nobio_CORECFS_daily/'; 
% ppsimdir = fullfile(moxdir, 'roms_for_public', 'B10K-K20nobio_CORECFS_daily'); % post-processed output from ROMS simulation
reloadb10kdata = false; % True to read data, false assumes this step has already been run and saved to the b10ksimdatafile .mat file
b10ksimdatafile = 'simdata_1980-2023_btemp_stemp.mat';

% ROMS grid file

grdfile = 'https://data.pmel.noaa.gov/aclim/thredds/dodsC/ancillary/Bering10K_extended_grid.nc';
% grdfile = fullfile(moxdir, 'roms_for_public', 'Bering10K_extended_grid.nc');

% ROMS survey-replicated file
% Note: This path currently points to the TDS-hosted copy of this file.  To
% recreate this file (i.e. rerunsrepflag = true), the path needs to be
% changed to a writable location.

srepfile = strrep(fullfile(ppsimdir, 'Level3', 'survey_replicates_B10K-K20nobio_CORECFS_daily.csv'), 'dodsC', 'fileServer/files');
srepsimname = 'model_bottom_temp';
rerunsrepflag = false; % True to redo survey-rep, false assumes its done already

% ETOPO 2022

etopotilename = 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/60s_surface_elev_netcdf/ETOPO_2022_v1_60s_N90W180_surface.nc';

%----------------------------
% Data included in repo
%----------------------------

% Popup float .csv files

pufdatafol = fullfile(wd, 'PUF_Data'); % PUF .csv files

% Mooring datasets: includes 1-hr gridded data with bering, kuskokwin, and
% innerfront subfolder

mooringfol = fullfile(wd, 'PMEL_Mooring_Data', '1hr Gridded');

% Bottom temperature rasters (created via coldpool R package, see details
% in folder)

srepyr = [1982:2019 2021:2023];

rasterfol = fullfile(wd, 'coldpool_rasters');

rasterF = [fullfile(rasterfol, "SEBS_Ste_gear_temperature_" + string(srepyr') + ".tif") ...
           fullfile(rasterfol, "SEBS_Ste_model_bottom_temp_" + string(srepyr') + ".tif")];

% rasterF = GetFullPath(rasterF); % uncomment if rasterfol includes relative path, b/c readgeoraster is annoyingly picky
                                  % GetFullPath available from https://www.mathworks.com/matlabcentral/fileexchange/28249-getfullpath

%----------------------------
% Data that users need to
% download separately
%----------------------------

% Pelland Data

pellandfol = '~/Documents/Data/PMEL_Pelland_Dataset/'; % Pelland EcoFOCI CTD dataset
pellandver = 'v0.97/';

%% Read validation data (observational datasets)

%% ... popup float (PUF) data:
%      
%       PUFData: cell array of tables, each holding info from PUF csvs
%       Id:      table with PUF coords, ID, and max depth

F = dir(fullfile(pufdatafol, '*.csv'));
nfile = length(F);

for ii = nfile:-1:1
    
    fname = fullfile(F(ii).folder, F(ii).name);
    
    if contains(F(ii).name, 'PopUP') % original IMEI-naming
        fid = fopen(fname);
        vname = fgetl(fid);
        fclose(fid);
        vname = regexp(vname, ',', 'split');
        Tmp = readtable(fname, 'headerlines', 2, 'readvariablenames', false);
        Tmp.Properties.VariableNames = vname;

        Tmp.time = datetime(Tmp.time, 'format', 'uuuu-MM-dd''T''HH:mm:ss''Z''');
    elseif contains(F(ii).name, 'Spring')
        fid = fopen(fname);
        vname = fgetl(fid);
        fclose(fid);
        vname = regexp(vname, ',', 'split');
        vname = strrep(vname, 'Temp_DegC_0', 'sr_temp'); % Standardize variables to older ones
        vname = strrep(vname, 'Temp_DegC_1', 'fr_temp');
        vname = strrep(vname, 'Pressure_Bar', 'pressure');

        Tmp = readtable(fname, 'headerlines', 2, 'readvariablenames', false);
        Tmp.Properties.VariableNames = vname;
        Tmp.time = datetime(Tmp.time, 'format', 'uuuu-MM-dd''T''HH:mm:ss''Z''');
        
    elseif contains(F(ii).name, '2023')

        fid = fopen(fname);
        vname = fgetl(fid);
        fclose(fid);

        vname = regexp(vname, ',', 'split');
        vname = strrep(vname, 'Temp_DegC_0', 'sr_temp'); % Standardize variables to older ones
        vname = strrep(vname, 'Temp_DegC_1', 'fr_temp');
        vname = strrep(vname, 'Pressure_Bar', 'pressure');
        vname = strrep(vname, 'time (UTC)', 'time');
        vname = strrep(vname, 'latitude (degrees_north)', 'latitude');
        vname = strrep(vname, 'longitude (degrees_east)', 'longitude');

        Tmp = readtable(fname, 'headerlines', 1, 'readvariablenames', false);
        Tmp.Properties.VariableNames = vname;
        Tmp.time = datetime(Tmp.time, 'format', 'uuuu-MM-dd''T''HH:mm:ss''Z''');

    else
        fid = fopen(fname);
        vname = fgetl(fid);
        fclose(fid);
        vname = regexp(vname, ',', 'split');
        vname = strrep(vname, 'Temp_DegC_0', 'sr_temp'); % Standardize variables to older ones
        vname = strrep(vname, 'Temp_DegC_1', 'fr_temp');
        vname = strrep(vname, 'Pressure_Bar', 'pressure');
        vname = strrep(vname, 'time (UTC)', 'time');
        
        Tmp = readtable(fname, 'headerlines', 2, 'readvariablenames', false, 'delimiter', ',');
%         tmp = fileread(fname);
%         tmp = splitlines(tmp);
%         tmp = regexp(tmp, ',', 'split');
        
        Tmp.Properties.VariableNames = vname;
        Tmp.time = datetime(Tmp.time, 'format', 'uuuu-MM-dd HH:mm:ss+00:00');
        
    end
    Tmp.time.Format = 'uuuu-MM-dd HH:mm:ss';
    PUFData{ii} = Tmp;
    
end

% Remove pre-at-depth data

isgood = true(size(PUFData));

for ii = 1:nfile
    [tf,s1] = ischange(PUFData{ii}.pressure);
    idx = find(tf);
    if s1(idx) > s1(idx-1)
        PUFData{ii} = PUFData{ii}(idx:end,:);
    end   
end

% ID info

Id = cellfun(@(X) X(1,{'latitude', 'longitude', 'timeseries_id'}), PUFData, 'uni', 0);
Id = cat(1, Id{:});
Id.timeseries_id = string(Id.timeseries_id);

% Toss obviously bad ones

badids = ["2010" "3610" "4860" "POPS0006"];
tf = ismember(Id.timeseries_id, badids);
Id = Id(~tf,:);
PUFData = PUFData(~tf);

nfloat = length(PUFData);

% Find assumed depth 

maxpres = cellfun(@(X) max(X.pressure), PUFData');
Id.maxdepth = -gsw_z_from_p(maxpres*10, Id.latitude);

%% ... PMEL Moorings, full depth profiles

%% ...... Long-term

mooringfiles = ff(dir(fullfile(mooringfol, 'bering/*.nc')));

for im = 1:length(mooringfiles)
    
    % Read data

    Tmp = ncstruct(mooringfiles{im});

    % Add NaNs to mark big time gaps (for plotting simplicity)

    ttmp = datetime(1970,1,1)+seconds(Tmp.time);
    gapidx = find(diff(ttmp) > days(1));
    tnew = sort([ttmp; mean(ttmp([gapidx gapidx+1]),2)]);
    [~,loc] = ismember(ttmp, tnew);
    nt = length(tnew);
    nz = length(Tmp.depth);

    Mprof(im) = struct('depth', Tmp.depth, ...
                       't', tnew, ...
                       'temperature', nan(nt,nz), ...
                       'latitude',    nan(nt,nz), ...
                       'longitude',   nan(nt,nz), ...
                       'stemp',       nan(nt,1), ...
                       'btemp',       nan(nt,1));
    Mprof(im).temperature(loc,:) = Tmp.temperature;
    Mprof(im).latitude(   loc,:) = Tmp.latitude;
    Mprof(im).longitude(  loc,:) = Tmp.longitude;

    % Surface-most and bottom-most temperatures

    Mprof(im).btemp = bottom(permute(Mprof(im).temperature, [3 1 2]));
    Mprof(im).stemp = bottom(permute(Mprof(im).temperature(:,end:-1:1), [3 1 2]));

end
Mprof = Mprof(:);

%% ...... Kuskokwim

mooringfiles = ff(dir(fullfile(mooringfol, 'kuskokwim/*.nc')));

for im = 1:length(mooringfiles)
    
    % Read data

    Tmp = ncstruct(mooringfiles{im});

    % Add NaNs to mark big time gaps (for plotting simplicity)

    ttmp = datetime(1970,1,1)+seconds(Tmp.time);
    gapidx = find(diff(ttmp) > days(1));
    tnew = sort([ttmp; mean(ttmp([gapidx gapidx+1]),2)]);
    [~,loc] = ismember(ttmp, tnew);
    nt = length(tnew);
    nz = length(Tmp.depth);

    Kprof(im) = struct('depth', Tmp.depth, ...
                       't', tnew, ...
                       'temperature', nan(nt,nz), ...
                       'latitude',    nan(nt,nz), ...
                       'longitude',   nan(nt,nz), ...
                       'stemp',       nan(nt,1), ...
                       'btemp',       nan(nt,1));
    Kprof(im).temperature(loc,:) = Tmp.temperature;
    Kprof(im).latitude(   loc,:) = Tmp.latitude;
    Kprof(im).longitude(  loc,:) = Tmp.longitude;

    % Surface-most and bottom-most temperatures

    Kprof(im).btemp = bottom(permute(Kprof(im).temperature, [3 1 2]));
    Kprof(im).stemp = bottom(permute(Kprof(im).temperature(:,end:-1:1), [3 1 2]));

end
Kprof = Kprof(:);

%% ...... Inner front

mooringfiles = ff(dir(fullfile(mooringfol, 'innerfront/*.nc')));

for im = 1:length(mooringfiles)
    
    % Read data

    Tmp = ncstruct(mooringfiles{im});

    % Add NaNs to mark big time gaps (for plotting simplicity)

    ttmp = datetime(1970,1,1)+seconds(Tmp.time);
    gapidx = find(diff(ttmp) > days(1));
    tnew = sort([ttmp; mean(ttmp([gapidx gapidx+1]),2)]);
    [~,loc] = ismember(ttmp, tnew);
    nt = length(tnew);
    nz = length(Tmp.depth);

    Iprof(im) = struct('depth', Tmp.depth, ...
                       't', tnew, ...
                       'temperature', nan(nt,nz), ...
                       'latitude',    nan(nt,nz), ...
                       'longitude',   nan(nt,nz), ...
                       'stemp',       nan(nt,1), ...
                       'btemp',       nan(nt,1));
    Iprof(im).temperature(loc,:) = Tmp.temperature;
    Iprof(im).latitude(   loc,:) = Tmp.latitude;
    Iprof(im).longitude(  loc,:) = Tmp.longitude;

    % Surface-most and bottom-most temperatures

    Iprof(im).btemp = bottom(permute(Iprof(im).temperature, [3 1 2]));
    Iprof(im).stemp = bottom(permute(Iprof(im).temperature(:,end:-1:1), [3 1 2]));

end
Iprof = Iprof(:);

%% ... PUF/mooring analysis "location" clustering
%      Group together any that are very close to each other (redeployments
%      of PUFs or PUF deployment alongside mooring)

lt = [Id.latitude; arrayfun(@(x) nanmean(x.latitude(:)), [Mprof; Kprof])];
ln = [Id.longitude; arrayfun(@(x) nanmean(wrapTo360(x.longitude(:))), [Mprof; Kprof])];

locidx = clusterdata([lt ln], ...
    'criterion' ,'distance', 'cutoff', 20, 'distance', dfun, ...
    'linkage', 'single');

Loc = table(...
    splitapply(@mean, lt, locidx), ...
    splitapply(@mean, wrapTo360(ln), locidx), ...
    'variablenames', {'lat', 'lon'});

nloc = height(Loc);

%% ... Pelland EcoFOCI data (v0.93)

% ctdfile = ff(dir(fullfile(pellandfol, pellandver, '*ctd_data_top_250*.nc')));
% ctdfile = ctdfile{1};
% 
% Pel = ncstruct(ctdfile);
% Pel.TEMP(Pel.TEMP>100) = NaN; % a few errant 1e35s instead of missing
% 
% Pel.btemp = bottom(permute(Pel.TEMP, [2 3 1]));
% Pel.t = ncdateread(ctdfile, 'PROF_TDAY');
% Pel = struct2table(rmfield(Pel, {'TEMP', 'SAL', 'SAL_QC', 'AX_P'}));
% Pel = Pel(Pel.PROF_REGION == 2,:);
% 
% dctd = pdist2([Pel.PROF_LAT Pel.PROF_LON], [Loc.lat Loc.lon], dfun);

%% ... Pelland EcoFOCI data (v0.97+)

ctdfile = ff(dir(fullfile(pellandfol, pellandver, '*CTD-250-CDF*.nc')));
ctdfile = ctdfile{1};

Pel = ncstruct(ctdfile);

[Pel.btemp, ibot] = bottom(permute(Pel.Temperature, [1 3 2]));
Pel.t = ncdateread(ctdfile, 'TIME');
Pel.deepest_pressure = Pel.PRESSURE(ibot);

pelpres = Pel.PRESSURE;
Pel = struct2table(rmfield(Pel, 'PRESSURE'));

Pel.REGION = strtrim(string(Pel.REGION));
Pel = Pel(Pel.REGION=="Bering Sea" & abs(Pel.deepest_pressure - Pel.BOTTOM_DEPTH)<15,:);

dctd = pdist2([Pel.LATITUDE Pel.LONGITUDE], [Loc.lat Loc.lon], dfun);

%% Read model data

%% ... B10K grid data

Grd = ncstruct(grdfile);

%% ...... Match PUF/mooring locations to nearest model grid cell

[distToLoc,gidx] = pdist2([Grd.lat_rho(:) Grd.lon_rho(:)], [Loc.lat Loc.lon], dfun, 'smallest', 1);
Loc.gidx = gidx';

%% ... Daily surface and bottom data from B10K physics-only daily hindcast
%      We only read the daily data from the grid cells in the
%      vicinity of the PUF and mooring locations.

ncell = 1; % grid cell buffer per side

if reloadb10kdata

    blockyr = (1980:5:2020)';

    files = [fullfile(ppsimdir, 'Level2', "B10K-K20nobio_CORECFS_daily_"+blockyr+"-"+(blockyr+4)+"_average_temp_bottom5m.nc"), ...
             fullfile(ppsimdir, 'Level2', "B10K-K20nobio_CORECFS_daily_"+blockyr+"-"+(blockyr+4)+"_average_temp_surface5m.nc")];
    
    Sim = struct;
    Sim.t = ncdateread(files(:,1), 'ocean_time');
    nt = length(Sim.t);
    
    Sim.btemp = nan(nt, (ncell*2+1).^2, nloc);
    Sim.stemp = nan(nt, (ncell*2+1).^2, nloc);
    
    [xi, eta] = ind2sub(size(Grd.h), Loc.gidx);
    for ii = 1:nloc
        fprintf('Reading model data, %d/%d\n', ii, nloc);
        scs = struct('xi_rho', [xi(ii)-ncell ncell*2+1 1], 'eta_rho', [eta(ii)-ncell ncell*2+1 1]);
        Stmp = ncstruct(files(:,1), 'temp', scs);
        Sim.btemp(:,:,ii) = reshape(Stmp.temp, (ncell*2+1).^2, nt)';    
        Stmp = ncstruct(files(:,2), 'temp', scs);
        Sim.stemp(:,:,ii) = reshape(Stmp.temp, (ncell*2+1).^2, nt)';   
    end
    save(b10ksimdatafile, '-struct', 'Sim');
else
    Sim = load(b10ksimdatafile);
end

%% ... Groundfish trawl and Survey-replicated model 

if rerunsrepflag && ~exist(srepfile, 'file')

    mask1 = ncread(grdfile, 'mask_rho');

    blockyr = (1980:5:2020)';

    files = fullfile(ppsimdir, 'Level2', "B10K-K20nobio_CORECFS_daily_"+blockyr+"-"+(blockyr+4)+"_average_temp_bottom5m.nc");
    
    coldpool_surveyreplicate('surveyfile', '', 'outputfile', srepfile, ...
        'label', srepsimname, ...
        'tempfiles', files, ...
        'vname', 'temp', ...
        'tname', 'ocean_time', ...
        'ltname', 'lat_rho', ...
        'mask', mask1);

end

Srep = readtable(srepfile);

[g, stationid] = findgroups(Srep.stationid);
S = splitapply(@(a,b) skillstats(a,b), Srep.gear_temperature, Srep.(srepsimname), g);

TrawlMeta = table(stationid, ...
                  splitapply(@mean, Srep.latitude, g), ...
                  splitapply(@mean, Srep.longitude, g), ...
                  splitapply(@std, doy(Srep.start_time), g), ...
                  splitapply(@range, doy(Srep.start_time), g), ...
                  'variablenames', {'id', 'lat', 'lon', 'doystd', 'doyrange'});

% Extract recent trawl data for ease of access in later plots

% isin = ismember(Srep.year, 2020:2022);
% Trawl = Srep(isin,:);
Trawl = Srep;

dtrawl = pdist2([Trawl.latitude Trawl.longitude], [Loc.lat Loc.lon], dfun);

% Variation in sampling date





%% ... interpolated rasters

[A,R] = readgeoraster(rasterF{1});
[x,y] = R.worldGrid;

Raster = struct;
[Raster.lt, Raster.ln] = projinv(R.ProjectedCRS, x, y);

Raster.survey = arrayfun(@readgeoraster, rasterF(:,1), 'uni', 0);
Raster.model  = arrayfun(@readgeoraster, rasterF(:,2), 'uni', 0);

Raster.survey = cat(3, Raster.survey{:});
Raster.model = cat(3, Raster.model{:});

hasdata = ~all(isnan(Raster.survey),3);
imask = any(hasdata,2);
jmask = any(hasdata,1);

Raster = structfun(@(x) x(imask,jmask,:), Raster, 'uni', 0);

w = warning('off');
for ii = size(Raster.survey,1):-1:1
    for jj = size(Raster.survey,2):-1:1
        Raster.Skill(ii,jj) = skillstats(...
            squeeze(Raster.survey(ii,jj,:)), ...
            squeeze(Raster.model(ii,jj,:)));
    end
end
warning(w);

extractfun = @(fld) arrayfun(@(x) x.(fld)(2), Raster.Skill);

%% ... ETOPO high-res bathymetry (here because extent is model-dependent)

% etop = '/Volumes/LaCie2023/ETOPO1/ETOPO1_Ice_g_gmt4.grd';
% E = ncstruct(etop, 'x', 'y');
% 
% latlim = minmax(Loc.lat, 'expand');
% lonlim = minmax(Loc.lon, 'expand');
% % latlim = minmax(Grd.lat_rho);
% 
% ymask = E.y >= latlim(1) & E.y <= latlim(2);
% 
% xmask = wrapTo360(E.x) >= lonlim(1) & wrapTo360(E.x) <= lonlim(2);
% scs = struct('x', [find(xmask,1) nnz(xmask) 1], ...
%              'y', [find(ymask,1) nnz(ymask) 1]);
% 
% E1 = ncstruct(etop, scs);

% ETOPO_2022 data

latlimbig = minmax(Grd.lat_rho);
lonlimbig = minmax(Grd.lon_rho);

E = ncstruct(etopotilename, 'lat', 'lon');
latmask = E.lat >= latlimbig(1) & E.lat <= latlimbig(2);
lonmask = wrapTo360(E.lon) >= lonlimbig(1) & wrapTo360(E.lon) <= lonlimbig(2);

latscs = [find(latmask, 1, 'first') nnz(latmask) 1];
[b, n, bi] = RunLength(lonmask);
lonscs = [bi(b==1) n(b==1) ones(nnz(b),1)];

Etop(1) = ncstruct(etopotilename, struct('lat', latscs, 'lon', lonscs(1,:)));
Etop(2) = ncstruct(etopotilename, struct('lat', latscs, 'lon', lonscs(2,:)));

Etop(1).lon = wrapTo360(cat(1, Etop([2 1]).lon));
Etop(1).z   = cat(1, Etop([2 1]).z);
Etop = Etop(1);

%% ********* PLOTS ***********

exportfigs = false;
return


%% ... color setup

% Color brewer qualitative set 2, 8 colors
Col = struct;
% Col.b10k  = 'ffd92f'; % gold
% Col.puf   = 'fc8d62'; % orange
% Col.moor  = 'e78ac3'; % pink
% Col.gfish = '8da0cb'; % blue
% Col.kusk  = 'e78ac3'; % tan
% Col.ctd   = 'a6d854'; % lime green
% Col.inf   = 'e78ac3'; % teal
% Col.crab  = '66c2a5'; % teal
% 
% Col = structfun(@hex2rgb, Col, 'uni', 0);

% Option 2: tsitsul.in "Normal" palette

Col.b10k  = [ 64  83 211]./255; % blue
Col.puf   = [251  73 176]./255; % pink
Col.moor  = [189  29  20]./255; % red
Col.gfish = [221 179  16]./255; % gold
Col.kusk  = Col.moor;
Col.ctd   = [  0 178  93]./255; % green  
Col.inf   = Col.moor;
Col.crab  = [  0 190 255]./255; % light blue
Col.gray   = [1 1 1]*202/255;


Col.connector = Col.b10k; %rgb('red');

adjustcol = @(c,x) interp1([-1 0 1]', [0 0 0; c; 1 1 1], x);

%% ... Map limits for most maps

latlim = minmax(Loc.lat, 'expand');
lonlim = minmax(Loc.lon, 'expand');

rlatlim = minmax(Raster.lt);
rlonlim = minmax(wrapTo360(Raster.ln));

%% ... Borders

[blat{1}, blon{1}] = borders('alaska');
[blat{2}, blon{2}] = borders('russia');

blon = cellfun(@wrapTo360, blon, 'uni', 0);

% Polygon (for shaded)

Border.poly = polyshape(blon, blat);

% Lines (for outline w/o 180 line)

is180 = abs(blon{2} - 180) < 0.01;
[b,n,bi] = RunLength(is180);
n = n(b);
bi = bi(b);
for ii = 1:length(bi)
    if n(ii)>2
        blon{2}(bi(ii)+(1:(n(ii)-2))) = NaN;
        blat{2}(bi(ii)+(1:(n(ii)-2))) = NaN;
    end
end
needextra = bi(n==2);
nper = diff([0 needextra length(blat{2})]);
blat(3+(0:length(needextra))) = mat2cell(blat{2}, 1, nper);
blon(3+(0:length(needextra))) = mat2cell(blon{2}, 1, nper);

Border.lon = [blon{:}];
Border.lat = [blat{:}];

%% Timeseries plot arranged spatially-ish with map in center

%% ... some setup parameters

nr = 8; % Number of rows in outer axis ring
myproj = 'mercator';

%% ... Axes setup (including order to minimize connector crossings)

% Create axes: map surrounded by ring of axes (in theta order)

nc = ceil((nloc + 4 -2*nr)/2);

h = plotgrid('size', [nr nc], 'sp', 0.01, 'mar', 0.04);
bij = {floor(nr/2):-1:1      nc
      1                     (nc-1):-1:1
      2:nr                  1
      nr                    2:(nc-1)
      nr:-1:(floor(nr/2)+1) nc};
for ii = 1:size(bij,1)
    if isscalar(bij{ii,1})
        bij{ii,1} = ones(size(bij{ii,2}))*bij{ii,1};
    end
    if isscalar(bij{ii,2})
        bij{ii,2} = ones(size(bij{ii,1}))*bij{ii,2};
    end
end
bij = [cat(2, bij{:,1}); cat(2, bij{:,2})];
bidx = sub2ind([nr nc], bij(1,:), bij(2,:));

anchorpoint = string;
anchorpoint(nr,nc) = "lt";
anchorpoint(nr, 1) = "rt";
anchorpoint( 1, 1) = "rb";
anchorpoint( 1,nc) = "lb";
anchorpoint(2:nr-1, 1) = "rc";
anchorpoint(2:nr-1,nc) = "lc";
anchorpoint( 1,2:nc-1) = "cb";
anchorpoint(nr,2:nc-1) = "ct";

h.bax = h.ax(bidx);
anchorpoint = anchorpoint(bidx);

% h.bax = [h.ax(floor(nr/2):-1:1,end)' h.ax(1,end-1:-1:1) h.ax(2:end,1)' h.ax(end,2:end-1) h.ax(end:-1:(floor(nr/2)+1),end)'];
% [baxrow, baxcol] = ndgrid(1:nr,1:nc);

delete(h.bax(nloc+1:end));
h.bax = h.bax(1:nloc);cd
bij = bij(:,1:nloc);

h.in = mergeaxes(h.ax(2:end-1,2:end-1));

h.fig.Position(3:4) = [797 529]; % (11in x 7.35in)
set(h.fig, 'color', 'w');

% Data point coordinates in figure units

axes(h.in);
worldmap(latlim, lonlim);
setm(h.in, 'mapprojection', myproj);
tightmap;
[xtmp, ytmp] = projfwd(getm(gca), Loc.lat, Loc.lon);
[x1, y1] = axescoord2figurecoord(xtmp, ytmp, h.in);

% Axis anchor point coordinates in figure units

pos = get(h.bax, 'Position');
pos = cat(1, pos{:}) .* [h.fig.Position([3 4 3 4])];

x2l = pos(:,1);
x2c = pos(:,1) + pos(:,3)./2;
x2r = pos(:,1) + pos(:,3);

y2b = pos(:,2);
y2c = pos(:,2) + pos(:,4)./2;
y2t = pos(:,2) + pos(:,4);

x2 = x2c;
y2 = y2c;
x2(startsWith(anchorpoint, 'l')) = x2l(startsWith(anchorpoint, 'l'));
x2(startsWith(anchorpoint, 'r')) = x2r(startsWith(anchorpoint, 'r'));
y2(endsWith(anchorpoint, 'b')) = y2b(endsWith(anchorpoint, 'b'));
y2(endsWith(anchorpoint, 't')) = y2t(endsWith(anchorpoint, 't'));

% Bipartite graph matching (fun new fact: minimizing total edge length
% results in no edge crossings, because math apparently does cool things
% sometimes.  :-)  Math terms to remember, if I ever want to remember what
% this is doing: euclidean bipartite matching problem, Hungarian algorthm)

costmat = pdist2([x1 y1], [x2 y2]);
[pairings, cost] = munkres(costmat);
[i1,i2] = find(pairings);

isrt = i1;

xc = [x1(isrt)'; x2'; nan(1,nloc)]; % Note; in figure coordinates
yc = [y1(isrt)'; y2'; nan(1,nloc)];

h.axan = axes('position', [0 0 1 1]);

h.cnct = patch(h.axan, xc, yc, 'r');
set(h.cnct, 'edgecolor', Col.connector, 'edgealpha', 0.2, 'linewidth', 4);
set(h.axan, 'xlim', [0 h.fig.Position(3)], 'ylim', [0 h.fig.Position(4)], 'visible', 'off');

% plot(h.axan, xc, yc, ':r', 'linewidth', 4);
% set(h.axan, 'xlim', [0 1], 'ylim',[0 1], 'visible', 'off');

%% ... Organize indices of which locations correspond to which axes

% Match up float and mooring data with primary locations

moorlat = arrayfun(@(x) nanmean(x.latitude(:)), Mprof);
moorlon = arrayfun(@(x) nanmean(wrapTo360(x.longitude(:))), Mprof);
kusklat = arrayfun(@(x) nanmean(x.latitude(:)), Kprof);
kusklon = arrayfun(@(x) nanmean(wrapTo360(x.longitude(:))), Kprof);
iflat = arrayfun(@(x) nanmean(x.latitude(:)), Iprof);
iflon = arrayfun(@(x) nanmean(wrapTo360(x.longitude(:))), Iprof);

[~,pufidx] = pdist2([Loc.lat Loc.lon], [Id.latitude Id.longitude], dfun, 'smallest', 1);
[~,moridx] = pdist2([Loc.lat Loc.lon], [moorlat moorlon], dfun, 'smallest', 1);
[~,kusidx] = pdist2([Loc.lat Loc.lon], [kusklat kusklon], dfun, 'smallest', 1);
[~,infidx] = pdist2([Loc.lat Loc.lon], [iflat iflon], dfun, 'smallest', 1);

% Table o' indices (so I can keep track of all this!)

lbls = string(char((1:length(bidx))+96)');
[~,lsrt] = ismember(bidx, sort(bidx));

Indices = struct('clockwise', (1:nloc)', ...
                'loc', isrt, ...
                'letters', lbls(lsrt));

Indices.moor = arrayfun(@(x) find(moridx==x), Indices.loc, 'uni', 0);
Indices.kusk = arrayfun(@(x) find(kusidx==x), Indices.loc, 'uni', 0);
Indices.puf = arrayfun(@(x) find(pufidx==x), Indices.loc, 'uni', 0);
Indices.inf = arrayfun(@(x) find(infidx==x), Indices.loc, 'uni', 0);

Indices = struct2table(Indices);

%% ... Plot timeseries (outer axes)

dmin = nm2km(sqrt(20^2 + 20^2)/2); % radius to match survey resolution

for ii = 1:nloc
    
    ip = isrt(ii);
    
    % B10K hindcast
    
    h.ts.b10k(ii,:) = plot(h.bax(ii), Sim.t, Sim.btemp(:,:,ip), 'color', adjustcol(Col.b10k, 0.5));
    set(h.ts.b10k(ii,5), 'color', Col.b10k);
    uistack(h.ts.b10k(ii,5),'top');
    
    % PUF
    
    hold(h.bax(ii),'on');
    for iplt = find(pufidx == ip)
        h.ts.puf(iplt) = plot(h.bax(ii), PUFData{iplt}.time, PUFData{iplt}.sr_temp, 'color', Col.puf);
    end
    
    % Moorings
    
    for iplt = find(moridx == ip)
        % h.moor(iplt) = plot(h.bax(ii), Moor.data{iplt}.time, Moor.data{iplt}.temperature, 'color', Col.moor);
        h.ts.moor(iplt) = plot(h.bax(ii), Mprof(iplt).t, Mprof(iplt).btemp, 'color', Col.moor);
    end
    
    % Groundfish trawl (recent years only)
    
    isclose = dtrawl(:,ip) < dmin;
    if any(isclose)
        h.ts.trawl(ii) = scatter(h.bax(ii), Trawl.start_time(isclose), Trawl.gear_temperature(isclose), dmin-dtrawl(isclose,ip), Col.gfish, 'filled');
    end
    
    % Groundfish trawl and survey-rep model
    
    % if any(isclose)
    %     sid = unique(Trawl.stationid(isclose));
    %     mask = ismember(Srep.stationid, sid);
    %     h.ts.trawl2(ii) = plot(h.bax(ii), Srep.start_time(mask), Srep.gear_temperature(mask), 'o');
    %     h.ts.srep(ii)   = plot(h.bax(ii), Srep.start_time(mask), Srep.(srepsimname)(mask), 'o');
    % end

    % EcoFOCI CTD data

    isclose = dctd(:,ip) < dmin;
    % isrecent = Pel.t >= datetime(2020,1,1);
    isrecent = true(size(Pel.t));
    if any(isclose)
        h.ts.ctd1(ii) = scatter(h.bax(ii), Pel.t(isclose&isrecent), Pel.btemp(isclose&isrecent), dmin-dctd(isclose&isrecent,ip), 'y', 'filled');
        % h.ts.ctd2(ii) = plot(h.bax(ii), Pel.t(isclose&~isrecent), Pel.btemp(isclose&~isrecent), 'o');
    end

    % Kuskokwim buoys

    for iplt = find(kusidx == ip)
        % h.kusk(iplt) = plot(h.bax(ii), Kusk{iplt}.t, Kusk{iplt}.temperature, 'color', Col.kusk);
        h.ts.kusk(iplt) = plot(h.bax(ii), Kprof(iplt).t, Kprof(iplt).btemp, 'color', Col.kusk);
    end

    % Inner front buoys

    % for iplt = find(infidx == ip)
    %     h.innerfront(iplt) = plot(h.bax(ii), IF{iplt}.t, IF{iplt}.temperature, 'color', Col.inf);
    % end

end

set(h.ts.trawl, 'markeredgecolor', adjustcol(Col.gfish, -0.3), 'markerfacecolor', Col.gfish);
set(h.ts.ctd1, 'markeredgecolor', adjustcol(Col.ctd, -0.3), 'markerfacecolor', Col.ctd, 'marker', 's');

% set(h.ts.trawl2, 'markerfacecolor', Col.gfish, 'markeredgecolor', 'none', 'markersize', 2);
% set(h.ts.srep,   'markerfacecolor', Col.b10k, 'markeredgecolor', 'none', 'markersize', 2);
% set(h.ts.ctd2, 'markerfacecolor', Col.ctd, 'markeredgecolor', 'none', 'markersize', 2);

% Axis tick labels

set(h.bax, 'ylim', [-2 12]);

Tk(1).tk = datetime(2020,1:3:48,1);
Tk(1).lim = datetime([2020 2024], 1, 1);
Tk(1).lbl = string(datestr(Tk(1).tk,'m'));
Tk(1).lbl(1:4:end) = string(datestr(Tk(1).tk(1:4:end), 'Jyy'));

Tk(2).tk = datetime(1982:5:2023,1,1);
Tk(2).lim = datetime([1982 2023],1,1);
Tk(2).lbl = datestr(datetime(1982:5:2023,1,1), 'yy');

set(h.bax, 'tickdir', 'out', 'fontsize', 6, 'tickdir', 'out');

set(h.bax(bij(1,:) == 1), 'xaxisloc', 'top');
set(h.bax(bij(2,:) == nc), 'yaxisloc', 'right');
set(h.bax(~ismember(bij(2,:), [1 nc])), 'yticklabel', '');
set(h.bax(~ismember(bij(1,:), [1 nr])), 'xticklabel', '');

% set(h.bax(baxnoxlbl), 'xticklabel', '');
% set(h.bax(baxyright), 'yaxisloc', 'right');

% Labels

ylabel(h.bax(bidx==1), 'Bottom temperature (\circC)')

labelaxes(h.bax, Indices.letters, 'northwest');

%% ... Plot map (inner axis)

latlim = minmax(Loc.lat, 'expand');
lonlim = minmax(Loc.lon, 'expand');

% Map setup and bathymetry

axes(h.in);
worldmap(latlim, lonlim);
% h.ras.bathy = pcolorm(E.y, E.x, E.z');
h.ras.bathy = pcolorm(Etop.lat, Etop.lon, Etop.z');
bordersm('alaska', 'k')
setm(h.in, 'flinewidth', 1, 'fontsize', 7, 'mapprojection', myproj, 'mlabelparallel', 54.8, 'plabelmeridian', -178.5, 'fontcolor', rgb('light gray'));
tightmap;
set(h.in, 'clim', [-200 0]);
% cmocean('gray');
cmocean('-turbid');

% Raster skill layers

h.ras.bias = pcolorm(Raster.lt, Raster.ln, arrayfun(@(X) X.bias(2),Raster.Skill));
h.ras.mef  = pcolorm(Raster.lt, Raster.ln, arrayfun(@(X) X.mef(2),Raster.Skill));

h.svy.bias = scatterm(TrawlMeta.lat, TrawlMeta.lon, 10, arrayfun(@(x) x.bias(2), S), 'filled', 'marker', 's', 'markeredgecolor', 'k');
h.svy.mef  = scatterm(TrawlMeta.lat, TrawlMeta.lon, 10, arrayfun(@(x) x.mef(2), S), 'filled', 'marker', 's', 'markeredgecolor', 'k');

% Groundfish trawl

h.loc.trawl = plotm(Trawl.latitude, Trawl.longitude, '.', 'color', Col.gfish);
set(h.loc.trawl, 'marker', 'o', 'markersize', 2, 'markerfacecolor', Col.gfish, 'markeredgecolor', adjustcol(Col.gfish, -0.3));

% PMEL CTDs

ctdmask = any(dctd<dmin,2);
% h.loc.ctd = plotm(Pel.PROF_LAT(ctdmask), Pel.PROF_LON(ctdmask), '.', 'color', Col.ctd);
h.loc.ctd = plotm(Pel.LATITUDE(ctdmask), Pel.LONGITUDE(ctdmask), 's', 'color', Col.ctd);
set(h.loc.ctd, 'marker', 's', 'markersize', 2, 'markerfacecolor', Col.ctd, 'markeredgecolor', adjustcol(Col.ctd, -0.3));


% PUFs

h.loc.puf = plotm(Id.latitude, Id.longitude, 'o');
set(h.loc.puf, 'markerfacecolor', Col.puf, 'color', adjustcol(Col.puf, -0.5));

% Moorings (long-term, kuskokwim, inner front)

h.loc.moor = plotm(moorlat, moorlon, '^');
set(h.loc.moor, 'markerfacecolor', Col.moor, 'color', adjustcol(Col.moor, -0.5), 'markersize', 10);

h.loc.kusk = plotm(kusklat, kusklon, 'v');
set(h.loc.kusk, 'markerfacecolor', Col.kusk, 'color', adjustcol(Col.kusk, -0.5), 'markersize', 5)

h.loc.if = plotm(iflat, iflon, 'v');
set(h.loc.if, 'markerfacecolor', Col.inf, 'color', adjustcol(Col.inf, -0.5), 'markersize', 5)

% % Crab mrPTs
% 
% for ic = 1:ncrab
%     if hasloc(ic) && hasdd(ic)
%         h.loc.crab(ic) = plotm(crabdaily{ic}.latestimate, crabdaily{ic}.lonestimate, 'marker', '+', 'color', Col.crab, 'markersize', 2);
%     end
%     % h.crab2 = plotm(crabdata.Latitude, crabdata.Longitude, 'marker', '+', 'color', Col.crab);
% end

% B10K grid cells

[xi, eta] = ind2sub(size(Grd.h), Loc.gidx);
[xbuf,ebuf] = ndgrid(-ncell:ncell, -ncell:ncell);
xall = xi +permute(xbuf, [3 1 2]);
eall = eta+permute(ebuf, [3 1 2]);
gall = sub2ind(size(Grd.h), xall, eall);

xpsi = xall(:) + [0 0 1 1 0] - 1;
epsi = eall(:) + [0 1 1 0 0] - 1;
gpsi = sub2ind(size(Grd.mask_psi), xpsi, epsi);
h.loc.b10k = plotm(Grd.lat_psi(gpsi'), Grd.lon_psi(gpsi'), 'color', Col.b10k);

% h.loc.b10k = plotm(Grd.lat_rho(gall(:)), Grd.lon_rho(gall(:)), 's');
% set(h.loc.b10k, 'color', Col.b10k, 'markerfacecolor', Col.b10k, 'markersize', 3);

% Add B10K bathymetry contours and ETOPO bathymetry colorbar

contourm(Grd.lat_rho, Grd.lon_rho, Grd.h, [50 100 200], 'w');

h.cb = colorbar('north');
pos = plotboxpos(h.in);
set(h.cb, 'position', [pos(1)+pos(3)-0.11 pos(2)+pos(4)-h.cb.Position(4)/2 0.1 h.cb.Position(4)/2]);
tk = -200:50:0;
set(h.cb, 'tickdir', 'out', 'fontsize', 6, 'ticks', tk, 'ticklabels', compose("%dm", -tk));

% Stacking

uistack(h.loc.puf, 'top');

uistack(h.axan, 'top');
arrayfun(@(x) uistack(x,'top'), h.bax);

%% ... Saveable versions of figure

% saveflag = false; 
plotopt = {'bathy'};
% plotopt = {'bathy', 'bias', 'mef', 'biasfull', 'meffull'};

for ii = 1:length(plotopt)

    structfun(@(x) set(x, 'visible', 'off'), h.ras);
    structfun(@(x) set(x, 'visible', 'off'), h.svy);

    % Switch focus of central map

    switch plotopt{ii}

        case {'bathy'}
            set(h.ras.bathy, 'visible', 'on');
            set(h.loc.trawl, 'visible', 'on');

            cmaptmp = linspace(0.4,1,22)'.*ones(1,3);
            set(h.in, 'clim', [-220 0], 'colormap', cmaptmp);
            cbtk = -200:50:0;
            set(h.cb, 'ticks', cbtk, 'ticklabels', compose("%dm", -cbtk));
            xlabel(h.cb, '');
            
            set(h.cnct, 'edgecolor', Col.b10k);
            
            set(h.ts.b10k, 'color', Col.gray);
            set(h.ts.b10k(:,5), 'color', Col.b10k);
            for il = 1:size(h.ts.b10k,1)
                uistack(h.ts.b10k(il,[1:4 6:9]), 'bottom');
            end

        case {'bias', 'biasfull'}
            set([h.svy.bias h.ras.bias], 'visible', 'on');
            set(h.loc.trawl, 'visible', 'off');

            set(h.cnct, 'edgecolor', rgb('dark gray'));

            set(h.in, 'colormap', cmocean('balance'), 'clim', [-2 2]);
            set(h.cb, 'TickLabelsMode', 'auto', 'Ticks', -2:2);
            xlabel(h.cb, 'Bias (\circC)');

        case {'mef', 'meffull'}
            set([h.svy.mef  h.ras.mef], 'visible', 'on');
            set(h.loc.trawl, 'visible', 'off');

            set(h.cnct, 'edgecolor', rgb('dark gray'));

            set(h.in, 'colormap', cmocean('curl'), 'clim', [-1 1]);
            set(h.cb, 'TickLabelsMode', 'auto', 'Ticks', -1:1);
            xlabel(h.cb, 'Model Efficiency');
    end

    % Switch limits of timeseries axes

    switch plotopt{ii}
        case {'bathy', 'bias', 'mef'}
            set(h.bax, 'xlim', [datetime(2020,3,1) datetime(2023,10,1)], ...
                       'xtick', datetime(2020:2024,1,1), ...
                       'XGrid', 'on', ...
                       'XMinorGrid', 'on', ...
                       'XMinorTick', 'on', ...
                       'GridColor', 'k', ...
                       'MinorGridColor', rgb('tan'));
            cellfun(@(x) set(x, 'MinorTickValues', datetime(2020,1:3:(12*5),1)), ...
                    get(h.bax, 'XAxis'));

            set(h.bax(~ismember(bij(2,:), [1 nc])), 'yticklabel', '');
            set(h.bax(~ismember(bij(1,:), [1 nr])), 'xticklabel', '');
        case {'biasfull', 'meffull'}
            % set(h.bax, 'xlim', datetime([1982 2023], 1, 1), ...
            set(h.bax, 'xlim', datetime([2000 2023], 1, 1), ...
                       'xtick', datetime(1980:5:2025,1,1), ...
                       'XGrid', 'off', ...
                       'XMinorGrid', 'off', ...
                       'XMinorTick', 'off');

            % cellfun(@(x) set(x, 'TickLabelFormat', 'yy'), ...
            %         get(h.bax, 'XAxis'));
            set(h.bax(~ismember(bij(2,:), [1 nc])), 'yticklabel', '');
            set(h.bax(~ismember(bij(1,:), [1 nr])), 'xticklabel', '');
    end

    arrayfun(@(x) uistack(x,'top'), h.bax);

    if exportfigs
        fname = sprintf('popup_vs_model_with_map_%s_%s', datetime('today', 'format', 'yyyyMM'), plotopt{ii});
        export_fig(fname, h.fig, '-png', '-r300', '-nocrop');
    end
end

%% ... Supplement: outer axes over longer time period

% TODO: add legends

legobj = [h.ts.b10k(1,5) h.ts.b10k(1,1) h.ts.trawl(1), h.ts.ctd1(1) h.ts.moor(1) h.ts.puf(1)];
legstr = {'Bering10K (center grid cell)' 'Bering10K (other grid cells)' 'Grounfish trawl' 'CTD' 'Mooring' 'PUF'};

[~,isrt] = sortrows(Indices, 'letters');
for ii = 1:height(Indices)
    hsup.fig = figure('color', 'w');
    hsup.fig.Position(3:4) = [612 216]; % 8.5in x 3in

    hsup.ax = copyobj(h.bax(isrt(ii)), hsup.fig);
    set(hsup.ax, ...
        'position', [0.05 0.1 0.93 0.85], ...
        'xlim', datetime([1985 2024],1,1), ...
        'xtickmode', 'auto', ...
        'xticklabelmode', 'auto', ...
        'xgrid', 'off', ...
        'xminortick', 'off', ...
        'xminorgrid', 'off', ...
        'xaxisloc', 'bottom', ...
        'yaxisloc', 'left', ...
        'fontsize', 8, ...
        'yticklabelmode', 'auto');

    ylabel(hsup.ax, 'Bottom temperature (\circC)')

    
    leg = legendflex(legobj, legstr, 'ref', hsup.ax, 'anchor', {'n','n'}, ...
        'buffer', [0 -5], 'box', 'off', 'nrow', 1, 'xscale', 0.5, 'fontsize', 8);

    if exportfigs
        fname = sprintf('supplement_timeseries_%02d', ii);
        export_fig(fname, hsup.fig, '-png', '-r300', '-nocrop');
    end
    close(hsup.fig);

end

%% ... rough comparison: model vs PUF

delta_mod_puf = cell(size(h.bax));
for ii = 1:length(h.bax)
    hmod = findall(h.bax(ii), 'color', Col.b10k);
    hpuf = [findall(h.bax(ii), 'color', Col.puf) findall(h.bax(ii), 'color', Col.moor)];

    if ~isempty(hpuf)
        x = [hpuf.XData];
        y = [hpuf.YData];
        [x,iunq] = unique(x);
        y = y(iunq);

        test = synchronize(timetable(x', y'), ...
                           timetable(hmod.XData', hmod.YData'), ...
                           'first', 'linear');
        delta_mod_puf{ii} = diff(test, 1, 2);
    end

end

figure;
axes;
hold on;

lettermiddle = ["i" "k" "h" "s" "a" "b" "f" "u"];
[tf,loc] = ismember(lettermiddle, Indices.letters);

plotmask = cellfun(@isempty, delta_mod_puf);

for ii = loc
    if ~isempty(delta_mod_puf{ii})
        
        % plot(doy(delta_mod_puf{ii}.Time), delta_mod_puf{ii}.Var1_1, '.');
        plot(delta_mod_puf{ii}.Time, delta_mod_puf{ii}.Var1_1, '.');

    end
end
% set(gca, 'xtick', doy(datetime(2023,3:3:12,1)), 'xlim', [0 365], 'xgrid', 'on');
legend(lettermiddle(~plotmask(loc)));


%% Raster plot

h = plotgrid('size', [3 3], 'sp', 0, 'mar', 0.01);

h.fig.Position(3:4) = [792 612]; % 11in x 8.5in
set(h.fig, 'color', 'w');

extractfun = @(fld) arrayfun(@(x) x.(fld)(2), Raster.Skill);

fld = setdiff(fieldnames(Raster.Skill),'ri');
Rplt = cell2struct(cellfun(extractfun, fld, 'uni', 0), fld, 1);

Rplt.mean = mean(Raster.survey, 3);
Rplt.std = std(Raster.survey, 1, 3);

data = {...
    mean(Raster.survey, 3) std(Raster.survey, 1, 3) []
    extractfun('cor') extractfun('bias') extractfun('stdnorm')
    extractfun('crmsd') extractfun('mef') []
    };

lbl = {...
    'Survey mean (\circC)', 'Survey standard deviation (\circC)', []
       'Correlation', 'Bias', 'Normalized standard deviation'
       'Centered RMSD', 'Model efficiency', []
       };

lims = {...
    [-2 3], [0 2], []
    [-1 1], [-2 2], [-0.5 2.5]
    [0 2], [-1 1] []
    };

cmap = {...
    cmocean('-dense') cmocean('matter') []
    crameri('broc') cmocean('balance') crameri('cork')
    cmocean('matter') cmocean('curl') []
    };


for ii = 1:numel(data)
    if ~isempty(data{ii})
        axes(h.ax(ii));
        worldmap(rlatlim, rlonlim);
        h.pc(ii) = pcolorm(Raster.lt, Raster.ln, data{ii});
        h.cb(ii) = colorbar('south');
        h.cb(ii).Position(3) = 0.1;
        xlabel(h.cb(ii), lbl{ii});
        set(h.ax(ii), 'clim', lims{ii}, 'colormap', cmap{ii});
        bordersm('alaska', 'k');

        contourm(Grd.lat_rho, Grd.lon_rho, Grd.h, [50 100 200], 'w');

    end
end

isemp = cellfun(@isempty, data);
set(h.ax(isemp), 'visible', 'off');
arrayfun(@(x) setm(x, 'frame', 'off', 'parallellabel', 'off', 'meridianlabel', 'off'), h.ax(~isemp));

set(h.cb(~isemp), 'tickdir', 'out');
set(h.cb(8), 'Limits', [0 2.5])

lbl2 = repmat("", 3, 3);
lbl2(~isemp') = string(char(96+(1:nnz(~isemp)))') + ")";
lbl2 = lbl2';
h.panellbl = labelaxes(h.ax(~isemp), lbl2(~isemp), 'northwest');

if exportfigs
    export_fig('surveyrep_raster_stats', h.fig, '-png', '-r300', '-nocrop');
end

%% ... add labels for reference (for analysis only, not for final figure)

for ii = 1:numel(data)
    if ~isempty(data{ii})
        axes(h.ax(ii));
        textm(Loc.lat(Indices.loc), Loc.lon(Indices.loc), Indices.letters, 'horiz', 'center');

        contourm(Grd.lat_rho, Grd.lon_rho, Grd.h, [50 100 200], 'g');
        
        if ii == 1
            [ln,lt] = ndgrid(Etop.lon, Etop.lat);
            [x,y] = projfwd(getm(gca),lt,ln);
        end
        contour(x, y, Etop.z, -[50 100 200], 'm');

    end
end


%% Region label plot (note: regions added in BioRender)

h = plotgrid('size', [1 1], 'mar', 0.02, 'mb', 0.1);
set(h.fig, 'color', 'w');
boxworldmap(rlatlim, rlonlim, 'latgrid', 50:5:60, 'longrid', -180:5:-160);

map1 = cmocean('balance', 7);
map2 = cmocean('curl', 7);

[bx,by] = projfwd(getm(h.ax), Border.poly.Vertices(:,2), Border.poly.Vertices(:,1));
plot(polyshape(bx,by), 'facecolor', rgb('gray'), 'edgecolor', 'none');
% bordersm('alaska', 'facecolor', rgb('light gray'), 'edgecolor', 'none');
contourm(Grd.lat_rho, Grd.lon_rho, Grd.h, [50 100 200], 'color', rgb('gray'), 'linewidth', 1.5, 'linestyle', '-.');

[ln,lt] = ndgrid(Etop.lon, Etop.lat);
[x,y] = projfwd(getm(gca),lt,ln);
contour(x, y, Etop.z, -[50 100 200], 'color', rgb('light brown'));

contourm(Raster.lt, Raster.ln, double(extractfun('bias')), [0 0], 'color', map1(3,:), 'linewidth', 1.0, 'linestyle', ':');
contourm(Raster.lt, Raster.ln, double(extractfun('mef')), [0 0], 'color', map2(3,:), 'linewidth', 1.0);

setm(h.ax, 'frame', 'off');

textm(Loc.lat(Indices.loc), Loc.lon(Indices.loc), Indices.letters, 'horiz', 'center');
% textm(Loc.lat(isrt), Loc.lon(isrt), lbls(lsrt), 'horiz', 'center');

if exportfigs
    % exportgraphics(h.fig, 'regions_map.png', 'resolution', 150);
    export_fig(h.fig, 'region_map_test', '-r150', '-png', '-nocrop');
end

%% ... now classify points based on those hand-drawn blobs

% Extract coordinates of blobs in normalized figure coordinates

im = imread('../figures_paper/regions_map_annotated_blobs.png');
L = imsegkmeans(im,7);

xfig = linspace(0, 1, size(L,2));
yfig = linspace(0, 1, size(L,1));
yfig = yfig(end:-1:1);

[xblob, yblob] = deal(cell(6));
for ii = 1:6
    bw = L == (ii+1);
    [b, ll] = bwboundaries(bw);
    [~,imax] = max(cellfun(@(x) size(x,1), b));

    xblob{ii} = xfig(b{imax}(:,2));
    yblob{ii} = yfig(b{imax}(:,1));

end

% Check math by plotting

% h.ax2 = axes('Position', [0 0 1 1], 'color', 'none');
% hold(h.ax2, 'on');
% cellfun(@(x,y) plot(x,y), xblob, yblob);
% 
% set(h.ax2, 'xlim', [0 1], 'ylim', [0 1], 'visible', 'off');

% Convert from figure coordinates to lat/lon

pos = plotboxpos(h.ax);
m = getm(h.ax);

Raster.blobmask = zeros(size(Raster.lt));

for ii = 1:6
    xnew = interp1([pos(1) pos(1)+pos(3)], h.ax.XLim, xblob{ii});
    ynew = interp1([pos(2) pos(2)+pos(4)], h.ax.YLim, yblob{ii});

    [xlat, xlon] = projinv(m, xnew, ynew);
    
    polymask = inpolygon(wrapTo360(Raster.ln), Raster.lt, wrapTo360(xlon), xlat);
    Raster.blobmask(polymask) = ii;

end

% Note:
% 1 = inner shelf
% 2 = inner front
% 3 = shelf break
% 4 = aleutians
% 5 = middle shelf
% 6 = northwest

%% Regional Taylor and target diagrams

col = Raster.blobmask;
clim = [-0.5 6.5];

cmap = ["efe645" "e935a1", "00e3ff" "e1562c" "537eff" "00cb85"]';
cmap = cell2mat(arrayfun(@hex2rgb, cmap, 'uni', 0));

cmap = [adjustcol(Col.gray, -0.2); cmap];
sz = 8;

% Note:
% 1 = inner shelf  #324068
% 2 = inner front  #8A4F21
% 3 = shelf break  #5A7B29
% 4 = aleutians    #874290
% 5 = middle shelf #226258
% 6 = northwest    #8A2A44

bloblabel = ["b) Inner Shelf", "c) Inner Front", "d) Shelf Break", "e) Alaska Peninsula", "f) Middle Shelf", "g) Northwest"];

% Axis setup

rtar = 2.5:-0.5:0;
thtar = linspace(0,2*pi,200);
xtar = rtar .* cos(thtar');
ytar = rtar .* sin(thtar');

h = plotgrid('size', [3 2], 'mar', 0.02, 'sh', 0.05, 'sv', 0.1, 'mt', 0.1);

% h.fig.Position(3:4) = [792 612];
h.ax = subgridaxes(h.ax, 1, 2);
h.ax = reshape(h.ax, 2, 6);

% Plotting data for Taylor and target

Plt.col = col(:);
Plt.r = real(acos(Rplt.cor(:))); % correlation
Plt.th = Rplt.stdnorm(:); % normalized standard dev.

fac = ones(size(Rplt.stdnorm));
fac(Rplt.stdnorm < 1) = -1;

Plt.x = Rplt.crmsd(:).*fac(:); % RMSD * sign(delta stdev)
Plt.y = Rplt.bias(:); % bias


% Plot

for ii = 1:6
    
    P = Plt;
    P.col(P.col ~= ii) = 0;
    [~,isrt] = sort(P.col);
    P = structfun(@(x) x(isrt), P, 'uni', 0);

    % Taylor plot
    
    h.tax(ii) = tayloraxis(h.ax(1,ii), 'stdmax', 2.5, 'npanel', 2);
    h.ps(ii) = polarscatter(P.r, P.th, sz, P.col, 'filled');
    
    % Target plot
    
    h.tar(ii) = patch(h.ax(2,ii), xtar, ytar, 'w');
    hold(h.ax(2,ii), 'on');
    
    fac = ones(size(data{2,3}));
    fac(data{2,3} < 1) = -1;
    
    h.sc(ii) = scatter(h.ax(2,ii), P.x, P.y, sz, P.col, 'filled');
    
end

% Prettify

set(h.fig, 'color', 'w');

set([h.tax.ax h.ax(2,:)], 'clim', clim, 'colormap', cmap, 'fontsize', 6);
set([h.tax.ax], 'RColor', 'k', 'thetalim', [0 pi*0.67]);
set([h.tax.rmsdline], 'color', rgb('gray'));
set(h.ax(2,:), 'color', 'none', 'layer', 'top', 'box', 'off', ...
    'xaxisloc', 'origin', 'yaxisloc', 'origin', ...
    'dataaspectratio', [1 1 1], 'xlim', [-2.5 2.5], 'ylim', [-2.5 2.5]);
set([h.ps h.sc], 'markerfacealpha', 0.1);
set(h.tar, 'linestyle', ':', 'edgecolor', rgb('gray'));

h.tax(1).ax.RAxis.Label.String = 'Normalized Standard Deviation';
h.tax(1).ax.ThetaAxis.Label.String = 'Correlation';
h.tax(1).ax.ThetaAxis.Label.Position = [0.8 3 0];
h.tax(1).ax.ThetaAxis.Label.Rotation = -45;

xlabel(h.ax(2,1), "RMSD*sgn(\Delta\sigma)");
h.ax(2,1).XAxis.Label.HorizontalAlignment = 'center';
ylabel(h.ax(2,1), "Bias");

% xlabel(h.tax(1).ax, 'Normalized Standard Deviation');
% thetalabel(h.tax(1).ax, 'Correlation');

h.lbl = labelaxes([h.tax.ax], bloblabel, 'northwestoutsideabove', 'vbuffer', 0.1);
set(h.lbl, {'color'}, num2cell(cmap(2:end,:),2));
h.lbl(1).Color = adjustcol(h.lbl(1).Color, -0.2);

% labelaxes([h.tax.ax], ["a)" "b)" "c)" "d)" "e)" "f)"], 'northwestoutside', 'hbuffer', 0.05)

if exportfigs
    export_fig('region_taylor_target', h.fig, '-png', '-r150', '-nocrop');
end

%% ... regional metric summaries

tmp = struct2array(structfun(@(x) x(:), Rplt, 'uni', 0));

Ravg = array2table(splitapply(@(x) mean(x,1,'omitmissing'), tmp, Raster.blobmask(:)+1), ...
    'variablenames', fieldnames(Rplt), 'rownames', compose("%d",0:6));

Rstd = array2table(splitapply(@(x) std(x,0,1,'omitmissing'), tmp, Raster.blobmask(:)+1), ...
    'variablenames', fieldnames(Rplt), 'rownames', compose("%d",0:6));




%% Compare moorings and models, surface and bottom (updated)

% Set order for plots and analysis (top to bottom)

nmoorplt = 4;

[~,idx] = ismember(('ikhs')', Indices.letters);
lbls = ["M8 (i)", "M5 (k)", "M4 (h)", "M2 (s)"];

% Match up with data arrays

midx = [Indices.moor{idx}]'; % index in Mprof array
lidx = Indices.loc(idx);     % index in Sim array (3rd dim)

% Set up axes, including copies for compressed early period and zoomed
% later period

tplt = [datetime(2005,1,1) datetime(2019,1,1), datetime(2023,9,1)];

ylim = [-2 15];
tedge = datetime(2005:2024,7,1);

h = plotgrid('size', [nmoorplt 1], 'sp', 0.08, 'ml', 0.05, 'mr', 0.01);
h.ax = subgridaxes(h.ax, 1, [0.4 0.6]);

for ii = 1:nmoorplt

    ttmp = cell(1,2);

    % Mooring data, converted to timetable

    ttmp{1} = timetable(Mprof(midx(ii)).btemp', Mprof(midx(ii)).stemp', ...
            'rowtimes', Mprof(midx(ii)).t, 'variablenames', {'btemp', 'stemp'});
    
    % Sim data, converted to timetable

    ttmp{2} = timetable(Sim.t, Sim.btemp(:,5,lidx(ii)), Sim.stemp(:,5,lidx(ii)), ...
        'variablenames', {'btemp', 'stemp'});
    ttmp{2} = ttmp{2}(ttmp{2}.Time >= tplt(1),:);

    % Synchronize for mooring skill stats

    tmskill = synchronize(ttmp{1}, ttmp{2}, 'daily', 'mean');
    moorskill(ii,1) = skillstats(tmskill.btemp_1, tmskill.btemp_2);
    moorskill(ii,2) = skillstats(tmskill.stemp_1, tmskill.stemp_2);

    % Retime mooring to daily for plot legibility and comparability

    ttmp{1} = retime(ttmp{1}, 'daily', 'mean');

    % Calculate mixing date for each year of each

    for itmp = 1:2
        ismixed = abs(ttmp{itmp}.stemp - ttmp{itmp}.btemp) < 0.2;
        [b,n,bi] = RunLength(ismixed);
        imix = bi(b==1 & n>20);
        ibin = findgroups(discretize(ttmp{itmp}.Time(imix), tedge));
        imix = splitapply(@(x) x(1), imix, ibin);

        ttmp{itmp}.mixmask = false(height(ttmp{itmp}),1);
        ttmp{itmp}.mixmask(imix) = true;
    end

    % Synchronize and save mixing event tables for analysis

    toctbin = datetime(2004:2024, 7, 1);
    toct = cell(1,2);
    for itmp = 1:2
        tmix = ttmp{itmp}.Time(ttmp{itmp}.mixmask);
        bn = discretize(tmix, toctbin);
        toct{itmp} = timetable(datetime(year(toctbin(bn))',1,1), days(tmix - toctbin(bn)'));
    end
    tt_mixing{ii} = synchronize(toct{:});
    tt_mixing{ii}.Properties.VariableNames = {'moor', 'model'};
    tt_mixing{ii}.delta = tt_mixing{ii}.model - tt_mixing{ii}.moor;

    % Plot

    for iax = 1:2

        % timeseries

        plot(h.ax(1,iax,ii), ttmp{1}, 'btemp', 'linewidth', 1.0, 'color', Col.moor);
        hold(h.ax(1,iax,ii), 'on');
        plot(h.ax(1,iax,ii), ttmp{2}, 'btemp', 'linewidth', 1.0, 'color', Col.b10k);
        plot(h.ax(1,iax,ii), ttmp{1}, 'stemp', 'linewidth', 1.0, 'color', adjustcol(Col.moor, 0.5));
        plot(h.ax(1,iax,ii), ttmp{2}, 'stemp', 'linewidth', 1.0, 'color', adjustcol(Col.b10k, 0.5));
        
        % mixing date option 1

        % plot(h.ax(1,iax,ii), repmat(ttmp{1}.Time(ttmp{1}.mixmask),1,2)', ylim, 'color', Col.moor, 'linestyle', '-.');
        % plot(h.ax(1,iax,ii), repmat(ttmp{2}.Time(ttmp{2}.mixmask),1,2)', ylim, 'color', Col.b10k, 'linestyle', '-.');

        % mixing date option 2

        plot(h.ax(1,iax,ii), ttmp{1}.Time(ttmp{1}.mixmask), ylim(2)-diff(ylim)*0.1, 'color', Col.moor, 'markerfacecolor', Col.moor, 'linestyle', 'none', 'marker', 'v', 'markersize', 4);
        plot(h.ax(1,iax,ii), ttmp{2}.Time(ttmp{2}.mixmask), ylim(2)-diff(ylim)*0.1, 'color', Col.b10k, 'markerfacecolor', 'none',   'linestyle', 'none', 'marker', 'v', 'markersize', 4);

    end

end

% Prettify

tks = datetime(year(tplt(1)):year(tplt(3)), 1, 1);
mtks = datetime(year(tplt(1)), 1:calmonths(between(tplt(1), tplt(3), 'months')), 1);

set(h.ax, 'ylim', ylim, 'box', 'off', 'xtick', tks, 'xminortick', 'on', 'tickdir', 'out');
cellfun(@(x) set(x, 'MinorTickValues', mtks), get(h.ax, 'XAxis'));
set(h.ax(1,1,1:end-1), 'xticklabel', '');
set(h.fig, 'color', 'w');

set(h.ax(1,1,:), 'xlim', tplt(1:2));
set(h.ax(1,2,:), 'xlim', tplt(2:3));
set(h.ax(1,2,:), 'ycolor', 'none');
set(h.ax(1,2,:), 'xtick', tks(tks>tplt(2)));

arrayfun(@(x) xlabel(x, ''), h.ax); % remove default x-axis label

h.lbl = labelaxes(h.ax(1,1,:), lbls, 'northwestoutsideabove');

h.fig.Position(3:4) = [800 400];

if exportfigs
    export_fig('mooring_temps', h.fig, '-r150', '-png');
end

%% ... double-check ice fraction cover at M8 vs SSMI

m8ltln = [Mprof(4).latitude(:) wrapTo360(Mprof(4).longitude(:))];
m8ltln = mean(m8ltln, 'omitnan');

dm8 = pdist2([Grd.lat_rho(:) Grd.lon_rho(:)], m8ltln, dfun);
[~,imin] = min(dm8);

[xi8, eta8] = ind2sub(size(Grd.h), imin);
rawb10k = dir(fullfile(moxdir, 'bering10k/output/hindcasts/phys_202206_30_popup/Out/*_avg_*.nc'));
rawb10k = fullfile({rawb10k.folder}, {rawb10k.name});
Scs = struct('xi_rho', [xi8 1 1], 'eta_rho', [eta8 1 1]);
I = cellfun(@(x) ncinfo(x, 'ocean_time'), rawb10k);
Aice = ncstruct(rawb10k, 'aice', 'temp', 'ocean_time', Scs);
Aice.t = ncdateread(rawb10k{1}, 'ocean_time', Aice.ocean_time);
Aice.aice = Aice.aice(:);


iceyr = 2005:2023;
for ii = 1:length(iceyr)
    Tmp = load("/Volumes/LaCie2023/SeaIce_SSMI/b10k_regrid/ssmr_" + iceyr(ii)+".mat");

    Ssmi.aice{ii} = Tmp.aice(xi8,eta8,:);
    Ssmi.t{ii} = Tmp.t;
end

Ssmi.t = cat(2, Ssmi.t{:});
Ssmi.aice = reshape(cat(3, Ssmi.aice{:}), 1, []);

%% ... plot 

hfig = figure;
yyaxis left;

ii = 1;
ttmp = cell(1,2);

% Mooring data, converted to timetable

ttmp{1} = timetable(Mprof(midx(ii)).btemp', Mprof(midx(ii)).stemp', ...
        'rowtimes', Mprof(midx(ii)).t, 'variablenames', {'btemp', 'stemp'});

% Sim data, converted to timetable

ttmp{2} = timetable(Sim.t, Sim.btemp(:,5,lidx(ii)), Sim.stemp(:,5,lidx(ii)), ...
    'variablenames', {'btemp', 'stemp'});
ttmp{2} = ttmp{2}(ttmp{2}.Time >= tplt(1),:);

% Plot ice

yyaxis left;
h.lnice(1) = area(Ssmi.t, Ssmi.aice);
hold on
h.lnice(2) = area(Aice.t, Aice.aice);

% Plot surface temps

yyaxis right;
h.lnsst(1) = plot(ttmp{1}, 'stemp');
hold on
h.lnsst(2) = plot(ttmp{2}, 'stemp');

tk = datetime(2019,1:1:48,1);
set(gca, 'xlim', minmax(Aice.t), 'xtick', tk);

% Prettify

set(h.lnice(2), 'facecolor', Col.b10k, 'facealpha', 0.2, 'edgecolor', 'none')
set(h.lnice(1), 'facecolor', Col.moor, 'facealpha', 0.2, 'edgecolor', 'none')
set(h.lnsst(2), 'color', Col.b10k, 'linestyle', '-');
set(h.lnsst(1), 'color', Col.moor, 'linestyle', '-');

yyaxis left;
set(gca, 'ycolor', 'k');
ylabel('Fraction ice cover');
yyaxis right;
set(gca, 'ycolor', 'k');
ylabel('SST (\circC)');


set(hfig, 'color', 'w');
hfig.Position(3:4) = [800 400];

legend([h.lnice' h.lnsst], {'Ice cover (SSMI)', 'Ice cover (Bering10K)', 'M8 SST', 'Bering10K SST'});

if exportfigs
    export_fig('sst_vs_aice', hfig, '-r150', '-png');
end

%% Compare moorings and model, surface and bottom

% % [~,idx] = ismember(('pokmglnv')', Indices.letters);
% [~,idx] = ismember(('ikhs')', Indices.letters);
% lbls = ["M8 (i)", "M5 (k)", "M4 (h)", "M2 (s)"];
% 
% tplt = [datetime(2005,1,1) datetime(2019,1,1), datetime(2023,9,1)];
% 
% h = plotgrid('size', [4 1], 'sp', 0.05, 'ml', 0.05, 'mr', 0.01);
% % h.ax = {h.ax};
% h.ax = subgridaxes(h.ax, 1, [0.4 0.6]);
% 
% nrun = 1;
% ylim = [-2 15];
% tedge = datetime(2005:2024,7,1);
% 
% for ii = 1:length(idx)
% 
%     midx = Indices.moor{idx(ii)};
%     pidx = Indices.puf{idx(ii)};
% 
%     for im = 1:length(midx)
%         ttmp = timetable(Mprof(midx(im)).btemp', Mprof(midx(im)).stemp', ...
%             'rowtimes', Mprof(midx(im)).t, 'variablenames', {'btemp', 'stemp'});
%         ttmp = retime(ttmp, 'daily', 'mean');
% 
%         plot(h.ax(1,1,ii), ttmp.Time, ttmp.btemp, 'k', ...
%                        ttmp.Time, ttmp.stemp, 'b');
%         hold(h.ax(1,1,ii),'on');
% 
%         plot(h.ax(1,2,ii), ttmp.Time, ttmp.btemp, 'k', ...
%                        ttmp.Time, ttmp.stemp, 'b');
%         hold(h.ax(1,2,ii),'on');
% 
%         ismixed = abs(ttmp.stemp - ttmp.btemp) < 0.2;
%         [b,n,bi] = RunLength(ismixed);
%         imix = bi(b==1 & n>20);
%         ibin = findgroups(discretize(ttmp.Time(imix), tedge));
%         imix = splitapply(@(x) x(1), imix, ibin);
% 
%         plot(h.ax(1,1,ii), [ttmp.Time(imix) ttmp.Time(imix)]', ylim, 'color', Col.moor, 'linestyle', '-.');
%         plot(h.ax(1,2,ii), [ttmp.Time(imix) ttmp.Time(imix)]', ylim, 'color', Col.moor, 'linestyle', '-.');
% 
%     end
% 
%     lidx = Indices.loc(idx(ii));
% 
%     plot(h.ax(1,1,ii), Sim.t, Sim.btemp(:,5,lidx), 'r', ...
%                    Sim.t, Sim.stemp(:,5,lidx), 'g');
% 
%     plot(h.ax(1,2,ii), Sim.t, Sim.btemp(:,5,lidx), 'r', ...
%                    Sim.t, Sim.stemp(:,5,lidx), 'g');
% 
% 
%     ismixed = abs(Sim.stemp(:,5,lidx) - Sim.btemp(:,5,lidx)) < 0.2;
%     [b,n,bi] = RunLength(ismixed);
%     imix = bi(b==1 & n>20);
%     ibin = findgroups(discretize(Sim.t(imix), tedge));
%     imix = splitapply(@(x) x(1), imix, ibin);
% 
%     plot(h.ax(1,1,ii), [Sim.t(imix) Sim.t(imix)]', ylim, 'color', Col.b10k, 'linestyle', '-.');
%     plot(h.ax(1,2,ii), [Sim.t(imix) Sim.t(imix)]', ylim, 'color', Col.b10k, 'linestyle', '-.');
% 
% 
%     % plot(h.ax{1}(ii), Sim.t(imix), Sim.btemp(imix,5,lidx), 's', ...
%     %     'markerfacecolor', Col.b10k, 'markeredgecolor', 'k');
%     % 
%     % plot(h.ax{2}(ii), Sim.t(imix), Sim.btemp(imix,5,lidx), 's', ...
%     %     'markerfacecolor', Col.b10k, 'markeredgecolor', 'k');
% 
% end
% 
% % cmap = cptcmap('d3cat20');
% 
% set(findall(h.ax, 'color', 'k'), 'color', adjustcol(Col.moor,  0.5), 'linewidth', 1.5); %cmap(1,:));
% set(findall(h.ax, 'color', 'b'), 'color', adjustcol(Col.moor, -0.0), 'linewidth', 1.5); %cmap(2,:));
% set(findall(h.ax, 'color', 'r'), 'color', adjustcol(Col.b10k,  0.5), 'linewidth', 1.5); %cmap(3,:));
% set(findall(h.ax, 'color', 'g'), 'color', adjustcol(Col.b10k, -0.0), 'linewidth', 1.5); %cmap(4,:));
% 
% tks = datetime(year(tplt(1)):year(tplt(3)), 1, 1);
% mtks = datetime(year(tplt(1)), 1:calmonths(between(tplt(1), tplt(3), 'months')), 1);
% 
% set(h.ax, 'ylim', ylim, 'box', 'off', 'xtick', tks, 'xminortick', 'on', 'tickdir', 'out');
% cellfun(@(x) set(x, 'MinorTickValues', mtks), get(h.ax, 'XAxis'));
% set(h.ax(1,1,1:end-1), 'xticklabel', '');
% set(h.fig, 'color', 'w');
% 
% set(h.ax(1,1,:), 'xlim', tplt(1:2));
% set(h.ax(1,2,:), 'xlim', tplt(2:3));
% set(h.ax(1,2,:), 'ycolor', 'none');
% set(h.ax(1,2,:), 'xtick', tks(tks>tplt(2)));
% 
% h.lbl = labelaxes(h.ax(1,1,:), lbls, 'northwestoutsideabove');
% 
% h.fig.Position(3:4) = [800 400];
% 
% if exportfigs
%     export_fig('mooring_temps', h.fig, '-r150', '-png');
% end

%% Inner front analysis
%

%% ... Read model data along transect

reload = false;

[~,imin] = arrayfun(@(x,y) min(dfun([x y], [Grd.lat_rho(:) Grd.lon_rho(:)])), iflat, iflon);
[ifxi, ifeta] = ind2sub(size(Grd.h), imin);
etaslice = mode(ifeta);

if reload
    Scs = struct('xi_rho', [80 40 1], 'eta_rho', [etaslice 1 1]);
    
    ftemp = fullfile(ppsimdir, 'Level1', 'B10K-K20nobio_CORECFS_daily_1995-1999_average_temp.nc');
    fsalt = fullfile(ppsimdir, 'Level1', 'B10K-K20nobio_CORECFS_daily_1995-1999_average_salt.nc');
    
    Slice = mergestruct(ncstruct(ftemp, Scs), ncstruct(fsalt, 'salt', Scs));
    Slice.t = ncdateread(F, 'ocean_time');
else
    Slice = load(fullfile(moxdir, 'kearney/Manuscripts/2022_pcodpopup/innerfront_model_slice.mat'), 'Slice');
    Slice = Slice.Slice;
end

% Some setup geometry calcs

[~,Slice.zw] = calcromsz(Grd.h(80:119, etaslice), Slice.zeta, 30);
dx = dfun([Slice.lat_rho(1:end-1) Slice.lon_rho(1:end-1)], [Slice.lat_rho(2:end) Slice.lon_rho(2:end)]);
Slice.x = cumsum([0; dx]);


%% ... Calculate stratification along model transect

dz = diff(Slice.zw, 1, 3);

zwrel = -(Slice.zw(:,:,end,:) - Slice.zw);
zmrel = (zwrel(:,:,1:end-1,:) + zwrel(:,:,2:end,:))./2;

sabs = max(Slice.salt,0); % GSW wants absolute salinity... close enough
tconc = gsw_CT_from_pt(sabs,Slice.temp); % convert from ROMS potential temp to conservative temp

% dens = sw_dens0(max(Slice.salt,0), Slice.temp);
dens = gsw_sigma0(sabs, tconc)+1000;

rhomean = sum(dens.*dz,3)./sum(dz,3);
g = -9.806; % m s^-2
penergy = (dens - rhomean).*g.*zmrel;  
  
Slice.strat = sum(penergy.*dz,3)./sum(dz,3);

%% ... Stratification for buoys

 sabs = 32;

for ii = 1:length(Iprof)
    
    % New

    z1 = -Iprof(ii).depth;
    x1 = Iprof(ii).t';
    v1 = Iprof(ii).temperature';

    if nnz(~isnan(v1)) == 0
        Iprof(ii).strat = nan(size(x1));
        fprintf('Skipping %d, no data\n', ii);
        continue;
    end

    if z1(1) < 0
        z1 = [0; z1];
        v1 = [v1(1,:); v1];
    end
    n = sum(~isnan(v1),1);

    z2 = (ceil(min(z1)):max(max(z1),0))';
    
    nz = length(z2);
    nx = length(x1);

    try
        v2 = interp1(z1, v1, z2, 'pchip', NaN);
    catch
        v2 = nan(nz, nx);
        try
            v2(:,n>1) = interp1(z1, v1(:,n>1), z2, 'pchip', NaN);
        catch ME
            rethrow(ME)
        end
    end

    p2 = gsw_p_from_z(z2, nanmean(Iprof(ii).latitude(:)));

    tc2 = gsw_CT_from_t(ones(nz,nx).*sabs, v2, p2);
    dens2 = gsw_sigma0(ones(nz,nx).*sabs, tc2)+1000;
    
    Iprof(ii).strat = nan(1,nx);
    for ix = 1:nx
        Iprof(ii).strat(ix) = trapz(z2, (dens2(:,ix) - nanmean(dens2(:,ix))).*g.*z2)./diff(minmax(z2));
    end



    % % Old
    % 
    % x1 = Iprof(ii).t;
    % z1 = Iprof(ii).depth;
    % v1 = Iprof(ii).temperature;
    % n = sum(~isnan(v1),2);
    % 
    % z2 = min(z1):max(z1);
    % try
    %     v2 = interp1(z1, v1', z2, 'pchip', NaN);
    % catch
    %     v2 = nan(length(z2), length(x1));
    %     try
    %         v2(:,n>1) = interp1(z1, v1(n>1,:)', z2, 'pchip', NaN);
    %     end
    % end
    % 
    % sabs = 32;
    % papprox = gsw_p_from_z(-z2, nanmean(Iprof(ii).latitude(:)));
    % tconc = gsw_CT_from_t(ones(size(v2)).*sabs, v2, papprox'.*ones(size(v2)));
    % dens = gsw_sigma0(ones(size(tconc)).*sabs, tconc)+1000;
    % 
    % % sabs = ones(size(v2))*32; % estimate salinity using mean-ish for this slice in the model
    % % pish = repmat(gsw_p_from_z(-z2, nanmean(Iprof(ii).latitude(:)))', size(v2,2));
    % % tconc = gsw_CT_from_t(sabs,v2,gsw_p_from_z(-z2,lattmp));
    % % dens2 = gsw_sigma0(sabs, tconc)+1000;
    % 
    % % dens = sw_dens0(ones(size(v2))*32, v2);
    % rhomean = mean(dens, 1);
    % penergy = (dens - rhomean).*g.*-z2';
    % Iprof(ii).strat = mean(penergy,1);

end

%% ... Supplement: mooring to model to real SI comparison

targetlat = nanmean(Iprof(1).latitude(:));
targetlon = nanmean(Iprof(1).longitude(:));
targett = mean(Iprof(1).t);

d = dfun([targetlat, targetlon], [Pel.LATITUDE, Pel.LONGITUDE]);
mask = d<1;
idx = find(mask);
[~,imin] = min(abs(Pel.t(mask) - targett));

idx = idx(17);

sabs = 32;

% The "real" CTD data

v0 = Pel.Temperature(idx,:)';
p0 = pelpres;
v0 = v0(end:-1:1);
p0 = p0(end:-1:1);
z0 = gsw_z_from_p(p0, targetlat);
mask = ~isnan(v0);
v0 = v0(mask);
z0 = z0(mask);
p0 = p0(mask);

tc0 = gsw_CT_from_t(ones(size(v0)).*sabs, v0, p0);
dens0 = gsw_sigma0(ones(size(tc0)).*sabs, tc0)+1000;
dens0 = gsw_sigma0(ones(size(tc0)).*sabs, tc0)+1000;
strat0 = mean((dens0 - mean(dens0)).*g.*-z0);

strat0 = trapz(z0, (dens0 - mean(dens0)).*g.*z0)./diff(minmax(z0));

% What the mooring at this location would sample

z1 = -Iprof(1).depth;
z1(end) = min(z0);
v1 = interp1(z0, v0, z1, 'linear', 'extrap');
if z1(1) < 0
    z1 = [0; z1];
    v1 = [v1(1); v1];
end

z2 = ceil(min(z1)):max(max(z1),0);
v2 = interp1(z1, v1', z2, 'pchip', NaN);

p2 = gsw_p_from_z(z2, targetlat);
tc2 = gsw_CT_from_t(ones(size(v2)).*sabs, v2, p2);
dens2 = gsw_sigma0(ones(size(tc2)).*sabs, tc2)+1000;

strat2 = trapz(z2, (dens2 - mean(dens2)).*g.*z2)./diff(minmax(z2));

% What a perfect model would look like

[zm3,zw3] = calcromsz(-min(z0), 0, 30);
zm3 = zm3(:); zw3 = zw3(:);
dz = diff(zw3);
z0hi = min(z0):0.1:0;
v0hi = interp1(z0,v0,z0hi);

ilayer = discretize(z0hi, zw3);
v3 = splitapply(@mean, v0hi, ilayer);

tc3 = gsw_CT_from_pt(ones(size(v3)).*sabs,v3); % convert from ROMS potential temp to conservative temp
dens3 = gsw_sigma0(ones(size(v3)).*sabs, tc3)+1000;

rhomean = sum(dens3'.*dz)./sum(dz);
strat3 = sum((dens3' - rhomean).*g.*zm3.*dz)./sum(dz);  

h = plot(v0, z0, 'k', ...
         v1, z1, 'ro', ...
         v2, z2, 'r', ...
         v3, zm3, 'bx'); 
gridxy([], zw3, 'color', adjustcol(Col.b10k, 0.5));

% set(h, {'color'}, num2cell(jet(length(h)),2));

%% ... Converting to along-transect space

% Projected coords for locs

latlim = minmax(Slice.lat_rho, 'expand', 0.2);
lonlim = minmax(Slice.lon_rho, 'expand', 1);

htmp = figure;
myproj = 'mercator';
worldmap(latlim, lonlim);
setm(gca, 'mapprojection', myproj);
mproj = getm(gca);
close(htmp);

% Model slice

[xslice, yslice] = projfwd(mproj, Slice.lat_rho, Slice.lon_rho);
[xsS, xsN] = xy2sn(xslice, yslice, xslice, yslice);

slim = minmax(xsS);

% Raster coords

[rx, ry] = projfwd(mproj, Raster.lt, Raster.ln);
[rS, rN, rL] = xy2sn(rx(:), ry(:), xslice, yslice);
rS = reshape(rS, size(rx));
rN = reshape(rN, size(rx));

% Groudfish stations

[tmx, tmy] = projfwd(mproj, TrawlMeta.lat, TrawlMeta.lon);
[tmS, tmN] = xy2sn(tmx, tmy, xslice, yslice);

[srx, sry] = projfwd(mproj, Srep.latitude, Srep.longitude);
[srS, srN] = xy2sn(srx, sry, xslice, yslice);

% Inner front survey

[ifx, ify] = projfwd(mproj, iflat, iflon);
[ifS, ifN] = xy2sn(ifx, ify, xslice, yslice);

nlim = minmax(ifN([1:16 18]), 'expand', 1);

% Borders

[blat, blon] = borders('alaska');
[bx, by] = projfwd(mproj, blat, blon);
[bS, bN] = xy2sn(bx', by', xslice, yslice);

% Conversions

km2s = @(s) interp1(Slice.x, xsS, s);
s2km = @(s) interp1(xsS, Slice.x, s);


%% ... Plot stratification analysis

extractfun = @(fld) arrayfun(@(x) x.(fld)(2), Raster.Skill);

% Note: 20 J/m contour is pretty good marker of front

% Figure setup

h = plotgrid('size', [1 1], 'mar', 0.05, 'ml', 0.1, 'mr', 0.1);
h.fig.Position(3:4) = [612 792]; % 8.5in x 11in 
set(h.fig, 'color', 'w');
h.ax = subgridaxes(h.ax, [0.1 0.7 0.2], 1);

%--------------------------
% Bathymetry along transect
%--------------------------

h.bathy = plot(h.ax(3), Slice.x, Grd.h(80:119, etaslice), 'k');
set(h.ax(3), 'ydir', 'reverse', 'ylim', [0 200], 'ytick', 0:50:200, 'ygrid', 'on', 'xlim', s2km(slim));
hold(h.ax(3), 'on');

tprop = {'rotation', 90, 'fontsize', 8, 'horiz', 'left', 'vert', 'base'};
aprop = {'color', adjustcol(Col.gray, -0.2), 'linewidth', 1.5};

[xi,yi] = polyxpoly(Slice.x, Grd.h(80:119,etaslice), minmax(Slice.x), [50 50]);
h.barrow(1) = plot(h.ax(3), [xi xi], [yi -50], aprop{:});
text(xi-4, yi-10, '50 m', tprop{:});

[xi,yi] = polyxpoly(Slice.x, Grd.h(80:119,etaslice), minmax(Slice.x), [100 100]);
h.barrow(2) = plot(h.ax(3), [xi xi], [yi -50], aprop{:});
text(xi-4, yi-10, '100 m', tprop{:});

[xi,yi] = polyxpoly(Slice.x, Grd.h(80:119,etaslice), minmax(Slice.x), [200 200]);
h.barrow(3) = plot(h.ax(3), [xi xi], [yi -50], aprop{:});
text(xi-4, yi-10, '200 m', tprop{:});

% Etopo

[xetop,yetop] = ndgrid(Etop.lon, Etop.lat);
fetop = griddedInterpolant(xetop, yetop, -Etop.z); %, Grd.lon_rho(:), Grd.lat_rho(:));

zetopslice = fetop(Slice.lon_rho, Slice.lat_rho);

plot(h.ax(3), Slice.x, zetopslice, ':k');

%--------------------------
% Along-transect map
%--------------------------

% MEF contours

h.mef = contour(h.ax(1), rS, rN, extractfun('mef'), -1:0.1:1);
shading(h.ax(1), 'flat');
hold(h.ax(1), 'on');
% plot(h.ax(3), bS, bN, 'k'); % TODO: clip bathy

% Model grid transect

plot(h.ax(1), xsS, xsN, 'k.'); 

% Inner front buoys

plot(h.ax(1), ifS, ifN, 'marker', 'v', 'linestyle', 'none', 'markeredgecolor', adjustcol(Col.moor,-0.3), 'markerfacecolor', Col.moor);

% Groundfish sampling stations

isnear = tmN >= nlim(1) & tmN <= nlim(2);
isnear2 = ismember(Srep.stationid, TrawlMeta.id(isnear));
isverynear = ismember(TrawlMeta.id, {'L-01', 'L-18', 'M-18', 'M-01', 'K-18', 'N-01', 'N-02'});
isverynear2 = ismember(Srep.stationid, {'L-01', 'L-18', 'M-18', 'M-01', 'K-18', 'N-01', 'N-02'});

plot(h.ax(1), tmS(isnear), tmN(isnear), 'marker', 'o', 'linestyle', 'none', ...
    'markeredgecolor', adjustcol(Col.gfish, -0.3), 'markerfacecolor', Col.gfish);
% text(h.ax(1), tmS(isverynear), tmN(isverynear), "  "+TrawlMeta.id(isverynear), 'fontsize', 8, 'color', adjustcol(Col.gfish, -0.3));

% Prettify

tk = [0:50:400 445.28];
set(h.ax(1), 'dataaspectratio', [1 1 1], 'xlim', slim, 'ylim', nlim, ...
    'colormap', cmocean('curl'), 'clim', [-1 1], 'xtick', [], 'ytick', []);

% km2s(tk), 'xticklabel', compose('%g',tk), 'tickdir', 'out', 

%--------------------------
% Stratification index
%--------------------------

whis = 2;
tmp = boxplotstats(Slice.strat(:), whis);
modrel = tmp.hi; %150;
tmp = boxplotstats(cat(2, Iprof.strat), whis)
obsrel = tmp.hi;

% Model stratification vs time and transect

% h.strat = contourf(h.ax(2), Slice.x, doy(Slice.t, 'decyear'), permute(Slice.strat, [4 1 2 3]), 0:10:150, 'edgecolor', 'none');
h.strat = contourf(h.ax(2), Slice.x, doy(Slice.t, 'decyear'), permute(Slice.strat, [4 1 2 3])./modrel, 0:0.1:1, 'edgecolor', 'none');
shading(h.ax(2), 'flat');
hold(h.ax(2), 'on');

% Buoy locs along this slice transect

h.ax2b = axes('position', h.ax(2).Position);

tbin = [Slice.t - days(3.5); Slice.t(end)+days(3.5)];

for ii = [1:7 10:16 18] %1:length(Iprof)
    tmp = timetable(Iprof(ii).t, Iprof(ii).strat');
   
    % tmp = retime(tmp, 'daily', 'mean');
    tmp = retime(tmp, tbin, 'mean'); 
    tmp = tmp(~isnan(tmp.Var1),:);

    % scatter(h.ax2b, ones(length(Iprof(ii).strat),1)*s2km(ifS(ii)), ...
    %     doy(Iprof(ii).t, 'decyear'), 5, max(0,Iprof(ii).strat'));
    h.scat(ii) = scatter(h.ax2b, ones(height(tmp),1)*s2km(ifS(ii)), ...
        doy(tmp.Time, 'decyear'), 15, max(0,tmp.Var1)./obsrel);
    hold(h.ax2b, 'on');
end

% [xy, dist, t] = distance2curve([Slice.lon_rho Slice.lat_rho], [iflon iflat]);
% xbuoy = t .* x(end);
% tbuoy = arrayfun(@(x) minmax(x.t), Iprof, 'uni', 0);
% tbuoy = cat(1, tbuoy{:});
% 
% plot(h.ax(2), [xbuoy xbuoy]', doy(tbuoy', 'decyear'), '-bx');

% L-01 sampling time

% l1lonlat = Srep{isl1, {'longitude', 'latitude'}};
% l1lonlat(:,1) = wrapTo360(l1lonlat(:,1));
% [xyl1, dtmp, tl1] = distance2curve([Slice.lon_rho Slice.lat_rho], l1lonlat);
% xl1 = tl1 .* x(end);

plot(h.ax2b, s2km(srS(isnear2)), doy(Srep.start_time(isnear2), 'decyear'),  'marker', 'o', 'linestyle', 'none', ...
    'markeredgecolor', adjustcol(Col.gfish, -0.3), 'markerfacecolor', Col.gfish);

% Prettify

set(h.ax(2), 'colormap', cmocean('-gray', 15), 'clim', [0 1], 'ytick', 1997:2000, 'ylim', [1997 2000], 'xlim', s2km(slim));
h.cb = colorbar(h.ax(2), 'east');
h.cb.Position = [sum(h.ax(2).Position([1 3]))+0.01 h.ax(2).Position(2) 0.015 0.3];

% h.cb.Position = [sum(h.ax(2).Position([1 3]))-0.32 h.ax(2).Position(2) 0.3 0.015]; 
% xlabel(h.cb, 'Model-based Stratification (J m^{-2})');
set(h.cb, 'AxisLoc', 'in', 'TickDir', 'in', 'box', 'off', 'ticks', 0:0.2:1, 'ticklabels', '');

set(h.ax2b, 'colormap', cmocean('amp',15), 'xlim', h.ax(2).XLim, 'ylim', h.ax(2).YLim, 'visible', 'off', 'clim', [0 1]);
h.cb2 = colorbar(h.ax2b, 'east');
h.cb2.Position = [sum(h.cb.Position([1 3])) h.ax(2).Position(2) 0.015 0.3];

% h.cb2 = colorbar(h.ax2b, 'south');
% h.cb2.Position = [sum(h.ax(2).Position([1 3]))-0.32 h.ax(2).Position(2)+h.cb.Position(4) 0.3 0.015]; 
ylabel(h.cb2, 'Relative Stratification', 'rotation', 270);
set(h.cb2, 'AxisLoc', 'out', 'TickDir', 'in', 'box', 'off', 'ticks', 0:0.2:1);

h.ax(1).Position = h.ax(1).Position + [0 0.02 0 0];
h.cb3 = colorbar(h.ax(1), 'east');
pos = plotboxpos(h.ax(1));
h.cb3.Position = [sum(pos([1 3]))+0.01 pos(2) 0.015 pos(4)];
ylabel(h.cb3, 'Model Efficiency (MEF)', 'rotation', 270);

set(h.ax(3), 'box', 'off', 'xcolor', 'none');
labelaxes(h.ax(3), "Model Transect Bottom Depth (m)", 'southeast');

set(h.ax(2:3), 'tickdir', 'out');
xlabel(h.ax(2), 'Distance along transect (km)');
ylabel(h.ax(2), 'Year');
set(h.ax(2), 'xaxisloc', 'top');

h.ax(3).Position = h.ax(3).Position - [0 0.03 0 0];
set(h.ax(3), 'visible', 'on');

h.lbl = labelaxes(h.ax, {'(a)','(b)','(c)'}, 'northwestoutside');
set(h.lbl, 'fontsize', 12, 'fontweight', 'b', 'vert', 'top');
for ii = 1:length(h.lbl)
    h.lbl(ii).Position(1) = -0.085;
end

pause(1);
for ii = 1:3
    h.barrow(ii).Clipping = 'off';
    line2arrow(h.barrow(ii));
end

% Export

if exportfigs
    export_fig('inner_front_stratification', h.fig, '-r300', '-png');
end

% Coordinates of transect box (to add to map later)

sbox = slim([0 0 1 1 0]+1) + [1 1 -1 -1 1].*10;
nbox = nlim([0 1 1 0 0]+1) + [1 -1 -1 1 1]*10;
dx = sqrt(diff(sbox).^2 + diff(nbox).^2);
tfrac = union(linspace(0,1,100), cumsum(dx)./sum(dx));
snbox = interparc(tfrac, sbox, nbox, 'linear');
[xbox, ybox] = sn2xy(snbox(:,1), snbox(:,2), xslice, yslice);
[ltbox, lnbox] = projinv(mproj, xbox, ybox);

%% Map of everything important (no data, just locations)

h = plotgrid('size', [1 1], 'margin', 0);
h.fig.Position(3:4) = [800 600];
set(h.fig, 'color', 'w');

% ETOPO shaded relief

h.ax(2) = copyobj(h.ax(1), h.fig);
axes(h.ax(1));
worldmap(latlimbig, lonlimbig);

[x,y] = ndgrid(Etop.lon, Etop.lat);
[x,y] = projfwd(getm(h.ax(1)), y,x);
h.topo = pcolor(h.ax(1), x,y,Etop.z);
shading flat;
% h.topo = pcolorm(Etop.lat(1:2:end), Etop.lon(1:2:end), Etop.z(1:2:end,1:2:end)'); % no idea why this is timing out
colormap('white');
shadem([90 45], 2);
setm(h.ax(1), 'flinewidth', 1, 'fedgecolor', rgb('gray'));


% B10K grid

axes(h.ax(2));
h.roms = plotromsrho(Grd, Grd.h);
setm(h.ax(2), 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off');
set(h.ax(2), 'colormap', cmocean('deep'), 'clim', minmax(Grd.h));
set(h.roms, 'facealpha', 0.5);

h.cb = colorbar('west');
h.cb.Position([1 2 4]) = [0.01 0.5 0.3];
tk = [0 200 1000:2000:7000];
set(h.cb, 'ticks', tk, 'ticklabels', compose("%d m", tk), 'tickdir', 'out');

% Borders

plotm(Border.lat, Border.lon, 'k');

% Focus region

latliml = minmax(Loc.lat, 'expand');
lonliml = minmax(Loc.lon, 'expand');

h.box = plotm(latliml([1 2 2 1 1]), lonliml([1 1 2 2 1]), ':k');

% Survey

Strata = shapeprjread('/Volumes/LaCie2023/AlaskaShapefiles/gis_updated/EBS_NBS_2019.shp');
isin = ~ismember([Strata.STRATUM], [70 71 81]);
Strata = Strata(isin);
pstrata = polyshape([Strata.Lon], [Strata.Lat]);
[ln,lt] = boundary(pstrata);
plotm(lt, ln, 'linewidth', 1.5);

% Bering Sea marine area

% Marea = shapeprjread('/Volumes/LaCie2023/NaturalEarth/ne_10m_geography_marine_polys/ne_10m_geography_marine_polys.shp');
% isber = ismember({Marea.name}, {'Bering Sea','Gulf of Anadyr''', 'Bristol Bay', 'Karaginskiy Gulf', 'Norton Sound'});
% Ber = Marea(isber);
% 
% Pber = polyshape(wrapTo360([Ber.Lon]), [Ber.Lat]);
% Pber = polybuffer(Pber, 0.01); % get rid of small gap at dateline
% Pber = rmholes(Pber);
% 
% h.ber = patchm(Pber.Vertices(:,2), Pber.Vertices(:,1), 'r', ...
%     'edgecolor', 'k', 'linewidth', 2, 'linestyle', '-.', 'facecolor', 'none');

% Moorings

h.m = plotm(moorlat, moorlon, '^');
set(h.m, 'color', adjustcol(Col.moor, -0.5), 'markerfacecolor', Col.moor, 'markersize', 8);

% h.mextra = plotm(iflat, iflon, 'v');
% set(h.mextra, 'color', adjustcol(Col.moor, -0.5), 'markerfacecolor', Col.moor, 'markersize', 4);

% PUFs

h.puf = plotm(Id.latitude, Id.longitude, 'o');
set(h.puf, 'markerfacecolor', Col.puf, 'color', adjustcol(Col.puf, -0.5), 'markersize', 6);

% Transect box

h.tsect = plotm(ltbox, lnbox, 'k--');

set([h.tsect h.box], 'linewidth', 1.5);

% Fix meridian labels

hmer = findobj(h.ax(1), 'tag', 'MLabel');
set(hmer, {'string'}, cellstr(angl2str(wrapTo180(-150:-10:-200), 'ew', 'degrees', 0)))

% Export

if exportfigs
    export_fig('intromap', h.fig, '-png', '-r300');
end

%% Particle Tracking
%
% Note: ROMS output prep for this was done with the
% modify_for_particle_tracking.m script in the anaysis/ folder for the
% model raw output.  It involves reordering time to run backward and
% reversing the signs of all velocity components, plus some other minor
% adjustments to coordinate variables to make them play nicely with the
% ROMSPath expectations.

%% ... create particle input file for ROMSPath

% Using trawl stations as reference to areas of interest.

akpen = ["Z-05", "A-04", "A-05", "A-06", "B-06", "B-07", "B-08", "C-08", "C-09", "D-09", "D-10", ...
    "E-10", "E-11", "E-12", "F-12", "F-13"];
nw = ["N-31", "N-30", "N-29", "O-31", "O-30", "O-29", "O-28", "P-31", "P-30", "P-29", "P-28", ...
    "Q-31", "Q-30", "Q-29", "Q-28", "R-32", "R-31", "R-30", "R-29", "S-31", "S-30", "S-29", ...
    "T-30", "T-29"];

[tf,loc] = ismember([akpen nw], TrawlMeta.id);

nperploc = 20;
nploc = length(loc);

rng(8);
Part = struct;
Part.lat = TrawlMeta.lat(loc) + rand(nploc,nperploc)*km2deg(20*1.852);
Part.lon = TrawlMeta.lon(loc) + rand(nploc,nperploc)*km2deg(20*1.852);
Part = struct2table(structfun(@(x) x(:), Part, 'uni', 0));
Part.lon = wrapTo360(Part.lon);

% TODO: get starting depths, reformat to input file syntax

ftmp = scatteredInterpolant(Grd.lon_rho(:), Grd.lat_rho(:), Grd.h(:), 'nearest');
Part.depth = -ftmp(Part.lon, Part.lat) + 5;
Part.time = zeros(height(Part),1); % time after tracking initialization

Part = Part(:, {'lon', 'lat', 'depth', 'time'});
writetable(Part, 'particle_tracking/init_nw_akpen.csv', 'WriteVariableNames', false);

npart = height(Part);
npartak = length(akpen);

isak = false(nploc, nperploc);
isak(1:npartak,:) = true;
isak = isak(:);

%% ... plot particle tracks

pfile = 'particle_tracking/particletracks_nw_akpen.nc';
P = ncstruct(pfile);

padforplt = @(x) cat(2, x, nan(npart,1))';

npt = length(P.model_time);

% Projection calcs (we need to use a non-map axis because map axes do weird
% things with transparency)

htmp = plotgrid('size', [1 1]);
worldmap(latlimbig, lonlimbig);

mproj = getm(htmp.ax);
close(htmp.fig);

[bx,by] = projfwd(mproj, Border.poly.Vertices(:,2), Border.poly.Vertices(:,1));
[px,py] = projfwd(mproj, padforplt(P.lat), padforplt(P.lon));
[gx,gy] = projfwd(mproj, Grd.lat_rho, Grd.lon_rho);


pz = padforplt(P.zp);
pc = padforplt((1:npt).*ones(npart,1));

% Plot

h = plotgrid('size', [1 1], 'mar', 0.02, 'mb', 0.1);
set(h.fig, 'color', 'w');
plot(polyshape(bx,by), 'facecolor', rgb('gray'), 'edgecolor', 'none');
hold on;

% [h.cc, h.ch] = contour3(gx, gy, -Grd.h, -[50 100 200], 'color', rgb('gray'));
contour3(gx, gy, -Grd.h, -[ 50  50], 'color', ones(1,3)*0.8);
contour3(gx, gy, -Grd.h, -[100 100], 'color', ones(1,3)*0.6);
contour3(gx, gy, -Grd.h, -[200 200], 'color', ones(1,3)*0.4);

% Tracklines

h.track(1) = patch(px(:, isak), py(:, isak), pz(:, isak), 0);
h.track(2) = patch(px(:,~isak), py(:,~isak), pz(:,~isak), 0);

h.end(1) = plot3(px(end-1, isak), py(end-1, isak), pz(end-1, isak), '.');
h.end(2) = plot3(px(end-1,~isak), py(end-1,~isak), pz(end-1,~isak), '.');

xtmp = double(px(1,isak)); ytmp = double(py(1,isak));
k1 = convhull(xtmp, ytmp);
h.start(1) = plot(xtmp(k1), ytmp(k1));
xtmp = double(px(1,~isak)); ytmp = double(py(1,~isak));
k2 = convhull(xtmp, ytmp);
h.start(2) = plot(xtmp(k2), ytmp(k2));

set(h.track(1), 'cdata', pc(:, isak), 'edgecolor','interp','facecolor','none', 'edgealpha', 0.1);
set(h.track(2), 'cdata', pc(:,~isak)*-1, 'edgecolor','interp','facecolor','none', 'edgealpha', 0.1);

cmap = crameri('bam');

% set(h.track(1), 'edgecolor', 'r', 'edgealpha', 0.1, 'facecolor', 'none');
% set(h.track(2), 'edgecolor', 'b', 'edgealpha', 0.1, 'facecolor', 'none');

% Prettify

set(h.ax, 'xlim', minmax(px, 'expand', 0.05), 'ylim', minmax(py(py>6e6),'expand', 0.08), 'zlim', [-1000 5], ...
    'dataaspectratio', [1 1 1/1000], ...
    'clim', [-366 366], 'colormap', cmap, ...
    'xtick', [], 'ytick', [], 'box', 'on');
set(h.end(1), 'color', cmap(end,:));
set(h.end(2), 'color', cmap(1,:));
set(h.start(1), 'color', cmap(end,:));
set(h.start(2), 'color', cmap(1,:));

h.cb(1) = colorbar('south');
h.cb.Position(3) = 0.2;
h.cb.Position(2) = 0.25;
set(h.cb, 'ticks', [-365 0 365], 'ticklabels', ["365" "0" "365"], 'tickdir', 'out', 'fontsize', 6, 'axisloc', 'out');
xlabel(h.cb, 'Northwest \leftarrow Days \rightarrow AK Peninsula');

% Grid lines and labels

ltgrid = 45:5:65;
lngrid = 160:10:210;

ltpts = latlimbig(1):latlimbig(2);
lnpts = lonlimbig(1):lonlimbig(2);

for iln = 1:length(lngrid)
    [xg,yg] = projfwd(mproj, ltpts, lngrid(iln).*ones(size(ltpts)));
    h.gln(iln) = patch([xg NaN], [yg NaN], 'r'); %, 'edgecolor', rgb('grey'), 'linestyle', ':', 'alpha', 0.3);
end
for ilt = 1:length(ltgrid)
    [xg,yg] = projfwd(mproj, ltgrid(ilt).*ones(size(lnpts)), lnpts);
    h.glt(ilt) = patch([xg NaN], [yg NaN], 'r'); %, 'edgecolor', rgb('grey'), 'linestyle', ':', 'alpha', 0.3);
end
set([h.gln h.glt],  'edgecolor', rgb('gray'), 'linestyle', ':', 'edgealpha', 0.5);

text(arrayfun(@(x) mean(x.XData(10:11)), h.gln), ...
     arrayfun(@(x) mean(x.YData(10:11)), h.gln), ...
     angl2str(wrapTo180(lngrid), 'ew', 'degrees', 0), ...
     'fontsize', 8, 'color', rgb('gray'), 'clipping', 'on', 'horiz', 'right');
text(arrayfun(@(x) x.XData(32), h.glt), ...
     arrayfun(@(x) x.YData(32), h.glt), ...
     angl2str(wrapTo180(ltgrid), 'ns', 'degrees', 0), ...
     'fontsize', 8, 'color', rgb('gray'), 'clipping', 'on', 'vert', 'base');



if exportfigs
    export_fig('particletracks', h.fig, '-png', '-r300');
end

%% ... Alternate colors (for analysis only, not for paper)

% Just by time, same color scheme for both

set(h.track(1), 'cdata', pc(:, isak), 'edgecolor','interp','facecolor','none', 'edgealpha', 0.1);
set(h.track(2), 'cdata', pc(:,~isak), 'edgecolor','interp','facecolor','none', 'edgealpha', 0.1);
cmocean('thermal');
set(h.ax, 'clim', [0 365]);

% By depth

set(h.track(1), 'cdata', pz(:, isak), 'edgecolor','interp','facecolor','none', 'edgealpha', 0.1);
set(h.track(2), 'cdata', pz(:,~isak), 'edgecolor','interp','facecolor','none', 'edgealpha', 0.1);
cmocean('thermal');
set(h.ax, 'clim', [-500 0]);

% isolate from west ones

iswest = px(end-1,:) < xy(1);
set(h.track(1), 'cdata', double(iswest( isak).*ones(npt+1,1)), 'edgecolor','interp','facecolor','none', 'edgealpha', 0.1);
set(h.track(2), 'cdata', double(iswest(~isak).*ones(npt+1,1)), 'edgecolor','interp','facecolor','none', 'edgealpha', 0.1);
set(h.ax, 'colormap', rgb(["red";"blue"]), 'clim', [-0.5 1.5]);

%% Pelland data distribution supplement

pltmask = Pel.t < datetime(2022,1,1);

h = plotgrid('size', [1 2], 'mar', 0.02, 'mb', 0.1, 'mt', 0.1, 'sp', 0.08, 'ml', 0.03);
set(h.fig, 'color', 'w');

axes(h.ax(1));
h.bax = boxworldmap(rlatlim, rlonlim, 'latgrid', 50:5:60, 'longrid', -180:5:-160);

m = getm(h.ax(1));

[bx,by] = projfwd(m, Border.poly.Vertices(:,2), Border.poly.Vertices(:,1));
plot(polyshape(bx,by), 'facecolor', rgb('gray'), 'edgecolor', 'none');

[sx,sy] = projfwd(m, Pel.LATITUDE(pltmask), Pel.LONGITUDE(pltmask));
h.s = scatter(sx, sy, 5, 'r', 'filled');
set(h.s, 'markerfacealpha', 0.1);

[ndata,yr,tmid] = reshapetimeseries(Pel.t(pltmask), Pel.t(pltmask), 'bin', 'month', 'func', @numel);
ndata(isnan(ndata)) = 0;

axes(h.ax(2));

[mn,yr] = ndgrid(month(tmid), yr);

h.s(2) = scatter(h.ax(2), mn(:), yr(:), 10, ndata(:), 'filled');
crameri('bilbao');
xlabel('Month');
ylabel('Year');
set(h.ax(2), 'ylim', [1970 2023], 'fontsize', 8);
h.cb = colorbar('southoutside');
xlabel(h.cb, 'Number of data points');

set([h.bax.lblpar; h.bax.lblmer], 'fontsize', 8);


if exportfigs
    exportgraphics(h.fig, 'pelland_data_density.png', 'resolution', 150);
end

%% Bathymetry mismatch supplement

[xetop,yetop] = ndgrid(Etop.lon, Etop.lat);
fetop = griddedInterpolant(xetop, yetop, -Etop.z); %, Grd.lon_rho(:), Grd.lat_rho(:));

hetop = fetop(Grd.lon_rho, Grd.lat_rho);
hetop(Grd.mask_rho == 0) = NaN;

h = plotgrid('size', [3 1], 'sv', 0.1);
h.ax(1) = mergeaxes(h.ax(1:2));
h.ax = h.ax([1 3]);
h.fig.Position(3:4) = [385 420]*1.1;

% h = plotgrid('size', [1 1]);
% h.ax = subgridaxes(h.ax, [0.2 0.8], 1);

axes(h.ax(1));
plotromsrho(Grd, Grd.h-hetop);
set(h.ax(1), 'clim', [-1500 1500], 'colormap', cmocean('balance', 30));

bordersm('alaska', 'facecolor', rgb('gray'), 'edgecolor', 'none');
bordersm('russia', 'facecolor', rgb('gray'), 'edgecolor', 'none');

h.cb = colorbar('east');
ylabel(h.cb, 'Bathymetry mismatch (m)');

pos = h.ax(1).Position;
h.cb.Position(1) = pos(1);
h.cb.Position(4) = h.cb.Position(4)*0.5;
h.cb.Position(2) = pos(2)+pos(4)-h.cb.Position(4);
set(h.cb, 'axisloc', 'out');

setm(h.ax(1), 'flinewidth', 1, 'fontsize', 8)

% slice bathymetry vs real world

ltslice = Grd.lat_rho(:,89);
lnslice = Grd.lon_rho(:,89);

plotm(ltslice, lnslice, ':k');

x = dfun([ltslice(1:end-1) lnslice(1:end-1)], [ltslice(2:end) lnslice(2:end)]);
x = cumsum([0; x]);

h.bathy = plot(h.ax(2), x, Grd.h(:,89), ...
                        x, fetop(lnslice,ltslice));
set(h.bathy(1), 'color', Col.b10k);
set(h.bathy(2), 'color', 'k');
set(h.ax(2), 'ydir', 'reverse', 'box', 'off');
xlabel(h.ax(2), 'Distance along transect (km)');
ylabel(h.ax(2), 'Depth (m)');

legend(h.ax(2), 'Bering10K', 'ETOPO 2022', 'location', 'best');
legend(h.ax(2), 'boxoff');

axis(h.ax(2), 'tight');
set(h.fig, 'color', 'w');

% Fix meridian labels

hmer = findobj(h.ax(1), 'tag', 'MLabel');
set(hmer, {'string'}, cellstr(angl2str(wrapTo180(-150:-10:-200), 'ew', 'degrees', 0)))



if exportfigs
    exportgraphics(gcf, 'bathymismatch.png', 'resolution', 150);
end
