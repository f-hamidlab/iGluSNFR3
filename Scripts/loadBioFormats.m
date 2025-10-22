%% [im_data, ops] = loadBioFormats(ops, varargin)
% Loads data using bioformats toolbox and reshapes data.
%
% Author: Jane Ling (yan_to.ling@kcl.ac.uk)
% Last updated: 2025-02-14 12:40
%
% USAGE:
% 1) [im_data, ops] = loadBioFormats(ops);
% 2) [im_data, ops] = loadBioFormats(ops, verbose);
%
% INPUTS:
%   - ops 
%       (struct) ops.filename defining the path to image file
%   - verbose 
%       (bool) whether to pring information about dataset or not.
%       (Default) false
%
% OUTPUTS:
%   - im_data
%       (numeric) 5D matrix if dimension order starts with 'XY', singular
%                 dimensions would be omitted.
%       (cell array) interleaved frames otherwise
%   - ops
%       (struct) parameters of the stack

function [im_data, ops] = loadBioFormats(ops, varargin)
    
    if nargin == 2
        verbose = varargin{1};
    else 
        verbose = 0;
    end

    % Bio-Formats Toolbox
    [~,~] = bfCheckJavaPath();   % added such that path to Bio-Formats Toolbox is know
    
    data = bfopen(ops.filename); % load data

    im_data = data{1,1};
    im_data = im_data(:,1); % image data
    
    omeMeta = data{1,4}; % metadata
    
    ops.imageIndex = 0;
    
    % size of stack in each dimension
    ops.Nx = omeMeta.getPixelsSizeX(ops.imageIndex).getValue(); % image width, pixels
    ops.Ny = omeMeta.getPixelsSizeY(ops.imageIndex).getValue(); % image height, pixels
    ops.Nz = omeMeta.getPixelsSizeZ(ops.imageIndex).getValue(); % number of Z slices
    ops.Nc = omeMeta.getPixelsSizeC(ops.imageIndex).getValue(); % number of channels
    ops.Nt = omeMeta.getPixelsSizeT(ops.imageIndex).getValue(); % number of time points
    
    % dimension order in data
    ops.dimOrder = char(omeMeta.getPixelsDimensionOrder(ops.imageIndex).getValue());
    
    % size of pixels (voxels)
    ops.sizeX = double(omeMeta.getPixelsPhysicalSizeX(ops.imageIndex).value(ome.units.UNITS.MICROMETER)); % [microns]
    ops.sizeY = double(omeMeta.getPixelsPhysicalSizeY(ops.imageIndex).value(ome.units.UNITS.MICROMETER)); % [microns]
    
    % axes
    ops.x = (0:ops.Nx-1)*ops.sizeX;
    ops.y = (0:ops.Ny-1)*ops.sizeY;
    
    if verbose
        fprintf('Size of image = %d x %d. \nNumber of Z slices = %d. \nNumber of channels = %d. \nNumber of time points = %d. \n', ops.Nx, ops.Ny, ops.Nz, ops.Nc, ops.Nt)
        fprintf('Dimension order = %s  \n', ops.dimOrder)
    end

    % for z-stack
    if ops.Nz > 1
        try
            ops.sizeZ = double(omeMeta.getPixelsPhysicalSizeZ(ops.imageIndex).value(ome.units.UNITS.MICROMETER)); % [microns]
            ops.z = (0:ops.Nz-1)*ops.sizeZ;
            if verbose
                fprintf('Voxel size = %.3f µm x %.3f µm x %.3f µm  \n', ops.sizeX, ops.sizeY, ops.sizeZ)
            end
        catch
            if verbose
                fprintf('Pixel size = %.3f µm x %.3f µm  \n', ops.sizeX, ops.sizeY)
            end
            warning('Could not read sizeZ from omeMeta.')
        end
    else
        if verbose
            fprintf('Pixel size = %.3f µm x %.3f µm  \n', ops.sizeX, ops.sizeY)
        end
    end
    
    % for time stack
    if ops.Nt > 1
        try
            ops.sizeT = omeMeta.getPixelsTimeIncrement(ops.imageIndex).value(ome.units.UNITS.SECOND);
            ops.t = (0:ops.Nt-1)*ops.sizeT;
            if verbose
                fprintf('Pixel time imcrement = %.3f s  \n', ops.sizeT)
            end
        catch
            warning('Dimension order may be incorrect.')
        end
    end
        
    % reshaping
    if strcmp(ops.dimOrder(1:2),'XY') 
        im_data = cell2mat(im_data);
        im_data = reshape(im_data,ops.Ny,ops.Nc*ops.Nz*ops.Nt,ops.Nx);
        im_data = permute(im_data, [1,3,2]);
        N = [ops.Ny,ops.Nx, 0, 0, 0];
    
        N(strfind(ops.dimOrder,'C')) = ops.Nc;
        N(strfind(ops.dimOrder,'Z')) = ops.Nz;
        N(strfind(ops.dimOrder,'T')) = ops.Nt;
        im_data = reshape(im_data,N(1),N(2),N(3),N(4),N(5));

        im_data = squeeze(im_data);
        ops.dimOrder = ops.dimOrder(N~=1);
        if verbose
            fprintf('Dimension order of im_data after reshaping = %s  \n', ops.dimOrder)
        end
    end
    
end
