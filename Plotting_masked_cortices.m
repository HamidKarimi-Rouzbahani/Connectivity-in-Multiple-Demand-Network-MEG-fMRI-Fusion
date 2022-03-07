clc;
clear all;
close all;

mask = niftiread('D:\Hamid\Postdoc\MD\Analyses\DuncanOwen2000_MDs\left_dorsal_lateral_PFc.nii');
mask2 = niftiread('D:\Hamid\Postdoc\MD\Analyses\spm12\canonical\avg152T1.nii');
mask3 = load('D:\Hamid\Postdoc\MD\Analyses\DuncanOwen2000_MDs\left_dorsal_lateral_PFc_roi.mat');
surf = gifti('D:\Hamid\Postdoc\MD\Analyses\spm12\canonical\cortex_5124.surf.gii');
% figure; plot(surf,mask);

% mask4 = niftiread('D:\Hamid\Postdoc\MD\Analyses\Results_temp_playing\dlpfc.nii');
% function createfigure(FaceVertexAlphaData1, Vertices1, Faces1, FaceVertexCData1)
% CREATEFIGURE(FaceVertexAlphaData1, Vertices1, Faces1, FaceVertexCData1)

Vertices=surf.vertices;
Faces1=surf.faces;
% Faces1=[surf.faces(1:5124,:);surf.faces(1:5116,:)];
% FaceVertexCData1=surf.vertices;
% face_color=ones(5124,3)*0.5;
% face_color=randi([1 64],[5124 3])./64;
face_color=ones(5124,3)*0.5;

% Create figure
figure1 = figure('InvertHardcopy','off','PaperUnits','normalized',...
    'Tag','Graphics',...
    'NumberTitle','off',...
    'Name','SPM12 (7771): Graphics',...
    'Color',[1 1 1]);
colormap(hot);

% Create axes
axes1 = axes('Tag','SPMMeshRenderBackground','Parent',figure1,...
    'Position',[-0.05 -0.05 1.05 0.555]);

% Set the remaining axes properties
set(axes1,'XTick',[],'YTick',[]);
% Create axes
axes2 = axes('Parent',figure1,'Position',[0.05 0.05 0.9 0.4]);
axis off

% Create patch
patch('Tag','SPMMeshRender','Parent',axes2,'FaceLighting','gouraud',...
    'Clipping','off',...
    'Vertices',Vertices,...
    'SpecularStrength',0,...
    'DiffuseStrength',0.8,...
    'Faces',Faces1,...
    'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceVertexCData',face_color);

% Create light
%     'Position',[820.981244982541 448.725194159879 828.276703732705],...
light('Parent',axes2,...
    'Position',[820.981244982541 800.725194159879 828.276703732705],...
    'Style','local');

view(axes2,[84.8259087175903 12.6760869565219]);
axis(axes2,'tight');
% Set the remaining axes properties
set(axes2,'DataAspectRatio',[1 1 1]);

