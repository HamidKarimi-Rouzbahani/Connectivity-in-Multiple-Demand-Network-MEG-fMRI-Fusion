% =========================================================================
function files = get_files(direc, filt)
% =========================================================================
% return a list of files
% filt = filter string
% direc = cell array of directory names

files = [];
%files = get_files(char(my_epi),'f*.nii');
for d=1:size(direc,1) % loop through each EPI session
    tmp = dir(fullfile(direc(d,:),filt)); % find all files matching f*.nii
    tmp = [repmat([direc(d,:) filesep],size(tmp,1),1) char(tmp.name)]; % build the full path name for these files
    files = char(files,tmp);
end

files = files(~all(files'==' ')',:);
end