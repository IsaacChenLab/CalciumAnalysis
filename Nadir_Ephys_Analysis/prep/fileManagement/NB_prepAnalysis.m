function [] = NB_prepAnalysis(numRecordings)

NB_createFolders(numRecordings);
disp 'Folders created'
NB_sortFiles(numRecordings);
disp 'Files sorted'
NB_batchRenameFiles(numRecordings);
disp 'Files renamed, ready for analysis'