function out = ExtractLickingBehavior(ThorSyncFolder)

ThorFile = fullfile(ThorSyncFolder, 'Episode001.h5');

x = h5info(ThorFile);
out = {x.Groups(1).Datasets.Name};