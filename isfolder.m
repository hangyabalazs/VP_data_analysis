function results = isfolder(foldername)

results = exist(foldername,'dir') ~= 0;