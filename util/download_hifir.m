function path_to_hifir = download_hifir(version)
% Download hifir

zip_file = fullfile(tempdir, 'hifir.zip');
zip_url = ['https://github.com/hifirworks/hifir/archive/refs/tags/v' ...
    version ...
    '.zip'];
urlwrite(zip_url, zip_file);
unzip(zip_file, tempdir);
delete(zip_file);
path_to_hifir = fullfile(tempdir, ['hifir-' version]);
end