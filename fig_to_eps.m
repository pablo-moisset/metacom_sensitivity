f = dir('*.fig') ;

for i=1:length(f)
    myfile = f(i) ;
    openfig([myfile.folder '/' myfile.name])
    saveas(gcf,[myfile.folder '/' myfile.name(1:end-3) 'eps'],'epsc')
end
