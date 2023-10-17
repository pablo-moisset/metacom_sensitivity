function [t, DParches, xP, yP] = gen_landscape(patch_qty, hurst, resolution)

%% La siguiente línea genera un paisaje fractal:

[z , ~, ~] = artificial_surf(1.0, hurst, 0.1, resolution, resolution);

%% La siguiente línea define un paisaje plano con equiprobabilidad en la colocación de todos los parches.

%% Preparacion de paisaje.
z = z - min(z(:)) ;
mat = 1 - z./max(z(:)) ;
mat_org = mat ;

%mat = 1000.^mat/1000 ;
mat = (tanh(10*(mat-0.5))+1)/2 ;
t = mat/sum(mat(:)) ;

mt = t(:) ;

bitmap_landscape = zeros(1,resolution*resolution) ;

chosen = 0 ;

while patch_qty - chosen > 0 
%    chosen
    R = mnrnd(patch_qty - chosen, mt);
    bitmap_landscape = (bitmap_landscape + R)>0 ;
    chosen = sum(bitmap_landscape) ; 
end

[xP,yP] = ind2sub([resolution,resolution],find(bitmap_landscape)) ;

% crear matriz de distancias:
DParches = squareform(pdist([xP',yP'],'euclidean')) ; 

end
