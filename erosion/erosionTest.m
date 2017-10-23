pkg load image

a = imread('E:\images\lenna_biLevel.bmp');
r = 2;
[X, Y] = meshgrid(-r:r, -r:r);
se = X.^2 + Y.^2;
se = se <= r^2;
ea = imerode(a, se);


imshow(ea)
imshow(a)
imshow(a-ea)
imwrite(mat2gray(ea), 'E:\images\octaveErode.bmp');