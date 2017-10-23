pkg load image

a = imread('E:\images\lenna_erosion.bmp');
r = 2;
[X, Y] = meshgrid(-r:r, -r:r);
se = X.^2 + Y.^2;
se = se <= r^2;
ea = imdilate(a, se);


imshow(ea)
imshow(a)
imshow(ea-a)
imwrite(mat2gray(ea), 'E:\images\octaveDilate.bmp');