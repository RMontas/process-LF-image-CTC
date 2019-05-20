function img_4DLF_VIEWS = deconstruct_lenslet_img10( img_4DLF_MI, mi_size )

siz = size(img_4DLF_MI);
height = (siz(1) - mod(siz(1),mi_size)) / mi_size;
width = (siz(2) - mod(siz(2),mi_size)) / mi_size;
    
img_4DLF_VIEWS = uint16(zeros(mi_size, mi_size, height, width, 3));

for mi = 1:mi_size
    for mj = 1:mi_size
        aux = img_4DLF_MI(mj + (0:height-1)*mi_size, mi + (0:width-1)*mi_size,:);
       % for i = 1:width
         %   for j = 1:height
               % img_4DLF_VIEWS(mj,mi,j,i,:) = img_4DLF_MI(mj + (j-1)*mi_size, mi + (i-1)*mi_size,:);
        %    end
        %end
        img_4DLF_VIEWS(mj,mi,:,:,:) = aux;
    end
end
