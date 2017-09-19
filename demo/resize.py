import cv2
import os
import sys

if len(sys.argv) != 4:
    print(sys.argv[0], '<source file> <dst file> <LxW>')
    print("E.g. ", sys.argv[0], 'test1.png test2.png 80x80')
    sys.exit(1)

src_filename = sys.argv[1];

if src_filename.endswith(".png"):

    img = cv2.imread(src_filename)
    img_size = sys.argv[3]
    index = img_size.find('x');
    if index == -1:
        print("you have wrong format for specifying size")
        sys.exit(1)
    
    img_length = int(img_size[:index])
    img_width = int(img_size[index+1:])
    
    
    img = cv2.resize(img, (img_length, img_width))
    
    dst_filename = sys.argv[2]
    cv2.imwrite(dst_filename, img)
