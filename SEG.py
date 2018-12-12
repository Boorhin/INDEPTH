import  os
Seg =  ["TIB01_56-175.sgy", "TIB01_176-275.sgy", "TIB01_276-375.sgy", "TIB01_376-475.sgy", "TIB01_476-581.sgy"]
EBDIC = 3200
Reel  = 400

with open(filename+'mod', 'wb') as new:
    with open(Seg[0], 'rb') as fo:
        new.write(fo.read())
    for filename in Seg[1:]:
        with open(filename, 'rb') as f:
            f.seek(EBDIC+Reel)
            new.write(f.read())
