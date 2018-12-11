import segyio, os
Seg =  ["TIB01_56-175.sgy", "TIB01_176-275.sgy", "TIB01_276-375.sgy", "TIB01_376-475.sgy", "TIB01_476-581.sgy"]
EBDIC = 3200
Reel  = 400

#TIB01_A,TIB01_B, TIB01_C, TIB01_D, TIB01_E = segyio.open(Seg[0], ignore_geometry=True), segyio.open(Seg[1], ignore_geometry=True),segyio.open(Seg[2], ignore_geometry=True),segyio.open(Seg[3], ignore_geometry=True),segyio.open(Seg[4], ignore_geometry=True)

for filename in Seg[1:]:
    with open(filename, 'rb') as f:
        Total = os.path.getsize(filename)
        f.seek(EBDIC+Reel)
        new =open(filename+'mod', 'wb')
        new.write(f.read())
                     
