# On Ubuntu, this requires the teem-apps package
teem-unu make -h -i dummy.raw -t uchar -s 5 6 > headers/test2duchar.nhdr
teem-unu make -h -i dummy.raw -t ushort -s 5 6 > headers/test2dushort.nhdr
teem-unu make -h -i dummy.raw -t float -s 5 6 > headers/test2dfloat.nhdr
teem-unu make -h -i dummy.raw -t uchar -s 3 5 6 -k "RGB-color" "space" "space" > headers/test2drgb.nhdr
teem-unu make -h -i dummy.raw -t uchar -s 5 6 3 -k "space" "space" "RGB-color" > headers/test2drgb_last.nhdr
teem-unu make -h -i dummy.raw -t uchar -s 3 5 6 -k "RGB-color" "space" "space" -spc 2 -spu "mm" "μm" --spacing nan 2 3.5 > headers/test2drgb_spu_spacing.nhdr
teem-unu make -h -i dummy.raw -t uchar -s 5 6 -k "space" "space" -spc 2 -spu "mm" "mm" --spacing 2 3 --label "first" "second" > headers/test2d_spu_spacing_labels.nhdr
teem-unu make -h -i dummy.raw -t uchar -s 3 5 6 -k "RGB-color" "space" "space" -spc 2 -spu "mm" "mm" --spacing nan 2 3 --label "rgb" "first" "second" > headers/test2drgb_spu_spacing_labels.nhdr
teem-unu make -h -i dummy.raw -t uchar -s 3 5 6 -k "RGB-color" "space" "space" -spc 2 -spu "mm" "mm" -dirs "none (1,0) (0,1)" > headers/test2drgb_spu_sd.nhdr

teem-unu make -h -i dummy.raw -t uchar -s 3 3 3 -k space space space -spc RAS -spu "mm" "mm" "mm" --spacing 2 2 5 -orig "(10,12,100)" > headers/test3d_units_range_axes.nhdr
teem-unu make -h -i dummy.raw -t uchar -s 3 3 3 2 -k space space space list -spc RAS -orig "(1,2,3)" > headers/test3dlist_origin.nhdr

teem-unu make -h -i dummy.raw -t ushort -s 8 10 3 800 -k space space space time -spc 3 > headers/test3d_time.nhdr
teem-unu make -h -i dummy.raw -t ushort -s 8 10 3 800 -k space space space time -spc 3 -l "x" "l" "s" "t"  > headers/test3d_time_labels.nhdr
teem-unu make -h -i dummy.raw -t ushort -s 8 10 3 800 -k space space space time -spc 3 -l "x" "l" "s" "time" -spu "μm" "μm" "μm" --spacing 0.25 0.25 2 nan > headers/test3d_time_labels_units.nhdr
teem-unu make -h -i dummy.raw -t ushort -s 8 10 3 800 -k space space space time -spc 3 -l "x" "l" "s" "time" -u "μm" "μm" "μm" "s" --spacing 0.25 0.25 2 0.8 > headers/test3d_time_labels_units_s.nhdr
