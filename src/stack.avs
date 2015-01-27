#clip1 = CoronaSequence("./2015-1-25_0-50_output_elliptic_lense/Ez_00*_rank-0_t0.png", sort=2,textmode=1)
#clip2 = CoronaSequence("./2015-1-26_17-22_output_wo_elliptic_lense/Ez_00*_rank-0_t0.png", sort=2,textmode=1)

clip1 = ImageSource("./2015-1-25_0-50_output_elliptic_lense/Ez_%05d_rank-0_t0.png"    , 1,401, 60)
clip2 = ImageSource("./2015-1-26_17-22_output_wo_elliptic_lense/Ez_%05d_rank-0_t0.png", 1,401, 60)

clip1 = AssumeFPS(clip1, 60).LanczosResize(720,360)
clip2 = AssumeFPS(clip2, 60).LanczosResize(720,360)

StackVertical(clip1,clip2)

#ffmpeg -i stack.avs -c:v libvpx -crf 4 -qmin 0 -qmax 9 -b:v 6M stacked2.webm
