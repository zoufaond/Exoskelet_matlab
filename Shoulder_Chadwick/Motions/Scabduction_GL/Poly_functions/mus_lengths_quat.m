function mus_lengths = mus_lengths_quat(t,q)
qsubs0 = q(1,:);
qsubs1 = q(2,:);
qsubs2 = q(3,:);
qsubs3 = q(4,:);
qsubs4 = q(5,:);
qsubs5 = q(6,:);
qsubs6 = q(7,:);
qsubs7 = q(8,:);
qsubs8 = q(9,:);
qsubs9 = q(10,:);
qsubs10 = q(11,:);
qsubs11 = q(12,:);
qsubs12 = q(13,:);
x0 = qsubs0.*qsubs3;
x1 = 2*x0;
x2 = 2*qsubs1.*qsubs2 - x1;
x3 = qsubs5.^2;
x4 = 2*x3;
x5 = qsubs7.^2;
x6 = 2*x5 - 1;
x7 = -x4 - x6;
x8 = 0.0173650793650794*x7;
x9 = qsubs2.^2;
x10 = 2*x9;
x11 = qsubs3.^2;
x12 = 2*x11;
x13 = -x10 - x12 + 1;
x14 = 2*qsubs4;
x15 = qsubs7.*x14;
x16 = 2*qsubs5.*qsubs6 - x15;
x17 = 0.0173650793650794*x16;
x18 = qsubs0.*qsubs2;
x19 = 2*x18;
x20 = qsubs1.*qsubs3;
x21 = x19 + 2*x20;
x22 = qsubs5.*x14;
x23 = 2*qsubs6;
x24 = qsubs7.*x23 + x22;
x25 = 0.0173650793650794*x24;
x26 = qsubs5.*x23 + x15;
x27 = 0.0518095238095238*x26;
x28 = qsubs6.*x14;
x29 = 2*qsubs5.*qsubs7 - x28;
x30 = 0.0518095238095238*x29;
x31 = qsubs6.^2;
x32 = 2*x31;
x33 = -x32 - x6;
x34 = 0.0518095238095238*x33;
x35 = 2*qsubs6.*qsubs7 - x22;
x36 = 0.0678730158730159*x35;
x37 = -x32 - x4 + 1;
x38 = 0.0678730158730159*x37;
x39 = qsubs5.*qsubs7;
x40 = x28 + 2*x39;
x41 = 0.0678730158730159*x40;
x42 = x11 + x9;
x43 = 2*qsubs1.*qsubs3 - x19;
x44 = qsubs0.*qsubs1;
x45 = 2*x44;
x46 = qsubs2.*qsubs3;
x47 = x45 + 2*x46;
x48 = qsubs1.^2;
x49 = 2*x48 - 1;
x50 = -x10 - x49;
x51 = x18 - x20;
x52 = 2*qsubs2.*qsubs3 - x45;
x53 = qsubs1.*qsubs2;
x54 = x1 + 2*x53;
x55 = -x12 - x49;
x56 = -x0 - x53;
x57 = 0.0223809523809524*x7;
x58 = 0.0223809523809524*x16;
x59 = 0.0223809523809524*x24;
x60 = 0.109873015873016*x26;
x61 = 0.109873015873016*x29;
x62 = 0.109873015873016*x33;
x63 = 0.104952380952381*x35;
x64 = 0.104952380952381*x37;
x65 = 0.104952380952381*x40;
x66 = 0.0315238095238095*x7;
x67 = 0.0315238095238095*x16;
x68 = 0.0315238095238095*x24;
x69 = 0.182444444444444*x26;
x70 = 0.182444444444444*x29;
x71 = 0.182444444444444*x33;
x72 = 0.110380952380952*x35;
x73 = 0.110380952380952*x37;
x74 = 0.110380952380952*x40;
x75 = x42 - 0.593015873015873;
x76 = 0.0346031746031746*x7;
x77 = 0.0346031746031746*x16;
x78 = 0.0346031746031746*x24;
x79 = 0.223492063492063*x26;
x80 = 0.223492063492063*x29;
x81 = 0.223492063492063*x33;
x82 = 0.110984126984127*x35;
x83 = 0.110984126984127*x37;
x84 = 0.110984126984127*x40;
x85 = 0.029047619047619*x7;
x86 = 0.029047619047619*x16;
x87 = 0.029047619047619*x24;
x88 = 0.246253968253968*x26;
x89 = 0.246253968253968*x29;
x90 = 0.246253968253968*x33;
x91 = 0.0925714285714286*x35;
x92 = 0.0925714285714286*x37;
x93 = 0.0925714285714286*x40;
x94 = x13.*x86 + x13.*x90 - x13.*x93 + x2.*x85 + x2.*x88 - x2.*x91 + x21.*x87 + x21.*x89 - x21.*x92 + x42;
x95 = x43.*x86 + x43.*x90 - x43.*x93 + x47.*x85 + x47.*x88 - x47.*x91 + x50.*x87 + x50.*x89 - x50.*x92 + x51;
x96 = x52.*x92;
x97 = x55.*x91;
x98 = x54.*x93;
x99 = 0.221301587301587*x26;
x100 = 0.221301587301587*x29;
x101 = 0.221301587301587*x33;
x102 = 0.0494285714285714*x7;
x103 = 0.0494285714285714*x16;
x104 = 0.0494285714285714*x24;
x105 = 0.122761904761905*x35;
x106 = 0.122761904761905*x37;
x107 = 0.122761904761905*x40;
x108 = x0 + x53;
x109 = 0.24352380952381*x26;
x110 = 0.24352380952381*x29;
x111 = 0.24352380952381*x33;
x112 = 0.0445714285714286*x7;
x113 = 0.0445714285714286*x16;
x114 = 0.0445714285714286*x24;
x115 = 0.116857142857143*x35;
x116 = 0.116857142857143*x37;
x117 = 0.116857142857143*x40;
x118 = 0.302666666666667*x26;
x119 = 0.302666666666667*x29;
x120 = 0.302666666666667*x33;
x121 = 0.0397777777777778*x7;
x122 = 0.0397777777777778*x16;
x123 = 0.0397777777777778*x24;
x124 = 0.0945714285714286*x35;
x125 = 0.0945714285714286*x37;
x126 = 0.0945714285714286*x40;
x127 = 0.0575555555555556*x7;
x128 = 0.0575555555555556*x16;
x129 = 0.0575555555555556*x24;
x130 = 0.287492063492064*x26;
x131 = 0.287492063492064*x29;
x132 = 0.287492063492064*x33;
x133 = 0.108444444444444*x35;
x134 = 0.108444444444444*x37;
x135 = 0.108444444444444*x40;
x136 = 0.0613015873015873*x7;
x137 = 0.0613015873015873*x16;
x138 = 0.0613015873015873*x24;
x139 = 0.267079365079365*x26;
x140 = 0.267079365079365*x29;
x141 = 0.267079365079365*x33;
x142 = 0.11568253968254*x35;
x143 = 0.11568253968254*x37;
x144 = 0.11568253968254*x40;
x145 = -x11 - x9;
x146 = 0.032952380952381*x7;
x147 = 0.032952380952381*x16;
x148 = 0.032952380952381*x24;
x149 = 0.371936507936508*x26;
x150 = 0.371936507936508*x29;
x151 = 0.371936507936508*x33;
x152 = 0.0712063492063492*x35;
x153 = 0.0712063492063492*x37;
x154 = 0.0712063492063492*x40;
x155 = 0.006*x7;
x156 = 0.006*x16;
x157 = 0.006*x24;
x158 = 0.0485396825396825*x35;
x159 = 0.0485396825396825*x37;
x160 = 0.0485396825396825*x40;
x161 = 0.140031746031746*x35;
x162 = 0.140031746031746*x37;
x163 = 0.140031746031746*x40;
x164 = 0.064031746031746*x7;
x165 = 0.064031746031746*x16;
x166 = 0.064031746031746*x24;
x167 = 0.0957142857142857*x26;
x168 = 0.0957142857142857*x29;
x169 = 0.0957142857142857*x33;
x170 = x13.*x163 + x13.*x165 + x13.*x169 + x161.*x2 + x162.*x21 + x164.*x2 + x166.*x21 + x167.*x2 + x168.*x21 + x42;
x171 = x161.*x47 + x162.*x50 + x163.*x43 + x164.*x47 + x165.*x43 + x166.*x50 + x167.*x47 + x168.*x50 + x169.*x43 + x51;
x172 = x161.*x55 + x162.*x52 + x163.*x54 + x164.*x55 + x165.*x54 + x166.*x52 + x167.*x55 + x168.*x52 + x169.*x54 + x56;
x173 = 0.0823174603174603*x26;
x174 = 0.0823174603174603*x29;
x175 = 0.0823174603174603*x33;
x176 = 0.148730158730159*x35;
x177 = 0.148730158730159*x37;
x178 = 0.148730158730159*x40;
x179 = 0.0833968253968254*x7;
x180 = 0.0833968253968254*x16;
x181 = 0.0833968253968254*x24;
x182 = 0.152539682539683*x35;
x183 = 0.152539682539683*x37;
x184 = 0.152539682539683*x40;
x185 = 0.0678730158730159*x26;
x186 = 0.0678730158730159*x29;
x187 = 0.0678730158730159*x33;
x188 = 0.109301587301587*x7;
x189 = 0.109301587301587*x16;
x190 = 0.109301587301587*x24;
x191 = 0.0857460317460317*x7;
x192 = 0.0857460317460317*x16;
x193 = 0.0857460317460317*x24;
x194 = 0.376444444444444*x26;
x195 = 0.376444444444444*x29;
x196 = 0.376444444444444*x33;
x197 = 0.0764444444444444*x35;
x198 = 0.0764444444444444*x37;
x199 = 0.0764444444444444*x40;
x200 = 0.311238095238095*x7;
x201 = 0.311238095238095*x16;
x202 = 0.311238095238095*x24;
x203 = 0.405269841269841*x26;
x204 = 0.405269841269841*x29;
x205 = 0.405269841269841*x33;
x206 = 0.0955238095238095*x35;
x207 = 0.0955238095238095*x37;
x208 = 0.0955238095238095*x40;
x209 = 0.328222222222222*x7;
x210 = 0.328222222222222*x16;
x211 = 0.328222222222222*x24;
x212 = 0.409301587301587*x26;
x213 = 0.409301587301587*x29;
x214 = 0.409301587301587*x33;
x215 = 0.0973968253968254*x35;
x216 = 0.0973968253968254*x37;
x217 = 0.0973968253968254*x40;
x218 = 0.366634920634921*x7;
x219 = 0.366634920634921*x16;
x220 = 0.366634920634921*x24;
x221 = 0.408031746031746*x26;
x222 = 0.408031746031746*x29;
x223 = 0.408031746031746*x33;
x224 = 0.098*x35;
x225 = 0.098*x37;
x226 = 0.098*x40;
x227 = 0.404571428571429*x26;
x228 = 0.404571428571429*x29;
x229 = 0.404571428571429*x33;
x230 = 0.404*x7;
x231 = 0.404*x16;
x232 = 0.404*x24;
x233 = 0.0912063492063492*x35;
x234 = 0.0912063492063492*x37;
x235 = 0.0912063492063492*x40;
x236 = qsubs2.*qsubs6;
x237 = qsubs6.*qsubs7;
x238 = qsubs1.*qsubs6;
x239 = qsubs6.^3;
x240 = qsubs7.*x239;
x241 = qsubs6.^4;
x242 = qsubs7.^4;
x243 = qsubs3.*qsubs7;
x244 = qsubs1.*x5;
x245 = qsubs5.*x31;
x246 = qsubs6.*x9;
x247 = x31.*x5;
x248 = qsubs1.^4;
x249 = qsubs7.^3;
x250 = qsubs1.*x249;
x251 = qsubs1.^3.*qsubs5;
x252 = qsubs2.^3;
x253 = qsubs3.*x252;
x254 = x243.*x9;
x255 = qsubs5.*qsubs6;
x256 = x255.*x5;
x257 = qsubs3.^3;
x258 = qsubs2.*x3;
x259 = qsubs10.*qsubs9;
x260 = qsubs11.*qsubs9;
x261 = qsubs10.^2;
x262 = qsubs9.^2;
x263 = qsubs10.*qsubs11;
x264 = qsubs11.^2;
x265 = qsubs11.^4;
x266 = qsubs10.*x264;
x267 = qsubs11.^3;
x268 = qsubs10.*x39;
x269 = qsubs10.*x31;
x270 = qsubs9.*x239;
x271 = x259.*x31;
x272 = qsubs9.*x261;
x273 = x260.*x31;
x274 = x261.*x31;
x275 = x264.*x5;
x276 = qsubs11.*qsubs6;
x277 = x276.*x5;
x278 = qsubs10.*x262;
x279 = qsubs9.^4;
x280 = qsubs11.*x261;
x281 = qsubs9.*x264;
x282 = x262.*x264;
x283 = qsubs10.*x267;
x284 = qsubs9.^3;
x285 = qsubs10.^3;
x286 = qsubs9.*x285;
x287 = x259.*x264;
x288 = qsubs12.*x260;
x289 = qsubs12.^2;
x290 = qsubs10.*x289;
x291 = qsubs9.*x289;
x292 = qsubs12.*qsubs9;
x293 = qsubs10.*qsubs12.^3;
x294 = qsubs11.*x262;
x295 = 0.158393294527867*qsubs11.*qsubs12 - 1.04542368081266*x260.*x289 + 0.228930997847512*x260 - 0.106942315781228*x290 + 0.380525055710162*x293 - 0.238965465949944*x294 + 0.298575760721053;
x296 = qsubs10.*qsubs7;
x297 = qsubs11.*qsubs7;
x298 = qsubs1.*qsubs5;
x299 = qsubs10.*x9;
x300 = qsubs6.*x252;
x301 = qsubs10.*qsubs2;
x302 = qsubs11.*x5;
x303 = x261.*x39;
x304 = qsubs1.*qsubs11;
x305 = qsubs1.*qsubs10;
x306 = qsubs1.*qsubs9.*x3;
x307 = x237.*x261;
x308 = x236.*x297;
x309 = qsubs2.*qsubs7;
x310 = qsubs11.*x309;
x311 = qsubs2.*x269;
x312 = qsubs5.*x9;
x313 = qsubs10.*x239;
x314 = qsubs1.*x269;
x315 = qsubs11.*x312;
x316 = x255.*x9;
mus_lengths = [0.315*sqrt((x13.*x17 + x13.*x34 - x13.*x41 + x2.*x27 - x2.*x36 + x2.*x8 + x21.*x25 + x21.*x30 - x21.*x38 + x42 - 0.605079365079365).^2 + (x17.*x43 + x25.*x50 + x27.*x47 + x30.*x50 + x34.*x43 - x36.*x47 - x38.*x50 - x41.*x43 + x47.*x8 + x51 + 0.342857142857143).^2 + (x17.*x54 + x25.*x52 + x27.*x55 + x30.*x52 + x34.*x54 - x36.*x55 - x38.*x52 - x41.*x54 + x55.*x8 + x56 + 0.147936507936508).^2); 0.315*sqrt((x13.*x58 + x13.*x62 - x13.*x65 + x2.*x57 + x2.*x60 - x2.*x63 + x21.*x59 + x21.*x61 - x21.*x64 + x42 - 0.604761904761905).^2 + (x43.*x58 + x43.*x62 - x43.*x65 + x47.*x57 + x47.*x60 - x47.*x63 + x50.*x59 + x50.*x61 - x50.*x64 + x51 + 0.346666666666667).^2 + (x52.*x59 + x52.*x61 - x52.*x64 + x54.*x58 + x54.*x62 - x54.*x65 + x55.*x57 + x55.*x60 - x55.*x63 + x56 + 0.138730158730159).^2); 0.315*sqrt((x13.*x67 + x13.*x71 - x13.*x74 + x2.*x66 + x2.*x69 - x2.*x72 + x21.*x68 + x21.*x70 - x21.*x73 + x75).^2 + (x43.*x67 + x43.*x71 - x43.*x74 + x47.*x66 + x47.*x69 - x47.*x72 + x50.*x68 + x50.*x70 - x50.*x73 + x51 + 0.354603174603175).^2 + (x52.*x68 + x52.*x70 - x52.*x73 + x54.*x67 + x54.*x71 - x54.*x74 + x55.*x66 + x55.*x69 - x55.*x72 + x56 + 0.113333333333333).^2); 0.315*sqrt((x13.*x77 + x13.*x81 - x13.*x84 + x2.*x76 + x2.*x79 - x2.*x82 + x21.*x78 + x21.*x80 - x21.*x83 + x42 - 0.591111111111111).^2 + (x43.*x77 + x43.*x81 - x43.*x84 + x47.*x76 + x47.*x79 - x47.*x82 + x50.*x78 + x50.*x80 - x50.*x83 + x51 + 0.366031746031746).^2 + (x52.*x78 + x52.*x80 - x52.*x83 + x54.*x77 + x54.*x81 - x54.*x84 + x55.*x76 + x55.*x79 - x55.*x82 + x56 + 0.0965079365079365).^2); 0.315*sqrt((x94 - 0.589206349206349).^2 + (x95 + 0.377460317460317).^2 + (x52.*x87 + x52.*x89 + x54.*x86 + x54.*x90 + x55.*x85 + x55.*x88 + x56 - x96 - x97 - x98 + 0.0612698412698413).^2); 0.315*sqrt((x100.*x21 + x101.*x13 + x102.*x2 + x103.*x13 + x104.*x21 - x105.*x2 - x106.*x21 - x107.*x13 + x2.*x99 + x75).^2 + (x100.*x50 + x101.*x43 + x102.*x47 + x103.*x43 + x104.*x50 - x105.*x47 - x106.*x50 - x107.*x43 + x47.*x99 + x51 + 0.413333333333333).^2 + (-x105.*x55 - x106.*x52 - x107.*x54 - x108 + 0.0494285714285714*x16.*x54 + 0.0494285714285714*x24.*x52 + 0.221301587301587*x26.*x55 + 0.221301587301587*x29.*x52 + 0.221301587301587*x33.*x54 + 0.0494285714285714*x55.*x7 - 0.0222222222222222).^2); 0.315*sqrt((x94 - 0.584444444444444).^2 + (x95 + 0.446984126984127).^2 + (-x108 + 0.029047619047619*x16.*x54 + 0.029047619047619*x24.*x52 + 0.246253968253968*x26.*x55 + 0.246253968253968*x29.*x52 + 0.246253968253968*x33.*x54 + 0.029047619047619*x55.*x7 - x96 - x97 - x98 - 0.129206349206349).^2); 0.315*sqrt((-x108 - x115.*x55 - x116.*x52 - x117.*x54 + 0.0445714285714286*x16.*x54 + 0.0445714285714286*x24.*x52 + 0.24352380952381*x26.*x55 + 0.24352380952381*x29.*x52 + 0.24352380952381*x33.*x54 + 0.0445714285714286*x55.*x7 - 0.221269841269841).^2 + (x109.*x2 + x110.*x21 + x111.*x13 + x112.*x2 + x113.*x13 + x114.*x21 - x115.*x2 - x116.*x21 - x117.*x13 + x42 - 0.567619047619048).^2 + (x109.*x47 + x110.*x50 + x111.*x43 + x112.*x47 + x113.*x43 + x114.*x50 - x115.*x47 - x116.*x50 - x117.*x43 + x51 + 0.46031746031746).^2); 0.315*sqrt((-x108 - x124.*x55 - x125.*x52 - x126.*x54 + 0.0397777777777778*x16.*x54 + 0.0397777777777778*x24.*x52 + 0.302666666666667*x26.*x55 + 0.302666666666667*x29.*x52 + 0.302666666666667*x33.*x54 + 0.0397777777777778*x55.*x7 - 0.373650793650794).^2 + (x118.*x2 + x119.*x21 + x120.*x13 + x121.*x2 + x122.*x13 + x123.*x21 - x124.*x2 - x125.*x21 - x126.*x13 + x42 - 0.558412698412698).^2 + (x118.*x47 + x119.*x50 + x120.*x43 + x121.*x47 + x122.*x43 + x123.*x50 - x124.*x47 - x125.*x50 - x126.*x43 + x51 + 0.463809523809524).^2); 0.315*sqrt((-x108 - x133.*x55 - x134.*x52 - x135.*x54 + 0.0575555555555556*x16.*x54 + 0.0575555555555556*x24.*x52 + 0.287492063492064*x26.*x55 + 0.287492063492064*x29.*x52 + 0.287492063492064*x33.*x54 + 0.0575555555555556*x55.*x7 - 0.491111111111111).^2 + (x127.*x2 + x128.*x13 + x129.*x21 + x13.*x132 - x13.*x135 + x130.*x2 + x131.*x21 - x133.*x2 - x134.*x21 + x42 - 0.546984126984127).^2 + (x127.*x47 + x128.*x43 + x129.*x50 + x130.*x47 + x131.*x50 + x132.*x43 - x133.*x47 - x134.*x50 - x135.*x43 + x51 + 0.46952380952381).^2); 0.315*sqrt((-x108 - x142.*x55 - x143.*x52 - x144.*x54 + 0.0613015873015873*x16.*x54 + 0.0613015873015873*x24.*x52 + 0.267079365079365*x26.*x55 + 0.267079365079365*x29.*x52 + 0.267079365079365*x33.*x54 + 0.0613015873015873*x55.*x7 - 0.591111111111111).^2 + (x13.*x137 + x13.*x141 - x13.*x144 + x136.*x2 + x138.*x21 + x139.*x2 + x140.*x21 - x142.*x2 - x143.*x21 + x42 - 0.536507936507936).^2 + (x136.*x47 + x137.*x43 + x138.*x50 + x139.*x47 + x140.*x50 + x141.*x43 - x142.*x47 - x143.*x50 - x144.*x43 + x51 + 0.44952380952381).^2); 0.22002*sqrt((0.044450504499591*x11 + 0.0798109262794291*x44 - 0.0798109262794291*x46 + 0.044450504499591*x48 + x56 + 0.534996818471048).^2 + (-0.044450504499591*x44 - 0.044450504499591*x46 + 0.0798109262794291*x48 + x51 + 0.0798109262794291*x9 + 0.443232433415144).^2 + (0.044450504499591*qsubs0.*qsubs3 - x145 - 0.0798109262794291*x18 - 0.0798109262794291*x20 - 0.044450504499591*x53 - 0.65771293518771).^2); 0.27294*sqrt((0.0210302630614787*x11 + 0.064776141276471*x44 - 0.064776141276471*x46 + 0.0210302630614787*x48 + x56 + 0.227266065802008).^2 + (-0.0210302630614787*x44 - 0.0210302630614787*x46 + 0.064776141276471*x48 + x51 + 0.064776141276471*x9 + 0.335824723382428).^2 + (0.0210302630614787*qsubs0.*qsubs3 - x145 - 0.064776141276471*x18 - 0.064776141276471*x20 - 0.0210302630614787*x53 - 0.630431596687917).^2); 0.315*sqrt((x13.*x147 + x13.*x151 - x13.*x154 + x146.*x2 + x148.*x21 + x149.*x2 + x150.*x21 - x152.*x2 - x153.*x21 + x42 - 0.456190476190476).^2 + (x146.*x47 + x147.*x43 + x148.*x50 + x149.*x47 + x150.*x50 + x151.*x43 - x152.*x47 - x153.*x50 - x154.*x43 + x51 + 0.162222222222222).^2 + (x146.*x55 + x147.*x54 + x148.*x52 + x149.*x55 + x150.*x52 + x151.*x54 - x152.*x55 - x153.*x52 - x154.*x54 + x56 + 0.447936507936508).^2); 0.315*sqrt((-x108 - x155.*x55 - x156.*x54 - x157.*x52 - x158.*x55 - x159.*x52 - x160.*x54 + 0.350603174603175*x26.*x55 + 0.350603174603175*x29.*x52 + 0.350603174603175*x33.*x54 + 0.376507936507937).^2 + (-x13.*x156 - x13.*x160 + 0.350603174603175*x13.*x33 - x145 - x155.*x2 - x157.*x21 - x158.*x2 - x159.*x21 + 0.350603174603175*x2.*x26 + 0.350603174603175*x21.*x29 - 0.486984126984127).^2 + (qsubs0.*qsubs2 - x155.*x47 - x156.*x43 - x157.*x50 - x158.*x47 - x159.*x50 - x160.*x43 - x20 + 0.350603174603175*x26.*x47 + 0.350603174603175*x29.*x50 + 0.350603174603175*x33.*x43 + 0.173333333333333).^2); 0.315*sqrt((x170 - 0.416507936507936).^2 + (x171 - 0.0765079365079365).^2 + (x172 - 0.213968253968254).^2); 0.315*sqrt((x170 - 0.393650793650794).^2 + (x171 - 0.104126984126984).^2 + (x172 - 0.317777777777778).^2); 0.315*sqrt((x13.*x175 + x13.*x178 + x13.*x180 + x173.*x2 + x174.*x21 + x176.*x2 + x177.*x21 + x179.*x2 + x181.*x21 + x42 - 0.273015873015873).^2 + (x173.*x47 + x174.*x50 + x175.*x43 + x176.*x47 + x177.*x50 + x178.*x43 + x179.*x47 + x180.*x43 + x181.*x50 + x51 - 0.0838095238095238).^2 + (x173.*x55 + x174.*x52 + x175.*x54 + x176.*x55 + x177.*x52 + x178.*x54 + x179.*x55 + x180.*x54 + x181.*x52 + x56 - 0.322222222222222).^2); 0.315*sqrt((x13.*x184 + x13.*x187 + x13.*x189 + x182.*x2 + x183.*x21 + x185.*x2 + x186.*x21 + x188.*x2 + x190.*x21 + x42 - 0.298095238095238).^2 + (x182.*x47 + x183.*x50 + x184.*x43 + x185.*x47 + x186.*x50 + x187.*x43 + x188.*x47 + x189.*x43 + x190.*x50 + x51 - 0.0965079365079365).^2 + (x182.*x55 + x183.*x52 + x184.*x54 + x185.*x55 + x186.*x52 + x187.*x54 + x188.*x55 + x189.*x54 + x190.*x52 + x56 - 0.382857142857143).^2); 0.315*sqrt((x13.*x192 + x13.*x196 - x13.*x199 + x191.*x2 + x193.*x21 + x194.*x2 + x195.*x21 - x197.*x2 - x198.*x21 + x42 - 0.608253968253968).^2 + (x191.*x47 + x192.*x43 + x193.*x50 + x194.*x47 + x195.*x50 + x196.*x43 - x197.*x47 - x198.*x50 - x199.*x43 + x51 + 0.323492063492063).^2 + (x191.*x55 + x192.*x54 + x193.*x52 + x194.*x55 + x195.*x52 + x196.*x54 - x197.*x55 - x198.*x52 - x199.*x54 + x56 + 0.177460317460317).^2); 0.315*sqrt((x13.*x201 + x13.*x205 - x13.*x208 + x2.*x200 + x2.*x203 - x2.*x206 + x202.*x21 + x204.*x21 - x207.*x21 + x42 - 0.603492063492063).^2 + (x200.*x47 + x201.*x43 + x202.*x50 + x203.*x47 + x204.*x50 + x205.*x43 - x206.*x47 - x207.*x50 - x208.*x43 + x51 + 0.346031746031746).^2 + (x200.*x55 + x201.*x54 + x202.*x52 + x203.*x55 + x204.*x52 + x205.*x54 - x206.*x55 - x207.*x52 - x208.*x54 + x56 + 0.139047619047619).^2); 0.315*sqrt((x13.*x210 + x13.*x214 - x13.*x217 + x2.*x209 + x2.*x212 - x2.*x215 + x21.*x211 + x21.*x213 - x21.*x216 + x42 - 0.562222222222222).^2 + (x209.*x47 + x210.*x43 + x211.*x50 + x212.*x47 + x213.*x50 + x214.*x43 - x215.*x47 - x216.*x50 - x217.*x43 + x51 + 0.378095238095238).^2 + (x209.*x55 + x210.*x54 + x211.*x52 + x212.*x55 + x213.*x52 + x214.*x54 - x215.*x55 - x216.*x52 - x217.*x54 + x56 + 0.079047619047619).^2); 0.315*sqrt((x218.*x55 + x219.*x54 + x220.*x52 + x221.*x55 + x222.*x52 + x223.*x54 - x224.*x55 - x225.*x52 - x226.*x54 + x56).^2 + (x13.*x219 + x13.*x223 - x13.*x226 + x2.*x218 + x2.*x221 - x2.*x224 + x21.*x220 + x21.*x222 - x21.*x225 + x42 - 0.588253968253968).^2 + (x218.*x47 + x219.*x43 + x220.*x50 + x221.*x47 + x222.*x50 + x223.*x43 - x224.*x47 - x225.*x50 - x226.*x43 + x51 + 0.408571428571429).^2); 0.315*sqrt((-x108 + 0.404*x16.*x54 - x233.*x55 - x234.*x52 - x235.*x54 + 0.404*x24.*x52 + 0.404571428571429*x26.*x55 + 0.404571428571429*x29.*x52 + 0.404571428571429*x33.*x54 + 0.404*x55.*x7 - 0.0850793650793651).^2 + (x13.*x229 + x13.*x231 - x13.*x235 + x2.*x227 + x2.*x230 - x2.*x233 + x21.*x228 + x21.*x232 - x21.*x234 + x42 - 0.577777777777778).^2 + (x227.*x47 + x228.*x50 + x229.*x43 + x230.*x47 + x231.*x43 + x232.*x50 - x233.*x47 - x234.*x50 - x235.*x43 + x51 + 0.440634920634921).^2); -0.399822570499435*qsubs1 + 2.08439558956446*qsubs3.*qsubs6.*qsubs7 - 0.337805912285226*qsubs5 - 0.180790241991161*qsubs7 - 0.487275968966565*x236 + 0.305616956876084*x5 + 0.130834291696229; -0.376462955285001*qsubs1 + 8.1167257326656*qsubs3.*qsubs6.*x5 - 0.329301245560537*qsubs5 - 0.55536750262959*x236 - 0.36946267645132*x237 + 1.22603768298756*x31.*x5 + 0.133536072817355; 34.0104684167836*qsubs1.*qsubs3.*x5 - 0.35268375963768*qsubs1 - 0.306290590298819*qsubs5 - 0.555385881056343*x236 - 0.349921929631479*x237 + 1.38649765126187*x31.*x5 + 0.137803415922912; 1.54971157330919*qsubs3.*x239 - 0.182077106822319*qsubs3 - 0.267519823420963*qsubs5 + 6.7828360514234*qsubs6.*qsubs7.*x11 - 0.536869621497297*x236 - 0.282938814081739*x237 - 0.608788078758737*x238 + 1.3986153803831*x31.*x5 + 0.136434505195217; 15.5955173621884*qsubs1.*qsubs3.*x5 + 1.62513766430436*qsubs3.*x239 - 0.23691067334513*qsubs3 + 1.3065827614352*qsubs5.*qsubs7.*x31 - 0.248972534646396*qsubs5 - 0.531621179777355*x236 - 0.465699759546125*x238 - 0.509174489059389*x240 + 0.127861053583364; 0.434766872368921*qsubs2.*qsubs3 - 0.142848998917806*qsubs5 - 0.48574913711943*x236 - 0.217916617815*x238 - 0.0795404564706214*x241 + 1.23894319516036*x242 + 0.123083864039469; -0.170433242647989*qsubs1.*x239 + 0.529619772207045*qsubs5.*x236 + 0.0472885934126085*qsubs7 - 0.4080223732204*x236 + 0.401474409773715*x46 + 0.0936457701983474; 0.453424983088349*qsubs2.*x245 - 1.26599876924795*qsubs3.*x246 + 0.0235210691377822*qsubs7 - 0.405939409530951*x236 - 0.433624330531867*x243 + 0.949085946132232*x244 + 0.485811355637738*x247 + 0.0697933662154756; -0.101189769155184*qsubs5.*x239 - 0.432388998829561*x236 + 3.50751758920994*x237.*x46 + 0.0169457669529891*x237 + 0.572372883148437*x247 + 6.41521199130295*x248 + 3.61164483793332*x250 + 6.15772047101551*x251 + 0.610571438617797*x253 + 0.0513885655348751; -2.98223336487227*qsubs1.*qsubs5.^3 - 0.950360833033398*x236 + 0.0568236450058926*x241 + 1.82179130932886*x244 + 0.671501705162839*x247 - 3.59384071716893*x250 + 10.1447865351693*x251 - 5.71617143856959*x254 + 0.466890522137556*x256 - 2.20479715165018*x31.*x9; 10.6447869171986*qsubs2.*x257 - 1.113145727487*qsubs2.*x5 - 0.998693268342873*qsubs6.*x249 + 0.0631927147806303*qsubs6 - 0.39278279025827*x236 - 0.370028206891712*x240 + 7.19301352356576*x248 + 4.37396748783788*x250 + 6.83986006706788*x251 + 0.721081562957841*x258 + 1.01021086009398*x3.*x31 + 0.225851495023584*x31.*x53 + 5.65253129235253*x46.*x5; 3.08333218443898*qsubs1.*qsubs3.*qsubs7 - 0.0143424537179752*qsubs1 + 0.12458674502311*qsubs2.*qsubs7 - 67.5269867445807*qsubs7.*x257 - 0.479737075098136*x236 - 3.06237386815939*x242 + 1.27124218843968*x249; 0.103119709721391*qsubs11 + 0.0163966296927972*qsubs9 + 0.0996925485901122*x259 - 0.351576540092049*x260 + 0.0275621766656404*x261 + 0.491853027011381*x262 + 0.173854334484687; 0.0575525196092037*qsubs11 + 0.0439620757076519*qsubs9 + 0.070475874894497*x263 - 0.0678402680064331*x264 + 0.152153544521287; 0.0440696501062466*qsubs9 + 0.0459234076908767*x259 + 0.124688385740236*x263 - 0.118618649653876*x265 - 0.177373036583537*x266 + 0.138707479919222; 0.0428588603097824*qsubs9 + 0.06671956781407*x259 + 0.0289089418484292*x261 + 0.0528831067517571*x263 - 0.060910813091805*x267 + 0.113540643294245; 0.0102911001944179*qsubs10 - 0.0201666410687433*qsubs11 + 0.0273740637472172*qsubs9 + 0.108277490626677*x259 - 0.0418579463780774*x264 + 0.107209134329072; 0.0111934918480468*qsubs10 + 0.0397471147998811*qsubs9 + 0.115889396453201*x259 - 0.0598083508099709*x264 + 0.143771673543336; 0.00912709852187035*qsubs10 + 0.038662662157145*qsubs9 + 0.123091570890802*x259 - 0.0614671657568094*x264 + 0.15904322614556; -0.0205700343247733*qsubs11 + 0.0240966956897555*qsubs9 + 0.148831619854703*x259 - 0.0530464821881656*x264 + 0.142282466915136; -0.0307343260929937*qsubs11 + 0.0143989254556171*qsubs9 + 0.15439526793052*x259 - 0.0416375068590483*x264 + 0.151213534137524; -0.138837172491663*qsubs11.*x259 - 0.0577540564750188*qsubs11 + 0.156985078427931*x259 + 0.0140852420262484*x261 + 0.141591032987101; 0.00484374139723913*qsubs10 - 0.0615111679418772*qsubs11 - 0.00611388193374257*qsubs9 - 0.0480170404843219*x262 + 0.131531202681082; -0.0226471267515631*qsubs5 - 0.136978804959263*qsubs7.*x264 + 0.0448736892316174*qsubs7 - 0.041612521298193*x264 - 0.450144112704636*x268 - 0.0648119005745192*x269 - 0.541476707759164*x270 + 0.147079173216751; 0.0967087296994761*qsubs10.*qsubs5.*qsubs6.*qsubs7 + 0.0593401696033345*qsubs11.*qsubs5.*x31 + 0.0150083394873738*qsubs11 - 0.456891358052375*qsubs6.*x272 + 0.0822920293068868*qsubs7 - 0.0565117274795606*x241 - 0.0929208618123594*x255 - 1.00562568014021*x271 - 0.981067387617282*x273 - 0.0927585537396005*x274 - 1.01249129658729*x275 + 0.178980641559182; 0.65348301747196*qsubs10.*qsubs5.*qsubs6.*qsubs7 + 0.943235664219779*qsubs5.*qsubs6.*x5 + 0.0774884389962967*qsubs7 - 0.0464825492604718*x241 - 0.104142482392191*x255 - 1.10625891879128*x271 - 0.827999768206373*x273 - 0.0971405053121815*x274 - 0.921858382866771*x277 + 0.156738082127525; 0.226317877930739*qsubs5.*x267 - 0.0423497574927131*qsubs5 + 0.090817308965499*qsubs7 - 0.866007485143719*x237.*x264 - 0.111359866267483*x268 - 0.150575962479161*x269 - 0.723007440620532*x270 + 0.251857218240734*x274 + 0.028487588639526*x276 + 0.170038762369198; 0.0188620756257299*qsubs11 - 0.063775914803762*qsubs9 - 0.0614301127259638*x263 - 0.026325961468536*x264 - 0.00422608702503588*x265 - 0.0974470873346288*x278 + 0.164797864939309; 0.0161868664922966*qsubs11 - 0.0617021986293006*qsubs9 - 1.11466708191827*x278 + 0.179463819129385; 0.224692885346393*qsubs10.*qsubs11 + 0.0659676379858262*qsubs10.*x262 + 1.57974089926596*qsubs9.*x261 - 0.0433655710727013*qsubs9 - 0.493641248484641*x259 - 0.261742210728811*x266 - 1.86262266502398*x279 - 0.888814674567964*x280 + 0.19399284158262; 0.0351600743224401*qsubs10 - 0.00548139887278397*qsubs9 - 0.0389969988161194*x261 - 0.0362120393795274*x264 + 0.0888796544503366; 0.00654382992935149*qsubs10.*qsubs9 + 0.0266356741616388*qsubs10 - 0.0131134760844284*qsubs11 - 0.00642196256592078*qsubs9 - 0.0391253793258046*x261 + 0.655096312832809*x262.*x264 - 0.259083222930599*x281 + 0.17373755616591; 0.028306893384337*qsubs10 - 0.015516448338057*qsubs11 - 0.00602156288292354*qsubs9 - 0.0391750470255718*x261 + 0.612311276474412*x262.*x264 - 0.246977598444628*x281 + 0.148811534845632; 0.0126226743839074*qsubs10 - 0.0380505714290533*qsubs11 - 0.00697715331593287*qsubs9 - 0.011714397057749*x260 + 0.0405702592492048*x263 + 0.142000685076213; 0.0112860882462805*qsubs10 - 0.0460808610822636*qsubs11 - 0.00461135114705151*qsubs9 + 0.0363385462925972*x263 + 0.0552032482961466*x272 + 0.0658762907510997*x282 + 0.128818690358278; 0.0185457526054213*qsubs10 - 0.0412718832714044*qsubs11 - 0.039660516307003*x259 + 0.0426699048603827*x263 + 0.120542673499911*x272 + 0.0181166417275067*x282 + 0.11167049423169; 0.0205380569914046*qsubs11 - 0.00264898816667054*qsubs9 + 0.96939410280136*x261.*x262 + 0.045930727922735*x261 + 0.0800263257673297*x263 - 0.0329982592913917*x267 + 0.0966260524804794; 0.186724644053863*qsubs10.*qsubs11 - 0.0146670125237523*x259 - 0.0173618745526012*x260 + 0.0418167581157975*x261 - 0.255591674215367*x280 - 0.578553529977569*x283 + 0.114411602909201; 1.91533221392964*qsubs10.*x284 + 0.0189316664652198*x259 + 0.045066101832991*x261 + 0.11891945259326*x263 - 0.0287097547473467*x283 - 0.277749763634791*x286 - 0.953878015385076*x287 + 0.0744684742354628; -0.0181890364076043*qsubs10 + 0.121186018887623*qsubs11 + 0.0425958280083888*qsubs9 + 0.095540283710137; -0.0170813438147603*qsubs10 + 0.151405139970193*qsubs11 + 0.0501686763910624*qsubs9 - 1.61158333748996*x286 + 0.112755578657816; -0.0168374479938973*qsubs10 + 0.129753793505414*qsubs11 + 0.0517624834535466*qsubs9 + 0.26005364650762*x278 + 0.107068558908372; -0.0196491066462149*qsubs10 + 0.091402507763326*qsubs11 + 0.0391313934303184*qsubs9 + 0.0904996665385804; 0.0167270265434543*qsubs10 - 0.0487200746084361*qsubs11 - 0.00533896889005948*qsubs9 + 0.00960816515992398*x262 - 0.0775024825257679*x285 + 0.0937577073077334; 0.0204394887054436*qsubs10 - 0.0483550558564048*qsubs11 - 0.00662992524573726*qsubs9 - 0.0542783189538078*x285 + 0.0974397225719821; 0.0543242069608756*qsubs10.*qsubs9 - 0.0228935928489545*qsubs10 - 0.0450001980765749*qsubs11 - 0.0214714110798854*qsubs9 + 0.118478898641868; 0.0365630636327083*qsubs10.*qsubs11 - 0.0278768327461594*qsubs10 - 0.0481554188595293*qsubs11 - 0.0133322033284317*qsubs9 + 0.0959675668749218; -0.0314832094589498*qsubs10 - 0.0168011736800951*qsubs11 + 0.0790457329800364*qsubs9.*x267 - 0.0319562253346934*x259 - 0.0762053892995415*x262 + 0.0352885551638182*x264 - 0.0717241898890914*x280 + 0.0742754476315841; -0.0144552707351456*qsubs10 + 0.133066081055831*qsubs9.*x261 - 0.0109083846308911*qsubs9 + 0.0699125771276063*x267 - 0.0931137882906683*x278 - 0.202743695710675*x279 - 0.150956938136127*x280 + 0.0974952560639379; -0.0187272117192585*qsubs10 - 0.0145886638531363*qsubs11 + 0.157779181945604*qsubs9.*x264 - 0.119179687649886*x260 - 0.0442721858085189*x261 - 0.0451702898597567*x263 + 0.114916159074208; -0.0439805161247504*qsubs10 - 0.00799137112146561*qsubs11 + 0.145774199795758*qsubs9.*x261 - 0.0972736426990291*x281 + 0.117606655555058; 0.0572806780620305*qsubs10.*qsubs9 - 0.0466644538309976*qsubs10 - 0.01281646306127*qsubs11 - 0.20135279054815*qsubs9.*x267 + 0.118995859083807; 0.103521455285466*qsubs10.*qsubs9 - 0.0524907627949378*qsubs10 - 0.0134128112630631*x263 - 0.0141913938955739*x264 + 0.143604437635507; -0.0472979505419058*qsubs10 - 0.453655755495786*qsubs11.*x284 - 0.00411974585487677*qsubs11 + 0.0187124692638796*qsubs9 + 0.0282957949412761*x259 + 0.138028129553961; 0.122340021125013*qsubs10.*qsubs9 - 0.0637430794194577*qsubs10 - 0.165424685054499*qsubs11.*x285 - 0.00284919143889254*qsubs11 + 0.0479446395455729*x261 - 0.0124439977856422*x265 + 0.0975732504406606; -0.0576395114293968*qsubs10 + 0.130113732408752*x259 + 0.0590247515039939*x260 + 0.0504668283735453*x261 - 0.0158642604245204*x263 - 0.193140261919918*x282 + 0.115591277378646; -0.029709879645705*qsubs10 + 0.0313659675100647*qsubs11 + 0.0210456019180538*qsubs9 + 0.0577377213139668; -0.0207079135285908*qsubs10 + 0.0437024029115847*qsubs11 + 0.0193533683310695*qsubs9 + 0.0383220730591165; -0.0201075385626289*qsubs11 + 0.134682886347591*qsubs12.*x267 - 0.165712714958471*x288 - 0.16496758403447*x290 - 0.228345187579548*x291 + 0.370306078988325; 0.0369816315835864*qsubs11.*qsubs12 - 0.0496884799460417*qsubs12.*x264 + 0.481737132473477*x261.*x289 - 0.307615678177299*x263.*x289 - 0.118479163191899*x282 - 0.24747346715704*x290 - 0.217790697409524*x292 + 0.34550552211701; 0.0764799758609488*qsubs11 + 0.400042426826882*qsubs12.*x285 - 0.591437926485703*x261.*x289 + 0.119237153409798*x288 + 0.0909396860450745*x292 + 0.422547631434471*x293 + 0.296194335879241; x295; 0.110796529808048*qsubs10.^4 + 0.0576583828079135*qsubs11 - 0.394580963511785*qsubs12.*x281 - 0.0329326418666452*x261 + 0.300540165014576*x288 + 0.129371364427049*x291 + 0.267485274425019*x293 + 0.299224219479144; x295; 0.10624889999299*qsubs1 + 0.0734035750239836*qsubs11 - 4.33429085931234*qsubs3.*qsubs6.*x259 + 0.17278317815387*qsubs3 + 0.109856256667747*qsubs5 + 0.14122619239113*qsubs7 - 1.48339107186098*qsubs9.*x298 - 0.769656246975948*x275 + 1.08311483113457*x287 - 0.0430468117130353*x296 - 0.0481745681811966*x297 - 0.468654843303337*x299 + 2.39682201788826*x300 + 0.206572844729267; 0.100761526529349*qsubs1 + 0.0824062870345599*qsubs11 + 0.244457497147725*qsubs3 + 0.126304340735247*qsubs5 + 0.0332290755392842*qsubs6.*qsubs9 + 0.120901852327823*qsubs7 - 3.02057880876778*x260.*x298 - 0.713628920145938*x275 + 0.0353624219782261*x297 - 0.441144182329852*x299 + 2.25272450745676*x300 + 0.234044347423462; 0.0956917776008449*qsubs1 + 0.0196629850648468*qsubs11 + 7.83523668935775*qsubs2.*qsubs9.*x20 + 0.292837892593824*qsubs3 + 4.32273589962962*qsubs5.*x272 + 0.127361069169863*qsubs5 + 0.0250804846965838*qsubs9 + 14.4935586675078*x242 + 12.4379707228384*x261.*x298 + 1.24358256083728*x297 + 1.97565779529429*x300 + 0.117929565753615*x301 - 5.15350455521096*x302 - 17.4240138195612*x303 + 0.273257304537178; -1.96763383614109*qsubs11.*x244 + 0.0316878165098613*qsubs11 + 0.324815100445844*qsubs3 + 0.118777220935895*qsubs5 + 0.0396789775933978*qsubs7 + 0.0233659961952198*qsubs9 - 2.45234713516217*x20.*x31 - 0.921028458409556*x275 + 0.297079362230543*x297 + 19.6589368664374*x3.*x305 + 1.57650747054983*x300 + 0.149245890879109*x301 + 0.295852161004074*x304 - 4.10038412943159*x306 + 1.37756039840569*x307 + 0.315686790005635; -0.0378144084612305*qsubs10 + 0.328513339450045*qsubs3 + 2.33682381277859*qsubs5.*x5 + 0.0621563274668142*qsubs5 - 4.51973779780459*qsubs7.*x272 + 0.0279495053813348*qsubs9 - 0.384821324425753*x246 + 0.812710449273725*x297 - 1.86955611584901*x302 + 0.208951990590126*x304 - 11.2204230285749*x306 + 2.30987478326477*x307 + 0.691113087493532*x308 - 1.46194708426307*x48.*x9 + 0.277317548908613; -0.310374070750152*qsubs1.*x260 - 0.0325703873366343*qsubs10 + 0.51700686572808*qsubs3 + 0.0644159315905401*qsubs5 - 0.315123798038869*x246 + 5.03844456663656*x256 + 0.0794661006289209*x260 + 4.4819345645131*x261.*x5 + 0.0857249008888518*x270 - 3.93261796806244*x277 + 0.7306342133818*x297 + 0.855174969714394*x298.*x31 - 12.3516206913975*x303 + 0.218830586438833*x304 - 0.0865181290920689*x305 - 14.3731503415791*x306 + 0.640298607511267*x46 + 0.250477771520652; 1.56855155499677*qsubs10.*x53 + 7.70580594891782*qsubs3.*x312 + 0.400114376915112*qsubs3 - 3.1418197671818*qsubs5.*x269 - 2.96896710609954*qsubs7.*qsubs9.*x9 - 13.0853450738482*x254 + 9.74144368874185*x263.*x48 + 0.30312161325448*x266 - 0.708602618868285*x269 - 4.24447223887906*x271 + 6.62072785057717*x286 + 11.0926016080849*x296.*x9 + 0.34207429952756*x298 + 7.20249095479694*x308 - 5.30186458951792*x310 - 1.02328755877549*x311 + 0.175141745778279; 8.75702512458436*qsubs1.*x284 + 3.09181556232702*qsubs11.*qsubs2.*x5 - 0.302515071877823*qsubs11.*qsubs2 + 2.77124123907308*qsubs2.*qsubs3.*qsubs6 + 0.600921844076446*qsubs3 - 7.67888797178209*qsubs5.*x299 + 0.194802757777472*qsubs7 - 4.06958158408515*x238.*x261 - 1.80419045344497*x261.*x276 - 1.26143131484465*x270 - 0.136817586936268*x271 - 0.585754402202252*x310 - 0.335123920635408*x313 + 0.179775656236991; -3.73044988969287*qsubs1.*x299 + 0.0836828796069849*qsubs11 - 0.625645881538184*qsubs2.*x260 + 0.35186246956345*qsubs3 + 3.5334489861314*qsubs5.*x252 - 0.884990315751177*x236.*x261 - 1.41705173716095*x258 - 1.57553617577551*x270 - 0.0468470972117773*x271 + 0.319181890488637*x298 - 0.834506801691598*x309 + 2.175847626812*x31.*x46 + 0.733813453472997*x310 + 1.26509823037655*x311 + 0.147080348497635*x48 + 0.194860970049359; 2.28465330510282*qsubs1.*x294 + 0.0777206069270714*qsubs11 + 0.23749959400274*qsubs3 - 0.103143782076465*x241 + 2.00052757800297*x253 - 0.989947409079768*x270 - 1.05704481641875*x271 - 0.822549977784942*x309 + 0.823767203096224*x310 + 0.578652791110293*x311 - 1.43356915421536*x314 + 3.55060977578907*x315 - 3.36008061457673*x316 + 0.188918738475118; 0.0563653367227722*qsubs11 + 0.0979960005652133*qsubs3 + 0.392385471088946*qsubs7.*x48 - 1.33362945192003*x268 - 0.993963618355393*x270 + 1.90717334085092*x308 - 0.671395069835252*x309 + 0.842348485834007*x311 - 0.750638035862421*x314 + 1.25896058581689*x315 - 2.16197517881284*x316 + 0.187818865639014; 3.96706905444155*qsubs1.*qsubs2.*x261 + 3.09558295905486*qsubs1.*x285 - 1.53087236176112*qsubs11.*x249 + 0.050163686046596*qsubs11 - 2.08536934661153*x237.*x9 - 0.271455508416888*x245 - 0.485085555955162*x264.*x9 - 2.18537155772031*x268 - 1.09572330141245*x270 - 1.18072302098663*x309.*x48 - 0.805727112123131*x309 - 0.45648540787015*x313 - 0.226978802545747*x46 + 0.191104486356178; 0.0434864586668707*qsubs11 - 0.170822370651655*qsubs6.*x263 - 0.10310945373823*qsubs6.*x264 + 0.125456699273232*qsubs7 - 0.179679739844359*x255 - 1.26185548574856*x268 - 1.2905974438352*x270 - 0.67570025658847*x275 - 0.272631535932667*x313 + 0.167933071538227; 4.83332786433184*qsubs10.*qsubs6.*x3 - 0.843719472113551*qsubs11.*x245 + 0.0844193382922012*qsubs7 - 2.62269906626424*x268 - 0.962371110651338*x270 - 0.378090052007643*x313 + 0.132728785338186; 0.0205525017445139*qsubs12 + 0.146154088427195; 0.0210747598991912*qsubs12 + 0.169710060635468; 0.0183286995272377*qsubs12 + 0.0558223257333852; 0.0210365931050275*qsubs12 + 0.0934988379050989; 0.021189549801054*qsubs12 + 0.136718812703705; 0.188349828999743 - 0.00935636496814369*qsubs12; 0.180735936766646 - 0.0115399436250548*qsubs12; 0.146695880367282 - 0.0125114621673938*qsubs12; 0.142856705518397 - 0.0101778654826064*qsubs12; 0.110531882969394 - 0.0112020979254651*qsubs12; 0.0981481135002244 - 0.0111856660891357*qsubs12; 0.0952952654336992 - 0.0097291513074297*qsubs12; 0.325325670237825 - 0.0278750450497314*x289; 0.312962084506815 - 0.0168445062433603*qsubs12; 0.298620425443565 - 0.0171791602867297*qsubs12; 0.156566379110345 - 0.00591328095539027*x289; 0.0841168136882597; 0.0214219226407107; 0.0192336228273652; 0.0185570263811613; 0.0180353484742324; 0.0186844961702131; 0.0363101932226377; 0.0291849168084541; 0.0279579165288706; 0.0219876579145918*qsubs12 + 0.221970762976743; 0.0226945032832469*qsubs12 + 0.224365826699157; 0.0197959668123318*qsubs12 + 0.229873924120311; 0.0196640729225794*qsubs12 + 0.190153580274778; 0.0193141915105747*qsubs12 + 0.185696914447041; 0.00794039835307144*qsubs12 + 0.0156191932975916; 0.00962939753164361*x289 + 0.0164687065275715; 0.00948511022761933*qsubs12 + 0.0175108488372233; 0.0478344240001351; 0.00525897556532031*x289 + 0.0315427261915965];