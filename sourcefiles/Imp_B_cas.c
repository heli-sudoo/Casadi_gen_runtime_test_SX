/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) Imp_B_cas_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_sq CASADI_PREFIX(sq)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

casadi_real casadi_sq(casadi_real x) { return x*x;}

static const casadi_int casadi_s0[18] = {14, 1, 0, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
static const casadi_int casadi_s1[8] = {4, 1, 0, 4, 0, 1, 2, 3};
static const casadi_int casadi_s2[6] = {2, 1, 0, 2, 0, 1};

/* Imp_B:(i0[14],i1[4])->(o0[14],o1[2]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a100, a101, a102, a103, a104, a105, a106, a107, a108, a109, a11, a110, a111, a112, a113, a114, a115, a116, a117, a118, a119, a12, a120, a121, a122, a123, a124, a125, a126, a127, a128, a129, a13, a130, a131, a132, a133, a134, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a4, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a5, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a6, a60, a61, a62, a63, a64, a65, a66, a67, a68, a69, a7, a70, a71, a72, a73, a74, a75, a76, a77, a78, a79, a8, a80, a81, a82, a83, a84, a85, a86, a87, a88, a89, a9, a90, a91, a92, a93, a94, a95, a96, a97, a98, a99;
  a0=arg[0]? arg[0][0] : 0;
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[0]? arg[0][1] : 0;
  if (res[0]!=0) res[0][1]=a0;
  a0=arg[0]? arg[0][2] : 0;
  if (res[0]!=0) res[0][2]=a0;
  a1=arg[0]? arg[0][3] : 0;
  if (res[0]!=0) res[0][3]=a1;
  a2=arg[0]? arg[0][4] : 0;
  if (res[0]!=0) res[0][4]=a2;
  a3=arg[0]? arg[0][5] : 0;
  if (res[0]!=0) res[0][5]=a3;
  a4=arg[0]? arg[0][6] : 0;
  if (res[0]!=0) res[0][6]=a4;
  a5=cos(a0);
  a6=5.4600000000000000e+000;
  a7=cos(a3);
  a8=1.2680000000000000e+000;
  a9=1.2800000000000000e-001;
  a10=cos(a4);
  a11=(a9*a10);
  a12=(a11*a10);
  a13=sin(a4);
  a14=(a9*a13);
  a15=(a14*a13);
  a12=(a12+a15);
  a12=(a8+a12);
  a15=(a7*a12);
  a16=sin(a3);
  a17=(a9*a10);
  a18=(a17*a13);
  a19=(a9*a13);
  a20=(a19*a10);
  a18=(a18-a20);
  a20=(a16*a18);
  a15=(a15+a20);
  a20=(a15*a7);
  a21=(a14*a10);
  a22=(a11*a13);
  a21=(a21-a22);
  a22=(a7*a21);
  a23=(a19*a13);
  a24=(a17*a10);
  a23=(a23+a24);
  a23=(a8+a23);
  a24=(a16*a23);
  a22=(a22+a24);
  a24=(a22*a16);
  a20=(a20+a24);
  a20=(a6+a20);
  a24=cos(a1);
  a25=cos(a2);
  a26=(a9*a25);
  a27=(a26*a25);
  a2=sin(a2);
  a28=(a9*a2);
  a29=(a28*a2);
  a27=(a27+a29);
  a27=(a8+a27);
  a29=(a24*a27);
  a1=sin(a1);
  a30=(a9*a25);
  a31=(a30*a2);
  a32=(a9*a2);
  a33=(a32*a25);
  a31=(a31-a33);
  a33=(a1*a31);
  a29=(a29+a33);
  a33=(a29*a24);
  a34=(a28*a25);
  a35=(a26*a2);
  a34=(a34-a35);
  a35=(a24*a34);
  a36=(a32*a2);
  a37=(a30*a25);
  a36=(a36+a37);
  a8=(a8+a36);
  a36=(a1*a8);
  a35=(a35+a36);
  a36=(a35*a1);
  a33=(a33+a36);
  a20=(a20+a33);
  a33=(a5*a20);
  a36=sin(a0);
  a37=(a7*a18);
  a38=(a16*a12);
  a37=(a37-a38);
  a38=(a37*a7);
  a39=(a7*a23);
  a40=(a16*a21);
  a39=(a39-a40);
  a40=(a39*a16);
  a38=(a38+a40);
  a40=(a24*a31);
  a41=(a1*a27);
  a40=(a40-a41);
  a41=(a40*a24);
  a42=(a24*a8);
  a43=(a1*a34);
  a42=(a42-a43);
  a43=(a42*a1);
  a41=(a41+a43);
  a38=(a38+a41);
  a41=(a36*a38);
  a33=(a33+a41);
  a41=(a33*a5);
  a43=(a22*a7);
  a44=(a15*a16);
  a43=(a43-a44);
  a44=(a35*a24);
  a45=(a29*a1);
  a44=(a44-a45);
  a43=(a43+a44);
  a44=(a5*a43);
  a45=(a39*a7);
  a46=(a37*a16);
  a45=(a45-a46);
  a6=(a6+a45);
  a45=(a42*a24);
  a46=(a40*a1);
  a45=(a45-a46);
  a6=(a6+a45);
  a45=(a36*a6);
  a44=(a44+a45);
  a45=(a44*a36);
  a41=(a41+a45);
  a45=1.;
  a46=-1.9500000000000001e-001;
  a47=cos(a4);
  a48=cos(a3);
  a49=cos(a0);
  a50=(a48*a49);
  a51=sin(a3);
  a52=sin(a0);
  a53=(a51*a52);
  a50=(a50-a53);
  a53=(a47*a50);
  a54=sin(a4);
  a55=cos(a3);
  a56=(a55*a52);
  a57=sin(a3);
  a49=(a57*a49);
  a56=(a56+a49);
  a56=(a54*a56);
  a53=(a53-a56);
  a53=(a46*a53);
  a56=-2.0899999999999999e-001;
  a50=(a56*a50);
  a49=-1.9000000000000000e-001;
  a52=(a49*a52);
  a50=(a50-a52);
  a53=(a53+a50);
  a50=casadi_sq(a53);
  a50=(a45+a50);
  a52=cos(a0);
  a58=cos(a3);
  a59=(a52*a58);
  a60=sin(a0);
  a61=sin(a3);
  a62=(a60*a61);
  a59=(a59-a62);
  a62=(a47*a59);
  a63=sin(a3);
  a64=(a52*a63);
  a3=cos(a3);
  a65=(a60*a3);
  a64=(a64+a65);
  a64=(a54*a64);
  a62=(a62-a64);
  a62=(a46*a62);
  a59=(a56*a59);
  a62=(a62+a59);
  a59=casadi_sq(a62);
  a50=(a50+a59);
  a59=(a52*a55);
  a64=(a60*a57);
  a59=(a59-a64);
  a64=cos(a4);
  a59=(a59*a64);
  a52=(a52*a51);
  a60=(a60*a48);
  a52=(a52+a60);
  a4=sin(a4);
  a52=(a52*a4);
  a59=(a59-a52);
  a59=(a46*a59);
  a52=casadi_sq(a59);
  a50=(a50+a52);
  a50=sqrt(a50);
  a52=(a41/a50);
  a60=-2.5360000000000001e-002;
  a65=-7.8079999999999998e-003;
  a66=(a65*a10);
  a67=(a56*a10);
  a11=(a11*a67);
  a66=(a66+a11);
  a11=(a56*a13);
  a14=(a14*a11);
  a66=(a66+a14);
  a66=(a60+a66);
  a14=(a7*a66);
  a17=(a17*a11);
  a68=(a65*a13);
  a19=(a19*a67);
  a68=(a68+a19);
  a17=(a17-a68);
  a68=(a16*a17);
  a14=(a14+a68);
  a68=1.9000000000000000e-001;
  a19=(a68*a16);
  a15=(a15*a19);
  a14=(a14-a15);
  a68=(a68*a7);
  a22=(a22*a68);
  a14=(a14+a22);
  a22=(a65*a25);
  a15=(a56*a25);
  a26=(a26*a15);
  a22=(a22+a26);
  a26=(a56*a2);
  a28=(a28*a26);
  a22=(a22+a28);
  a22=(a60+a22);
  a28=(a24*a22);
  a30=(a30*a26);
  a69=(a65*a2);
  a32=(a32*a15);
  a69=(a69+a32);
  a30=(a30-a69);
  a69=(a1*a30);
  a28=(a28+a69);
  a69=(a49*a1);
  a29=(a29*a69);
  a28=(a28-a29);
  a29=(a49*a24);
  a35=(a35*a29);
  a28=(a28+a35);
  a14=(a14+a28);
  a28=(a5*a14);
  a35=(a7*a17);
  a32=(a16*a66);
  a35=(a35-a32);
  a37=(a37*a19);
  a35=(a35-a37);
  a39=(a39*a68);
  a35=(a35+a39);
  a39=(a24*a30);
  a37=(a1*a22);
  a39=(a39-a37);
  a40=(a40*a69);
  a39=(a39-a40);
  a42=(a42*a29);
  a39=(a39+a42);
  a35=(a35+a39);
  a39=(a36*a35);
  a28=(a28+a39);
  a39=(a53/a50);
  a42=(a28*a39);
  a52=(a52+a42);
  a42=(a7*a66);
  a40=(a16*a17);
  a42=(a42+a40);
  a40=(a5*a42);
  a37=(a7*a17);
  a32=(a16*a66);
  a37=(a37-a32);
  a32=(a36*a37);
  a40=(a40+a32);
  a32=(a62/a50);
  a70=(a40*a32);
  a52=(a52+a70);
  a70=(a65*a10);
  a71=(a7*a70);
  a72=(a65*a13);
  a73=(a16*a72);
  a71=(a71-a73);
  a73=(a5*a71);
  a16=(a16*a70);
  a7=(a7*a72);
  a16=(a16+a7);
  a7=(a36*a16);
  a73=(a73-a7);
  a7=(a59/a50);
  a74=(a73*a7);
  a52=(a52+a74);
  a74=(a52/a50);
  a74=(a41-a74);
  a75=sin(a0);
  a76=(a48*a75);
  a77=cos(a0);
  a78=(a51*a77);
  a76=(a76+a78);
  a78=(a47*a76);
  a79=(a55*a77);
  a75=(a57*a75);
  a79=(a79-a75);
  a79=(a54*a79);
  a78=(a78+a79);
  a78=(a46*a78);
  a76=(a56*a76);
  a49=(a49*a77);
  a76=(a76+a49);
  a78=(a78+a76);
  a76=(a78*a39);
  a49=cos(a0);
  a61=(a49*a61);
  a0=sin(a0);
  a58=(a0*a58);
  a61=(a61+a58);
  a47=(a47*a61);
  a3=(a49*a3);
  a63=(a0*a63);
  a3=(a3-a63);
  a54=(a54*a3);
  a47=(a47+a54);
  a47=(a46*a47);
  a56=(a56*a61);
  a47=(a47+a56);
  a56=(a47*a32);
  a76=(a76+a56);
  a48=(a49*a48);
  a51=(a0*a51);
  a48=(a48-a51);
  a48=(a48*a4);
  a0=(a0*a55);
  a49=(a49*a57);
  a0=(a0+a49);
  a0=(a0*a64);
  a48=(a48+a0);
  a46=(a46*a48);
  a48=(a46*a7);
  a76=(a76+a48);
  a48=(a76/a50);
  a0=casadi_sq(a48);
  a0=(a0+a45);
  a45=(a76*a39);
  a45=(a45-a78);
  a64=casadi_sq(a45);
  a0=(a0+a64);
  a64=(a76*a32);
  a64=(a64-a47);
  a49=casadi_sq(a64);
  a0=(a0+a49);
  a49=(a76*a7);
  a49=(a49-a46);
  a57=casadi_sq(a49);
  a0=(a0+a57);
  a0=sqrt(a0);
  a48=(a48/a0);
  a57=(a74*a48);
  a44=(a44*a5);
  a33=(a33*a36);
  a44=(a44-a33);
  a33=(a44/a0);
  a57=(a57+a33);
  a33=(a52*a39);
  a33=(a28-a33);
  a45=(a45/a0);
  a55=(a33*a45);
  a57=(a57+a55);
  a55=(a52*a32);
  a55=(a40-a55);
  a64=(a64/a0);
  a4=(a55*a64);
  a57=(a57+a4);
  a4=(a52*a7);
  a4=(a73-a4);
  a49=(a49/a0);
  a51=(a4*a49);
  a57=(a57+a51);
  a51=(a57*a48);
  a74=(a74-a51);
  a51=(a28/a50);
  a56=1.1641900000000001e-001;
  a61=4.7131999999999999e-003;
  a54=9.7228799999999997e-004;
  a3=(a65*a67);
  a3=(a54+a3);
  a63=(a9*a67);
  a63=(a65+a63);
  a58=(a63*a67);
  a3=(a3+a58);
  a58=(a9*a11);
  a11=(a58*a11);
  a3=(a3+a11);
  a3=(a61+a3);
  a11=(a19*a66);
  a11=(a3-a11);
  a77=(a68*a17);
  a11=(a11+a77);
  a77=(a63*a10);
  a79=(a58*a13);
  a77=(a77+a79);
  a77=(a60+a77);
  a12=(a19*a12);
  a77=(a77-a12);
  a18=(a68*a18);
  a77=(a77+a18);
  a77=(a77*a19);
  a11=(a11-a77);
  a58=(a58*a10);
  a63=(a63*a13);
  a58=(a58-a63);
  a21=(a19*a21);
  a58=(a58-a21);
  a23=(a68*a23);
  a58=(a58+a23);
  a58=(a58*a68);
  a11=(a11+a58);
  a56=(a56+a11);
  a11=(a65*a15);
  a11=(a54+a11);
  a58=(a9*a15);
  a58=(a65+a58);
  a23=(a58*a15);
  a11=(a11+a23);
  a9=(a9*a26);
  a26=(a9*a26);
  a11=(a11+a26);
  a61=(a61+a11);
  a11=(a69*a22);
  a11=(a61-a11);
  a26=(a29*a30);
  a11=(a11+a26);
  a26=(a58*a25);
  a23=(a9*a2);
  a26=(a26+a23);
  a60=(a60+a26);
  a27=(a69*a27);
  a60=(a60-a27);
  a31=(a29*a31);
  a60=(a60+a31);
  a60=(a60*a69);
  a11=(a11-a60);
  a9=(a9*a25);
  a58=(a58*a2);
  a9=(a9-a58);
  a34=(a69*a34);
  a9=(a9-a34);
  a8=(a29*a8);
  a9=(a9+a8);
  a9=(a9*a29);
  a11=(a11+a9);
  a56=(a56+a11);
  a11=(a56*a39);
  a51=(a51+a11);
  a66=(a19*a66);
  a66=(a3-a66);
  a17=(a68*a17);
  a66=(a66+a17);
  a17=(a66*a32);
  a51=(a51+a17);
  a67=(a65*a67);
  a67=(a54+a67);
  a19=(a19*a70);
  a19=(a67-a19);
  a68=(a68*a72);
  a19=(a19-a68);
  a68=(a19*a7);
  a51=(a51+a68);
  a68=(a51/a50);
  a68=(a28-a68);
  a72=(a68*a48);
  a35=(a5*a35);
  a14=(a36*a14);
  a35=(a35-a14);
  a14=(a35/a0);
  a72=(a72+a14);
  a14=(a51*a39);
  a14=(a56-a14);
  a70=(a14*a45);
  a72=(a72+a70);
  a70=(a51*a32);
  a70=(a66-a70);
  a17=(a70*a64);
  a72=(a72+a17);
  a17=(a51*a7);
  a17=(a19-a17);
  a11=(a17*a49);
  a72=(a72+a11);
  a11=(a72*a48);
  a68=(a68-a11);
  a11=casadi_sq(a68);
  a9=(a72/a0);
  a9=(a35-a9);
  a8=casadi_sq(a9);
  a11=(a11+a8);
  a8=(a72*a45);
  a14=(a14-a8);
  a8=casadi_sq(a14);
  a11=(a11+a8);
  a8=(a69*a22);
  a8=(a61-a8);
  a34=(a29*a30);
  a8=(a8+a34);
  a34=casadi_sq(a8);
  a11=(a11+a34);
  a15=(a65*a15);
  a15=(a54+a15);
  a25=(a65*a25);
  a69=(a69*a25);
  a69=(a15-a69);
  a65=(a65*a2);
  a29=(a29*a65);
  a69=(a69-a29);
  a29=casadi_sq(a69);
  a11=(a11+a29);
  a29=(a72*a64);
  a70=(a70-a29);
  a29=casadi_sq(a70);
  a11=(a11+a29);
  a29=(a72*a49);
  a17=(a17-a29);
  a29=casadi_sq(a17);
  a11=(a11+a29);
  a29=casadi_sq(a53);
  a11=(a11+a29);
  a29=casadi_sq(a78);
  a11=(a11+a29);
  a11=sqrt(a11);
  a68=(a68/a11);
  a29=(a74*a68);
  a2=(a57/a0);
  a2=(a44-a2);
  a9=(a9/a11);
  a34=(a2*a9);
  a29=(a29+a34);
  a34=(a57*a45);
  a33=(a33-a34);
  a14=(a14/a11);
  a34=(a33*a14);
  a29=(a29+a34);
  a34=(a24*a22);
  a58=(a1*a30);
  a34=(a34+a58);
  a58=(a5*a34);
  a30=(a24*a30);
  a22=(a1*a22);
  a30=(a30-a22);
  a22=(a36*a30);
  a58=(a58+a22);
  a22=(a8/a11);
  a60=(a58*a22);
  a29=(a29+a60);
  a60=(a24*a25);
  a31=(a1*a65);
  a60=(a60-a31);
  a31=(a5*a60);
  a1=(a1*a25);
  a24=(a24*a65);
  a1=(a1+a24);
  a24=(a36*a1);
  a31=(a31-a24);
  a24=(a69/a11);
  a65=(a31*a24);
  a29=(a29+a65);
  a65=(a57*a64);
  a55=(a55-a65);
  a70=(a70/a11);
  a65=(a55*a70);
  a29=(a29+a65);
  a65=(a57*a49);
  a4=(a4-a65);
  a17=(a17/a11);
  a65=(a4*a17);
  a29=(a29+a65);
  a53=(a53/a11);
  a29=(a29+a53);
  a65=(a29*a68);
  a74=(a74-a65);
  a65=(a58/a50);
  a25=(a8*a39);
  a65=(a65+a25);
  a25=(a65/a50);
  a25=(a58-a25);
  a27=(a25*a48);
  a30=(a5*a30);
  a34=(a36*a34);
  a30=(a30-a34);
  a34=(a30/a0);
  a27=(a27+a34);
  a34=(a65*a39);
  a34=(a8-a34);
  a26=(a34*a45);
  a27=(a27+a26);
  a26=(a65*a32);
  a23=(a26*a64);
  a27=(a27-a23);
  a23=(a65*a7);
  a21=(a23*a49);
  a27=(a27-a21);
  a21=(a27*a48);
  a25=(a25-a21);
  a21=(a25*a68);
  a63=(a27/a0);
  a63=(a30-a63);
  a13=(a63*a9);
  a21=(a21+a13);
  a13=(a27*a45);
  a34=(a34-a13);
  a13=(a34*a14);
  a21=(a21+a13);
  a13=(a61*a22);
  a21=(a21+a13);
  a13=(a15*a24);
  a21=(a21+a13);
  a13=(a27*a64);
  a26=(a26+a13);
  a13=(a26*a70);
  a21=(a21-a13);
  a13=(a27*a49);
  a23=(a23+a13);
  a13=(a23*a17);
  a21=(a21-a13);
  a13=(a21*a68);
  a25=(a25-a13);
  a13=casadi_sq(a25);
  a10=(a21*a9);
  a63=(a63-a10);
  a10=casadi_sq(a63);
  a13=(a13+a10);
  a10=(a21*a14);
  a34=(a34-a10);
  a10=casadi_sq(a34);
  a13=(a13+a10);
  a10=(a21*a22);
  a10=(a61-a10);
  a77=casadi_sq(a10);
  a13=(a13+a77);
  a77=(a21*a24);
  a77=(a15-a77);
  a18=casadi_sq(a77);
  a13=(a13+a18);
  a18=(a21*a70);
  a26=(a26+a18);
  a18=casadi_sq(a26);
  a13=(a13+a18);
  a18=(a21*a17);
  a23=(a23+a18);
  a18=casadi_sq(a23);
  a13=(a13+a18);
  a18=(a21*a53);
  a12=casadi_sq(a18);
  a13=(a13+a12);
  a78=(a78/a11);
  a12=(a21*a78);
  a79=casadi_sq(a12);
  a13=(a13+a79);
  a13=sqrt(a13);
  a25=(a25/a13);
  a79=(a74*a25);
  a75=(a29*a9);
  a2=(a2-a75);
  a63=(a63/a13);
  a75=(a2*a63);
  a79=(a79+a75);
  a75=(a29*a14);
  a33=(a33-a75);
  a34=(a34/a13);
  a75=(a33*a34);
  a79=(a79+a75);
  a75=(a29*a22);
  a75=(a58-a75);
  a10=(a10/a13);
  a80=(a75*a10);
  a79=(a79+a80);
  a80=(a29*a24);
  a80=(a31-a80);
  a77=(a77/a13);
  a81=(a80*a77);
  a79=(a79+a81);
  a81=(a29*a70);
  a55=(a55-a81);
  a26=(a26/a13);
  a81=(a55*a26);
  a79=(a79-a81);
  a81=(a29*a17);
  a4=(a4-a81);
  a23=(a23/a13);
  a81=(a4*a23);
  a79=(a79-a81);
  a81=-1.;
  a82=(a29*a53);
  a82=(a81+a82);
  a18=(a18/a13);
  a83=(a82*a18);
  a79=(a79+a83);
  a83=(a29*a78);
  a12=(a12/a13);
  a84=(a83*a12);
  a79=(a79+a84);
  a84=(a79*a25);
  a74=(a74-a84);
  a84=(a31/a50);
  a85=(a69*a39);
  a84=(a84+a85);
  a85=(a84/a50);
  a85=(a31-a85);
  a86=(a85*a48);
  a60=(a36*a60);
  a1=(a5*a1);
  a60=(a60+a1);
  a1=(a60/a0);
  a86=(a86-a1);
  a1=(a84*a39);
  a1=(a69-a1);
  a87=(a1*a45);
  a86=(a86+a87);
  a87=(a84*a32);
  a88=(a87*a64);
  a86=(a86-a88);
  a88=(a84*a7);
  a89=(a88*a49);
  a86=(a86-a89);
  a89=(a86*a48);
  a85=(a85-a89);
  a89=(a85*a68);
  a90=(a86/a0);
  a90=(a60+a90);
  a91=(a90*a9);
  a89=(a89-a91);
  a91=(a86*a45);
  a1=(a1-a91);
  a91=(a1*a14);
  a89=(a89+a91);
  a91=(a15*a22);
  a89=(a89+a91);
  a91=(a54*a24);
  a89=(a89+a91);
  a91=(a86*a64);
  a87=(a87+a91);
  a91=(a87*a70);
  a89=(a89-a91);
  a91=(a86*a49);
  a88=(a88+a91);
  a91=(a88*a17);
  a89=(a89-a91);
  a91=(a89*a68);
  a85=(a85-a91);
  a91=(a85*a25);
  a92=(a89*a9);
  a90=(a90+a92);
  a92=(a90*a63);
  a91=(a91-a92);
  a92=(a89*a14);
  a1=(a1-a92);
  a92=(a1*a34);
  a91=(a91+a92);
  a92=(a89*a22);
  a92=(a15-a92);
  a93=(a92*a10);
  a91=(a91+a93);
  a93=(a89*a24);
  a93=(a54-a93);
  a94=(a93*a77);
  a91=(a91+a94);
  a94=(a89*a70);
  a87=(a87+a94);
  a94=(a87*a26);
  a91=(a91+a94);
  a94=(a89*a17);
  a88=(a88+a94);
  a94=(a88*a23);
  a91=(a91+a94);
  a94=(a89*a53);
  a95=(a94*a18);
  a91=(a91+a95);
  a95=(a89*a78);
  a96=(a95*a12);
  a91=(a91+a96);
  a96=(a91*a25);
  a85=(a85-a96);
  a96=casadi_sq(a85);
  a97=(a91*a63);
  a90=(a90+a97);
  a97=casadi_sq(a90);
  a96=(a96+a97);
  a97=(a91*a34);
  a1=(a1-a97);
  a97=casadi_sq(a1);
  a96=(a96+a97);
  a97=(a91*a10);
  a92=(a92-a97);
  a97=casadi_sq(a92);
  a96=(a96+a97);
  a97=(a91*a77);
  a93=(a93-a97);
  a97=casadi_sq(a93);
  a96=(a96+a97);
  a97=(a91*a26);
  a97=(a97-a87);
  a87=casadi_sq(a97);
  a96=(a96+a87);
  a87=(a91*a23);
  a87=(a87-a88);
  a88=casadi_sq(a87);
  a96=(a96+a88);
  a88=(a91*a18);
  a94=(a94-a88);
  a88=casadi_sq(a94);
  a96=(a96+a88);
  a88=(a91*a12);
  a88=(a88-a95);
  a95=casadi_sq(a88);
  a96=(a96+a95);
  a96=sqrt(a96);
  a85=(a85/a96);
  a95=(a74*a85);
  a98=(a79*a63);
  a2=(a2-a98);
  a90=(a90/a96);
  a98=(a2*a90);
  a95=(a95-a98);
  a98=(a79*a34);
  a33=(a33-a98);
  a1=(a1/a96);
  a98=(a33*a1);
  a95=(a95+a98);
  a98=(a79*a10);
  a75=(a75-a98);
  a92=(a92/a96);
  a98=(a75*a92);
  a95=(a95+a98);
  a98=(a79*a77);
  a80=(a80-a98);
  a93=(a93/a96);
  a98=(a80*a93);
  a95=(a95+a98);
  a98=(a79*a26);
  a55=(a55+a98);
  a97=(a97/a96);
  a98=(a55*a97);
  a95=(a95+a98);
  a98=(a79*a23);
  a4=(a4+a98);
  a87=(a87/a96);
  a98=(a4*a87);
  a95=(a95+a98);
  a98=(a79*a18);
  a82=(a82-a98);
  a94=(a94/a96);
  a98=(a82*a94);
  a95=(a95+a98);
  a98=(a79*a12);
  a98=(a98-a83);
  a88=(a88/a96);
  a83=(a98*a88);
  a95=(a95+a83);
  a83=(a95*a85);
  a74=(a74-a83);
  a83=(a40/a50);
  a99=(a66*a39);
  a83=(a83+a99);
  a99=(a3*a32);
  a83=(a83+a99);
  a99=(a67*a7);
  a83=(a83+a99);
  a99=(a83/a50);
  a99=(a40-a99);
  a100=(a99*a48);
  a37=(a5*a37);
  a42=(a36*a42);
  a37=(a37-a42);
  a42=(a37/a0);
  a100=(a100+a42);
  a42=(a83*a39);
  a42=(a66-a42);
  a101=(a42*a45);
  a100=(a100+a101);
  a101=(a83*a32);
  a101=(a3-a101);
  a102=(a101*a64);
  a100=(a100+a102);
  a102=(a83*a7);
  a102=(a67-a102);
  a103=(a102*a49);
  a100=(a100+a103);
  a103=(a100*a48);
  a99=(a99-a103);
  a103=(a99*a68);
  a104=(a100/a0);
  a104=(a37-a104);
  a105=(a104*a9);
  a103=(a103+a105);
  a105=(a100*a45);
  a42=(a42-a105);
  a105=(a42*a14);
  a103=(a103+a105);
  a105=(a100*a64);
  a101=(a101-a105);
  a105=(a101*a70);
  a103=(a103+a105);
  a105=(a100*a49);
  a102=(a102-a105);
  a105=(a102*a17);
  a103=(a103+a105);
  a105=(a62*a53);
  a103=(a103+a105);
  a105=(a47*a78);
  a103=(a103+a105);
  a105=(a103*a68);
  a99=(a99-a105);
  a105=(a99*a25);
  a106=(a103*a9);
  a104=(a104-a106);
  a106=(a104*a63);
  a105=(a105+a106);
  a106=(a103*a14);
  a42=(a42-a106);
  a106=(a42*a34);
  a105=(a105+a106);
  a106=(a103*a22);
  a107=(a106*a10);
  a105=(a105-a107);
  a107=(a103*a24);
  a108=(a107*a77);
  a105=(a105-a108);
  a108=(a103*a70);
  a101=(a101-a108);
  a108=(a101*a26);
  a105=(a105-a108);
  a108=(a103*a17);
  a102=(a102-a108);
  a108=(a102*a23);
  a105=(a105-a108);
  a108=(a103*a53);
  a108=(a108-a62);
  a62=(a108*a18);
  a105=(a105+a62);
  a62=(a103*a78);
  a47=(a47-a62);
  a62=(a47*a12);
  a105=(a105-a62);
  a62=(a105*a25);
  a99=(a99-a62);
  a62=(a99*a85);
  a109=(a105*a63);
  a104=(a104-a109);
  a109=(a104*a90);
  a62=(a62-a109);
  a109=(a105*a34);
  a42=(a42-a109);
  a109=(a42*a1);
  a62=(a62+a109);
  a109=(a105*a10);
  a106=(a106+a109);
  a109=(a106*a92);
  a62=(a62-a109);
  a109=(a105*a77);
  a107=(a107+a109);
  a109=(a107*a93);
  a62=(a62-a109);
  a109=(a105*a26);
  a101=(a101+a109);
  a109=(a101*a97);
  a62=(a62+a109);
  a109=(a105*a23);
  a102=(a102+a109);
  a109=(a102*a87);
  a62=(a62+a109);
  a109=(a105*a18);
  a108=(a108-a109);
  a109=(a108*a94);
  a62=(a62+a109);
  a109=(a105*a12);
  a47=(a47+a109);
  a109=(a47*a88);
  a62=(a62+a109);
  a109=(a62*a85);
  a99=(a99-a109);
  a109=casadi_sq(a99);
  a110=(a62*a90);
  a104=(a104+a110);
  a110=casadi_sq(a104);
  a109=(a109+a110);
  a110=(a62*a1);
  a42=(a42-a110);
  a110=casadi_sq(a42);
  a109=(a109+a110);
  a110=(a62*a92);
  a106=(a106+a110);
  a110=casadi_sq(a106);
  a109=(a109+a110);
  a110=(a62*a93);
  a107=(a107+a110);
  a110=casadi_sq(a107);
  a109=(a109+a110);
  a110=(a62*a97);
  a101=(a101-a110);
  a110=casadi_sq(a101);
  a109=(a109+a110);
  a110=(a62*a87);
  a102=(a102-a110);
  a110=casadi_sq(a102);
  a109=(a109+a110);
  a110=(a62*a94);
  a108=(a108-a110);
  a110=casadi_sq(a108);
  a109=(a109+a110);
  a110=(a62*a88);
  a47=(a47-a110);
  a110=casadi_sq(a47);
  a109=(a109+a110);
  a109=sqrt(a109);
  a99=(a99/a109);
  a110=(a74*a99);
  a111=(a95*a90);
  a2=(a2+a111);
  a104=(a104/a109);
  a111=(a2*a104);
  a110=(a110+a111);
  a111=(a95*a1);
  a33=(a33-a111);
  a42=(a42/a109);
  a111=(a33*a42);
  a110=(a110+a111);
  a111=(a95*a92);
  a75=(a75-a111);
  a106=(a106/a109);
  a111=(a75*a106);
  a110=(a110-a111);
  a111=(a95*a93);
  a80=(a80-a111);
  a107=(a107/a109);
  a111=(a80*a107);
  a110=(a110-a111);
  a111=(a95*a97);
  a55=(a55-a111);
  a101=(a101/a109);
  a111=(a55*a101);
  a110=(a110+a111);
  a111=(a95*a87);
  a4=(a4-a111);
  a102=(a102/a109);
  a111=(a4*a102);
  a110=(a110+a111);
  a111=(a95*a94);
  a82=(a82-a111);
  a108=(a108/a109);
  a111=(a82*a108);
  a110=(a110+a111);
  a111=(a95*a88);
  a98=(a98-a111);
  a47=(a47/a109);
  a111=(a98*a47);
  a110=(a110+a111);
  a111=(a110*a99);
  a74=(a74-a111);
  a111=(a73/a50);
  a112=(a19*a39);
  a111=(a111+a112);
  a112=(a67*a32);
  a111=(a111+a112);
  a112=(a54*a7);
  a111=(a111+a112);
  a112=(a111/a50);
  a112=(a73-a112);
  a113=(a112*a48);
  a71=(a36*a71);
  a16=(a5*a16);
  a71=(a71+a16);
  a16=(a71/a0);
  a113=(a113-a16);
  a16=(a111*a39);
  a16=(a19-a16);
  a114=(a16*a45);
  a113=(a113+a114);
  a114=(a111*a32);
  a114=(a67-a114);
  a115=(a114*a64);
  a113=(a113+a115);
  a115=(a111*a7);
  a115=(a54-a115);
  a116=(a115*a49);
  a113=(a113+a116);
  a116=(a113*a48);
  a112=(a112-a116);
  a116=(a112*a68);
  a117=(a113/a0);
  a117=(a71+a117);
  a118=(a117*a9);
  a116=(a116-a118);
  a118=(a113*a45);
  a16=(a16-a118);
  a118=(a16*a14);
  a116=(a116+a118);
  a118=(a113*a64);
  a114=(a114-a118);
  a118=(a114*a70);
  a116=(a116+a118);
  a118=(a113*a49);
  a115=(a115-a118);
  a118=(a115*a17);
  a116=(a116+a118);
  a118=(a59*a53);
  a116=(a116+a118);
  a118=(a46*a78);
  a116=(a116+a118);
  a118=(a116*a68);
  a112=(a112-a118);
  a118=(a112*a25);
  a119=(a116*a9);
  a117=(a117+a119);
  a119=(a117*a63);
  a118=(a118-a119);
  a119=(a116*a14);
  a16=(a16-a119);
  a119=(a16*a34);
  a118=(a118+a119);
  a119=(a116*a22);
  a120=(a119*a10);
  a118=(a118-a120);
  a120=(a116*a24);
  a121=(a120*a77);
  a118=(a118-a121);
  a121=(a116*a70);
  a114=(a114-a121);
  a121=(a114*a26);
  a118=(a118-a121);
  a121=(a116*a17);
  a115=(a115-a121);
  a121=(a115*a23);
  a118=(a118-a121);
  a121=(a116*a53);
  a121=(a121-a59);
  a59=(a121*a18);
  a118=(a118+a59);
  a59=(a116*a78);
  a46=(a46-a59);
  a59=(a46*a12);
  a118=(a118-a59);
  a59=(a118*a25);
  a112=(a112-a59);
  a59=(a112*a85);
  a122=(a118*a63);
  a117=(a117+a122);
  a122=(a117*a90);
  a59=(a59+a122);
  a122=(a118*a34);
  a16=(a16-a122);
  a122=(a16*a1);
  a59=(a59+a122);
  a122=(a118*a10);
  a119=(a119+a122);
  a122=(a119*a92);
  a59=(a59-a122);
  a122=(a118*a77);
  a120=(a120+a122);
  a122=(a120*a93);
  a59=(a59-a122);
  a122=(a118*a26);
  a114=(a114+a122);
  a122=(a114*a97);
  a59=(a59+a122);
  a122=(a118*a23);
  a115=(a115+a122);
  a122=(a115*a87);
  a59=(a59+a122);
  a122=(a118*a18);
  a121=(a121-a122);
  a122=(a121*a94);
  a59=(a59+a122);
  a122=(a118*a12);
  a46=(a46+a122);
  a122=(a46*a88);
  a59=(a59+a122);
  a122=(a59*a85);
  a112=(a112-a122);
  a122=(a112*a99);
  a123=(a59*a90);
  a123=(a123-a117);
  a117=(a123*a104);
  a122=(a122+a117);
  a117=(a59*a1);
  a16=(a16-a117);
  a117=(a16*a42);
  a122=(a122+a117);
  a117=(a59*a92);
  a119=(a119+a117);
  a117=(a119*a106);
  a122=(a122+a117);
  a117=(a59*a93);
  a120=(a120+a117);
  a117=(a120*a107);
  a122=(a122+a117);
  a117=(a59*a97);
  a114=(a114-a117);
  a117=(a114*a101);
  a122=(a122+a117);
  a117=(a59*a87);
  a115=(a115-a117);
  a117=(a115*a102);
  a122=(a122+a117);
  a117=(a59*a94);
  a121=(a121-a117);
  a117=(a121*a108);
  a122=(a122+a117);
  a117=(a59*a88);
  a46=(a46-a117);
  a117=(a46*a47);
  a122=(a122+a117);
  a117=(a122*a99);
  a112=(a112-a117);
  a117=casadi_sq(a112);
  a124=(a122*a104);
  a123=(a123-a124);
  a124=casadi_sq(a123);
  a117=(a117+a124);
  a124=(a122*a42);
  a16=(a16-a124);
  a124=casadi_sq(a16);
  a117=(a117+a124);
  a124=(a122*a106);
  a124=(a124-a119);
  a119=casadi_sq(a124);
  a117=(a117+a119);
  a119=(a122*a107);
  a119=(a119-a120);
  a120=casadi_sq(a119);
  a117=(a117+a120);
  a120=(a122*a101);
  a114=(a114-a120);
  a120=casadi_sq(a114);
  a117=(a117+a120);
  a120=(a122*a102);
  a115=(a115-a120);
  a120=casadi_sq(a115);
  a117=(a117+a120);
  a120=(a122*a108);
  a121=(a121-a120);
  a120=casadi_sq(a121);
  a117=(a117+a120);
  a120=(a122*a47);
  a46=(a46-a120);
  a120=casadi_sq(a46);
  a117=(a117+a120);
  a117=sqrt(a117);
  a112=(a112/a117);
  a120=(a74*a112);
  a125=(a110*a104);
  a2=(a2-a125);
  a123=(a123/a117);
  a125=(a2*a123);
  a120=(a120+a125);
  a125=(a110*a42);
  a33=(a33-a125);
  a16=(a16/a117);
  a125=(a33*a16);
  a120=(a120+a125);
  a125=(a110*a106);
  a75=(a75+a125);
  a124=(a124/a117);
  a125=(a75*a124);
  a120=(a120+a125);
  a125=(a110*a107);
  a80=(a80+a125);
  a119=(a119/a117);
  a125=(a80*a119);
  a120=(a120+a125);
  a125=(a110*a101);
  a55=(a55-a125);
  a114=(a114/a117);
  a125=(a55*a114);
  a120=(a120+a125);
  a125=(a110*a102);
  a4=(a4-a125);
  a115=(a115/a117);
  a125=(a4*a115);
  a120=(a120+a125);
  a125=(a110*a108);
  a82=(a82-a125);
  a121=(a121/a117);
  a125=(a82*a121);
  a120=(a120+a125);
  a125=(a110*a47);
  a98=(a98-a125);
  a46=(a46/a117);
  a125=(a98*a46);
  a120=(a120+a125);
  a125=(a120*a112);
  a74=(a74-a125);
  a125=casadi_sq(a74);
  a126=(a120*a123);
  a2=(a2-a126);
  a126=casadi_sq(a2);
  a125=(a125+a126);
  a126=(a120*a16);
  a33=(a33-a126);
  a126=casadi_sq(a33);
  a125=(a125+a126);
  a126=(a120*a124);
  a75=(a75-a126);
  a126=casadi_sq(a75);
  a125=(a125+a126);
  a126=(a120*a119);
  a80=(a80-a126);
  a126=casadi_sq(a80);
  a125=(a125+a126);
  a126=(a120*a114);
  a55=(a55-a126);
  a126=casadi_sq(a55);
  a125=(a125+a126);
  a126=(a120*a115);
  a4=(a4-a126);
  a126=casadi_sq(a4);
  a125=(a125+a126);
  a126=(a120*a121);
  a82=(a82-a126);
  a126=casadi_sq(a82);
  a125=(a125+a126);
  a126=(a120*a46);
  a98=(a98-a126);
  a126=casadi_sq(a98);
  a125=(a125+a126);
  a125=sqrt(a125);
  a74=(a74/a125);
  a126=arg[0]? arg[0][7] : 0;
  a41=(a41*a126);
  a127=arg[0]? arg[0][8] : 0;
  a128=(a44*a127);
  a41=(a41+a128);
  a128=arg[0]? arg[0][9] : 0;
  a129=(a28*a128);
  a41=(a41+a129);
  a129=arg[0]? arg[0][10] : 0;
  a130=(a58*a129);
  a41=(a41+a130);
  a130=arg[0]? arg[0][11] : 0;
  a131=(a31*a130);
  a41=(a41+a131);
  a131=arg[0]? arg[0][12] : 0;
  a132=(a40*a131);
  a41=(a41+a132);
  a132=arg[0]? arg[0][13] : 0;
  a133=(a73*a132);
  a41=(a41+a133);
  a133=(a74*a41);
  a2=(a2/a125);
  a134=(a44*a126);
  a6=(a5*a6);
  a43=(a36*a43);
  a6=(a6-a43);
  a6=(a6*a5);
  a5=(a5*a38);
  a20=(a36*a20);
  a5=(a5-a20);
  a5=(a5*a36);
  a6=(a6-a5);
  a5=(a6*a127);
  a134=(a134+a5);
  a5=(a35*a128);
  a134=(a134+a5);
  a5=(a30*a129);
  a134=(a134+a5);
  a5=(a60*a130);
  a134=(a134-a5);
  a5=(a37*a131);
  a134=(a134+a5);
  a5=(a71*a132);
  a134=(a134-a5);
  a5=(a2*a134);
  a133=(a133+a5);
  a33=(a33/a125);
  a28=(a28*a126);
  a5=(a35*a127);
  a28=(a28+a5);
  a56=(a56*a128);
  a28=(a28+a56);
  a56=(a8*a129);
  a28=(a28+a56);
  a56=(a69*a130);
  a28=(a28+a56);
  a56=(a66*a131);
  a28=(a28+a56);
  a56=(a19*a132);
  a28=(a28+a56);
  a56=(a33*a28);
  a133=(a133+a56);
  a75=(a75/a125);
  a58=(a58*a126);
  a56=(a30*a127);
  a58=(a58+a56);
  a8=(a8*a128);
  a58=(a58+a8);
  a61=(a61*a129);
  a58=(a58+a61);
  a61=(a15*a130);
  a58=(a58+a61);
  a61=(a75*a58);
  a133=(a133+a61);
  a80=(a80/a125);
  a31=(a31*a126);
  a61=(a60*a127);
  a31=(a31-a61);
  a69=(a69*a128);
  a31=(a31+a69);
  a15=(a15*a129);
  a31=(a31+a15);
  a130=(a54*a130);
  a31=(a31+a130);
  a130=(a80*a31);
  a133=(a133+a130);
  a55=(a55/a125);
  a40=(a40*a126);
  a130=(a37*a127);
  a40=(a40+a130);
  a66=(a66*a128);
  a40=(a40+a66);
  a3=(a3*a131);
  a40=(a40+a3);
  a3=(a67*a132);
  a40=(a40+a3);
  a3=(a55*a40);
  a133=(a133+a3);
  a4=(a4/a125);
  a73=(a73*a126);
  a127=(a71*a127);
  a73=(a73-a127);
  a19=(a19*a128);
  a73=(a73+a19);
  a67=(a67*a131);
  a73=(a73+a67);
  a54=(a54*a132);
  a73=(a73+a54);
  a54=(a4*a73);
  a133=(a133+a54);
  a54=(a44/a50);
  a132=(a35*a39);
  a54=(a54+a132);
  a132=(a37*a32);
  a54=(a54+a132);
  a132=(a71*a7);
  a54=(a54-a132);
  a132=(a54/a50);
  a44=(a44-a132);
  a132=(a44*a48);
  a67=(a6/a0);
  a132=(a132+a67);
  a67=(a54*a39);
  a35=(a35-a67);
  a67=(a35*a45);
  a132=(a132+a67);
  a67=(a54*a32);
  a37=(a37-a67);
  a67=(a37*a64);
  a132=(a132+a67);
  a67=(a54*a7);
  a71=(a71+a67);
  a67=(a71*a49);
  a132=(a132-a67);
  a67=(a132*a48);
  a44=(a44-a67);
  a67=(a44*a68);
  a131=(a132/a0);
  a6=(a6-a131);
  a131=(a6*a9);
  a67=(a67+a131);
  a131=(a132*a45);
  a35=(a35-a131);
  a131=(a35*a14);
  a67=(a67+a131);
  a131=(a30*a22);
  a67=(a67+a131);
  a131=(a60*a24);
  a67=(a67-a131);
  a131=(a132*a64);
  a37=(a37-a131);
  a131=(a37*a70);
  a67=(a67+a131);
  a131=(a132*a49);
  a71=(a71+a131);
  a131=(a71*a17);
  a67=(a67-a131);
  a67=(a67-a78);
  a131=(a67*a68);
  a44=(a44-a131);
  a131=(a44*a25);
  a19=(a67*a9);
  a6=(a6-a19);
  a19=(a6*a63);
  a131=(a131+a19);
  a19=(a67*a14);
  a35=(a35-a19);
  a19=(a35*a34);
  a131=(a131+a19);
  a19=(a67*a22);
  a30=(a30-a19);
  a19=(a30*a10);
  a131=(a131+a19);
  a19=(a67*a24);
  a60=(a60+a19);
  a19=(a60*a77);
  a131=(a131-a19);
  a19=(a67*a70);
  a37=(a37-a19);
  a19=(a37*a26);
  a131=(a131-a19);
  a19=(a67*a17);
  a71=(a71+a19);
  a19=(a71*a23);
  a131=(a131+a19);
  a53=(a67*a53);
  a19=(a53*a18);
  a131=(a131+a19);
  a78=(a67*a78);
  a81=(a81-a78);
  a78=(a81*a12);
  a131=(a131-a78);
  a78=(a131*a25);
  a44=(a44-a78);
  a78=(a44*a85);
  a19=(a131*a63);
  a6=(a6-a19);
  a19=(a6*a90);
  a78=(a78-a19);
  a19=(a131*a34);
  a35=(a35-a19);
  a19=(a35*a1);
  a78=(a78+a19);
  a19=(a131*a10);
  a30=(a30-a19);
  a19=(a30*a92);
  a78=(a78+a19);
  a19=(a131*a77);
  a60=(a60+a19);
  a19=(a60*a93);
  a78=(a78-a19);
  a19=(a131*a26);
  a37=(a37+a19);
  a19=(a37*a97);
  a78=(a78+a19);
  a19=(a131*a23);
  a19=(a19-a71);
  a71=(a19*a87);
  a78=(a78+a71);
  a18=(a131*a18);
  a53=(a53-a18);
  a18=(a53*a94);
  a78=(a78+a18);
  a12=(a131*a12);
  a81=(a81+a12);
  a12=(a81*a88);
  a78=(a78+a12);
  a12=(a78*a85);
  a44=(a44-a12);
  a12=(a44*a99);
  a18=(a78*a90);
  a6=(a6+a18);
  a18=(a6*a104);
  a12=(a12+a18);
  a18=(a78*a1);
  a35=(a35-a18);
  a18=(a35*a42);
  a12=(a12+a18);
  a18=(a78*a92);
  a30=(a30-a18);
  a18=(a30*a106);
  a12=(a12-a18);
  a18=(a78*a93);
  a60=(a60+a18);
  a18=(a60*a107);
  a12=(a12+a18);
  a18=(a78*a97);
  a37=(a37-a18);
  a18=(a37*a101);
  a12=(a12+a18);
  a18=(a78*a87);
  a19=(a19-a18);
  a18=(a19*a102);
  a12=(a12+a18);
  a94=(a78*a94);
  a53=(a53-a94);
  a94=(a53*a108);
  a12=(a12+a94);
  a88=(a78*a88);
  a81=(a81-a88);
  a88=(a81*a47);
  a12=(a12+a88);
  a88=(a12*a99);
  a44=(a44-a88);
  a88=(a44*a112);
  a94=(a12*a104);
  a6=(a6-a94);
  a94=(a6*a123);
  a88=(a88+a94);
  a94=(a12*a42);
  a35=(a35-a94);
  a94=(a35*a16);
  a88=(a88+a94);
  a94=(a12*a106);
  a30=(a30+a94);
  a94=(a30*a124);
  a88=(a88+a94);
  a94=(a12*a107);
  a94=(a94-a60);
  a60=(a94*a119);
  a88=(a88+a60);
  a60=(a12*a101);
  a37=(a37-a60);
  a60=(a37*a114);
  a88=(a88+a60);
  a60=(a12*a102);
  a19=(a19-a60);
  a60=(a19*a115);
  a88=(a88+a60);
  a108=(a12*a108);
  a53=(a53-a108);
  a108=(a53*a121);
  a88=(a88+a108);
  a47=(a12*a47);
  a81=(a81-a47);
  a47=(a81*a46);
  a88=(a88+a47);
  a47=(a88*a112);
  a44=(a44-a47);
  a47=(a44*a74);
  a108=(a88*a123);
  a6=(a6-a108);
  a108=(a6*a2);
  a47=(a47+a108);
  a108=(a88*a16);
  a35=(a35-a108);
  a108=(a35*a33);
  a47=(a47+a108);
  a108=(a88*a124);
  a30=(a30-a108);
  a108=(a30*a75);
  a47=(a47+a108);
  a108=(a88*a119);
  a94=(a94-a108);
  a108=(a94*a80);
  a47=(a47+a108);
  a108=(a88*a114);
  a37=(a37-a108);
  a108=(a37*a55);
  a47=(a47+a108);
  a108=(a88*a115);
  a19=(a19-a108);
  a108=(a19*a4);
  a47=(a47+a108);
  a121=(a88*a121);
  a53=(a53-a121);
  a82=(a82/a125);
  a121=(a53*a82);
  a47=(a47+a121);
  a46=(a88*a46);
  a81=(a81-a46);
  a98=(a98/a125);
  a46=(a81*a98);
  a47=(a47+a46);
  a74=(a47*a74);
  a44=(a44-a74);
  a74=casadi_sq(a44);
  a2=(a47*a2);
  a6=(a6-a2);
  a2=casadi_sq(a6);
  a74=(a74+a2);
  a33=(a47*a33);
  a35=(a35-a33);
  a33=casadi_sq(a35);
  a74=(a74+a33);
  a75=(a47*a75);
  a30=(a30-a75);
  a75=casadi_sq(a30);
  a74=(a74+a75);
  a80=(a47*a80);
  a94=(a94-a80);
  a80=casadi_sq(a94);
  a74=(a74+a80);
  a55=(a47*a55);
  a37=(a37-a55);
  a55=casadi_sq(a37);
  a74=(a74+a55);
  a4=(a47*a4);
  a19=(a19-a4);
  a4=casadi_sq(a19);
  a74=(a74+a4);
  a82=(a47*a82);
  a53=(a53-a82);
  a53=casadi_sq(a53);
  a74=(a74+a53);
  a98=(a47*a98);
  a81=(a81-a98);
  a81=casadi_sq(a81);
  a74=(a74+a81);
  a74=sqrt(a74);
  a44=(a44/a74);
  a44=(a44*a41);
  a6=(a6/a74);
  a6=(a6*a134);
  a44=(a44+a6);
  a35=(a35/a74);
  a35=(a35*a28);
  a44=(a44+a35);
  a30=(a30/a74);
  a30=(a30*a58);
  a44=(a44+a30);
  a94=(a94/a74);
  a94=(a94*a31);
  a44=(a44+a94);
  a37=(a37/a74);
  a37=(a37*a40);
  a44=(a44+a37);
  a19=(a19/a74);
  a19=(a19*a73);
  a44=(a44+a19);
  a44=(a44/a74);
  a47=(a47*a44);
  a133=(a133-a47);
  a133=(a133/a125);
  if (res[0]!=0) res[0][7]=a133;
  if (res[0]!=0) res[0][8]=a44;
  a68=(a68*a41);
  a9=(a9*a134);
  a68=(a68+a9);
  a14=(a14*a28);
  a68=(a68+a14);
  a22=(a22*a58);
  a68=(a68+a22);
  a24=(a24*a31);
  a68=(a68+a24);
  a70=(a70*a40);
  a68=(a68+a70);
  a17=(a17*a73);
  a68=(a68+a17);
  a67=(a67*a44);
  a68=(a68-a67);
  a29=(a29*a133);
  a68=(a68-a29);
  a112=(a112*a41);
  a123=(a123*a134);
  a112=(a112+a123);
  a16=(a16*a28);
  a112=(a112+a16);
  a124=(a124*a58);
  a112=(a112+a124);
  a119=(a119*a31);
  a112=(a112+a119);
  a114=(a114*a40);
  a112=(a112+a114);
  a115=(a115*a73);
  a112=(a112+a115);
  a88=(a88*a44);
  a112=(a112-a88);
  a120=(a120*a133);
  a112=(a112-a120);
  a112=(a112/a117);
  a116=(a116*a112);
  a68=(a68-a116);
  a99=(a99*a41);
  a104=(a104*a134);
  a99=(a99+a104);
  a42=(a42*a28);
  a99=(a99+a42);
  a106=(a106*a58);
  a99=(a99-a106);
  a107=(a107*a31);
  a99=(a99-a107);
  a101=(a101*a40);
  a99=(a99+a101);
  a102=(a102*a73);
  a99=(a99+a102);
  a12=(a12*a44);
  a99=(a99-a12);
  a110=(a110*a133);
  a99=(a99-a110);
  a122=(a122*a112);
  a99=(a99-a122);
  a99=(a99/a109);
  a103=(a103*a99);
  a68=(a68-a103);
  a85=(a85*a41);
  a90=(a90*a134);
  a85=(a85-a90);
  a1=(a1*a28);
  a85=(a85+a1);
  a92=(a92*a58);
  a85=(a85+a92);
  a93=(a93*a31);
  a85=(a85+a93);
  a97=(a97*a40);
  a85=(a85+a97);
  a87=(a87*a73);
  a85=(a85+a87);
  a78=(a78*a44);
  a85=(a85-a78);
  a95=(a95*a133);
  a85=(a85-a95);
  a59=(a59*a112);
  a85=(a85-a59);
  a62=(a62*a99);
  a85=(a85-a62);
  a85=(a85/a96);
  a89=(a89*a85);
  a68=(a68-a89);
  a25=(a25*a41);
  a63=(a63*a134);
  a25=(a25+a63);
  a34=(a34*a28);
  a25=(a25+a34);
  a10=(a10*a58);
  a25=(a25+a10);
  a77=(a77*a31);
  a25=(a25+a77);
  a26=(a26*a40);
  a25=(a25-a26);
  a23=(a23*a73);
  a25=(a25-a23);
  a131=(a131*a44);
  a25=(a25-a131);
  a79=(a79*a133);
  a25=(a25-a79);
  a118=(a118*a112);
  a25=(a25-a118);
  a105=(a105*a99);
  a25=(a25-a105);
  a91=(a91*a85);
  a25=(a25-a91);
  a25=(a25/a13);
  a21=(a21*a25);
  a68=(a68-a21);
  a68=(a68/a11);
  if (res[0]!=0) res[0][9]=a68;
  if (res[0]!=0) res[0][10]=a25;
  if (res[0]!=0) res[0][11]=a85;
  if (res[0]!=0) res[0][12]=a99;
  if (res[0]!=0) res[0][13]=a112;
  a11=(a41/a50);
  a39=(a39*a28);
  a11=(a11+a39);
  a32=(a32*a40);
  a11=(a11+a32);
  a7=(a7*a73);
  a11=(a11+a7);
  a54=(a54*a44);
  a11=(a11-a54);
  a52=(a52*a133);
  a11=(a11-a52);
  a111=(a111*a112);
  a11=(a11-a111);
  a83=(a83*a99);
  a11=(a11-a83);
  a84=(a84*a85);
  a11=(a11-a84);
  a65=(a65*a25);
  a11=(a11-a65);
  a51=(a51*a68);
  a11=(a11-a51);
  a48=(a48*a41);
  a134=(a134/a0);
  a48=(a48+a134);
  a45=(a45*a28);
  a48=(a48+a45);
  a64=(a64*a40);
  a48=(a48+a64);
  a49=(a49*a73);
  a48=(a48+a49);
  a132=(a132*a44);
  a48=(a48-a132);
  a57=(a57*a133);
  a48=(a48-a57);
  a113=(a113*a112);
  a48=(a48-a113);
  a100=(a100*a99);
  a48=(a48-a100);
  a86=(a86*a85);
  a48=(a48-a86);
  a27=(a27*a25);
  a48=(a48-a27);
  a72=(a72*a68);
  a48=(a48-a72);
  a48=(a48/a0);
  a76=(a76*a48);
  a11=(a11+a76);
  a11=(a11/a50);
  if (res[1]!=0) res[1][0]=a11;
  if (res[1]!=0) res[1][1]=a48;
  return 0;
}

CASADI_SYMBOL_EXPORT int Imp_B(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int Imp_B_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int Imp_B_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void Imp_B_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int Imp_B_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void Imp_B_release(int mem) {
}

CASADI_SYMBOL_EXPORT void Imp_B_incref(void) {
}

CASADI_SYMBOL_EXPORT void Imp_B_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int Imp_B_n_in(void) { return 2;}

CASADI_SYMBOL_EXPORT casadi_int Imp_B_n_out(void) { return 2;}

CASADI_SYMBOL_EXPORT casadi_real Imp_B_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* Imp_B_name_in(casadi_int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* Imp_B_name_out(casadi_int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* Imp_B_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* Imp_B_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int Imp_B_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 2;
  if (sz_res) *sz_res = 2;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
