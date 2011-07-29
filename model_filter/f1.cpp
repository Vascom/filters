/**************************************************************************** \
*                   Copyright (c) 2006 by Javad GNSS.                        *
*                          All Rights Reserved.                              *
******************************************************************************

Data : 01.03.2006

Author : Glazov Vasiliy

Contents: main program of filter modeling

Description: these program contein model of FIR-filters and signal generator
for testing

\****************************************************************************/
// Filter.cpp : Defines the entry point for the console application.
//#define allow_iir_filter 0    // 0 - no, 1 - yes
#define CIC_filter 1             // 0 - ordinary FIR-filter, 1 - CIC-filter
#define adc_agc_mode 0          // 0 - wide, 1 - narrow
#define anti_jamm_on 1           // 0 - off, 1 - on
#define spectrum_analyser_on 0  // 0 - off, 1 - on
#define ca_p_code 1             // 1 - C/A, 10 - P
#define pwm_dac 1               // 0 - PWM, 1 - DAC
#define ajm_real 1              // 0 - complex, 1 - real
//#define linux
#define TRESHOLD_SWITCH_ALGORITHM 0
#define LOCATA 0
#define AVT2 0

#define PI 3.1415926538
#define V5_COMPAT
//#define linux
//#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
//#include <string.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef linux
#include "slstrn.h"
//#include "slrnrc.h"
#include "slru.h"
#else
#include <slrn.h>
#include <slrnrc.h>
#include <slru.h>
#endif
#include <math.h>
#include <algorithm>
#include "class_Filter.h"
//#include "class_ajm.h"

using namespace std;

double *FFTw( int *dIn, int nn )
{
    int i, j, n, m, mmax, istep;
    double tempr, tempi, wtemp, theta, wpr, wpi, wr, wi;

    int isign = -1;
    double *data = new double [ nn * 2 + 1 ];

    for( i = 0; i < nn; i++ )
    {
        data[ i * 2 ] = 0;
        data[ i * 2 + 1 ] = dIn[ i ]; //cout<<data[ i*2 ]<<' '<<data[ i*2+1 ]<<endl;
    }
	data[512] = 0;
for( short i = 0; i <= nn*2; i++ )
    {
        //cout<<i<<' '<<data[ i ]<<endl;
    }
    n = nn << 1;
    j = 1;
    i = 1;
    while( i < n )
    {
        if( j > i )
        {
            tempr = data[ i ]; data[ i ] = data[ j ]; data[ j ] = tempr;
            tempr = data[ i + 1 ]; data[ i + 1 ] = data[ j + 1 ]; data[ j + 1 ] = tempr;
        }
        m = n >> 1;
        while( ( m >= 2 ) && ( j > m ) )
        {
            j = j - m;
            m = m >> 1;
        }
        j = j + m;
        i = i + 2;//cout<<i<<' '<<n<<' '<<data[ i ]<<endl;
    }
	for( short i = 0; i <= nn*2; i++ )
    {
        //cout<<i<<' '<<data[ i ]<<endl;
    }
	
    mmax = 2;
    while( n > mmax )
    {
        istep = 2 * mmax;
        theta = 2.0 * PI / ( isign * mmax );
        wtemp = sin( 0.5 * theta );
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin( theta );
        wr = 1.0;
        wi = 0.0;
        m = 1;
		for( short  i = 0; i < nn; i++ )
    {
        //cout<<istep<<' '<<theta<<' '<<wtemp<<' '<<wpr<<' ' <<wpi<<endl;
    }
        while( m < mmax )
        {
            i = m;
            while( i < n )
            {
                j = i + mmax;
                tempr = wr * data[ j ] - wi * data[ j + 1 ];
                tempi = wr * data[ j + 1 ] + wi * data[ j ];
                data[ j ] = data[ i ] - tempr;
                data[ j + 1 ] = data[ i + 1 ] - tempi;
                data[ i ] = data[ i ] + tempr;
                data[ i + 1 ] = data[ i + 1 ] + tempi;//cout<<m<<' '<<i<<' '<<j<<' '<<data[ i ]<<' '<<data[ i+1 ]<<' '<<data[ j ]<<' '<<data[ j+1 ]<<' '<<tempr<<' '<<tempi<<endl;
                i = i + istep;
				
            }
            wtemp = wr;
            wr = wtemp * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
            m = m + 2;
        }
        mmax = istep;
    }
    double *dOut = new double [ nn / 2 ];

    for( i = 0; i < ( nn / 2 ); i++ )
    {
        dOut[ i ] = sqrt( data[ (i+1) * 2 ] * data[ (i+1) * 2 ] + data[ i * 2 + 1 ] * data[ i * 2 + 1 ] );//cout<<data[ i ]<<endl;
    }

    delete []data;
    return dOut;
}
// keep sign bits
int keep_sign(int num, int width)
{
    for (int i = width; i < 32; i++){ num = (num | (((num >> (i-1)) & 0x1) << i)); }
    return num;
}
//Generator E5aI PRN
void GALE5aIgen(short *e5a)
{
 short geng1[14]={1,1,1,1,1,1,1,1,1,1,1,1,1,1};
 short geng2[14]={0,1,1,0,1,1,0,0,1,0,1,0,1,0};
 short G=0,g1=0,g2=0,g1s=0,g2s=0,x=0;
 for(short n=0;n<10000;n++)
 {
  g1=geng1[13];
  g1s=geng1[0]^geng1[5]^geng1[7]^geng1[13];
  //for(x=0;x<14;x++){geng1[14-1-x]=geng1[14-2-x];}
  rotate(&geng1[0], &geng1[13], &geng1[14]);
  geng1[0]=g1s;
  g2=geng2[13];
  g2s=geng2[3]^geng2[4]^geng2[6]^geng2[7]^geng2[11]^geng2[13];
  //for(x=0;x<14;x++){geng2[14-1-x]=geng2[14-2-x];}
  rotate(&geng2[0], &geng2[13], &geng2[14]);
  geng2[0]=g2s;
  G=g1^g2;
  if(G==0) G=-1;
  e5a[n]=G;
 }
}
//Generator E5bI PRN
void GALE5bIgen(short *e5b)
{
 short geng1[14]={1,1,1,1,1,1,1,1,1,1,1,1,1,1};
 short geng2[14]={1,1,1,0,0,0,1,0,1,0,0,0,1,0};
 short G=0,g1=0,g2=0,g1s=0,g2s=0,x=0;
 for(short n=0;n<10000;n++)
 {
  g1=geng1[13];
  g1s=geng1[3]^geng1[10]^geng1[12]^geng1[13];
  //for(x=0;x<14;x++){geng1[14-1-x]=geng1[14-2-x];}
  rotate(&geng1[0], &geng1[13], &geng1[14]);
  geng1[0]=g1s;
  g2=geng2[13];
  g2s=geng2[1]^geng2[4]^geng2[7]^geng2[8]^geng2[11]^geng2[13];
  //for(x=0;x<14;x++){geng2[14-1-x]=geng2[14-2-x];}
  rotate(&geng2[0], &geng2[13], &geng2[14]);
  geng2[0]=g2s;
  G=g1^g2;
  if(G==0) G=-1;
  e5b[n]=G;
 }
}

//Make two impulses for multipath error for I channel
short RadioImpI(unsigned long t)
{
 if(t<120)//if(t<120)
 {
  return 0;
  //double a=-cos(4e6*2*pi*t/180e6);
  //return short(floor(a*(1<<(8-1))+0.5));
 } else
 if(t>=120 && t<300)//if(t>=540 && t<1260)
 {
  return 1;
  //double a=cos(4e6*2*pi*t/180e6);
  //return short(floor(a*(1<<(8-1))+0.5));
 } else
 if(t>=300 && t<480)
 {
  return -1;
  //double a=-cos(4e6*2*pi*t/180e6);
  //return short(floor(a*(1<<(8-1))+0.5));
 } else
  if(t>=480)return 0; else return 0;
}

//Make two impulses for multipath error for Q channel
short RadioImpQ(unsigned long t)
{
 if(t<120)//if(t<540)
 {
  return 0;
  //double a=sin(4e6*2*pi*t/180e6);
  //return short(floor(a*(1<<(8-1))+0.5));
 } else
 if(t>=120 && t<138)//if(t>=540 && t<1260)
 {
  return 1;
  //double a=sin(4e6*2*pi*t/180e6);
  //return short(floor(a*(1<<(8-1))+0.5));
 } else
 if(t>=138 && t<156)
 {
  return -1;
  //double a=-sin(4e6*2*pi*t/180e6);
  //return short(floor(a*(1<<(8-1))+0.5));
 } else
  if(t>=156)return 0; else return 0;
}

//Input PRN-code from file
void InputPosled(short inp[])
{
 ifstream input1;
 input1.open("data_in/g.txt");//("posl_1.txt");
 if(input1.fail())cout<<"file g.txt not open";
 for(short i=0;i<1023;i++)
 {
  input1>>inp[i];
 }
 input1.close();
}
void InputPosledGlo(short inp[])
{
 ifstream input1;
 input1.open("data_in/glo_ca.txt");//("posl_1.txt");
 if(input1.fail())cout<<"file not open";
 for(short i=0;i<511;i++)
 {
  input1>>inp[i];
 }
 input1.close();
}

void InputVerilogPosled(int inp[])
{
 ifstream input1;
 input1.open("data_out/rf_in.txt");//("posl_1.txt");
 if(input1.fail())cout<<"file not open";
 for(short i=0;i<1700;i++)
 {
  input1>>inp[i];
 }
 input1.close();
}

void InputVerilogPosled2(int inp[])
{
 ifstream input1;
 input1.open("data_out/nco_out2.txt");//("posl_1.txt");
 if(input1.fail())cout<<"file not open";
 for(short i=0;i<1700;i++)
 {
  input1>>inp[i];
 }
 input1.close();
}

void InputAVT2data(short inp[])
{
    ifstream input1;
    input1.open("data_in/avt2.txt");//("posl_1.txt");
    if(input1.fail())cout<<"file not open";
    for(int i=0;i<65536;i++)
    {
        input1>>inp[i];
    }
    input1.close();
}
// cos generator
//#define nco_bit_out 8
short NCO_gen_cos(double freq, long long t)
{
 //return short(floor(cos(double(t) * freq)*127 + 0.5));
 double phase = (((long long)floor(double(t) * freq / 2.0 / pi * 16384 + 0.5) % 16384) >> 4) * 2.0 * pi / 1024.0;
 if (((cos(phase) * sin(phase)) >= 0) & ((cos(phase + 2.0 * pi / 1024.0)  * sin(phase + 2.0 * pi / 1024.0)) >= 0)) {phase += 2.0 * pi / 1024.0;}
 return short(floor(cos(phase)*127 + 0.5));
}
short NCO_gen_cos_phase(double freq, double phase_in, long long t)
{
	return short(floor(cos(double(t) * freq + phase_in)*127 + 0.5));
 double phase = (((long long)floor(double(t) * freq / 2.0 / pi * 16384 + 0.5) % 16384) >> 4) * 2.0 * pi / 1024.0+phase_in;
 if (((cos(phase) * sin(phase)) >= 0) & ((cos(phase + 2.0 * pi / 1024.0)  * sin(phase + 2.0 * pi / 1024.0)) >= 0)) {phase += 2.0 * pi / 1024.0;}
 return short(floor(cos(phase)*127 + 0.5));
}
short NCO_gen_sin_phase(double freq, double phase_in, long long t)
{
	return short(floor(sin(double(t) * freq + phase_in)*127 + 0.5));
 double phase = (((long long)floor(double(t) * freq / 2.0 / pi * 16384 + 0.5) % 16384) >> 4) * 2.0 * pi / 1024.0 + phase_in;
 if (((cos(phase) * sin(phase)) <= 0) & ((cos(phase + 2.0 * pi / 1024.0)  * sin(phase + 2.0 * pi / 1024.0)) <= 0)) {phase += 2.0 * pi / 1024.0;}
 return short(floor(sin(phase)*127 + 0.5));
}
// sin generator
short NCO_gen_sin(double freq, long long t)
{
    return short(floor(sin(double(t) * freq)*127 + 0.5));
 double phase = (((long long)floor(double(t) * freq / 2.0 / pi * 16384 + 0.5) % 16384) >> 4) * 2.0 * pi / 1024.0;
 if (((cos(phase) * sin(phase)) <= 0) & ((cos(phase + 2.0 * pi / 1024.0)  * sin(phase + 2.0 * pi / 1024.0)) <= 0)) {phase += 2.0 * pi / 1024.0;}
 return short(floor(sin(phase)*127 + 0.5));
}


//! Main program

int main(int argc, char *argv[]){//cout<<argc<<endl;
//printf("%s", argv[1]);
short z=0;
ofstream norm, str, viv, viv2;
ofstream test_out_0,test_out_1,test_out_2,test_out_4,test_out_5,test_out_6,test_out_7,test_out_8,test_out_9,test_out_dbg_0;
norm.open("data_out/norm.txt");//,ios::app
str.open("data_out/str.txt");
viv.open("data_out/viv.txt");
viv2.open("data_out/viv2.txt");

test_out_0.open("data_out/test_out_0.txt");
test_out_1.open("data_out/test_out_1.txt");
test_out_2.open("data_out/test_out_2.txt");
//test_out_3.open("data_out/test_out_3.txt");
test_out_4.open("data_out/test_out_4.txt");
test_out_5.open("data_out/test_out_5.txt");
test_out_6.open("data_out/test_out_6.txt");
test_out_7.open("data_out/test_out_7.txt");
test_out_8.open("data_out/test_out_8.txt");
test_out_9.open("data_out/test_out_9.txt");
test_out_dbg_0.open("data_out/test_out_dbg_0.txt");
//test_out_dbg_1.open("data_out/test_out_dbg_1.txt");
double fr0 = 0, fr1 = 0, fr;
cout<<"Enter fr: ";
//cin>>fr;
fr = 10;
fr0 = 13+fr*4;
fr1 = fr0+3.5;

short input_bw = short(fr);
short input_mult = 1<<(short(fr)-input_bw);
short lms_period = 64;
short filter_length = 64;
double JA = 0;
int parameter_0 = 4096;//32768;//256;
int parameter_1 = 131072;
//short agc_algorithm = 5;
short params_file_number = 0;
int* max_point = 0;
int max_point_value = 0;
int del_filt = 0;
short k_cordic = 2;
/*
for(short s = 0;s<64;s++)
{
    short f1 = s/2-s/8;
    short f2 = s/2-s/8;
    short f3_pre_1 = (s>>2)&0x1 && !(s&0x1);
    short f3_pre_2 = !((s>>2)&0x1) && s&0x1;
    short f3 = 0;
    if(f3_pre_1) f3 = s/4 - 1;
    else
        if(f3_pre_2) f3 = s/4 + 1;
        else f3 = s/4;
    cout<<s<<' '<<s-(f1+f2+f3)<<' '<<f1<<' '<<f2<<' '<<f3<<endl;
}
*/
/*
for(short s = 0;s<64;s++)
{
    short f1 = s*3/8;
    short f2 = s*3/8;
    //short f3_pre_1 = (s>>2)&0x1 && !(s&0x1);
    //short f3_pre_2 = !((s>>2)&0x1) && s&0x1;
    short f3 = s*2/8+((((((s*3)<<1)&0x3)+((s*2)&0x3))>>2)&0x3);
        cout<<s<<' '<<s-(f1+f2+f3)<<' '<<f1<<' '<<f2<<' '<<f3<<endl;
}
*/
for(short s = 0;s<128;s++)
{
    short f = ((s<<1) + (s<<2) + (s<<4))>>6;
    short err = ((s&0x3) + (f&0x3))&0x3;

    short f1 = f + ((err&0x1) || ((err>>1)&0x1));
    short f2 = f + ((err>>1)&0x1);
    short f3 = f + ((err&0x1) && ((err>>1)&0x1));
    //short f3_pre_1 = (s>>2)&0x1 && !(s&0x1);
    //short f3_pre_2 = !((s>>2)&0x1) && s&0x1;
    //short f3 = s*2/8+((((((s*3)<<1)&0x3)+((s*2)&0x3))>>2)&0x3);
    /*if(f3_pre_1) f3 = s/4 - 1;
    else
        if(f3_pre_2) f3 = s/4 + 1;
        else f3 = s/4;
        */
    cout<<s<<' '<<f<<' '<<' '<<f1<<' '<<f2<<' '<<f3<<' '<<s-(f1+f2+f3)<<' '<<err<<endl;
}
/*
for(short s = 0;s<64;s++)
{
    short f1 = s/4+s/16+((s>>3)&0x1);
    short f2 = s/4+s/16+(s&0x3);
    //short f3_pre_1 = (s>>2)&0x1 && !(s&0x1);
    //short f3_pre_2 = !((s>>2)&0x1) && s&0x1;
    short f3 = s/4+s/8+((s>>2)&0x1);
    //if(f3_pre_1) f3 = s/4 - 1;
   // else
     //   if(f3_pre_2) f3 = s/4 + 1;
        //else
            //f3 = s/4+s/8;
        cout<<s<<' '<<s-(f1+f2+f3)<<' '<<f1<<' '<<f2<<' '<<f3<<' '<<s%16<<endl;
}
*/
/*
for(short s = 0;s<64;s++)
{
    short f1 = s/4;
    short f2 = s/4+s/8;
    //short f3_pre_1 = (s>>2)&0x1 && !(s&0x1);
    //short f3_pre_2 = !((s>>2)&0x1) && s&0x1;
    short f3 = s/4+s/8;
    //if(f3_pre_1) f3 = s/4 - 1;
    // else
    //   if(f3_pre_2) f3 = s/4 + 1;
    //else
    //f3 = s/4+s/8;
    cout<<s<<' '<<s-(f1+f2+f3)<<' '<<f1<<' '<<f2<<' '<<f3<<' '<<s%16<<endl;
}
*/
double locata_pot_pre = 6;
short pos_filt_en = 0;
double k_cycle = 0;
short fil_pos_length = 128;
short fil_pos_number = 2;
//for(pos_filt_en = 1; pos_filt_en<=1; pos_filt_en++)
//for(params_file_number = 0; params_file_number<=6; params_file_number++)
//for(k_cycle = 0; k_cycle<=12; k_cycle = k_cycle +1)
    //for(locata_pot_pre = 1; locata_pot_pre<=10; locata_pot_pre++)
//for(agc_algorithm = 0; agc_algorithm<3; agc_algorithm++)
//for(del_filt = -8; del_filt<17; del_filt++)
{
    double locata_pot = locata_pot_pre*20;
    //if (parameter_0==65536) agc_algorithm = 0;
    //else agc_algorithm = 1;
    //int tty1=-799;
    //int tty2=799;
    //cout<<tty1/8<<' '<<tty2/8<<' '<<(tty1>>3)<<' '<<(tty2>>3)<<endl;
    short tf=0, j=0, k=0, l=0, l1=0, l2=0, filter_delay=0, kk=0;
    int input_signal=0, input_signal2=0, base_cosine=0, base_sine=0, gln_opc=0,gln_ops=0, inp_prn_ajm[6000], inp_prn_ajm2[6000];
    int complex_signal_i_L5=0, complex_signal_q_L5=0,base_cosine_L5=0, base_sine_L5=0;
    short aj_i_0=0,aj_i_1=0,aj_q_0=0,aj_q_1=0;
    short  sts=0, stc=0, ct=0, sign_cosine=1, sign_sine=1,outre=0, outim=0, max_coeff=0;
    short *outp1=NULL, *outp2=NULL, sign_cosine_semicycle=3, inp_prn[1023], inp_prn_glo[511], *e5a=NULL, *e5b=NULL;
    short f_sec_agc=0,f_sec_agc_tr=0,f_sec_agc_tr2=0,f_sec_agc_q=0,f_sec_agc_tr_q=0,f_sec_agc_tr2_q=0;
    int complex_signal_i=0, complex_signal_q=0;
    long long loe=0;
    int v=0, tt=0, outm=0, outm2=0;
#if adc_agc_mode==0
    short sampling_frequency=200;
#else
    short sampling_frequency=50;//43
#endif
#if LOCATA==1
    sampling_frequency=500;
#endif
    double sampling_interval=1/(sampling_frequency*1e6), sh1=1, sh2=1, sh=1;
    double *nmetod1=NULL, *nmetod2=NULL;
    char *f1=NULL, *f2=NULL;
    nmetod1=new double[24000];
    nmetod2=new double[24000];
    const short input_sh_delay=2000;
    outp1=new short[input_sh_delay];
    outp2=new short[input_sh_delay];
    e5a=new short[10000];
    e5b=new short[10000];
    for(short g=0;g<input_sh_delay;g++)
    {
        outp1[g]=0;
        outp2[g]=0;
    }
    //copy(&outp1[0], &outp1[64], ostream_iterator<short>(cout, " "));
    f1=new char[30];
    f2=new char[30];
    StNormDev gaussian_noise, code, L5_noise, locata_RF_noise, test0, test1;
    double signal_amplitude=0;//0.1;//0.0707;//0.0577;
    GALE5aIgen(e5a);
    GALE5bIgen(e5b);
LB1:
    l1=22;//cin>>l1;//l1=12
    l2=12;//cin>>l2;
    j=10;//cin>>j;
    long pulse_width = 0;
    cout<<"Filter length: 64 "<<endl;
    cout<<"Enter pulse width: ";
    //cin>>pulse_width;
    l = 21;
    short switch_filter_use=1;
    f2="data_in/f10_2_60.fcf";
    switch(l)
    {
        case 0: f1="data_in/f10_1_180_new.fcf"; f2="data_in/f10_2_60_new.fcf"; l1=26; filter_delay=0; switch_filter_use=0; break;
        case 1: f1="data_in/f10_01_150.fcf"; f2="data_in/f10_02_50.fcf"; l1=12; l2=23; filter_delay=27; switch_filter_use=1; break;
        case 2: f1="data_in/f20_1_150.fcf"; f2="data_in/f20_2_50.fcf"; l1=22; l2=22; filter_delay=29; switch_filter_use=1; break;
        case 3: f1="data_in/f45_06.fcf"; f2="data_in/f10_02_50.fcf"; l1=26; l2=23; filter_delay=13; switch_filter_use=1; break;
        case 4: f1="data_in/f45_04.fcf"; l1=26; filter_delay=13; break;
        case 5: f1="data_in/f45_05.fcf"; l1=21; filter_delay=10; break;
        case 6: f1="data_in/f20_1.fcf"; f2="data_in/f20_2.fcf"; l1=22; l2=22; filter_delay=29; switch_filter_use=1; break;
        case 7: f1="data_in/f10_1_180.fcf"; f2="data_in/f10_2_60.fcf"; l1=12; l2=23; filter_delay=27; switch_filter_use=1; break;
        case 8: f1="data_in/f45_08.fcf"; l1=11; filter_delay=5; break;
        case 9: f1="data_in/f45_09.fcf"; l1=11; filter_delay=5; break;
        case 10: f1="data_in/f45_10.fcf"; l1=4; filter_delay=2; break;
        case 11: f1="data_in/f45_11.fcf"; l1=2; filter_delay=1; break;
        case 12: f1="data_in/f45_12.fcf"; l1=4; filter_delay=2; break;
        case 13: f1="data_in/f12_1_180.fcf"; f2="data_in/f12_2_60.fcf"; l1=16; l2=23; filter_delay=28; switch_filter_use=1; break;
        case 14: f1="data_in/f10_1_180_new.fcf"; f2="data_in/f10_2_60_new.fcf"; l1=10; l2=18; filter_delay=21; switch_filter_use=1; break;
        case 15: f1="data_in/f20_1_180_new.fcf"; f2="data_in/f20_2_60_new.fcf"; l1=19; l2=18; filter_delay=24; switch_filter_use=1; break;
        case 16: f1="data_in/f12_1_180_new.fcf"; f2="data_in/f12_2_60_new.fcf"; l1=11; l2=20; filter_delay=24; switch_filter_use=1; break;
        case 20: f1="data_in/f45_12.fcf"; f2="data_in/f12_2_60_new.fcf"; l1=20; l2=20; filter_delay=20; switch_filter_use=2; break;
        //case 21: f1="data_in/equalizing_180_10_v03.fcf"; f2="data_in/equalizing_27_v.fcf";max_coeff=2;l1=20;filter_delay=23*1+0; switch_filter_use=2; break;
        //! For filt_div = 8 : filter_delay = 265
        //! For filt_div = 8 : filter_delay = 425 new real signal
        case 21:    f1="data_in/equalizing_180_10_v04_80.fcf";
                    f2="data_in/equalizing_27_v.fcf";
                    max_coeff=2;l1=20;
                    #if ajm_real==1
                    filter_delay=420 + 5-8 +8 -4 + 8 -232 -0 +64+16       +fil_pos_length*4*pos_filt_en*fil_pos_number;//-80
                    #else
                    filter_delay=265;//-80
                    #endif
                    //filter_delay=0;//-4
                    switch_filter_use=2;
                    break;
        //case 21: f1="data_in/equalizing_180_20_v01.fcf"; f2="data_in/equalizing_27_v.fcf";max_coeff=2;l1=20;filter_delay=22*1+0; switch_filter_use=2; break;
        //case 21: f1="data_in/equalizing_180_10_v05_97.fcf"; f2="data_in/equalizing_27_v.fcf";max_coeff=2;l1=20;filter_delay=25*1+0; switch_filter_use=2; break;
        default: cout<<" this filter does not exist!"<<endl; goto LB1;												//22+64	//8	//256
    }
    double jam_amplitude,jam_frequency_offset;
    short jam_on;
    short cic_del_m[6] = {0,0,0,0,0,0};
    short nco_freq = 0; //7281
    short filt_cic_out_shift,mid_agc_lev_n,mid_agc_lev_p,mid_agc_count_n,mid_agc_count_p,mid_agc_man,lms_mju;
    short h1[30];
    short h2[10];
    short *h_pos;
    h_pos = new short[fil_pos_length];
    short rf_real_0_en,rf_I_en,rf_Q_en,comb_out_en,d6_en,co_rnd_en,cic_out_en,sfir_en,aj_out_en,qnt_en;
    long log_length;

    double jamming_frequency=0,jamming_frequency_2=0,jamming_frequency_3=0;
    double jamming_frequency_4=0,jamming_frequency_5=0,jamming_frequency_6=0;
    double jamming_frequency_7=0,jamming_frequency_8=0,jamming_frequency_9=0,jamming_frequency_10=0;
    double jamming_amplitude=0,jamming_amplitude_2=0,jamming_amplitude_3=0;
    double jamming_amplitude_4=0,jamming_amplitude_5=0,jamming_amplitude_6=0;
    double jamming_amplitude_7=0,jamming_amplitude_8=0,jamming_amplitude_9=0,jamming_amplitude_10=0;

    double nco_carr_frec;
    int jam_pulse_on;
    int jam_pulse_width;
    int jam_pulse_period;
    int jam_fm_on;
    double jam_fm_freq;
    double jam_fm_speed;
    int jam_am_on;
    double jam_am_freq;
    double jam_am_depth;
    int filt_div;
    //! Reading parameters from configuration file
    char shift1[64];
    FILE *config;
    if(argv[1] == NULL)
        switch(params_file_number)
        {
            #if LOCATA==0
            case 0:    config = fopen( "fparams.txt", "r+" ); break;
            #else
            case 0:    config = fopen( "fparams_locata.txt", "r+" ); break;
            #endif
            case 1:    config = fopen( "fparams02.txt", "r+" ); break;
            case 2:    config = fopen( "fparams03.txt", "r+" ); break;
            case 3:    config = fopen( "fparams04.txt", "r+" ); break;
            case 4:    config = fopen( "fparams05.txt", "r+" ); break;
            case 5:    config = fopen( "fparams06.txt", "r+" ); break;
            case 6:    config = fopen( "fparams07.txt", "r+" ); break;
            case 7:    config = fopen( "fparams08.txt", "r+" ); break;
            case 8:    config = fopen( "fparams09.txt", "r+" ); break;
            case 9:    config = fopen( "fparams10.txt", "r+" ); break;
            case 10:   config = fopen( "fparams11.txt", "r+" ); break;
            case 11:   config = fopen( "fparams12.txt", "r+" ); break;
            case 12:   config = fopen( "fparams13.txt", "r+" ); break;
            case 13:   config = fopen( "fparams14.txt", "r+" ); break;
            case 14:   config = fopen( "fparams15.txt", "r+" ); break;
            case 15:   config = fopen( "fparams16.txt", "r+" ); break;
            case 16:   config = fopen( "fparams17.txt", "r+" ); break;
            case 17:   config = fopen( "fparams18.txt", "r+" ); break;
            case 18:   config = fopen( "fparams19.txt", "r+" ); break;
            default:   config = fopen( "fparams.txt", "r+" ); break;
        }
    else config = fopen( argv[1], "r+" );
    //config = fopen( "fscanf.txt", "r+" );
    #if ajm_real==0
    config = fopen( "fparamsc.txt", "r+" );
    #endif
    if( config == NULL )
        printf( "The file fparams.txt was not opened\n" );
    else
    {
        fseek( config, 0L, SEEK_SET );
        for(short i=0;i<3;i++)fscanf( config, "%s", &shift1 );
        fscanf( config, "%s", &shift1 );fscanf( config, "%lf", &signal_amplitude );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &nco_carr_frec );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency_2 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency_3 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency_4 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency_5 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency_6 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency_7 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency_8 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency_9 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_frequency_10 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude_2 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude_3 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude_4 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude_5 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude_6 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude_7 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude_8 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude_9 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jamming_amplitude_10 );
        fscanf (config, "%s", &shift1 );fscanf( config, "%ld", &jam_pulse_on );
        fscanf (config, "%s", &shift1 );fscanf( config, "%ld", &jam_pulse_width );
        fscanf (config, "%s", &shift1 );fscanf( config, "%ld", &jam_pulse_period );
        fscanf (config, "%s", &shift1 );fscanf( config, "%ld", &jam_fm_on );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jam_fm_freq );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jam_fm_speed );
        fscanf (config, "%s", &shift1 );fscanf( config, "%ld", &jam_am_on );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jam_am_freq );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jam_am_depth );

        for(short i=0;i<2;i++)fscanf( config, "%s", &shift1 );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &nco_freq );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &filt_div );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &cic_del_m[0] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &cic_del_m[1] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &cic_del_m[2] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &cic_del_m[3] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &cic_del_m[4] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &cic_del_m[5] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &filt_cic_out_shift );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &mid_agc_lev_n );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &mid_agc_lev_p );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &mid_agc_count_n );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &mid_agc_count_p );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &mid_agc_man );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &lms_mju );
        for(short i=0;i<30;i++)
        {
            fscanf( config, "%s", &shift1 );
            fscanf( config, "%ld", &h1[i] );
        }
        for(short i=0;i<10;i++)
        {
            fscanf( config, "%s", &shift1 );
            fscanf( config, "%ld", &h2[i] );
        }
        for(short i=0;i<fil_pos_length;i++)
        {
            fscanf( config, "%s", &shift1 );
            fscanf( config, "%ld", &h_pos[i] );
        }

        for(short i=0;i<22;i++)fscanf( config, "%s", &shift1 );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &log_length );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &rf_real_0_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &rf_I_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &rf_Q_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &comb_out_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &d6_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &co_rnd_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &cic_out_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &sfir_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &aj_out_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &qnt_en );

        fclose( config );
    }
    //jamming_amplitude = jamming_amplitude * pow(1.7783,k_cycle);
    //! Maiking classes for filters and filters initialization
    Filter Flt,Flt_Q;//,Flt_L5,Flt_L5_Q;
    Filter FltPos,FltPos2,FltPos_Q,FltPos2_Q;
    AGC inp_agc;
    //L5_amp_ctrl gain_L5;
    double *hsat,*hsat2;
    short lsat=90,lsat2=90;
    hsat=new double[lsat];
    hsat2=new double[lsat2];

    //lms_mju = 5;
    cout<<"................."<<endl;
    cout<<5/2<<' '<<-5/2<<endl;
    //f1="data_in/equalizing_180_20_v01.fcf";
    Flt.LoadCoeff(switch_filter_use,cic_del_m,filt_cic_out_shift,mid_agc_lev_n,mid_agc_lev_p,
                  mid_agc_count_n,mid_agc_count_p,mid_agc_man,lms_mju,h1,h2,h_pos);
    Flt_Q.LoadCoeff(switch_filter_use,cic_del_m,filt_cic_out_shift,mid_agc_lev_n,mid_agc_lev_p,
                  mid_agc_count_n,mid_agc_count_p,mid_agc_man,lms_mju,h1,h2,h_pos);
    FltPos.LoadCoeff(switch_filter_use,cic_del_m,filt_cic_out_shift,mid_agc_lev_n,mid_agc_lev_p,
                  mid_agc_count_n,mid_agc_count_p,mid_agc_man,lms_mju,h1,h2,h_pos);
    FltPos2.LoadCoeff(switch_filter_use,cic_del_m,filt_cic_out_shift,mid_agc_lev_n,mid_agc_lev_p,
                  mid_agc_count_n,mid_agc_count_p,mid_agc_man,lms_mju,h1,h2,h_pos);
    /*FltPos_Q.LoadCoeff(switch_filter_use,cic_del_m,filt_cic_out_shift,mid_agc_lev_n,mid_agc_lev_p,
                  mid_agc_count_n,mid_agc_count_p,mid_agc_man,lms_mju,h3,h2);
    FltPos2_Q.LoadCoeff(switch_filter_use,cic_del_m,filt_cic_out_shift,mid_agc_lev_n,mid_agc_lev_p,
                  mid_agc_count_n,mid_agc_count_p,mid_agc_man,lms_mju,h3,h2);*/
    f1="data_in/equalizing_180_aj_v01.fcf";//f1="data_in/equalizing_180_12.5_v02.fcf";
    max_coeff=1;
    UnifDev jaja(1,80), locata_sh(-1,1), locata_amp_jump(8,16);
    int point[80];
    for(int i=1; i<=80; i++)
    {
        LB_POINT:
        point[i-1]=int(jaja.RandValue());
        for(int j=2; j<i; j++) if(point[i-1]==point[j-1]) goto LB_POINT;
    }
    InputPosled(inp_prn);
    Flt.LoadAJFilter(80,255,short(6+4+0), 10);
    lms_period = 5;
    Flt_Q.LoadAJFilter(80,255,short(6+4+0), 10);

#if adc_agc_mode==0
    double carrier_frequency_pre=nco_carr_frec+0.0, base_frequency=(nco_carr_frec + 0.0)*1e6*2*pi*sampling_interval;
#else
    double carrier_frequency_pre=16, base_frequency=16*1e6*2*pi*sampling_interval;
#endif
    double carrier_frequency=carrier_frequency_pre*1e6*2*pi*sampling_interval;

    cout << "carrier_frequency_pre: " << carrier_frequency_pre <<endl;
    cout << "sampling_interval: " << sampling_interval <<endl;
    cout << "jamming_frequency: " << jamming_frequency<<' '<<jamming_amplitude<<endl;

    cout << "lms_mju: " << lms_mju <<endl;

    jamming_frequency=(carrier_frequency_pre+jamming_frequency - 0.0)*1e6*2*pi*sampling_interval;
    jamming_frequency_2=(carrier_frequency_pre+jamming_frequency_2)*1e6*2*pi*sampling_interval;
    jamming_frequency_3=(carrier_frequency_pre+jamming_frequency_3)*1e6*2*pi*sampling_interval;
    jamming_frequency_4=(carrier_frequency_pre+jamming_frequency_4)*1e6*2*pi*sampling_interval;
    jamming_frequency_5=(carrier_frequency_pre+jamming_frequency_5)*1e6*2*pi*sampling_interval;
    jamming_frequency_6=(carrier_frequency_pre+jamming_frequency_6)*1e6*2*pi*sampling_interval;
    jamming_frequency_7=(carrier_frequency_pre+jamming_frequency_7)*1e6*2*pi*sampling_interval;
    jamming_frequency_8=(carrier_frequency_pre+jamming_frequency_8)*1e6*2*pi*sampling_interval;
    jamming_frequency_9=(carrier_frequency_pre+jamming_frequency_9)*1e6*2*pi*sampling_interval;
    jamming_frequency_10=(carrier_frequency_pre+jamming_frequency_10)*1e6*2*pi*sampling_interval;

    double w=124.58;

    double pomeha=0;

    tt=0;outm=0;outm2=0;sts=0;stc=0;kk=0;v=0;ct=0;sign_cosine=1;sign_sine=1;
    sign_cosine_semicycle=sampling_frequency/60;
    for(short g=0;g<input_sh_delay;g++)
    {
        outp1[g]=0;
        outp2[g]=0;
    }
    const short correlate_offset=50;
    //loe=(((100+correlate_offset)*1000)*sampling_frequency+1);//loe=50000;//21600000001
    loe=log_length;

    sh1=double(e5a[kk]+e5b[kk]);
    sh2=double(e5a[kk]-e5b[kk]);
    sh=1;
    double d_input_signal=0, d_sat_signal=0, L5_amp=0;
    short jamming_situation=0;
    if (jam_pulse_on==1) jamming_situation=5;

    cout << "signal_amplitude: " << signal_amplitude <<endl;
    cout << "nco_carr_frec: " << nco_carr_frec <<endl;
    cout << "jamming_frequency: " << jamming_frequency <<endl;


    cout << "jam_pulse_on: " << jam_pulse_on << endl;
    cout << "--jamming_situation: " << jamming_situation << endl;

    cout << "jam_am_on: " << jam_am_on << endl;
    cout << "jam_fm_on: " << jam_fm_on << endl;

    int outren=0,outimn=0,outren_d=0,outimn_d=0,outren_L5=0,outimn_L5=0;
    tt=0;outm=0;outm2=0;sts=0;stc=0;kk=0;v=0;ct=0;sign_cosine=1;sign_sine=1;
    sign_cosine_semicycle=sampling_frequency/60;
    int adc_of_bit=0, agc_mid_level=0, adc_high_bit=0, adc_low_bit=0;
    int adc_of_bitp=0, adc_high_bitp=0, adc_low_bitp=0;
    int adc_of_bitn=0, adc_high_bitn=0, adc_low_bitn=0;
    long ref_mid_level=320*512; //vk
    k_cordic = 2;
    int jam_detector=0;
    long agc_mid_level_av=0;
    short adc_of_bit_av=0;
    short adc_high_bit_av=0;
    short adc_low_bit_av=0;
    short adc_of_bit_av_temp=0;
    short adc_low_bit_av_temp=0;
    short adc_high_bit_av_temp=0;

    short input_signal_tr_temp[1024];
    short input_signal_tr_temp2[1024];
    short input_signal_tr_temp3[1024];
    short input_signal_tr_temp4[1024];
    for(short i=0; i<1024; i++)
    {
        input_signal_tr_temp[i]=0;
        input_signal_tr_temp2[i]=0;
        input_signal_tr_temp3[i]=0;
        input_signal_tr_temp4[i]=0;
    }

    int input_signal_mem[128];
    short input_signal_tr_1[128];
    short input_signal_tr_2[128];
    short input_signal_tr_3[128];
    short input_signal_tr_4[128];
    for(int i=0; i<128; i++)
    {
        input_signal_mem[i]=0;
        input_signal_tr_1[i]=0;
        input_signal_tr_2[i]=0;
        input_signal_tr_3[i]=0;
        input_signal_tr_4[i]=0;
    }

    int agc_out=0;
    int sum_out_fir=0, sum_out_fir_L5=0, L5_agc_out=0;
    double real_gain_coeff=1,pomeha_FM=0,pomeha_AM=0;
    StNormDev imp_jam_amp; //0,20,40,...,600
    UnifDev imp_jam_dur(-1,1); //20mks-200mks
    StUnifDev4::SetKernel(35189);
    int imp_jam_d=int(imp_jam_dur.RandValue()), imp_jam_count=0;
    double imp_jam_a=floor(abs(imp_jam_amp.RandValue()))*20, gaussian_noise_2=0;
    //int pulse_period=100;
    //int pulse width=20;
    int spectr_an_outre=0,spectr_an_outim=0,sa_outren=0,sa_outimn=0;
    short spectran_cnt1 = 0;
    short spectran_cnt2 = 0;
    int spectran_sum[512], spectran_sum_av = 0, sa_y = 0, sa_y1 = 4095, sa_max = 0, sa_min = 0;
    for(short g=0;g<512;g++) spectran_sum[g] = 0;

    short complex_signal_i_i=32;
    short k_agc_deay[3]={0,0,0};
    short sh_fm = 0;
    short sinc_inp_i[18],sinc_inp_q[18];
    short sinc_inp_si=0,sinc_inp_sq=0;
    for(short g=0;g<18;g++)
    {
        sinc_inp_i[g]=0;
        sinc_inp_q[g]=0;
    }
    short aj_in_on=0;
    short sig_en = 1;
    double phase_in = 0;
    double jamm_noise_comp = 1;
    double jnc_discr = 0.5;
    short input_max = 1<<(input_bw-1);
    int sec_agc_sum = 1, sec_agc_sum_pre = 1;
    short f_sec_agc_old = 0;
    f_sec_agc=0,f_sec_agc_tr=0,f_sec_agc_tr2=0,f_sec_agc_q=0,f_sec_agc_tr_q=0,f_sec_agc_tr2_q=0;
    //! Delays for bit true correspondence to verilog
    int delay_line_inp[64];
    int delay_line_cos[64];
    int delay_line_sin[64];
    int delay_line_sign_i[64];
    int delay_line_sign_q[64];
    int delay_line_out_1[64];
    int delay_line_out_2[64];
    int delay_line_out_4[64];
    int delay_line_out_4q[64];
    int delay_line_out_6[64];
    int delay_line_out_7[64];
    int delay_line_out_6q[64];
    int delay_line_out_7q[64];
    int delay_line_out_8[64];
    int delay_line_out_9[64];
    int delay_line_out_8q[64];
    int delay_line_out_9q[64];
    int delay_line_cics[64];
    for(short i=0;i<64;i++)
    {
        delay_line_inp [i]=0;
        delay_line_cos [i]=2047;
        delay_line_sin [i]=0;
        delay_line_sign_i [i]=0;
        delay_line_sign_q [i]=0;
        delay_line_out_1 [i]=0;
        delay_line_out_2 [i]=0;
        delay_line_out_4 [i]=0;
        delay_line_out_4q [i]=0;
        delay_line_out_6 [i]=0;
        delay_line_out_7 [i]=0;
        delay_line_out_6q [i]=0;
        delay_line_out_7q [i]=0;
        delay_line_out_8 [i]=0;
        delay_line_out_9 [i]=0;
        delay_line_out_8q [i]=0;
        delay_line_out_9q [i]=0;
        delay_line_cics[i] = 0;
    }

    int cos_tmp = 0;
    int sin_tmp = 0;
    int inp_tmp = 0;

    int d_x = 0;
    double int_d_x = 0;
    double int_d_x2 = 0;
    int int_d_x_pre = 0;
    int outren_pre = 0;
    int outren_pre2 = 0;
    int outren_pre3 = 0;
    int outren_pre4 = 0;
    int slow_mode = 0;
    int slow_mode_cnt = 0;
    int t_freq_move = 0;
    int t_freq_move_fast = 0;

    double PLL_BW = 3.5e5;
    double Alpha = 2*PLL_BW;
    double Beta = Alpha*PLL_BW*double(sampling_interval)/0.6*8;
    double Gamma = PLL_BW/2.4*double(sampling_interval)*8;
    double Discr = 0;
    double PLLOut = 0;
    double Sum1 = 0;
    double Sum2 = 0;

    Flt.NCO(nco_freq,0);
    input_signal=0;
    short kkk = 0;
    short kkk2 = 1022*0;
    short new_epoch = 0;
    short new_epoch2 = 0;
    short new_epoch_pre = 0;
    short new_chip = 0;
    short new_chip2 = 0;
    int outren_tmp_agc = 0;
    int code_shift = 8600000*200*0;
    long long code_shift_all = 0;
    short accum_period = 0;
    short accum_period2 = 0;
    short accumulation_time_ms = 1;
    #if LOCATA==1
        //double locata_noise_sigma = sqrt(1.0/(pow(10.0,locata_pot/10.0)*2/50e6));
        //locata_pot = floor(locata_amp_jump.RandValue()+0.5)*20;
        double locata_noise_sigma = 0.05/128.;
        double locata_amp = sqrt(locata_noise_sigma*locata_noise_sigma*pow(10.0,locata_pot/10.0)*2/50e6);
        double agc_k_locata = 0;
        double locata_in = 0;
        short locata_in_dig = 0;
    #endif
    double lfm = 1;
    short *avt2_data = NULL;
    avt2_data = new short[65536];
    InputAVT2data(avt2_data);
    short gain_mode = 0;
    short jam_detect_cnt = 255;
    short jam_detect_number = 0;
    int aj_coeff_sum_0 = 0;
    int aj_coeff_sum_1 = 0;
    short abs_outren_0 = 0;
    short abs_outren_pre_0 = 0;
    short abs_outren_1 = 0;
    short abs_outren_pre_1 = 0;
    short w_aj_abs_max_0 = 0;
    short w_aj_abs_max_1 = 0;

    short jam_osn_thr_int_0 = 0;
    short jam_osn_thr_int_1 = 0;
    int jam_osn_thr_sum_0 = 0;
    int jam_osn_thr_sum_1 = 0;
    int jam_osn_thr_mulsum_0 = 0;
    int jam_osn_thr_mulsum_1 = 0;
    short jam_first_thr_0 = 0;
    short jam_first_thr_1 = 0;
    short jam_exist_0 = 0;
    short jam_exist_1 = 0;
    //for(long cort = 1; cort<=2048*1024; cort++)
        //norm<<inp_agc.Cordic_log2(555)<<endl;
    //cin>>new_chip;
    short Nfft = 64;
    int a_fft[64],b_fft[64];
     //Initialising
    /*a_fft[0]  =    -89;
    a_fft[1]  =     22;
    a_fft[2]  =    226;
    a_fft[3]  =    -58;
    a_fft[4]  =    252;
    a_fft[5]  =     95;
    a_fft[6]  =    235;
    a_fft[7]  =     25;
    a_fft[8]  =   -131;
    a_fft[9]  =   -296;
    a_fft[10] =     31;
    a_fft[11] =    164;
    a_fft[12] =    -59;
    a_fft[13] =   -108;
    a_fft[14] =    -62;
    a_fft[15] =   -219;
    a_fft[16] =    -99;
    a_fft[17] =    -36;
    a_fft[18] =      9;
    a_fft[19] =    -13;
    a_fft[20] =    122;
    a_fft[21] =     22;
    a_fft[22] =    363;
    a_fft[23] =     62;
    a_fft[24] =    361;
    a_fft[25] =   -145;
    a_fft[26] =    105;
    a_fft[27] =    -52;
    a_fft[28] =    120;
    a_fft[29] =    119;
    a_fft[30] =   -437;
    a_fft[31] =   -265;
    a_fft[32] =   -288;
    a_fft[33] =     80;
    a_fft[34] =    294;
    a_fft[35] =    -65;
    a_fft[36] =    162;
    a_fft[37] =    109;
    a_fft[38] =   -210;
    a_fft[39] =     79;
    a_fft[40] =   -150;
    a_fft[41] =    303;
    a_fft[42] =     -7;
    a_fft[43] =    327;
    a_fft[44] =    -85;
    a_fft[45] =    118;
    a_fft[46] =    -13;
    a_fft[47] =   -404;
    a_fft[48] =   -196;
    a_fft[49] =    123;
    a_fft[50] =    -11;
    a_fft[51] =   -224;
    a_fft[52] =   -125;
    a_fft[53] =     50;
    a_fft[54] =   -199;
    a_fft[55] =    195;
    a_fft[56] =   -128;
    a_fft[57] =    362;
    a_fft[58] =   -216;
    a_fft[59] =     40;
    a_fft[60] =   -304;
    a_fft[61] =   -145;
    a_fft[62] =   -119;
    a_fft[63] =     80;
*/
    for(short i = 0; i < Nfft; i++) a_fft[i] = 0;
    for(short i = 0; i < Nfft; i++) b_fft[i] = 0;
    int *fft_data_in = NULL;
    int *fft_data0 = NULL;
    int *fft_data1 = NULL;
    short *fft_filt = NULL;
    fft_data_in = new int[Nfft];
    fft_data0 = new int[Nfft];
    fft_data1 = new int[Nfft];
    fft_filt = new short[Nfft];
    for(short i = 0; i < Nfft; i++) fft_data0[i] = 0;
    for(short i = 0; i < Nfft; i++) fft_data1[i] = 0;
    /*Flt.FFT(a_fft,b_fft,Nfft);
    for(short i=0;i<Nfft;i++)
    {
        //cout<<Flt.GetFFT(i,0)<<' '<<Flt.GetFFT(i,1)<<endl;
        viv2<<a_fft[i]<<' '<<b_fft[i]<<endl;
        norm<<Flt.GetFFT(i,0)<<' '<<Flt.GetFFT(i,1)<<endl;
        fft_data0[i] = Flt.GetFFT(i,0);
        fft_data1[i] = Flt.GetFFT(i,1);
        a_fft[i] = Flt.GetFFT(i,0);
        b_fft[i] = Flt.GetFFT(i,1);
        if(i==15  || i==49) fft_filt[i] = 0;
        else fft_filt[i] = 1;
    }
    for(short i=0;i<Nfft;i++)
    {
        fft_data0[i] = fft_data0[i] * fft_filt[i];
        fft_data1[i] = fft_data1[i] * fft_filt[i];
        str<<fft_data0[i]<<' '<<fft_data1[i]<<endl;
        a_fft[i] = fft_data0[i];
        b_fft[i] = fft_data1[i];
    }
    Flt.FFT(b_fft,a_fft,Nfft);
    for(short i=0;i<Nfft;i++)
    {
        //cout<<Flt.GetFFT(i,0)<<' '<<Flt.GetFFT(i,1)<<endl;
        viv<<Flt.GetFFT(i,0)<<' '<<Flt.GetFFT(i,1)<<endl;
        fft_data0[i] = Flt.GetFFT(i,0);
        fft_data1[i] = Flt.GetFFT(i,1);
        a_fft[i] = Flt.GetFFT(i,0);
        b_fft[i] = Flt.GetFFT(i,1);
        if(i>15 && i==16) fft_filt[i] = 0;
        else fft_filt[i] = 1;
    }
    cin>>new_chip;*/
	int *fft_in;
	fft_in = new int[512];
    double *fft_out;
	fft_out = new double[512];
    //! Main cycle start
    cout<<loe<<endl;
    for(long long t=0;t<loe;t++) //116000000
    {//cout<<t<<endl;
    #if AVT2==1
    if(t==65536) break;
    #endif
    if(ca_p_code==1) new_chip = Flt.NCO_code_gen_phase_gate(21968758,0,short(sh),0);//87875031  200 - 21968758; 50 - 87875031 150 - 29291676
        else if(ca_p_code==10) new_chip = Flt.NCO_code_gen(219687577);
        kkk+=new_chip;
        if(kkk==(1023*ca_p_code))
        //if(t%50000==0)
        {
            new_epoch++;
            accum_period++;
            if(accum_period==accumulation_time_ms)
            {
                if(v>(correlate_offset-1))
                {

                    nmetod1[v-correlate_offset]=double(outm);
                }
                str<<double(outm)<<' '<<double(outm2)<<' '<<filter_delay<<' '<<jamming_frequency/2/pi*200e6<<' '<<new_epoch<<endl;
                cout<<double(outm)<<' '<<double(outm2)<<' '<<filter_delay<<' '<<jamming_frequency/2/pi*200e6<<' '<<new_epoch<<endl;
                /*if(outm2<0)
                {
                    cout<<"code_shift_all = "<<code_shift_all<<endl;
                    goto LB3;
                }*/
                outm = 0;
                outm2 = 0;
                accum_period = 0;
                v++;
            }
            kkk = 0;
        }
        //else new_epoch = 0;
        if(ca_p_code==1) sh = double(inp_prn[kkk]);
        else if(ca_p_code==10)
            if(new_chip)
            {
                sh=code.RandValue();
                if(sh>0)sh=1;
                else sh=-1;
            }
        //!LOCATA code generator
        //if(t==1)code_shift = 8600000*200*1;
        //if(t==2)code_shift = 8600000*70*1;
        //if(t==3)code_shift = 500000*365*1+50000*12*1+5000*22*1+500*58*1;

        code_shift_all += code_shift;
        //if(ca_p_code==1) new_chip2 = Flt_Q.NCO_code_gen_phase_gate(87875031,1,short(sh2),code_shift);//87875031
            //else if(ca_p_code==10) new_chip2 = Flt_Q.NCO_code_gen(219687577);
            kkk2+=new_chip2;
        code_shift = 0;
        if(kkk2==(1023*ca_p_code))
        //if(t%50000==0)
        {
            new_epoch2++;
            accum_period2++;
            kkk2 = 0;
            if(accum_period2==accumulation_time_ms)
            {
                code_shift = 8600000*1;//8600000;//500000 5*0
                accum_period2 = 0;
            }
        }
        //else new_epoch = 0;
        if(ca_p_code==1) sh2 = double(inp_prn[kkk2]);
        else if(ca_p_code==10)
            if(new_chip2)
            {
                sh2=code.RandValue();
                if(sh2>0)sh2=1;
                else sh2=-1;
            }
        #if LOCATA==0
        //! Jamming situation setting
        switch(jamming_situation)
        {
            case 0: jamming_frequency=jamming_frequency + 0*4*0.25*0.125*1e0*2*pi*sampling_interval;
                    //if(lfm==0) lfm = -1;
                    //if(jamming_frequency>(150e6*2*pi*sampling_interval) && lfm==1) lfm = 0;
                    /*if(jamming_frequency>=(150e6*2*pi*sampling_interval) && lfm==1)
                    {
                        //jamming_frequency = 150e6*2*pi*sampling_interval;
                        lfm = -1;
                        //jamming_frequency += 4*0.25*0.125*1e0*2*pi*sampling_interval*t;
                        //t = -1;
                    }*/
                    //if(jamming_frequency>(143.58e6*2*pi*sampling_interval)) lfm = 1;
                        break;
            case 1: if (t>200000) {jamming_amplitude=50;}if (t>400000) jamming_amplitude=0;
                    if (t>600000) jamming_amplitude=50;
                    if (t>800000) jamming_frequency=141e6*2*pi*sampling_interval;
                    if (t>1000000) jamming_amplitude=0; break;
            case 2: if (t<00000) jamming_amplitude=0; else jamming_amplitude=50; break;//sig_en = 0; else sig_en = 1; break;
            case 3: jamming_amplitude=0; if (t>1000000 && t<1200000) jamming_amplitude=138;if (t>1400000 && t<1600000) jamming_amplitude=15; break;
            case 4: if (t>1000000) jamming_amplitude=double(138*sin(t*1e4*2*pi*sampling_interval));if (t>2800000) jamming_amplitude=0; break;
            case 5: if ((t - pulse_width)%jam_pulse_width==0 && t>0) jamming_amplitude=0;
                    if(t%jam_pulse_period==0  && t>0) jamming_amplitude=100; break;
            case 6: if (t%150000==0 && t>0) jamming_amplitude*=1.12201845430196; break;
            case 7: if(imp_jam_count>=imp_jam_d)
                    {
                        imp_jam_d=int(imp_jam_dur.RandValue());
                        imp_jam_a=floor(abs(imp_jam_amp.RandValue()))*20;
                        imp_jam_count=0;
                    }
                    imp_jam_count++;
                    jamming_amplitude=imp_jam_a; break;
            default: jamming_amplitude=0; break;
        }
        //pomeha_FM = jamming_amplitude*cos(jamming_frequency*double(t)+jam_fm_freq*2*pi*sampling_interval*sin(jam_fm_speed*2*pi*sampling_interval*t));
        pomeha_FM = jamming_amplitude*cos(jamming_frequency*double(t)+jam_fm_freq*sin(jam_fm_speed*2*pi*sampling_interval*t));
        pomeha_AM = jamming_amplitude*(1+jam_am_depth*cos(jam_am_freq*1e3*2*pi*sampling_interval*t))*cos(double(t)*jamming_frequency);

        //! Generate input radio signal
        d_input_signal = gaussian_noise.RandValue()*1+1*sh*signal_amplitude*cos(double(t)*carrier_frequency);
        d_input_signal += jam_am_on*pomeha_AM+jam_fm_on*pomeha_FM+(jam_am_on-1)*(jam_fm_on-1)*jamming_amplitude*cos(double(t)*jamming_frequency);
        d_input_signal += jamming_amplitude_2*cos(double(t)*jamming_frequency_2);
        d_input_signal += jamming_amplitude_3*cos(double(t)*jamming_frequency_3);
		d_input_signal += jamming_amplitude_4*cos(double(t)*jamming_frequency_4);
        d_input_signal += jamming_amplitude_5*cos(double(t)*jamming_frequency_5);
		d_input_signal += jamming_amplitude_6*cos(double(t)*jamming_frequency_6);
        //d_input_signal += jamming_amplitude*cos(double(t)*2*pi*(jamming_frequency-11e6)*sampling_interval+pi/8);
        //d_input_signal += jamming_amplitude*cos(double(t)*2*pi*(jamming_frequency-3*1e6)*sampling_interval+pi/3);
        //d_input_signal += jamming_amplitude*cos(double(t)*2*pi*(jamming_frequency+9*1e6)*sampling_interval+pi/2);
        //d_input_signal += jam_am_on*pomeha_AM+jam_fm_on*pomeha_FM;
        //d_input_signal += (jam_am_on-1)*(jam_fm_on-1)*jamming_amplitude*cos(double(t)*jamming_frequency+25e6*2*pi*sampling_interval*sampling_interval*t*t);
        //viv2<<d_input_signal<<endl;
        if(rf_real_0_en==1) test_out_0<<input_signal<<endl;
        //if(new_epoch>=30 && new_epoch<=33) str<<d_input_signal<<' ';
        //real_gain_coeff = 1;
        d_input_signal *= real_gain_coeff;
        //viv<<d_input_signal<<endl;
        if(new_epoch>=575 && new_epoch<=579) norm<<d_input_signal<<endl;
        if(new_epoch>=575 && new_epoch<=579) viv2<<double(t)*jamming_frequency<<' '<<t<<' '<<jamming_frequency<<endl;
        //! ADC quantization of input radio signal
#if adc_agc_mode==0
        input_signal=short(floor(d_input_signal*(input_max-1)+0.5));
        if(input_signal>=(input_max-1))
        {
            input_signal = (input_max-1);
        }
        else if(input_signal<=(-input_max))
        {
            input_signal = -input_max;
        }
        #if AVT2==1
        input_signal = avt2_data[t];
        #endif
        //input_signal += 6;
        //viv<<input_signal<<' '<<input_signal/128<<endl;
        //input_signal = (input_signal & 0xfffffff0);
        //int input_signal_test = (input_signal & 0xfffffff0);
        short input_signal_sign = 0;
        if(input_signal<0) input_signal_sign = 1;
        input_signal = (abs(input_signal) & 0xffffffff);
        if(input_signal_sign==1) input_signal = -input_signal;
        //norm<<input_signal<<' '<<input_signal_test<<endl;
        agc_mid_level+=abs(input_signal<<2);
        //viv2<<input_signal<<endl;//<<' '<<real_gain_coeff<<' '<<agc_mid_level<<endl;
        input_signal*=input_mult;
        /*
        for(short i = 63; i>0; i--) delay_line_inp[i] = delay_line_inp[i-1];
        delay_line_inp[0] = input_signal;
        inp_tmp = delay_line_inp[2];
        */
        //str<<input_signal<<endl;
    #else
        /*input_signal=short(floor(d_input_signal*(128-1)+0.5));
        if(input_signal>=(128-1))
        {
            input_signal = (128-1);
        }
        else if(input_signal<=(-128))
        {
            input_signal = -128;
        }
        */
        input_signal = 0;
                if(d_input_signal>0)
            {
                if(d_input_signal>0.25) input_signal = 1;
                if(d_input_signal>0.5) input_signal = 2;
                if(d_input_signal>1) input_signal = 3;
            } else
            {
                if(d_input_signal<-0.25) input_signal = -1;
                if(d_input_signal<-0.5) input_signal = -2;
                if(d_input_signal<-1) input_signal = -3;
            }
        //norm<<d_input_signal<<' '<<input_signal<<endl;
        /*
        agc_mid_level=10000;
        if(d_input_signal>=0 && d_input_signal<0.000000002) {input_signal=1;}//viv<<0<<endl;}
        else if(d_input_signal>=0.000000002) {input_signal=1; adc_of_bit++;}//viv<<1<<endl;} //0.46755 //0.12
            else if(d_input_signal<0 && d_input_signal>-0.000000002) {input_signal=-1;}//viv<<0<<endl;}
            else if(d_input_signal<=-0.000000002) {input_signal=-1;}// adc_of_bit++;}//viv<<1<<endl;} //0.46755
        */
#endif

#if adc_agc_mode==0
        //! Generation of signals for first complex multiplier
        base_cosine=int(floor(cos(double(t)*base_frequency)*((1<<11) - 1)+0.5));
        base_sine=int(floor(sin(double(t)*base_frequency)*((1<<11) - 1)+0.5));

        for(short i = 15; i>0; i--) delay_line_cos[i] = delay_line_cos[i-1];
        delay_line_cos[0] = Flt.NCO(nco_freq,1);
        cos_tmp = delay_line_cos[10];

        for(short i = 15; i>0; i--) delay_line_sin[i] = delay_line_sin[i-1];
        delay_line_sin[0] = Flt.sin_out;
        sin_tmp = delay_line_sin[10];
        //viv2<<d_input_signal<<endl;
        //! Complex multiplying and rounding
        //complex_signal_i=int(floor(double(input_signal*cos_tmp)/(1<<10)+0.5));
        complex_signal_i=int(floor(double(input_signal*base_cosine)/(1<<10)+0.5));//(input_signal*base_cosine)/8;//
        //complex_signal_q=int(floor(double(input_signal*sin_tmp)/(1<<10)+0.5));
        complex_signal_q=int(floor(double(input_signal*base_sine)/(1<<10)+0.5));

        if(complex_signal_i == 1024) complex_signal_i = 1024-1;
        if(complex_signal_q == 1024) complex_signal_q = 1024-1;
        //if(new_epoch>=1 && new_epoch<=1) str<<input_signal<<' ';

        //! Complex filtering for AGC
        int complex_signal_agc = 0;
        int outren_tmp_agc = 0;
        if(t_freq_move_fast==4) t_freq_move_fast = 0;
        int half_band_signal_i = 0;
        int half_band_signal_q = 0;
        #if LOCATA==0
        switch(t_freq_move_fast)
        {
            case 0  : half_band_signal_i = Flt.FirH(input_signal,0,11,0); break;
            case 1  : half_band_signal_i = Flt.FirH(0,0,11,1); break;
            case 2  : half_band_signal_i = Flt.FirH(-input_signal,0,11,0); break;
            case 3  : half_band_signal_i = Flt.FirH(0,0,11,1); break;
        }
        switch(t_freq_move_fast)
        {
            case 0  : half_band_signal_q = Flt_Q.FirH(0,0,11,1); break;
            case 1  : half_band_signal_q = Flt_Q.FirH(input_signal,0,11,0); break;
            case 2  : half_band_signal_q = Flt_Q.FirH(0,0,11,1); break;
            case 3  : half_band_signal_q = Flt_Q.FirH(-input_signal,0,11,0); break;
        }
        outren_tmp_agc = half_band_signal_i*half_band_signal_i+half_band_signal_q*half_band_signal_q;
        #endif
        //norm<<half_band_signal_i<<endl;

        t_freq_move_fast++;

        //if(outren_tmp_agc>=(511*511)) adc_of_bit++;
        //if(outren_tmp_agc>=(490*490)) adc_high_bit++;
        //if(outren_tmp_agc>=(470*470)) adc_low_bit++;
        if(outren_tmp_agc>=(511*511)) adc_of_bit++;//511
        if(outren_tmp_agc>=(250*250)) adc_high_bit++;//470
        if(outren_tmp_agc>=(200*200)) adc_low_bit++;//300

        #if pwm_dac==0
        real_gain_coeff=inp_agc.PWM(k_agc_deay[2],16,128);//160 - wide, 50 - narrow
        #else
        //real_gain_coeff=inp_agc.DAC(k_agc_deay[2],0);

        short porog_1 = 10;
        short porog_2 = 126;
        short value_1 = 70<<7;
        short value_2 = 16<<7;
        short bit_zero = 0;
        short bit_one = 1;
        short trm_enable = 1;

        if(trm_enable == 1)
        {
            short gain_mode_pre = gain_mode;

            if(k_agc_deay[2] < porog_1) gain_mode = bit_zero;
            if(k_agc_deay[2] > porog_2) gain_mode = bit_one;

            if((gain_mode - gain_mode_pre) < 0) inp_agc.gain_coeff = value_1;
            if((gain_mode - gain_mode_pre) > 0) inp_agc.gain_coeff = value_2;
        }

        real_gain_coeff=inp_agc.DACnew(k_agc_deay[2],gain_mode);
        //viv2<<k_agc_deay[2]<<' '<<gain_mode<<' '<<d_input_signal<<' '<<input_signal<<endl;
        #endif
        k_agc_deay[2]=k_agc_deay[1];
        k_agc_deay[1]=k_agc_deay[0];
        k_agc_deay[0]=agc_out;

        if(abs(input_signal)>max_point_value) max_point_value = abs(input_signal);
        //! AGC============================================================================
        if(t%128==0)
        {
            if(adc_high_bit==0) ref_mid_level += 512;
            else if(adc_of_bit>0) ref_mid_level -= 512;
            else ref_mid_level = agc_mid_level;
            if(ref_mid_level<250*512) ref_mid_level = 250*512;
            if(ref_mid_level>410*512) ref_mid_level = 410*512;

            agc_out = inp_agc.GainCoefficient_CordicNew(adc_of_bit,adc_high_bit,adc_low_bit,agc_mid_level,ref_mid_level,0,k_cordic);
            //agc_out = inp_agc.GainCoefficient(adc_of_bit,adc_high_bit,agc_mid_level,0);
            //agc_out = inp_agc.GainCoefficient_CordicNew2(adc_of_bit,adc_high_bit,adc_low_bit,agc_mid_level,ref_mid_level,0,k_cordic,max_point_value);
            //cout<<(agc_out/1024)<<endl;//<<' '<<adc_of_bit_av<<' '<<adc_high_bit_av<<' '<<' '<<agc_mid_level_av<<' '<<ref_mid_level<<endl;
            //if(new_epoch>=124 && new_epoch<=128)
            //str<<agc_out<<' '<<adc_of_bit<<' '<<adc_high_bit<<' '<<adc_low_bit<<' '<<agc_mid_level<<' '<<ref_mid_level<<endl;
            //agc_out = 1590;
            adc_of_bit = 0;
            agc_mid_level = 0;
            adc_high_bit = 0;
            adc_low_bit = 0;

            max_point_value = 0;
        }
        //viv2<<d_input_signal<<' '<<input_signal<<' '<<complex_signal_i<<endl;

        for(short i = 15; i>0; i--) delay_line_sign_i[i] = delay_line_sign_i[i-1];
        delay_line_sign_i[0] = complex_signal_i;
        complex_signal_i = delay_line_sign_i[4];

        for(short i = 15; i>0; i--) delay_line_sign_q[i] = delay_line_sign_q[i-1];
        delay_line_sign_q[0] = complex_signal_q;
        complex_signal_q = delay_line_sign_q[4];

        if(rf_I_en==1) test_out_1<<complex_signal_i<<endl;
        //norm<<complex_signal_i<<" "<<complex_signal_q<<endl;
        //if(rf_Q_en==1) test_out_2<<complex_signal_q<<endl;
        //complex_signal_i=int(floor(double(input_signal*Flt.cos_out)/(1<<10)+0.5));//(input_signal*base_cosine)/8;//
        //complex_signal_q=int(floor(double(input_signal*Flt.sin_out)/(1<<10)+0.5));
#else
        base_cosine=short(floor(cos(double(t)*base_frequency)*3+0.5));
        base_sine=short(floor(sin(double(t)*base_frequency)*3+0.5));
        complex_signal_i=(input_signal*base_cosine);
        complex_signal_q=(input_signal*base_sine);
#endif
        #endif
#if CIC_filter==1
        #if LOCATA==0
        //! CIC filter
        for(short i = 63; i>0; i--) delay_line_cics[i] = delay_line_cics[i-1];
        delay_line_cics[0] = complex_signal_i;

        //if(new_epoch>=1 && new_epoch<=1) str<<complex_signal_i<<' ';
        int outren_tmp = Flt.CICFilter(complex_signal_i);
        //int outren_tmp = 0;

        for(short i = 35; i>0; i--) delay_line_out_4[i] = delay_line_out_4[i-1];
        delay_line_out_4[0] = outren_tmp;
        outren_tmp = delay_line_out_4[12];
        outren_d = delay_line_out_4[33];

        if(d6_en==1) test_out_4<<outren_tmp<<endl;

        int outimn_tmp = Flt_Q.CICFilter(complex_signal_q);
        //int outimn_tmp = 0;

        for(short i = 35; i>0; i--) delay_line_out_4q[i] = delay_line_out_4q[i-1];
        delay_line_out_4q[0] = outimn_tmp;
        outimn_tmp = delay_line_out_4q[12];
        outimn_d = delay_line_out_4q[33];

        //if(d6_en==1) test_out_4<<outimn_tmp<<endl;
        #if AVT2==1
        //viv2<<outren_tmp<<" "<<outimn_tmp<<endl;
        #endif
        //! First middle digital AGC
        outren = Flt.AGCMid2(outren_tmp,outimn_tmp,outren_d,outimn_d,t);
        outimn = Flt.agcmidout_q;
        if(co_rnd_en==1) test_out_5<<outren_tmp<<endl;
        //str<<outren<<endl;
        /*
        viv<<complex_signal_i<<' '<<outren<<' ';
        outren_tmp = FltPos.CICFilter(delay_line_cics[63]*3-outren*2/3);
        outimn_tmp = FltPos_Q.CICFilter(delay_line_cics[63]*3-outimn*2/3);
        outren = FltPos.AGCMid2(outren_tmp,outimn_tmp,outren_d,outimn_d,t);
        outimn = FltPos.agcmidout_q;
        viv<<outren<<' ';
        outren_tmp = FltPos2.CICFilter(outren/3);
        outimn_tmp = FltPos2_Q.CICFilter(outimn/3);
        outren = FltPos2.AGCMid2(outren_tmp,outimn_tmp,outren_d,outimn_d,t);
        outimn = FltPos2.agcmidout_q;
        viv<<outren<<endl;*/
        #endif
#else
        //outre=Flt.FirI(complex_signal_i);
        //outim=Flt.FirQ(complex_signal_q);
#endif
#if adc_agc_mode==0
        #if LOCATA==0
        rotate(&outp1[0],&outp1[input_sh_delay-1],&outp1[input_sh_delay]);
        outp1[0]=short(sh);
        //for(short i = 1; i<input_sh_delay; i++) outp2[i] = outp1[i-1];
        //str<<outren<<endl;
        //! Signal decimation
        if(((t+2)%filt_div)==0)
        #endif
        //if(0)]
        {
            #if LOCATA==0
            for(short i = 7; i>0; i--) delay_line_out_6[i] = delay_line_out_6[i-1];
            delay_line_out_6[0] = outren;
            outren = delay_line_out_6[1];

            //! FIR filtering
            if(cic_out_en==1) test_out_6<<outren<<endl;
            outren=Flt.FirI(outren,0,0); //8 3

            for(short i = 25; i>0; i--) delay_line_out_7[i] = delay_line_out_7[i-1];
            delay_line_out_7[0] = outren;
            outren = delay_line_out_7[5];
            outren_d = delay_line_out_7[21];

            if(sfir_en==1) test_out_7<<outren<<endl;

            for(short i = 7; i>0; i--) delay_line_out_6q[i] = delay_line_out_6q[i-1];
            delay_line_out_6q[0] = outimn;
            outimn = delay_line_out_6q[1];

            if(cic_out_en==1) test_out_6<<outimn<<endl;
            outimn=Flt_Q.FirI(outimn,0,0);

            for(short i = 25; i>0; i--) delay_line_out_7q[i] = delay_line_out_7q[i-1];
            delay_line_out_7q[0] = outimn;
            outimn = delay_line_out_7q[5];
            outimn_d = delay_line_out_7q[21];
            #if AVT2==1
            //viv2<<outren<<" "<<outimn<<endl;
            #endif
            //str<<outren<<endl;
#if anti_jamm_on==1
            //! Second middle digital AGC
            //if(new_epoch>=1 && new_epoch<=6) viv2<<outren<<' '<<outimn<<' ';
            outren = Flt_Q.AGCMids(outren,outimn,outren_d,outimn_d,t/filt_div);
            outimn = Flt_Q.agcmidout_q;

            for(short i = 7; i>0; i--) delay_line_out_8[i] = delay_line_out_8[i-1];
            delay_line_out_8[0] = outren;
            outren = delay_line_out_8[1];
            if (aj_out_en == 1) test_out_8 << outren << endl;

            for(short i = 7; i>0; i--) delay_line_out_8q[i] = delay_line_out_8q[i-1];
            delay_line_out_8q[0] = outimn;
            outimn = delay_line_out_8q[1];

            for(short i = 7; i>0; i--) delay_line_out_9[i] = delay_line_out_9[i-1];
            delay_line_out_9[0] = outren;
            outren = delay_line_out_9[1];

            if (qnt_en==1) test_out_9 << outren << endl;

            for(short i = 7; i>0; i--) delay_line_out_9q[i] = delay_line_out_9q[i-1];
            delay_line_out_9q[0] = outimn;
            outimn = delay_line_out_9q[1];
            //str<<outren<<endl;
#else
            outren=inp_agc.AGC_second(outren,0);//outimn
#endif
            //! Complex to real Converter
            //norm<<outren<<' '<<outimn<<endl;
            //outren_pre = outren;
            if(t_freq_move==4) t_freq_move = 0;
            #if ajm_real==1
            switch(t_freq_move)
            {
                case 0  : outren = outren; break;
                case 1  : outren = outimn; break;
                case 2  : outren = -outren; break;
                case 3  : outren = -outimn; break;
            }
            #endif
            #endif
            //int outrensh = outren*1 + test0.RandValue()*500;
            //outren = outren*1 + test1.RandValue()*500;
            //viv2<<outren<<endl;
            if(new_epoch>=318 && new_epoch<=322)viv<<outren<<' ';
            //outren=Flt.AJFilterNew2(outren,16,0,point);
            //rotate(&a_fft[0],&a_fft[63],&a_fft[64]);
            //a_fft[0] = outren;
            //if(t>100000)
            /*if((((t+2)%(filt_div*64)))==0)
            {

                Flt.FFT(a_fft,b_fft,Nfft);
                for(short i=0;i<Nfft;i++)
                {
                    //norm
                    fft_data0[i] = Flt.GetFFT(i,0);
                    fft_data1[i] = Flt.GetFFT(i,1);
                    //a[i] = Flt.GetFFT(i,0);
                    //b[i] = Flt.GetFFT(i,1);
                    if(i>15 && i<49) fft_filt[i] = 0;
                            else fft_filt[i] = 1;
                    fft_data0[i] = fft_data0[i] * fft_filt[i];
                    fft_data1[i] = fft_data1[i] * fft_filt[i];
                }
                Flt.FFT(fft_data1,fft_data0,Nfft);
                for(short i=0;i<Nfft;i++)
                {
                    fft_data0[i] = Flt.GetFFT(i,0)/64;
                    fft_data1[i] = Flt.GetFFT(i,1)/64;
                }
            }*/
            //norm<<outren<<endl;
            //short t_fft = 63 - (((t+2)/filt_div) & 0x3F);
            //str<<t_fft<<endl;
            //outren = fft_data1[t_fft];
            //norm<<outren<<' '<<fft_data0[t_fft]<<endl;
            //norm<<outren<<' ';
            //norm<<FltPos.IirI(outren,0,0)<<endl;
			if(((t+2)/filt_div)>=50000 && ((t+2)/filt_div)<50512) fft_in[((t+2)/filt_div)-50000] = outren;
			if(((t+2)/filt_div)==50512)
			{
				double *dfourier = FFTw(fft_in,256);
				double max_el_fft = *max_element(&dfourier[0],&dfourier[128]);
				double max_k = 128;
				for(short k=0; k<128; k++)
				{
					fft_out[k] = 20*log10(dfourier[k]/max_el_fft);
					viv<<fft_in[k]<<' '<<dfourier[k]<<endl;
				}
				for(short k=0; k<128; k++)
				{
					if(abs(fft_out[k]) < 0.01) max_k = k;
				}
				cout<<"max_k = "<<max_k<<endl;
				for(short k=0; k<128; k++)
				{
					if((fft_out[short(max_k)] - fft_out[k]) < 15) cout<<"k2 = "<<k<<endl;
				}
			}
            if(pos_filt_en == 1)
            {
                //norm<<outren<<' ';
                outren=FltPos.FirIpos(outren,0,0);

                outren = FltPos.AGCMidsPos(outren,0,outren,0,t/filt_div);

                //outren=FltPos2.FirIpos(outren,0,0);

                //outren = FltPos2.AGCMidsPos(outren,0,outren,0,t/filt_div);
                //norm<<outren<<' ';
                //outren=FltPos2.FirIpos(outren,0,0);
                //norm<<outren<<endl;
                //outren = FltPos2.AGCMidsPos(outren,0,outren,0,t/filt_div);
            }

            abs_outren_0 = abs(outren);
            aj_coeff_sum_0 += abs_outren_0;
            if(abs_outren_0 > w_aj_abs_max_0) w_aj_abs_max_0 = abs_outren_0;

            //if(((t+2)%100)==0)
            //	outren=Flt.AJFilterNew2(outren,16,1,point);
            //else
                outren=Flt.AJFilterNew2(outren,16,0,point);
                //Flt_Q.AJFilterNew2(outren,16,1,point);
            abs_outren_1 = abs(outren);
            aj_coeff_sum_1 += abs_outren_1;
            if(abs_outren_1 > w_aj_abs_max_1) w_aj_abs_max_1 = abs_outren_1;

            jam_detect_cnt--;
            if(jam_detect_cnt==-1)
            {
                jam_osn_thr_int_0 = int(floor(double(aj_coeff_sum_0/256)/double(w_aj_abs_max_0)*256+0.5));
                jam_osn_thr_int_1 = int(floor(double(aj_coeff_sum_1/256)/double(w_aj_abs_max_1)*256+0.5));
                jam_osn_thr_sum_0 += jam_osn_thr_int_0;
                jam_osn_thr_sum_1 += jam_osn_thr_int_1;
                jam_osn_thr_mulsum_0 += jam_osn_thr_int_0*jam_osn_thr_int_0;
                jam_osn_thr_mulsum_1 += jam_osn_thr_int_1*jam_osn_thr_int_1;
                //viv2<<(double(aj_coeff_sum_0/256)/double(w_aj_abs_max_0))<<' ';
                //viv2<<(double(aj_coeff_sum_1/256)/double(w_aj_abs_max_1))<<endl;
                if(jam_osn_thr_int_0 > 102) jam_first_thr_0++;
                if(jam_osn_thr_int_1 > 102) jam_first_thr_1++;

                jam_detect_cnt = 255;
                aj_coeff_sum_0 = 0;
                w_aj_abs_max_0 = 0;
                aj_coeff_sum_1 = 0;
                w_aj_abs_max_1 = 0;
                jam_detect_number++;
            }

            //if(new_epoch_pre != new_epoch)
			if(jam_detect_number==128 || new_epoch_pre != new_epoch)
            {
                int m_x_0 = jam_osn_thr_sum_0/jam_detect_number;
                int m_x_1 = jam_osn_thr_sum_1/jam_detect_number;
                int m_x2_0 = jam_osn_thr_mulsum_0/jam_detect_number;
                int m_x2_1 = jam_osn_thr_mulsum_1/jam_detect_number;
                int d_x_0 = m_x2_0 - m_x_0*m_x_0;
                int d_x_1 = m_x2_1 - m_x_1*m_x_1;
                //viv2<<double(m_x_0)/256.0<<' ';
                //viv2<<double(m_x_1)/256.0<<' ';
                //viv2<<double(d_x_0)/256.0<<' ';
                //viv2<<double(d_x_1)/256.0<<endl;
                if(jam_first_thr_0>0) jam_exist_0 = 1;
                else if(d_x_0<63) jam_exist_0 = 1;
                else jam_exist_0 = 0;
                if(jam_first_thr_1>0) jam_exist_1 = 1;
                else if(d_x_1<63) jam_exist_1 = 1;
                else jam_exist_1 = 0;
				if(jam_detect_number==128)
				{
				viv2<<jam_first_thr_0<<' ';
				viv2<<d_x_0<<' ';
                viv2<<jam_first_thr_1<<' ';
                viv2<<d_x_1<<endl;
				}
				
					jam_detect_number = 0;
					jam_osn_thr_sum_0 = 0;
					jam_osn_thr_sum_1 = 0;
					jam_osn_thr_mulsum_0 = 0;
					jam_osn_thr_mulsum_1 = 0;
					jam_first_thr_0 = 0;
					jam_first_thr_1 = 0;

					jam_detect_cnt = 255;
					aj_coeff_sum_0 = 0;
					w_aj_abs_max_0 = 0;
					aj_coeff_sum_1 = 0;
					w_aj_abs_max_1 = 0;
				
            }
            new_epoch_pre = new_epoch;
            //norm<<outren<<endl;
            //outren=Flt.AJFilterNew3(outren,outrensh,16,0,point);
            if(new_epoch>=318 && new_epoch<=322)viv<<outren<<endl;
            #if AVT2==1
            //viv<<outren<<endl;
            #endif
            //cout<<((outren>>3)<<0)<<endl;
            #if LOCATA==1
            rotate(&outp1[0],&outp1[input_sh_delay-1],&outp1[input_sh_delay]);
            outp1[0]=Flt_Q.code_p ? 1 : -1;
            rotate(&outp2[0],&outp2[input_sh_delay-1],&outp2[input_sh_delay]);
            outp2[0]=Flt_Q.strb_h;
            //outren = inp_agc.AGC_second((sh*cos(double(t)*12.5e6*2*pi*sampling_interval))*8);
            locata_in = (Flt.code_p==1) ? cos(double(t)*10.5e6*2*pi*sampling_interval*1) : -cos(double(t)*10.5e6*2*pi*sampling_interval*1);
            //norm << Flt.IirIn(locata_in) << endl;
            //norm<<locata_in<<' ';
            //str<<locata_in<<' ';
            //locata_in = (double(Flt.FirIn((floor(locata_in*2047)+0.5),0,0))/512/2 + 1*Flt.NoiseGenLocata())/2047;
            //locata_in = (double(Flt.FirI((floor(locata_in*2047)+0.5),0,0))/512/2 + 0*locata_sh.RandValue())/1;
            //locata_in = (double(Flt_Q.FirI((floor(locata_in)+0.5),0,0))/512/2 + 0*locata_sh.RandValue())/1;
            //if((t%200)==0)norm<<locata_in<<' ';

            locata_noise_sigma = 1;
            double noise_RF = locata_RF_noise.RandValue();
            locata_in = Flt_Q.IirIn(noise_RF)*locata_noise_sigma;
            //norm << Flt.IirIn(outp1[0]) << endl;
            //norm<<locata_in<<endl;
            //locata_in = Flt.locata_out;
            //if(t%450==0) locata_pot = floor(locata_amp_jump.RandValue()+0.5)*10*0+100;
            //locata_amp = sqrt(locata_noise_sigma*locata_noise_sigma*pow(10.0,locata_pot/10.0)*2/50e6);
            locata_in += locata_amp * Flt.locata_out * cos(double(t)*10.5e6*2*pi*sampling_interval)*1 + 0*cos(double(t)*10.4e6*2*pi*sampling_interval);
            locata_in = locata_in/1;
            //norm<<Flt.locata_out<<' '<<Flt_Q.IirIn(locata_in)<<endl;
            //if(((t+2)%3)==0)
            //{
            //norm<<locata_in<<' ';
            locata_in /= pow(10.,3-agc_k_locata);
            //norm<<locata_in<<' '<<agc_k_locata<<endl;
            //!LOCATA ADC
            //locata_in = floor(locata_in*8191)+0.5 + 0*locata_sh.RandValue();
            locata_in_dig=short(floor(locata_in*8191)+0.5);// + 0*locata_sh.RandValue());
            if(locata_in_dig>=(8192-1))
            {
                locata_in_dig = (8192-1);
            }
            else if(locata_in_dig<=-8192)
            {
                locata_in_dig = -8192;
            }
            agc_k_locata = Flt.AGCLocata(abs(locata_in_dig),t);

            //norm<<locata_in<<' '<<locata_in_dig<<' '<<agc_k_locata<<endl;
            //locata_in = double((floor(locata_in*2047)+0.5)+0*locata_sh.RandValue())/2048;
            //locata_in = (double(Flt.FirIn((floor(locata_in*2047)+0.5),0,0))/512/2*0+0.6597*2047 + 1*locata_sh.RandValue())/2047;
            //norm<<Flt.NoiseGenLocata()<<endl;
            //norm<<locata_sh.RandValue()<<endl;
            //locata_in = locata_in * 0.3;
            //if(t%2==0) locata_in = double(Flt.FirH((floor(locata_in*2048)+0.5),0,11,0))/2048;
            //else locata_in = double(Flt.FirH(0,0,11,1))/2048;
            //norm<<locata_in<<' ';
            //str<<locata_in<<' '<<outp1[0]<<' '<<outp2[0]<<endl;
            outren = 0;
            if(locata_in>=0)
            {
                if(locata_in>293) outren = 1;//0.14285 293
                if(locata_in>878) outren = 2;//0.42855 878
                if(locata_in>1463) outren = 3;//0.71425 1463
                //outren = 3;
            } else
            {
                if(locata_in<-293) outren = -1;
                if(locata_in<-878) outren = -2;
                if(locata_in<-1463) outren = -3;
                //outren = -3;
            }
            /*if(locata_in>=0)
            {
                outren = 1;
                if(locata_in>0.33) outren = 2;
                if(locata_in>0.66) outren = 3;
            } else
            {
                outren = -1;
                if(locata_in<-0.33) outren = -2;
                if(locata_in<-0.66) outren = -3;
            }*/
            //norm<<outren<<endl;
            //norm<<outren<<' '<<locata_in<<endl;
            //if(outren<3 && outren>-3) cout<<outren<<endl;
            //cout<< (512*sh*cos(double(t)*25e6*2*pi*sampling_interval))*8<< endl;
            outren = Flt.AGCSecondLocata(locata_in_dig, t, 1, 1, locata_sh.RandValue(), 10);
            #else
            //norm<<outren<<endl;
            //outren = outren/256;
            //if(outren>3) outren = 3;
            //if(outren<-3) outren = 0;
            //norm<<outren<<endl;
            outren=inp_agc.AGC_second(((outren>>0)<<3));
            #if AVT2==1
            viv2<<outren<<endl;
            #endif
            #endif
            //outren=outren>>3;

            //! Finish bit-true model
            #if ajm_real==1
            switch(t_freq_move)
            {
                case 0  : complex_signal_i = outren; break;
                case 1  : complex_signal_i = 0; break;
                case 2  : complex_signal_i = -outren; break;
                case 3  : complex_signal_i = 0; break;
            }
            switch(t_freq_move)
            {
                case 0  : complex_signal_q = 0; break;
                case 1  : complex_signal_q = outren; break;
                case 2  : complex_signal_q = 0; break;
                case 3  : complex_signal_q = -outren; break;
            }
            #else
            complex_signal_i = outren;
            complex_signal_q = outimn;
            #endif
            //! Signal accumulating
            //filter_delay = 20;
            #if LOCATA==1
            //if(new_epoch>=500 && outp2[30]!=0) norm<<outren<<' '<<locata_in<<endl;
            //norm<<int(Flt_Q.nco_reg&0xffffffff)<<' '<<outp1[0]<<' '<<outp2[0]<<endl;
            //if(outp2[0]!=0) str<<(int(Flt_Q.nco_reg&0xff000000)>>24)<<endl;

            double loc_opor_cos = cos(double(t)*10.5e6*2*pi*sampling_interval+double(-1)/100*1);//-64 -13
            short locata_opor_cos = 0;
            if(loc_opor_cos>0)
            {
                if(loc_opor_cos>0.33) locata_opor_cos = 1;
                if(loc_opor_cos>0.66) locata_opor_cos = 2;
                //if(loc_opor_cos>0.75) locata_opor_cos = 3;
            } else
            {
                if(loc_opor_cos<-0.33) locata_opor_cos = -1;
                if(loc_opor_cos<-0.66) locata_opor_cos = -2;
                //if(loc_opor_cos<-0.75) locata_opor_cos = -3;
            }
            double loc_opor_sin = sin(double(t)*10.5e6*2*pi*sampling_interval+double(-1)/100*1);
            short locata_opor_sin = 0;
            if(loc_opor_sin>0)
            {
                if(loc_opor_sin>0.33) locata_opor_sin = 1;
                if(loc_opor_sin>0.66) locata_opor_sin = 2;
                //if(loc_opor_cos>0.75) locata_opor_cos = 3;
            } else
            {
                if(loc_opor_sin<-0.33) locata_opor_sin = -1;
                if(loc_opor_sin<-0.66) locata_opor_sin = -2;
                //if(loc_opor_cos<-0.75) locata_opor_cos = -3;
            }
            //complex_signal_i = int(floor(double(outren)*loc_opor_cos+0.5));
            //complex_signal_q = int(floor(double(outren)*loc_opor_sin+0.5));
            complex_signal_i = outren*locata_opor_cos;
            complex_signal_q = outren*locata_opor_sin;

            locata_in = (Flt.code_p==1) ? cos(double(t)*10.5e6*2*pi*sampling_interval*1) : -cos(double(t)*10.5e6*2*pi*sampling_interval*1);
            //if(outp2[0] !=0)
                //norm << Flt.IirIn(locata_in)*locata_opor_cos << endl;
                outp1[0]=Flt.code_p ? 1 : -1;
                double test_norm = Flt.IirIn(outp1[0]);
                if(t%10==0) norm << test_norm << ' ' << outp1[0] << endl;
                //norm << Flt.locata_out << endl;
            //complex_signal_i = outren;
            //complex_signal_q = outren;
            //norm<<outp1[30]<<' '<<outp2[30]<<endl;
            //complex_signal_i = Flt.FirIn((complex_signal_i+0.5),0,0)/512/2;
            //str<<complex_signal_i<<endl;
            //complex_signal_i = int(floor(locata_in*2048)+0.5)*locata_opor_cos;
            //complex_signal_i = int(floor(loc_opor_cos*8192)+0.5)*outren;
            outm+=complex_signal_i*outp1[2];//30
            outm2+=complex_signal_i*outp2[2];
            //outm2+=complex_signal_i*outp2[32];
            //outm2+=complex_signal_q*Flt_Q.code_p;strb_h
            //str<<complex_signal_i<<' '<<sh2<<' '<<outren<<' '<<sh<<endl;
            //norm<<Flt_Q.code_p<<' '<<Flt_Q.strb_h<<endl;
            //}
            #else
            outm+=complex_signal_i*outp1[filter_delay];
            outm2+=complex_signal_q*outp1[filter_delay];
    //viv<<complex_signal_i<<' '<<outp1[filter_delay]<<endl;
            #endif
            //base_cosine=int(floor(cos(double(t/(filt_div))*0.05*1e6*2*pi*sampling_interval*4)*((1<<11) - 1)+0.5));
            //outm+=complex_signal_i;
            //outm2+=complex_signal_q;
            //outm+=complex_signal_i*(outp2[filter_delay] - outp1[filter_delay])/2;
            //outm2+=complex_signal_q*(outp2[filter_delay] - outp1[filter_delay])/2;
            t_freq_move++;
            //if(t%(50000*20)==0 && t!=0)
            //if(new_epoch%20==0) filter_delay = 117+4*(new_epoch/20);
            //if(t>300*2e5 && t<308*2e5)
            //norm<<outp1[filter_delay]<<' '<<outp2[filter_delay]<<' '<<outp1[filter_delay] - outp2[filter_delay]<<endl;
            //norm<<complex_signal_i<<' '<<outp1[filter_delay]<<endl;
        }
#else
        //input_signal = input_signal / 32;
        //if(input_signal<-3) input_signal = -4;
        complex_signal_i=short(floor(double(input_signal*NCO_gen_cos_phase(base_frequency,phase_in,t))/32.0+0.0));//(input_signal*base_cosine)/8;//
        complex_signal_q=short(floor(double(input_signal*(-NCO_gen_sin_phase(base_frequency,phase_in,t)))/32.0+0.0));
        //cout<<input_signal<<' '<<NCO_gen_cos_phase(base_frequency,phase_in,t)<<endl;
        outm+=complex_signal_i*short(sh);//+23
        outm2+=complex_signal_q*short(sh);
#endif
        sts++;
        stc++;
        ct++;
        //sinc_inp_si=0;
        //sinc_inp_sq=0;
    }

    v=v-correlate_offset;
    cout<<"Poteri "<<endl;
    viv<<"w     "<<w<<endl;
    viv<<"App   "<<jamming_amplitude<<endl;
    viv<<"v     "<<v<<endl;
    viv<<"N fil "<<l<<endl;
    cout<<"w     "<<w<<endl;
    cout<<"App   "<<jamming_amplitude<<endl;
    cout<<"v     "<<v<<endl;
    cout<<"ADC "<<double(adc_of_bit)/(sampling_frequency*1000*loe)*100<<endl;

    //! Calculating SNR
    double xmid1=0,dmid1=0,smid1=0,xmid2=0,dmid2=0,smid2=0,dmid3=0;

    xmid1=accumulate<double *__w64,double>(&nmetod1[0],&nmetod1[v],0)/double(v);
    xmid2=accumulate<double *__w64,double>(&nmetod2[0],&nmetod2[v],0)/double(v);

    transform(&nmetod1[0], &nmetod1[v], &nmetod1[0], bind2nd(minus<double>(), xmid1));
    transform(&nmetod1[0], &nmetod1[v], &nmetod1[0], &nmetod1[0], multiplies<double>());
    dmid1=accumulate<double *__w64,double>(&nmetod1[0],&nmetod1[v],0)/double(v-1);
    transform(&nmetod2[0], &nmetod2[v], &nmetod2[0], bind2nd(minus<double>(), xmid2));
    transform(&nmetod2[0], &nmetod2[v], &nmetod2[0], &nmetod2[0], multiplies<double>());
    dmid2=accumulate<double *__w64,double>(&nmetod2[0],&nmetod2[v],0)/double(v-1);

    smid1=xmid1*xmid1/dmid1;
    smid2=xmid2*xmid2/dmid2;

    double signal_power = pow(signal_amplitude,2);
#if adc_agc_mode==0
    double ideal_snr=200e6*signal_power/2/1e3;
#else
    double ideal_snr=50e6*signal_power/2/1e3; //600
#endif
    cout<<"loss I "<<10*log10(smid1/ideal_snr)<<endl;
    viv<<"loss I "<<10*log10(smid1/ideal_snr)<<" parameter_0 = "<<parameter_0<<endl;
    norm<<-10*log10(smid1/ideal_snr)<<endl;
    if(10*log10(smid1/ideal_snr)<-100) goto LB3;
}
LB3:;
//LB2:
//cin>>z;
return 0;
}
