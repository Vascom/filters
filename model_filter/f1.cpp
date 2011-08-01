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

#define TRESHOLD_SWITCH_ALGORITHM 0
#define LOCATA 0
#define AVT2 0

#define V5_COMPAT
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include "class_Filter.h"
#include <limits.h>

using namespace std;

double gaussian_noise() 
{
    double tmp = 0;
    for(short i=0; i<12; i++) tmp += double(random()) - INT_MAX/2;
    return tmp/INT_MAX; 
}

double *FFTw( int *dIn, int nn )
{
    int i, j, n, m, mmax, istep;
    double tempr, tempi, wtemp, theta, wpr, wpi, wr, wi;

    int isign = -1;
    double *data = new double [ nn * 2 + 1 ];

    for( i = 0; i < nn; i++ )
    {
        data[ i * 2 ] = 0;
        data[ i * 2 + 1 ] = dIn[ i ];
    }
    data[512] = 0;
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
        i = i + 2;
    }

    mmax = 2;
    while( n > mmax )
    {
        istep = 2 * mmax;
        theta = 2.0 * pi / ( isign * mmax );
        wtemp = sin( 0.5 * theta );
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin( theta );
        wr = 1.0;
        wi = 0.0;
        m = 1;

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
                data[ i + 1 ] = data[ i + 1 ] + tempi;
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

//Generator E5aI PRN
void GALE5aIgen(short *e5a)
{
 short geng1[14]={1,1,1,1,1,1,1,1,1,1,1,1,1,1};
 short geng2[14]={0,1,1,0,1,1,0,0,1,0,1,0,1,0};
 short G=0,g1=0,g2=0,g1s=0,g2s=0;
 for(short n=0;n<10000;n++)
 {
  g1=geng1[13];
  g1s=geng1[0]^geng1[5]^geng1[7]^geng1[13];
  rotate(&geng1[0], &geng1[13], &geng1[14]);
  geng1[0]=g1s;
  g2=geng2[13];
  g2s=geng2[3]^geng2[4]^geng2[6]^geng2[7]^geng2[11]^geng2[13];
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
 short G=0,g1=0,g2=0,g1s=0,g2s=0;
 for(short n=0;n<10000;n++)
 {
  g1=geng1[13];
  g1s=geng1[3]^geng1[10]^geng1[12]^geng1[13];
  rotate(&geng1[0], &geng1[13], &geng1[14]);
  geng1[0]=g1s;
  g2=geng2[13];
  g2s=geng2[1]^geng2[4]^geng2[7]^geng2[8]^geng2[11]^geng2[13];
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

short NCO_gen_cos_phase(double freq, double phase_in, long long t)
{
    return short(floor(cos(double(t) * freq + phase_in)*127 + 0.5));
 //double phase = (((long long)floor(double(t) * freq / 2.0 / pi * 16384 + 0.5) % 16384) >> 4) * 2.0 * pi / 1024.0+phase_in;
 //if (((cos(phase) * sin(phase)) >= 0) & ((cos(phase + 2.0 * pi / 1024.0)  * sin(phase + 2.0 * pi / 1024.0)) >= 0)) {phase += 2.0 * pi / 1024.0;}
 //return short(floor(cos(phase)*127 + 0.5));
}
short NCO_gen_sin_phase(double freq, double phase_in, long long t)
{
    return short(floor(sin(double(t) * freq + phase_in)*127 + 0.5));
 //double phase = (((long long)floor(double(t) * freq / 2.0 / pi * 16384 + 0.5) % 16384) >> 4) * 2.0 * pi / 1024.0 + phase_in;
 //if (((cos(phase) * sin(phase)) <= 0) & ((cos(phase + 2.0 * pi / 1024.0)  * sin(phase + 2.0 * pi / 1024.0)) <= 0)) {phase += 2.0 * pi / 1024.0;}
 //return short(floor(sin(phase)*127 + 0.5));
}


//! Main program

int main(int argc, char *argv[]){
srandom( time( NULL ) );
short z=0;
ofstream norm, str, viv, viv2;
ofstream test_out_0,test_out_1,test_out_2,test_out_4,test_out_5,test_out_6,test_out_7,test_out_8,test_out_9;
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

unsigned int inp_prn_hex[32] = {
    0xC83949E5,0x13EAD115,0x591E9FB7,0x37CAA100,0xEA44DE0F,0x5CCF602F,0x3EA62DC6,0xF5158201,
    0x31D81C6, 0xFFA74B61,0x56272DD8,0xEEF0D864,0x906D2DE2,0xE0527E02,0xB9F5F331,0xC6D56C6E,
    0xE002CD9D,0xA0ABAE94,0x7389452D,0xADAD8E7, 0xB21F9688,0x7D5CC925,0xFF87DE37,0x2C3950A5,
    0x7E3DA767,0xEFA31F01,0x28B444D8,0x1DA3448E,0x2CC9E6FC,0xCA69AF36,0xA778D442,0x24E1CA20};

short inp_prn[1023];
for(short k = 0; k<32; k++)
    {
        for(short n = 31; n>=0; n--)
        {
            if(((inp_prn_hex[k]>>n)&0x1) > 0) inp_prn[k*32+(31-n)] = 1;
            else inp_prn[k*32+(31-n)] = -1;
        }
    }

double fr0 = 0, fr1 = 0, fr;
cout<<"Enter fr: ";
//cin>>fr;
fr = 10;
fr0 = 13+fr*4;
fr1 = fr0+3.5;

short input_bw = short(fr);
short input_mult = 1<<(short(fr)-input_bw);
int parameter_0 = 4096;//32768;//256;
int parameter_1 = 131072;
//short agc_algorithm = 5;
short params_file_number = 0;
int del_filt = 0;
short k_cordic = 2;

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
    //cout<<s<<' '<<f<<' '<<' '<<f1<<' '<<f2<<' '<<f3<<' '<<s-(f1+f2+f3)<<' '<<err<<endl;
}

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
    short k=0, filter_delay=0, kk=0;
    int input_signal=0, base_cosine=0, base_sine=0, inp_prn_ajm[6000], inp_prn_ajm2[6000];
    short  sts=0, stc=0, ct=0;
    short *outp1=NULL, *outp2=NULL, *e5a=NULL, *e5b=NULL;
    int complex_signal_i=0, complex_signal_q=0;
    long long loe=0;
    int v=0, outm=0, outm2=0;
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
    //!!StNormDev gaussian_noise, code, locata_RF_noise, test0, test1;
    double signal_amplitude=0;
    GALE5aIgen(e5a);
    GALE5bIgen(e5b);

    long pulse_width = 0;
    cout<<"Filter length: 64 "<<endl;
    cout<<"Enter pulse width: ";
    //cin>>pulse_width;
    short switch_filter_use=1;
    //! For filt_div = 8 : filter_delay = 265
    //! For filt_div = 8 : filter_delay = 425 new real signal
    #if ajm_real==1
    filter_delay=420 + 5-8 +8 -4 + 8 -232 -0 +64+16       +fil_pos_length*4*pos_filt_en*fil_pos_number;//-80
    #else
    filter_delay=265;//-80
    #endif
    //filter_delay=0;//-4
    switch_filter_use=2;

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
    short int filt_div;
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
        fscanf (config, "%s", &shift1 );fscanf( config, "%d", &jam_pulse_on );
        fscanf (config, "%s", &shift1 );fscanf( config, "%d", &jam_pulse_width );
        fscanf (config, "%s", &shift1 );fscanf( config, "%d", &jam_pulse_period );
        fscanf (config, "%s", &shift1 );fscanf( config, "%d", &jam_fm_on );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jam_fm_freq );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jam_fm_speed );
        fscanf (config, "%s", &shift1 );fscanf( config, "%d", &jam_am_on );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jam_am_freq );
        fscanf (config, "%s", &shift1 );fscanf( config, "%lf", &jam_am_depth );

        for(short i=0;i<2;i++)fscanf( config, "%s", &shift1 );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &nco_freq );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &filt_div );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &cic_del_m[0] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &cic_del_m[1] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &cic_del_m[2] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &cic_del_m[3] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &cic_del_m[4] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &cic_del_m[5] );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &filt_cic_out_shift );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &mid_agc_lev_n );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &mid_agc_lev_p );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &mid_agc_count_n );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &mid_agc_count_p );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &mid_agc_man );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &lms_mju );
        for(short i=0;i<30;i++)
        {
            fscanf( config, "%s", &shift1 );
            fscanf( config, "%hd", &h1[i] );
        }
        for(short i=0;i<10;i++)
        {
            fscanf( config, "%s", &shift1 );
            fscanf( config, "%hd", &h2[i] );
        }
        for(short i=0;i<fil_pos_length;i++)
        {
            fscanf( config, "%s", &shift1 );
            fscanf( config, "%hd", &h_pos[i] );
        }

        for(short i=0;i<22;i++)fscanf( config, "%s", &shift1 );
        fscanf( config, "%s", &shift1 );fscanf( config, "%ld", &log_length );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &rf_real_0_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &rf_I_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &rf_Q_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &comb_out_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &d6_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &co_rnd_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &cic_out_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &sfir_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &aj_out_en );
        fscanf( config, "%s", &shift1 );fscanf( config, "%hd", &qnt_en );

        fclose( config );
    }
    //jamming_amplitude = jamming_amplitude * pow(1.7783,k_cycle);
    //! Maiking classes for filters and filters initialization
    Filter Flt,Flt_Q;
    Filter FltPos,FltPos2;//,FltPos_Q,FltPos2_Q;
    AGC inp_agc;

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
    //!!UnifDev jaja(1,80), locata_sh(-1,1), locata_amp_jump(8,16);
    //InputPosled(inp_prn);
    Flt.LoadAJFilter(80,255,short(6+4+0), 10);
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

    double pomeha=0;
        outm=0;outm2=0;sts=0;stc=0;kk=0;v=0;ct=0;
    for(short g=0;g<input_sh_delay;g++)
    {
        outp1[g]=0;
        outp2[g]=0;
    }
    const short correlate_offset=50;
    loe=log_length;

    sh1=double(e5a[kk]+e5b[kk]);
    sh2=double(e5a[kk]-e5b[kk]);
    sh=1;
    double d_input_signal=0;
    short jamming_situation=0;
    if (jam_pulse_on==1) jamming_situation=5;

    cout << "signal_amplitude: " << signal_amplitude <<endl;
    cout << "nco_carr_frec: " << nco_carr_frec <<endl;
    cout << "jamming_frequency: " << jamming_frequency <<endl;

    cout << "jam_pulse_on: " << jam_pulse_on << endl;
    cout << "--jamming_situation: " << jamming_situation << endl;

    cout << "jam_am_on: " << jam_am_on << endl;
    cout << "jam_fm_on: " << jam_fm_on << endl;

    int outren=0,outimn=0,outren_d=0,outimn_d=0;
    outm=0;outm2=0;sts=0;stc=0;kk=0;v=0;ct=0;
    int adc_of_bit=0, agc_mid_level=0, adc_high_bit=0, adc_low_bit=0;
    long ref_mid_level=320*512; //vk
    k_cordic = 2;

    int agc_out=0;
    double real_gain_coeff=1,pomeha_FM=0,pomeha_AM=0;
    //!!StNormDev imp_jam_amp; //0,20,40,...,600
    //UnifDev imp_jam_dur(-1,1); //20mks-200mks
    //StUnifDev4::SetKernel(35189);
    //int imp_jam_d=int(imp_jam_dur.RandValue()), imp_jam_count=0;
    //double imp_jam_a=floor(fabs(imp_jam_amp.RandValue()))*20, gaussian_noise_2=0;

    short k_agc_deay[3]={0,0,0};
    short input_max = 1<<(input_bw-1);

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

    int t_freq_move = 0;
    int t_freq_move_fast = 0;

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
    short *avt2_data = NULL;
    avt2_data = new short[65536];
    //InputAVT2data(avt2_data);
    short gain_mode = 0;
    short jam_detect_cnt = 255;
    short jam_detect_number = 0;
    int aj_coeff_sum_0 = 0;
    int aj_coeff_sum_1 = 0;
    short abs_outren_0 = 0;
    short abs_outren_1 = 0;
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

    short Nfft = 64;
    int *fft_data0 = NULL;
    int *fft_data1 = NULL;
    fft_data0 = new int[Nfft];
    fft_data1 = new int[Nfft];
    for(short i = 0; i < Nfft; i++) fft_data0[i] = 0;
    for(short i = 0; i < Nfft; i++) fft_data1[i] = 0;

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
                //sh = double(rand());
                if(random() > INT_MAX/2)sh = 1;
                else sh = -1;             
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
                //sh2=double(rand());
                if(random() > INT_MAX/2) sh2 = 1;
                else sh2 = -1;
            }
        #if LOCATA==0
        //! Jamming situation setting
        switch(jamming_situation)
        {
            case 0: jamming_frequency=jamming_frequency + 0*4*0.25*0.125*1e0*2*pi*sampling_interval;
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
            /*case 7: if(imp_jam_count>=imp_jam_d)
                    {
                        imp_jam_d=int(imp_jam_dur.RandValue());
                        imp_jam_a=floor(fabs(imp_jam_amp.RandValue()))*20;
                        imp_jam_count=0;
                    }
                    imp_jam_count++;
                    jamming_amplitude=imp_jam_a; break;*/
            default: jamming_amplitude=0; break;
        }
        //pomeha_FM = jamming_amplitude*cos(jamming_frequency*double(t)+jam_fm_freq*2*pi*sampling_interval*sin(jam_fm_speed*2*pi*sampling_interval*t));
        pomeha_FM = jamming_amplitude*cos(jamming_frequency*double(t)+jam_fm_freq*sin(jam_fm_speed*2*pi*sampling_interval*t));
        pomeha_AM = jamming_amplitude*(1+jam_am_depth*cos(jam_am_freq*1e3*2*pi*sampling_interval*t))*cos(double(t)*jamming_frequency);
        pomeha = jam_am_on*pomeha_AM + jam_fm_on*pomeha_FM;
        //! Generate input radio signal
        d_input_signal = gaussian_noise()*1+1*sh*signal_amplitude*cos(double(t)*carrier_frequency);
        d_input_signal += pomeha + (jam_am_on-1)*(jam_fm_on-1)*jamming_amplitude*cos(double(t)*jamming_frequency);
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
                    if(fabs(fft_out[k]) < 0.01) max_k = k;
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
                outren=Flt.AJFilterNew2(outren,16,0,0);
                //Flt_Q.AJFilterNew2(outren,16,1,0);
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
            //outren=Flt.AJFilterNew3(outren,outrensh,16,0,0);
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
    }

    v=v-correlate_offset;
    cout<<"Poteri "<<endl;
    viv<<"App   "<<jamming_amplitude<<endl;
    viv<<"v     "<<v<<endl;
    cout<<"App   "<<jamming_amplitude<<endl;
    cout<<"v     "<<v<<endl;
    cout<<"ADC "<<double(adc_of_bit)/(sampling_frequency*1000*loe)*100<<endl;

    //! Calculating SNR
    double xmid1=0,dmid1=0,smid1=0,xmid2=0,dmid2=0,smid2=0;

    xmid1=accumulate<double *,double>(&nmetod1[0],&nmetod1[v],0)/double(v);
    xmid2=accumulate<double *,double>(&nmetod2[0],&nmetod2[v],0)/double(v);

    transform(&nmetod1[0], &nmetod1[v], &nmetod1[0], bind2nd(minus<double>(), xmid1));
    transform(&nmetod1[0], &nmetod1[v], &nmetod1[0], &nmetod1[0], multiplies<double>());
    dmid1=accumulate<double *,double>(&nmetod1[0],&nmetod1[v],0)/double(v-1);
    transform(&nmetod2[0], &nmetod2[v], &nmetod2[0], bind2nd(minus<double>(), xmid2));
    transform(&nmetod2[0], &nmetod2[v], &nmetod2[0], &nmetod2[0], multiplies<double>());
    dmid2=accumulate<double *,double>(&nmetod2[0],&nmetod2[v],0)/double(v-1);

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
//cin>>z;
return 0;
}
