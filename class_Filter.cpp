/**************************************************************************** \
*                   Copyright (c) 2006 by Javad GNSS.                        *
*                          All Rights Reserved.                              *
******************************************************************************

Data : 01.03.2006

Author : Glazov Vasiliy

Contents: implementation of classes for filter modeling

Description: these classes contain interface data and methods for filter modeling

\****************************************************************************/
// class_Filter.cpp : class for FIR filter model.

#include <fstream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <math.h>
#include "class_Filter.h"
#define complex_aj_filter 0

using namespace std;

//! Expand sign bit of number if it bit width was reduced
int Filter::keep_sign(int num, int width)
{
    for (int i = width; i < 32; i++){ num = (num | (((num >> (i-1)) & 0x1) << i)); }
    return num;
}

//! Class Filter contain functions for bit-true filtering and processing
Filter::Filter()
{
    vivod1.open("data_out/f1.txt",ios::out);
    vivod2.open("data_out/f2.txt",ios::out);
    vivod3.open("data_out/f3.txt",ios::out);

    test_out_dbg_1.open("data_out/test_out_dbg_1.txt",ios::out);
    test_out_dbg_2.open("data_out/test_out_dbg_2.txt",ios::out);
    test_out_dbg_3.open("data_out/test_out_dbg_3.txt",ios::out);
    test_out_dbg_4.open("data_out/test_out_dbg_4.txt",ios::out);
    test_out_dbg_5.open("data_out/test_out_dbg_5.txt",ios::out);

    test_out_3.open("data_out/test_out_3.txt",ios::out);
    N1=0;
    Nsat=0;
    i=0;
    p1=0;
    psat=0;
    for(i=0;i<cic_stage_number;i++)
    {
        stage_int[i]=0;
        stage_comb[i]=0;
        for(short j=0;j<cic_stage_max_length;j++) stage_comb_z[i][j]=0;
    }
    for(i=0;i<3;i++)
    {
        sa_stage_int[i]=0;
        sa_stage_comb[i]=0;
        for(short j=0;j<513;j++) sa_stage_comb_z[i][j]=0;
    }
    sa_cic_out=0;
    for(i=0;i<cic_stage_number;i++) cic_delay_matrix[i]=0;
    cic_out=0;
    N_aj=0; aj_delay=0; aj_d_delay=0; eps_d=0; yk_d=0; mju_d=0;aj_in_z=0;
    eps_d_i=0;
    eps_d_q=0;
    eps_i=0;
    eps_q=0;

    sec_agc_sum = 1; sec_agc_sum_pre = 1;
    f_sec_agc_old = 0;
    f_sec_agc=0;f_sec_agc_tr=0;f_sec_agc_tr2=0;
    agcmidout = 0;
    nco_acc = 0;

    ifstream input;
    input.open("data_in/flt_tbl.txt");
    for(i=0;i<4096;i++) input>>flt_tbl[i];
    sas1 = 0;sas2=0;sas3=0;sas4=0;
    nco_reg = 0;
    code_e = 0; code_p = 0; code_l = 0;
    strb_h = 0;
    for(i=0;i<14;i++) geng1[i] = 1;
    nco_out_d1 = 0;
    nco_out_d2 = 0;
    k_nco_locata = 0;
    k_nco_locata2 = 0;
    k_nco_locata_lock = 0;
    for(i=0;i<256;i++)
    {
        signal_filtered_line[i] = 0;
        signal_filtered_line_real[i] = 0;
    }
    f_sec_agc_locata = 0;
    f_sec_agc_locata_tr = 0;
    f_sec_agc_locata_tr2 = 0;
    L1 = 0;
    L2 = 0;
    L3 = 0;
    max_locata = 0;
}

Filter::~Filter()
{
    //delete inputI;
}

//! Noise Generator for LOCATA
short Filter::NoiseGenLocata()
{
    short g1s=0;
    short noise=0;

    g1s=geng1[0]^geng1[5]^geng1[7]^geng1[13];
    rotate(&geng1[0], &geng1[13], &geng1[14]);
    geng1[0]=g1s;
    if(geng1[9]==0) noise = (geng1[9]<<9) |(geng1[8]<<8) |
                             (geng1[7]<<7) | (geng1[6]<<6) | (geng1[5]<<5) | (geng1[4]<<4) | (geng1[3]<<3) | (geng1[2]<<2) |
                             (geng1[1]<<1) | (geng1[0]<<0);
    else noise = (0xFF<<10) | (geng1[9]<<9)  | (geng1[8]<<8) |
        (geng1[7]<<7) | (geng1[6]<<6) | (geng1[5]<<5) | (geng1[4]<<4) | (geng1[3]<<3) | (geng1[2]<<2) |
        (geng1[1]<<1) | (geng1[0]<<0);

    return noise;
}
//! Code NCO
short Filter::NCO_code_gen(int freq)
{
    short out = 0;
    nco_reg += freq;
    if((nco_reg>>32) == 0)
        out = 0;
    else
    {
        nco_reg = (nco_reg<<32)>>32;
        out = 1;
    }
    return out;
}

short Filter::NCO_code_gen_phase_gate(int freq, int thr, short code, int phase_sh)
{
    short out = 0;
    short strob;
    short front;
    short sign;
    short magn;

    nco_reg += freq + phase_sh;
    if((nco_reg>>32) == 0)
        out = 0;
    else
    {
        nco_reg = nco_reg&0xffffffff;
        out = 1;
    }
    if(out==1) code_l = code_p;
    if(out==1) code_p = code_e;
    if(out==1) code_e = (code==1) ? 1 : 0;
    if(int(nco_reg&0xffc00000)>=(-thr<<22) && int(nco_reg&0xffc00000)<(thr<<22)) strob = 1;
    else strob = 0;
    if((nco_reg>>31)==1) front = code_e ^ code_p;
    else                 front = code_p ^ code_l;
    if((nco_reg>>31)==1) sign = (!code_e) & 0x1;
    else                 sign = (!code_p) & 0x1;
    magn = (front&0x1) & (strob&0x1);
    strb_h = sign ? -magn : magn;

    short cycle = (code_p ^ code_e) ;

    if(out==1)
    {
		if(cycle==0) k_nco_locata++;
		else k_nco_locata = 0;
    }
	double signal_filtered = locata_pattern[((nco_reg>>19)&0x1fff)/10+819*k_nco_locata];

    double locata_sign = code_e ? signal_filtered : -signal_filtered;
    locata_out = locata_sign;

    return out;
}

short Filter::NCO(short nco_freq, short init)
{
    if(init!=0)
    {
        short nco_phase, cos_phase, sin_phase, cos_sign, sin_sign, a, b;
        nco_acc += nco_freq;
        //vivod1<<nco_acc<<endl;
        //vivod3<<(16384 - nco_acc)<<" "<<nco_acc<<endl;
        if((16384 - nco_acc) <= 0) nco_acc = nco_acc - 16384;
        //nco_phase  = ((nco_acc >> 11) & 0x1) | ((nco_acc >> 10) & 0x1) | ((nco_acc >> 9) & 0x1) | ((nco_acc >> 8) & 0x1)
        //		 | ((nco_acc >> 7) & 0x1) | ((nco_acc >> 6) & 0x1) | ((nco_acc >> 5) & 0x1) | ((nco_acc >> 4) & 0x1)
        //		 | ((nco_acc >> 3) & 0x1) | ((nco_acc >> 2) & 0x1) | ((nco_acc >> 1) & 0x1) | ((nco_acc >> 0) & 0x1);
        nco_phase = nco_acc & 0x3fff;
        a = (nco_acc >> 12) & 0x1;
        b = (a | (a << 1) | (a << 2) | (a << 3) | (a << 4) | (a << 5) | (a << 6) | (a << 7) | (a << 8) | (a << 9) | (a << 10) | (a << 11) | (a << 12));
        //cos_phase = nco_phase ^ (nco_acc & 0xfff);
        cos_phase = (nco_phase ^ b) & 0xfff;
        sin_phase = (~cos_phase) & 0xfff;

        cos_out = flt_tbl[cos_phase];
        sin_out = flt_tbl[sin_phase];
        //vivod1<<cos_phase<<" "<<sin_phase<<endl;
        //vivod2<<nco_acc<<" "<<nco_phase<<endl;

        cos_sign = ((nco_acc >> 13) & 0x1) ^ ((nco_acc >> 12) & 0x1);;
        sin_sign = ((nco_acc >> 13) & 0x1);
        if(cos_sign == 1) cos_out = -cos_out;
        if(sin_sign == 1) sin_out = -sin_out;
        //vivod3<<cos_out<<" "<<sin_out<<cos_sign<<" "<<sin_sign<<" "<<endl;
    }
    else
    {
        nco_acc = 0;
    }
    return cos_out;
}

short Filter::NCO_real_cos(short freq_code, long long t)
{
    nco_acc += freq_code;
    //if (((cos(phase) * sin(phase)) >= 0) & ((cos(phase + 2.0 * pi / 1024.0)  * sin(phase + 2.0 * pi / 1024.0)) >= 0)) {phase += 2.0 * pi / 1024.0;}
    //vivod2<<short(floor(cos(double(nco_acc>>20))*127 + 0.5))<<endl;
    return short(floor(cos(double(nco_acc>>20))*127 + 0.5));
}
void Filter::FilterReset(short s, short cic_del_m[])
{
    sw1=s;
    for(i=0;i<cic_stage_number;i++) cic_delay_matrix[i]=cic_del_m[i];
}
//! Load coefficients of FIR-filter and some parameters
//void Filter::LoadCoeff(short k,short l1,short max_coeff, char *f1)
void Filter::LoadCoeff(short s,short cic_del_m[],
                        short filt_cic_out_shift_pre, short mid_agc_lev_n_pre, short mid_agc_lev_p_pre,
                        short mid_agc_count_n_pre, short mid_agc_count_p_pre, short mid_agc_man_pre, short lms_mju_pre,short h[],short h_2[],short h_pos[])
{
    lms_mju = lms_mju_pre;
    sw1=s;
    for(i=0;i<cic_stage_number;i++) cic_delay_matrix[i]=cic_del_m[i];
    N1=30;//30 20
    N2=10;
    N1_pos = 128;
    h1=new short[N1];
    h1_pos = new short[N1_pos];
    double *h1d;
    h1d=new double[N1];
    h1 = h;
    h1[0] = -5;//-5 1
    h2=new short[N2];
    h2 = h_2;
    h2_dop = 2048;
    h1_pos = h_pos;
    filt_cic_out_shift = filt_cic_out_shift_pre;
    mid_agc_lev_n = mid_agc_lev_n_pre;
    mid_agc_lev_p = mid_agc_lev_p_pre;
    mid_agc_count_n = mid_agc_count_n_pre+1;
    mid_agc_count_p = mid_agc_count_p_pre+1;
    mid_agc_man = mid_agc_man_pre;
    //for(i=0;i<N1;i++)vivod1<<cic_delay_matrix[i]<<endl;
    /*ifstream input;
    input.open(f1);
    //if(input.fail())cout<<"file not open "<<f1;
    // Чтение точных коэффициентов фильтра из файла
    for(i=0;i<N1;i++) input>>h1d[i];
    // Квантование коэффициентов фильтра
    for(i=0;i<N1;i++)
    {
    h1[i]=short(floor(h1d[i]/abs(h1d[N1-max_coeff])*((1<<(k-1))-0)+0.5));
    //vivod1<<h1[i]<<endl;
    //if(h1[i]>=127) h1[i]=127;
    //else if(h1[i]<=-127) h1[i]=-127;
    }
    input.close();
    */
    /*
    for (short i = 0; i < N1; i++)
    {
        cout<<h[i]<<' ';

    }*/
    inputI=new int[N1*2+1];
    input1_rev=new int[N1];
    inputH=new int[N2*2+1];
    input1_revH=new int[N2];
    inputI_pos=new int[N1_pos*2+1];
    input1_rev_pos=new int[N1_pos];
    for(i=0;i<N1*2+1;i++)
    {
        inputI[i]=0;
        input1_rev[i/2]=0;
    }
    for(i=0;i<N2*2+1;i++)
    {
        inputH[i]=0;
        input1_revH[i/2]=0;
    }
    for(i=0;i<N1_pos*2+1;i++)
    {
        inputI_pos[i]=0;
        input1_rev_pos[i/2]=0;
        //vivod1<<h1_pos[i]<<endl;
    }
    /*for(i=0;i<N1_pos;i++)
    {
        h1_pos[i] = h_pos[i]/16;
    }*/
/*
    iira[0] = int(floor( 0.94818115234375   *(int(1)<<22)+0.5));
    iira[1] = int(floor(-2.8444824218750027 *(int(1)<<22)+0.5));
    iira[2] = int(floor( 2.8444824218750044 *(int(1)<<22)+0.5));
    iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb[0] = int(floor(1.0*(1<<22))+0.5);
    iirb[1] = int(floor(-2.8935241699218754 *(int(1)<<22)+0.5));
    iirb[2] = int(floor( 2.7927551269531268 *(int(1)<<22)+0.5));
    iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    //cout<<"iira: "<<iira[0]<<endl;
*/
    /*
    iira[0] = 3977019;
    iira[1] = -11930586;
    iira[2] = 11930586;
    iira[3] = -3977019;
    iirb[0] = floor(1.0*(1<<22))+0.5;
    iirb[1] = -12136303;
    iirb[2] = 11713613;
    iirb[3] = -3770991;
    */
    //nnum = round([ 1.0000592209384946 -2   1.0000592209384951]*2^38);
    //nden = round([ 1 -1.9453434189141721   0.94821382880963578]*2^38);
    //nnum = round([ 2 -2  0]*2^38);
    //nden = round([ 1 -0.94817664023846659  0]*2^38);
    /*
    iirb1[0] = floor( 1.0000592209384946 *(int(1)<<22)+0.5);
    iirb1[1] = floor(-2.0                *(int(1)<<22)+0.5);
    iirb1[2] = floor( 1.0000592209384951 *(int(1)<<22)+0.5);
    iira1[0] = floor(1.0                 *(1<<22))+0.5;
    iira1[1] = floor(-1.9453434189141721 *(int(1)<<22)+0.5);
    iira1[2] = floor( 0.94821382880963578*(int(1)<<22)+0.5);

    iirb2[0] = floor( 2.0                *(int(1)<<22)+0.5);
    iirb2[1] = floor(-2.0                *(int(1)<<22)+0.5);
    iirb2[2] = floor( 0.0                *(int(1)<<22)+0.5);
    iira2[0] = floor(1.0                 *(1<<22))+0.5;
    iira2[1] = floor(-0.94817664023846659*(int(1)<<22)+0.5);
    iira2[2] = floor( 0.0                *(int(1)<<22)+0.5);
    */
    /*!
    iirb1[0] = floor( 1.0000592209384949        *(int(1)<<22)+0.5);
    iirb1[1] = floor(-2.0                *(int(1)<<22)+0.5);
    iirb1[2] = floor( 1.0000592209384951        *(int(1)<<22)+0.5);
    iira1[0] = floor(1.0                 *(1<<22))+0.5;
    iira1[1] = floor(-1.9183673391484659        *(int(1)<<22)+0.5);
    iira1[2] = floor( 0.92442617582521869      *(int(1)<<22)+0.5);

    iirb2[0] = floor( 1.9999999999999998                       *(int(1)<<22)+0.5);
    iirb2[1] = floor(-1.9999999999999998                       *(int(1)<<22)+0.5);
    iirb2[2] = floor( 0.0                *(int(1)<<22)+0.5);
    iira2[0] = floor(1.0                 *(1<<22))+0.5;
    iira2[1] = floor(-0.9243115298054132       *(int(1)<<22)+0.5);
    iira2[2] = floor( 0.0                *(int(1)<<22)+0.5);
    */
/*
    iirb1[0] = int(floor(1.0002696402279485        *(int(1)<<22)+0.5));
    iirb1[1] = int(floor(2.0                *(int(1)<<22)+0.5));
    iirb1[2] = int(floor(1.0002696402279481        *(int(1)<<22)+0.5));
    iira1[0] = int(floor(1.0                 *(1<<22))+0.5);
    iira1[1] = int(floor(1.9130330778968239        *(int(1)<<22)+0.5));
    iira1[2] = int(floor(0.92386597886434962       *(int(1)<<22)+0.5));

    iirb2[0] = int(floor(1.0000462577827243        *(int(1)<<22)+0.5));
    iirb2[1] = int(floor(2*(int(1)<<22)+0.5));
    iirb2[2] = int(floor(1.0000462577827243        *(int(1)<<22)+0.5));
    iira2[0] = int(floor(1.0                 *(1<<22))+0.5);
    iira2[1] = int(floor(1.8156941346856157        *(int(1)<<22)+0.5));
    iira2[2] = int(floor(0.82556805190590121       *(int(1)<<22)+0.5));
*/
    for(short bk = 0; bk<7; bk++)
    {
        iira[bk][0] = 2048*4;
        iirb[bk][0] = 2048*4;
        iira[bk][2] = 2048*4;
    }
    short stepen = 13;
    iira[0][1] = int(floor(-1.496985581573328   *(int(1)<<stepen)+0.5));
    iirb[0][1] = int(floor(-1.3164437465072689  *(int(1)<<stepen)+0.5));
    iirb[0][2] = int(floor(0.92762729803735011  *(int(1)<<stepen)+0.5));
    iira[1][1] = int(floor(-1.5897327128841141  *(int(1)<<stepen)+0.5));
    iirb[1][1] = int(floor(-1.6333592412614473  *(int(1)<<stepen)+0.5));
    iirb[1][2] = int(floor(0.94566083652593935  *(int(1)<<stepen)+0.5));
    iira[2][1] = int(floor(-1.4676325861473418  *(int(1)<<stepen)+0.5));
    iirb[2][1] = int(floor(-1.4165081360362106  *(int(1)<<stepen)+0.5));
    iirb[2][2] = int(floor(0.98606532340181319  *(int(1)<<stepen)+0.5));
    iira[3][1] = int(floor(-1.613104641435855   *(int(1)<<stepen)+0.5));
    iirb[3][1] = int(floor(-1.632957496874897   *(int(1)<<stepen)+0.5));
    iirb[3][2] = int(floor(0.98863929354696367  *(int(1)<<stepen)+0.5));
    iira[4][1] = int(floor(-1.4566419370024952  *(int(1)<<stepen)+0.5));
    iirb[4][1] = int(floor(-1.4380109921935293  *(int(1)<<stepen)+0.5));
    iirb[4][2] = int(floor(0.99724935645233137  *(int(1)<<stepen)+0.5));
    iira[5][1] = int(floor(-1.621277226426741   *(int(1)<<stepen)+0.5));
    iirb[5][1] = int(floor(-1.631243752263059   *(int(1)<<stepen)+0.5));
    iirb[5][2] = int(floor(0.99771142765011933  *(int(1)<<stepen)+0.5));
    iira[6][1] = int(floor(-1.5454171698701393  *(int(1)<<stepen)+0.5));
    iirb[6][1] = int(floor(-1.1448125406843981  *(int(1)<<stepen)+0.5));
    iirb[6][2] = int(floor(0.48155794176998312  *(int(1)<<stepen)+0.5));

    for(short j = 0; j<8; j++)
    for(short i=0;i<4;i++)
    {
        rein[j][i] = 0;
        reind[i] = 0;
        reout[i] = 0;
        reout1[i] = 0;
        reout2[i] = 0;
        max_fast[i] = 0;
        max_fast2[i] = 0;

        reind[i] = 0;
        reoutd[i] = 0;
        reout1d[i] = 0;
        reout2d[i] = 0;
        reout3d[i] = 0;
    }

    max_slow = 0;
    max_slow2 = 0;

    locata_pattern=new double[8192];
    ifstream loc_pattern;
    loc_pattern.open("data_in/locata_pattern.txt");//("data_in/sat10_filt1.fcf");//("data_in/h1_03.fcf");
    if(loc_pattern.fail())cout<<"file not open sat_filt1.fcf";
    for(i=0;i<8192;i++)
    {
        loc_pattern>>locata_pattern[i];
        //cout<<locata_pattern[i]<<endl;
    }
    loc_pattern.close();
}

void Filter::LoadSatCoeff(short lsat, char *f1)
{
    Nsat=lsat;
    hsat=new double[Nsat];
    ifstream inputsat;
    inputsat.open(f1);//("data_in/sat10_filt1.fcf");//("data_in/h1_03.fcf");
    //if(inputsat.fail())cout<<"file not open sat_filt1.fcf";
    for(i=0;i<lsat;i++) inputsat>>hsat[i];
    inputsat.close();
    inputSat=new double[Nsat*2];
    inputSat_rev=new double[Nsat];
    for(i=0;i<Nsat*2;i++)
    {
        inputSat[i]=0;
        inputSat_rev[i/2]=0;
    }
}

double Filter::SatFilter(double signal)
{
    rotate(&inputSat[0], &inputSat[2*Nsat-1], &inputSat[2*Nsat]);
    inputSat[0]=signal;
    reverse_copy(&inputSat[Nsat], &inputSat[2*Nsat], &inputSat_rev[0]);
    transform(&inputSat[0], &inputSat[Nsat], &inputSat_rev[0], &inputSat_rev[0], plus<double>());
    transform(&inputSat_rev[0], &inputSat_rev[Nsat], &hsat[0], &inputSat_rev[0], multiplies<double>());
    psat=accumulate<double *__w64,double>(&inputSat_rev[0],&inputSat_rev[Nsat],0);
    return psat;
}

//! IIR-filter
int Filter::IirI(int re, short k1, short k2)
{
    int iir_out[7] = {0,0,0,0,0,0,0};
    //long long iir_out1 = 0;
    //long long iir_out2 = 0;
    rotate(&rein[0][0], &rein[0][3], &rein[0][4]);
    rein[0][0]=re;

    for(short i = 0;i<7;i++)
    {
        iir_out[i] = (rein[i][0]*iirb[i][0])-(rein[i+1][0]*iira[i][1]);
        iir_out[i] += (rein[i][1]*iirb[i][1])-(rein[i+1][1]*iira[i][2]);
        iir_out[i] += (rein[i][2]*iirb[i][2]);

        rotate(&rein[i+1][0], &rein[i+1][3], &rein[i+1][4]);
        rein[i+1][0]=int(floor(double(iir_out[i])/(1<<14)+0.5));
        vivod1<<iir_out[i]<<' ';
    }
        vivod1<<endl;
/*
    rotate(&reout1[0], &reout1[3], &reout1[4]);
    reout1[0]=int(floor(double(iir_out1)/(int(1)<<22)+0.5));

    iir_out[1] = (reout1[0]*iirb2[0])-(reout2[0]*iira2[1]);
    iir_out[1] += (reout1[1]*iirb2[1])-(reout2[1]*iira2[2]);
    iir_out[1] += (reout1[2]*iirb2[2]);

    rotate(&reout2[0], &reout2[3], &reout2[4]);
    reout2[0]=int(floor(double(iir_out2)/(int(1)<<22)+0.5));
    //if(iir_out%(1<<22)==0 && iir_out<0) reout[0] = reout[0]-1;
*/
    //vivod1<<re<<' '<<rein[7][0]<<endl;

    return int(rein[7][0]);
}

double Filter::IirIn(double re)
{
    double iir_out = 0;
    double iir_out1 = 0;
    double iir_out2 = 0;
    double iir_out3 = 0;
    double iira1d[3],iira2d[3],iira3d[3],iirb1d[3],iirb2d[3],iirb3d[3];
/*
iira1d[0] = 1                        ;
    iira1d[1] =-1.982891822585457                     ;
    iira1d[2] =  0.98297418042823004            ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb1d[0] =  0.017455157327397271                ;
    iirb1d[1] =-0.034880597155728242         ;
    iirb1d[2] =  0.017455157327397271           ;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira2d[0] =   1                        ;
    iira2d[1] =-1.9971311259106153                     ;
    iira2d[2] =   0.9971819987139805                 ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb2d[0] =      0.77058782238070112                        ;
    iirb2d[1] = -1.5410814152050689                 ;
    iirb2d[2] =      0.770587822380701                  ;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira3d[0] =   1                        ;
    iira3d[1] =-1.9904248795668047                      ;
    iira3d[2] =   0.99048773454400563                  ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb3d[0] =      1.0001141013056154                            ;
    iirb3d[1] =-2       ;
    iirb3d[2] =   1.0001141013056156                     ;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
*/
// 50MHz carr
    iira1d[0] = 1;
    iira1d[1] =-1.9042720018419801       ;
    iira1d[2] = 0.91286274399206757      ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb1d[0] =1.000265024764184;
    iirb1d[1] =-2;
    iirb1d[2] = 1.0002650247641838       ;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira2d[0] =   1                        ;
    iira2d[1] =1.7623844749443052        ;
    iira2d[2] =0.80730682166183831       ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb2d[0] =1.0014828231747273        ;
    iirb2d[1] =2;
    iirb2d[2] =1.0014828231747268        ;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira3d[0] =   1                        ;
    iira3d[1] =-0.10905166816419065      ;
    iira3d[2] =-0.73112623229031626      ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb3d[0] =2;
    iirb3d[1] =0;
    iirb3d[2] =-2;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    rotate(&reind[0], &reind[3], &reind[4]);
    reind[0]=re*1              ;
// 50MHz zero
    iira1d[0] = 1;
    iira1d[1] =0.159816466824174       ;
    iira1d[2] =0.039430316894020226      ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb1d[0] =0.16415207142390692;
    iirb1d[1] =0.29952166245226869;
    iirb1d[2] =0.16415207142390684;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira2d[0] =   1                        ;
    iira2d[1] =-0.23674308486882578        ;
    iira2d[2] =0.68843585111042194       ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb2d[0] =0.87961420419836778        ;
    iirb2d[1] =0.38809134697680892;
    iirb2d[2] =0.87961420419836778         ;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira3d[0] =   1                        ;
    iira3d[1] =-0.018620703344422965      ;
    iira3d[2] =0.26763354820599561     ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb3d[0] =0.61528462593029842;
    iirb3d[1] =0.60310521849742349;
    iirb3d[2] =0.61528462593029842;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    //rotate(&reind[0], &reind[3], &reind[4]);
    //reind[0]=re*1              ;

// 500MHz zero
    iira1d[0] = 1;
    iira1d[1] =-1.6500805113598851       ;
    iira1d[2] =0.68511375559292909      ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb1d[0] =0.020645293963506332;
    iirb1d[1] =-0.032699622400934734;
    iirb1d[2] =0.020645293963506329;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira2d[0] =   1                        ;
    iira2d[1] =-1.9358006496806299        ;
    iira2d[2] =0.9546363633013335       ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb2d[0] =1.0168127282314094        ;
    iirb2d[1] =-2;
    iirb2d[2] =1.0168127282314092         ;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira3d[0] =   1                        ;
    iira3d[1] =-1.8162832897033074      ;
    iira3d[2] =0.84052087195989489     ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb3d[0] =1.031603052152922;
    iirb3d[1] =-2;
    iirb3d[2] =1.031603052152922;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));             ;
//150 Mhz
/*
    iira1d[0] = 1;
    iira1d[1] =-1.9714511547752931       ;
    iira1d[2] = 0.97824593318139885      ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb1d[0] = 1.0018651015711693       ;
    iirb1d[1] =-2;
    iirb1d[2] = 1.0018651015711695       ;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira2d[0] =1;
    iira2d[1] =-1.2067260711781589       ;
    iira2d[2] = 0.81921928650471942      ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb2d[0] = 1.9999999999999998       ;
    iirb2d[1] =-1.823465568428378        ;
    iirb2d[2] = 1.9999999999999998       ;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    iira3d[0] =   1                        ;
    iira3d[1] =-1.5347646653599734       ;
    iira3d[2] = 0.59313314904164915      ;
    //iira[3] = int(floor(-0.94818115234375244*(int(1)<<22)+0.5));
    iirb3d[0] =2;
    iirb3d[1] =0;
    iirb3d[2] =-2;
    //iirb[3] = int(floor(-0.899078369140626  *(int(1)<<22)+0.5));
    rotate(&reind[0], &reind[3], &reind[4]);
    reind[0]=re*1              ;
*/
    iir_out1 = (reind[0]*iirb1d[0])-(reout1d[0]*iira1d[1]);
    iir_out1 += (reind[1]*iirb1d[1])-(reout1d[1]*iira1d[2]);
    iir_out1 += (reind[2]*iirb1d[2]);

    rotate(&reout1d[0], &reout1d[3], &reout1d[4]);
    reout1d[0]=iir_out1;

    iir_out2 = (reout1d[0]*iirb2d[0])-(reout2d[0]*iira2d[1]);
    iir_out2 += (reout1d[1]*iirb2d[1])-(reout2d[1]*iira2d[2]);
    iir_out2 += (reout1d[2]*iirb2d[2]);

    rotate(&reout2d[0], &reout2d[3], &reout2d[4]);
    reout2d[0]=iir_out2;

	iir_out3 = (reout2d[0]*iirb3d[0])-(reout3d[0]*iira3d[1]);
    iir_out3 += (reout2d[1]*iirb3d[1])-(reout3d[1]*iira3d[2]);
    iir_out3 += (reout2d[2]*iirb3d[2]);

    rotate(&reout3d[0], &reout3d[3], &reout3d[4]);
    reout3d[0]=iir_out3;
    //if(iir_out%(1<<22)==0 && iir_out<0) reout[0] = reout[0]-1;

    //vivod1<<re<<' '<<reout2[0]<<' '<<iir_out2<<endl;

    return reout3d[0]/2.727;//2.427/0.9998 /1.121        ;2.727 39
}
/*
int Filter::IirI(int re, short k1, short k2)
{
    double iir_out = 0;
    rotate(&rein[0], &rein[3], &rein[4]);
    rein[0]=double(re);


    iir_out = (rein[0]*iira[0])-(reout[0]*iirb[1]);
    iir_out += (rein[1]*iira[1])-(reout[1]*iirb[2]);
    iir_out += (rein[2]*iira[2])-(reout[2]*iirb[3]);
    iir_out += (rein[3]*iira[3]);

    rotate(&reout[0], &reout[3], &reout[4]);
    reout[0]=iir_out;

    vivod1<<re<<' '<<reout[0]<<' '<<iir_out<<endl;

    return rein[3];
}
*/
//! FIR-filter
int Filter::FirI(int re, short k1, short k2)
{
    //if(sw1==0)return re;
    rotate(&inputI[0], &inputI[2*N1-1], &inputI[2*N1]);
    inputI[0]=re;
    reverse_copy(&inputI[N1], &inputI[2*N1], &input1_rev[0]);
    //! Суммирование парных отсчётов сигнала
    transform(&inputI[0], &inputI[N1], &input1_rev[0], &input1_rev[0], plus<int>());
    //! Умножение на коэффициенты фильтра
    transform(&input1_rev[0], &input1_rev[N1], &h1[0], &input1_rev[0], multiplies<int>());
    //! Уменьшение разрядности после перемножения
    //transform(&input1_rev[0], &input1_rev[N1], &input1_rev[0], bind2nd(divides<int>(), (1<<k1)));//
    //! Суммирование результатов перемножения
    p1 = int(accumulate(&input1_rev[0],&input1_rev[N1],0));
    //cout<<((p1>>(k2-1)) & 0x1)<<endl;
    //p1 = (p1>>k2) + ((p1>>(k2-1)) & 0x1);
    ///if (p1 > 32767) p1 = 32767;
    ///if (p1 < -32767) p1 = -32767;
    return p1;
}
int Filter::FirIpos(int re, short k1, short k2)
{
    //if(sw1==0)return re;
    rotate(&inputI_pos[0], &inputI_pos[2*N1_pos], &inputI_pos[2*N1_pos+1]);
    inputI_pos[0]=re;
    reverse_copy(&inputI_pos[N1_pos+1], &inputI_pos[2*N1_pos+1], &input1_rev_pos[0]);
    //! Суммирование парных отсчётов сигнала
    transform(&inputI_pos[0], &inputI_pos[N1_pos], &input1_rev_pos[0], &input1_rev_pos[0], plus<int>());
    //! Умножение на коэффициенты фильтра
    transform(&input1_rev_pos[0], &input1_rev_pos[N1_pos], &h1_pos[0], &input1_rev_pos[0], multiplies<int>());
    //! Уменьшение разрядности после перемножения
    //transform(&input1_rev_pos[0], &input1_rev_pos[N1], &input1_rev_pos[0], bind2nd(divides<int>(), (1<<k1)));//
    //! Суммирование результатов перемножения
    p1 = int(accumulate(&input1_rev_pos[0],&input1_rev_pos[N1_pos],0))+inputI_pos[N1_pos]*8192;
    //vivod1<<int(accumulate(&input1_rev_pos[0],&input1_rev_pos[N1_pos],0))<<' '<<inputI_pos[N1_pos]*511<<endl;
    //cout<<((p1>>(k2-1)) & 0x1)<<endl;
    //p1 = (p1>>k2) + ((p1>>(k2-1)) & 0x1);
    ///if (p1 > 32767) p1 = 32767;
    ///if (p1 < -32767) p1 = -32767;
    return p1;
}

int Filter::FFT(int a_in[], int b_in[], short Nfft)
{
    int a[128],b[128];

    for(short i = 0; i < Nfft; i++)
    {
        a[i] = a_in[i];
        b[i] = b_in[i];
    }


    short Nfft2 = Nfft/2;
    short Nfft4 = Nfft/4;
    short Nfft8 = Nfft/8;
    short Nfft16 = Nfft/16;
    short Nfft32 = Nfft/32;
    short Nfftlog = short(log(double(Nfft))/log(2.0));

    fft_mux_a1_pre = new int[Nfft];
    fft_mux_a2_pre = new int[Nfft];
    fft_mux_a1 = new int[Nfft];
    fft_mux_a2 = new int[Nfft];
    fft_mux_sw = new int*[Nfft];
    fft_mux_sw_rev = new int[Nfft];
    fft_mux_p1 = new int*[Nfft4];
    fft_mux_p2 = new int*[Nfft4];
    for(int i = 0; i < Nfft; i++)
    {
        fft_mux_sw[i] = new int[Nfftlog-2];
        fft_mux_p1[i] = new int[Nfftlog-2];
        fft_mux_p2[i] = new int[Nfftlog-2];
    }
    fft_mux_a1_s = new int[Nfft];
    fft_mux_a2_s = new int[Nfft];

    short *j_seq;
    j_seq = new short[Nfftlog];
    short j_rev;

    double cos_t_0 = 0;
    double sin_t_0 = 0;
    double cos_t = 0;
    double sin_t = 0;

    //! First stage
	//cout<<"zero"<<endl;
    //for(short i=0;i<Nfft;i++)
    //{
    //    cout<<a[i]<<' '<<b[i]<<' '<<i<<endl;
    //}
    for(short i=0;i<Nfft2;i++)
    {
        fft_mux_a1_pre[i] = a[i] + a[i+Nfft2];
        fft_mux_a1_pre[i+Nfft2] = a[i] - a[i+Nfft2];
        fft_mux_a2_pre[i] = b[i] + b[i+Nfft2];
        fft_mux_a2_pre[i+Nfft2] = b[i] - b[i+Nfft2];
    }
    //cout<<"first"<<endl;
    //for(short i=0;i<Nfft4;i++)
    //{
    //    cout<<a[i]<<' '<<b[i]<<' '<<a[i+Nfft2]<<' '<<b[i+Nfft2]<<' '<<fft_mux_a1_pre[i]<<' '<<fft_mux_a2_pre[i]<<' '<<i<<endl;
    //}
    //! Second stage
    for(short i=0;i<Nfft4;i++)
    {
        fft_mux_a1[i] = fft_mux_a1_pre[i] + fft_mux_a1_pre[i+Nfft4];
        fft_mux_a2[i] = fft_mux_a2_pre[i] + fft_mux_a2_pre[i+Nfft4];//0

        fft_mux_a1[i+Nfft4] = fft_mux_a1_pre[i] - fft_mux_a1_pre[i+Nfft4];
        fft_mux_a2[i+Nfft4] = fft_mux_a2_pre[i] - fft_mux_a2_pre[i+Nfft4];//0

        fft_mux_a1[i+Nfft2] = fft_mux_a1_pre[i+Nfft2]+fft_mux_a2_pre[i+Nfft2+Nfft4];
        fft_mux_a2[i+Nfft2] = fft_mux_a2_pre[i+Nfft2]-fft_mux_a1_pre[i+Nfft2+Nfft4];

        fft_mux_a1[i+Nfft2+Nfft4] = fft_mux_a1_pre[i+Nfft2]-fft_mux_a2_pre[i+Nfft2+Nfft4];
        fft_mux_a2[i+Nfft2+Nfft4] = fft_mux_a2_pre[i+Nfft2]+fft_mux_a1_pre[i+Nfft2+Nfft4];
    }
    //cout<<"second"<<endl;
    //for(short i=0;i<Nfft4;i++)
    //{
    //    cout<<fft_mux_a1_pre[i]<<' '<<fft_mux_a1_pre[i+Nfft4]<<' '<<fft_mux_a1[i]<<' '<<fft_mux_a2[i]<<' '<<i<<endl;
    //}
    //! Bit-reverse sequence
    for(short i=0;i<Nfft;i++)
    {
        short j = 0;
        {
            fft_mux_sw[i][j]   = i;
            j_rev = 0;
            for(short n=0;n<Nfftlog;n++) j_seq[n] = (i>>n) & 0x1;
            for(short n=0;n<Nfftlog;n++) j_rev = j_rev | (j_seq[n]<<(Nfftlog-1-n));
            fft_mux_sw_rev[i] = j_rev;
            //cout<<fft_mux_sw[i][j]<<' '<<i<<' '<<j_rev<<' '<<fft_mux_sw_rev[i]<<' '<<j_seq[0]<<' '<<j_seq[1]<<' '<<j_seq[2]<<endl;
        }
    }
    //! Ordinary sequences
    for(short i=0;i<Nfft2;i++)
    {
        short j = 1;
        fft_mux_sw[i][j] = fft_mux_sw[i*2+0][0];
        fft_mux_sw[i+Nfft2][j] = fft_mux_sw[i*2+1][0];
    }
    for(short i=0;i<Nfft4;i++)
    {
        short j = 2;
        fft_mux_sw[i+Nfft4*0][j] = fft_mux_sw[i*4+0][0];
        fft_mux_sw[i+Nfft4*1][j] = fft_mux_sw[i*4+2][0];
        fft_mux_sw[i+Nfft4*2][j] = fft_mux_sw[i*4+1][0];
        fft_mux_sw[i+Nfft4*3][j] = fft_mux_sw[i*4+3][0];
    }
    for(short i=0;i<Nfft8;i++)
    {
        short j = 3;
        fft_mux_sw[i+Nfft8*0][j]    = fft_mux_sw[i*8+0][0];
        fft_mux_sw[i+Nfft8*1][j]    = fft_mux_sw[i*8+4][0];
        fft_mux_sw[i+Nfft8*2][j]    = fft_mux_sw[i*8+2][0];
        fft_mux_sw[i+Nfft8*3][j]    = fft_mux_sw[i*8+6][0];
        fft_mux_sw[i+Nfft8*4][j]    = fft_mux_sw[i*8+1][0];
        fft_mux_sw[i+Nfft8*5][j]    = fft_mux_sw[i*8+5][0];
        fft_mux_sw[i+Nfft8*6][j]    = fft_mux_sw[i*8+3][0];
        fft_mux_sw[i+Nfft8*7][j]    = fft_mux_sw[i*8+7][0];
    }
    for(short i=0;i<Nfft16;i++)
    {
        short j = 4;
        fft_mux_sw[i+Nfft16*0][j]     = fft_mux_sw[i*8+0][0];
        fft_mux_sw[i+Nfft16*1][j]     = fft_mux_sw[i*8+8][0];
        fft_mux_sw[i+Nfft16*2][j]     = fft_mux_sw[i*8+4][0];
        fft_mux_sw[i+Nfft16*3][j]     = fft_mux_sw[i*8+12][0];
        fft_mux_sw[i+Nfft16*4][j]     = fft_mux_sw[i*8+2][0];
        fft_mux_sw[i+Nfft16*5][j]     = fft_mux_sw[i*8+10][0];
        fft_mux_sw[i+Nfft16*6][j]     = fft_mux_sw[i*8+6][0];
        fft_mux_sw[i+Nfft16*7][j]     = fft_mux_sw[i*8+14][0];
        fft_mux_sw[i+Nfft16*8][j]     = fft_mux_sw[i*8+1][0];
        fft_mux_sw[i+Nfft16*9][j]     = fft_mux_sw[i*8+9][0];
        fft_mux_sw[i+Nfft16*10][j]    = fft_mux_sw[i*8+5][0];
        fft_mux_sw[i+Nfft16*11][j]    = fft_mux_sw[i*8+13][0];
        fft_mux_sw[i+Nfft16*12][j]    = fft_mux_sw[i*8+3][0];
        fft_mux_sw[i+Nfft16*13][j]    = fft_mux_sw[i*8+1][0];
        fft_mux_sw[i+Nfft16*14][j]    = fft_mux_sw[i*8+7][0];
        fft_mux_sw[i+Nfft16*15][j]    = fft_mux_sw[i*8+15][0];
    }
    //for(short j=0;j<4;j++){
    //    for(short i=0;i<Nfft;i++){
    //        cout<<fft_mux_sw[i][j]<<' '<<i<<' '<<j<<' '<<endl;}}

    for(short j=0;j<Nfftlog-2;j++)
    {
        for(short i=0;i<Nfft2;i++)
        {
            cos_t_0 =  cos(2*pi/Nfft);
            sin_t_0 = -sin(2*pi/Nfft);
            cos_t = cos_t_0;
            sin_t = sin_t_0;

            for(short n=2;n<=fft_mux_sw_rev[i*2];n++)
            {
                cos_sin_cplx(cos_t,sin_t,cos_t_0,sin_t_0);
                cos_t = Get_cos_sin_cplx(0);
                sin_t = Get_cos_sin_cplx(1);
            }

            if(i==0)
            {
                cos_t = 1;
                sin_t = 0;
            }
            short cos_int = int(floor(cos_t*511 + 0.5));
            short sin_int = int(floor(sin_t*511 + 0.5));

            if(j==0)
            {
                fft_mux_p1[i][j] = cos_int;
                fft_mux_p2[i][j] = sin_int;
            }
            if(j==1 && i<Nfft4)
            {
                for(short n=0;n<2;n++)
                {
                    fft_mux_p1[i+Nfft4*n][j] = cos_int;
                    fft_mux_p2[i+Nfft4*n][j] = sin_int;
                }
            }
            if(j==2 && i<Nfft8)
            {
                for(short n=0;n<4;n++)
                {
                    fft_mux_p1[i+Nfft8*n][j] = cos_int;
                    fft_mux_p2[i+Nfft8*n][j] = sin_int;
                }
            }
            if(j==3 && i<Nfft16)
            {
                for(short n=0;n<8;n++)
                {
                    fft_mux_p1[i+Nfft16*n][j] = cos_int;
                    fft_mux_p2[i+Nfft16*n][j] = sin_int;
                }
            }
            if(j==4 && i<Nfft32)
            {
                for(short n=0;n<16;n++)
                {
                    fft_mux_p1[i+Nfft32*n][j] = cos_int;
                    fft_mux_p2[i+Nfft32*n][j] = sin_int;
                }
            }
        }
    }

    //for(short j=0;j<4;j++){
    //    for(short i=0;i<Nfft2;i++){
    //        vivod1<<fft_mux_p1[i][j]<<' '<<fft_mux_p2[i][j]<<' '<<i<<' '<<j<<' '<<endl;}}
	//for(short i=0;i<Nfft;i++){
	//	vivod3<<"*"<<i<<endl;
	//	vivod2<<"*"<<i+1<<endl;}

    for(short j=Nfftlog-2-1;j>=0;j--)
    {
        if(j<Nfftlog-2-1)
        {
            for(short k=0;k<Nfft;k++)
            {
                fft_mux_a1[k] = fft_mux_a1_s[k];
                fft_mux_a2[k] = fft_mux_a2_s[k];
            }
        }
        short swa;
        short swb;
        //cout<<"A1"<<endl;
        //for(short i=0;i<8;i++) cout<<fft_mux_a1[i]<<' '<<fft_mux_a2[i]<<endl;
        for(short r=0;r<2;r++)
        {
            for(short i=0;i<Nfft4;i++)
            {
                short ir = i + r*Nfft4;
                if(j==0)
                {
                    swa = fft_mux_sw_rev[ir*2];
                    swb = fft_mux_sw_rev[ir*2+1];
                }else
                {
                    swa = fft_mux_sw[ir*2][j];
                    swb = fft_mux_sw[ir*2+1][j];
                }

                int a1 = fft_mux_a1[fft_mux_sw[ir*2][j]];
                int a2 = fft_mux_a2[fft_mux_sw[ir*2][j]];
                int b1 = fft_mux_a1[fft_mux_sw[ir*2+1][j]];
                int b2 = fft_mux_a2[fft_mux_sw[ir*2+1][j]];
                int p1 = fft_mux_p1[ir][j];
                int p2 = fft_mux_p2[ir][j];

                FFTb(a1,a2,b1,b2,p1,p2);
                fft_mux_a1_s[swa] = GetFFTb(0);
                fft_mux_a2_s[swa] = GetFFTb(1);
                fft_mux_a1_s[swb] = GetFFTb(2);
                fft_mux_a2_s[swb] = GetFFTb(3);
				//cout<<a1<<' '<<a2<<' '<<b1<<' '<<b2<<' '<<p1<<' '<<p2<<' '<<j<<' '<<ir<<' '<<swa<<' '<<swb<<' ';
				//cout<<fft_mux_a1_s[swa]<<' '<<fft_mux_a2_s[swa]<<' '<<fft_mux_a1_s[swb]<<' '<<fft_mux_a2_s[swb]<<endl;
            }
        }
    }

    return 1;
}

int Filter::FFTb(int a1, int a2, int b1, int b2, int p1, int p2)
{
    int k1 = b1*(p1+p2);
    int k2 = p2*(b1+b2);
    int k3 = p1*(b2-b1);
    int t1 = k1 - k2;
    int t2 = k1 + k3;

    c1_fft = a1 + int(floor(double(t1)/512+0.5));
    c2_fft = a2 + int(floor(double(t2)/512+0.5));
    d1_fft = a1 - int(floor(double(t1)/512+0.5));
    d2_fft = a2 - int(floor(double(t2)/512+0.5));

    return 1;
}

int Filter::GetFFT(short a, short b)
{
    int c;
    switch(b)
    {
        case 0: c = fft_mux_a1_s[a]; break;
        case 1: c = fft_mux_a2_s[a]; break;
    }
    return c;
}

int Filter::cos_sin_cplx(double a1, double a2, double b1, double b2)
{
    double k1 = a1*(b1+b2);
    double k2 = b2*(a1+a2);
    double k3 = b1*(a2-a1);
    cos_sin_cplx_re = k1 - k2;
    cos_sin_cplx_im = k1 + k3;

    return 1;
}

double Filter::Get_cos_sin_cplx(short a)
{
    double c;
    switch(a)
    {
        case 0: c = cos_sin_cplx_re; break;
        case 1: c = cos_sin_cplx_im; break;
    }
    return c;
}

int Filter::GetFFTb(short a)
{
    int b;
    switch(a)
    {
        case 0: b = c1_fft; break;
        case 1: b = c2_fft; break;
        case 2: b = d1_fft; break;
        case 3: b = d2_fft; break;
    }
    return b;
}

int Filter::FirIn(int re, short k1, short k2)
{
    //if(sw1==0)return re;
    rotate(&inputI[0], &inputI[2*N1_pos-1+1], &inputI[2*N1_pos+1]);
    inputI[0]=re;
    reverse_copy(&inputI[N1_pos+1], &inputI[2*N1_pos+1], &input1_rev[0]);
    //! Суммирование парных отсчётов сигнала
    transform(&inputI[0], &inputI[N1_pos], &input1_rev[0], &input1_rev[0], plus<int>());
    //! Умножение на коэффициенты фильтра
    transform(&input1_rev[0], &input1_rev[N1_pos], &h1[0], &input1_rev[0], multiplies<int>());
    //! Уменьшение разрядности после перемножения
    //transform(&input1_rev[0], &input1_rev[N1_pos], &input1_rev[0], bind2nd(divides<int>(), (1<<k1)));//
    //! Суммирование результатов перемножения
    p1 = int(accumulate(&input1_rev[0],&input1_rev[N1_pos],0)) + inputI[N1_pos]*511;
    //cout<<((p1>>(k2-1)) & 0x1)<<endl;
    //p1 = (p1>>k2) + ((p1>>(k2-1)) & 0x1);
    ///if (p1 > 32767) p1 = 32767;
    ///if (p1 < -32767) p1 = -32767;
    return p1;
}
//! Half-band FIR-filter
int Filter::FirH(int re, short k1, short k2, short sw)
{
    if(sw == 1)
    {
        p1 = inputH[N2-1] * h2_dop;
        p1 = (p1>>k2) + ((p1>>(k2-1)) & 0x1);
    }
    else
    {
        rotate(&inputH[0], &inputH[2*N2-1], &inputH[2*N2]);
        inputH[0]=re;
        reverse_copy(&inputH[N2], &inputH[2*N2], &input1_revH[0]);
        //! Суммирование парных отсчётов сигнала
        transform(&inputH[0], &inputH[N2], &input1_revH[0], &input1_revH[0], plus<int>());
        //! Умножение на коэффициенты фильтра
        transform(&input1_revH[0], &input1_revH[N2], &h2[0], &input1_revH[0], multiplies<int>());
        //! Уменьшение разрядности после перемножения
        //for (short i = 0; i < N2; i++) cout<<h2[i]<<' ';
        //cout<<endl;
        transform(&input1_revH[0], &input1_revH[N2], &input1_revH[0], bind2nd(divides<int>(), (1<<k1)));//
        //! Суммирование результатов перемножения
        p1 = int(accumulate(&input1_revH[0],&input1_revH[N2],0));
        p1 = (p1>>k2) + ((p1>>(k2-1)) & 0x1);
        //if (p1 > 32767) p1 = 32767;
        //if (p1 < -32767) p1 = -32767;
    }
    //vivod1<<re<<' '<<p1<<endl;
    return p1;
}
//! CIC-filter
int Filter::CICFilter(int signal)
{
    //! Comb part
    rotate(&stage_comb_z[0][0], &stage_comb_z[0][cic_delay_matrix[0]], &stage_comb_z[0][cic_delay_matrix[0]+1]);
    stage_comb_z[0][0]=signal;
    stage_comb[0]=signal-stage_comb_z[0][cic_delay_matrix[0]];
    for(i=1;i<cic_stage_number;i++)
    {
        rotate(&stage_comb_z[i][0], &stage_comb_z[i][cic_delay_matrix[i]], &stage_comb_z[i][cic_delay_matrix[i]+1]);
        stage_comb_z[i][0]=stage_comb[i-1];
        stage_comb[i]=stage_comb[i-1]-stage_comb_z[i][cic_delay_matrix[i]];
    }
    //! Int part
    stage_int[0]+=stage_comb[cic_stage_number-1];
    cic_out=stage_int[0];
    for(i=1;i<cic_stage_number;i++)
    {
        stage_int[i]+=cic_out;
        cic_out=stage_int[i];
    }
    return cic_out;
}

short Filter::AGCMid(int cicout, long long t)
{
 if((t-9)%mid_agc_count_n == 0)
 {
  f_sec_agc_old = f_sec_agc;
  if (f_sec_agc_tr<mid_agc_lev_n) f_sec_agc--;
  f_sec_agc_tr = 0;

  //test_out_dbg_1 << sas2 <<" "<<sec_agc_sum_pre<<endl;
  if(sas2-sec_agc_sum_pre>=0){
  if(sas2<sec_agc_sum_pre*2) f_sec_agc = f_sec_agc_old;}
  else if(sec_agc_sum_pre<sas2*2) f_sec_agc = f_sec_agc_old;
  //sec_agc_sum_pre = sec_agc_sum;
  //sas4 = sas3;
  //sas3 = sas2;
  //sas2 = sas1;
  //sas1 = sec_agc_sum;
  //sec_agc_sum = 0;
 }
  //sas4 = sas3;
  sas3 = sas2;
  sas2 = sas1;
  sas1 = sec_agc_sum_pre;
 if((t-7)%mid_agc_count_n == 0)
 {
  sec_agc_sum_pre = sec_agc_sum;
  sec_agc_sum = 0;
 }

 if (f_sec_agc<0) f_sec_agc = 0;
 if ((t-9)%mid_agc_count_p == 0)
 {
  if (f_sec_agc_tr2 > mid_agc_lev_p) f_sec_agc++;
  f_sec_agc_tr2 = 0;
 }
 if ( f_sec_agc > 26) f_sec_agc = 26;
 /*
 if((t+17)%mid_agc_count_n == 0)
 {
  if(sec_agc_sum_pre-sec_agc_sum>=0){
  if(sec_agc_sum_pre/sec_agc_sum<2) f_sec_agc = f_sec_agc_old;}
  else if(sec_agc_sum/sec_agc_sum_pre<2) f_sec_agc = f_sec_agc_old;
  sec_agc_sum_pre = sec_agc_sum;
  sec_agc_sum = 1;
 }
 */
 if(mid_agc_man==1)f_sec_agc = filt_cic_out_shift;
 agcmidout = cicout >> f_sec_agc;

 //f_sec_agc = 13;
 //vivod3<<f_sec_agc<<endl;
 if (((cicout>>31) & 0x1)^((agcmidout>>12) & 0x1)) f_sec_agc_tr++;//15
 if (((cicout>>31) & 0x1)^((agcmidout>>13) & 0x1)) f_sec_agc_tr2++;//16
 agcmidout = cicout >> f_sec_agc;

 //vivod1 << f_sec_agc <<endl;
 //agcmidout = keep_sign(agcmidout, 13);
 agcmidout = (agcmidout>>1)  + (agcmidout & 0x1);
 agcmidout = keep_sign(agcmidout, 13);

 if(agcmidout>4095 || agcmidout<-4096)sas4 = 8192 - abs(agcmidout);
 else sas4 = abs(agcmidout);
 sec_agc_sum += abs(sas4);

 //test_out_dbg_1 << sas4 <<endl;
 //test_out_dbg_2 << sec_agc_sum_pre <<endl;
 //test_out_dbg_3 << sec_agc_sum <<endl;
if(t%200000==0){
	 test_out_dbg_1 << f_sec_agc <<endl;
     //test_out_dbg_3 << sec_agc_sum_pre <<endl;
 }
 return (agcmidout & 0x1FFF);
}
//! Second Middle AGC
int Filter::AGCMids(int cicout_i, int cicout_q,int cicout_id, int cicout_qd, long long t)
{
    int cicout = 0;
    if(abs(cicout_i) > abs(cicout_q)) cicout = abs(cicout_i) + abs(cicout_q)/4;
    else cicout = abs(cicout_q) + abs(cicout_i)/4;

    //vivod1<<f_sec_agc_tr<<' '<<f_sec_agc_tr2<<endl;
    if(cicout>max_slow) max_slow = cicout;
    if(cicout>max_slow2) max_slow2 = cicout;
    f_sec_agc_old = f_sec_agc;
    if((t-9)%(mid_agc_count_n/1) == 0)
    {
        //f_sec_agc_old = f_sec_agc;
        f_sec_agc = int(floor(log(double(max_slow))/log(2.))+2-15);
        if(f_sec_agc<0) f_sec_agc = 0;
        //if((f_sec_agc_old - f_sec_agc)==1) f_sec_agc = f_sec_agc_old;
        //test_out_dbg_4<<f_sec_agc<<' '<<f_sec_agc_old<<' '<<max_slow<<' '<<log(double(max_slow))/log(2.)<<endl;
        max_slow = 0;
    }

    if (f_sec_agc<0) f_sec_agc = 0;
    if ((t-9)%(mid_agc_count_p) == 0)
    {
        f_sec_agc_old = f_sec_agc;
        f_sec_agc = int(floor(log(double(max_slow2))/log(2.))+2-15);
        if(f_sec_agc_old > f_sec_agc) f_sec_agc = f_sec_agc_old;
            //test_out_dbg_4<<f_sec_agc<<' '<<f_sec_agc_old<<' '<<max_slow2<<' '<<log(double(max_slow2))/log(2.)<<endl;
        max_slow2 = 0;
    }
    if ( f_sec_agc > 26) f_sec_agc = 26;
	//f_sec_agc = 6;
    agcmidout = cicout_id >> f_sec_agc;
    agcmidout_q = cicout_qd >> f_sec_agc;

    //cout << f_sec_agc <<endl;
    //vivod1 << cicout_id <<endl;
    //agcmidout = keep_sign(agcmidout, 13);
    agcmidout = (agcmidout>>1)  + (agcmidout & 0x1);
    agcmidout = keep_sign(agcmidout, 14);
    agcmidout_q = (agcmidout_q>>1)  + (agcmidout_q & 0x1);
    agcmidout_q = keep_sign(agcmidout_q, 14);
	//vivod1 << agcmidout <<endl;
    //test_out_dbg_1 << f_sec_agc <<endl;
    //test_out_dbg_2 << sec_agc_sum_pre <<endl;
    //test_out_dbg_3 << sec_agc_sum <<endl;
    if(t%50000==0)
    //if(f_sec_agc_old != f_sec_agc)
    {
        //cout<<f_sec_agc<<' '<<f_sec_agc_old<<' '<<t<<endl;
        test_out_dbg_2<<f_sec_agc<<endl;//' '<<f_sec_agc_old<<' '<<t<<endl;
    }
    //vivod2<<cicout<<endl;
    //vivod3<<f_sec_agc<<endl;
    //if(t>=296*2.5e4 && t< 304*2.5e4) test_out_dbg_4 << (((cicout>>31) & 0x1)^((agcmidout>>15) & 0x1))<<' '<<(((cicout>>31) & 0x1)^((agcmidout>>16) & 0x1)) <<endl;
    //if(t>=180*2.5e4 && t< 184*2.5e4) test_out_dbg_4 << (((cicout>>31) & 0x1)^((agcmidout>>15) & 0x1))<<' '<<(((cicout>>31) & 0x1)^((agcmidout>>16) & 0x1)) <<endl;
    return agcmidout;
}

int Filter::AGCMidsPos(int cicout_i, int cicout_q,int cicout_id, int cicout_qd, long long t)
{
    int cicout = 0;
    if(abs(cicout_i) > abs(cicout_q)) cicout = abs(cicout_i) + abs(cicout_q)/4;
    else cicout = abs(cicout_q) + abs(cicout_i)/4;

    //vivod1<<f_sec_agc_tr<<' '<<f_sec_agc_tr2<<endl;
    if(cicout>max_slow) max_slow = cicout;
    if(cicout>max_slow2) max_slow2 = cicout;
    f_sec_agc_old = f_sec_agc;
    if((t-9)%(mid_agc_count_n/1) == 0)
    {
        //f_sec_agc_old = f_sec_agc;
        f_sec_agc = int(floor(log(double(max_slow))/log(2.))+2-15-0);
        if(f_sec_agc<0) f_sec_agc = 0;
        //if((f_sec_agc_old - f_sec_agc)==1) f_sec_agc = f_sec_agc_old;
        //test_out_dbg_4<<f_sec_agc<<' '<<f_sec_agc_old<<' '<<max_slow<<' '<<log(double(max_slow))/log(2.)<<endl;
        max_slow = 0;
    }

    if (f_sec_agc<0) f_sec_agc = 0;
    if ((t-9)%(mid_agc_count_p) == 0)
    {
        f_sec_agc_old = f_sec_agc;
        f_sec_agc = int(floor(log(double(max_slow2))/log(2.))+2-15-0);
        if(f_sec_agc_old > f_sec_agc) f_sec_agc = f_sec_agc_old;
            //test_out_dbg_4<<f_sec_agc<<' '<<f_sec_agc_old<<' '<<max_slow2<<' '<<log(double(max_slow2))/log(2.)<<endl;
        max_slow2 = 0;
    }
    if ( f_sec_agc > 26) f_sec_agc = 26;

    agcmidout = cicout_id >> f_sec_agc;
    agcmidout_q = cicout_qd >> f_sec_agc;

    //cout << f_sec_agc <<endl;
    //vivod1 << cicout_id <<endl;
    //agcmidout = keep_sign(agcmidout, 13);
    agcmidout = (agcmidout>>1)  + (agcmidout & 0x1);
    agcmidout = keep_sign(agcmidout, 14);
    agcmidout_q = (agcmidout_q>>1)  + (agcmidout_q & 0x1);
    agcmidout_q = keep_sign(agcmidout_q, 14);
	//vivod1 << agcmidout <<endl;
    //test_out_dbg_1 << f_sec_agc <<endl;
    //test_out_dbg_2 << sec_agc_sum_pre <<endl;
    //test_out_dbg_3 << sec_agc_sum <<endl;
    if(t%50000==0)
    //if(f_sec_agc_old != f_sec_agc)
    {
        //cout<<f_sec_agc<<' '<<f_sec_agc_old<<' '<<t<<endl;
        test_out_dbg_2<<f_sec_agc<<endl;//' '<<f_sec_agc_old<<' '<<t<<endl;
    }
    //vivod2<<cicout<<endl;
    //vivod3<<f_sec_agc<<endl;
    //if(t>=296*2.5e4 && t< 304*2.5e4) test_out_dbg_4 << (((cicout>>31) & 0x1)^((agcmidout>>15) & 0x1))<<' '<<(((cicout>>31) & 0x1)^((agcmidout>>16) & 0x1)) <<endl;
    //if(t>=180*2.5e4 && t< 184*2.5e4) test_out_dbg_4 << (((cicout>>31) & 0x1)^((agcmidout>>15) & 0x1))<<' '<<(((cicout>>31) & 0x1)^((agcmidout>>16) & 0x1)) <<endl;
    return agcmidout;
}
//! Middle AGC
short Filter::AGCMid2(int cicout_i, int cicout_q,int cicout_id, int cicout_qd, long long t)
{
    int cicout = 0;
    int abs_cicout_i = abs(cicout_i);
    int abs_cicout_q = abs(cicout_q);
    if(abs_cicout_i > abs_cicout_q) cicout = abs_cicout_i + abs_cicout_q/4;
    else cicout = abs_cicout_q + abs_cicout_i/4;

    max_fast[4] = max_fast[3];
    max_fast[3] = max_fast[2];
    max_fast[2] = max_fast[1];
    max_fast[1] = max_fast[0];
    if(cicout>max_fast[0]) max_fast[0] = cicout;
    max_fast2[4] = max_fast2[3];
    max_fast2[3] = max_fast2[2];
    max_fast2[2] = max_fast2[1];
    max_fast2[1] = max_fast2[0];
    if(cicout>max_fast2[0]) max_fast2[0] = cicout;

    //vivod1<<cicout_i<<' '<<cicout<<' '<<agcmidout<<' '<<f_sec_agc<<' '<<((t-7)%mid_agc_count_p == 0)<<' '<<((t-7)%mid_agc_count_n == 0)<<' ';
    //vivod1<<int(floor(log(double(max_fast[4]))/log(2.))+2-14)<<' '<<int(floor(log(double(max_fast2[4]))/log(2.))+2-14)<<' '<<max_fast[4]<<' '<<max_fast2[4]<<endl;

    f_sec_agc_old = f_sec_agc;
    if((t-7)%mid_agc_count_n == 0)
    {
        f_sec_agc_old = f_sec_agc;
        f_sec_agc = int(floor(log(double(max_fast[4]))/log(2.))+2-14);
        if(f_sec_agc<0) f_sec_agc = 0;

        //test_out_dbg_5<<f_sec_agc<<' '<<f_sec_agc_old<<' '<<max_fast<<' '<<log(double(max_fast))/log(2.)<<endl;
        if((f_sec_agc_old - f_sec_agc)==1) f_sec_agc = f_sec_agc_old;
        max_fast[0] = cicout;
    }

    if (f_sec_agc<0) f_sec_agc = 0;
    if ((t-7)%mid_agc_count_p == 0)
    {

        f_sec_agc_old = f_sec_agc;
        f_sec_agc = int(floor(log(double(max_fast2[4]))/log(2.))+2-14);
        if(f_sec_agc_old > f_sec_agc) f_sec_agc = f_sec_agc_old;
        max_fast2[0] = cicout;
    }
    if ( f_sec_agc > 26) f_sec_agc = 26;

    if(mid_agc_man==1)f_sec_agc = filt_cic_out_shift;
	//f_sec_agc = 12;
    agcmidout = cicout >> f_sec_agc;

    agcmidout = cicout_id >> f_sec_agc;
    agcmidout = (agcmidout>>1)  + (agcmidout & 0x1);
    if(agcmidout >4095) cout << "PROBLEM" <<endl;
    agcmidout = keep_sign(agcmidout, 13);

    agcmidout_q = cicout_qd >> f_sec_agc;
    agcmidout_q = (agcmidout_q>>1)  + (agcmidout_q & 0x1);
    agcmidout_q = keep_sign(agcmidout_q, 13);


    //test_out_dbg_1 << sas4 <<endl;
    //test_out_dbg_2 << sec_agc_sum_pre <<endl;
    //test_out_dbg_3 << sec_agc_sum <<endl;
    if(t%200000==0)
    {
        test_out_dbg_1 << f_sec_agc <<endl;
        //test_out_dbg_3 << sec_agc_sum_pre <<endl;
    }
    return agcmidout;
}
//! Loading parameters for anti-jaming filter
void Filter::LoadAJFilter(short N, short adapt_delay, short mult_rnd, short input_b)
{
    N_aj=N;
    w_aj=new short[N_aj];
	w_aj_abs=new short[256];
    w_ajj=new short[N_aj];
    w_aj1=new short[N_aj];
    w_aj2=new short[N_aj];
    w_aj11=new short[N_aj];
    w_aj22=new short[N_aj];
    aj_delay=1*adapt_delay;//30+0*N/2;
    aj_d_delay=0;
    inputAJ=new int[N_aj+aj_delay+aj_d_delay+3];
    inputAJ2=new int[N_aj+aj_delay+aj_d_delay+3];
    inputAJ_rev=new int[N_aj];

    waj_res=new int[N_aj];
    waj_res1=new int[N_aj];
    waj_res2=new int[N_aj];
    w_lms=new int[N_aj];
    w_lms1=new int[N_aj];
    w_lms2=new int[N_aj];

    input_delay_d=new double[aj_d_delay];
    w_aj_d=new double[N_aj];
    inputAJ_d=new double[N_aj+aj_delay+aj_d_delay+3];
    inputAJ_rev_d=new double[N_aj];

    w_aj_d_cplx_i=new double[N_aj];
    inputAJ_d_cplx_i=new double[N_aj+aj_delay];
    inputAJ_rev_d_cplx_i=new double[N_aj];
    w_aj_d_cplx_q=new double[N_aj];
    inputAJ_d_cplx_q=new double[N_aj+aj_delay];
    inputAJ_rev_d_cplx_q=new double[N_aj];

    for(i=0;i<N_aj;i++)
    {
        w_aj[i] = 0;
        w_ajj[i] = 0;
        w_aj1[i] = 0;
        w_aj2[i] = 0;
        w_aj11[i] = 0;
        w_aj22[i] = 0;
        inputAJ_rev[i] = 0;

        waj_res[i] =  0;
        waj_res1[i] =  0;
        waj_res2[i] =  0;
        w_lms[i] = 0;
        w_lms1[i] = 0;
        w_lms2[i] = 0;

        w_aj_d[i]=0;
        w_aj_d_cplx_i[i]=0;
        w_aj_d_cplx_q[i]=0;
        inputAJ_rev_d[i]=0;
        inputAJ_rev_d_cplx_i[i]=0;
        inputAJ_rev_d_cplx_q[i]=0;
    }
    /*w_aj[0]=-3219;
    w_aj[1]=-2994;
    w_aj[2]=-2706;
    w_aj[3]=-2429;
    w_aj[4]=-2119;
    w_aj[5]=-1789;
    w_aj[6]=-1441;
    w_aj[7]=-1067;
    w_aj[8]=-685;
    w_aj[9]=-307;
    w_aj[10]=55;
    w_aj[11]=436;
    w_aj[12]=817;
    w_aj[13]=1180;
    w_aj[14]=1508;
    w_aj[15]=1832;
    //w_aj[16]=1426;
    //w_aj[17]=266;
    //w_aj[18]=-1125;
    //w_aj[19]=-1513;Columns 1 through 7
    */
    yk=0;
    yk1=0;
    yk2=0;
    yk3=0;
    yk4=0;
    ykk=0;
    ykk1=0;
    ykk2=0;
    ykk3=0;
    ykk4=0;
    eps1=0;
    eps2=0;
    eps3=0;
    eps4=0;
    eps5=0;
    for(i=0;i<(N_aj+aj_delay+aj_d_delay+3);i++)
    {
        inputAJ[i]=0;
        inputAJ2[i]=0;
        inputAJ_d[i]=0;
    }
    mju=1<<21;//131072
    eps=0;
    epss=0;

    for(i=0;i<(N_aj+aj_delay);i++)
    {
        inputAJ_d_cplx_i[i]=0;
        inputAJ_d_cplx_q[i]=0;
    }
    for(i=0;i<(aj_d_delay);i++) input_delay_d[i]=0;
    mju_d=1.0/(1048576.0*1.0);//1.0/(1000000.0);
    eps_d=0;
    yk_d=0;
    aj_in_z=0;
    yk_d_cplx_i=0;
    yk_d_cplx_q=0;
    eps_d_cplx_i=0;
    eps_d_cplx_q=0;

    // r1_i=0;
    //r1_q=0;
    input_bw = input_b+2; //8
    fltcoeff_bw = 11; //10
    multround = mult_rnd-1+4+1   -2  -2 ; //10
   // cout<<"dsfasdfsadfsadfdsaf"<<multround<<endl;
    multround_bw = (input_bw + fltcoeff_bw - 1 - multround);
    sumround = 2;
    sumround_bw = (multround_bw + 6); //multround_bw + log2(N_aj)
    coeffcompmult_bw = (sumround_bw + input_bw + 0 - 1);
    coeffcompldel = input_bw*2+2-0   +lms_mju-5; // -lms_mju+1;     // -2;//14
    coeffcompdiv_bw = coeffcompmult_bw - coeffcompldel;

    /*
    input_bw = input_b+2; //8
    fltcoeff_bw = 11; //10
    multround = mult_rnd-1+4+1   -2  -2; //10
    multround_bw = (input_bw + fltcoeff_bw - 1 - multround);
    sumround = 2;
    sumround_bw = (multround_bw + 6); //multround_bw + log2(N_aj)
    coeffcompmult_bw = (sumround_bw + input_bw + 0 - 1);
    coeffcompldel = input_bw*2+2-0   +lms_mju-5; // -lms_mju+1;     // -2;//14
    coeffcompdiv_bw = coeffcompmult_bw - coeffcompldel;
    */

    ave_input1_z = 0;
    ave_input2_z = 0;
}

short Filter::AJFilter(short aj_in,short aj_in2,short k,short compute_coeff_en, short log_en)
{
#if complex_aj_filter==1
 for(i=(N_aj+aj_delay-1);i>0;i--)
 {
  inputAJ[i]=inputAJ[i-1];
  inputAJ2[i]=inputAJ2[i-1];
 }
 inputAJ[0]=aj_in;
 inputAJ2[0]=aj_in2;
 yk=0;
 yk2=0;
 for(i=0;i<N_aj;i++)
 {
  yk+=int(floor(double(inputAJ[aj_delay+i]*w_aj[i]-inputAJ2[aj_delay+i]*w_aj2[i])/16384.0+0.5));
  yk2+=int(floor(double(inputAJ[aj_delay+i]*w_aj2[i]+inputAJ2[aj_delay+i]*w_aj[i])/16384.0+0.5));
 }
 eps=inputAJ[0]-yk;
 eps2=inputAJ2[0]-yk2;
 //if(compute_coeff_en==1) {copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<__int16>(vivod1, " ")); vivod1<<endl;}
 //vivod1<<aj_in<<endl;
 //vivod2<<yk<<endl;
 //vivod3<<yk<<endl;
 //if(compute_coeff_en==1)
 for(i=0;i<N_aj;i++)
 {
	 w_aj[i]=w_aj[i]+int(floor(double(2*(eps*inputAJ[aj_delay+i]+eps2*inputAJ2[aj_delay+i]))/double(mju)+0.5));//-int(floor(double(w_aj[i])/double(16536)+0.5));
	 w_aj2[i]=w_aj2[i]+int(floor(double(2*(-eps*inputAJ2[aj_delay+i]+eps2*inputAJ[aj_delay+i]))/double(mju)+0.5));
	 //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<__int16>(vivod1, " ")); vivod1<<endl;
 }

#else
 //Delay line
 for(i = (N_aj+aj_delay+aj_d_delay+2); i>0; i--) inputAJ[i] = inputAJ[i-1];//2 (N_aj+aj_delay+aj_d_delay-1+3)
 //inputAJ[0] = ((aj_in<<16)>>16);
 inputAJ[0] = (aj_in<<19)>>19;//2 3 inputAJ[0] = ((aj_in/(1<<0))<<19)>>19;
 //vivod1<<inputAJ[aj_delay+i]<<endl;
 //Filter multipliers
 yk = 0;
 //{copy(&inputAJ[aj_delay+0], &inputAJ[aj_delay+N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;}
 //{copy(&w_aj2[0], &w_aj2[N_aj], ostream_iterator<double>(vivod3, " ")); vivod3<<endl;}
 //vivod3<<inputAJ[aj_delay+0]<<"   "<<w_aj2[0]<<endl;
 for(i = 0; i<N_aj; i++) {yk += int(floor(double((inputAJ[aj_delay+i]*w_aj2[i]))/double(1024)+0.5));//(16384/4))<<21)>>21 - 32 mult///w_aj1[i]
 //vivod2<<yk<<"   "<<aj_in<<endl;
 //vivod2<<((((inputAJ[aj_delay+0]*w_aj2[0])/(16384/2))<<0)>>0)<<endl;
 //vivod2<<int(floor(double((inputAJ[aj_delay+i]*w_aj2[i]))/double(16384/8)+0.5))<<"  ";
 }
 //vivod2<<inputAJ[0]<<endl;
 //vivod3<<yk3<<endl;
 yk = (yk<<19)>>19;
 eps=inputAJ[0]*8-yk3;//yk2
 yk4=yk3;
 yk3=yk2;
 yk2=yk1;
 yk1=yk;
 eps=(eps<<19)>>19;
 //if(compute_coeff_en==1)
 //{
  //vivod2<<yk4<<endl;
  //vivod3<<inputAJ[0]*8<<endl;
 //}
/*
 for(i=0;i<N_aj;i++)
 {
  w_aj[i] += eps*inputAJ[aj_delay+i]/(1<<18);
  w_aj[i] = ((w_aj[i]<<20)>>20);
 }
 */
 //LMS compute coefficients
 //{copy(&inputAJ[aj_delay+3], &inputAJ[aj_delay+N_aj+3], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;}
 //vivod1<<eps1<<'   '<<inputAJ[aj_delay+2]<<endl;
 for(i=0;i<N_aj;i++) w_aj2[i] = w_aj1[i];
 if(compute_coeff_en==1){
 int ldel = 1<<14;//16 15
 //{copy(&inputAJ[aj_delay+3], &inputAJ[aj_delay+N_aj+3], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;}
 //vivod1<<eps<<"   "<<inputAJ[aj_delay+3]<<endl;
 for(i=0;i<N_aj;i++)
 {
 //w_aj2[i] = w_aj1[i];
 w_aj1[i] = w_aj[i];

  w_lms[i] = ((eps*(inputAJ[aj_delay+i+3]*8))<<5)>>5;//eps1


  //Накопление остатка целочисленного деления
  /*w_aj[i] += (((w_lms[i]+waj_res[i])/(ldel))<<19)>>19;
  waj_res[i] = w_aj[i]%(ldel);
  */
  //waj_res[i] += w_lms[i]%(ldel);//Накопление остатка целочисленного деления
  w_aj[i] += ((w_lms[i]/(ldel))<<19)>>19;

  //учёт остатка от целочисленного деления
  if(waj_res[i]>(ldel) || waj_res[i]<-(ldel))
  {
   w_aj[i] += waj_res[i]/(ldel);
   waj_res[i] = w_lms[i]%(ldel) + waj_res[i]%(ldel);
  }//-=(sdvig/4);}
  else waj_res[i] += w_lms[i]%(ldel);

  //недопускание переполнения коэффициентов МНК
  if(w_aj[i]>511 || w_aj[i]<-512)
  {
   //vivod1<<i<<endl;
   if(w_aj[i]>0) w_aj[i] = 511;
   else w_aj[i] = -511;
  }
  w_aj[i] = ((w_aj[i]<<22)>>22);
 }
 }
 eps5=eps4;
 eps4=eps3;
 eps3=eps2;
 eps2=eps1;
 eps1=eps;
 //if(log_en==1)
 //	{copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;}
#endif
 //vivod2<<eps<<endl;
 //if(compute_coeff_en==1)
 //{
  //vivod1<<eps<<endl;
  //vivod3<<inputAJ[0]<<endl;
  //vivod3<<((aj_in<<16)>>16)<<" "<<inputAJ[0]<<endl;
 //}
 return eps;
}

short Filter::AJFilterNew(short aj_in,short aj_in2,short k,short compute_coeff_en, short log_en)
{
 //Delay line
 for(i = (N_aj+aj_delay+aj_d_delay+2); i>0; i--) inputAJ[i] = inputAJ[i-1];//2 (N_aj+aj_delay+aj_d_delay-1+3)
 inputAJ[0] = (aj_in<<(32 - input_bw))>>(32 - input_bw);//2 3 inputAJ[0] = ((aj_in/(1<<0))<<19)>>19;
 yk = 0;
 //vivod1<<inputAJ[0]<<endl;
 //Multipliers and accumulator
 for(i = 0; i<N_aj; i++) {yk += int(floor(double(inputAJ[aj_delay+i]*w_aj2[i])/double(1<<multround)+0.5));//(16384/4))<<21)>>21 - 32 mult///w_aj1[i]
 }
 //vivod2<<yk<<endl;
 //yk /= (1<<sumround);
 //vivod2<<yk<<endl;
 yk = (yk<<(32 - sumround_bw))>>(32 - sumround_bw);
 //Difference
 eps=inputAJ[0]*16-yk3;//8
 yk3=yk2;
 yk2=yk1;
 yk1=yk;
 //vivod1<<inputAJ[0]*8<<endl;
 eps=(eps<<(32 - sumround_bw))>>(32 - sumround_bw);
 //vivod2<<yk3<<endl;
 //vivod3<<eps<<endl;
 //Coefficient computing
 for(i=0;i<N_aj;i++) w_aj2[i] = w_aj1[i];
 if(compute_coeff_en==1)
 {
	int ldel = 1<<coeffcompldel;//16 15
	for(i=0;i<N_aj;i++)
	{
	w_aj1[i] = w_aj[i];
	//w_lms[i] = ((eps*(inputAJ[aj_delay+i+3]*16))<<(32 - coeffcompmult_bw))>>(32 - coeffcompmult_bw);//eps1
	w_lms[i] = eps*(inputAJ[aj_delay+i+3]*1);//eps1
	//vivod3<<(eps*(inputAJ[aj_delay+i+3]*8))<<endl;
	//vivod2<<inputAJ[aj_delay+i+3]<<endl;
	//Накопление остатка целочисленного деления
	 w_aj[i] += ((w_lms[i]/(ldel))<<(32 - coeffcompdiv_bw))>>(32 - coeffcompdiv_bw);
	 //if(w_lms[i]/(ldel)>0) w_aj[i] +=1;
	 //else if(w_lms[i]/(ldel)<0) w_aj[i] -=1;
	 //учёт остатка от целочисленного деления
	 if(waj_res[i]>(ldel) || waj_res[i]<-(ldel))
	 {
	  w_aj[i] += waj_res[i]/(ldel);
	  waj_res[i] = w_lms[i]%(ldel) + waj_res[i]%(ldel);
	 }
	 else waj_res[i] += w_lms[i]%(ldel);

	 //недопускание переполнения коэффициентов МНК
	 if(w_aj[i]>((1<<(fltcoeff_bw-1)) - 1) || w_aj[i]<-(1<<(fltcoeff_bw-1)))
	 {
	  if(w_aj[i]>0) w_aj[i] = ((1<<(fltcoeff_bw-1)) - 1);
	  else w_aj[i] = -((1<<(fltcoeff_bw-1)) - 1);
	 }
	 w_aj[i] = ((w_aj[i]<<(32 - fltcoeff_bw))>>(32 - fltcoeff_bw));

	}

	//vivod2<<(w_lms[0]/(ldel))<<endl;
	//vivod3<<waj_res[0]/(ldel)<<endl;
 }
 eps5=eps4;
 eps4=eps3;
 eps3=eps2;
 eps2=eps1;
 eps1=eps;

 if(log_en==1)
 {
	//copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
	//copy(&waj_res[0], &waj_res[N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;
 }
 //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
 //vivod1<<w_aj[0]<<endl;
 //vivod2<<waj_res[0]<<endl;
 //vivod3<<waj_res[0]/(1<<coeffcompldel)<<endl;
 //copy(&waj_res[0], &waj_res[N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;
 //vivod3<<eps<<endl;
 //if(compute_coeff_en==1)
 //{
  //vivod1<<eps<<endl;
  //vivod3<<inputAJ[0]<<endl;
  //vivod3<<((aj_in<<16)>>16)<<" "<<inputAJ[0]<<endl;
 //}
 return eps;
}
//! Anti-Jamming filter
int Filter::AJFilterNew2(int aj_in,short k, short log_en,int point[])
{
    //! Delay line
    for(i = (N_aj+aj_delay+aj_d_delay+2); i>0; i--) inputAJ[i] = inputAJ[i-1];//2 (N_aj+aj_delay+aj_d_delay-1+3)
    inputAJ[0] = aj_in;//2 3 inputAJ[0] = ((aj_in/(1<<0))<<19)>>19;
    yk = 0;
    yk_d = 0;
    //! Multipliers and accumulator
    for(i = 0; i<N_aj; i++)
    {
        yk += int(floor(double(inputAJ[aj_delay+i]*w_aj1[i])/double(1<<multround)+0.5));
    }
    yk = int(floor(double(yk)/double(1<<(sumround))+0.5));
    //Difference
    eps=inputAJ[0]*1-yk3;//8
    //vivod2<<inputAJ[0]<<' '<<yk3<<' '<<eps<<endl;
    if(log_en==1)
    {
        //vivod2<<inputAJ[0]<<' '<<yk3<<endl;
    }
    yk3=yk2;
    yk2=yk1;
    yk1=yk;

    //! Coefficient computing
    for(i=0;i<N_aj;i++) w_aj1[i] = w_aj[i];
    int ldel = 1<<(coeffcompldel);//16 15
    for(i=0;i<N_aj;i++)
    {
        w_lms[i] = eps*(inputAJ[aj_delay+i+3]*1);//eps1  +3
        //Накопление остатка целочисленного деления
        w_aj[i] += w_lms[i]/ldel;
        //учёт остатка от целочисленного деления
        if(waj_res[i]>(ldel-1) || waj_res[i]<-(ldel))
        {
            w_aj[i] += waj_res[i]/(ldel);
            waj_res[i] = (w_lms[i]%(ldel)) + waj_res[i]%(ldel);
        }
        else waj_res[i] = waj_res[i] + w_lms[i]%(ldel);
        //недопускание переполнения коэффициентов МНК
        if(w_aj[i]>((1<<(fltcoeff_bw-1)) - 1) || w_aj[i]<-(1<<(fltcoeff_bw-1)))
        {
            if(w_aj[i]>0) w_aj[i] = ((1<<(fltcoeff_bw-1)) - 1);
            else w_aj[i] = -((1<<(fltcoeff_bw-1)) - 1);
        }
    }

    if(log_en==1)
    {
        //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod3, " ")); vivod3<<endl;
        //vivod3<<w_lms[0]<<' '<<eps<<' '<<(inputAJ[aj_delay+0+3]*1)<<' '<<w_lms[0]/(1<<coeffcompldel)<<' '<<w_aj[0]<<endl;
        //copy(&waj_res[0], &waj_res[N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;

		//copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
        //vivod2<<aj_in<<endl;
    }

	/*int aj_coeff_sum = 0;
	for(short k=0;k<255;k++)
	{
		w_aj_abs[k] = abs(inputAJ[k]);//abs(inputAJ[k]);
		aj_coeff_sum += w_aj_abs[k];
	}
	short w_aj_abs_max = *max_element(&w_aj_abs[0],&w_aj_abs[256]);
	if(w_aj_abs_max!=0)
	{
		//vivod2<<(double(aj_coeff_sum/80)/double(w_aj_abs_max))<<' '<<aj_coeff_sum/80<<' '<<w_aj_abs_max<<endl;
		//if(log_en==0)
		//	vivod2<<(double(aj_coeff_sum/256)/double(w_aj_abs_max))<<endl;
		//else
		//	vivod3<<(double(aj_coeff_sum/256)/double(w_aj_abs_max))<<endl;
		//vivod1<<inputAJ[0]<<endl;
		//copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
		//if((aj_coeff_sum/256/double(w_aj_abs_max)) > 0.41)
		//{vivod2<<'1'<<endl; cout<<'1'<<endl;}
		//else
		//	vivod2<<'0'<<endl;
	}
    */
    //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
    /*
    w_aj[0] = 87;
    w_aj[1] = -4;
    w_aj[2] = -87;
    w_aj[3] = 6;
    w_aj[4] = 86;
    w_aj[5] = -7;
    w_aj[6] = -86;
    w_aj[7] = 8;
    w_aj[8] = 86;
    w_aj[9] = -9;
    w_aj[10] = -85;
    w_aj[11] = 10;
    w_aj[12] = 86;
    w_aj[13] = -11;
    w_aj[14] = -86;
    w_aj[15] = 12;
    w_aj[16] = 86;
    w_aj[17] = -13;
    w_aj[18] = -86;
    w_aj[19] = 14;
    w_aj[20] = 87;
    w_aj[21] = -16;
    w_aj[22] = -87;
    w_aj[23] = 17;
    w_aj[24] = 88;
    w_aj[25] = -18;
    w_aj[26] = -88;
    w_aj[27] = 20;
    w_aj[28] = 89;
    w_aj[29] = -21;
    w_aj[30] = -90;
    w_aj[31] = 22;
    w_aj[32] = 91;
    w_aj[33] = -24;
    w_aj[34] = -92;
    w_aj[35] = 25;
    w_aj[36] = 92;
    w_aj[37] = -27;
    w_aj[38] = -94;
    w_aj[39] = 29;
    w_aj[40] = 94;
    w_aj[41] = -30;
    w_aj[42] = -95;
    w_aj[43] = 32;
    w_aj[44] = 97;
    w_aj[45] = -33;
    w_aj[46] = -97;
    w_aj[47] = 35;
    w_aj[48] = 99;
    w_aj[49] = -37;
    w_aj[50] = -100;
    w_aj[51] = 39;
    w_aj[52] = 101;
    w_aj[53] = -41;
    w_aj[54] = -102;
    w_aj[55] = 43;
    w_aj[56] = 102;
    w_aj[57] = -45;
    w_aj[58] = -104;
    w_aj[59] = 47;
    w_aj[60] = 105;
    w_aj[61] = -49;
    w_aj[62] = -105;
    w_aj[63] = 51;
    w_aj[64] = 106;
    w_aj[65] = -53;
    w_aj[66] = -107;
    w_aj[67] = 55;
    w_aj[68] = 108;
    w_aj[69] = -57;
    w_aj[70] = -109;
    w_aj[71] = 60;
    w_aj[72] = 109;
    w_aj[73] = -61;
    w_aj[74] = -110;
    w_aj[75] = 64;
    w_aj[76] = 111;
    w_aj[77] = -66;
    w_aj[78] = -111;
    w_aj[79] = 68;
    */
    //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
    //vivod1<<w_aj[0]<<endl;
    //vivod2<<waj_res[0]<<endl;
    //vivod3<<waj_res[0]/(1<<coeffcompldel)<<endl;
    //copy(&waj_res[0], &waj_res[N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;
    //vivod3<<eps<<endl;
    //if(compute_coeff_en==1)
    //{
    //vivod1<<eps<<endl;
    //vivod3<<inputAJ[0]<<endl;
    //vivod3<<((aj_in<<16)>>16)<<" "<<inputAJ[0]<<endl;
    //}
    return eps;
}

int Filter::AJFilterNew3(int aj_in,int aj_in2,short k, short log_en,int point[])
{
    //! Delay line
    for(i = (N_aj+aj_delay+aj_d_delay+2); i>0; i--) inputAJ[i] = inputAJ[i-1];//2 (N_aj+aj_delay+aj_d_delay-1+3)
    inputAJ[0] = aj_in;//2 3 inputAJ[0] = ((aj_in/(1<<0))<<19)>>19;
    yk = 0;
    yk_d = 0;
    //! Multipliers and accumulator
    for(i = 0; i<N_aj; i++)
    {
        yk += int(floor(double(inputAJ[aj_delay+i]*w_aj1[i])/double(1<<multround)+0.5));
    }
    yk = int(floor(double(yk)/double(1<<(sumround))+0.5));
    //Difference
    eps = aj_in*1-yk3;//8
	vivod1<<eps<<endl;
	//vivod3<<aj_in<<' '<<inputAJ[80+256]<<' '<<yk3<<' '<<aj_in2<<endl;
    //vivod2<<inputAJ[0]<<' '<<yk3<<' '<<eps<<endl;
    if(log_en==1)
    {
        vivod2<<inputAJ[0]<<' '<<yk3<<endl;
    }
    yk3=yk2;
    yk2=yk1;
    yk1=yk;

    //! Coefficient computing
    for(i=0;i<N_aj;i++) w_aj1[i] = w_aj[i];
    int ldel = 1<<(coeffcompldel);//16 15
    for(i=0;i<N_aj;i++)
    {
        w_lms[i] = eps*(inputAJ[aj_delay+i+3]*1);//eps1  +3
        //Накопление остатка целочисленного деления
        w_aj[i] += w_lms[i]/ldel;
        //учёт остатка от целочисленного деления
        if(waj_res[i]>(ldel-1) || waj_res[i]<-(ldel))
        {
            w_aj[i] += waj_res[i]/(ldel);
            waj_res[i] = (w_lms[i]%(ldel)) + waj_res[i]%(ldel);
        }
        else waj_res[i] = waj_res[i] + w_lms[i]%(ldel);
        //недопускание переполнения коэффициентов МНК
        if(w_aj[i]>((1<<(fltcoeff_bw-1)) - 1) || w_aj[i]<-(1<<(fltcoeff_bw-1)))
        {
            if(w_aj[i]>0) w_aj[i] = ((1<<(fltcoeff_bw-1)) - 1);
            else w_aj[i] = -((1<<(fltcoeff_bw-1)) - 1);
        }
    }

    if(log_en==1)
    {
        //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod3, " ")); vivod3<<endl;
        //vivod3<<w_lms[0]<<' '<<eps<<' '<<(inputAJ[aj_delay+0+3]*1)<<' '<<w_lms[0]/(1<<coeffcompldel)<<' '<<w_aj[0]<<endl;
        //copy(&waj_res[0], &waj_res[N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;
        copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
        //vivod2<<aj_in<<endl;
    }
    //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;

    /*
    w_aj[0]  = -63 ;
    w_aj[1]  = 110 ;
    w_aj[2]  = -61 ;
    w_aj[3]  = -38 ;
    w_aj[4]  = 105 ;
    w_aj[5]  = -90 ;
    w_aj[6]  = -2  ;
    w_aj[7]  = 85  ;
    w_aj[8]  = -101;
    w_aj[9]  = 32  ;
    w_aj[10] = 63  ;
    w_aj[11] = -104;
    w_aj[12] = 61  ;
    w_aj[13] = 31  ;
    w_aj[14] = -102;
    w_aj[15] = 86  ;
    w_aj[16] = 0   ;
    w_aj[17] = -84 ;
    w_aj[18] = 97  ;
    w_aj[19] = -34 ;
    w_aj[20] = -60 ;
    w_aj[21] = 104 ;
    w_aj[22] = -59 ;
    w_aj[23] = -34 ;
    w_aj[24] = 95  ;
    w_aj[25] = -86 ;
    w_aj[26] = 3   ;
    w_aj[27] = 87  ;
    w_aj[28] = -97 ;
    w_aj[29] = 29  ;
    w_aj[30] = 61  ;
    w_aj[31] = -104;
    w_aj[32] = 62  ;
    w_aj[33] = 35  ;
    w_aj[34] = -94 ;
    w_aj[35] = 85  ;
    w_aj[36] = 3   ;
    w_aj[37] = -81 ;
    w_aj[38] = 100 ;
    w_aj[39] = -32 ;
    w_aj[40] = -59 ;
    w_aj[41] = 102 ;
    w_aj[42] = -60 ;
    w_aj[43] = -35 ;
    w_aj[44] = 96  ;
    w_aj[45] = -79 ;
    w_aj[46] = 0   ;
    w_aj[47] = 85  ;
    w_aj[48] = -95 ;
    w_aj[49] = 34  ;
    w_aj[50] = 62  ;
    w_aj[51] = -101;
    w_aj[52] = 58  ;
    w_aj[53] = 33  ;
    w_aj[54] = -96 ;
    w_aj[55] = 83  ;
    w_aj[56] = 2   ;
    w_aj[57] = -82 ;
    w_aj[58] = 96  ;
    w_aj[59] = -30 ;
    w_aj[60] = -62 ;
    w_aj[61] = 103 ;
    w_aj[62] = -58 ;
    w_aj[63] = -35 ;
    w_aj[64] = 96  ;
    w_aj[65] = -82 ;
    w_aj[66] = -5  ;
    w_aj[67] = 81  ;
    w_aj[68] = -97 ;
    w_aj[69] = 30  ;
    w_aj[70] = 57  ;
    w_aj[71] = -100;
    w_aj[72] = 56  ;
    w_aj[73] = 33  ;
    w_aj[74] = -97 ;
    w_aj[75] = 82  ;
    w_aj[76] = 3   ;
    w_aj[77] = -81 ;
    w_aj[78] = 97  ;
    w_aj[79] = -29 ;
	*/
    //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
    //vivod1<<w_aj[0]<<endl;
    //vivod2<<waj_res[0]<<endl;
    //vivod3<<waj_res[0]/(1<<coeffcompldel)<<endl;
    //copy(&waj_res[0], &waj_res[N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;
    //vivod3<<eps<<endl;
    //if(compute_coeff_en==1)
    //{
    //vivod1<<eps<<endl;
    //vivod3<<inputAJ[0]<<endl;
    //vivod3<<((aj_in<<16)>>16)<<" "<<inputAJ[0]<<endl;
    //}
    return eps;
}

//! Complex Anti-Jamming filter
int Filter::AJFilterNew2_cplx(int aj_in,int aj_in2,short k,short compute_coeff_en, short log_en,int point[])
{
    //! Delay line
    for(i = (N_aj+aj_delay+aj_d_delay+2); i>0; i--) inputAJ[i] = inputAJ[i-1];
    inputAJ[0] = aj_in;
    for(i = (N_aj+aj_delay+aj_d_delay+2); i>0; i--) inputAJ2[i] = inputAJ2[i-1];
    inputAJ2[0] = aj_in2;

    yk = 0;
    ykk = 0;

    //! Multipliers and accumulator
    for(i = 0; i<N_aj; i++)
    {
        yk += int(floor(double(inputAJ[aj_delay+i]*w_aj1[i] - inputAJ2[aj_delay+i]*w_aj22[i])/double(1<<multround)+0.5));
        ykk += int(floor(double(inputAJ2[aj_delay+i]*w_aj1[i] + inputAJ[aj_delay+i]*w_aj22[i])/double(1<<multround)+0.5));
    }
    yk = int(floor(double(yk)/double(1<<(sumround))+0.5));
    ykk = int(floor(double(ykk)/double(1<<(sumround))+0.5));
    //! Difference
    eps=inputAJ[0]*1-yk3;//8
    yk3=yk2;
    yk2=yk1;
    yk1=yk;
    eps2=inputAJ2[0]*1-ykk3;//8
    ykk3=ykk2;
    ykk2=ykk1;
    ykk1=ykk;

    //vivod1<<inputAJ[0]<<' '<<yk3<<' '<<inputAJ2[0]<<' '<<ykk3<<endl;
    //! Coefficient computing
    //if(compute_coeff_en==0)eps1=eps;
    int i_f,i_l;
    for(i=0;i<N_aj;i++) w_aj1[i] = w_aj[i];
    for(i=0;i<N_aj;i++) w_aj22[i] = w_aj2[i];
    //if(compute_coeff_en==1)
    {
        int ldel = 1<<(coeffcompldel);//16 15
        switch(compute_coeff_en)
        {
            case 0:    i_f = 0; i_l=20; break;
            case 1:    i_f = 20; i_l=40; break;
            case 2:    i_f = 40; i_l=60; break;
            case 3:    i_f = 60; i_l=80; break;
            case 4:    i_f = 64; i_l=80; break;
            default:   i_f = 0; i_l=16; break;
        }
        //for(i=i_f;i<i_l;i=i+1)
        //for(i=16*compute_coeff_en;i<(N_aj/5*(compute_coeff_en+1));i=i+1)
        //for(i=40*compute_coeff_en;i<(N_aj/2*(compute_coeff_en+1));i=i+1)
        //for(i=10*compute_coeff_en;i<(N_aj/8*(compute_coeff_en+1));i=i+1)
        //for(i=1*compute_coeff_en;i<(N_aj/80*(compute_coeff_en+1));i=i+1)
        for(i=0;i<N_aj;i++)
        //for(int j=16*compute_coeff_en;j<(N_aj/5*(compute_coeff_en+1));j=j+1)
        {//i = point[j];
            //w_aj1[i] = w_aj[i];
            //w_lms[i] = ((eps*(inputAJ[aj_delay+i+3]*16))<<(32 - coeffcompmult_bw))>>(32 - coeffcompmult_bw);//eps1
            w_lms[i] = eps*(inputAJ[aj_delay+i+3]*1) + eps2*(inputAJ2[aj_delay+i+3]*1);
            w_lms2[i] = eps*(-inputAJ2[aj_delay+i+3]*1) + eps2*(inputAJ[aj_delay+i+3]*1);
            //Накопление остатка целочисленного деления
            w_aj[i] += w_lms[i]/ldel;
            w_aj2[i] += w_lms2[i]/ldel;
            //if(w_lms[i]/(ldel)>0) w_aj[i] +=1;
            //else if(w_lms[i]/(ldel)<0) w_aj[i] -=1;
            //учёт остатка от целочисленного деления
            if(waj_res[i]>(ldel-1) || waj_res[i]<-(ldel))
            {
                w_aj[i] += waj_res[i]/(ldel);
                waj_res[i] = (w_lms[i]%(ldel)) + waj_res[i]%(ldel);
            }
            else waj_res[i] = waj_res[i] + w_lms[i]%(ldel);
            //недопускание переполнения коэффициентов МНК
            if(w_aj[i]>((1<<(fltcoeff_bw-1)) - 1) || w_aj[i]<-(1<<(fltcoeff_bw-1)))
            {
                if(w_aj[i]>0) w_aj[i] = ((1<<(fltcoeff_bw-1)) - 1);
                else w_aj[i] = -((1<<(fltcoeff_bw-1)) - 1);
            }

            if(waj_res2[i]>(ldel-1) || waj_res2[i]<-(ldel))
            {
                w_aj2[i] += waj_res2[i]/(ldel);
                waj_res2[i] = (w_lms2[i]%(ldel)) + waj_res2[i]%(ldel);
            }
            else waj_res2[i] = waj_res2[i] + w_lms2[i]%(ldel);
            //недопускание переполнения коэффициентов МНК
            if(w_aj2[i]>((1<<(fltcoeff_bw-1)) - 1) || w_aj2[i]<-(1<<(fltcoeff_bw-1)))
            {
                if(w_aj2[i]>0) w_aj2[i] = ((1<<(fltcoeff_bw-1)) - 1);
                else w_aj2[i] = -((1<<(fltcoeff_bw-1)) - 1);
            }
        }
    }

    if(log_en==1)
    {
        //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod3, " ")); vivod3<<endl;
        //vivod3<<w_lms[0]<<' '<<eps<<' '<<(inputAJ[aj_delay+0+3]*1)<<' '<<w_lms[0]/(1<<coeffcompldel)<<' '<<w_aj[0]<<endl;
        //copy(&waj_res[0], &waj_res[N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;
    }
    return eps;
}

//! Simple complex Anti-Jamming filter (for slow jams)
int Filter::AJFilterNew2_cplxs(int aj_in,int aj_in2,short k,short compute_coeff_en, short log_en,int point[])
{
    //! Delay line
    for(i = (N_aj+aj_delay+aj_d_delay+2); i>0; i--) inputAJ[i] = inputAJ[i-1];
    inputAJ[0] = aj_in;
    for(i = (N_aj+aj_delay+aj_d_delay+2); i>0; i--) inputAJ2[i] = inputAJ2[i-1];
    inputAJ2[0] = aj_in2;

    yk = 0;

    //! Multipliers and accumulator
    yk += int(floor(double(inputAJ[aj_delay+0]*w_aj1[0])/double(1<<multround)+0.5));
    //yk += int(floor(double(inputAJ2[aj_delay+0]*w_aj22[0])/double(1<<multround)+0.5));
    yk = int(floor(double(yk)/double(1<<(sumround))+0.5));
    //! Difference
    eps=inputAJ[0]*1-yk3;//8
    yk3=yk2;
    yk2=yk1;
    yk1=yk;

    //vivod1<<inputAJ[0]<<' '<<yk3<<' '<<inputAJ2[0]<<' '<<ykk3<<endl;
    //! Coefficient computing
    //if(compute_coeff_en==0)eps1=eps;
    int i_f,i_l;
    for(i=0;i<N_aj;i++) w_aj1[i] = w_aj[i];
    for(i=0;i<N_aj;i++) w_aj22[i] = w_aj2[i];
    //if(compute_coeff_en==1)
    {
        int ldel = 1<<(coeffcompldel);//16 15
        switch(compute_coeff_en)
        {
            case 0:    i_f = 0; i_l=20; break;
            case 1:    i_f = 20; i_l=40; break;
            case 2:    i_f = 40; i_l=60; break;
            case 3:    i_f = 60; i_l=80; break;
            case 4:    i_f = 64; i_l=80; break;
            default:   i_f = 0; i_l=16; break;
        }
        //for(i=i_f;i<i_l;i=i+1)
        //for(i=16*compute_coeff_en;i<(N_aj/5*(compute_coeff_en+1));i=i+1)
        //for(i=40*compute_coeff_en;i<(N_aj/2*(compute_coeff_en+1));i=i+1)
        //for(i=10*compute_coeff_en;i<(N_aj/8*(compute_coeff_en+1));i=i+1)
        //for(i=1*compute_coeff_en;i<(N_aj/80*(compute_coeff_en+1));i=i+1)
        for(i=0;i<1;i++)
        //for(int j=16*compute_coeff_en;j<(N_aj/5*(compute_coeff_en+1));j=j+1)
        {//i = point[j];
            //w_aj1[i] = w_aj[i];
            //w_lms[i] = ((eps*(inputAJ[aj_delay+i+3]*16))<<(32 - coeffcompmult_bw))>>(32 - coeffcompmult_bw);//eps1
            w_lms[i] = eps*(inputAJ[aj_delay+i+3]*1);
            w_lms2[i] = eps*(inputAJ2[aj_delay+i+3]*1);
            //Накопление остатка целочисленного деления
            w_aj[i] += w_lms[i]/ldel;
            w_aj2[i] += w_lms2[i]/ldel;
            //if(w_lms[i]/(ldel)>0) w_aj[i] +=1;
            //else if(w_lms[i]/(ldel)<0) w_aj[i] -=1;
            //учёт остатка от целочисленного деления
            if(waj_res[i]>(ldel-1) || waj_res[i]<-(ldel))
            {
                w_aj[i] += waj_res[i]/(ldel);
                waj_res[i] = (w_lms[i]%(ldel)) + waj_res[i]%(ldel);
            }
            else waj_res[i] = waj_res[i] + w_lms[i]%(ldel);
            //недопускание переполнения коэффициентов МНК
            if(w_aj[i]>((1<<(fltcoeff_bw-1)) - 1) || w_aj[i]<-(1<<(fltcoeff_bw-1)))
            {
                if(w_aj[i]>0) w_aj[i] = ((1<<(fltcoeff_bw-1)) - 1);
                else w_aj[i] = -((1<<(fltcoeff_bw-1)) - 1);
            }

            if(waj_res2[i]>(ldel-1) || waj_res2[i]<-(ldel))
            {
                w_aj2[i] += waj_res2[i]/(ldel);
                waj_res2[i] = (w_lms2[i]%(ldel)) + waj_res2[i]%(ldel);
            }
            else waj_res2[i] = waj_res2[i] + w_lms2[i]%(ldel);
            //недопускание переполнения коэффициентов МНК
            if(w_aj2[i]>((1<<(fltcoeff_bw-1)) - 1) || w_aj2[i]<-(1<<(fltcoeff_bw-1)))
            {
                if(w_aj2[i]>0) w_aj2[i] = ((1<<(fltcoeff_bw-1)) - 1);
                else w_aj2[i] = -((1<<(fltcoeff_bw-1)) - 1);
            }
        }
    }

    if(log_en==1)
    {
        //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod3, " ")); vivod3<<endl;
        //vivod3<<w_lms[0]<<' '<<eps<<' '<<(inputAJ[aj_delay+0+3]*1)<<' '<<w_lms[0]/(1<<coeffcompldel)<<' '<<w_aj[0]<<endl;
        //copy(&waj_res[0], &waj_res[N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;
    }
    //copy(&w_aj[0], &w_aj[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
    //vivod1<<w_aj[0]<<endl;
    //vivod2<<waj_res[0]<<endl;
    //vivod3<<waj_res[0]/(1<<coeffcompldel)<<endl;
    //copy(&waj_res[0], &waj_res[N_aj], ostream_iterator<double>(vivod2, " ")); vivod2<<endl;
    //vivod3<<eps<<endl;
    //if(compute_coeff_en==1)
    //{
    //vivod1<<eps<<endl;
    //vivod3<<inputAJ[0]<<endl;
    //vivod3<<((aj_in<<16)>>16)<<" "<<inputAJ[0]<<endl;
    //}
    return eps;
}
double Filter::AJFilter_double(short aj_in_opor,short aj_in,short k,short compute_coeff_en, short log_en)
{
 //if(sw1==0)return re;
 rotate(&inputAJ_d[0], &inputAJ_d[N_aj+aj_delay-1], &inputAJ_d[N_aj+aj_delay]);
 inputAJ_d[0]=double(aj_in)/2000;//-aj_in_z;
 //aj_in_z=(floor(double(aj_in_opor)/1.0+0.0))*1.0;
 //rotate(&input_delay_d[0], &input_delay_d[aj_d_delay-1], &input_delay_d[aj_d_delay]);
 //input_delay_d[0]=double(aj_in);
 // Умножение на коэффициенты фильтра
 transform(&inputAJ_d[aj_delay], &inputAJ_d[N_aj+aj_delay], &w_aj_d[0], &inputAJ_rev_d[0], multiplies<double>());
 //vivod1<<inputAJ_rev[0]<<endl;
 // Уменьшение разрядности после перемножения
 transform(&inputAJ_rev_d[0], &inputAJ_rev_d[N_aj], &inputAJ_rev_d[0], bind2nd(divides<double>(), double(1<<0)));//(1<<k)
 //for(i=0;i<N_aj;i++) inputAJ_rev_d[i]=floor(inputAJ_rev_d[i]+0.5);
 // Суммирование результатов перемножения
 yk_d=accumulate<double *,double>(&inputAJ_rev_d[0],&inputAJ_rev_d[N_aj],0);
 //yk_d=(floor(yk_d/128.0+0.5))*128.0;
 eps_d=inputAJ_d[0]-yk_d;
 //vivod1<<w_aj_d[0]<<endl;
 //vivod2<<aj_in<<endl;
 //vivod3<<eps_d<<endl;
 //if(compute_coeff_en==1)copy(&w_aj_d[0], &w_aj_d[N_aj], ostream_iterator<double>(vivod1, " "));
 for(i=0;i<N_aj;i++) w_aj_d[i]=w_aj_d[i]+2*inputAJ_d[aj_delay+i]*eps_d*mju_d;
 if(log_en==1)
 {
	copy(&w_aj_d[0], &w_aj_d[N_aj], ostream_iterator<double>(vivod1, " ")); vivod1<<endl;
 }
 return eps_d;
}


short Filter::AJFilter_cplx(short aj_in_i,short aj_in_q,short k,short compute_coeff_en)
{
 //rotate(&inputAJ_d_cplx_i[0], &inputAJ_d_cplx_i[N_aj+aj_delay-1], &inputAJ_d_cplx_i[N_aj+aj_delay]);
 //inputAJ_d_cplx_i[0]=double(aj_in_i);
 //rotate(&inputAJ_d_cplx_q[0], &inputAJ_d_cplx_q[N_aj+aj_delay-1], &inputAJ_d_cplx_q[N_aj+aj_delay]);
 //inputAJ_d_cplx_q[0]=double(aj_in_q);
 for(i=(N_aj+aj_delay-1);i>0;i--)
 {
  inputAJ_d_cplx_i[i]=inputAJ_d_cplx_i[i-1];
  inputAJ_d_cplx_q[i]=inputAJ_d_cplx_q[i-1];
 }
 inputAJ_d_cplx_i[0]=double(aj_in_i);
 inputAJ_d_cplx_q[0]=double(aj_in_q);
 yk_d_cplx_i=0;
 yk_d_cplx_q=0;
 for(i=0;i<N_aj;i++)
 {
  //inputAJ_rev_d_cplx_i[i]=inputAJ_d_cplx_i[aj_delay+i]*w_aj_d_cplx_i[i]-inputAJ_d_cplx_q[aj_delay+i]*w_aj_d_cplx_q[i];
  //inputAJ_rev_d_cplx_q[i]=inputAJ_d_cplx_q[aj_delay+i]*w_aj_d_cplx_i[i]+inputAJ_d_cplx_i[aj_delay+i]*w_aj_d_cplx_q[i];
  yk_d_cplx_i+=(inputAJ_d_cplx_i[aj_delay+i]*w_aj_d_cplx_i[i]-inputAJ_d_cplx_q[aj_delay+i]*w_aj_d_cplx_q[i])/double(65536);
  yk_d_cplx_q+=(inputAJ_d_cplx_q[aj_delay+i]*w_aj_d_cplx_i[i]+inputAJ_d_cplx_i[aj_delay+i]*w_aj_d_cplx_q[i])/double(65536);
 }
 // Уменьшение разрядности после перемножения
 //transform(&inputAJ_rev_d_cplx_i[0], &inputAJ_rev_d_cplx_i[N_aj], &inputAJ_rev_d_cplx_i[0], bind2nd(divides<double>(), double(32768)));//(1<<k)
 //transform(&inputAJ_rev_d_cplx_q[0], &inputAJ_rev_d_cplx_q[N_aj], &inputAJ_rev_d_cplx_q[0], bind2nd(divides<double>(), double(32768)));//(1<<k)
 // Суммирование результатов перемножения
 //yk_d_cplx_i=accumulate<double *,double>(&inputAJ_rev_d_cplx_i[0],&inputAJ_rev_d_cplx_i[N_aj],0);
 //yk_d_cplx_q=accumulate<double *,double>(&inputAJ_rev_d_cplx_q[0],&inputAJ_rev_d_cplx_q[N_aj],0);
 eps_d_cplx_i=inputAJ_d_cplx_i[0]-yk_d_cplx_i;
 eps_d_cplx_q=inputAJ_d_cplx_q[0]-yk_d_cplx_q;
 if(compute_coeff_en==1)
 {
  for(i=0;i<N_aj;i++)
  {
   w_aj_d_cplx_i[i]=w_aj_d_cplx_i[i]+2.0*(inputAJ_d_cplx_i[aj_delay+i]*eps_d_cplx_i+inputAJ_d_cplx_q[aj_delay+i]*eps_d_cplx_q)*mju_d;
   w_aj_d_cplx_q[i]=w_aj_d_cplx_q[i]+2.0*(-inputAJ_d_cplx_q[aj_delay+i]*eps_d_cplx_i+inputAJ_d_cplx_i[aj_delay+i]*eps_d_cplx_q)*mju_d;
  }
 }
 //if(compute_coeff_en==1)copy(&w_aj_d[0], &w_aj_d[N_aj], ostream_iterator<double>(vivod1, " "));
 //if(compute_coeff_en==1){copy(&w_aj_d_cplx_i[0], &w_aj_d_cplx_i[N_aj], ostream_iterator<double>(vivod3, " "));vivod3<<endl;}
 //copy(&inputAJ_d_cplx_q[aj_delay], &inputAJ_d_cplx_q[aj_delay+N_aj], ostream_iterator<double>(vivod3, " "));vivod3<<endl;
 //vivod1<<w_aj_d_cplx_i[0]<<endl;
 //vivod2<<yk_d_cplx_i<<endl;
 //vivod3<<eps_d_cplx_i<<endl;
 return short(eps_d_cplx_i);
}

void Filter::LoadAJFilter_adv(short N, short adapt_delay, short rez_aj[])
{
 Naj=N;
 Naj_delay=0;
 aj_comb_i[0]=0;
 aj_comb_q[0]=0;
 aj_comb_d_i=0;
 aj_comb_d_q=0;
 //r1_i=0;
 //r1_q=0;
// r1_i_pos=0;
 //r1_q_pos=0;
 rn=0.99999;
 //rnn=pow(rn,double(Naj));
 //rnn=0.87489824456875;//512 (9) //0.9928375;//rnn=0.8666;//rnn=0.96875;
 //rnn=0.25;//256 (5) 0.73979227953412 0.734375
 rnn[0]=0.25;
 for(i=1;i<8;i++) rnn[i] = rnn[i-1] + 0.0625;
 //0.99489305938955004392597107192554
 i=0;

 sdvig=1048576/64*64;//*64;//  /4
 hc_aj=new short[Naj];
 aj_cplx_const_i=new short[Naj];
 aj_cplx_const_q=new short[Naj];
 w_aj_d_i=new double[Naj];
 w_aj_d_q=new double[Naj];
 r_i=new double[Naj];
 r_q=new double[Naj];
 r_i_n=new double[Naj];
 r_q_n=new double[Naj];
 input_d_i=new double[Naj+1];
 input_d_q=new double[Naj+1];

 w_aj_i=new int[Naj];
 w_aj_q=new int[Naj];
 r_i_s=new int[Naj];
 r_q_s=new int[Naj];
 r_i_n_s=new int[Naj];
 r_q_n_s=new int[Naj];
 input_i=new short[Naj+1];
 input_q=new short[Naj+1];
 waj_res_i=new int[Naj];
 waj_res_q=new int[Naj];
 w_lms_i=new int[Naj];
 w_lms_q=new int[Naj];

 r_i_s_46=0;
 r_q_s_46=0;
 eps_i1=0;
 eps_q1=0;
 eps_i2=0;
 eps_q2=0;

 for(i=0;i<(Naj+1);i++)
 {
  input_d_i[i]=0;
  input_d_q[i]=0;
  input_i[i]=0;
  input_q[i]=0;
 }
 for(i=0;i<Naj;i++)
 {
  w_aj_d_i[i]=0;
  w_aj_d_q[i]=0;
  w_aj_i[i]=0;
  w_aj_q[i]=0;
  hc_aj[i]=rez_aj[i];
  r_i[i]=0;
  r_q[i]=0;
  r_i_n[i]=0;
  r_q_n[i]=0;
  r_i_s[i]=0;
  r_q_s[i]=0;
  r_i_n_s[i]=0;
  r_q_n_s[i]=0;
  //aj_cplx_const_i[i]=short(floor(cos(2*pi*i/double(Naj))*2047+0.5));
  //aj_cplx_const_q[i]=short(floor(sin(2*pi*i/double(Naj))*2047+0.5));
  //if(double(aj_cplx_const_i[i]*aj_cplx_const_i[i]/2048/2048+aj_cplx_const_q[i]*aj_cplx_const_q[i]/2048/2048)>1 && hc_aj[i]==1)
  // cout<<"AJ coeff singularity!"<<" N = "<<Naj<<"; n = "<<i<<endl;
  aj_cplx_const_i[i]=short(floor(cos(2*pi*i/double(Naj))*127+0.5));
  aj_cplx_const_q[i]=short(floor(sin(2*pi*i/double(Naj))*127+0.5));
  waj_res_i[i]=0;
  waj_res_q[i]=0;
  w_lms_i[i]=0;
  w_lms_q[i]=0;
 }
 //aj_cplx_const_i[8]=2038;
 //aj_cplx_const_q[8]=201;

}

double Filter::ComplexConst_i(short n)
{
 return cos(2*pi*n/double(Naj));
}
double Filter::ComplexConst_q(short n)
{
 return sin(2*pi*n/double(Naj));
}
double Filter::AJFilter_adv_double(double aj_in_i,double aj_in_q,short k,short compute_coeff_en)
{
 for(i=Naj;i>0;i--)
 {
  input_d_i[i]=input_d_i[i-1];
  input_d_q[i]=input_d_q[i-1];
 }
 input_d_i[0]=aj_in_i;
 input_d_q[0]=aj_in_q;
 aj_comb_d_i=input_d_i[0]-input_d_i[Naj]*rnn[0];
 aj_comb_d_q=input_d_q[0]-input_d_q[Naj]*rnn[0];

 for(i=0;i<Naj;i++)//hc_aj[i] - массив включения резонатов
  if(hc_aj[i]==1)
  {//w_aj_d_i[504]=w_aj_d_i[8];
   //w_aj_d_i[8]=1.27726;
   //w_aj_d_i[504]=1.27726;
   r_i_n[i]=aj_comb_d_i+rn*(r_i[i]*ComplexConst_i(i)-r_q[i]*ComplexConst_q(i));
   r_q_n[i]=aj_comb_d_q+rn*(r_i[i]*ComplexConst_q(i)+r_q[i]*ComplexConst_i(i));
   r_i[i]=r_i_n[i];
   r_q[i]=r_q_n[i];
   r_i_n[i]=r_i[i]*w_aj_d_i[i]-1*r_q[i]*w_aj_d_q[i];
   r_q_n[i]=r_i[i]*w_aj_d_q[i]*1+r_q[i]*w_aj_d_i[i];
  }
 //eps_d_i=
 eps_d_i=input_d_i[Naj/2]-accumulate<double *,double>(&r_i_n[0],&r_i_n[Naj],0)/double(Naj);
 eps_d_q=input_d_q[Naj/2]-accumulate<double *,double>(&r_q_n[0],&r_q_n[Naj],0)/double(Naj);
 //if(compute_coeff_en==1)
 {
  for(i=0;i<Naj;i++)
   if(hc_aj[i]==1)
   {
    w_aj_d_i[i]=w_aj_d_i[i]+2.0*(eps_d_i*r_i[i]+eps_d_q*r_q[i])*0.0001;
	w_aj_d_q[i]=w_aj_d_q[i]+2.0*(-eps_d_i*r_q[i]+eps_d_q*r_i[i])*0.0001;
   }
 }
 //eps_d_i=r_i_n[8]/double(Naj);
 //eps_d_q=r_q_n[8]/double(Naj);
 //if(compute_coeff_en==1)copy(&w_aj_d[0], &w_aj_d[N_aj], ostream_iterator<double>(vivod1, " "));
 //if(compute_coeff_en==1){copy(&w_aj_d_cplx_i[0], &w_aj_d_cplx_i[N_aj], ostream_iterator<double>(vivod3, " "));vivod3<<endl;}
 //copy(&inputAJ_d_cplx_q[aj_delay], &inputAJ_d_cplx_q[aj_delay+N_aj], ostream_iterator<double>(vivod3, " "));vivod3<<endl;
 if(compute_coeff_en==1)
 {
  vivod1<<w_aj_d_i[4]<<' '<<w_aj_d_q[4]<<endl;
  //vivod2<<w_aj_d_i[77]<<' '<<w_aj_d_q[77]<<endl;
  //vivod3<<w_aj_d_i[78]<<' '<<w_aj_d_q[78]<<endl;
 }
 //vivod3<<r_i[504]*w_aj_d_i[8]/double(Naj-1)<<endl;
 return 0;
}
//! Целочисленный модифицированный адапривный фильтр
short Filter::AJFilter_adv(short aj_in_i,short aj_in_q,short k,short compute_coeff_en,short log_en, short ha_aj[])
{
 //! Линия задержки адаптивного фильтра, Naj-длина линии задержки
 for(i=Naj;i>0;i--)
 {
  input_i[i]=input_i[i-1];
  input_q[i]=input_q[i-1];
 }
 input_i[0]=aj_in_i/(1<<2);//cout<<input_i[0]<<endl;
 input_q[0]=aj_in_q/(1<<2);
 //! Умножение задержанного сигнала на коэффициент rnn и вычитание из основного сигнала.
 //! Коэффициент rnn рассчитывается, как sqrt(aj_cplx_const_i[i]*aj_cplx_const_i[i]+aj_cplx_const_i[i]*aj_cplx_const_i[i])
 //! возвердённое в степень Naj
 for(i=0;i<8;i++)
 {
  aj_comb_i[i]=input_i[0]-short(floor(input_i[Naj]*rnn[i]+0.5));
  aj_comb_q[i]=input_q[0]-short(floor(input_q[Naj]*rnn[i]+0.5));
 }
 //! Обработка сигнала резонаторами (интеграторами)
 //! hc_aj[i] - массив включения резонатов
 for(i=0;i<Naj;i++)
 if(ha_aj[i]==1)
 {
	 short j;
	 if(i==46) j=5;
	 //else j=0;
	 else j=2;
  //! Округление сигнала r_i_s[i] на число разрядов Naj/2, умножение на комплексную константу, округление результата
  //! на число разрядов, зависящее от разрядности целочисленной комплексной константы 2048/(Naj/2)
  r_i_n_s[i]=aj_comb_i[j]+int(((r_i_s[i])/(Naj>>1))*aj_cplx_const_i[i]-((r_q_s[i])/(Naj>>1))*aj_cplx_const_q[i])/(128/(Naj>>1));//деление в двух точках с целью получения верной комплексной константы
  r_q_n_s[i]=aj_comb_q[j]+int(((r_i_s[i])/(Naj>>1))*aj_cplx_const_q[i]+((r_q_s[i])/(Naj>>1))*aj_cplx_const_i[i])/(128/(Naj>>1));//уменьшает размер перемножителя при длинной линии задержки																								                  //первое деление на 1 меньше Nj, второе = (размер aj_cplx_const_i)/первое
  //cout<<r_i_n_s[46]<<endl;
  r_i_s[i]=((r_i_n_s[i])<<12)>>12;//cout<<(r_i_s[46]>>17)<<endl;
  r_q_s[i]=((r_q_n_s[i])<<12)>>12;
  //! Умножение на комплексный коэффициент LMS (МНК) для компенсации амплитуды и фазы помехи
  r_i_n_s[i]=(r_i_s[i]/Naj)*w_aj_i[i]-(r_q_s[i]/Naj)*w_aj_q[i];//if(i==89)cout<<w_aj_i[i]<<" "<<w_aj_q[i]<<endl;
  r_q_n_s[i]=(r_i_s[i]/Naj)*w_aj_q[i]+(r_q_s[i]/Naj)*w_aj_i[i];
 }
 //! Округление суммы выходов всех резонаторов и вычитание её из основного сигнала, задержанного на Naj/2
 //eps_i=input_i[Naj>>1]-accumulate<int *,int>(&r_i_n_s[0],&r_i_n_s[Naj],0)/(16384/4);//16384 - приблизительная разрядность коэффициентов
 //eps_q=input_q[Naj>>1]-accumulate<int *,int>(&r_q_n_s[0],&r_q_n_s[Naj],0)/(16384/4);
 eps_i=0;
 eps_q=0;
 for(i=0;i<Naj;i++)
 if(ha_aj[i]==1)
 {
	eps_i+=((r_i_n_s[i]/2)<<1)>>1;//4096
	eps_q+=((r_q_n_s[i]/2)<<1)>>1;//4096
 }
 //vivod3<<input_i[Naj>>1]<<' '<<eps_i<<endl;
 eps_i=input_i[Naj>>1]-eps_i/1;// /(16384/4);
 eps_q=input_q[Naj>>1]-eps_q/1;// /(16384/4);
 //cout<<r_i_s[46]<<endl;
 //! Kоэффициент скорости сходимости МНК int sdvig=1048576;
 //! Вычисление коэффициентов целочисленным методом наименьших квадратов

 if(compute_coeff_en==1)
 {
  for(i=0;i<Naj;i++)
   if(ha_aj[i]==1)
   {


    w_lms_i[i]=eps_i2*(r_i_s[i]/Naj)+eps_q2*(r_q_s[i]/Naj);//cout<<w_lms_i[4]/sdvig<<endl;
    w_lms_q[i]=-eps_i2*(r_q_s[i]/Naj)+eps_q2*(r_i_s[i]/Naj);

	w_aj_i[i] += (w_lms_i[i]+waj_res_i[i])/sdvig;
    waj_res_i[i] = (w_lms_i[i]+waj_res_i[i])%sdvig;

	w_aj_q[i] += (w_lms_q[i]+waj_res_q[i])/sdvig;
    waj_res_q[i] = (w_lms_q[i]+waj_res_q[i])%sdvig;
	/*
    waj_res_i[i]+=w_lms_i[i]%sdvig;
    waj_res_q[i]+=w_lms_q[i]%sdvig;

    w_aj_i[i]+=w_lms_i[i]/sdvig;
	w_aj_q[i]+=w_lms_q[i]/sdvig;
    //учёт остатка от целочисленного деления
	if(waj_res_i[i]>=sdvig || waj_res_i[i]<=-sdvig) {w_aj_i[i]+=waj_res_i[i]/(sdvig); waj_res_i[i]=0;}//-=(sdvig/4);}
    if(waj_res_q[i]>=sdvig || waj_res_q[i]<=-sdvig) {w_aj_q[i]+=waj_res_q[i]/(sdvig); waj_res_q[i]=0;}//-=(sdvig/4);}
	*/
	if(w_aj_i[i]>1023 || w_aj_i[i]<-1024)
    if(w_aj_i[i]>0) w_aj_i[i] = 1023; else w_aj_i[i] = -1024;
	if(w_aj_q[i]>1023 || w_aj_q[i]<-1024)
    if(w_aj_q[i]>0) w_aj_q[i] = 1023; else w_aj_q[i] = -1024;
   }
 r_i_s_46=r_i_s[46];
 r_q_s_46=r_q_s[46];

 eps_i2=eps_i1;
 eps_q2=eps_q1;
 eps_i1=eps_i;
 eps_q1=eps_q;

 }
 //Debuging
 if(log_en==1)
 {
	 //vivod1<<aj_comb_i[13]<<' '<<aj_comb_q[13]<<endl;
	 //vivod1<<w_aj_i[46]<<' '<<w_aj_q[46]<<endl;
	 //vivod3<<aj_in_i<<endl;
	 //vivod1<<eps_i<<endl;
	 //vivod2<<eps_q<<endl;
    //vivod1<<waj_res_i[4]<<' '<<waj_res_q[4]<<endl;
  //vivod2<<w_aj_i[9]<<' '<<w_aj_q[9]<<endl;
  //vivod3<<eps_i<<' '<<accumulate<int *,int>(&r_i_n_s[0],&r_i_n_s[Naj],0)/16384<<endl;
 }
 return 0;
}

//! Спектральный анализатор
int Filter::SpectrAn(int sa_signal, short sa_reset)
{
 // Comb part
 rotate(&sa_stage_comb_z[0][0], &sa_stage_comb_z[0][sa_st0], &sa_stage_comb_z[0][sa_st0+1]);
 sa_stage_comb_z[0][0]=sa_signal;
 sa_stage_comb[0]=sa_signal-sa_stage_comb_z[0][sa_st0];

 rotate(&sa_stage_comb_z[1][0], &sa_stage_comb_z[1][sa_st1], &sa_stage_comb_z[1][sa_st1+1]);
 sa_stage_comb_z[1][0]=sa_stage_comb[1-1];
 sa_stage_comb[1]=sa_stage_comb[1-1]-sa_stage_comb_z[1][sa_st1];

 rotate(&sa_stage_comb_z[2][0], &sa_stage_comb_z[2][sa_st2], &sa_stage_comb_z[2][sa_st2+1]);
 sa_stage_comb_z[2][0]=sa_stage_comb[2-1];
 sa_stage_comb[2]=sa_stage_comb[2-1]-sa_stage_comb_z[2][sa_st2];

 // Int part
 sa_stage_int[0]+=sa_stage_comb[3-1];
 sa_cic_out=sa_stage_int[0];

 sa_stage_int[1]+=sa_cic_out;
 sa_cic_out=sa_stage_int[1];

 sa_stage_int[2]+=sa_cic_out;
 sa_cic_out=sa_stage_int[2];

 //vivod1<<stage_comb[cic_stage_number-1]<<endl;
 //vivod2<<stage_int[1]<<endl;
 //vivod3<<stage_int[0]<<endl;
 //vivod1<<cic_out<<' '<<(cic_out<<(19-k))<<' '<<((cic_out<<(19-k))>>19)<<endl;
 //cic_out=cic_out<<(19-k);
 //cic_out=cic_out>>19;
 if(sa_reset==1)
 {
  for(i=0;i<3;i++)
  {
   sa_stage_int[i]=0;
   sa_stage_comb[i]=0;
   for(short j=0;j<513;j++) sa_stage_comb_z[i][j]=0;
  }
 }
 return (sa_cic_out/(65536/2));//  /(1<<k));
}

short Filter::AGCLocata(short cicout, long long t)
{
    if (t%10 == 0)
    {
        if (f_sec_agc_locata_tr > 1) f_sec_agc_locata--;
        f_sec_agc_locata_tr = 0;
    }

    if (t%10 == 0)
    {
        if (f_sec_agc_locata_tr2 == 0) f_sec_agc_locata++;
        f_sec_agc_locata_tr2 = 0;
    }

    if (f_sec_agc_locata < 0) f_sec_agc_locata = 0;
    if (f_sec_agc_locata > 3) f_sec_agc_locata = 3;

    if (cicout >= 8000) f_sec_agc_locata_tr++;//15
    if (cicout >= 100) f_sec_agc_locata_tr2++;//16
    //vivod1<<cicout<<' '<<f_sec_agc_locata_tr<<' '<<f_sec_agc_locata_tr2<<endl;

    return f_sec_agc_locata;
}
short Filter::AGCSecondLocata(short cicout, long long t, short sh_en, short zero_en, double sh, short mult_width)
{
    if(abs(cicout)>max_locata) max_locata = abs(cicout);
    if(zero_en==0)
    {
        if (t%10 == 0)
        {
            L3 = (max_locata*(short(floor((1<<mult_width)/7+0.5))))>>mult_width;
            L2 = (max_locata*(short(floor((3<<mult_width)/7+0.5))))>>mult_width;
            L1 = (max_locata*(short(floor((5<<mult_width)/7+0.5))))>>mult_width;
            Ln = (max_locata*(short(floor((1<<mult_width)/7+0.5))))>>mult_width;
            max_locata = 0;
        }
    }
    else
    {
        if (t%10 == 0)
        {
            L3 = 0;
            L2 = max_locata*1/3;
            L1 = max_locata*2/3;
            Ln = max_locata*1/3;
            max_locata = 0;
        }
    }
    //vivod1<<L3<<' '<<L2<<' '<<L1<<endl;
    //vivod2<<cicout<<endl;
    //!Noise generator
    short g1s=0;
    short noise=0;

    g1s=geng1[0]^geng1[5]^geng1[7]^geng1[13];
    rotate(&geng1[0], &geng1[13], &geng1[14]);
    geng1[0]=g1s;
    if(geng1[13]==0) noise = (geng1[13]<<13) |(geng1[12]<<12) | (geng1[11]<<11) |(geng1[10]<<10) | (geng1[9]<<9) |(geng1[8]<<8) |
            (geng1[7]<<7) | (geng1[6]<<6) | (geng1[5]<<5) | (geng1[4]<<4) | (geng1[3]<<3) | (geng1[2]<<2) |
                    (geng1[1]<<1) | (geng1[0]<<0);
    else noise = (0xF<<14) | (geng1[13]<<13) |(geng1[12]<<12) | (geng1[11]<<11) |(geng1[10]<<10) | (geng1[9]<<9) |(geng1[8]<<8) |
        (geng1[7]<<7) | (geng1[6]<<6) | (geng1[5]<<5) | (geng1[4]<<4) | (geng1[3]<<3) | (geng1[2]<<2) |
                            (geng1[1]<<1) | (geng1[0]<<0);

    //! Quantization
    short last_out=0;
    //vivod1<<cicout<<' ';
    if(sh_en==1) cicout += ((noise*Ln)>>13);
    //vivod1<<((noise*Ln)>>13)<<' '<<noise<<' '<<Ln<<endl;
    //if(sh_en==1) cicout += short(sh*Ln);
    if(cicout>0)
    {
        if(cicout>=L3) last_out=1;
        if(cicout>=L2) last_out=2;
        if(cicout>=L1) last_out=3;
    }
    else if(cicout<0)
        {
            if(cicout<=(-1*L3)) last_out=-1;
            if(cicout<=(-1*L2)) last_out=-2;
            if(cicout<=(-1*L1)) last_out=-3;
        }

    return last_out;
}

//! Class contain AGC, PWM and DAC algorithms
AGC::AGC()
{
    vivod1.open("data_out/f1.txt",ios::out);
    vivod2.open("data_out/f2.txt",ios::out);
    vivod3.open("data_out/f3.txt",ios::out);

    overflow_adc=0;
    gain_coeff=0;	//начальное значение коэффициента усиления рассчитанное для As/sigma=0.1;
    gain_coeff3=0;
    count1=0;
    s_zero=0;

    overflow_bit=0;
    agc_mid_level=0;

    x = 0;
    y = 0;
    y1 = 0;
    y2 = 0;

    k_agc = 1;
    pwm_count = 0;
    pwm_out = 0;
    pwm_filter_out = 0;
    pwm_z = 0;
    last_out = 0;
    L1 = 0;L2 = 0;L3 = 0;
    delay_mult = 0;
    gain_mult = 10;
    pwm_count_sec = 0;
    pwm_number = 0;
    pwm_remainder = 0;
    pwm_number_2 = 0;
    pwm_remainder_2 = 0;
    count_fix = 10;//14
    gain_coeff2 = 0;
    //filter_alpha = 0.002/8*66/10/16;// 10 - 1 KOm
    filter_alpha = 1.0/(0.01e-6*4e3*200e6*2*pi);// 10 - 1 KOm filter_alpha = 1.0/(0.01e-6*4e3*60e6*2*pi);
}
AGC::~AGC()
{
// test_out_dbg_1.close();
}/*
//AGC (Вычисление коэффициента усиления входного сигнала)
int AGC::GainCoefficient(short overflow_bit_in, short agc_mid_level_in, short overflow_percent)
{//vivod1<<overflow_bit<<endl;
 if(overflow_bit>0)			//-17 - порог включения режима "помеха пропала" можно изменять -19...-10 (1,7%)
 {								//20 - означает 2% превышения верхнего уровня по модулю при подсчёте за 1000 тактов
  overflow_adc=overflow_percent-overflow_bit;//-2
  gain_coeff+=overflow_adc*2;//*(-2)
  count1=10;
  //cout<<overflow_bit<<endl;
 }
 else
 {//cout<<overflow_bit<<" "<<agc_mid_level<<endl;
  s_zero=6;//6- для широкополосной фильтрации //8 - для узкополосной
  gain_coeff+=s_zero;
  if(agc_mid_level<=1500)
  { //cout<<count1<<endl;
   if(count1==9)		//count1==9 - для широкополосной фильтрации; count1=15 - для узкополосной
   {//cout<<"Jump"<<endl;
    if(agc_mid_level<=2) {gain_mult=16; delay_mult=11;}
    if(agc_mid_level>2 && agc_mid_level<=35) {gain_mult=12; delay_mult=11;}
    if(agc_mid_level>35 && agc_mid_level<=160) {gain_mult=8; delay_mult=11;}
    if(agc_mid_level>160 && agc_mid_level<=600) {gain_mult=6; delay_mult=11;}
    if(agc_mid_level>600) {gain_mult=5; delay_mult=11;}
    gain_coeff=(gain_coeff*gain_mult)>>2;
   }
  count1--;
  if(count1==-1) count1=10;//count1==10 - для широкополосной фильтрации;
  }
 }//vivod2<<(gain_coeff)<<endl;
 if(gain_coeff<63) gain_coeff=63;
 else if(gain_coeff>10240) gain_coeff=10240; //10240 - шир 6400 - узк
 //vivod1<<(gain_coeff>>6)<<endl;

 overflow_bit=overflow_bit_in;
 agc_mid_level=agc_mid_level_in;

 //vivod1<<(gain_coeff>>7)<<endl;
 return gain_coeff>>6; //6 - шир 7 - узк
}
*/
#define pwm_dac 1
int AGC::GainCoefficient(short overflow_bit_in, int adc_high_bit, int agc_mid_level_in, short overflow_percent)
{//vivod2<<agc_mid_level<<endl;
 if(adc_high_bit>0 && overflow_bit==0);
 else
 if(overflow_bit>0)			//-17 - порог включения режима "помеха пропала" можно изменять -19...-10 (1,7%)
 {								//20 - означает 2% превышения верхнего уровня по модулю при подсчёте за 1000 тактов
  overflow_adc=overflow_percent-overflow_bit;//-2

  //if(abs(overflow_adc)<30 && (overflow_adc)<0) overflow_adc = 0;
  //vivod1<<overflow_adc<<endl;
  //if(abs(overflow_adc)<30) overflow_adc = 0;

  gain_coeff+=overflow_adc*1;//*(-2)
  count1=count_fix;
  //cout<<overflow_bit<<endl;
 }
 else
 {//cout<<overflow_bit<<" "<<agc_mid_level<<endl;
  s_zero=6;//6- для широкополосной фильтрации //8 - для узкополосной
  gain_coeff+=s_zero;
 }
  //vivod1<<agc_mid_level<<' '<<-1500*8<<' '<<(agc_mid_level<=1500*8)<<endl;
  if(agc_mid_level<=1500*8)
  { //cout<<count1<<endl;
   if(count1 == count_fix-1)		//count1==9 - для широкополосной фильтрации; count1=15 - для узкополосной
   {//cout<<"Jump"<<endl;
    /*
	if(agc_mid_level<=2) {gain_mult=16; delay_mult=11;}
    if(agc_mid_level>2 && agc_mid_level<=35) {gain_mult=12; delay_mult=11;}
    if(agc_mid_level>35 && agc_mid_level<=160) {gain_mult=8; delay_mult=11;}
    if(agc_mid_level>160 && agc_mid_level<=600) {gain_mult=6; delay_mult=11;}
    if(agc_mid_level>600) {gain_mult=5; delay_mult=11;}
	gain_coeff=(gain_coeff*gain_mult)>>2;
    */

	if(agc_mid_level<=2*8) {gain_mult=24; delay_mult=11;}
    if(agc_mid_level>2*8 && agc_mid_level<=35*8) {gain_mult=16; delay_mult=11;}
    if(agc_mid_level>35*8 && agc_mid_level<=160*8) {gain_mult=12; delay_mult=11;}
    if(agc_mid_level>160*8 && agc_mid_level<=600*8) {gain_mult=10; delay_mult=11;}
    if(agc_mid_level>600*8) {gain_mult=9; delay_mult=11;}
	gain_coeff=(gain_coeff*gain_mult)>>3;

   }
  count1--;
  if(count1==-1) count1=count_fix;//count1==10 - для широкополосной фильтрации;
  //vivod1<<(gain_coeff)<<endl;
  }
 //vivod2<<(gain_coeff)<<endl;
 //if(gain_coeff<(64/1-1)) gain_coeff=(64/1-1);
 //if(gain_coeff<0) gain_coeff=0;
  //vivod1<<(gain_coeff)<<endl;
#if pwm_dac==0
 if(gain_coeff>(16383/1)) gain_coeff=16383/1; //10240 - шир 6400 - узк 16383 7232 10368
 if(gain_coeff<(64/1-1)) gain_coeff=(64/1-1);
#else
 if(gain_coeff>(16383/1)) gain_coeff=16383/1;
 if(gain_coeff<0) gain_coeff=0;
#endif
 //if(gain_coeff<63) gain_coeff=63;
 //else if(gain_coeff>10240) gain_coeff=10240; //10240 - шир 6400 - узк
 //vivod1<<(gain_coeff>>6)<<endl;

 overflow_bit=overflow_bit_in;
 agc_mid_level=agc_mid_level_in;

 //vivod1<<agc_mid_level<<endl;
 //vivod2<<gain_coeff<<endl;
#if pwm_dac==0
 return gain_coeff>>7; //6 - шир 7 - узк
#else
 return gain_coeff>>0; //6 - шир 7 - узк
#endif
}

int AGC::GainCoefficientAdv(short overflow_bit_in, int agc_mid_level_in, short overflow_percent, int parameter_0, int parameter_1)
{//vivod2<<agc_mid_level<<endl;
 if(overflow_bit>0)         //-17 - порог включения режима "помеха пропала" можно изменять -19...-10 (1,7%)
 {                              //20 - означает 2% превышения верхнего уровня по модулю при подсчёте за 1000 тактов
  overflow_adc=overflow_percent-overflow_bit;//-2

  //if(abs(overflow_adc)<30 && (overflow_adc)<0) overflow_adc = 0;
  //vivod1<<overflow_adc<<endl;
  //if(abs(overflow_adc)<30) overflow_adc = 0;

  gain_coeff+=overflow_adc*1;//*(-2)
  count1=count_fix;
  //cout<<overflow_bit<<endl;
 }
 else
 {//cout<<overflow_bit<<" "<<agc_mid_level<<endl;
  s_zero=6;//6- для широкополосной фильтрации //8 - для узкополосной
  //gain_coeff+=s_zero;
  gain_coeff+=(parameter_1-agc_mid_level_in)/parameter_0;//parameter_0 = 4096 parameter_1 = 131072
 }
  //vivod1<<agc_mid_level<<' '<<-1500*4<<' '<<(agc_mid_level<=-1500*4)<<endl;
  if(agc_mid_level<=1500*4)
  { //cout<<count1<<endl;
   if(count1 == count_fix-1)        //count1==9 - для широкополосной фильтрации; count1=15 - для узкополосной
   {//cout<<"Jump"<<endl;
    /*
    if(agc_mid_level<=2) {gain_mult=16; delay_mult=11;}
    if(agc_mid_level>2 && agc_mid_level<=35) {gain_mult=12; delay_mult=11;}
    if(agc_mid_level>35 && agc_mid_level<=160) {gain_mult=8; delay_mult=11;}
    if(agc_mid_level>160 && agc_mid_level<=600) {gain_mult=6; delay_mult=11;}
    if(agc_mid_level>600) {gain_mult=5; delay_mult=11;}
    gain_coeff=(gain_coeff*gain_mult)>>2;
    */
    if(agc_mid_level<=2*4) {gain_mult=24; delay_mult=11;}
    if(agc_mid_level>2*4 && agc_mid_level<=35*4) {gain_mult=16; delay_mult=11;}
    if(agc_mid_level>35*4 && agc_mid_level<=160*4) {gain_mult=12; delay_mult=11;}
    if(agc_mid_level>160*4 && agc_mid_level<=600*4) {gain_mult=10; delay_mult=11;}
    if(agc_mid_level>600*4) {gain_mult=9; delay_mult=11;}
    gain_coeff=(gain_coeff*gain_mult)>>3;
   }
  count1--;
  if(count1==-1) count1=count_fix;//count1==10 - для широкополосной фильтрации;
  //vivod1<<(gain_coeff)<<endl;
  }
 //vivod2<<(gain_coeff)<<endl;
 //if(gain_coeff<(64/1-1)) gain_coeff=(64/1-1);
 //if(gain_coeff<0) gain_coeff=0;
  //vivod1<<(gain_coeff)<<endl;
    //vivod1<<overflow_adc<<' '<<(5000-agc_mid_level_in)<<' '<<agc_mid_level_in<<endl;
#if pwm_dac==0
 if(gain_coeff>(16383/1)) gain_coeff=16383/1; //10240 - шир 6400 - узк 16383 7232 10368
 if(gain_coeff<(64/1-1)) gain_coeff=(64/1-1);
#else
 if(gain_coeff>(8191/1)) gain_coeff=8191/1;
 if(gain_coeff<0) gain_coeff=0;
#endif
 //if(gain_coeff<63) gain_coeff=63;
 //else if(gain_coeff>10240) gain_coeff=10240; //10240 - шир 6400 - узк
 //vivod1<<(gain_coeff>>6)<<endl;

 overflow_bit=overflow_bit_in;
 agc_mid_level=agc_mid_level_in;

 //vivod1<<agc_mid_level<<endl;
 //vivod2<<gain_coeff<<endl;
#if pwm_dac==0
 return gain_coeff>>7; //6 - шир 7 - узк
#else
 return gain_coeff>>1; //6 - шир 7 - узк
#endif
}
double AGC::log2(double a)
{
    return log(a) / log(double(2.));
}
int AGC::GainCoefficientLog(short overflow_bit_in, int agc_mid_level_in, short overflow_percent, int parameter_0, int parameter_1)
{//vivod2<<agc_mid_level<<endl;
    if(agc_mid_level_in==0)agc_mid_level_in = 1;
    gain_coeff2+=(log2(double(parameter_1/64))-log2(double(agc_mid_level_in)))/256;//parameter_0 = 4096 parameter_1 = 131072
    //vivod1<<log2(double(parameter_1))<<' '<<log2(double(agc_mid_level_in))<<' '<<log2(double(parameter_0))<<endl;
    gain_coeff=int(pow(2.0, gain_coeff2));
#if pwm_dac==0
    if(gain_coeff>(16383/1)) gain_coeff=16383/1; //10240 - шир 6400 - узк 16383 7232 10368
    if(gain_coeff<(64/1-1)) gain_coeff=(64/1-1);
#else
    if(gain_coeff>(8191/1)) gain_coeff=8191/1;
    if(gain_coeff<0) gain_coeff=0;
#endif
    overflow_bit=overflow_bit_in;
    agc_mid_level=agc_mid_level_in;
#if pwm_dac==0
    return gain_coeff>>7; //6 - шир 7 - узк
#else
    return gain_coeff>>1; //6 - шир 7 - узк
#endif
}
//------------------------------------------------------------------------------------------------------------------
int AGC::GainCoefficient_ALT3(short overflow_bit_in, int agc_mid_level_pre,short adc_low_in, long agc_mid_level_in, long ref_mid_level, int overflow_percent)
{//vivod2<<agc_mid_level<<endl;
    //double a = 1;
    double Dscr_out = 0;

    if(agc_mid_level_in==0) agc_mid_level_in=ref_mid_level/2;
    //double Dscr_filt_out_in = *Dscr_filt_out;
    overflow_bit=overflow_bit_in;
    //if((max_point < 450*4) && (ref_mid_level == 255*512)) ref_mid_level = 320*512;
    //if((max_point > 510*4) && (ref_mid_level == 320*512)) ref_mid_level = 255*512;
    //vivod1<<ref_mid_level<<' '<<max_point/4<<endl;
    //vivod1<<overflow_bit_in<<' '<<agc_mid_level_pre<<' '<<adc_low_in<<endl;
    if(overflow_bit_in==0 && agc_mid_level_pre>0);
    //else if((agc_mid_level_pre-overflow_bit_in)>=1 && overflow_bit_in>0)
    /*else if(overflow_bit_in<32 && overflow_bit_in>0)
    {
        gain_coeff += overflow_percent-overflow_bit_in;
        gain_coeff2 = log2(gain_coeff);//cout<<overflow_bit_in<<endl;
    }
    else if(adc_low_in>0 && overflow_bit_in==0)
    {
        s_zero=6;//6- для широкополосной фильтрации //8 - для узкополосной
        gain_coeff += s_zero;
        gain_coeff2 = log2(gain_coeff);//cout<<adc_low_in<<endl;
    }
    */else
    {
        Dscr_out=(log2(ref_mid_level)-log2(agc_mid_level_in))/4;     //250
        //Dscr_out=(log2(ref_mid_level)*8-log2(agc_mid_level_in))/4;
        //Dscr_out=(log2(ref_mid_level)-agc_mid_level_in);
        //vivod1<<overflow_bit_in<<' '<<agc_mid_level_pre<<endl;
        //vivod1<<agc_mid_level_in<<' '<<log2(ref_mid_level)*8<<endl;
        //vivod2<<Dscr_out<<endl;
        //vivod3<<gain_coeff2<<' '<<gain_coeff<<endl;
        //cout<<adc_low_in<<' '<<agc_mid_level_pre<<endl;
        //*Dscr_filt_out=a*Dscr_out + (1-a)*Dscr_filt_out_in;
        gain_coeff2+=Dscr_out;

        //if(gain_coeff2>13*512) gain_coeff2=13*512;
        if(gain_coeff2>13) gain_coeff2=13;
        //if(overflow_bit>0) gain_coeff+=overflow_adc<<loop_amplif;
        //else gain_coeff=pow(gain_coeff2,2.0);
        gain_coeff=int(pow(2.0, gain_coeff2));
    }
    //vivod1<<adc_low_in<<' '<<agc_mid_level_pre<<' '<<overflow_bit_in<<endl;
    //vivod2<<gain_coeff<<endl;
#if pwm_dac==0

    if(gain_coeff>16383) gain_coeff=16383;
    if(gain_coeff<64) gain_coeff=64;

#else
    if(gain_coeff>(8191/1)) gain_coeff=8191/1;
    if(gain_coeff<0) gain_coeff=0;
#endif
    //*gain_coeff_out = gain_coeff;
#if pwm_dac==0
    return gain_coeff>>7; //6
#else
    return gain_coeff>>1; //1
#endif
}

//----------------------------------------------------------------------------------------
const int N = 29;
const int INT_FACTOR = 10;
const double INT_MULT = static_cast<double>( 1 << INT_FACTOR );

const int TABLE_SIZE = 30;
int PosSteps[TABLE_SIZE];
int NegSteps[TABLE_SIZE];

inline int sign(int val)
{
    return (val >= 0) ? 1 : -1;
}

inline int ToInteger(double val)
{
    return static_cast<int>( floor( INT_MULT * val + 0.5f) );
}

inline double ToDouble(int val)
{
    return static_cast<double>(val) / INT_MULT;
}

//----------------------------------------------------------------------------------------

int AGC::Cordic_log2(int InValue_int)
{
    // Table forming
    for (int step=0; step<TABLE_SIZE; step++)
    {
        PosSteps[step] = ToInteger(log(1.0f + pow(2.0f,-step)) / log(2.0f));
        //vivod1<<PosSteps[step]<<endl;

        if ( step != 0)
        {
            NegSteps[step] = ToInteger(log(1.0f - pow(2.0f,-step)) / log(2.0f));
            //vivod2<<NegSteps[step]<<endl;
        }
    }

    // Reference
    //double ResValue = log( InValue ) / log(2.0f);

    // Initialization
    int X = InValue_int<<INT_FACTOR;
    int Y = 0;
    int Z = (1<<INT_FACTOR) - (InValue_int<<INT_FACTOR);
    int p = 0;
    int E = 1;
    //vivod3<<Z<<endl;
    // Calculation
    for (int step=0; step<=N; step++)
    {
        int Epr = E;
        E = sign(Z);
        if ( sign(Epr) != sign(E) ) { p++; }

        Z = Z - E*( X>>p );
        X = X + E*( X>>p );

        int index = p ;
        if ( E > 0 ) { Y = Y + PosSteps[index]; }
        else         { Y = Y + NegSteps[index]; }
        /*vivod1 << "p: " << p
             << "  Y: " << Y << "  (" << ToDouble(Y) << ")  "
             << "  Z: " << Z << "  (" << ToDouble(Z) << ")  "
             << "  X: " << X << "  (" << ToDouble(X) << ")  "
             << "  E: " << E << endl;*/
    }
    //vivod2 <<" stop "<<endl;
    return -Y;
}

int AGC::Cordic_pow2(int InValue_int)
{
    // Table forming
    for (int step=0; step<TABLE_SIZE; step++)
    {
        PosSteps[step] = ToInteger(log(1.0f + pow(2.0f,-step)) / log(2.0f));

        if ( step != 0)
        {
            NegSteps[step] = ToInteger(log(1.0f - pow(2.0f,-step)) / log(2.0f));
        }
        //vivod1<<PosSteps[step]<<' '<<NegSteps[step]<<endl;
    }

    // Reference
    //double ResValue = log( InValue ) / log(2.0f);

    // Initialization
    int X = 1;
    int Y = InValue_int;
    //int Z = static_cast<int>( floor(INT_MULT * 1 + 0.5f) ) - INT_MULT * InValue_int;
    int p = 0;
    int E = 1;
    //vivod3<<Z<<endl;
    // Calculation
    for (int step=0; step<=N; step++)
    {
        int Epr = E;
        E = sign(Y);
        if (E == 0) E = 1;
        if ( sign(Epr) != sign(E) ) { p++; }

        //Z = Z - E*( X>>p );
        X = X + E*( X>>p );

        int index = static_cast<int>( p );
        if ( E > 0 ) { Y = Y - PosSteps[index]; }
        else         { Y = Y - NegSteps[index]; }
        /*vivod2 << "p: " << p
             << "  i: " << InValue_int << "  (" << ToDouble(InValue_int) << ")  "
             << "  Y: " << Y << "  (" << ToDouble(Y) << ")  "
            //<< "  Z: " << Z << "  (" << ToDouble(Z) << ")  "
             << "  X: " << X << "  (" << ToDouble(X) << ")  "
             << "  E: " << E << endl;*/
    }//vivod2 <<" stop "<<endl;
    return X;
}

int AGC::GainCoefficient_Cordic(short overflow_bit_in, int agc_mid_level_pre,short adc_low_in, long agc_mid_level_in, long ref_mid_level, int overflow_percent, short k_cordic)
{
    int Dscr_out = 0;

    if(agc_mid_level_in==0) agc_mid_level_in=2;//agc_mid_level_in=ref_mid_level/2;
    overflow_bit=overflow_bit_in;
    if(overflow_bit_in==0 && agc_mid_level_pre>0);
    else
    {
        //Dscr_out=(log2(ref_mid_level)-log2(agc_mid_level_in))/4;     //250
        //Dscr_out=(ToDouble(Cordic_log2(ref_mid_level)) - ToDouble(Cordic_log2(agc_mid_level_in)))/4;
        //Dscr_out=ToDouble((Cordic_log2(ref_mid_level) - Cordic_log2(agc_mid_level_in))>>2);
        Dscr_out=(Cordic_log2(ref_mid_level) - Cordic_log2(agc_mid_level_in))>>k_cordic;
        //gain_coeff2+=double(Dscr_out);
        gain_coeff3+=Dscr_out;
        //vivod1<<ToDouble(Cordic_log2(agc_mid_level_in))<<' '<<log2(agc_mid_level_in)<<endl;
        //if(gain_coeff2>13*512) gain_coeff2=13*512;
        //if(gain_coeff2>13) gain_coeff2=13;
        if(gain_coeff3>(pow(2.,INT_FACTOR)*13)) gain_coeff3=int(pow(2.,INT_FACTOR)*13);
        //if(overflow_bit>0) gain_coeff+=overflow_adc<<loop_amplif;
        //else gain_coeff=pow(gain_coeff2,2.0);
        //gain_coeff=pow(2.0, gain_coeff2);
        //vivod1<<gain_coeff<<' '<<Cordic_pow2(ToInteger(gain_coeff2))<<endl;
        gain_coeff=Cordic_pow2(gain_coeff3);
    }
    //vivod1<<adc_low_in<<' '<<agc_mid_level_pre<<' '<<overflow_bit_in<<endl;
    //vivod3<<gain_coeff<<endl;
	//vivod1<<ref_mid_level<<' '<<agc_mid_level_in<<endl;
    //vivod2<<overflow_bit_in<<' '<<agc_mid_level_pre<<endl;
	//vivod3<<Cordic_log2(agc_mid_level_in)<<' '<<gain_coeff<<endl;
#if pwm_dac==0
    if(gain_coeff>16383) gain_coeff=16383;
    if(gain_coeff<64) gain_coeff=64;
#else
    if(gain_coeff>(8191/1)) gain_coeff=8191/1;
    if(gain_coeff<0) gain_coeff=0;
#endif
#if pwm_dac==0
    return gain_coeff>>7; //6
#else
    return gain_coeff>>1; //1
#endif
}

int AGC::GainCoefficient_CordicNew(short overflow_bit_in, int agc_mid_level_pre,short adc_low_in, long agc_mid_level_in, long ref_mid_level, int overflow_percent, short k_cordic)
{
    int Dscr_out = 0;

    if(agc_mid_level_in==0) agc_mid_level_in=2;//agc_mid_level_in=ref_mid_level/2;
    overflow_bit=overflow_bit_in;

    Dscr_out = Cordic_log2(ref_mid_level) - Cordic_log2(agc_mid_level_in);

    if(overflow_bit_in==0 && agc_mid_level_pre>0);
    else if(overflow_bit_in<5 && overflow_bit_in>0)
    {
        gain_coeff += overflow_percent-overflow_bit_in;
    }
    else if(adc_low_in>0 && overflow_bit_in==0)
    {
        s_zero=4;//6- для широкополосной фильтрации //8 - для узкополосной
        gain_coeff += s_zero;
    }
    else
    {
        if(overflow_bit_in>0 && overflow_bit_in<5) k_cordic = 32;
        else if(overflow_bit_in>=5 && overflow_bit_in<=37) k_cordic = (32-(overflow_bit_in-5))/2;
        else k_cordic = 0;
        if(overflow_bit_in==0) k_cordic = 2;

        //if(overflow_bit_in>0)
            gain_coeff += (Dscr_out>>k_cordic)*1 - overflow_bit_in;//72
        //else
            //gain_coeff += Dscr_out>>k_cordic;
    }
    //vivod1<<ref_mid_level<<' '<<agc_mid_level_in<<' '<<Dscr_out/4<<' '<<overflow_bit_in<<endl;
    //if(log_en==1) vivod1<<gain_coeff<<' '<<gain_coeff3<<' '<<Dscr_out/4<<' '<<overflow_bit_in<<endl;

    #if pwm_dac==0
    if(gain_coeff>16383) gain_coeff=16383;
    if(gain_coeff<64) gain_coeff=64;
    #else
    if(gain_coeff>(16383/1)) gain_coeff=16383/1;
    //if(gain_coeff>(16383/1)) gain_coeff=16383/1;
    if(gain_coeff<128) gain_coeff=128;
    #endif
    #if pwm_dac==0
    return gain_coeff>>7; //6
    #else
    return gain_coeff>>7; //1
    #endif
}

int AGC::GainCoefficient_CordicNew2(short overflow_bit_in, int agc_mid_level_pre,short adc_low_in, long agc_mid_level_in, long ref_mid_level, int overflow_percent, short k_cordic, short max_inp)
{
    int Dscr_out = 0;

    if(agc_mid_level_in==0) agc_mid_level_in=2;
    overflow_bit=overflow_bit_in;

    Dscr_out = Cordic_log2(ref_mid_level) - Cordic_log2(agc_mid_level_in);

    if(overflow_bit_in==0 && agc_mid_level_pre>0) gain_coeff2 = 0;
    else if(overflow_bit_in<5 && overflow_bit_in>0)
    {
        gain_coeff += (overflow_percent-overflow_bit_in);
        //gain_coeff2 += (overflow_percent-overflow_bit_in)/2;
        //gain_coeff += gain_coeff2;
    }
    else if(adc_low_in>0 && overflow_bit_in==0)
    {
        s_zero=4;
        gain_coeff += s_zero;
        //gain_coeff2 += 16;
        //gain_coeff += gain_coeff2;
    }
    else
    {
        if(overflow_bit_in>0 && overflow_bit_in<5) k_cordic = 32;
        else if(overflow_bit_in>=5 && overflow_bit_in<=37) k_cordic = 18-(overflow_bit_in-5)/2;
        else k_cordic = 2;
        if(overflow_bit_in==0) k_cordic = 4;

        if(overflow_bit_in>0)
        {
            gain_coeff += (Dscr_out>>k_cordic)*1 - overflow_bit_in;//72
            gain_coeff2 += ((Dscr_out>>k_cordic)*1 - overflow_bit_in)/2;
            gain_coeff += gain_coeff2;
        }
        //else if(adc_low_in>0
        else
        {
            //short max_discr = 9 - int(floor(log(double(max_inp))/log(2.))+1);
            //gain_coeff += 340*max_discr;
            //gain_coeff += ((Cordic_log2(510)-Cordic_log2(max_inp))>>3)*1;

            gain_coeff += (Dscr_out>>k_cordic)*1 - overflow_bit_in;//72
            gain_coeff2 += ((Dscr_out>>k_cordic)*1 - overflow_bit_in)/2;
            gain_coeff += gain_coeff2;

            /*
            gain_coeff += (Cordic_log2(490)-Cordic_log2(max_inp))>>5;//72
            gain_coeff2 += ((Cordic_log2(490)-Cordic_log2(max_inp))>>5)/8;
            gain_coeff += gain_coeff2;
            */
            //vivod1<<((Cordic_log2(490)-Cordic_log2(max_inp))>>5)<<' '<<(Dscr_out>>k_cordic)<<endl;
        }
    }

    //vivod1<<ref_mid_level<<' '<<agc_mid_level_in<<' '<<Dscr_out/4<<' '<<overflow_bit_in<<endl;
    //if(log_en==1) vivod1<<gain_coeff<<' '<<gain_coeff3<<' '<<Dscr_out/4<<' '<<overflow_bit_in<<endl;

    #if pwm_dac==0
    if(gain_coeff>16383) gain_coeff=16383;
    if(gain_coeff<64) gain_coeff=64;
    #else
    if(gain_coeff>(8191/1)) gain_coeff=8191/1;
    //if(gain_coeff>(16383/1)) gain_coeff=16383/1;
    if(gain_coeff<0) gain_coeff=0;
    #endif
    #if pwm_dac==0
    return gain_coeff>>7; //6
    #else
    return gain_coeff>>1; //1
    #endif
}

int AGC::GainCoefficient_ALT4(short overflow_bit_in, int agc_mid_level_pre, long agc_mid_level_in, long ref_mid_level, int max_point)
{//vivod2<<agc_mid_level<<endl;
    double a = 1;
    double Dscr_out = 0;

    //double Dscr_filt_out_in = *Dscr_filt_out;
    //vivod2<<agc_mid_level_in<<' ';
    if(agc_mid_level_in<=17 && agc_mid_level_in>1) {agc_mid_level_in = agc_mid_level_in*32-32; goto LB;}
    if(agc_mid_level_in<=92 && agc_mid_level_in>1) {agc_mid_level_in = agc_mid_level_in*4+466; goto LB;}
    if(agc_mid_level_in<=960 && agc_mid_level_in>1) {agc_mid_level_in = agc_mid_level_in/2+788; goto LB;}
    if(agc_mid_level_in<=18420 && agc_mid_level_in>1) {agc_mid_level_in = agc_mid_level_in/32+1238; goto LB;}
    if(agc_mid_level_in<=262144 && agc_mid_level_in>1) {agc_mid_level_in = agc_mid_level_in/512+1778; goto LB;}
    LB:
    if(agc_mid_level_in==0) agc_mid_level_in=1;
    ref_mid_level = 2033;
    //vivod2<<agc_mid_level_in<<endl;
    overflow_bit=overflow_bit_in;
    if(overflow_bit_in==0 && agc_mid_level_pre>0);
    else
    //Dscr_out=(log2(ref_mid_level)-log2(agc_mid_level_in));     //250
    //Dscr_out=(log2(ref_mid_level)*8-log2(agc_mid_level_in))/4;
    Dscr_out=(ref_mid_level-agc_mid_level_in);
    //vivod1<<ref_mid_level<<' '<<agc_mid_level_in<<endl;
    //vivod1<<Dscr_out<<' '<<gain_coeff2<<endl;
    gain_coeff2+=Dscr_out;
    //vivod3<<gain_coeff2<<' '<<gain_coeff<<endl;
    //if(gain_coeff2>13*512) gain_coeff2=13*512;
    if(gain_coeff2>13*512) gain_coeff2=13*512;
    //if(overflow_bit>0) gain_coeff+=overflow_adc<<loop_amplif;
    //else gain_coeff=pow(gain_coeff2,2.0);
    //gain_coeff=pow(2.0, gain_coeff2/512);
    if(gain_coeff2<=1362) {gain_coeff = int(gain_coeff2)/256+1; goto LBG;}
    if(gain_coeff2<=3040) {gain_coeff = int(gain_coeff2)/32-34; goto LBG;}
    if(gain_coeff2<=5199) {gain_coeff = int(gain_coeff2)/2-1459; goto LBG;}
    if(gain_coeff2<=6656) {gain_coeff = int(gain_coeff2)*5-24859; goto LBG;}
    LBG:


#if pwm_dac==0

 if(gain_coeff>16383) gain_coeff=16383;
 if(gain_coeff<64) gain_coeff=64;

#else
 if(gain_coeff>(8191/1)) gain_coeff=8191/1;
 if(gain_coeff<0) gain_coeff=0;
#endif

 //*gain_coeff_out = gain_coeff;

#if pwm_dac==0
 return gain_coeff>>7; //6
#else
 return gain_coeff>>6; //1
#endif
}
//! PWM (ШИМ)
#define PWM_advanced 1 //Enable advanced PWM scheme
//#define pwm_divide 16
//#define pwm_period 160
double AGC::PWM(int agc_gain_coeff,short pwm_divide, short pwm_period)
{
    //PWM algorithm
#if PWM_advanced == 1
    pwm_number=agc_gain_coeff/pwm_divide;
    pwm_remainder = agc_gain_coeff%pwm_divide;
    if(pwm_count == pwm_period/pwm_divide)
    {
        pwm_count = 0;
        pwm_count_sec++;
    }
    if(pwm_remainder > pwm_count_sec)  pwm_number += 1;
    if(pwm_count_sec == pwm_divide) pwm_count_sec = 0;
#else
    if(pwm_count == pwm_period) pwm_count = 0; //160 - шир 50 - узк
    pwm_number = agc_gain_coeff;
#endif
    if(pwm_count < pwm_number) pwm_out = 1;
    else pwm_out = 0;
    pwm_count++;//pwm_out=1;
    //vivod2<<pwm_out<<endl;
    //Smooth filter
    //pwm_filter_out=pwm_z*0.998+double(pwm_out)*0.002; //(0.998 0.002) - для широкополосной фильтрации	(0.9994 0.0006) - для узкополосной
    pwm_filter_out = pwm_z*(1-filter_alpha)+double(pwm_out)*filter_alpha;
    pwm_z = pwm_filter_out;
    //Nonlinear transformer
    //vivod1<<pwm_filter_out<<endl;
    pwm_filter_out = pow(28.0,pwm_filter_out*2.5)/4000;  //    (...)/4000 = (...)*5/20000
    //pwm_filter_out = pow(40.0,pwm_filter_out*10)/1.049e16;  //    (...)/4000 = (...)*5/20000
    //vivod1<<pwm_filter_out<<endl;
    return pwm_filter_out;
}

//! DAC
double AGC::DAC(int agc_gain_coeff,short mode)
{
    if(mode == 0) pwm_filter_out = pow(28.0,(double(agc_gain_coeff)/4095)*2.5)/4000/1;
    else pwm_filter_out = pow(28.0,(double(agc_gain_coeff/1024)/16)*1.6)/400/1;
	//vivod2<<pow(28.0,(double(agc_gain_coeff)/4095)*2.5)<<endl;
    //else pwm_filter_out = double(agc_gain_coeff/1024)*1.44/80;
    //pwm_filter_out = pwm_z*(1-filter_alpha*10)+pwm_filter_out*filter_alpha*10;
    //pwm_z = pwm_filter_out;
    //pwm_filter_out = pow(28.0,(double(agc_gain_coeff)/4095)*2.5)/4000;
    //vivod2<<pwm_filter_out<<endl;

    return pwm_filter_out;
}

double AGC::DACnew(int agc_gain_coeff,short mode)
{
    if(mode == 0) pwm_filter_out = agc_gain_coeff*0.0546+0.06;
        else pwm_filter_out = agc_gain_coeff*0.3906+0.4;

    return pwm_filter_out/200;
}

//! Secondary quantization
short AGC::AGC_second(int filtered_signal_I)
{
    //! Filter  RC_tau=5, mult_coeff=2;
    y1=y;
    x=int(abs(filtered_signal_I))<<mult_coeff;
    y2 += (y1&0x1F);
    y=y1-(y1>>RC_tau)+(x>>RC_tau);// - (y2>>RC_tau);
    if(y2>=32) y2 -= 32;
    //y = 120;
    //Tuning of AGC levels
    //y=y>>mult_coeff;
    L1=y<<1;
    // ***************** y - dbg
    //test_out_dbg_1 << y << endl;

    L2=y+(y>>2)-(y>>4)+(y>>6);
    L3=(y>>1)-(y>>3);

    L1=L1>>mult_coeff;
    L2=L2>>mult_coeff;
    L3=L3>>mult_coeff;

    //test_out_dbg_2 << L1 << endl;

    //test_out_dbg_3 << L2 << endl;
    //test_out_dbg_4 << L3 << endl;


    /*
    L1=y<<1;
    // ***************** y - dbg
    test_out_dbg_1 << y << endl;
    L2=y+(y>>2)-(y>>4)+(y>>6);
    L3=(y>>1)-(y>>3);
    */

    //! Quantization
    last_out=0;
    if(filtered_signal_I>=0)
    {
        if(filtered_signal_I>L3) last_out=1;
        if(filtered_signal_I>L2) last_out=2;
        if(filtered_signal_I>L1) last_out=3;
        //last_out = 1;
    }
    else
    {
        if(filtered_signal_I<(-1*L3)) last_out=-1;
        if(filtered_signal_I<(-1*L2)) last_out=-2;
        if(filtered_signal_I<(-1*L1)) last_out=-3;
        //last_out = -1;
    }

    //vivod1<<filtered_signal<<'\t'<<y<<'\t'<<last_out<<endl;
    //vivod1<<filtered_signal_I<<endl;
    //vivod2<<filtered_signal_I<<endl;
    //vivod3<<L1<<' '<<L2<<' '<<L3<<endl;
    return last_out;
}

L5_amp_ctrl::L5_amp_ctrl()
{
 vivod1.open("data_out/f1.txt",ios::out);
 vivod2.open("data_out/f2.txt",ios::out);
 vivod3.open("data_out/f3.txt",ios::out);

 amp_coeff=2000;
 discrim=0;

 k_agc=1;
 pwm_count=0;
 pwm_out=0;
 pwm_filter_out=0;
 pwm_z=0;
 last_out=0;
 delay_mult=0;
 gain_mult=1;
 pwm_count_sec=0;
 pwm_number=0;
 pwm_remainder=0;
 pwm_number_2=0;
 pwm_remainder_2=0;
}
L5_amp_ctrl::~L5_amp_ctrl()
{
}
//AGC (Вычисление коэффициента усиления входного сигнала)
int L5_amp_ctrl::GainCoefficient(int L0_sum, int L5_sum)
{
 L0_sum/=16384;
 L5_sum/=16384;
 discrim=L0_sum*2-L5_sum;
 amp_coeff+=(discrim/1);

 //if(amp_coeff<63) amp_coeff=63;
 //else if(amp_coeff>16240) amp_coeff=16240; //9600 - шир 6400 - узк
 //vivod1<<amp_coeff<<endl;
 return amp_coeff;
}

// PWM (ШИМ)
//#define PWM_advanced 1 //Enable advanced PWM scheme
//#define pwm_divide 16
//#define pwm_period 160
double L5_amp_ctrl::PWM(int agc_gain_coeff,short pwm_divide, short pwm_period)
{
 //PWM algorithm
#ifdef PWM_advanced
 pwm_number=((agc_gain_coeff>>2)/pwm_divide);
 pwm_remainder=((agc_gain_coeff>>2)%pwm_divide);

 //pwm_number_2=pwm_remainder/8;
 //pwm_remainder_2=pwm_remainder%8;
 //if(pwm_remainder<8) pwm_remainder_2=pwm_remainder;
 if(pwm_count==pwm_period/pwm_divide)
 {
  pwm_count=0;
  pwm_count_sec++;
 }

 //if(pwm_count_sec%2==0 && pwm_number_2==1) pwm_number+=1;
 //if(pwm_count_sec%2==0 && pwm_remainder_2*2>pwm_count_sec ) pwm_number+=1;

 if(pwm_remainder>pwm_count_sec)  pwm_number+=1;
 if(pwm_count_sec==pwm_divide) pwm_count_sec=0;

#else
 if(pwm_count==pwm_period) pwm_count=0; //150 - шир 50 - узк
 pwm_number=((agc_gain_coeff>>6));
#endif
 //vivod1<<pwm_number<<endl;
 //vivod2<<pwm_remainder<<endl;
 //vivod3<<(agc_gain_coeff>>6)<<endl;
 if(pwm_count<pwm_number) pwm_out=1; //6 - шир 7 - узк
 else pwm_out=0;
 pwm_count++;//pwm_out=1;
 //Smooth filter
 pwm_filter_out=pwm_z*0.999+double(pwm_out)*0.001; //(0.998 0.002) - для широкополосной фильтрации	(0.9994 0.0006) - для узкополосной
 pwm_z=pwm_filter_out;
 //Nonlinear transformer
 pwm_filter_out=pwm_filter_out*7;  //    (...)/4000 = (...)*5/20000
 //vivod2<<pwm_z<<endl;
 return pwm_filter_out;
}
