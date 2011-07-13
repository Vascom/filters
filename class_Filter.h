/**************************************************************************** \
*                   Copyright (c) 2006 by Javad GNSS.                        *
*                          All Rights Reserved.                              *
******************************************************************************

Data : 01.03.2006

Author : Glazov Vasiliy

Contents: declaration of classes for filter modeling

Description: these classes contain interface data and methods for filter modeling

\****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <numeric>
#define V5_COMPAT

using namespace std;
//const double pi=3.14159265358979323846;//26433832795;
const double pi=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067;
//const double pi=3.141592653589792;
const double pi_2=6.283185307179586476925286766559;

class Filter
{
    short *h1,*h1_pos,*h2,h2_dop,N1,N1_pos,N2,Nsat,i,sw1,N_aj,aj_delay,aj_d_delay,N_aj_adv,aj_delay_adv;
    int p1,*input1_rev,*input1_rev_pos,*inputI,*inputI_pos,*input1_revH,*inputH,*inputAJ_rev,*inputAJ,*inputAJ2,yk,yk1,yk2,yk3,yk4,ykk,ykk1,ykk2,ykk3,ykk4,eps,eps1,eps2,eps3,eps4,eps5,mju,epss,epss1,epss2,epss3,epss4,epss5;
    short *w_aj,*w_aj_abs,*w_ajj,*w_aj1,*w_aj2,*w_aj11,*w_aj22;
    double *inputAJ_rev_d,eps_d,*inputAJ_d,yk_d,*w_aj_d,mju_d,aj_in_z,*input_delay_d;
    double *inputAJ_rev_d_cplx_i,*inputAJ_d_cplx_i,yk_d_cplx_i,*w_aj_d_cplx_i;
    double *inputAJ_rev_d_cplx_q,*inputAJ_d_cplx_q,yk_d_cplx_q,*w_aj_d_cplx_q,eps_d_cplx_i,eps_d_cplx_q;
    double psat,*inputSat,*inputSat_rev,*hsat;

    short  Naj,Naj_delay;
    double aj_comb_d_i,aj_comb_d_q,rn,rnn[8];
    double *w_aj_d_i,*w_aj_d_q,*input_d_i,*input_d_q;
    double *r_i,*r_q,*r_i_n,*r_q_n;
    short *hc_aj;
    short  input_bw, fltcoeff_bw, multround, multround_bw, sumround, sumround_bw, coeffcompmult_bw, coeffcompldel, coeffcompdiv_bw;

    short *input_i,*input_q,aj_comb_i[8],aj_comb_q[8];
    int *r_i_s,*r_q_s,*r_i_n_s,*r_q_n_s,*w_aj_i,*w_aj_q,*w_lms_i,*w_lms_q,*waj_res_i,*waj_res_q,*w_lms,*w_lms1,*w_lms2,*waj_res,*waj_res1,*waj_res2;
    short *aj_cplx_const_i,*aj_cplx_const_q;

    int r_i_s_46,r_q_s_46,eps_i1,eps_q1,eps_i2,eps_q2;
    int sas1,sas2,sas3,sas4;

    double ComplexConst_i(short);
    double ComplexConst_q(short);

    static const short cic_stage_number=6, cic_stage_max_length=17;

    int cic_out,cic_delay_matrix[cic_stage_number];
    int stage_int[cic_stage_number],stage_comb[cic_stage_number],stage_comb_z[cic_stage_number][cic_stage_max_length];
    int sdvig;
    //For spectrum analysers CIC-filter
    int sa_cic_out, sa_stage_int[3], sa_stage_comb[3], sa_stage_comb_z[3][513];
    static const short sa_st0=256, sa_st1=197, sa_st2=121;// 512 395 245

    int sec_agc_sum, sec_agc_sum_pre;
    int f_sec_agc_old, f_sec_agc, f_sec_agc_tr, f_sec_agc_tr2, f_sec_agc_tr3, agcmidout, max_fast[5], max_fast2[5], max_slow, max_slow2;
    int nco_acc;

    short flt_tbl[4096];
    short filt_cic_out_shift,mid_agc_lev_n,mid_agc_lev_p,
        mid_agc_count_n,mid_agc_count_p,mid_agc_man;

    short lms_mju;
    //---------------------------------
    int iira[7][3], iirb[7][3], iira1[3], iirb1[3], iira2[3], iirb2[3];
    int rein[8][4], reout[4], reout1[4], reout2[4];
    double iirad[4], iirbd[4], iira1d[3], iirb1d[3], iira2d[3], iirb2d[3], iira3d[3], iirb3d[3];
    double reind[4], reoutd[4], reout1d[4], reout2d[4], reout3d[4];
    double *locata_pattern;

    double ave_input1_z, ave_input2_z;

    short geng1[14];
    short nco_out_d1, nco_out_d2;
    short k_nco_locata,k_nco_locata2, k_nco_locata_lock;
    double signal_filtered_line[256],signal_filtered_line_real[256];
    short f_sec_agc_locata,f_sec_agc_locata_tr,f_sec_agc_locata_tr2;
    int L1,L2,L3,Ln, max_locata;

    int *fft_mux_a1_pre,*fft_mux_a2_pre,*fft_mux_a1,*fft_mux_a2,**fft_mux_sw,**fft_mux_p1,**fft_mux_p2;
    int c1_fft,c2_fft,d1_fft,d2_fft;
    int *fft_mux_a1_s,*fft_mux_a2_s,*fft_mux_b1_s,*fft_mux_b2_s,*fft_mux_sw_rev;
    double cos_sin_cplx_re,cos_sin_cplx_im;

public:
    double eps_d_i,eps_d_q;
    int eps_i,eps_q, agcmidout_q;
    short cos_out, sin_out;
    short code_e, code_p, code_l, strb_h;
    long long nco_reg;
    double locata_out;

    fstream test_out_dbg_1,test_out_dbg_2,test_out_dbg_3,test_out_dbg_4;
    fstream test_out_dbg_5;//,test_out_dbg_6,test_out_dbg_7,test_out_dbg_8
    fstream test_out_3;
    fstream vivod1,vivod2,vivod3;

    int keep_sign(int num, int width);
    Filter();
    ~Filter();
    short NoiseGenLocata();
    short NCO_code_gen(int);
    short NCO_code_gen_phase_gate(int,int,short,int);
    short NCO(short,short);
    short NCO_real_cos(short,long long);
    void FilterReset(short,short[]);
    void LoadCoeff(short,short[],short,short,short,short,short,short,short,short[],short[],short[]);
    int IirI(int,short,short);
    double IirIn(double);
    int FirI(int,short,short);
    int FirIpos(int,short,short);
    int FirIn(int,short,short);
    int FirH(int,short,short,short);
    double SatFilter(double);
    void LoadSatCoeff(short,char *);
    int CICFilter(int);
    short AGCMid(int, long long);
    int AGCMids(int,int,int,int, long long);
	int AGCMidsPos(int,int,int,int, long long);
    short AGCMid2(int,int,int,int, long long);
    void LoadAJFilter(short,short,short,short);
    short AJFilter(short,short,short,short,short);
    short AJFilterNew(short,short,short,short,short);
    int AJFilterNew2(int,short,short,int[]);
    int AJFilterNew3(int,int,short,short,int[]);
    int AJFilterNew2_cplx(int,int,short,short,short,int[]);
    int AJFilterNew2_cplxs(int,int,short,short,short,int[]);
    double AJFilter_double(short,short,short,short,short);
    double AJFilter_double_cplx(short,short,short,short);
    short AJFilter_cplx(short,short,short,short);

    double AJFilter_adv_double(double,double,short,short);
    short AJFilter_adv(short,short,short,short,short,short[]);
    void LoadAJFilter_adv(short,short,short[]);
    int SpectrAn(int, short);

    short AGCLocata(short, long long);
    short AGCSecondLocata(short, long long, short, short, double, short);

	//double *FFTw(int[],int);
    int FFT(int[],int[],short);
    int FFTb(int,int,int,int, int,int);
    int GetFFT(short,short);
    int GetFFTb(short);
    int cos_sin_cplx(double,double,double,double);
    double Get_cos_sin_cplx(short);
};


class AGC
{
    short overflow_adc,count1,s_zero,gain_coeff_second,k_agc,pwm_count,pwm_out,last_out,delay_mult,gain_mult,count_fix;
    short pwm_count_sec,pwm_number,pwm_remainder,pwm_number_2,pwm_remainder_2;
    int x,y,y1,y2,L1,L2,L3,gain_coeff3;
    double pwm_filter_out,pwm_z, filter_alpha, gain_coeff2;

    int overflow_bit,agc_mid_level;
    //Static parameters
    static const short RC_tau=3, mult_coeff=2;
public:
    fstream vivod1,vivod2,vivod3;
    int gain_coeff;
    short last_out_Q;
    AGC();
    ~AGC();
    int GainCoefficient(short,int,int,short);
    int GainCoefficientAdv(short,int,short,int,int);
    int GainCoefficientLog(short,int,short,int,int);
    int GainCoefficient_ALT3(short,int,short,long,long,int);
    int GainCoefficient_ALT4(short,int,long,long,int);
    int GainCoefficient_Cordic(short,int,short,long,long,int,short);
    int GainCoefficient_CordicNew(short,int,short,long,long,int,short);
    int GainCoefficient_CordicNew2(short,int,short,long,long,int,short,short);
    int Cordic_log2(int);
    int Cordic_pow2(int);

    double log2(double);
    double PWM(int,short,short);
    double DAC(int,short);
    double DACnew(int,short);
    short AGC_second(int);
};

class L5_amp_ctrl
{
    short a,count1,s_zero,gain_coeff_second,k_agc,pwm_count,pwm_out,last_out,delay_mult,gain_mult;
    short pwm_count_sec,pwm_number,pwm_remainder,pwm_number_2,pwm_remainder_2;
    int amp_coeff,discrim;
    double pwm_filter_out,pwm_z;
    public:
    fstream vivod1,vivod2,vivod3;
    L5_amp_ctrl();
    ~L5_amp_ctrl();
    int GainCoefficient(int,int);
    double PWM(int,short,short);
};