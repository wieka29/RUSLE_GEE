# -*- coding: utf-8 -*-
import numpy as np
from netCDF4 import Dataset
from sys import argv

#code Adjusted RUSLE model as presented in Naipal et al., 2015., 2016

def K_factor(fname,llx,urx,lly,ury,lat_h,lon_h): 

    '''
    This function calculates the K factor based on the method of Torri et al., 1997 
    at 1 km spatial resolution and then upscales the result to 30m

    '''
    from Input_File import input_K
    from scaling_func import upscaling
    #read data
    sand_t,silt_t,clay_t,oc_t,su_code90,lat_c,lon_c = input_K(fname,llx,urx,lly,ury) 
    #calculate Dg & K
    Dg = (sand_t*np.log((2.*0.05)**0.5))+(silt_t*np.log((0.05*0.002)**0.5))+(clay_t*np.log((0.002*0.00005)**0.5))
    K = 0.0293*(0.65-Dg+0.24*Dg**2.)*np.exp((-0.0021*(oc_t/clay_t))-(0.00037*(oc_t/clay_t)**2.)-(4.02*clay_t)+(1.72*clay_t**2.))  
    #Volcanic soils K = 0.08 (Andisols)
    #volcanic soils are andisols numbers 20-27 from Taxousda
    su_code90[su_code90>28]=0;su_code90[su_code90<20]=0
    k_volc=np.zeros((sand_t.shape))
    k_volc[su_code90!=0]=0.08
    k_volc[su_code90==0]=K[su_code90==0]
    #rescale K factor to 30m
    r_sel,k_sel=upscaling(lat_h,lon_h,lat_c,lon_c)
    K_out1=k_volc[r_sel[:]];K_out=K_out1[:,k_sel[:]]
    
    return K_out
  
def C_factor(fname,year_start,year_end,llx,urx,lly,ury):
    
    '''
    This function calculates the C factor based on ndvi at monthly resolution 
    and adjusts maximum values for forest and grassland 
    '''
    from Input_File import input_lc,input_C
    from scaling_func import upscaling
    #read data
    ndvi=input_C(fname,year_start,year_end,llx,urx,lly,ury)
    treefrac,grassfrac,cropfrac=input_lc(fname,llx,urx,lly,ury)
    rows=treefrac.shape[0];kols=treefrac.shape[1]   
    #calculate monthly C factor 
    ndvi[ndvi<0.] = np.nan
    ndvi=ndvi.reshape(12,rows*kols)
    c = np.exp(-2*(ndvi/(1-ndvi)))
    tree_h = np.ravel(treefrac)
    grass_h = np.ravel(grassfrac)
    #define maximum c values for grass and tree
    for i in range(12):
        c_temp_tree=np.zeros((rows*kols))
        c_temp_tree[tree_h>0.7]=c[i,tree_h>0.7]
        c_temp_tree[c_temp_tree>0.01]=0.01
     
        c_temp_grass=np.zeros((rows*kols))
        c_temp_grass[grass_h>0.7]=c[i,grass_h>0.7]
        c_temp_grass[c_temp_grass>0.05]=0.05
        
        c[i,tree_h>0.7]=c_temp_tree[tree_h>0.7]
        c[i,grass_h>0.7]=c_temp_grass[grass_h>0.7]
               
    C_all = np.asarray(c).reshape(12,rows,kols)
    C_all[C_all==0.]=np.nan
    C_out=np.nanmean(C_all,axis=0)
    C_out[np.isnan(C_out)==True]=0.
    
    return C_out

def LS_factor(fname,llx,urx,lly,ury):
    
    '''
    This function calculates the S factor based on Nearing (1997) and 
    the L factor is derived with the method of Desmet and Govers (1996)
    '''
    
    from Input_File import input_LS
    #read data
    slope,aspect,A_contr=input_LS(fname,llx,urx,lly,ury) #in degrees and m2
    rows=slope.shape[0]
    kols=slope.shape[1]
    #make correction for A_contr shape
    A_contr=A_contr[:-1,:]
    A_contr[A_contr>10000.]=10000. #max slope length for RUSLE is 305m
    #calculate slope length
    slope_rad=np.radians(slope)
    aspect_rad=np.radians(aspect)
    F=(np.sin(slope_rad)/0.0896)/(3*np.sin(slope_rad)**0.8+0.56) #slope in degrees
    m=F/(1+F)
    X=np.sin(aspect_rad)+np.cos(aspect_rad)
    D=30. 
    L=((A_contr+D**2)**(m+1)-A_contr**(m+1))/(D**(m+2)*X**m*22.13**m)
    L[L<0.] = 0.;
    L[np.isnan(L)==True] = 0.
    # calculate S factor
    S=-1.5+(17./(1+np.exp(2.3-6.1*np.sin(slope_rad)))) #slope in degrees
    S[S<0.] = 0.
    S[np.isnan(S)==True] = 0.
    LS_out=S*L
    
    return LS_out
            
def R_factor(fname,llx,urx,lly,ury,lat_h,lon_h):
    
    '''
    This function calculates the erosivity factor at 5 arcmin resolution using the erosivity.f90 script
    erosivity.f90 has to be wrapped to be able to import as a function in python using f2py
    Before executing this function run in your terminal: f2py -c -m erosivity erosivity.f90
    '''
    
    from Input_File import input_R
    import sys
    sys.path.insert(0,'/home/wieka20/Dokumente/Adj_RUSLE')
    import erosivity     
    from scaling_func import upscaling
    #read data input for erosivity
    z,p,k,sdii,rows_sel,kols_sel,lat_c,lon_c=input_R(fname,llx,urx,lly,ury)  
    #calculate R factor at 10km res
    r=erosivity.erosivity.main(z,p,k,sdii,len(rows_sel),len(kols_sel)) #shape(kols*rows)
    r[r==-9999.]=np.nan
    R=r.reshape(len(rows_sel),len(kols_sel))
    #rescale R factor to 30m
    r_sel,k_sel=upscaling(lat_h,lon_h,lat_c,lon_c)
    R_out1=R[r_sel[:]];R_out=R_out1[:,k_sel[:]]
    
    return R_out

def P_factor(fname,llx,urx,lly,ury):
    
    '''
    This function calculates the P factor for areas under terracing according to Haidong et al, 2016
    '''
    from Input_File import input_P
        
    terrace_frac=input_P(fname,llx,urx,lly,ury)
    P_out =0.12*terrace_frac
    P_out[P_out==0.]=1. #if no terracing then no reduction in erosion
    
    return P_out
    

