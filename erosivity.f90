module erosivity
  implicit none
  
  contains
  
  subroutine main(z,p,k,sdii,rows,kols,r_out) !input should be p[t],sdii[t],k[t]; -9999 for r_out should be defined as nan in python after
    !use ieee_arithmetic
    
    implicit none
    
    real(8),dimension(kols*rows),intent(in)::sdii,p,z
    integer,intent(in):: kols, rows
    integer:: i   
    integer,dimension(kols*rows),intent(in)::k
    real(8),dimension(kols*rows),intent(out):: r_out 
    
    do i=1,size(z)
      if (p(i)==-9999.) then
          r_out(i)=-9999.
      else
          if (k(i)==5) then!BWk
              if (sdii(i)==-9999.) then
                  if (p(i)<=850.) then
                      r_out(i)=0.0483*p(i)**1.61
                  else
                      r_out(i)=587.8-1.219*p(i)+0.004105*p(i)**2
                  endif
              else
                  r_out(i)=0.809*p(i)**0.957 + 0.000189*sdii(i)**6.285
              endif
          elseif (k(i)==6) then!BSh
              if (sdii(i)==-9999.) then
                  r_out(i)=exp(-8.164+2.455*log(p(i)))
              else
                  r_out(i)=exp(-7.72+1.595*log(p(i))+2.068*log(sdii(i)))
              endif
          elseif (k(i)==7) then!BSk
              if (z(i)>=200.) then
                  if (sdii(i)==-9999.) then
                      r_out(i)=exp(5.52+1.33*log(p(i))-(0.977*log(z(i))))
                  else
                      r_out(i)=exp(0.0793+0.887*log(p(i)) + 1.892*log(sdii(i))-(0.429*log(z(i))))
                  endif
              else
                  if (sdii(i)==-9999.) then
                      r_out(i)=exp(5.52+1.33*log(p(i))-(0.977*log(200.)))
                  else
                      r_out(i)=exp(0.0793+0.887*log(p(i)) + 1.892*log(sdii(i))-(0.429*log(200.)))
                  endif
              endif
          elseif (k(i)==9) then!Csb
              r_out(i)=98.35+0.0003549*p(i)**1.987
          elseif (k(i) == 14) then!Cfa
              if (z(i) >= 4.5) then
                  if (sdii(i)==-9999.) then
                      r_out(i)=exp(3.378+0.852*log(p(i))-(0.191*log(z(i))))
                  else
                      r_out(i)=exp(0.524+0.462*log(p(i))-(0.106*log(z(i)))+1.97*log(sdii(i)))
                  endif
              else
                  if (sdii(i)==-9999.) then
                      r_out(i)=exp(3.378+0.852*log(p(i))-(0.191*log(4.5)))
                  else
                      r_out(i)=exp(0.524+0.462*log(p(i))-(0.106*log(4.5))+1.97*log(sdii(i)))
                  endif
              endif
          elseif (k(i)==15) then!Cfb
              if (sdii(i)==-9999.) then
                  if (z(i) >= 75.) then
                      r_out(i)=exp(5.267+0.839*log(p(i))-(0.635*log(z(i))))
                  else
                      r_out(i)=exp(5.267+0.839*log(p(i))-(0.635*log(75.)))
                  endif
              else
                  r_out(i)=exp(-4.853+0.676*log(p(i))+3.34*log(sdii(i)))
              endif
          elseif (k(i) == 17) then!Dsa
              if (sdii(i)==-9999.) then
                  r_out(i)=exp(7.49+-0.0512*log(p(i))-(0.272*log(z(i))))
              else
                  r_out(i)=exp(8.602+-0.963*log(sdii(i))-(0.247*log(z(i))))
              endif
          elseif (k(i) == 18) then!Dsb
              r_out(i)=exp(2.166+0.494*log(p(i)))
          elseif (k(i) == 19) then!Dsc
              if (sdii(i)==-9999.) then
                  r_out(i)=exp(4.416+0.0594*log(p(i)))
              else
                  r_out(i)=exp(6.236-(0.869*log(sdii(i))))
              endif
          elseif (k(i) == 21) then!Dwa
              r_out(i)=exp(-0.572+1.238*log(p(i)))
          elseif (k(i) == 22) then!Dwb
              if (sdii(i)==-9999.) then
                  r_out(i)=exp(1.882+0.819*log(p(i)))
              else
                  r_out(i)=exp(-1.7+0.788*log(p(i))+1.824*log(sdii(i)))
              endif
          elseif (k(i) == 25) then!Dfa
              if (sdii(i)==-9999.) then
                  r_out(i)=exp(-2.396+1.5*log(p(i)))
              else
                  r_out(i)=exp(-1.99+0.737*log(p(i))+2.033*log(sdii(i)))
              endif
          elseif (k(i) == 26) then!Dfb
              if (z(i) >= 65.) then
                  if (sdii(i)==-9999.) then
                      r_out(i)=exp(1.96+1.084*log(p(i))-(0.34*log(z(i))))
                  else
                      r_out(i)=exp(-0.5+0.266*log(p(i))-(0.131*log(z(i)))+3.1*log(sdii(i)))
                  endif
              else
                  if (sdii(i)==-9999.) then
                      r_out(i)=exp(1.96+1.084*log(p(i))-(0.34*log(65.)))
                  else
                      r_out(i)=exp(-0.5+0.266*log(p(i))-(0.131*log(65.))+3.1*log(sdii(i)))
                  endif
              endif
          elseif (k(i) == 27) then!Dfc
              if (sdii(i)==-9999.) then
                  r_out(i)=exp(-3.263+1.576*log(p(i)))
              else
                  r_out(i)=exp(-1.259+3.862*log(sdii(i)))
              endif
          elseif (k(i) == 29) then!ET
              r_out(i)=exp(-3.945+1.54*log(p(i)))
          elseif (k(i) == 30) then!EF
              r_out(i)=exp(16.39-1.286*log(p(i)))
          elseif (k(i) == 31) then!ETH
              if (sdii(i)==-9999.) then
                  r_out(i)=exp(-10.66+2.43*log(p(i)))
              else
                  r_out(i)=exp(21.44+1.293*log(p(i))-(10.579*log(sdii(i))))
              endif
          elseif (k(i) == 32) then!EFH
              r_out(i)=exp(16.39-1.286*log(p(i)))
          else !Af(1),Am(2),Aw(3),BWh(4),Csa(8),Csc(10),Cwa(11),Cwb(12),Cwc(13),Cfc(16),Dsd(20),Dwc(23),Dwd(24),Dfd(28)
              if (p(i) <= 850.) then
                  r_out(i)=0.0483*p(i)**1.61
              else
                  r_out(i)=587.8-1.219*p(i)+0.004105*p(i)**2
              endif
          endif
      endif
    enddo
  end subroutine main
end module erosivity

 !f2py -c -m routing_sub routing_sub.f90
