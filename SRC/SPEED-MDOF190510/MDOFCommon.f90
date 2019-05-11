!!!Some useful subroutines for time integration computation. by Group of XZ LU
subroutine Central_Difference(M, C, M_Inv, dT, U, U1, U0, P, NDOF, NewOrNot) ! Central_Difference
	implicit none;
	integer NDOF; 
	integer NewOrNot;
	real*8 M(NDOF,NDOF),C(NDOF,NDOF); 
	real*8 M_Inv(NDOF, NDOF)
	real*8 U(NDOF), U1(NDOF), U0(NDOF), P(NDOF), dT; 

	integer I, J; 
	real*8 M1(NDOF,NDOF), M2(NDOF,NDOF), M3(NDOF,NDOF) 
	real*8 M4(NDOF,NDOF), x, y, T(NDOF) 

	M1=M/dT/dT+C/2./dT; 
	M2=-2.*M/dT/dT;     
	M3=M/dT/dT-C/2./dT;

	if(NewOrNot .eq. 0) then
	    call matinv(M1, M4, NDOF); 
	    M_Inv = M4
		NewOrNot=1;
	else
	    M4 = M_Inv
	end if

	do I=1, NDOF;
		x=0.;	y=0.;
		do J=1, NDOF;
			U(J)=0;	x=x+M2(I,J)*U1(J);	y=y+M3(I,J)*U0(J);
		end do; 
		T(I)=P(I)-x-y;
	end do
	do I=1, NDOF
		do J=1, NDOF; U(I)=U(I)+M4(I,J)*T(J); end do;
	end do
	
	return;
end subroutine Central_Difference



subroutine matinv(A, B, N) 
	implicit none
	integer N
	real(8) ,intent (in)::A(N,N)
	real(8) ::B(N,N)
	integer(4):: I,J,K;	real(8)::D,T
	integer(4)::IS(N),JS(N)		
	B=A
	do 	K=1,N
		D=0.0D0
		do I=K,N
			do J=K,N
				if(abs(B(I,J))>D) then
					D=abs(B(I,J)); IS(K)=I; JS(K)=J
				end if
			end do
		end do
		do J=1,N
			T=B(K,J); B(K,J)=B(IS(K),J); B(IS(K),J)=T
		end do
		do I=1,N
			T=B(I,K);	B(I,K)=B(I,JS(K));	B(I,JS(K))=T
		end do
		B(K,K)=1.d0/B(K,K);
		do J=1,N; if(J.NE.K) B(K,J)=B(K,J)*B(K,K);  end do
		do I=1,N
			if(I.NE.K) then
				do J=1,N; if(J.NE.K) B(I,J)=B(I,J)-B(I,K)*B(K,J); end do
			end if
		end do
		do I=1,N; if(I.NE.K) B(I,K)=-B(I,K)*B(K,K);	end do
	end do
	do K=N,1,-1
		do J=1,N
			T=B(K,J); B(K,J)=B(JS(K),J); B(JS(K),J)=T
		end do
		do I=1,N
			T=B(I,K); B(I,K)=B(I,IS(K)); B(I,IS(K))=T
		end do
	end do
	return
end subroutine matinv


!!----------------------------------------------------------------------------------------

subroutine ksteel02(props,s,e,de,Et,statev,spd, yield, IDeath,M,ndof)
!
	implicit none
	real*8 E0,sy0,eta,mu,gama,esoft,alpha,beta,a_k,Omega
    real*8 emax,emin,ert,srt,erc,src,Ehc,Eh1,dt,dc,eu
    real*8 de,s,e,s0,et,e_unload,sign,sy,evs,eve,epeak,smax,max
    real*8 sres,eres,x,e_slip,s_slip,e_close,s_close,srel,ET1
    real*8 smin,spd,strain_end
    
    !real*8 mu
	integer kon, yield, IDeath,ndof
    real*8 M(NDOF,NDOF)!
    real*8 props(10), statev(11)

      E0  = props(1) 
      sy0 = props(2) 
      eta = props(3) 
      mu  = props(4) 
      gama= props(5) 
      esoft=props(6) 
      alpha= props(7) 
      beta = props(8) 
      a_k= props(9) 
      Omega= props(10) 
      emax  = statev(1) !maximum strain
      emin  = statev(2) !minimum strain
      ert   = statev(3) !stain at load reversal toward tension
      srt   = statev(4)
      erc   = statev(5) !stain at load reversal toward compression
      src   = statev(6)
      kon   = nint(statev(7)) !
      Ehc   = statev(8) !effective cummulative hysteresis energy
      Eh1   = statev(9) !hysteresis energy in a half cycle
      dt    = statev(10) !damage index for tension
      dc    = statev(11) !damage index for compression
      
      eu    = mu * sy0/E0 !characteristic ultimate strain


	if(Omega<=0.) Omega=0.5  
	if(a_k<0.) a_k=0. 
	if(eta<=0.) eta=1.d-6;      
	if(esoft>=0.) esoft=-1.d-6  
      
      if (kon.eq.0) then
        emax =  sy0/E0
        emin = -beta*sy0/E0
        if (de.ge.0.0) then
            kon = 1
        else
            kon = 2
        end if
      else if ((kon.eq.1).and.(de.lt.0.0)) then !Load reversal
            kon = 2
            if (s.gt.0.0) then
                erc = e
                src = s
            end if
            Ehc = Ehc + Eh1 * (erc / eu ) ** 2.0
            Eh1 = 0.0 !a new half cycle is to begin
            if (e.gt.emax) emax = e
      else if ((kon.eq.2).and.(de.gt.0.0)) then !Load reversal
            kon = 1
            if (s.lt.0.0) then
                ert = e
                srt = s
            endif
            Ehc = Ehc + Eh1 * (ert / eu ) ** 2.0
            Eh1 = 0.0 !a new half cycle is to begin
            if (e.lt.emin) emin = e
      end if
!
	s0=s
	s = s + E0 * de
	Et = E0
	
	if(a_k>0.) then
		if(s0>0.) E_unload=E0*(abs(emax/(sy0/E0)))**(-a_k)
		if(s0<0.) E_unload=E0*(abs(emin/(sy0/E0)))**(-a_k) 
		if(E_unload<0.1*E0) E_unload=0.1*E0
	else
		E_unload=E0
	end if

	if(s0*de<0.) then 
		s = s0 + E_unload * de
		Et = E_unload
		if(s*s0<0.) then 
			de=de-s0/E_unload
			s0=1.D-6*sy0*sign(1.d0,s) 
			Et=E0
		end if
	end if

		
      if ( de .ge. 0.0 .and. s0>=0.) then
          sy = (1.0 - dt) * sy0
          !loading envelope
		! Hardening
		if(e+de>sy/E0) then
			evs = max( sy + ( e + de - sy/E0) * eta * E0, 0.)
			evE = eta * E0
	       if (s .ge. evs) then
		      s = evs
		      Et = evE
				yield=1
		   end if
		end if
		! Softening
		epeak=sy/E0+(alpha-1.)*sy/E0/eta
		if(e+0.5*de>epeak) then
			evs=max(sy*alpha+esoft*E0*(e+de-epeak),0.0*sy)
			if(sy*alpha+esoft*E0*(e+de-epeak)<=0.) IDeath=1 ! Complete damage
			evE=esoft*E0
		    if (s .ge. evs) then
			   s = evs
			   Et = evE
			   yield=1
		  end if
		end if

          !reloading envelope
          smax = max(sy, sy + (emax - sy/E0) * eta * E0)  
		if(emax>epeak)	then
			smax=max(sy*alpha+esoft*E0*(emax-epeak),0.0*sy)
		end if	
          sres = 0.02 * smax                         
          eres = ert - (srt - sres) / E_unload                 

		x=emax-smax/E0  
		e_slip=gama*emax+(1.-gama)*x
		s_slip=smax*gama 
		e_close=e_slip*Omega 
		s_close=(e_close-eres)/(e_slip-eres) * (s_slip-sres) + sres

          if (eres .le. emax - smax / E0) then   
			if(e+0.5*de<e_close)  then  
				srel = (e+de-eres)/(e_slip-eres) * (s_slip-sres) + sres
				Et1=(s_slip-sres)/(e_slip-eres)
			else                 
				srel=(e+de-e_close)/(emax-e_close)*(smax-s_close)+s_close
				Et1 = (smax - s_close) / (emax - e_close)
			end if
            if (s .gt. srel) then
               s = max( srel, 0.)
               Et = Et1
            end if
          end if

      elseif ( de .lt. 0.0 .and. s0<0. ) then
          sy = (1.0 - dc) * sy0 *beta
          !loading envelope
		! Hardening
		if(e+de<-sy/E0) then
			evs =  min(-sy + ( e + de + sy/E0) * eta * E0,-0.0*sy)
			evE = eta * E0
			if (s .le. evs) then
				s = evs
				Et = evE
				yield=1
			end if
		end if
		! Softening
		epeak=-sy/E0-(alpha-1.)*sy/E0/eta
		if(e+0.5*de<epeak) then
			evs=min(-sy*alpha+esoft*E0*(e+de-epeak),-0.*sy)
			if(-sy*alpha+esoft*E0*(e+de-epeak)>=0.) IDeath=1 
			evE=esoft*E0
		    if (s .le. evs) then
			    s = evs
			    Et = evE
				yield=1
		  end if
		end if

          !reloading envelope
          smin = min(-sy, -sy + (emin + sy/E0) * eta * E0)
 		if(emin<epeak)	then 
			smin=min(-sy*alpha+esoft*E0*(emin-epeak),0.)
		end if	
          sres = 0.02 * smin 
          eres = erc - (src - sres) /  E_unload

		x=emin-smin/E0 
		e_slip=gama*emin+(1.-gama)*x 
		s_slip=smin*gama 
		e_close=e_slip*Omega 
		s_close=(e_close-eres)/(e_slip-eres) * (s_slip-sres) + sres

          if (eres .ge. emin - smin / E0) then    
			if(e+0.5*de>e_close) then 
				srel = (e+de-eres)/(e_slip-eres) * (s_slip-sres) + sres
				Et1=(s_slip-sres)/(e_slip-eres)
			else         
				srel=(e+de-e_close)/(emin-e_close)*(smin-s_close)+s_close
				Et1 = (smin - s_close) / (emin - e_close)
			end if
            if (s .lt. srel) then
                s = min (srel, 0.)
                Et = Et1
            end if
          end if
      end if


      if (Et.ne.E0 .and. Et.ne. E_unload) then 
            spd = spd + s * de
            Eh1 = Eh1 + s * de
            if ( s .ge. 0.0 ) then
                dc = min(Ehc /(3.0 * beta* sy0 * eu), 0.7)
            else
                dt = min(Ehc /(3.0 * sy0 * eu), 0.7)
            end if
      end if

      	x=max(sy0, beta*sy0)
	epeak=x/E0+(alpha-1.)*x/E0/eta
	Strain_End=epeak+abs(x/(esoft*E0))
	x=max(abs(emax),abs(emin),abs(e+de))

	if (IDeath==1) then
		s=0; 
		Et=1E-6*E0
	end if
!
	statev(1)   = emax
	statev(2)   = emin
	statev(3)   = ert
	statev(4)   = srt
	statev(5)   = erc
	statev(6)   = src
	statev(7)   = kon
	statev(8)   = Ehc
	statev(9)   = Eh1
	statev(10)  = dt
	statev(11)  = dc
    return
end subroutine ksteel02
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	