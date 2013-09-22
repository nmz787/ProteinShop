C/***************************************************************************************
C*cr	1. FUNCTIONS OF THE PROGRAM
C*cr	---------------------------
C*cr	ANALYZE processes the dihehral-angle outo.* files obtained from 
C*cr	calculations (usually global conformational analysis with the EDMC
C*cr	method) using ECEPPAK. These functions include the following:
C*cr	
C*cr	1. Calculations of conformational characteristics, such as hydrogen bonds,
C*cr	   turn position and types, RMS deviation from a reference conformation,
C*cr	   interchromophore distances, interproton distances, etc.
C*cr	
C*cr	2. Calculation of Boltzmann-averaged properties of the conformational 
C*cr	   ensemble.
C*cr	
C*cr	3. Calculate the dihedral angles from supplied Cartesian coordinates.
C*cr	
C*cr	4. Cluster analysis of the conformational ensemble by the minimal spanning
C*cr	   tree or minimum-variance method.
C*cr	
C*cr	5. Fitting the statistical weights of the conformations so as to achieve
C*cr	   the best agreement between the calculated average and experimental NOE
C*cr	   spectra and coupling constants. 
C*cr	
C*cr	Details of the functions can be found in INPUT.DOC
C*cr	
C*cr	2. CONTENTS 
C*cr	-----------
C*cr	The package includes two versions of ANALYZE: one for fitting NMR spectra
C*cr	and one for cluster analysis and other purposes. The source files are 
C*cr	contained in source_clust and source_nmr, respectively. This division is 
C*cr	caused by practical reasons: the NMR and clustering parts are very memory 
C*cr	consuming, which practically excludes their incorporation into one program. 
C*cr	
C*cr	This NMR part of ANALYZE includes the core part of the MORASS package from 
C*cr	Dr. D. Gorenstein laboratory, Sealy Center for Structural Biology, 
C*cr	University of Texas Medical Branch, Galveston, TX (http://www.nmr.utmb.edu/).
C*cr	
C*cr	The authors of ANALYZE are grateful to the MORASS developer group for the 
C*cr	permission of including the MORASS modules in the publicly available version 
C*cr	of ANALYZE. Note that this permission was granted subject to the condition 
C*cr	that the software remains free for academic and educational use,
C*cr	as is the MORASS package.
C*cr	
C*cr	For more information about MORASS see http://www.nmr.utmb.edu/#mrass.
C*cr	
C*cr	For the program to work, some database files are required. The location
C*cr	of the database directory must be specified by the DBASDIR environmental
C*cr	variable. See INSTALL.DOC and COMMAND.DOC for details.
C*cr	
C*cr	3. AVAILABILITY AND CITING
C*cr	--------------------------
C*cr	The package is available for free to academic users. When publishing the 
C*cr	work based on this program please cite the following references:
C*cr	
C*cr	ECEPP/ECEPPAK:
C*cr	
C*cr	1. F.A. Momany, R.F. McGuire, A.W. Burgess and H.A. Scheraga,
C*cr	   J. Phys. Chem., 79, 2361-2381 (1975).
C*cr	
C*cr	2. G. Nemethy and H.A. Scheraga, J. Phys. Chem. 87, 1883-1891 (1983).
C*cr	
C*cr	3. M.J. Sippl, G. Nemethy, and H.A. Scheraga, J. Phys. Chem. 88, 6231-6233 
C*cr	   (1984).
C*cr	
C*cr	4. G. Nemethy, K.D. Gibson, K.A. Palmer, C.N. Yoon, C.N., G. Paterlini, 
C*cr	   A. Zagari, S. Rumsey, S., and H.A. Scheraga, H.A. J. Phys. Chem., 96, 
C*cr	   6472-6484 (1992).
C*cr	
C*cr	5. D.R. Ripoll, M.S. Pottle, K.D. Gibson, A. Liwo and H.A. Scheraga.
C*cr	   J. Comput. Chem., 16, 1153-1163 (1995).
C*cr	
C*cr	6. D.R. Ripoll, A. Liwo and C. Czaplewski. TASK Quart., 3, 313-323 (1999).
C*cr	
C*cr	NMR fitting:
C*cr	
C*cr	1. J. Malicka, M. Groth, C. Czaplewski, A. Liwo, W. Wiczk, and L. Lankiewicz.
C*cr	   Lett. Pept. Sci., 5, 445-447 (1998).
C*cr	
C*cr	2. M. Groth, J. Malicka, C.  Czaplewski, S. Oldziej, L. Lankiewicz, 
C*cr	   W. Wiczk and A. Liwo, J. Biomol. NMR, 4, 315-330 (1999).
C*cr	
C*cr	3. D.R. Ripoll, A. Liwo and C. Czaplewski. TASK Quart., 3, 313-323 (1999).
C*cr	
C*cr	MORASS:
C*cr	
C*cr	1. Robert P. Meadows, Carol Beth Post, Bruce A. Luxon, and David G. Gorenstein,
C*cr	   MORASS 2.1, Purdue University, W. Lafayette (1994).
C*cr	
C*cr	2. C.B. Post, R.P. Meadows and D.G. Gorenstein, J. Am. Chem. Soc., 112, 6796 
C*cr	   (1990)
C*cr	
C*cr	4. SUPPORT 
C*cr	----------
C*cr	User's feedback is greatly appreciated. Please send questions and comments
C*cr	to:
C*cr	
C*cr	Dr. Adam Liwo
C*cr	Faculty of Chemistry, University of Gdansk
C*cr	Sobieskiego 18, 80-952 Gdansk, Poland
C*cr	
C*cr	e-mail: adam@rutyl.chem.univ.gda.pl
C*cr	
C*cr	or
C*cr	
C*cr	Dr. Cezary Czaplewski
C*cr	Baker Laboratory of Chemistry, Cornell University
C*cr	Ithaca, NY 14853-1301
C*cr	
C*cr	e-mail: czarek@scheraga2.chem.cornell.edu
C*cr	
C*cr
C***************************************************************************************/

C		DOUBLE FUNCTION FITSQ(RMS,X,Y,NN,T,B)
      SUBROUTINE FITSQ(RMS,X,Y,NN,T,B)
C  X AND Y ARE THE VECTORS OF COORDINATES (DIMENSIONED (3,N)) OF THE TWOSUP00860
C  STRUCTURES TO BE SUPERIMPOSED.  NN IS 3*N, WHERE N IS THE NUMBER OF  SUP00870
C  POINTS.   T AND B ARE RESPECTIVELY THE TRANSLATION VECTOR AND THE    SUP00880
C  ROTATION MATRIX THAT TRANSFORMS THE SECOND SET OF COORDINATES TO THE SUP00890
C  FRAME OF THE FIRST SET.                                              SUP00900
C  ETA =  MACHINE-SPECIFIC VARIABLE                                     SUP00910
                                                                        SUP00920
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(*),Y(*),T(3)                                          SUP00930
      DIMENSION B(3,3),Q(3,3),R(3,3),V(3),XAV(3),YAV(3),E(3),C(3,3)     SUP00940
C     SMALL=25.0*RMDCON(3)                                              SUP00960
C     SMALL=25.0*ETA                                                    SUP00970
C     SMALL=25.0*10.E-10                                                SUP00980
C THE FOLLOWING IS A VERY LENIENT VALUE FOR 'SMALL'                     SUP00990
      SMALL = 0.0001                                                    SUP01000
      FN=NN                                                             SUP01010
      DO 10 I=1,3                                                       SUP01020
      XAV(I)=0.0                                                        SUP01030
      YAV(I)=0.0                                                        SUP01040
      DO 10 J=1,3                                                       SUP01050
   10 B(J,I)=0.0                                                        SUP01060
      NC=0                                                              SUP01070
C                                                                       SUP01080
      DO 30 N=1,NN                                                      SUP01090
      DO 20 I=1,3                                                       SUP01100
C       PRINT*,'X = ',X(NC+I),'  Y = ',Y(NC+I)                           SUP01110
C		DIFF =  X(NC+I)-Y(NC+I)       
C		PRINT*,'Diff = ', (X(NC+I)-Y(NC+I))                			     
      XAV(I)=XAV(I)+X(NC+I)/FN                                          SUP01120
   20 YAV(I)=YAV(I)+Y(NC+I)/FN                                          SUP01130
   30 NC=NC+3                                                           SUP01140
C                                                                       SUP01150
                                                                        SUP01160
c      PRINT*,'XAV = ',(XAV(J),J=1,3)                                    SUP01170
c      PRINT*,'YAV = ',(YAV(J),J=1,3)                                    SUP01180
                                                                        SUP01190
      NC=0                                                              SUP01200
      RMS=0.0                                                           SUP01210
      DO 50 N=1,NN                                                      SUP01220
      DO 40 I=1,3                                                       SUP01230
      RMS=RMS+((X(NC+I)-XAV(I))**2+(Y(NC+I)-YAV(I))**2)/FN              SUP01240
      DO 40 J=1,3                                                       SUP01250
      B(J,I)=B(J,I)+(X(NC+I)-XAV(I))*(Y(NC+J)-YAV(J))/FN                SUP01260
   40 C(J,I)=B(J,I)                                                     SUP01270
   50 NC=NC+3                                                           SUP01280
      CALL SIVADE(B,Q,R,D)                                              SUP01290
      SN3=DSIGN(1.0D0,D)                                                SUP01300
      DO 120 I=1,3                                                      SUP01310
      DO 120 J=1,3                                                      SUP01320
  120 B(J,I)=-Q(J,1)*R(I,1)-Q(J,2)*R(I,2)-SN3*Q(J,3)*R(I,3)             SUP01330
      CALL MVVAD(B,XAV,YAV,T)                                           SUP01340
      DO 130 I=1,3                                                      SUP01350
      DO 130 J=1,3                                                      SUP01360
      RMS=RMS+2.0*C(J,I)*B(J,I)                                         SUP01370
  130 B(J,I)=-B(J,I)                                                    SUP01380
      IF (DABS(RMS).GT.SMALL) GO TO 140                                  SUP01390
*     	WRITE (6,301)                                                     SUP01400
      RETURN                                                            SUP01410
  140 IF (RMS.GT.0.0) GO TO 150                                         SUP01420
*     	WRITE (6,303) RMS                                                 SUP01430
      rms=0.0
*     STOP                                                              SUP01440
*  150 WRITE (6,302) DSQRT(RMS)                                           SUP01450
  150 continue
      RETURN                                                            SUP01460
  301 FORMAT (5X,'RMS DEVIATION NEGLIGIBLE')                            SUP01470
  302 FORMAT (5X,'RMS DEVIATION ',F14.6)                                SUP01480
  303 FORMAT (//,5X,'NEGATIVE MS DEVIATION - ',F14.6)                   SUP01490
      END                                                               SUP01500

      SUBROUTINE SIVADE(X,Q,R,DT)                                       SUP01510
C  COMPUTES Q,E AND R SUCH THAT Q(T)XR = DIAG(E)                        SUP01520
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(3,3),Q(3,3),R(3,3),E(3)                               SUP01530
      DIMENSION H(3,3),P(3,3),U(3,3),D(3)                               SUP01540
      SMALL=25.0*10.E-10                                                SUP01560
C     SMALL=25.0*ETA                                                    SUP01570
C     SMALL=2.0*RMDCON(3)                                               SUP01580
      XNRM=0.0                                                          SUP01590
      DO 20 I=1,3                                                       SUP01600
      DO 10 J=1,3                                                       SUP01610
      XNRM=XNRM+X(J,I)*X(J,I)                                           SUP01620
      U(J,I)=0.0                                                        SUP01630
      R(J,I)=0.0                                                        SUP01640
   10 H(J,I)=0.0                                                        SUP01650
      U(I,I)=1.0                                                        SUP01660
   20 R(I,I)=1.0                                                        SUP01670
      XNRM=DSQRT(XNRM)                                                   SUP01680
      DO 110 N=1,2                                                      SUP01690
      XMAX=0.0                                                          SUP01700
      DO 30 J=N,3                                                       SUP01710
   30 IF (DABS(X(J,N)).GT.XMAX) XMAX=DABS(X(J,N))                         SUP01720
      A=0.0                                                             SUP01730
      DO 40 J=N,3                                                       SUP01740
      H(J,N)=X(J,N)/XMAX                                                SUP01750
   40 A=A+H(J,N)*H(J,N)                                                 SUP01760
      A=DSQRT(A)                                                         SUP01770
      DEN=A*(A+DABS(H(N,N)))                                             SUP01780
      D(N)=1.0/DEN                                                      SUP01790
      H(N,N)=H(N,N)+DSIGN(A,H(N,N))                                      SUP01800
      DO 70 I=N,3                                                       SUP01810
      S=0.0                                                             SUP01820
      DO 50 J=N,3                                                       SUP01830
   50 S=S+H(J,N)*X(J,I)                                                 SUP01840
      S=D(N)*S                                                          SUP01850
      DO 60 J=N,3                                                       SUP01860
   60 X(J,I)=X(J,I)-S*H(J,N)                                            SUP01870
   70 CONTINUE                                                          SUP01880
      IF (N.GT.1) GO TO 110                                             SUP01890
      XMAX=DMAX1(DABS(X(1,2)),DABS(X(1,3)))                               SUP01900
      H(2,3)=X(1,2)/XMAX                                                SUP01910
      H(3,3)=X(1,3)/XMAX                                                SUP01920
      A=DSQRT(H(2,3)*H(2,3)+H(3,3)*H(3,3))                               SUP01930
      DEN=A*(A+DABS(H(2,3)))                                             SUP01940
      D(3)=1.0/DEN                                                      SUP01950
      H(2,3)=H(2,3)+DSIGN(A,H(2,3))                                      SUP01960
      DO 100 I=1,3                                                      SUP01970
      S=0.0                                                             SUP01980
      DO 80 J=2,3                                                       SUP01990
   80 S=S+H(J,3)*X(I,J)                                                 SUP02000
      S=D(3)*S                                                          SUP02010
      DO 90 J=2,3                                                       SUP02020
   90 X(I,J)=X(I,J)-S*H(J,3)                                            SUP02030
  100 CONTINUE                                                          SUP02040
  110 CONTINUE                                                          SUP02050
      DO 130 I=1,3                                                      SUP02060
      DO 120 J=1,3                                                      SUP02070
  120 P(J,I)=-D(1)*H(J,1)*H(I,1)                                        SUP02080
  130 P(I,I)=1.0+P(I,I)                                                 SUP02090
      DO 140 I=2,3                                                      SUP02100
      DO 140 J=2,3                                                      SUP02110
      U(J,I)=U(J,I)-D(2)*H(J,2)*H(I,2)                                  SUP02120
  140 R(J,I)=R(J,I)-D(3)*H(J,3)*H(I,3)                                  SUP02130
      CALL MMMUL(P,U,Q)                                                 SUP02140
  150 NP=1                                                              SUP02150
      NQ=1                                                              SUP02160
      IF (DABS(X(2,3)).GT.SMALL*(DABS(X(2,2))+DABS(X(3,3)))) GO TO 160     SUP02170
      X(2,3)=0.0                                                        SUP02180
      NQ=NQ+1                                                           SUP02190
  160 IF (DABS(X(1,2)).GT.SMALL*(DABS(X(1,1))+DABS(X(2,2)))) GO TO 180     SUP02200
      X(1,2)=0.0                                                        SUP02210
      IF (X(2,3).NE.0.0) GO TO 170                                      SUP02220
      NQ=NQ+1                                                           SUP02230
      GO TO 180                                                         SUP02240
  170 NP=NP+1                                                           SUP02250
  180 IF (NQ.EQ.3) GO TO 310                                            SUP02260
      NPQ=4-NP-NQ                                                       SUP02270
      IF (NP.GT.NPQ) GO TO 230                                          SUP02280
      N0=0                                                              SUP02290
      DO 220 N=NP,NPQ                                                   SUP02300
      NN=N+NP-1                                                         SUP02310
      IF (DABS(X(NN,NN)).GT.SMALL*XNRM) GO TO 220                        SUP02320
      X(NN,NN)=0.0                                                      SUP02330
      IF (X(NN,NN+1).EQ.0.0) GO TO 220                                  SUP02340
      N0=N0+1                                                           SUP02350
      GO TO (190,210,220),NN                                            SUP02360
  190 DO 200 J=2,3                                                      SUP02370
  200 CALL GIVNS(X,Q,1,J)                                               SUP02380
      GO TO 220                                                         SUP02390
  210 CALL GIVNS(X,Q,2,3)                                               SUP02400
  220 CONTINUE                                                          SUP02410
      IF (N0.NE.0) GO TO 150                                            SUP02420
  230 NN=3-NQ                                                           SUP02430
      A=X(NN,NN)*X(NN,NN)                                               SUP02440
      IF (NN.GT.1) A=A+X(NN-1,NN)*X(NN-1,NN)                            SUP02450
      B=X(NN+1,NN+1)*X(NN+1,NN+1)+X(NN,NN+1)*X(NN,NN+1)                 SUP02460
      C=X(NN,NN)*X(NN,NN+1)                                             SUP02470
      DD=0.5*(A-B)                                                      SUP02480
      XN2=C*C                                                           SUP02490
      RT=B-XN2/(DD+DSIGN(DSQRT(DD*DD+XN2),DD))                            SUP02500
      Y=X(NP,NP)*X(NP,NP)-RT                                            SUP02510
      Z=X(NP,NP)*X(NP,NP+1)                                             SUP02520
      DO 300 N=NP,NN                                                    SUP02530
      IF (DABS(Y).LT.DABS(Z)) GO TO 240                                   SUP02540
      T=Z/Y                                                             SUP02550
      C=1.0/DSQRT(1.0+T*T)                                               SUP02560
      S=C*T                                                             SUP02570
      GO TO 250                                                         SUP02580
  240 T=Y/Z                                                             SUP02590
      S=1.0/DSQRT(1.0+T*T)                                               SUP02600
      C=S*T                                                             SUP02610
  250 DO 260 J=1,3                                                      SUP02620
      V=X(J,N)                                                          SUP02630
      W=X(J,N+1)                                                        SUP02640
      X(J,N)=C*V+S*W                                                    SUP02650
      X(J,N+1)=-S*V+C*W                                                 SUP02660
      A=R(J,N)                                                          SUP02670
      B=R(J,N+1)                                                        SUP02680
      R(J,N)=C*A+S*B                                                    SUP02690
  260 R(J,N+1)=-S*A+C*B                                                 SUP02700
      Y=X(N,N)                                                          SUP02710
      Z=X(N+1,N)                                                        SUP02720
      IF (DABS(Y).LT.DABS(Z)) GO TO 270                                   SUP02730
      T=Z/Y                                                             SUP02740
      C=1.0/DSQRT(1.0+T*T)                                               SUP02750
      S=C*T                                                             SUP02760
      GO TO 280                                                         SUP02770
  270 T=Y/Z                                                             SUP02780
      S=1.0/DSQRT(1.0+T*T)                                               SUP02790
      C=S*T                                                             SUP02800
  280 DO 290 J=1,3                                                      SUP02810
      V=X(N,J)                                                          SUP02820
      W=X(N+1,J)                                                        SUP02830
      A=Q(J,N)                                                          SUP02840
      B=Q(J,N+1)                                                        SUP02850
      X(N,J)=C*V+S*W                                                    SUP02860
      X(N+1,J)=-S*V+C*W                                                 SUP02870
      Q(J,N)=C*A+S*B                                                    SUP02880
  290 Q(J,N+1)=-S*A+C*B                                                 SUP02890
      IF (N.GE.NN) GO TO 300                                            SUP02900
      Y=X(N,N+1)                                                        SUP02910
      Z=X(N,N+2)                                                        SUP02920
  300 CONTINUE                                                          SUP02930
      GO TO 150                                                         SUP02940
  310 DO 320 I=1,3                                                      SUP02950
  320 E(I)=X(I,I)                                                       SUP02960
  330 N0=0                                                              SUP02970
      DO 360 I=1,3                                                      SUP02980
      IF (E(I).GE.0.0) GO TO 350                                        SUP02990
      E(I)=-E(I)                                                        SUP03000
      DO 340 J=1,3                                                      SUP03010
  340 Q(J,I)=-Q(J,I)                                                    SUP03020
  350 IF (I.EQ.1) GO TO 360                                             SUP03030
      IF (DABS(E(I)).LT.DABS(E(I-1))) GO TO 360                           SUP03040
      CALL SWITCH(I,1,Q,R,E)                                            SUP03050
      N0=N0+1                                                           SUP03060
  360 CONTINUE                                                          SUP03070
      IF (N0.NE.0) GO TO 330                                            SUP03080
      IF (DABS(E(3)).GT.SMALL*XNRM) GO TO 370                            SUP03090
      E(3)=0.0                                                          SUP03100
      IF (DABS(E(2)).GT.SMALL*XNRM) GO TO 370                            SUP03110
      E(2)=0.0                                                          SUP03120
  370 DT=DET(Q(1,1),Q(1,2),Q(1,3))*DET(R(1,1),R(1,2),R(1,3))            SUP03130
*     WRITE (1,501) (E(I),I=1,3)                                        SUP03140
      RETURN                                                            SUP03150
  501 FORMAT (/,5X,'SINGULAR VALUES - ',3E15.5)                         SUP03160
      END                                                               SUP03170

      SUBROUTINE GIVNS(A,B,M,N)                                         SUP03180
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(3,3),B(3,3)                                           SUP03190
      IF (DABS(A(M,N)).LT.DABS(A(N,N))) GO TO 10                          SUP03200
      T=A(N,N)/A(M,N)                                                   SUP03210
      S=1.0/DSQRT(1.0+T*T)                                               SUP03220
      C=S*T                                                             SUP03230
      GO TO 20                                                          SUP03240
   10 T=A(M,N)/A(N,N)                                                   SUP03250
      C=1.0/DSQRT(1.0+T*T)                                               SUP03260
      S=C*T                                                             SUP03270
   20 DO 30 J=1,3                                                       SUP03280
      V=A(M,J)                                                          SUP03290
      W=A(N,J)                                                          SUP03300
      X=B(J,M)                                                          SUP03310
      Y=B(J,N)                                                          SUP03320
      A(M,J)=C*V-S*W                                                    SUP03330
      A(N,J)=S*V+C*W                                                    SUP03340
      B(J,M)=C*X-S*Y                                                    SUP03350
   30 B(J,N)=S*X+C*Y                                                    SUP03360
      RETURN                                                            SUP03370
      END                                                               SUP03380

      SUBROUTINE SWITCH(N,M,U,V,D)                                      SUP03390
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION U(3,3),V(3,3),D(3)                                      SUP03400
      DO 10 I=1,3                                                       SUP03410
      TEM=U(I,N)                                                        SUP03420
      U(I,N)=U(I,N-1)                                                   SUP03430
      U(I,N-1)=TEM                                                      SUP03440
      IF (M.EQ.0) GO TO 10                                              SUP03450
      TEM=V(I,N)                                                        SUP03460
      V(I,N)=V(I,N-1)                                                   SUP03470
      V(I,N-1)=TEM                                                      SUP03480
   10 CONTINUE                                                          SUP03490
      TEM=D(N)                                                          SUP03500
      D(N)=D(N-1)                                                       SUP03510
      D(N-1)=TEM                                                        SUP03520
      RETURN                                                            SUP03530
      END                                                               SUP03540
C     SUBROUTINE MVVAD(A,B,C,D)                                         SUP03550

      SUBROUTINE MVVAD(B,XAV,YAV,T)                                     SUP03560
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION B(3,3),XAV(3),YAV(3),T(3)                               SUP03570
C     DIMENSION A(3,3),B(3),C(3),D(3)                                   SUP03580
C     DO 10 J=1,3                                                       SUP03590
C     D(J)=C(J)                                                         SUP03600
C     DO 10 I=1,3                                                       SUP03610
C  10 D(J)=D(J)+A(J,I)*B(I)                                             SUP03620
      DO 10 J=1,3                                                       SUP03630
      T(J)=YAV(J)                                                       SUP03640
      DO 10 I=1,3                                                       SUP03650
   10 T(J)=T(J)+B(J,I)*XAV(I)                                           SUP03660
      RETURN                                                            SUP03670
      END                                                               SUP03680

      DOUBLE PRECISION FUNCTION DET (A,B,C)                             SUP03690
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(3),B(3),C(3)                                          SUP03700
      DET=A(1)*(B(2)*C(3)-B(3)*C(2))+A(2)*(B(3)*C(1)-B(1)*C(3))         SUP03710
     1  +A(3)*(B(1)*C(2)-B(2)*C(1))                                     SUP03720
      RETURN                                                            SUP03730
      END                                                               SUP03740

      SUBROUTINE MMMUL(A,B,C)                                           SUP03750
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(3,3),B(3,3),C(3,3)                                    SUP03760
      DO 10 I=1,3                                                       SUP03770
      DO 10 J=1,3                                                       SUP03780
      C(I,J)=0.0                                                        SUP03790
      DO 10 K=1,3                                                       SUP03800
   10 C(I,J)=C(I,J)+A(I,K)*B(K,J)                                       SUP03810
      RETURN                                                            SUP03820
      END                                                               SUP03830

      SUBROUTINE MATVEC(UVEC,TMAT,PVEC,NBACK)                           SUP03840
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION TMAT(3,3),UVEC(3,NBACK), PVEC(3,NBACK)           SUP03850
      DO 2 J=1,NBACK                                                    SUP03870
         DO 1 I=1,3                                                     SUP03880
         UVEC(I,J) = 0.0                                                SUP03890
         DO 1 K=1,3                                                     SUP03900
    1    UVEC(I,J)=UVEC(I,J)+TMAT(I,K)*PVEC(K,J)                        SUP03910
    2 CONTINUE                                                          SUP03920
      RETURN                                                            SUP03930
      END                                                               SUP03940

