    
      subroutine spanr(x,y,freq,n,m, long,combA, combAc,hc,hg,copt,
     &          betain,status,npt,qo,po,gopt,gfreq,binw,gamma, surv)
c!dir$ ATTRIBUTES DLLEXPORT :: SPANR
c!dir$ ATTRIBUTES C, REFERENCE, ALIAS : 'spanr_' :: SPANR
c-------------------------------------------------------------
C                      S P A N R
c---------------------------------------------------------------
c   subroutine SPANR for conversion to R
c   Input   X  n*m matrix of binary (assumed) positive predictors
c              of if surv=T is an n*(m+1) matrix, the m+1th being event
c             indicator 
c           n  rows of X
c           m  (or m+1) columns of X  
c           freq  n * 1  frequency count of  each row of X
c           y   outcome n * 1 column 
c           long 
c  output  combA, combAc are returned integer arrays of resultant Boolean
c          combinations. Each of the form
c          0  2 5 7  0  5 3  0  7 0...
c         where 0 represents a bracket or pair of brackets
c           ( 2 5 7)( 5 3) (7)....          
c-------------------------------------------------------------

C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS
      PARAMETER (MAXP=100, MAXADD=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================


 
      integer combAc(maxp2), combA(maxp2)
c
      integer LPX(MAXADD,MAXP2),
     *LPXC(MAXADD,MAXP2),LP(MAXP2),LPP(MAXP2),LPC(MAXP2),
     *LPNEW(MAXP2),LPNEWC(MAXP2)
c   complexities up to MAXP      
      real  Gcmplx(MAXCX), hg(MAXCX), hc(MAXCX)
      
      integer Qcmpx(MAXCX), Pcmplx(MAXCX,maxp),
     *  LPcx(MAXCX,MAXP2), IDcx(MAXCX)
      
      integer Qccx(MAXCX), Pccx(MAXCX,maxp),
     *  LPccx(MAXCX,MAXP2)
      
      integer c, cx(MAXCX), nbn(MAXCX), cmax, copt
      
      integer qx(maxadd), qxc(maxadd), px(maxadd,maxp), pxc(maxadd,maxp)

      integer qopt, popt(maxp),lpopt(maxp2)
      integer qcopt, pcopt(maxp),lpcopt(maxp2)
      integer qoptv, poptv(maxp),lpoptv(maxp2)
      
      integer gfreq(2000)
      
           
      REAL y(n), ytemp(n)
      integer event(n)
      integer mrank(n) 
      integer x(n,m+1)

      integer f(long)
      integer freq(n)

c  need to KAEON this m+2 for GNU compilation
      integer naddx(m+2) 
c      integer naddx(100) 
      integer status, npt
      integer p(maxp),po(maxp), pp(maxp), u(maxp), pc(maxp)
      integer ip(maxp2),ipc(maxp2), ipx(maxp2)
      integer pnew(maxp), qnew, vnew
      integer pnewc(maxp), qnewc, vnewc
      integer lpwork(maxp2)
      integer q, qq, qo, qc
      logical fa, tr
      integer prod,v, vv, vc

      logical fini, begn,noemb, more
      logical idpart
       real ghold(m)
      
      double precision sy0(maxp), ssy0(maxp), sy1(maxp), ssy1(maxp)
      double precision SSAG0,SSA0,SSAG1,SSA1, TSS
      
      logical  subg, skip
      logical surv
            common/partix/ Dwork(10),Cwork(10),CUMNS(11),KEY(10)
            integer Dwork,Cwork,CUMNS,KEY
c       write(*,*)'Entering SPANR'
      status=0
c     status is return indicator:
c     1  =  no convergence  
c     2  =  overpenalised, input beta too large   
c     3  =  Boolean transposition error
c     4  =  Boolean expansion error  
      do 13 i=1,m+2
       naddx(i)=0
13    continue
      do 14 i=1,maxcx
      Gcmplx(i)=0
14    continue       
      skip=.true.
c  set up required arrays
      lc=maxadd+ m
     
c      write(*,*)'long =',long
c      long=lc*n
C      allocate(f(long))
c   make a long array of X data as required for calculation
c   of objective function G   
      
c       write(*,*)'long,maxadd,m,lc,n=',long,maxadd,m,lc,n
c       write(*,*)'lc*n=',lc*n
      do 555 i=1,n
          ipos=(i-1)*lc
          do 556 j=1,m
          l=ipos+j
          f(l)=x(i,j)
 556      continue
 555      continue
      fa=.false.
      tr=.true.
      nhold=0
 
      nptcheck=0
      noemb=tr
      subg=fa
      cmax=20
      maxc=0
c  testing gamma
c      gamma=0
c----------------------------------------      
c survival settings, last column of x is assumed event indicator
      if(surv) then
      do 321 i=1,n
      event(i)=x(i,m+1)
321   continue           
      ncall=0
      
c      write(*,*)"at surv=T?"
      end if 
c--------------------------------------      
c  initialise iteration no.
      its=0
      gopt=0
c  betain <1  indicates beta to be automatically set
c  otherwise  actual input value is used      
      if(betain.lt.0) then
      beta=0.0
      else
      beta=betain
      end if
c  try building g distribution - assume betain is input known      
       do 444 i=1,2000
       gfreq(i)=0
 444  continue
c     binw=beta*2
c      write(*,*)'binw=',binw
      npt=0
      k0=m
      gsum=0

      csum=0      
      csum2=0 
      gcsum=0
      nsum=0
c&&&      write(*,*) 'Size q=', qo, ' p(1)..p(q)=', (po(j),j=1,qo)
C-----------------------------------------------------------------------
C   Begin the search, or new iteration in iterative mode 
C--------------------------------------------------------------------- 

1557  continue   
c     starts a new iteration
        u(1)=1
          q=1
          p(1)=1
          gmax=0
        
           
8151      continue
c&&         WRITE(*,*) '**ITERATION ** its +1',ITS + 1
c&&         call intpr ("iteration", -1,its,1)
         

C********************************************************************
C  Key label 155 marks start of loop for increments of q, p1,...,pq *
C********************************************************************
155   continue
              id0=0
              if(its.eq.0.and.Q.eq.1.and.p(1).eq.1)id0=1
              id1=1
c              if(q.eq.1.and.p(1).eq.1) id0=1


      V=0
      DO 401 J=1,Q
401   V=V+ P(J)
c
C    ****************************************************
C    *  FROM HERE TO LABEL 200 PARTITIONS ARE GENERATED  *
C    *****************************************************
C  start partition for current Q, p1,p2,...,pq.
c
c&&      WRITE(*,8156) (P(J),J=1,Q)
c&&8156  FORMAT(10x,'Class p1...pq=',10i3)
C

c

c---------------------------------------------------------------
      
      BEGN=TR
      FINI=FA 
C  Skip automatic evalation of partitions...
151   CONTINUE
c

C----------------------------------------------------------------------
C  here lies the heart of the program. PARTI does the regular
C  partition generation.

c        write(*,*)'enter parti'
        CALL PARTI(K0,Q,P,IP,BEGN,FINI,NOEMB)
 
        IF(FINI) then
 
      GO TO 201
      end if 

C-----------------------------------------------------------------------
C  Compute the partition. First set up the vector LP which indexes
C the columns of the binary data array F to be accessed for the
C partition.
C
C  if in iterative mode skip unnecessary partitions which dont 
C  include current added variable i.e the K0th element
C
             IF(ITS.GE.1) THEN
             DO  557 J=1,V
             IF(IP(J).EQ.K0) GO TO 6760
557          continue
C
C exit loop , K0th element not included
C
             go to 151  
c             skip evaluation, go to next PARTI
             END IF 
 
6760         continue
c   note  that IP is combination of numbers range 1,m
C   but LP vector in SPAN must preclude lp()=1 which designates
c   null set. need to set  lp=IP+1

      do 558 j=1,v
       lp(j)=ip(j)+1
558    continue
              
c       if(its.eq.1.and.numb.lt.100) then
c      if(q.eq.3.and.p(1).eq.2.and.p(2).eq.2.and.p(3).eq.1) then
c          write(*,*)(lp(kunt),kunt=1,v)
c      end if
c      end if
           
c
            do 1166 id=id0,id1 
                 
             call trans(q,p,lp,v,qc,pc,lpc,vc,ifail)

c  need to do same to LPC to convert back to IPC
             
            if(ifail.ne.0) then
             status=3
             end if
             
             do 559 j=1,vc
             ipc(j)=lpc(j)-1
559          continue
c-------------------------------------------------------
           if(its.ge.1) then
         
c         call outit('lp',q,p,lp)
c   iteration >0, must be  an added attribute         
         if(id.eq.1) then 
       call expand(q,p,lp,qx,px,lpx,qxc,pxc,lpxc,naddx,qnew,
     *            pnew,lpnew,vnew,ifail)
            if(ifail.ne.0) then
             status=4
c             go to 151
             end if
         call expand(qc,pc,lpc,qxc,pxc,lpxc,qx,px,lpx,naddx,qnewc,
     *            pnewc,lpnewc,vnewc,ifail)
            if(ifail.ne.0) then
             status=4
c             go to 151
             end if
          else
             
       call expand(q,p,lp,qxc,pxc,lpxc,qx,px,lpx,naddx,qnew,
     *             pnew,lpnew,vnew,ifail)
            if(ifail.ne.0) then
             status=4
c             go to 151
             end if
        call expand(qc,pc,lpc,qx,px,lpx,qxc,pxc,lpxc,naddx,qnewc,
     *             pnewc,lpnewc,vnewc,ifail)
            if(ifail.ne.0) then
             status=4
c             go to 151
             end if
          end if
c         call outit('lpnew',qnew,pnew,lpnew)
          
          
          else  
c       else of if(its.ge.1
c  copy original index vector to new ones for 0th iteration
        CALL EXCHNG(Q, P ,LP ,QNEW ,VNEW, PNEW ,LPNEW)
        CALL EXCHNG(QC,PC,LPC,QNEWC,VNEWC,PNEWC,LPNEWC)
               
       
      end if 

      c= qnew+qnewc-1
c   cmax is maximum allowable complexity. =20 
c  1166 is do id=id0,id1 loop. Effectively ignores partition
c      write(*,*)'c,cmax',c,cmax
      if(c.gt.cmax) go to 1166  
      
c  need to account for LP vectors starting at 2
      do 674 j=1,vnew
          ipx(j)=lpnew(j)-1
674   continue
      
      do 675 j=1,vnewc
          ipc(j)=lpnewc(j)-1
675   continue
c--- COMPUTE  CRITERION STATISTIC G -------------------------      
      if(surv) then

       call logrnk(f,y,event,n,IPx,Qnew,Pnew,ID
     *,freq,skip,na0,na1,G,ncall,mrank,ytemp)
       
cc       if(G.gt.gopt-1) write(*,*) G
       ncall=ncall+1
  
       else

c  establish MSE criterion for a split, returned value G           
           CALL CRNCH(f,Y,N,IPx,Qnew,Pnew,IPC,QnewC,PnewC,
     *     ID,FREQ, skip,SSAG0,SSA0,SSAG1,SSA1,
     *     SUBG, sy0,sy1, ssy0,ssy1, na0,na1, tss, 
     *     nag0,nag1, G)
      end if     
c   balance if gamma>0
               if(gamma.gt.0) then
               pa=float(na1)/float(na0+na1)
               G = G * (pa*(1-pa))**gamma
               end if 
           
c------------------------------------------------------------------           
c   increment best at this complexity      

      if(c.ge.1) then 
        if(G.gt.Gcmplx(c)) then
c        if(g.ge.Gcmplx(c)) then
          Gcmplx(c)=G
         Qcmpx(c)=qnew
         Qccx(c)=qnewc
         IDcx(c)=id
          cx(c)=c
         l=0
         do 554 i=1, qnew
         Pcmplx(c,i)=pnew(i)
            do 553 j=1,pnew(i)
                l=l+1
                LPcx(c,l)=lpnew(l)
553          continue
554       continue
         
                 l=0
         do 552 i=1, qnewc
         Pccx(c,i)=pnewc(i)
            do 551 j=1,pnewc(i)
                l=l+1
                LPccx(c,l)=lpnewc(l)
551          continue
552       continue

c         if(its.eq.0.and.c.ge.2) then
c         gsum=gsum+g
c         csum=csum+c
c         csum2=csum2+c*c
c         gcsum=gcsum+c*g
c         nsum=nsum+1
c         end if
         if(c.gt.maxc) maxc=c
         
 
        end if
        end if 
        
        npt=npt +1 

c  cumulate gfreq bins - long array in  100 segment blocks
c   c=1, 1,...100
c   c=2, 101,....200 
c        
         if(its.eq.0.and.Q.eq.1.and.p(1).eq.1) then
          nhold=nhold+1
          ghold(nhold)=g
        else
c        if(c.gt.1) then
        ibinc=min(int(g/binw)+1,100)
        ibin=(c-1)*100 + ibinc
        if(c.le.20) gfreq(ibin)=gfreq(ibin) +1 
        end if 
c   complexity penalise 
           
         g=g-beta*c
         
c
c    ge or gt ? I.e. first or last found?
            if(g.ge.gopt) then

c            if(g.gt.gopt) then
                gopt=g
                qopt=qnew
                popt=pnew
                lpopt=lpnew
                
                qcopt=qnewc
                pcopt=pnewc
                lpcopt=lpnewc

                idopt=id
           end if
 

1166  CONTINUE
C end of loop over ID
c----------------------------------------------------------------
C
C   continue generating for different cut levels unless all
C   variables are binary, or levels are fixed.
C
156   GO TO 151
C  now go back to 151 to continue  generation.
C------------------------------------------------------------------
201   CONTINUE
c  at the start iteration=0 set beta (up to this point working beta =0)
         if(its.eq.0.and.Q.eq.1.and.p(1).eq.1) then
c  arbitrarily set bin widhth to cumulative save frequencies
c  of partitions at each complexity             
             binw=gopt/20
             

             if(betain.lt.0) then 
c  arbitrary  fraction of MSE of a single attribute to set beta
c   with betain=NA  (set to -1 in fortran)                  
             beta=0.03*gopt
             gopt=gopt-beta
c             write(*,*)'beta set on iteration 0=', beta
             else
c   input beta giving overpenalised, abort run                 
             if(gopt-beta.lt.0) then
                 status=2
                 return
             end if
             end if
             end if 
c
C  Test whether complexity reached or limits on P1,..,PQ
c
         CALL WINDUP(MORE,Q,QO,P,PO,U)


c&&          write(*,*) 'New Size q=', q, ' p(1)..p(q)=', (p(j),j=1,q)
C-----------------------------------------------------------------------



         IF(MORE) GO TO 155

C-----------------------------------------------------------------------
C   **************************************************
C   *   200 MARKS END OF SEARCHES FOR THIS ITERATION*
C   **************************************************
C  Scanning; output and skip over ranking partitions etc
200   continue

c--------------------------------------------------------------------
              
       its=its+1
       
 
c----------------------------------------------------------------------       
c  need to save the optimal partition 
c   note m+2 in naddx. naddx is index of the expanded partitions
c   it is addressed in EXPAND as  naddx(lp()). As lp()>=2
c   and m+1th is an added attributes. Need to add 1, giving m+1+1=m+2       
        naddx(m+2)=its
c         write(*,*)'m+its+1', m+its+1  
        if(idopt.eq.1)  then 
        qx(its)=qopt
         l=0
         do  548 i=1,qopt
           px(its,i)=popt(i)
           do 549 j=1,popt(i)
           l=l+1
          lpx(its,l)=lpopt(l)
549       continue
548       continue
         
          qxc(its)=qcopt
          l=0
         do  547 i=1,qcopt
         pxc(its,i)=pcopt(i)
           do 546 j=1,pcopt(i)
           l=l+1
            lpxc(its,l)=lpcopt(l)
546       continue
547       continue
c       call outit('optimal idopt=1',qopt,popt,lpopt)         
      else
          
        qx(its)=qcopt
         l=0
         do  545 i=1,qcopt
           px(its,i)=pcopt(i)
           do 543 j=1,pcopt(i)
           l=l+1
          lpx(its,l)=lpcopt(l)
543       continue
545       continue
          qxc(its)=qopt
          l=0
         do  541 i=1,qopt
         pxc(its,i)=popt(i)
           do 544 j=1,popt(i)
           l=l+1
            lpxc(its,l)=lpopt(l)
544       continue
541      continue
         
c       call outit('optimal idopt=0',qcopt,pcopt,lpcopt)
       end if
c&&      write(*,*)'gopt=',gopt
c---------------------------------------------------------------
c   test for same partition as previous iteration, convergence        
       if(idpart(IDopt,qopt,popt,lpopt,
     & IDoptv,qoptv, poptv,lpoptv)) then
c&&         write(*,*)'Converged on iteration ', its
c**       call intpr ("converged iteration", -1,its,1)
 
         go to 7777
         end if
c  or 10 iterations  reached   
       if(its.gt.10) then
c&&         write(*,*)'No convergence by iteration ', its
        status=1
         go to 7777
          end if
c
c  keep current optimal to test for convergence on next iteration

       call exchng(Qopt, Popt ,LPopt ,
     & Qoptv, idummy, Poptv ,LPoptv)
       IDoptv=IDopt
c   ensure attribute added to search set       
       k0=m+1
c&&       write(*,*)'new iteration', its
c  start a new iteration
       go to 1557   
c--------------------------------------------------------------------------------
 
7777   continue
c       write(*,*)'idopt G opt=', idopt, Gopt
c&&       call outit('optimal',qopt,popt,lpopt)
c&&        write(*,*)'Total No. of partitions generated =',npt

       do 539 c=1, maxc
        cx(c)=c
c&&        write(*,*)'c,G', c,gcmplx(c)
539    continue
       
c   compute complexity hull  
c&&&        call intpr ("enter hull maxc", -1,maxc,1)
  
       
         call hull(cx,Gcmplx,maxc,hc,hg,nh,nbn)
         gopt=0
c         write(*,*)'Beta penalising =', beta
c         write(*,*) 'Convex hull:'

c         write(*,*)'c   G           Gbeta'
c&&&        call intpr ("exit hull nh", -1,nh,1)
 
          do  538 j=1,nh
             
             gbeta=hg(j)-beta*hc(j)
           if(gbeta.gt.gopt) then
             gopt=gbeta    
             jopt=j    
             copt=nbn(j)
             end if 
                          
538       continue
 
         l=0
c  extract the optimal from the hull  to be output 
         
         If(copt.le.0)  go to 123
         qopt=Qcmpx(copt)
         do 536 i=1,qopt
             popt(i)=Pcmplx(copt,i)
             do 537 j=1,popt(i)
            l=l+1
            lpopt(l)=LPcx(copt,l) 

537       continue
536       continue
         
         l=0
         
                 qcopt=Qccx(copt)
             do 535 i=1,qcopt
             pcopt(i)=Pccx(copt,i)
             do 534 j=1,pcopt(i)
            l=l+1
            lpcopt(l)=lPccx(copt,l) 

534       continue
535       continue

123       continue
c          write(*,*)'Optimal is complexity c=',copt
c&&        call intpr ("Optimal c", -1,copt,1)
 
         if(idopt.eq.1) then
             call outtoA(combA, qopt,popt,lpopt, MA, lA)
             
             call outtoA(combAc,qcopt,pcopt,lpcopt, MAC,lAc)
             
             else
           
            
             call outtoA(combA,qcopt,pcopt,lpcopt, MA,lA)
             call outtoA(combAc,qopt,popt,lpopt, MAC,lAc)
          end if 
        betain=beta
c      write(*,*) 'gopt exit=',gopt
c                       write(*,*)'frequency bins'
c need to add back in  first iteration p(1)=1,q=1 
cc                       write(*,*)'m, nhold=',m, nhold
                       do 222 im=1,m
                      ibin=min(int(ghold(im)/binw)+1,100)
                      gfreq(ibin)=gfreq(ibin) +1
222                   continue                      
                           
c        sum=0
c        do 333 ic=1,10
c        write(*,*)'c=',ic
c        do 333 jc=1,100
c        if(gfreq((ic-1)*100+ jc).gt.0) then
c            write(*,*)jc, binw*jc, gfreq((ic-1)*100+ jc)
c            sum=sum +gfreq((ic-1)*100+ jc)
c            end if
c333   continue  

       return
      end
      
c-------------------------------------------------------------
      SUBROUTINE  WINDUP(MORE,Q,QO,P,PO,U)
      IMPLICIT INTEGER (A-Z)
      DIMENSION P(*),PO(*),U(*)
      LOGICAL MORE
      MORE=.FALSE.
C
C  a subroutine to wind "up to" an order QO subordinate PO
C  partition
C  
      R=Q 
4902  iq=Q+1-R
      IF(P(iq).EQ.U(iq)) THEN
        IF(R.EQ.1) GO TO 211
        R=R-1
        GO TO 4902
C     
        ELSE
          IF(R.EQ.Q) THEN
          P(iq)=P(iq)+1
          ELSE            
          JP=P(iq)+1
            DO 735 J=R,Q
735         P(Q+1-J)=JP
          END IF
        END IF
        MORE=.TRUE.
        GO TO 789
211   IF(Q.GE.QO) GO TO 789
      Q=Q+1
c      CQ=C+1-Q  
c
c   U is upper "upto" on suborder
c
        DO 215 R=1,Q
        U(R)=PO(R)
215     P(R)=1
      MORE=.TRUE.
789   CONTINUE
      END
C-----------------------------------------------------------------
      SUBROUTINE PARTI(N, Q, NS, SEQ, START, FINI, NOEMB)
C
C  An algorithm to generate the  ways of selecting Q groups
C  of objects from a collection of N objects.
C
      IMPLICIT INTEGER(A-Z)
C
C  There are two vectors passed as subroutine arguments :
C
      DIMENSION SEQ(*), NS(*)
C
C There are 4 locally dimensioned vectors :
      common/partix/ D(10),C(10),CUMNS(11),KEY(10)
      LOGICAL START ,FINI ,EQUL ,NOEMB ,EMBED
      INTEGER NCM
c      character*80 quest
      data one/1/
      IF(.NOT.START) GO TO 200
      START=.FALSE. 
C  First call to PARTI.
C Form cumulative array CUMNS for NS.
C     

      CUMNS(1)=0
      DO 80 I=1,Q
80    CUMNS(I+ 1)=CUMNS(I) + NS(I)
C
      IF(NOEMB) GO TO 84
C
C Embedding is allowed. Set upper key values in C.
C
      DO 85 I=1,Q
      IF(N.LT.0. OR. NS(I).LT.0. OR. N.LT.NS(I)) GO TO 700
      C(I)=NCM(N, NS(I))
      D(I)=0
85    CONTINUE
      GO TO 102
C
C  Embedding is not allowed. Determine binary indicator D such that
C  D(I)=0 if NS(I) NE NS(I-1); =1 if NS(I)=NS(I-1). Also EQUL=.TRUE.
C  only if all NS's are equal.
C
84    I=1
      J=0
      D(1)=0
      EQUL=.TRUE.
99    IJ=I+J + 1
      IF(IJ.GT.Q) GO TO 100
      IF(NS(I).NE.NS(IJ)) GO TO 109
      D(IJ)=1
      J=J + 1
      IF(IJ.EQ.Q) GO TO 100
      GO TO 99
109   IF( NS(I). LT. NS(IJ) ) GO TO 700
      EQUL=.FALSE.
      D(IJ)=0
100   IF(N.LT.0. OR. NS(I).LT.0. OR. N.LT.NS(I)) GO TO 700
      NC=NCM(N, NS(I))
        DO 101 LS=0, J
101     C(I+LS)=NC - J + LS
      IF(IJ.GT.Q) GO TO 102
      I=IJ
      J=0
      GO TO 99
C
C Initialise the combination key, using the binary indicator D.
C  If the keys exceed upper bound, there are no combination, return
C  with FINI=.TRUE.
C
102   KEY(1)=1
c**      write(*,*) 'D(I)',(D(I),I=1,Q)
c**      write(*,*) 'C(I)',(C(I),I=1,Q)
      DO 103 I=2,Q
      IF(D(I). EQ. 0) GO TO 119
      KEY(I)=KEY(I-1) + 1
      IF(KEY(I). LE. C(I)) GO TO 103
      FINI=.TRUE.

      RETURN
119   KEY(I)=1
103   CONTINUE                

C  Determine the combinations corresponding to the initial keys.
C
      DO 478 I=1,Q
      CALL UNDO(N, NS(I), KEY(I), SEQ, CUMNS(I) + one)
478   CONTINUE
      IF(.NOT.NOEMB. OR. EQUL. OR. Q.EQ.1) RETURN
C
C  Test for embedding when none of above conditions hold.
C       
C*      write(*,*)'initial SEQ',(SEQ(I),I=1,7)
      DO 480 I=2,Q
481   CONTINUE

      IF(.NOT.EMBED(I, NS,Q, SEQ, CUMNS)) GO TO 480
      KEY(I)=KEY(I)+ 1
         IF(KEY(I).GT.C(I)) GO TO 700
      CALL UNDO(N, NS(I), KEY(I), SEQ, CUMNS(I)+ one)
      GO TO 481
480   CONTINUE   
c      write(*,*)'initial SEQ',(SEQ(I),I=1,7)
c     write(*,*)'equl noemb', equl ,noemb
      RETURN
C-----------------------------------------------------------------------
C
C  Past first call to PARTI. Increment the keys and unlock the
C  combinations.
C
200   IF(NOEMB) GO TO 250
      I=Q
C
C Embedding is allowable. Simulate nested DO loops subject only
C to the constraint 1<= KEY(I) <=C(I) for all I.
C  All combinations have been generated when KEY(1)=C(1).
C
261   IF(KEY(I).NE.C(I)) GO TO 600
          IF(I.NE.1) GO TO 601
          GO TO 700
601   I=I-1
      GO TO 261
600   KEY(I)=KEY(I) + 1
      CALL UNDO(N, NS(I), KEY(I), SEQ, CUMNS(I)+ one)
      IF(I.EQ.Q) RETURN
          DO 606 II=I+ 1,Q
          KEY(II)=1
         CALL UNDO(N, NS(II), KEY(II), SEQ, CUMNS(II) + one )
606      CONTINUE
       RETURN
C
C  Embedding is not allowable. Separate th cases when EQUL=.TRUE. and
C  EQUL=.FALSE.
C
250   IF(.NOT.EQUL) GO TO 270
      I=Q
C
C  All NS are equal and no embedding permitted. Simulate nested
C  DO loops with KEY(1)<KEY(2)<...<KEY(Q).
C
251   IF(KEY(I).NE.C(I)) GO TO 500
          IF(I.NE.1) GO TO 501
          FINI=.TRUE.
          RETURN
501   I=I-1
      GO TO 251
500   KEY(I)=KEY(I) + 1
      CALL UNDO(N, NS(I), KEY(I), SEQ, CUMNS(I)+ one)
      IF(I.EQ.Q) RETURN
          KY=KEY(I)
          DO 505 II=I+ 1,Q
          KY=KY + 1
          KEY(II)=KY
         CALL UNDO(N, NS(II), KEY(II), SEQ, CUMNS(II) + one)
505      CONTINUE
      RETURN
C
C  Different NS's and no embedding. Simulate nested DO loops with
C  the constraint KEY(1) < D(2)*KEY(2) < D(3)*KEY(3)... < D(Q)*KEY(Q).
C  Test for embedding at each incremented KEY
C
270   I=Q
201   IF(KEY(I).NE.C(I)) GO TO 300
          IF(I.NE.1) GO TO 301
          FINI=.TRUE.
          RETURN
301   I=I-1
      GO TO 201
300   KEY(I)=KEY(I) + 1
      CALL UNDO(N, NS(I), KEY(I), SEQ, CUMNS(I)+ one)
      IF(EMBED(I, NS,Q, SEQ, CUMNS)) GO TO 270
      IF(I.EQ.Q) RETURN
         KY=KEY(I)
         DO 303 II=I+ 1,Q
         IF(D(II). EQ. 0) GO TO 56
         KY=KY + 1
         KEY(II)=KY
         GO TO 306
56       KEY(II)=1
         KY=1
306      CALL UNDO(N, NS(II), KEY(II), SEQ, CUMNS(II) + one)
         IF(EMBED(II, NS, Q, SEQ, CUMNS)) GO TO 270
303      CONTINUE
      RETURN
700   FINI=.TRUE.
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION EMBED(I, NS, Q, SEQ, CUMNS)
C
C  To test whether the Ith combination of SEQ is embedded in at least
C  one of the first, second,.., I-1 th combinations held in SEQ.
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL EMBED
      DIMENSION NS(*), SEQ(*), CUMNS(*)
      I1=CUMNS(I)
      PI=NS(I)
      TOP=SEQ(I1+PI)
      BOT=SEQ(I1+ 1)
      DO 501 J=1,I-1
C
C  First  test for overlapping. If they don't no further testing is
C  needed.
C
        J1=CUMNS(J)
        IF( BOT .GT. SEQ(J1 + NS(J))) GO TO 501
        IF(SEQ(J1+ 1).GT. TOP )  GO TO 501
C
C  Sequences overlap and are possibly embedded. Do an exhaustive
C  to see whether sequence I is contained in J.
C                        
           DO 10 S=1,NS(I)
           SQ=SEQ(I1 + S)
           DO 11 R=1,NS(J)
11         IF(SQ . EQ. SEQ( J1 + R)) GO TO 10
           GO TO 501
10         CONTINUE
C
C  Exit double loop implies that embedding is .TRUE.
C
      EMBED=.TRUE.
      GO TO 506
501   CONTINUE
C  Exit loop implies no embedding found
502   EMBED=.FALSE.
506   RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE UNDO(N, M, K, COMBI, J1)
      implicit integer(a-z)
C
C To unlock the key K to obtain the combination COMBI.
C
      DIMENSION COMBI(*)
      MD=M-1
      ND=N
      KD=K
      J=J1
      JUP=J1 + MD
      CUMI=0
C
102   I=0
      C1=0
      R=ND-MD
      R0=NCM(ND, MD)
C
101   CONTINUE
      R1=(R0*( R - I))/(ND - I )
      C=R1 + C1

      IF(KD.GT.C) GO TO 100
      CUMI=CUMI+I+1
      COMBI(J)=CUMI
      IF( J . EQ. JUP) RETURN
         J= J + 1
         MD= MD - 1
         ND= ND - I - 1
         KD= KD - C1
      GO TO 102
100   I=I + 1
      C1=C
      R0=R1
      GO TO 101
      END
C---------------------------------------------------------------------
      FUNCTION NCM(N, M)
      IMPLICIT INTEGER (A-Z)
      INTEGER NCM
C
C  To compute  the number of ways of selecting M objects from N.
C
      NCM=1
      D=N-M
      DO 101 I=1,M
101   NCM = (NCM *(D + I))/I
      RETURN
      END
c-----------------------------------------------------------
      subroutine exchng(q,p,lp,qn,vn,pn,lpn)
      implicit integer (a-z)
      dimension p(*),lp(*),pn(*),lpn(*)
c
c  to copy  partition definition q,p, lp into qn, pn, lpn.
c
       v=0
       do  539 i=1,q
        v=v+p(i)
539       continue
       qn=q
       vn=v
       do 1 i=1,q
1      pn(i)=p(i)
       do 2 i=1,v
2      lpn(i)=lp(i)
       return
      end  
c------------------------------------------------------------------      
       subroutine trans(q,p,lp,v,qc,pc,lpc,vc,ifail)

  
C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS

      PARAMETER (MAXP=100, maxadd=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================



      character*50 error


      dimension p(*), lp(*), pc(*), lpc(*)
      integer p, pc, q , qc,vc,v
      integer ncall, ncallx
	logical temp

C  local array jp...
      dimension jp(maxp)
	common/debug/bug
	logical bug
c  FUDGE FACTOR TO DEAL WITH LP()=1 REPRESENTING  NULL
      DO 539 j=1,v
      LP(j)=LP(j)+1
539       continue
      

     
	temp=.false.
c	temp  temporary debugger
    
      ifail=0

c  a subroutine to transpose a partition and reduce it
c
c  if single attribute - trivial 
c
      if(q.eq.1.and.p(1).eq.1) then
      qc=1
      pc(1)=1
      lpc(1)=lp(1)
      vc=1
c	write(*,*)'returning trans lpc(1)', lpc(1)
	go to 9999
      end if

      if(q.le.0)then
c  null input partition return fail=3

	qc=1
	lpc(1)=1
	pc(1)=1
      ifail=3
c	write(*,*)'returning null lpc(1)', lpc(1)

      go to 9999
      end if

	jp=1 
c
c  simulate nested do loops to obtain "cross-products"
c  reduce to simplest form as each new cross-product group
C  is created.
c 
      vc=0              
c initialise with one group  qc=1...
      ncall=0
      ncallx=1
1000      k=0
      do 210 i=1,q
      ncallx=ncallx*p(i)
      vc=vc+1
      lpc(vc)=lp(k+jp(i))
210   k=k+p(i)

	if(ncallx.gt.50000) then
c&& 	write(*,*)'Too many transpose operations :', ncallx
	error='ERROR: prohibitive transpose'
c	call help(error)
	ifail=1
	go to 9999
	end if 

      qc=1
      pc(qc)=q 
C
C   now begin simulating nested loops..
c
262   i=q
261   if(jp(i).ne.p(i)) go to 600
        if(i.ne.1) go to 601
        go to 700
601   i=i-1  
      go to 261
600   k=0 
         jp(i)=jp(i)+1
         do 300 ii=i+1,q
300      jp(ii)=1
c
c  form new group
c


         do 200 j=1,q
         vc=vc+1
C  test for array bounds....
            if(vc.gt.maxp2) then 
            ifail=1
            go to 9999
            end if
         lpc(vc)=lp(k+jp(j))
200      k=k+p(j)


         qc=qc+1
C  test for array bounds....
C
           if(qc.gt.maxp) then
           ifail=1
           go to 9999
           end if
         pc(qc)=q
c
C  invoke reduction routine

c	if(bug)	call outit('into reduce in trans qc',qc,pc,lpc)

      call reduce(qc,pc,lpc,vc,ifail)

      
c	if(bug)	call outit('out of reduce in trans qc',qc,pc,lpc)

	if(ifail.eq.1) go to 9999

	if(ifail.lt.0) then
c   return from reduce with null set. But need to
c  fix up by resetting qc=0
c	write(*,*) 'exit null/unic',lpc(1)	
      qc=0
	vc=0
	go to 262
	end if
	        IF(IFAIL.EQ.1.OR.IFAIL.EQ.10)THEN

			go to 9999
			END IF
 
        if(ifail.eq.1) RETURN  
c        indicates q too big, q=maxp
        ncall=ncall+1
C
        if(ncall.gt.ncallx) then
        error='ERROR: problem in transposition algorithm'
c        call help(error) 
        ifail=1
        go to 9999
        end if

      go to 262
c
700    continue


      call reduce(qc,pc,lpc,vc,ifail)


c	if(bug)      call outit('in qc exit trans',qc,pc,lpc) 
      
9999  CONTINUE
c   DEFUDGE!
      do 531 j=1,v
          LP(j)=LP(j)-1
531       continue
      do 529 j=1,vc
          LPC(j)=LPC(j)-1
529       continue

      return
      end
c-----------------------------------------------------------------

      subroutine reduce(q,p,lp,v,ifail)

 
C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS
C

      PARAMETER (MAXP=100, maxadd=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================

      character*80 error

      common/wk1/pold(maxp)
      common/wk2/lpold(maxp2)
      common/wk3/drop(maxp)
      common/wk4/dropsb(maxp2)
      dimension p(*), lp(*)
      logical drop, dropsb, change ,  uni
	logical nulls(maxp), allnull
      logical coef, same, resu, temp
      integer q,p,qold,pold,e,f,v, ae, af
	logical changed

      ifail=0
c  temp is temporary debugging output switch
	
      temp=.false.
	if(temp) then
	
	call outit('entering reduce',q,p,lp)
	end if 

	nulls=.false.

c
c  return immediately if q=1, p(1)=1 no reduction possible
      if(q.eq.1.and.p(1).eq.1) then
      ifail=0

c	if(lp(1).eq.1) ifail=-1   null set
c     if(lp(1).eq.-1)ifail=-2   universal set

	v=1
      return
      end if

c--------------------------------------------------------------

c  initialise drop arrays
      its=0 
      uni=.false. 
1000  continue
	change=.false.
	if(temp)then
c&& 	write(*,*)'cycle ', its 
      call outit('in reduce',q,p,lp)
	end if  
      its=its+1
      if(its.gt.100) then
      error=
     *'Problem in Boolean algebra reduction (Search:End will terminate)'
      ifail=10
c	write(*,*)error
      return
      end if                                
      if(q.gt.maxp)  then
      ifail=1
      return
      end if
c  initialise dropping vectors
	l=0
      do 139 i=1,q
      drop(i)=.false.                                          
        do 139 j=1,p(i)
        l=l+1  
	dropsb(l)=.false.

139   continue  
c----------------------------------------------------
c  begin loop over   group i
c-------------------------------------------------------------
      k=0
      do 100 i=1,q
      If(drop(i)) go to 100
c  i indexes "test" string


c-------------------------------------------------------------
c  begin sub-loop over   group ii
c-------------------------------------------------------------

       kk=0
       do 101 ii=1,q  
   
       if(drop(ii)) go to 101
c  ii  indexes "comparison" string
c  test whether all elements of test string are in comparison
c  string. If so, comparison string is dropped.

c=============================================================== 
          if(i.ne.ii) then


c-------------------------------------------------------------
c  begin sub-sub-loop over   group j
c-------------------------------------------------------------
            
          do 102 j=1,p(i)
c          if(dropsb(k+j)) go to 102
          e=lp(k+j)
c          ae=abs(e)

             if(p(ii).eq.0) go to 102   
c             shouldn't occur   

c-------------------------------------------------------------
c  begin sub-sub-sub-loop over   group jj
c-------------------------------------------------------------


             do 103 jj=1,p(ii)
c             if(dropsb(jj+kk)) go to 103
             f=lp(kk+jj)
c             af=abs(f)


            if(e.eq.f) go to 102
c
c------------------------------------------------------
c      data driven reductions  of string i ne ii
c------------------------------------------------------

103          continue     
c  of              jj=1,p(ii)              
c  exit inner loop over jj implies e, in test string i, is not
c  found in comparison string ii. Comparispon string not nested

         drop(ii)=.false.

          go to 111
102       continue        
c  of          j=1,p(i)


c  exit from outer loop implies every element of test string 
c  in comparison string
c				 write(*,*)'drop(ii) ii exit', ii    

          drop(ii)=.true.
	go to 2030

111       continue
c=============================================================== 

          else 
c              else of if(i.ne.ii)
c=============================================================== 
c i=ii, test for repeats within a string
c
c         write(*,*) 'i,ii at else', i,ii
         do 104 j=1,p(i)
         if(dropsb(k+j)) go to 104
           e=lp(k+j)
c           ae=abs(e)


              do 105 jj=j+1,p(ii)
              jjkk=jj+kk
              if(dropsb(jjkk)) go to 105
              f=lp(jjkk) 
c              af=abs(f)


              if(e.eq.f) then
c              write(*,*)'e.eq.f dropsb(jjkk)=.true.',e,f,jjkk
              dropsb(jjkk)=.true.
			go to 2030
c
c              else  
c   else of else of if(e.eq.f)

         
c                 if(ae.eq.af) then
c   e ne f but abs(e)=abs(f) i.e. attribute and complement
c                  drop(ii)=.true. 
c				nulls(ii)=.true.
c                  go to 2030
c                 end if  
c------------------------------------------------------
c      data driven reductions  within string i=ii
c------------------------------------------------------
c       if(temp) then
c&&       write(*,*)'i,ii', i,ii
c&&       write(*,*)' j,jj',j,jj
c      end if
 
           end if       
c        if(e.eq.f
105       continue
104       continue		
       end if         
c  end of       if(i.ne.ii)
101     kk=kk+p(ii)
100     k=k+p(i) 
c  

c       write(*,*)'at 2030 i,ii,', i,ii
2030   continue

c  now carry out the reduction, first retain q, p, lp
c   option to print out the droppings
	if(temp) then 
c&& 	write(*,*)'big drop sequence'
	ll=0
	do 539 i=1,q
c&& 	write(*,*)'drops :',drop(i),(dropsb(l),l=ll+1,ll+p(i))

	ll=ll+p(i)
539       continue
	end if

c
	   qold=q
       l=0
       do 203 i=1,q
       pold(i)=p(i)
         do 203 j=1,p(i)
         l=l+1
	
203      lpold(l)=lp(l)
c
      q=0
      l=0
      k=0 
	allnull=.true.
      do 200 i=1,qold
cwrite(*,*)'drop(i)', drop(i), i
	if(.not.nulls(i))allnull=.false.
      if(drop(i))go to 200
      q=q+1
      p(q)=0
        do 201 j=1,pold(i)
cwrite(*,*)'dropsb(k+j), k,j',dropsb(k+j), k,j
        if(dropsb(k+j))go to 201
        p(q)=p(q)+1
        l=l+1
        lp(l)=lpold(k+j) 
201     continue
c  correct for possibly empty p( )
      if(p(q).eq.0) then
	q=q-1

	if(q.eq.0) go to 2045

	end if
200   k=k+pold(i)
c   recycle if no change



2045  continue

	

        if(q.ne.qold) then
        change=.true.
	  go to 1000
        else
         do 341 i=1,q
341      if(p(i).ne.pold(i)) change=.true.
         end if

        if(change) go to 1000

c  possibility of change to lp ? really, but is not reduced?
c  i.e.  annulus rule may change lp, but really not worth
c  doing since 
c	call outit('lp',q,p,lp)
c	call outit('lpold',qold,pold,lpold)
	L=0
	do 518 i=1,q
	do 519 j=1,p(i)
	l=l+1
	if(lp(l).ne.lpold(l))change=.true.
519       continue
518       continue
	if(change) go to 1000

c  no change, exit
c  if one cycle only, use IFAIL=2 to indicate input=output
c      if(its.eq.1) ifail=2
      v=0
      do 823 i=1,q
823   v=v+p(i)

8888	continue
      return
      end      
c-------------------------------------------------      
      subroutine outtoA(combo, q,p,lp, Lend,lc)
      integer q,p(*),lp(*),v

      integer combo(*)

c  utility to output partition index strings
      
       L=1
      lx=0
      LL=1
      combo(1)=0
      lc=1
      do 1 i=1,q
      M=L+p(i)-1
      LU=LL+3*p(i)
      
        do 549 j=1,p(i)
            lc=lc+1
            lx=lx+1
            combo(lc)=lp(lx) - 1
549       continue
        
        lc=lc+1
        combo(lc)=0


      LL=LU+2
      l=l+p(i)
1     continue

      Lend=LU+1 
      return
      end
c---------------------------------------------------
      subroutine outit(string,q,p,lp)
      integer q,p(*),lp(*)
      character*(*) string

c  utility to output partition index strings
c      write(*,*)string
       l=1
	if(q.ge.1) then
      do 1 i=1,q
c     write(*,*)'lp=',(lp(j), j=l,l+p(i)-1)

      l=l+p(i)
1     continue
	else
c&& 	write(*,*)'q=0'
	end if
      return
      end
c----------------------------------------------------------------------

      subroutine CRNCH(f,y,n,lp,q,p,lpc,qc,pc, id,freq,
     *skip,ssag,ssa,ssagc,ssac,subg,
     *sy,syc,ssy,ssyc,na,nac,tss,nag,nagc, G)

C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS
C

      PARAMETER (MAXP=100, maxadd=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================


C
C  subroutine to evaluate subgroups and write out diversities
C  on one "side" of a partition.
C
      common/debug/bug
      logical bug
      integer p(*),q,freq(*),lp(*),qc,pc(*),lpc(*)
      integer f(*)
      logical skip,subg
      real y(*)
      integer noin(maxp),noinc(maxp)  
      double precision sy(*), ssy(*),syc(*), ssyc(*), ss, s
      double precision ssa, ssac, ssag, ssagc,tss

c	call outit('lp',q,p,lp)
c	call outit('lpc',qc,pc,lpc)

C  dtaloo loops over data to create partition statistics
        call dtaloo(f,y,n,lp,q,p,lpc,qc,pc, id,freq,
     *  skip,noin,sy,ssy,noinc,syc,ssyc,subg)
       myq=(q+1)
      myqc=(qc+1)
      no=noin(myq)+noinc(myqc)

c   call divi to set diversities from sums of squares
c   on each side of partition
      call divi(q, sy, ssy, noin, subg,ssag, ssa, na, nag)
      call divi(qc,syc,ssyc,noinc,subg,ssagc,ssac,nac,nagc)

c
c  compute overall diversity by adding ss and sums in A and A-
c
      myq=(q+1)
      myqc=(qc+1)
      ss=ssy(myq)+ssyc(myqc)
      s=sy(myq)+syc(myqc)
      no=noin(myq)+noinc(myqc)
      

      xno=float(no)
      if(xno.gt.0) then
      tss=(ss-s*s/xno)
      xm=s/xno
      dv=tss/xno
      else
      tss=0
      xm=0
      dv=0
      end if
      if(dv.le.0) dv=0

      
      nobsin=na+nac         
           IF(NA.eq.0.or.NAc.eq.0) THEN
           WSS=TSS
           BSS=0.0
           G=0.0
           ELSE

           WSS=SSA +SSAc
  
           BSS=TSS-WSS
           G=BSS/NOBSIN

           END IF
           
      return
      end
c--------------------------------------------------------------------
      subroutine divi(q,sy,ssy,noin,subg,ssg,ss,na,nag)

C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS

      PARAMETER (MAXP=100, maxadd=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================


C
C  subroutine to work  out diversities
C  on one "side" of a partition.
C
      integer q,qq
      logical subg
      integer noin(*)
      double precision sy(*), ssy(*)
      double precision ssg, ss, xm
c
      ssg=0
      qq=q+1
      nag=0
C  elements 1,..,q of subgroup arrays not needed if subgroups are
C  not
      l1=1
      if(.not.subg)l1=qq
c
c  qq th element contains the statistics for the q groups combined
c 
      do 1 imy=l1,qq
      na=noin(imy)
      xm=0.0
      if(na.gt.0) then
      xm=sy(imy)/na
      ss=ssy(imy)-xm*sy(imy)
          if(imy.le.q) then
          ssg=ssg+ss   
          nag=nag+na
          end if
      end if
1     continue
      return                     
      end
C----------------------------------------------------------------
      subroutine dtaloo(x,y,n,lp,q,p,lpc,qc,pc,id
     *,freq,skip,noin,sy,ssy,noinc,syc,ssyc,subg)

C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS

      PARAMETER (MAXP=100, maxadd=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================



      integer x(*)
      integer noin(*),noinc(*)
      double precision sy(*), ssy(*), syc(*),ssyc(*)
      double precision yyf
      integer z
      integer freq(*)
      integer q,qq,p(*),lp(*),qc,qqc,pc(*),lpc(*)
      integer t,f
      integer ll, lp4
      logical skip ,subg  , ina
      dimension y(*)

872   format(30i2)
c   initialise over A
 
      qq=q+1
         do 88 i=1,qq
         noin(i)=0
         ssy(i)=0.
88       sy(i)=0
c   initialise over A-
      qqc=qc+1
         do 880 i=1,qqc
         noinc(i)=0
         ssyc(i)=0.
880      syc(i)=0
c
C  loop over observations...
C
c      mn=n*(myx-1)
         mn=0
      ll=0
      do 1 m=1,n     
      ifreq=freq(m)
c      write(*,*)'m,ll=',m,ll
 
      if(ifreq.eq.0) go to 1


c      if(.not.inc(m)) go to 1

c
      ym=y(m)
      if(ym.lt.0) then
c      write(*,*)'observation is .not.inc or freq=0', m
      go to 1
      end if
c      if(skip) then
c
c  missing values handled by skipping ? skip if any missing
      l=0
         do 576 i=1,q
         do 577 j=1,p(i)
         lp4=lp(l+j)
         z=x(ll+lp4)
         if(z.eq.-1) then
             
             
          go to 1
         end if
577      continue
576      l=l+p(i)
c      end if  
c            
c   sort out subgroup membership and increment...
c  is a multiplication worth avoiding? Not sure. Do anyway.
        if(ifreq.eq.1) then
             yf=ym
             else
             yf=ym*ifreq
             end if
      yyf=yf*ym
c  id=1  when  value 1 indicates in a subgroup. 
c  id=0  when value 0 indicates in a subgroup       
      t=id
      f=1-id
      l=0
      a=1-t
      
c      if(ll.eq.0) write(*,*)'z=',Z
c

         do 10 iq=1,q
         b=t
            do 12 j=1,p(iq)
            lp4=lp(l+j)
            z=x(ll+lp4)
            if(z.eq.-1) then
                           
               b=-1 
            else
               if(z.eq.f) go to 10
            end if
12          continue
c           exit inner loop with b either -1 or =t
            if(b.eq.t) then 
             a=t
             if(.not.subg) go to 200
c  increment numbers, sum and ssq  in the subgroups
             noin(iq)=noin(iq)+ifreq
             sy(iq)=sy(iq)+yf
             ssy(iq)=ssy(iq)+yyf
            go to 10
            end if
         if(a.ne.t) a=-1                   
10       l=l+p(iq)
c        exit outer loop with a=-1 or =f or =t
c                     
        
200     ina=.false.
        if(a.ne.-1) then
         if(a.eq.t) then
             noin(qq)=noin(qq)+ifreq
             ssy(qq)=ssy(qq)+yyf
             sy(qq)=sy(qq)+yf
c   in this side of partition, return to skip recall of increm
             ina=.true.  
         end if       
        
      end if          
      if(ina) go to 1


      l=0
c  switch t, f roles
      t=1-id
      f=1-t
c
      a=1-t
c

         do 100 iq=1,qc
         b=t
            do 120 j=1,pc(iq)
            lp4=lpc(l+j)
            z=x(ll+lp4)
            if(z.eq.-1) then

               b=-1 
            else
               if(z.eq.f) go to 100
            end if
120          continue
c           exit inner loop with b either -1 or =t
            if(b.eq.t) then 
             a=t
             if(.not.subg) go to 201
c  increment numbers, sum and ssq  in the subgroups
             noinc(iq)=noinc(iq)+ifreq
             syc(iq)=syc(iq)+yf
             ssyc(iq)=ssyc(iq)+yyf
            go to 100
            end if
         if(a.ne.t) a=-1                   
100       l=l+pc(iq)
c        exit outer loop with a=-1 or =f or =t
c                     
       
201    ina=.false.
        if(a.ne.-1) then
         if(a.eq.t) then
             noinc(qqc)=noinc(qqc)+ifreq
             ssyc(qqc)=ssyc(qqc)+yyf
             syc(qqc)=syc(qqc)+yf
c   in this side of partition, return to skip recall of increm
         end if      
         
 
      end if          
       
1     ll=ll+lc 
c      lc spacing of between observations 1,2,3...n of long data x
      return
      end
c------------------------------------------------------------
      subroutine expand(q,p,lp,qx,px,lpx,qxc,pxc,lpxc,naddx,qnew,
     *pnew,lpnew,vnew,ifail)

      implicit integer (a-z)
 
C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS

      PARAMETER (MAXP=100, maxadd=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================
      integer p(maxp),lp(maxp2), pnew(maxp),lpnew(maxp2),naddx(lc)
      dimension px(maxadd,maxp),lpx(maxadd,maxp2),qx(maxadd)
      dimension pxc(maxadd,maxp),lpxc(maxadd,maxp2),qxc(maxadd)


C  local arrays
      common/wk1/pold(maxp)
      common/wk2/lpold(maxp2)
      
      logical change ,itself

c      logical restrict
	logical pos

      ifail=0
      loop=0
       qold=q         
      l=0
      do 1 i=1,q
       pold(i)=p(i)
       do 1 j=1,p(i)
       l=l+1

1      lpold(l)=lp(l)        
       lold=l

500   l=0
      change=.false.

      qnew=0
      k=0
      do 100 i=1,qold


        do 200 j=1,pold(i)
c  is this element an added factor ?
            
        nox=naddx(lpold(k+j))
c      if(abs(lpold(k+j)).gt.11)write(*,*)   abs(lpold(k+j)),'>11'
c   is added attribute -ve?

c      write(*,*) 'nox,lpold(k+j) ', nox,lpold(k+j)
	pos=.true.
	if(lpold(k+j).lt.0) pos=.false.

c
c  is this an added attribute nox>0.
c  is it also not a data derived added attribute for which
c  save combination is same as itself.
         if(nox.ge.1) then

c       write(*,*)'abs(lpold(k+j))', abs(lpold(k+j))

c	if(restrict) then
	
c	if(nox.eq.n1.or.nox.eq.n2) go to 300

c	else
c	 write(*,*)'lpx(nox,1)', lpx(nox,1)
c	if(pos) then
         itself=( qx(nox).eq.1. .and. px(nox,1).eq.1.and. 
     *          lpx(nox,1).eq.lpold(k+j) )
c	else
c	    itself=( qxc(nox).eq.1. .and. pxc(nox,1).eq.1.and. 
c     *          lpxc(nox,1).eq.lpold(k+j) )
c	end if

        if(.not.itself) go to 300

c      end if

	end if 
200    continue 

c	write(*,*)'no expansion itself=', itself
                                                                   
c exit this loop means element not added, no expansion
      qnew=qnew+1

      if(qnew.gt.maxp) then
      ifail=1
      return
      end if


      pnew(qnew)=pold(i)
         do 301 j=1,pold(i)
         l=l+1

301      lpnew(l)=lpold(k+j)
         go to 100



300   mj=j
      change=.true.
      kx=0  

c	if(pos) then
	qxx=qx(nox)
c	else
c	qxx=qxc(nox)
c	end if
	      
      if(qnew+qxx.gt.maxp) then
c	quest='Array insufficiency in Boolean algebra expansion'//
c     *' Partition ignored'
c	call mssg(quest,1)
c&&       write(*,5712) qnew +qxx, maxp
c5712  format(' WARNING: array insufficiency in EXPAND. Unable to '/
c     *' expand out partition. New no. of subgroups at least ',i4,
c     */' exceeds allowable MAXP=',i4)
      ifail=1
c	call outit('entry',q,p,lp)
c	call outit('warn:', qold,pold,lpold)

	return
      end if


         do 302 ix=1,qxx
         qnew=qnew+1  
c	   if(pos) then   
         pxx=px(nox,ix)
c		else
c	   pxx=pxc(nox,ix)
c		end if


         pnew(qnew)=pold(i)-1 + pxx

           do 303 j=1,pold(i)
           if(j.eq.mj) go to 303
           l=l+1
           lpnew(l)=lpold(k+j)
303        continue


           do 304 jx=1,pxx
           l=l+1
c		if(pos) then
           lpnew(l)=lpx(nox,jx+kx) 
c	else
c           lpnew(l)=-lpxc(nox,jx+kx) 
c	end if
c	write(*,*)'pos=',pos
c	write(*,*)' lpnew(l),lpx(nox,jx+kx)',  lpnew(l),lpx(nox,jx+kx)
c	write(*,*)' lpnew(l),lpxc(nox,jx+kx)', lpnew(l),lpxc(nox,jx+kx)


304        continue

302      kx=kx+pxx


100   k=k+pold(i)
c

      if(.not.change) go to 501 

c  reduce new partition to simplest form 
c       if(qproc)call outit('in reduce from expand',qnew,pnew,lpnew)
       call reduce(qnew,pnew,lpnew,vnew,ifailr)

c       if(qproc)call outit('out reduce from expand',qnew,pnew,lpnew)
cc ifail<0 for null or universal set allowed. Use different "ifail" for reduce


			if(lpnew(1).eq.1) then
c  null or universal set, which??
c			ifail=-1
c	write(*,*) 'id lpnew(l)', id,lpnew(1)
c  this is a questionable fix...?
c			if(id.eq.0)lpnew(1)=-lpnew(1)

			return
			end if

			IF(ifailr.EQ.1.OR.ifailr.EQ.10)THEN

			ifail=ifailr
			return
			END IF
 
c  if new length of element list is unchanged, all
c  expansions are done, else repeat.
      if(loop.gt.20) then
c         write(*,*)'loops >20'
c      call pause('loop>20 in EXPAND exit expansion$')
          ifail=1
      return
      end if
      lold=l
      qold=qnew           
      l=0
       do 4 i=1,qold
       pold(i)=pnew(i)
         do 4 j=1,pold(i)
         l=l+1

4        lpold(l)=lpnew(l)    
c	call outit('at loop',qnew,pnew,lpnew)   
      loop=loop+1
      go to 500
501   continue   



       call reduce(qnew,pnew,lpnew,vnew,ifailr)

		
	ifail=ifailr

c	call outit('exit expand',qnew,pnew,lpnew)    
      return
      end   
c   end of EXPAND
c-----------------------------------------------------------------------      
      subroutine hull(compx,gk,nbest,hc,hg,nh,nbn)
      integer compx(*)
      real gk(*)
      real hc(*)
      real hg(*)
      integer nbn(*)
c  determine pointer to envelope positions on the hull
c
c  include (0,0) as point 1 on the hull

c      eps=abs(  0.001*gk(1) )
      eps=1.0e-20
c      admit points pretty close to hull to be on hull
      i=1

c  revised

      clow=compx(1)
      hlow=gk(1)
      nh=1
      hg(1)=hlow
      hc(1)=clow
      g0=hlow
      c0=clow    
      nbn(1)=1
      
      
2     b=0
      n=0
c  pick out point with steepest slope from position i on      
      do 1 j=i+1,nbest

c	write(*,*)'j,g0,gk(j),c0,compx(j)',j,g0,gk(j),c0,compx(j) 
      bj=(gk(j)-g0)/(compx(j)-c0)
      if((bj-b).gt.eps) then
      n=j
      b=bj
      end if
1     continue
c
c  exit loop, if n=0, must have hit last item
c  
      if(n.eq.0) then
       return
      end if
c   
      nh=nh+1
      nbn(nh)=n   
c      index of which maxima is on hull
      
c      if(n.eq.0) write(*,*)'n=0 nh,compx(n)',nh,compx(n)
      hc(nh)=compx(n)
      hg(nh)=gk(n)
      i=n
      g0=gk(n)
      c0=compx(n)
      go to 2
      end   
c----------------------------------------------------      
      function idpart(id,q,p,lp,idd,qq,pp,lpp)
c
c  logical function to test for identiacl partition combination
c
      implicit integer (a-z)
 
C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS

      PARAMETER (MAXP=100, maxadd=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================


      logical idpart,same2,same
      dimension p(*),pp(*),lp(*),lpp(*)
      logical done(maxp)
      idpart=.false.
      if(idd.ne.id) return 
c**      write(*,*)'q,qq',q,qq
      if(q.ne.qq) return
     
      do 1 j=1,q
1     done(j)=.false.
       l=0
       do 3 j=1,q
           ll=0
c  check each p(j) string elements against those of pp(jj), crossing
c  off (done) each pp(jj) string once it matches
           do 4 jj=1,qq
c           write(*,*)'j,jj,p(j),p(jj)',j,jj,p(j),p(jj)
           if(done(jj)) go to 4
           if(p(j).eq.pp(jj)) then
c   p and pp same lenght, but are same elements ?
           same=same2(lp,l+1,l+p(j),lpp,ll+1,ll+p(j))
c               write(*,*)'lp',(lp(i),i=l+1,l+p(j))
c               write(*,*)'lpp',(lpp(i),i=ll+1,ll+pp(j))
              if(same) then
              done(jj)=.true.
              go to 3
              end if
c	write(*,*)'not same'
           end if
4          ll=ll+pp(jj)
c  exit this loop means no match on remaining element
           idpart=.false.
           return
3       l=l+p(j)
        idpart=.true.
        return
        end
c----------------------------------------------------------------
       function same2(x,l1,l2,y,m1,m2)
c  to test that elements of strings X Y beginning and end locations
c  l1,l2 for x  and m1,m2 for y.
       implicit integer (a-z)

C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS
C
      PARAMETER (MAXP=100, maxadd=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================


       logical done(maxp2)
       logical same2
       integer x(*),y(*)
c
       do 10 j=m1,m2
10     done(j)=.false.
c
       do 1 j=l1,l2
       xj=x(j)
          do 2 jj=m1,m2
          if(done(jj)) go to 2
          if(xj.eq.y(jj)) then
          done(jj)=.true.
          go to 1
          end if
2         continue
c  
c  exit loop 2 means xj.ne. y(jj) for the not done y
c
          same2=.false.
          return
1       continue
         same2=.true.
         return
         end

c==================================================================
      subroutine logrnk(x,y,event,n,lp,q,p,t
     *,freq,skip,na0,na1,G,ncall,mrank,ytemp)

   
      
c          call logrnk(f,y,event,n,IPx,Qnew,Pnew,ID
c     *,freq,skip,na0,na1,G,ncall,mrank,ytemp)
  
C an adaption of FOURF to set up 4 log-rank test
c
c  asssumes data sorted by ascending Y
c
C============================================================
C  PARAMETERS SETTING ARRAY BOUNDS
C
      PARAMETER (MAXP=100, maxadd=10,  MAXCX=100)
      PARAMETER( MAXP2=10000)
      common/dims/lc
c==================================================================



      integer event(*)
      integer x(*)

      integer ll, lp4
      
      integer a,b,z
      integer freq(*)
      integer q,p(*),t,f,lp(*)
      logical skip

      dimension y(*)
      integer mrank(*)
      real ytemp(*)

	double precision ta0, ta1
	double precision E,V, d , risk, o, var
c
c   rank data on first call

      if(ncall.eq.0) then
	n_omit=0
c      write(7,765)
c765    formAt('Note: Ranking data for log-rank calculations')
      do 543 m=1,n
      mrank(m)=m

      ytemp(m)=y(m)

      if(ytemp(m).eq.0) then
	n_omit=n_omit+freq(m)
	end if


543   continue

c	if(n_omit.gt.0) then
c	write(7,674)n_omitted
c674	format('NOTE: ',i4,' zero survival times omitted in log-rank')
c	end if
cccc      call pause('calling sort2$')
      call sort2(n,ytemp,mrank)
      end if
   

      o=0  ! observed
      E=0  ! expected
      V=0  ! variance
      na1=0
      na0=0
c
      n1=0
      n0=0
c	ta0=0
c	ta1=0
c	da0=0
c	da1=0
      ll=0
	sum=0
	d0=0
	d1=0
c      do 104 ii=n,1,-1 
c  make explicit for standard fortran      
      do 104 iidum=1,n
          ii=n-iidum+ 1

	ym=ytemp(ii)
      m=mrank(ii)    
      ll=lc*(m-1)
c??      write(7,*)'i,m,ll,ytemp(m)',i,m,ll,ytemp(m)
c      ms=n*(my-1)+m
c whether ym.lt or ym.le depnds on whether throw out zero survival
c   as n_omitted calculated above. 
c
c      if(.not.inc(m).or.ym.le.0.or.freq(m).eq.0) go to 104
      if(ym.lt.0.or.freq(m).eq.0) go to 104
c	if(event(m).eq.-1) then
c  missing censored variable value?
c	if(x(ll+idc).eq.-1) go to 104
c          go to 104
c	end if

       if(skip) then
c  missing values handled by skipping ? skip if any missing
      l=0
         do 5760 ix=1,q
         do 5770 j=1,p(ix)
         lp4=lp(l+j)
         a=x(ll+lp4)
5770      if(a.eq.-1) go to 200
5760      l=l+p(ix)
      end if


c
	      ifreq=freq(m)
      l=0
      a=1-t
         do 10 i=1,q 
         b=t
            do 12 j=1,p(i)
      lp4=lp(l+j)

      f=1-t
c      if((ip4.gt.1 .and. bit(lp4).eq.-1).or.(ip4.le.0 .and. 
c     *      bit(lp4).eq.1)) f=t

c	if(ip4.eq.-1)f=1  !universal set
c	if(ip4.eq.1)f=0	   !null set


         z=x(ll+lp4)
            if(z.eq.-1) then
               b=-1 
            else
               if(z.eq.f) go to 10
            end if
12          continue
c           exit inner loop with b either -1 or =t
            if(b.eq.t) then 
             a=t
            go to 200
            end if
         a=-1                   
10       l=l+p(i)
c        exit outer loop with a=-1 or =f or =t(when div=true)
c                     
200       continue



c  locate next included observation
	jump=0
8899	jump=jump+1
  

	next=ii-jump
	if(next.eq.0) then
	ynext=ytemp(ii)+87864  !what does 87864 mean???
c  ans: nothing, just ensures ym.ne.ynext in if() below
	go to 8898
	end if
	mnext=mrank(next)
	ynext=ytemp(next)
	if(freq(mnext).eq.0.or.ynext.lt.0) go to 8899



8898	continue

c	write(7,*)'ii,next,jump, ym,ynext',ii,next,jump, ym,ynext

	if(ym.ne.ynext) then

	if(a.eq.0) then
	n0=n0+ifreq
	end if

	if(a.eq.1) then
	n1=n1+ifreq
	end if


      if(event(m).eq.1) then
       
             if(a.eq.0) then
c in  A-
             
             d0=d0+ifreq
c			da0=da0+ifreq

             end if
c  in A
             if(a.eq.1) then
		   
             d1=d1+ifreq
c			da1=da1+ifreq

             end if

	end if !if(event
c
          d=d1+d0
          risk=n0+n1

		if(ym.eq.0) go to 999

          o=o+d1

          if(risk.ge.1) then
          E=E+d*n1/risk
          if(risk.gt.1) then 
          Var=(risk-d)/(risk-1)
     *    * (dfloat(n1)/risk) * (dfloat(n0)/risk) * d
c  following computation gives errors for large n
c  so have broken down computation 15.0.9.09
c          Var=(d*(risk-d)/(risk-1))*((float(n1*n0))/(risk*risk))
          V=V+max(0.0,var)  ! ensure V>0

c	if(var.lt.0) then 
c	write(7,*)'d,risk,n1,n0,var',d,risk,n1,n0,var 
c	end if
c	V=V+var

	          end if
          end if
c
	d0=0
	d1=0
	else  !if(ym.ne.ynext

			if(a.eq.0)then
			n0=n0+ifreq
			end if
			if(a.eq.1)then
			n1=n1+ifreq
			end if

c		if(.not.( censor.and.x(ll+idc).eq.ival) ) then
			
		if( event(m).eq.1 ) then
          if(a.eq.0)then
			d0=d0+ifreq
c			da0=da0+ifreq
			end if
			if(a.eq.1)then
			d1=d1+ifreq
c			da1=da1+ifreq
			end if

		end if

	end if !if(ym.ne.ynext).
999	continue

c  cumulate person-time
c             if(a.eq.0) then
c in  A-
c			ta0=ta0+ym*ifreq
c             end if
c  in A
c			if(a.eq.1) then
c			ta1=ta1+ym*ifreq
c             end if
c
104        ll=ll+lc
	na1=n1
	na0=n0

c	write(7,*)'d, n0, n1, risk',d, n0, n1, risk
c	write(7,*)' o, E, v', o, E, V 
	G=0.0
	if(V.gt.0) G = ((o-E)/V)*(o-E)	
      return
      end !logrank

c-------------------------------------------------
      SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      REAL arr(*)
      INTEGER brr(*),btemp
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,atemp
      if(n.le.0)return 
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        atemp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=atemp
        btemp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=btemp
        if(arr(l+1).gt.arr(ir))then
          atemp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l).gt.arr(ir))then
          atemp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l+1).gt.arr(l))then
          atemp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=btemp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        atemp=arr(i)
        arr(i)=arr(j)
        arr(j)=atemp
        btemp=brr(i)
        brr(i)=brr(j)
        brr(j)=btemp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
c        if(jstack.gt.NSTACK)call pause('NSTACK too small in sort2$')
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END


