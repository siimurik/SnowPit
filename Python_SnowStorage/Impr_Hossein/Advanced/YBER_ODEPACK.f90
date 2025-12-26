!=======================================================================
!                       START OF ODEPACK                                
!=======================================================================
!DECK DLSODA                                                            
      SUBROUTINE DLSODA (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,   &
     &            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT,        &
     &            common_data)                                          
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_interface, only: DINTDY, DROOTS, DSTODA 
      use odepack_common 
      EXTERNAL F, JAC 
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT 
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK 
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW) 
!-----------------------------------------------------------------------
! This is the 12 November 2003 version of                               
! DLSODA: Livermore Solver for Ordinary Differential Equations, with    
!         Automatic method switching for stiff and nonstiff problems.   
!                                                                       
! This version is in double precision.                                  
!                                                                       
! DLSODA solves the initial value problem for stiff or nonstiff         
! systems of first order ODEs,                                          
!     dy/dt = f(t,y) ,  or, in component form,                          
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).    
!                                                                       
! This a variant version of the DLSODE package.                         
! It switches automatically between stiff and nonstiff methods.         
! This means that the user does not have to determine whether the       
! problem is stiff or not, and the solver will automatically choose the 
! appropriate method.  It always starts with the nonstiff method.       
!                                                                       
! Authors:       Alan C. Hindmarsh                                      
!                Center for Applied Scientific Computing, L-561         
!                Lawrence Livermore National Laboratory                 
!                Livermore, CA 94551                                    
! and                                                                   
!                Linda R. Petzold                                       
!                Univ. of California at Santa Barbara                   
!                Dept. of Computer Science                              
!                Santa Barbara, CA 93106                                
!                                                                       
! References:                                                           
! 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE     
!     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),  
!     North-Holland, Amsterdam, 1983, pp. 55-64.                        
! 2.  Linda R. Petzold, Automatic Selection of Methods for Solving      
!     Stiff and Nonstiff Systems of Ordinary Differential Equations,    
!     Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.                 
!-----------------------------------------------------------------------
! Summary of Usage.                                                     
!                                                                       
! Communication between the user and the DLSODA package, for normal     
! situations, is summarized here.  This summary describes only a subset 
! of the full set of options available.  See the full description for   
! details, including alternative treatment of the Jacobian matrix,      
! optional inputs and outputs, nonstandard options, and                 
! instructions for special situations.  See also the example            
! problem (with program and output) following this summary.             
!                                                                       
! A. First provide a subroutine of the form:                            
!               SUBROUTINE F (NEQ, T, Y, YDOT)                          
!               DOUBLE PRECISION T, Y(*), YDOT(*)                       
! which supplies the vector function f by loading YDOT(i) with f(i).    
!                                                                       
! B. Write a main program which calls Subroutine DLSODA once for        
! each point at which answers are desired.  This should also provide    
! for possible use of logical unit 6 for output of error messages       
! by DLSODA.  On the first call to DLSODA, supply arguments as follows: 
! F      = name of subroutine for right-hand side vector f.             
!          This name must be declared External in calling program.      
! NEQ    = number of first order ODEs.                                  
! Y      = array of initial values, of length NEQ.                      
! T      = the initial value of the independent variable.               
! TOUT   = first point where output is desired (.ne. T).                
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.       
! RTOL   = relative tolerance parameter (scalar).                       
! ATOL   = absolute tolerance parameter (scalar or array).              
!          the estimated local error in y(i) will be controlled so as   
!          to be less than                                              
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or        
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.           
!          Thus the local error test passes if, in each component,      
!          either the absolute error is less than ATOL (or ATOL(i)),    
!          or the relative error is less than RTOL.                     
!          Use RTOL = 0.0 for pure absolute error control, and          
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error    
!          control.  Caution: actual (global) errors may exceed these   
!          local tolerances, so choose them conservatively.             
! ITASK  = 1 for normal computation of output values of y at t = TOUT.  
! ISTATE = integer flag (input and output).  Set ISTATE = 1.            
! IOPT   = 0 to indicate no optional inputs used.                       
! RWORK  = real work array of length at least:                          
!             22 + NEQ * MAX(16, NEQ + 9).                              
!          See also Paragraph E below.                                  
! LRW    = declared length of RWORK (in user's dimension).              
! IWORK  = integer work array of length at least  20 + NEQ.             
! LIW    = declared length of IWORK (in user's dimension).              
! JAC    = name of subroutine for Jacobian matrix.                      
!          Use a dummy name.  See also Paragraph E below.               
! JT     = Jacobian type indicator.  Set JT = 2.                        
!          See also Paragraph E below.                                  
! Note that the main program must declare arrays Y, RWORK, IWORK,       
! and possibly ATOL.                                                    
!                                                                       
! C. The output from the first call (or any call) is:                   
!      Y = array of computed values of y(t) vector.                     
!      T = corresponding value of independent variable (normally TOUT). 
! ISTATE = 2  if DLSODA was successful, negative otherwise.             
!          -1 means excess work done on this call (perhaps wrong JT).   
!          -2 means excess accuracy requested (tolerances too small).   
!          -3 means illegal input detected (see printed message).       
!          -4 means repeated error test failures (check all inputs).    
!          -5 means repeated convergence failures (perhaps bad Jacobian 
!             supplied or wrong choice of JT or tolerances).            
!          -6 means error weight became zero during problem. (Solution  
!             component i vanished, and ATOL or ATOL(i) = 0.)           
!          -7 means work space insufficient to finish (see messages).   
!                                                                       
! D. To continue the integration after a successful return, simply      
! reset TOUT and call DLSODA again.  No other parameters need be reset. 
!                                                                       
! E. Note: If and when DLSODA regards the problem as stiff, and         
! switches methods accordingly, it must make use of the NEQ by NEQ      
! Jacobian matrix, J = df/dy.  For the sake of simplicity, the          
! inputs to DLSODA recommended in Paragraph B above cause DLSODA to     
! treat J as a full matrix, and to approximate it internally by         
! difference quotients.  Alternatively, J can be treated as a band      
! matrix (with great potential reduction in the size of the RWORK       
! array).  Also, in either the full or banded case, the user can supply 
! J in closed form, with a routine whose name is passed as the JAC      
! argument.  These alternatives are described in the paragraphs on      
! RWORK, JAC, and JT in the full description of the call sequence below.
!                                                                       
!-----------------------------------------------------------------------
! Example Problem.                                                      
!                                                                       
! The following is a simple example problem, with the coding            
! needed for its solution by DLSODA.  The problem is from chemical      
! kinetics, and consists of the following three rate equations:         
!     dy1/dt = -.04*y1 + 1.e4*y2*y3                                     
!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2                         
!     dy3/dt = 3.e7*y2**2                                               
! on the interval from t = 0.0 to t = 4.e10, with initial conditions    
! y1 = 1.0, y2 = y3 = 0.  The problem is stiff.                         
!                                                                       
! The following coding solves this problem with DLSODA,                 
! printing results at t = .4, 4., ..., 4.e10.  It uses                  
! ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because           
! y2 has much smaller values.                                           
! At the end of the run, statistical quantities of interest are         
! printed (see optional outputs in the full description below).         
!                                                                       
!     EXTERNAL FEX                                                      
!     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y                    
!     DIMENSION Y(3), ATOL(3), RWORK(70), IWORK(23)                     
!     NEQ = 3                                                           
!     Y(1) = 1.                                                         
!     Y(2) = 0.                                                         
!     Y(3) = 0.                                                         
!     T = 0.                                                            
!     TOUT = .4                                                         
!     ITOL = 2                                                          
!     RTOL = 1.D-4                                                      
!     ATOL(1) = 1.D-6                                                   
!     ATOL(2) = 1.D-10                                                  
!     ATOL(3) = 1.D-6                                                   
!     ITASK = 1                                                         
!     ISTATE = 1                                                        
!     IOPT = 0                                                          
!     LRW = 70                                                          
!     LIW = 23                                                          
!     JT = 2                                                            
!     DO 40 IOUT = 1,12                                                 
!       CALL DLSODA(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,       
!    1     IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT)                            
!       WRITE(6,20)T,Y(1),Y(2),Y(3)                                     
! 20    FORMAT(' At t =',D12.4,'   Y =',3D14.6)                         
!       IF (ISTATE .LT. 0) GO TO 80                                     
! 40    TOUT = TOUT*10.                                                 
!     WRITE(6,60)IWORK(11),IWORK(12),IWORK(13),IWORK(19),RWORK(15)      
! 60  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4/      
!    1   ' Method last used =',I2,'   Last switch was at t =',D12.4)    
!     STOP                                                              
! 80  WRITE(6,90)ISTATE                                                 
! 90  FORMAT(///' Error halt.. ISTATE =',I3)                            
!     STOP                                                              
!     END                                                               
!                                                                       
!     SUBROUTINE FEX (NEQ, T, Y, YDOT)                                  
!     DOUBLE PRECISION T, Y, YDOT                                       
!     DIMENSION Y(3), YDOT(3)                                           
!     YDOT(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)                              
!     YDOT(3) = 3.D7*Y(2)*Y(2)                                          
!     YDOT(2) = -YDOT(1) - YDOT(3)                                      
!     RETURN                                                            
!     END                                                               
!                                                                       
! The output of this program (on a CDC-7600 in single precision)        
! is as follows:                                                        
!                                                                       
!   At t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02  
!   At t =  4.0000e+00   Y =  9.055333e-01  2.240655e-05  9.444430e-02  
!   At t =  4.0000e+01   Y =  7.158403e-01  9.186334e-06  2.841505e-01  
!   At t =  4.0000e+02   Y =  4.505250e-01  3.222964e-06  5.494717e-01  
!   At t =  4.0000e+03   Y =  1.831975e-01  8.941774e-07  8.168016e-01  
!   At t =  4.0000e+04   Y =  3.898730e-02  1.621940e-07  9.610125e-01  
!   At t =  4.0000e+05   Y =  4.936363e-03  1.984221e-08  9.950636e-01  
!   At t =  4.0000e+06   Y =  5.161831e-04  2.065786e-09  9.994838e-01  
!   At t =  4.0000e+07   Y =  5.179817e-05  2.072032e-10  9.999482e-01  
!   At t =  4.0000e+08   Y =  5.283401e-06  2.113371e-11  9.999947e-01  
!   At t =  4.0000e+09   Y =  4.659031e-07  1.863613e-12  9.999995e-01  
!   At t =  4.0000e+10   Y =  1.404280e-08  5.617126e-14  1.000000e+00  
!                                                                       
!   No. steps = 361  No. f-s = 693  No. J-s =  64                       
!   Method last used = 2   Last switch was at t =  6.0092e-03           
!-----------------------------------------------------------------------
! Full description of user interface to DLSODA.                         
!                                                                       
! The user interface to DLSODA consists of the following parts.         
!                                                                       
! 1.   The call sequence to Subroutine DLSODA, which is a driver        
!      routine for the solver.  This includes descriptions of both      
!      the call sequence arguments and of user-supplied routines.       
!      following these descriptions is a description of                 
!      optional inputs available through the call sequence, and then    
!      a description of optional outputs (in the work arrays).          
!                                                                       
! 2.   Descriptions of other routines in the DLSODA package that may be 
!      (optionally) called by the user.  These provide the ability to   
!      alter error message handling, save and restore the internal      
!      Common, and obtain specified derivatives of the solution y(t).   
!                                                                       
! 3.   Descriptions of Common blocks to be declared in overlay          
!      or similar environments, or to be saved when doing an interrupt  
!      of the problem and continued solution later.                     
!                                                                       
! 4.   Description of a subroutine in the DLSODA package,               
!      which the user may replace with his/her own version, if desired. 
!      this relates to the measurement of errors.                       
!                                                                       
!-----------------------------------------------------------------------
! Part 1.  Call Sequence.                                               
!                                                                       
! The call sequence parameters used for input only are                  
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, JT,   
! and those used for both input and output are                          
!     Y, T, ISTATE.                                                     
! The work arrays RWORK and IWORK are also used for conditional and     
! optional inputs and optional outputs.  (The term output here refers   
! to the return from Subroutine DLSODA to the user's calling program.)  
!                                                                       
! The legality of input parameters will be thoroughly checked on the    
! initial call for the problem, but not checked thereafter unless a     
! change in input parameters is flagged by ISTATE = 3 on input.         
!                                                                       
! The descriptions of the call arguments are as follows.                
!                                                                       
! F      = the name of the user-supplied subroutine defining the        
!          ODE system.  The system must be put in the first-order       
!          form dy/dt = f(t,y), where f is a vector-valued function     
!          of the scalar t and the vector y.  Subroutine F is to        
!          compute the function f.  It is to have the form              
!               SUBROUTINE F (NEQ, T, Y, YDOT)                          
!               DOUBLE PRECISION T, Y(*), YDOT(*)                       
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)   
!          is output.  Y and YDOT are arrays of length NEQ.             
!          Subroutine F should not alter Y(1),...,Y(NEQ).               
!          F must be declared External in the calling program.          
!                                                                       
!          Subroutine F may access user-defined quantities in           
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array      
!          (dimensioned in F) and/or Y has length exceeding NEQ(1).     
!          See the descriptions of NEQ and Y below.                     
!                                                                       
!          If quantities computed in the F routine are needed           
!          externally to DLSODA, an extra call to F should be made      
!          for this purpose, for consistent and accurate results.       
!          If only the derivative dy/dt is needed, use DINTDY instead.  
!                                                                       
! NEQ    = the size of the ODE system (number of first order            
!          ordinary differential equations).  Used only for input.      
!          NEQ may be decreased, but not increased, during the problem. 
!          If NEQ is decreased (with ISTATE = 3 on input), the          
!          remaining components of Y should be left undisturbed, if     
!          these are to be accessed in F and/or JAC.                    
!                                                                       
!          Normally, NEQ is a scalar, and it is generally referred to   
!          as a scalar in this user interface description.  However,    
!          NEQ may be an array, with NEQ(1) set to the system size.     
!          (The DLSODA package accesses only NEQ(1).)  In either case,  
!          this parameter is passed as the NEQ argument in all calls    
!          to F and JAC.  Hence, if it is an array, locations           
!          NEQ(2),... may be used to store other integer data and pass  
!          it to F and/or JAC.  Subroutines F and/or JAC must include   
!          NEQ in a Dimension statement in that case.                   
!                                                                       
! Y      = a real array for the vector of dependent variables, of       
!          length NEQ or more.  Used for both input and output on the   
!          first call (ISTATE = 1), and only for output on other calls. 
!          On the first call, Y must contain the vector of initial      
!          values.  On output, Y contains the computed solution vector, 
!          evaluated at T.  If desired, the Y array may be used         
!          for other purposes between calls to the solver.              
!                                                                       
!          This array is passed as the Y argument in all calls to       
!          F and JAC.  Hence its length may exceed NEQ, and locations   
!          Y(NEQ+1),... may be used to store other real data and        
!          pass it to F and/or JAC.  (The DLSODA package accesses only  
!          Y(1),...,Y(NEQ).)                                            
!                                                                       
! T      = the independent variable.  On input, T is used only on the   
!          first call, as the initial point of the integration.         
!          on output, after each call, T is the value at which a        
!          computed solution Y is evaluated (usually the same as TOUT). 
!          on an error return, T is the farthest point reached.         
!                                                                       
! TOUT   = the next value of t at which a computed solution is desired. 
!          Used only for input.                                         
!                                                                       
!          When starting the problem (ISTATE = 1), TOUT may be equal    
!          to T for one call, then should .ne. T for the next call.     
!          For the initial t, an input value of TOUT .ne. T is used     
!          in order to determine the direction of the integration       
!          (i.e. the algebraic sign of the step sizes) and the rough    
!          scale of the problem.  Integration in either direction       
!          (forward or backward in t) is permitted.                     
!                                                                       
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after    
!          the first call (i.e. the first call with TOUT .ne. T).       
!          Otherwise, TOUT is required on every call.                   
!                                                                       
!          If ITASK = 1, 3, or 4, the values of TOUT need not be        
!          monotone, but a value of TOUT which backs up is limited      
!          to the current internal T interval, whose endpoints are      
!          TCUR - HU and TCUR (see optional outputs, below, for         
!          TCUR and HU).                                                
!                                                                       
! ITOL   = an indicator for the type of error control.  See             
!          description below under ATOL.  Used only for input.          
!                                                                       
! RTOL   = a relative error tolerance parameter, either a scalar or     
!          an array of length NEQ.  See description below under ATOL.   
!          Input only.                                                  
!                                                                       
! ATOL   = an absolute error tolerance parameter, either a scalar or    
!          an array of length NEQ.  Input only.                         
!                                                                       
!             The input parameters ITOL, RTOL, and ATOL determine       
!          the error control performed by the solver.  The solver will  
!          control the vector E = (E(i)) of estimated local errors      
!          in y, according to an inequality of the form                 
!                      max-norm of ( E(i)/EWT(i) )   .le.   1,          
!          where EWT = (EWT(i)) is a vector of positive error weights.  
!          The values of RTOL and ATOL should all be non-negative.      
!          The following table gives the types (scalar/array) of        
!          RTOL and ATOL, and the corresponding form of EWT(i).         
!                                                                       
!             ITOL    RTOL       ATOL          EWT(i)                   
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL        
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)     
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL     
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)  
!                                                                       
!          When either of these parameters is a scalar, it need not     
!          be dimensioned in the user's calling program.                
!                                                                       
!          If none of the above choices (with ITOL, RTOL, and ATOL      
!          fixed throughout the problem) is suitable, more general      
!          error controls can be obtained by substituting a             
!          user-supplied routine for the setting of EWT.                
!          See Part 4 below.                                            
!                                                                       
!          If global errors are to be estimated by making a repeated    
!          run on the same problem with smaller tolerances, then all    
!          components of RTOL and ATOL (i.e. of EWT) should be scaled   
!          down uniformly.                                              
!                                                                       
! ITASK  = an index specifying the task to be performed.                
!          Input only.  ITASK has the following values and meanings.    
!          1  means normal computation of output values of y(t) at      
!             t = TOUT (by overshooting and interpolating).             
!          2  means take one step only and return.                      
!          3  means stop at the first internal mesh point at or         
!             beyond t = TOUT and return.                               
!          4  means normal computation of output values of y(t) at      
!             t = TOUT but without overshooting t = TCRIT.              
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to   
!             or beyond TOUT, but not behind it in the direction of     
!             integration.  This option is useful if the problem        
!             has a singularity at or beyond t = TCRIT.                 
!          5  means take one step, without passing TCRIT, and return.   
!             TCRIT must be input as RWORK(1).                          
!                                                                       
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT        
!          (within roundoff), it will return T = TCRIT (exactly) to     
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT, 
!          in which case answers at t = TOUT are returned first).       
!                                                                       
! ISTATE = an index used for input and output to specify the            
!          the state of the calculation.                                
!                                                                       
!          On input, the values of ISTATE are as follows.               
!          1  means this is the first call for the problem              
!             (initializations will be done).  See note below.          
!          2  means this is not the first call, and the calculation     
!             is to continue normally, with no change in any input      
!             parameters except possibly TOUT and ITASK.                
!             (If ITOL, RTOL, and/or ATOL are changed between calls     
!             with ISTATE = 2, the new values will be used but not      
!             tested for legality.)                                     
!          3  means this is not the first call, and the                 
!             calculation is to continue normally, but with             
!             a change in input parameters other than                   
!             TOUT and ITASK.  Changes are allowed in                   
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,        
!             and any optional inputs except H0, MXORDN, and MXORDS.    
!             (See IWORK description for ML and MU.)                    
!          Note:  A preliminary call with TOUT = T is not counted       
!          as a first call here, as no initialization or checking of    
!          input is done.  (Such a call is sometimes useful for the     
!          purpose of outputting the initial conditions.)               
!          Thus the first call for which TOUT .ne. T requires           
!          ISTATE = 1 on input.                                         
!                                                                       
!          On output, ISTATE has the following values and meanings.     
!           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
!           2  means the integration was performed successfully.        
!          -1  means an excessive amount of work (more than MXSTEP      
!              steps) was done on this call, before completing the      
!              requested task, but the integration was otherwise        
!              successful as far as T.  (MXSTEP is an optional input    
!              and is normally 500.)  To continue, the user may         
!              simply reset ISTATE to a value .gt. 1 and call again     
!              (the excess work step counter will be reset to 0).       
!              In addition, the user may increase MXSTEP to avoid       
!              this error return (see below on optional inputs).        
!          -2  means too much accuracy was requested for the precision  
!              of the machine being used.  This was detected before     
!              completing the requested task, but the integration       
!              was successful as far as T.  To continue, the tolerance  
!              parameters must be reset, and ISTATE must be set         
!              to 3.  The optional output TOLSF may be used for this    
!              purpose.  (Note: If this condition is detected before    
!              taking any steps, then an illegal input return           
!              (ISTATE = -3) occurs instead.)                           
!          -3  means illegal input was detected, before taking any      
!              integration steps.  See written message for details.     
!              Note:  If the solver detects an infinite loop of calls   
!              to the solver with illegal input, it will cause          
!              the run to stop.                                         
!          -4  means there were repeated error test failures on         
!              one attempted step, before completing the requested      
!              task, but the integration was successful as far as T.    
!              The problem may have a singularity, or the input         
!              may be inappropriate.                                    
!          -5  means there were repeated convergence test failures on   
!              one attempted step, before completing the requested      
!              task, but the integration was successful as far as T.    
!              This may be caused by an inaccurate Jacobian matrix,     
!              if one is being used.                                    
!          -6  means EWT(i) became zero for some i during the           
!              integration.  Pure relative error control (ATOL(i)=0.0)  
!              was requested on a variable which has now vanished.      
!              The integration was successful as far as T.              
!          -7  means the length of RWORK and/or IWORK was too small to  
!              proceed, but the integration was successful as far as T. 
!              This happens when DLSODA chooses to switch methods       
!              but LRW and/or LIW is too small for the new method.      
!                                                                       
!          Note:  Since the normal output value of ISTATE is 2,         
!          it does not need to be reset for normal continuation.        
!          Also, since a negative input value of ISTATE will be         
!          regarded as illegal, a negative output value requires the    
!          user to change it, and possibly other inputs, before         
!          calling the solver again.                                    
!                                                                       
! IOPT   = an integer flag to specify whether or not any optional       
!          inputs are being used on this call.  Input only.             
!          The optional inputs are listed separately below.             
!          IOPT = 0 means no optional inputs are being used.            
!                   default values will be used in all cases.           
!          IOPT = 1 means one or more optional inputs are being used.   
!                                                                       
! RWORK  = a real array (double precision) for work space, and (in the  
!          first 20 words) for conditional and optional inputs and      
!          optional outputs.                                            
!          As DLSODA switches automatically between stiff and nonstiff  
!          methods, the required length of RWORK can change during the  
!          problem.  Thus the RWORK array passed to DLSODA can either   
!          have a static (fixed) length large enough for both methods,  
!          or have a dynamic (changing) length altered by the calling   
!          program in response to output from DLSODA.                   
!                                                                       
!                       --- Fixed Length Case ---                       
!          If the RWORK length is to be fixed, it should be at least    
!               MAX (LRN, LRS),                                         
!          where LRN and LRS are the RWORK lengths required when the    
!          current method is nonstiff or stiff, respectively.           
!                                                                       
!          The separate RWORK length requirements LRN and LRS are       
!          as follows:                                                  
!          IF NEQ is constant and the maximum method orders have        
!          their default values, then                                   
!             LRN = 20 + 16*NEQ,                                        
!             LRS = 22 + 9*NEQ + NEQ**2           if JT = 1 or 2,       
!             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ   if JT = 4 or 5.       
!          Under any other conditions, LRN and LRS are given by:        
!             LRN = 20 + NYH*(MXORDN+1) + 3*NEQ,                        
!             LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT,                 
!          where                                                        
!             NYH    = the initial value of NEQ,                        
!             MXORDN = 12, unless a smaller value is given as an        
!                      optional input,                                  
!             MXORDS = 5, unless a smaller value is given as an         
!                      optional input,                                  
!             LMAT   = length of matrix work space:                     
!             LMAT   = NEQ**2 + 2              if JT = 1 or 2,          
!             LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5.          
!                                                                       
!                       --- Dynamic Length Case ---                     
!          If the length of RWORK is to be dynamic, then it should      
!          be at least LRN or LRS, as defined above, depending on the   
!          current method.  Initially, it must be at least LRN (since   
!          DLSODA starts with the nonstiff method).  On any return      
!          from DLSODA, the optional output MCUR indicates the current  
!          method.  If MCUR differs from the value it had on the        
!          previous return, or if there has only been one call to       
!          DLSODA and MCUR is now 2, then DLSODA has switched           
!          methods during the last call, and the length of RWORK        
!          should be reset (to LRN if MCUR = 1, or to LRS if            
!          MCUR = 2).  (An increase in the RWORK length is required     
!          if DLSODA returned ISTATE = -7, but not otherwise.)          
!          After resetting the length, call DLSODA with ISTATE = 3      
!          to signal that change.                                       
!                                                                       
! LRW    = the length of the array RWORK, as declared by the user.      
!          (This will be checked by the solver.)                        
!                                                                       
! IWORK  = an integer array for work space.                             
!          As DLSODA switches automatically between stiff and nonstiff  
!          methods, the required length of IWORK can change during      
!          problem, between                                             
!             LIS = 20 + NEQ   and   LIN = 20,                          
!          respectively.  Thus the IWORK array passed to DLSODA can     
!          either have a fixed length of at least 20 + NEQ, or have a   
!          dynamic length of at least LIN or LIS, depending on the      
!          current method.  The comments on dynamic length under        
!          RWORK above apply here.  Initially, this length need         
!          only be at least LIN = 20.                                   
!                                                                       
!          The first few words of IWORK are used for conditional and    
!          optional inputs and optional outputs.                        
!                                                                       
!          The following 2 words in IWORK are conditional inputs:       
!            IWORK(1) = ML     these are the lower and upper            
!            IWORK(2) = MU     half-bandwidths, respectively, of the    
!                       banded Jacobian, excluding the main diagonal.   
!                       The band is defined by the matrix locations     
!                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU    
!                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.        
!                       These are required if JT is 4 or 5, and         
!                       ignored otherwise.  ML and MU may in fact be    
!                       the band parameters for a matrix to which       
!                       df/dy is only approximately equal.              
!                                                                       
! LIW    = the length of the array IWORK, as declared by the user.      
!          (This will be checked by the solver.)                        
!                                                                       
! Note: The base addresses of the work arrays must not be               
! altered between calls to DLSODA for the same problem.                 
! The contents of the work arrays must not be altered                   
! between calls, except possibly for the conditional and                
! optional inputs, and except for the last 3*NEQ words of RWORK.        
! The latter space is used for internal scratch space, and so is        
! available for use by the user outside DLSODA between calls, if        
! desired (but not for use by F or JAC).                                
!                                                                       
! JAC    = the name of the user-supplied routine to compute the         
!          Jacobian matrix, df/dy, if JT = 1 or 4.  The JAC routine     
!          is optional, but if the problem is expected to be stiff much 
!          of the time, you are encouraged to supply JAC, for the sake  
!          of efficiency.  (Alternatively, set JT = 2 or 5 to have      
!          DLSODA compute df/dy internally by difference quotients.)    
!          If and when DLSODA uses df/dy, it treats this NEQ by NEQ     
!          matrix either as full (JT = 1 or 2), or as banded (JT =      
!          4 or 5) with half-bandwidths ML and MU (discussed under      
!          IWORK above).  In either case, if JT = 1 or 4, the JAC       
!          routine must compute df/dy as a function of the scalar t     
!          and the vector y.  It is to have the form                    
!               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)          
!               DOUBLE PRECISION T, Y(*), PD(NROWPD,*)                  
!          where NEQ, T, Y, ML, MU, and NROWPD are input and the array  
!          PD is to be loaded with partial derivatives (elements of     
!          the Jacobian matrix) on output.  PD must be given a first    
!          dimension of NROWPD.  T and Y have the same meaning as in    
!          Subroutine F.                                                
!               In the full matrix case (JT = 1), ML and MU are         
!          ignored, and the Jacobian is to be loaded into PD in         
!          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).     
!               In the band matrix case (JT = 4), the elements          
!          within the band are to be loaded into PD in columnwise       
!          manner, with diagonal lines of df/dy loaded into the rows    
!          of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!          ML and MU are the half-bandwidth parameters (see IWORK).     
!          The locations in PD in the two triangular areas which        
!          correspond to nonexistent matrix elements can be ignored     
!          or loaded arbitrarily, as they are overwritten by DLSODA.    
!               JAC need not provide df/dy exactly.  A crude            
!          approximation (possibly with a smaller bandwidth) will do.   
!               In either case, PD is preset to zero by the solver,     
!          so that only the nonzero elements need be loaded by JAC.     
!          Each call to JAC is preceded by a call to F with the same    
!          arguments NEQ, T, and Y.  Thus to gain some efficiency,      
!          intermediate quantities shared by both calculations may be   
!          saved in a user Common block by F and not recomputed by JAC, 
!          if desired.  Also, JAC may alter the Y array, if desired.    
!          JAC must be declared External in the calling program.        
!               Subroutine JAC may access user-defined quantities in    
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array      
!          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).   
!          See the descriptions of NEQ and Y above.                     
!                                                                       
! JT     = Jacobian type indicator.  Used only for input.               
!          JT specifies how the Jacobian matrix df/dy will be           
!          treated, if and when DLSODA requires this matrix.            
!          JT has the following values and meanings:                    
!           1 means a user-supplied full (NEQ by NEQ) Jacobian.         
!           2 means an internally generated (difference quotient) full  
!             Jacobian (using NEQ extra calls to F per df/dy value).    
!           4 means a user-supplied banded Jacobian.                    
!           5 means an internally generated banded Jacobian (using      
!             ML+MU+1 extra calls to F per df/dy evaluation).           
!          If JT = 1 or 4, the user must supply a Subroutine JAC        
!          (the name is arbitrary) as described above under JAC.        
!          If JT = 2 or 5, a dummy argument can be used.                
!-----------------------------------------------------------------------
! Optional Inputs.                                                      
!                                                                       
! The following is a list of the optional inputs provided for in the    
! call sequence.  (See also Part 2.)  For each such input variable,     
! this table lists its name as used in this documentation, its          
! location in the call sequence, its meaning, and the default value.    
! The use of any of these inputs requires IOPT = 1, and in that         
! case all of these inputs are examined.  A value of zero for any       
! of these optional inputs will cause the default value to be used.     
! Thus to use a subset of the optional inputs, simply preload           
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and   
! then set those of interest to nonzero values.                         
!                                                                       
! Name    Location      Meaning and Default Value                       
!                                                                       
! H0      RWORK(5)  the step size to be attempted on the first step.    
!                   The default value is determined by the solver.      
!                                                                       
! HMAX    RWORK(6)  the maximum absolute step size allowed.             
!                   The default value is infinite.                      
!                                                                       
! HMIN    RWORK(7)  the minimum absolute step size allowed.             
!                   The default value is 0.  (This lower bound is not   
!                   enforced on the final step before reaching TCRIT    
!                   when ITASK = 4 or 5.)                               
!                                                                       
! IXPR    IWORK(5)  flag to generate extra printing at method switches. 
!                   IXPR = 0 means no extra printing (the default).     
!                   IXPR = 1 means print data on each switch.           
!                   T, H, and NST will be printed on the same logical   
!                   unit as used for error messages.                    
!                                                                       
! MXSTEP  IWORK(6)  maximum number of (internally defined) steps        
!                   allowed during one call to the solver.              
!                   The default value is 500.                           
!                                                                       
! MXHNIL  IWORK(7)  maximum number of messages printed (per problem)    
!                   warning that T + H = T on a step (H = step size).   
!                   This must be positive to result in a non-default    
!                   value.  The default value is 10.                    
!                                                                       
! MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff    
!                   (Adams) method.  the default value is 12.           
!                   if MXORDN exceeds the default value, it will        
!                   be reduced to the default value.                    
!                   MXORDN is held constant during the problem.         
!                                                                       
! MXORDS  IWORK(9)  the maximum order to be allowed for the stiff       
!                   (BDF) method.  The default value is 5.              
!                   If MXORDS exceeds the default value, it will        
!                   be reduced to the default value.                    
!                   MXORDS is held constant during the problem.         
!-----------------------------------------------------------------------
! Optional Outputs.                                                     
!                                                                       
! As optional additional output from DLSODA, the variables listed       
! below are quantities related to the performance of DLSODA             
! which are available to the user.  These are communicated by way of    
! the work arrays, but also have internal mnemonic names as shown.      
! except where stated otherwise, all of these outputs are defined       
! on any successful return from DLSODA, and on any return with          
! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return           
! (ISTATE = -3), they will be unchanged from their existing values      
! (if any), except possibly for TOLSF, LENRW, and LENIW.                
! On any error return, outputs relevant to the error will be defined,   
! as noted below.                                                       
!                                                                       
! Name    Location      Meaning                                         
!                                                                       
! HU      RWORK(11) the step size in t last used (successfully).        
!                                                                       
! HCUR    RWORK(12) the step size to be attempted on the next step.     
!                                                                       
! TCUR    RWORK(13) the current value of the independent variable       
!                   which the solver has actually reached, i.e. the     
!                   current internal mesh point in t.  On output, TCUR  
!                   will always be at least as far as the argument      
!                   T, but may be farther (if interpolation was done).  
!                                                                       
! TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,         
!                   computed when a request for too much accuracy was   
!                   detected (ISTATE = -3 if detected at the start of   
!                   the problem, ISTATE = -2 otherwise).  If ITOL is    
!                   left unaltered but RTOL and ATOL are uniformly      
!                   scaled up by a factor of TOLSF for the next call,   
!                   then the solver is deemed likely to succeed.        
!                   (The user may also ignore TOLSF and alter the       
!                   tolerance parameters in any other way appropriate.) 
!                                                                       
! TSW     RWORK(15) the value of t at the time of the last method       
!                   switch, if any.                                     
!                                                                       
! NST     IWORK(11) the number of steps taken for the problem so far.   
!                                                                       
! NFE     IWORK(12) the number of f evaluations for the problem so far. 
!                                                                       
! NJE     IWORK(13) the number of Jacobian evaluations (and of matrix   
!                   LU decompositions) for the problem so far.          
!                                                                       
! NQU     IWORK(14) the method order last used (successfully).          
!                                                                       
! NQCUR   IWORK(15) the order to be attempted on the next step.         
!                                                                       
! IMXER   IWORK(16) the index of the component of largest magnitude in  
!                   the weighted local error vector ( E(i)/EWT(i) ),    
!                   on an error return with ISTATE = -4 or -5.          
!                                                                       
! LENRW   IWORK(17) the length of RWORK actually required, assuming     
!                   that the length of RWORK is to be fixed for the     
!                   rest of the problem, and that switching may occur.  
!                   This is defined on normal returns and on an illegal 
!                   input return for insufficient storage.              
!                                                                       
! LENIW   IWORK(18) the length of IWORK actually required, assuming     
!                   that the length of IWORK is to be fixed for the     
!                   rest of the problem, and that switching may occur.  
!                   This is defined on normal returns and on an illegal 
!                   input return for insufficient storage.              
!                                                                       
! MUSED   IWORK(19) the method indicator for the last successful step:  
!                   1 means Adams (nonstiff), 2 means BDF (stiff).      
!                                                                       
! MCUR    IWORK(20) the current method indicator:                       
!                   1 means Adams (nonstiff), 2 means BDF (stiff).      
!                   This is the method to be attempted                  
!                   on the next step.  Thus it differs from MUSED       
!                   only if a method switch has just been made.         
!                                                                       
! The following two arrays are segments of the RWORK array which        
! may also be of interest to the user as optional outputs.              
! For each array, the table below gives its internal name,              
! its base address in RWORK, and its description.                       
!                                                                       
! Name    Base Address      Description                                 
!                                                                       
! YH      21             the Nordsieck history array, of size NYH by    
!                        (NQCUR + 1), where NYH is the initial value    
!                        of NEQ.  For j = 0,1,...,NQCUR, column j+1     
!                        of YH contains HCUR**j/factorial(j) times      
!                        the j-th derivative of the interpolating       
!                        polynomial currently representing the solution,
!                        evaluated at T = TCUR.                         
!                                                                       
! ACOR     LACOR         array of size NEQ used for the accumulated     
!         (from Common   corrections on each step, scaled on output     
!           as noted)    to represent the estimated local error in y    
!                        on the last step.  This is the vector E in     
!                        the description of the error control.  It is   
!                        defined only on a successful return from       
!                        DLSODA.  The base address LACOR is obtained by 
!                        including in the user's program the            
!                        following 2 lines:                             
!                           COMMON /DLS001/ RLS(218), ILS(37)           
!                           LACOR = ILS(22)                             
!                                                                       
!-----------------------------------------------------------------------
! Part 2.  Other Routines Callable.                                     
!                                                                       
! The following are optional calls which the user may make to           
! gain additional capabilities in conjunction with DLSODA.              
! (The routines XSETUN and XSETF are designed to conform to the         
! SLATEC error handling package.)                                       
!                                                                       
!     Form of Call                  Function                            
!   CALL XSETUN(LUN)          set the logical unit number, LUN, for     
!                             output of messages from DLSODA, if        
!                             the default is not desired.               
!                             The default value of LUN is 6.            
!                                                                       
!   CALL XSETF(MFLAG)         set a flag to control the printing of     
!                             messages by DLSODA.                       
!                             MFLAG = 0 means do not print. (Danger:    
!                             This risks losing valuable information.)  
!                             MFLAG = 1 means print (the default).      
!                                                                       
!                             Either of the above calls may be made at  
!                             any time and will take effect immediately.
!                                                                       
!   CALL DSRCMA(RSAV,ISAV,JOB) saves and restores the contents of       
!                             the internal Common blocks used by        
!                             DLSODA (see Part 3 below).                
      type(odepack_common_data), target, intent(inout) :: common_data 
!                             RSAV must be a real array of length 240   
!                             or more, and ISAV must be an integer      
!                             array of length 46 or more.               
!                             JOB=1 means save Common into RSAV/ISAV.   
!                             JOB=2 means restore Common from RSAV/ISAV.
!                                DSRCMA is useful if one is             
!                             interrupting a run and restarting         
!                             later, or alternating between two or      
!                             more problems solved with DLSODA.         
!                                                                       
!   CALL DINTDY(,,,,,)        provide derivatives of y, of various      
!        (see below)          orders, at a specified point t, if        
!                             desired.  It may be called only after     
!                             a successful return from DLSODA.          
!                                                                       
! The detailed instructions for using DINTDY are as follows.            
! The form of the call is:                                              
!                                                                       
!   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)                      
!                                                                       
! The input parameters are:                                             
!                                                                       
! T         = value of independent variable where answers are desired   
!             (normally the same as the T last returned by DLSODA).     
!             For valid results, T must lie between TCUR - HU and TCUR. 
!             (See optional outputs for TCUR and HU.)                   
! K         = integer order of the derivative desired.  K must satisfy  
!             0 .le. K .le. NQCUR, where NQCUR is the current order     
!             (see optional outputs).  The capability corresponding     
!             to K = 0, i.e. computing y(T), is already provided        
!             by DLSODA directly.  Since NQCUR .ge. 1, the first        
!             derivative dy/dt is always available with DINTDY.         
! RWORK(21) = the base address of the history array YH.                 
! NYH       = column length of YH, equal to the initial value of NEQ.   
!                                                                       
! The output parameters are:                                            
!                                                                       
! DKY       = a real array of length NEQ containing the computed value  
!             of the K-th derivative of y(t).                           
! IFLAG     = integer flag, returned as 0 if K and T were legal,        
!             -1 if K was illegal, and -2 if T was illegal.             
!             On an error return, a message is also written.            
!-----------------------------------------------------------------------
! Part 3.  Common Blocks.                                               
!                                                                       
! If DLSODA is to be used in an overlay situation, the user             
! must declare, in the primary overlay, the variables in:               
!   (1) the call sequence to DLSODA, and                                
!   (2) the two internal Common blocks                                  
!         /DLS001/  of length  255  (218 double precision words         
!                      followed by 37 integer words),                   
!         /DLSA01/  of length  31    (22 double precision words         
!                      followed by  9 integer words).                   
!                                                                       
! If DLSODA is used on a system in which the contents of internal       
! Common blocks are not preserved between calls, the user should        
! declare the above Common blocks in the calling program to insure      
! that their contents are preserved.                                    
!                                                                       
! If the solution of a given problem by DLSODA is to be interrupted     
! and then later continued, such as when restarting an interrupted run  
! or alternating between two or more problems, the user should save,    
! following the return from the last DLSODA call prior to the           
! interruption, the contents of the call sequence variables and the     
! internal Common blocks, and later restore these values before the     
! next DLSODA call for that problem.  To save and restore the Common    
! blocks, use Subroutine DSRCMA (see Part 2 above).                     
!                                                                       
!-----------------------------------------------------------------------
! Part 4.  Optionally Replaceable Solver Routines.                      
!                                                                       
! Below is a description of a routine in the DLSODA package which       
! relates to the measurement of errors, and can be                      
! replaced by a user-supplied version, if desired.  However, since such 
! a replacement may have a major impact on performance, it should be    
! done only when absolutely necessary, and only with great caution.     
! (Note: The means by which the package version of a routine is         
! superseded by the user's version may be system-dependent.)            
!                                                                       
! (a) DEWSET.                                                           
! The following subroutine is called just before each internal          
! integration step, and sets the array of error weights, EWT, as        
! described under ITOL/RTOL/ATOL above:                                 
!     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)              
! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODA call sequence,   
! YCUR contains the current dependent variable vector, and              
! EWT is the array of weights set by DEWSET.                            
!                                                                       
! If the user supplies this subroutine, it must return in EWT(i)        
! (i = 1,...,NEQ) a positive quantity suitable for comparing errors     
! in y(i) to.  The EWT array returned by DEWSET is passed to the        
! DMNORM routine, and also used by DLSODA in the computation            
! of the optional output IMXER, and the increments for difference       
! quotient Jacobians.                                                   
!                                                                       
! In the user-supplied version of DEWSET, it may be desirable to use    
! the current values of derivatives of y.  Derivatives up to order NQ   
! are available from the history array YH, described above under        
! optional outputs.  In DEWSET, YH is identical to the YCUR array,      
! extended to NQ + 1 columns with a column length of NYH and scale      
! factors of H**j/factorial(j).  On the first call for the problem,     
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.            
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST       
! can be obtained by including in DEWSET the statements:                
!     DOUBLE PRECISION RLS                                              
!     COMMON /DLS001/ RLS(218),ILS(37)                                  
!     NQ = ILS(33)                                                      
!     NST = ILS(34)                                                     
!     H = RLS(212)                                                      
! Thus, for example, the current value of dy/dt can be obtained as      
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is               
! unnecessary when NST = 0).                                            
!-----------------------------------------------------------------------
!                                                                       
!***REVISION HISTORY  (YYYYMMDD)                                        
! 19811102  DATE WRITTEN                                                
! 19820126  Fixed bug in tests of work space lengths;                   
!           minor corrections in main prologue and comments.            
! 19870330  Major update: corrected comments throughout;                
!           removed TRET from Common; rewrote EWSET with 4 loops;       
!           fixed t test in INTDY; added Cray directives in STODA;      
!           in STODA, fixed DELP init. and logic around PJAC call;      
!           combined routines to save/restore Common;                   
!           passed LEVEL = 0 in error message calls (except run abort). 
! 19970225  Fixed lines setting JSTART = -2 in Subroutine LSODA.        
! 20010425  Major update: convert source lines to upper case;           
!           added *DECK lines; changed from 1 to * in dummy dimensions; 
!           changed names R1MACH/D1MACH to RUMACH/DUMACH;               
!           renamed routines for uniqueness across single/double prec.; 
!           converted intrinsic names to generic form;                  
!           removed ILLIN and NTREP (data loaded) from Common;          
!           removed all 'own' variables from Common;                    
!           changed error messages to quoted strings;                   
!           replaced XERRWV/XERRWD with 1993 revised version;           
!           converted prologues, comments, error messages to mixed case;
!           numerous corrections to prologues and internal comments.    
! 20010507  Converted single precision source to double precision.      
! 20010613  Revised excess accuracy test (to match rest of ODEPACK).    
! 20010808  Fixed bug in DPRJA (matrix in DBNORM call).                 
! 20020502  Corrected declarations in descriptions of user routines.    
! 20031105  Restored 'own' variables to Common blocks, to enable        
!           interrupt/restart feature.                                  
! 20031112  Added SAVE statements for data-loaded constants.            
!                                                                       
!-----------------------------------------------------------------------
! Other routines in the DLSODA package.                                 
!                                                                       
! In addition to Subroutine DLSODA, the DLSODA package includes the     
! following subroutines and function routines:                          
!  DINTDY   computes an interpolated value of the y vector at t = TOUT. 
!  DSTODA   is the core integrator, which does one step of the          
!           integration and the associated error control.               
!  DCFODE   sets all method coefficients and test constants.            
!  DPRJA    computes and preprocesses the Jacobian matrix J = df/dy     
!           and the Newton iteration matrix P = I - h*l0*J.             
!  DSOLSY   manages solution of linear system in chord iteration.       
!  DEWSET   sets the error weight vector EWT before each step.          
!  DMNORM   computes the weighted max-norm of a vector.                 
!  DFNORM   computes the norm of a full matrix consistent with the      
!           weighted max-norm on vectors.                               
!  DBNORM   computes the norm of a band matrix consistent with the      
!           weighted max-norm on vectors.                               
!  DSRCMA   is a user-callable routine to save and restore              
!           the contents of the internal Common blocks.                 
!  DGEFA and DGESL   are routines from LINPACK for solving full         
!           systems of linear algebraic equations.                      
!  DGBFA and DGBSL   are routines from LINPACK for solving banded       
!           linear systems.                                             
!  DUMACH   computes the unit roundoff in a machine-independent manner. 
!  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all 
!           error messages and warnings.  XERRWD is machine-dependent.  
! Note:  DMNORM, DFNORM, DBNORM, DUMACH, IXSAV, and IUMACH are          
! function routines.  All the others are subroutines.                   
!                                                                       
!-----------------------------------------------------------------------
      EXTERNAL DPRJA, DSOLSY 
      DOUBLE PRECISION DUMACH, DMNORM 
      INTEGER, pointer ::                                               &
     &   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(:),            &
     &   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                     &
     &   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,               &
     &   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU         
      INTEGER, pointer ::                                               &
     &   INSUFR, INSUFI, IXPR, IOWNS2(:), JTYP, MUSED, MXORDN, MXORDS   
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LF0,                        &
     &   LENIW, LENRW, LENWM, ML, MORD, MU, MXHNL0, MXSTP0              
      INTEGER LEN1, LEN1C, LEN1N, LEN1S, LEN2, LENIWC, LENRWC 
      DOUBLE PRECISION, pointer :: ROWNS(:),                            &
     &   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND                  
      DOUBLE PRECISION, pointer :: TSW, ROWNS2(:), PDNORM 
      DOUBLE PRECISION ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI, &
     &   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0             
      DIMENSION MORD(2) 
      LOGICAL IHIT 
      CHARACTER*60 MSG 
      SAVE MORD, MXSTP0, MXHNL0 
!-----------------------------------------------------------------------
! The following two internal Common blocks contain                      
! (a) variables which are local to any subroutine but whose values must 
!     be preserved between calls to the routine ("own" variables), and  
! (b) variables which are communicated between subroutines.             
! The block DLS001 is declared in subroutines DLSODA, DINTDY, DSTODA,   
! DPRJA, and DSOLSY.                                                    
! The block DLSA01 is declared in subroutines DLSODA, DSTODA, and DPRJA.
! Groups of variables are replaced by dummy arrays in the Common        
! declarations in routines where those variables are not used.          
!-----------------------------------------------------------------------
!      COMMON /DLS001/ ROWNS(209),                                      
!     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,                
!     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),           
!     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                    
!     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,              
!     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU        
!                                                                       
!      COMMON /DLSA01/ TSW, ROWNS2(20), PDNORM,                         
!     1   INSUFR, INSUFI, IXPR, IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS  
!                                                                       
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/ 
!                                                                       
!     Common block pointers                                             
      type(DLS001_type), pointer :: DLS001 
      type(DLSA01_type), pointer :: DLSA01 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      DLS001 => common_data%DLS001 
                                                                        
      tmp_ptr = c_loc(DLS001%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWNS,[209]) 
      tmp_ptr = c_loc(DLS001%reals(210)) 
      call c_f_pointer(tmp_ptr,CCMAX) 
      tmp_ptr = c_loc(DLS001%reals(211)) 
      call c_f_pointer(tmp_ptr,EL0) 
      tmp_ptr = c_loc(DLS001%reals(212)) 
      call c_f_pointer(tmp_ptr,H) 
      tmp_ptr = c_loc(DLS001%reals(213)) 
      call c_f_pointer(tmp_ptr,HMIN) 
      tmp_ptr = c_loc(DLS001%reals(214)) 
      call c_f_pointer(tmp_ptr,HMXI) 
      tmp_ptr = c_loc(DLS001%reals(215)) 
      call c_f_pointer(tmp_ptr,HU) 
      tmp_ptr = c_loc(DLS001%reals(216)) 
      call c_f_pointer(tmp_ptr,RC) 
      tmp_ptr = c_loc(DLS001%reals(217)) 
      call c_f_pointer(tmp_ptr,TN) 
      tmp_ptr = c_loc(DLS001%reals(218)) 
      call c_f_pointer(tmp_ptr,UROUND) 
                                                                        
      tmp_ptr = c_loc(DLS001%ints(1)) 
      call c_f_pointer(tmp_ptr,INIT) 
      tmp_ptr = c_loc(DLS001%ints(2)) 
      call c_f_pointer(tmp_ptr,MXSTEP) 
      tmp_ptr = c_loc(DLS001%ints(3)) 
      call c_f_pointer(tmp_ptr,MXHNIL) 
      tmp_ptr = c_loc(DLS001%ints(4)) 
      call c_f_pointer(tmp_ptr,NHNIL) 
      tmp_ptr = c_loc(DLS001%ints(5)) 
      call c_f_pointer(tmp_ptr,NSLAST) 
      tmp_ptr = c_loc(DLS001%ints(6)) 
      call c_f_pointer(tmp_ptr,NYH) 
      tmp_ptr = c_loc(DLS001%ints(7)) 
      call c_f_pointer(tmp_ptr,IOWNS,[6]) 
      tmp_ptr = c_loc(DLS001%ints(13)) 
      call c_f_pointer(tmp_ptr,ICF) 
      tmp_ptr = c_loc(DLS001%ints(14)) 
      call c_f_pointer(tmp_ptr,IERPJ) 
      tmp_ptr = c_loc(DLS001%ints(15)) 
      call c_f_pointer(tmp_ptr,IERSL) 
      tmp_ptr = c_loc(DLS001%ints(16)) 
      call c_f_pointer(tmp_ptr,JCUR) 
      tmp_ptr = c_loc(DLS001%ints(17)) 
      call c_f_pointer(tmp_ptr,JSTART) 
      tmp_ptr = c_loc(DLS001%ints(18)) 
      call c_f_pointer(tmp_ptr,KFLAG) 
      tmp_ptr = c_loc(DLS001%ints(19)) 
      call c_f_pointer(tmp_ptr,L) 
      tmp_ptr = c_loc(DLS001%ints(20)) 
      call c_f_pointer(tmp_ptr,LYH) 
      tmp_ptr = c_loc(DLS001%ints(21)) 
      call c_f_pointer(tmp_ptr,LEWT) 
      tmp_ptr = c_loc(DLS001%ints(22)) 
      call c_f_pointer(tmp_ptr,LACOR) 
      tmp_ptr = c_loc(DLS001%ints(23)) 
      call c_f_pointer(tmp_ptr,LSAVF) 
      tmp_ptr = c_loc(DLS001%ints(24)) 
      call c_f_pointer(tmp_ptr,LWM) 
      tmp_ptr = c_loc(DLS001%ints(25)) 
      call c_f_pointer(tmp_ptr,LIWM) 
      tmp_ptr = c_loc(DLS001%ints(26)) 
      call c_f_pointer(tmp_ptr,METH) 
      tmp_ptr = c_loc(DLS001%ints(27)) 
      call c_f_pointer(tmp_ptr,MITER) 
      tmp_ptr = c_loc(DLS001%ints(28)) 
      call c_f_pointer(tmp_ptr,MAXORD) 
      tmp_ptr = c_loc(DLS001%ints(29)) 
      call c_f_pointer(tmp_ptr,MAXCOR) 
      tmp_ptr = c_loc(DLS001%ints(30)) 
      call c_f_pointer(tmp_ptr,MSBP) 
      tmp_ptr = c_loc(DLS001%ints(31)) 
      call c_f_pointer(tmp_ptr,MXNCF) 
      tmp_ptr = c_loc(DLS001%ints(32)) 
      call c_f_pointer(tmp_ptr,N) 
      tmp_ptr = c_loc(DLS001%ints(33)) 
      call c_f_pointer(tmp_ptr,NQ) 
      tmp_ptr = c_loc(DLS001%ints(34)) 
      call c_f_pointer(tmp_ptr,NST) 
      tmp_ptr = c_loc(DLS001%ints(35)) 
      call c_f_pointer(tmp_ptr,NFE) 
      tmp_ptr = c_loc(DLS001%ints(36)) 
      call c_f_pointer(tmp_ptr,NJE) 
      tmp_ptr = c_loc(DLS001%ints(37)) 
      call c_f_pointer(tmp_ptr,NQU) 
                                                                        
      DLSA01 => common_data%DLSA01 
                                                                        
      tmp_ptr = c_loc(DLSA01%reals(1)) 
      call c_f_pointer(tmp_ptr,TSW) 
      tmp_ptr = c_loc(DLSA01%reals(2)) 
      call c_f_pointer(tmp_ptr,ROWNS2,[20]) 
      tmp_ptr = c_loc(DLSA01%reals(22)) 
      call c_f_pointer(tmp_ptr,PDNORM) 
                                                                        
      tmp_ptr = c_loc(DLSA01%ints(1)) 
      call c_f_pointer(tmp_ptr,INSUFR) 
      tmp_ptr = c_loc(DLSA01%ints(2)) 
      call c_f_pointer(tmp_ptr,INSUFI) 
      tmp_ptr = c_loc(DLSA01%ints(3)) 
      call c_f_pointer(tmp_ptr,IXPR) 
      tmp_ptr = c_loc(DLSA01%ints(4)) 
      call c_f_pointer(tmp_ptr,IOWNS2,[2]) 
      tmp_ptr = c_loc(DLSA01%ints(6)) 
      call c_f_pointer(tmp_ptr,JTYP) 
      tmp_ptr = c_loc(DLSA01%ints(7)) 
      call c_f_pointer(tmp_ptr,MUSED) 
      tmp_ptr = c_loc(DLSA01%ints(8)) 
      call c_f_pointer(tmp_ptr,MXORDN) 
      tmp_ptr = c_loc(DLSA01%ints(9)) 
      call c_f_pointer(tmp_ptr,MXORDS) 
!-----------------------------------------------------------------------
! Block A.                                                              
! This code block is executed on every call.                            
! It tests ISTATE and ITASK for legality and branches appropriately.    
! If ISTATE .gt. 1 but the flag INIT shows that initialization has      
! not yet been done, an error return occurs.                            
! If ISTATE = 1 and TOUT = T, return immediately.                       
!-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601 
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602 
      IF (ISTATE .EQ. 1) GO TO 10 
      IF (INIT .EQ. 0) GO TO 603 
      IF (ISTATE .EQ. 2) GO TO 200 
      GO TO 20 
   10 INIT = 0 
      IF (TOUT .EQ. T) RETURN 
!-----------------------------------------------------------------------
! Block B.                                                              
! The next code block is executed for the initial call (ISTATE = 1),    
! or for a continuation call with parameter changes (ISTATE = 3).       
! It contains checking of all inputs and various initializations.       
!                                                                       
! First check legality of the non-optional inputs NEQ, ITOL, IOPT,      
! JT, ML, and MU.                                                       
!-----------------------------------------------------------------------
   20 IF (NEQ(1) .LE. 0) GO TO 604 
      IF (ISTATE .EQ. 1) GO TO 25 
      IF (NEQ(1) .GT. N) GO TO 605 
   25 N = NEQ(1) 
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606 
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607 
      IF (JT .EQ. 3 .OR. JT .LT. 1 .OR. JT .GT. 5) GO TO 608 
      JTYP = JT 
      IF (JT .LE. 2) GO TO 30 
      ML = IWORK(1) 
      MU = IWORK(2) 
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609 
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610 
   30 CONTINUE 
! Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40 
      IXPR = 0 
      MXSTEP = MXSTP0 
      MXHNIL = MXHNL0 
      HMXI = 0.0D0 
      HMIN = 0.0D0 
      IF (ISTATE .NE. 1) GO TO 60 
      H0 = 0.0D0 
      MXORDN = MORD(1) 
      MXORDS = MORD(2) 
      GO TO 60 
   40 IXPR = IWORK(5) 
      IF (IXPR .LT. 0 .OR. IXPR .GT. 1) GO TO 611 
      MXSTEP = IWORK(6) 
      IF (MXSTEP .LT. 0) GO TO 612 
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0 
      MXHNIL = IWORK(7) 
      IF (MXHNIL .LT. 0) GO TO 613 
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0 
      IF (ISTATE .NE. 1) GO TO 50 
      H0 = RWORK(5) 
      MXORDN = IWORK(8) 
      IF (MXORDN .LT. 0) GO TO 628 
      IF (MXORDN .EQ. 0) MXORDN = 100 
      MXORDN = MIN(MXORDN,MORD(1)) 
      MXORDS = IWORK(9) 
      IF (MXORDS .LT. 0) GO TO 629 
      IF (MXORDS .EQ. 0) MXORDS = 100 
      MXORDS = MIN(MXORDS,MORD(2)) 
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 614 
   50 HMAX = RWORK(6) 
      IF (HMAX .LT. 0.0D0) GO TO 615 
      HMXI = 0.0D0 
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX 
      HMIN = RWORK(7) 
      IF (HMIN .LT. 0.0D0) GO TO 616 
!-----------------------------------------------------------------------
! Set work array pointers and check lengths LRW and LIW.                
! If ISTATE = 1, METH is initialized to 1 here to facilitate the        
! checking of work space lengths.                                       
! Pointers to segments of RWORK and IWORK are named by prefixing L to   
! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).  
! Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.    
! If the lengths provided are insufficient for the current method,      
! an error return occurs.  This is treated as illegal input on the      
! first call, but as a problem interruption with ISTATE = -7 on a       
! continuation call.  If the lengths are sufficient for the current     
! method but not for both methods, a warning message is sent.           
!-----------------------------------------------------------------------
   60 IF (ISTATE .EQ. 1) METH = 1 
      IF (ISTATE .EQ. 1) NYH = N 
      LYH = 21 
      LEN1N = 20 + (MXORDN + 1)*NYH 
      LEN1S = 20 + (MXORDS + 1)*NYH 
      LWM = LEN1S + 1 
      IF (JT .LE. 2) LENWM = N*N + 2 
      IF (JT .GE. 4) LENWM = (2*ML + MU + 1)*N + 2 
      LEN1S = LEN1S + LENWM 
      LEN1C = LEN1N 
      IF (METH .EQ. 2) LEN1C = LEN1S 
      LEN1 = MAX(LEN1N,LEN1S) 
      LEN2 = 3*N 
      LENRW = LEN1 + LEN2 
      LENRWC = LEN1C + LEN2 
      IWORK(17) = LENRW 
      LIWM = 1 
      LENIW = 20 + N 
      LENIWC = 20 
      IF (METH .EQ. 2) LENIWC = LENIW 
      IWORK(18) = LENIW 
      IF (ISTATE .EQ. 1 .AND. LRW .LT. LENRWC) GO TO 617 
      IF (ISTATE .EQ. 1 .AND. LIW .LT. LENIWC) GO TO 618 
      IF (ISTATE .EQ. 3 .AND. LRW .LT. LENRWC) GO TO 550 
      IF (ISTATE .EQ. 3 .AND. LIW .LT. LENIWC) GO TO 555 
      LEWT = LEN1 + 1 
      INSUFR = 0 
      IF (LRW .GE. LENRW) GO TO 65 
      INSUFR = 2 
      LEWT = LEN1C + 1 
      if (common_data%iprint == 1) then 
      MSG='DLSODA-  Warning.. RWORK length is sufficient for now, but  ' 
      CALL XERRWD (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG='      may not be later.  Integration will proceed anyway.   ' 
      CALL XERRWD (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG = '      Length needed is LENRW = I1, while LRW = I2.' 
      CALL XERRWD (MSG, 50, 103, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0) 
      endif 
   65 LSAVF = LEWT + N 
      LACOR = LSAVF + N 
      INSUFI = 0 
      IF (LIW .GE. LENIW) GO TO 70 
      INSUFI = 2 
      if (common_data%iprint == 1) then 
      MSG='DLSODA-  Warning.. IWORK length is sufficient for now, but  ' 
      CALL XERRWD (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG='      may not be later.  Integration will proceed anyway.   ' 
      CALL XERRWD (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG = '      Length needed is LENIW = I1, while LIW = I2.' 
      CALL XERRWD (MSG, 50, 104, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0) 
      endif 
   70 CONTINUE 
! Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1) 
      ATOLI = ATOL(1) 
      DO 75 I = 1,N 
        IF (ITOL .GE. 3) RTOLI = RTOL(I) 
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I) 
        IF (RTOLI .LT. 0.0D0) GO TO 619 
        IF (ATOLI .LT. 0.0D0) GO TO 620 
   75   CONTINUE 
      IF (ISTATE .EQ. 1) GO TO 100 
! If ISTATE = 3, set flag to signal parameter changes to DSTODA. -------
      JSTART = -1 
      IF (N .EQ. NYH) GO TO 200 
! NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH 
      I2 = LYH + (MAXORD + 1)*NYH - 1 
      IF (I1 .GT. I2) GO TO 200 
      DO I = I1,I2 
        RWORK(I) = 0.0D0 
      END DO
      GO TO 200 
!-----------------------------------------------------------------------
! Block C.                                                              
! The next block is for the initial call only (ISTATE = 1).             
! It contains all remaining initializations, the initial call to F,     
! and the calculation of the initial step size.                         
! The error weights in EWT are inverted after being loaded.             
!-----------------------------------------------------------------------
  100 UROUND = DUMACH() 
      TN = T 
      TSW = T 
      MAXORD = MXORDN 
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110 
      TCRIT = RWORK(1) 
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GO TO 625 
      IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0)           &
     &   H0 = TCRIT - T                                                 
  110 JSTART = 0 
      NHNIL = 0 
      NST = 0 
      NJE = 0 
      NSLAST = 0 
      HU = 0.0D0 
      NQU = 0 
      MUSED = 0 
      MITER = 0 
      CCMAX = 0.3D0 
      MAXCOR = 3 
      MSBP = 20 
      MXNCF = 10 
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH 
      CALL F (NEQ, T, Y, RWORK(LF0), common_data%ierr) 
      if (common_data%ierr < 0) then; istate = -8; return; endif 
      NFE = 1 
! Load the initial value vector in YH. ---------------------------------
      DO I = 1,N 
        RWORK(I+LYH-1) = Y(I)
      END DO
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1 
      H = 1.0D0 
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT)) 
      DO I = 1,N 
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621 
        RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
      END DO
!-----------------------------------------------------------------------
! The coding below computes the step size, H0, to be attempted on the   
! first step, unless the user has supplied a value for this.            
! First check that TOUT - T differs significantly from zero.            
! A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))          
! if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted    
! so as to be between 100*UROUND and 1.0E-3.                            
! Then the computed value H0 is given by:                               
!                                                                       
!   H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2                
!                                                                       
! where   w0     = MAX ( ABS(T), ABS(TOUT) ),                           
!         F      = the initial value of the vector f(t,y), and          
!         norm() = the weighted vector norm used throughout, given by   
!                  the DMNORM function routine, and weighted by the     
!                  tolerances initially loaded into the EWT array.      
! The sign of H0 is inferred from the initial values of TOUT and T.     
! ABS(H0) is made .le. ABS(TOUT-T) in any case.                         
!-----------------------------------------------------------------------
      IF (H0 .NE. 0.0D0) GO TO 180 
      TDIST = ABS(TOUT - T) 
      W0 = MAX(ABS(T),ABS(TOUT)) 
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622 
      TOL = RTOL(1) 
      IF (ITOL .LE. 2) GO TO 140 
      DO I = 1,N 
        TOL = MAX(TOL,RTOL(I))
      END DO
  140 IF (TOL .GT. 0.0D0) GO TO 160 
      ATOLI = ATOL(1) 
      DO 150 I = 1,N 
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I) 
        AYI = ABS(Y(I)) 
        IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI) 
  150   CONTINUE 
  160 TOL = MAX(TOL,100.0D0*UROUND) 
      TOL = MIN(TOL,0.001D0) 
      SUM = DMNORM (N, RWORK(LF0), RWORK(LEWT)) 
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2 
      H0 = 1.0D0/SQRT(SUM) 
      H0 = MIN(H0,TDIST) 
      H0 = SIGN(H0,TOUT-T) 
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
  180 RH = ABS(H0)*HMXI 
      IF (RH .GT. 1.0D0) H0 = H0/RH 
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0 
      DO I = 1,N 
        RWORK(I+LF0-1) = H0*RWORK(I+LF0-1) 
      END DO
      GO TO 270 
!-----------------------------------------------------------------------
! Block D.                                                              
! The next code block is for continuation calls only (ISTATE = 2 or 3)  
! and is to check stop conditions before taking a step.                 
!-----------------------------------------------------------------------
  200 NSLAST = NST 
      GO TO (210, 250, 220, 230, 240), ITASK 
  210 IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250 
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG, common_data) 
      IF (IFLAG .NE. 0) GO TO 627 
      T = TOUT 
      GO TO 420 
  220 TP = TN - HU*(1.0D0 + 100.0D0*UROUND) 
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623 
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250 
      T = TN 
      GO TO 400 
  230 TCRIT = RWORK(1) 
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624 
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625 
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245 
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG, common_data) 
      IF (IFLAG .NE. 0) GO TO 627 
      T = TOUT 
      GO TO 420 
  240 TCRIT = RWORK(1) 
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624 
  245 HMX = ABS(TN) + ABS(H) 
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX 
      IF (IHIT) T = TCRIT 
      IF (IHIT) GO TO 400 
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND) 
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250 
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND) 
      IF (ISTATE .EQ. 2 .AND. JSTART .GE. 0) JSTART = -2 
!-----------------------------------------------------------------------
! Block E.                                                              
! The next block is normally executed for all calls and contains        
! the call to the one-step core integrator DSTODA.                      
!                                                                       
! This is a looping point for the integration steps.                    
!                                                                       
! First check for too many steps being taken, update EWT (if not at     
! start of problem), check for too much accuracy being requested, and   
! check for H below the roundoff level in T.                            
!-----------------------------------------------------------------------
  250 CONTINUE 
      IF (METH .EQ. MUSED) GO TO 255 
      IF (INSUFR .EQ. 1) GO TO 550 
      IF (INSUFI .EQ. 1) GO TO 555 
  255 IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500 
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT)) 
      DO I = 1,N 
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510 
        RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1) 
      END DO
  270 TOLSF = UROUND*DMNORM (N, RWORK(LYH), RWORK(LEWT)) 
      IF (TOLSF .LE. 1.0D0) GO TO 280 
      TOLSF = TOLSF*2.0D0 
      IF (NST .EQ. 0) GO TO 626 
      GO TO 520 
  280 IF ((TN + H) .NE. TN) GO TO 290 
      NHNIL = NHNIL + 1 
      IF (NHNIL .GT. MXHNIL) GO TO 290 
      if (common_data%iprint == 1) then 
      MSG = 'DLSODA-  Warning..Internal T (=R1) and H (=R2) are' 
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG='      such that in the machine, T + H = T on the next step  ' 
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG = '     (H = step size). Solver will continue anyway.' 
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H) 
      IF (NHNIL .LT. MXHNIL) GO TO 290 
      MSG = 'DLSODA-  Above warning has been issued I1 times.  ' 
      CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG = '     It will not be issued again for this problem.' 
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0) 
      endif 
  290 CONTINUE 
!-----------------------------------------------------------------------
!   CALL DSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPRJA,DSOLSY)
!-----------------------------------------------------------------------
      CALL DSTODA (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),    &
     &   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM),           &
     &   F, JAC, DPRJA, DSOLSY, common_data)                            
      if (common_data%ierr < 0) then; istate = -8; return; endif 
      KGO = 1 - KFLAG 
      GO TO (300, 530, 540), KGO 
!-----------------------------------------------------------------------
! Block F.                                                              
! The following block handles the case of a successful return from the  
! core integrator (KFLAG = 0).                                          
! If a method switch was just made, record TSW, reset MAXORD,           
! set JSTART to -1 to signal DSTODA to complete the switch,             
! and do extra printing of data if IXPR = 1.                            
! Then, in any case, check for stop conditions.                         
!-----------------------------------------------------------------------
  300 INIT = 1 
      IF (METH .EQ. MUSED) GO TO 310 
      TSW = TN 
      MAXORD = MXORDN 
      IF (METH .EQ. 2) MAXORD = MXORDS 
      IF (METH .EQ. 2) RWORK(LWM) = SQRT(UROUND) 
      INSUFR = MIN(INSUFR,1) 
      INSUFI = MIN(INSUFI,1) 
      JSTART = -1 
      IF (IXPR .EQ. 0) GO TO 310 
      IF (METH .EQ. 2) THEN 
      MSG='DLSODA- A switch to the BDF (stiff) method has occurred     ' 
      CALL XERRWD (MSG, 60, 105, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      ENDIF 
      IF (METH .EQ. 1) THEN 
      MSG='DLSODA- A switch to the Adams (nonstiff) method has occurred' 
      CALL XERRWD (MSG, 60, 106, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      ENDIF 
      MSG='     at T = R1,  tentative step size H = R2,  step NST = I1 ' 
      CALL XERRWD (MSG, 60, 107, 0, 1, NST, 0, 2, TN, H) 
  310 GO TO (320, 400, 330, 340, 350), ITASK 
! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
  320 IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250 
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG, common_data) 
      T = TOUT 
      GO TO 420 
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
  330 IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400 
      GO TO 250 
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary. 
  340 IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345 
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG, common_data) 
      T = TOUT 
      GO TO 420 
  345 HMX = ABS(TN) + ABS(H) 
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX 
      IF (IHIT) GO TO 400 
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND) 
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250 
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND) 
      IF (JSTART .GE. 0) JSTART = -2 
      GO TO 250 
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
  350 HMX = ABS(TN) + ABS(H) 
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX 
!-----------------------------------------------------------------------
! Block G.                                                              
! The following block handles all successful returns from DLSODA.       
! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.        
! ISTATE is set to 2, and the optional outputs are loaded into the      
! work arrays before returning.                                         
!-----------------------------------------------------------------------
  400 DO I = 1,N 
        Y(I) = RWORK(I+LYH-1) 
      END DO
      T = TN 
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420 
      IF (IHIT) T = TCRIT 
  420 ISTATE = 2 
      RWORK(11) = HU 
      RWORK(12) = H 
      RWORK(13) = TN 
      RWORK(15) = TSW 
      IWORK(11) = NST 
      IWORK(12) = NFE 
      IWORK(13) = NJE 
      IWORK(14) = NQU 
      IWORK(15) = NQ 
      IWORK(19) = MUSED 
      IWORK(20) = METH 
      RETURN 
!-----------------------------------------------------------------------
! Block H.                                                              
! The following block handles all unsuccessful returns other than       
! those for illegal input.  First the error message routine is called.  
! If there was an error test or convergence test failure, IMXER is set. 
! Then Y is loaded from YH and T is set to TN.                          
! The optional outputs are loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
! 500  MSG = 'DLSODA-  At current T (=R1), MXSTEP (=I1) steps   '       
!      CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      taken on this call before reaching TOUT     '       
!      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)        
  500 common_data%error_message =                                       &
     & 'MXSTEP steps taken on this call before reaching TOUT.'          
      ISTATE = -1 
      GO TO 580 
! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
! 510  EWTI = RWORK(LEWT+I-1)                                           
!      MSG = 'DLSODA-  At T (=R1), EWT(I1) has become R2 .le. 0.'       
!      CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)              
  510 common_data%error_message =                                       &
     & 'EWT(i) became zero for some i during the integration.'          
      ISTATE = -6 
      GO TO 580 
! Too much accuracy requested for machine precision. -------------------
! 520  MSG = 'DLSODA-  At T (=R1), too much accuracy requested  '       
!      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      for precision of machine..  See TOLSF (=R2) '       
!      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)             
  520 common_data%error_message =                                       &
     & 'Too much accuracy requested for machine precision.'             
      RWORK(14) = TOLSF 
      ISTATE = -2 
      GO TO 580 
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
! 530  MSG = 'DLSODA-  At T(=R1) and step size H(=R2), the error'       
!      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      test failed repeatedly or with ABS(H) = HMIN'       
!      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)                 
  530 common_data%error_message =                                       &
     & 'Error test failed repeatedly or with ABS(H) = HMIN.'            
      ISTATE = -4 
      GO TO 560 
! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
! 540  MSG = 'DLSODA-  At T (=R1) and step size H (=R2), the    '       
!      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      corrector convergence failed repeatedly     '       
!      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      or with ABS(H) = HMIN   '                           
!      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)                 
  540 common_data%error_message =                                       &
     & 'Convergence failed repeatedly or with ABS(H) = HMIN.'           
      ISTATE = -5 
      GO TO 560 
! RWORK length too small to proceed. -----------------------------------
! 550  MSG = 'DLSODA-  At current T(=R1), RWORK length too small'       
!      CALL XERRWD (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG='      to proceed.  The integration was otherwise successful.
!      CALL XERRWD (MSG, 60, 206, 0, 0, 0, 0, 1, TN, 0.0D0)             
  550 common_data%error_message =                                       &
     & 'RWORK length too small to proceed.'                             
      ISTATE = -7 
      GO TO 580 
! IWORK length too small to proceed. -----------------------------------
! 555  MSG = 'DLSODA-  At current T(=R1), IWORK length too small'       
!      CALL XERRWD (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG='      to proceed.  The integration was otherwise successful.
!      CALL XERRWD (MSG, 60, 207, 0, 0, 0, 0, 1, TN, 0.0D0)             
  555 common_data%error_message =                                       &
     & 'IWORK length too small to proceed.'                             
      ISTATE = -7 
      GO TO 580 
! Compute IMXER if relevant. -------------------------------------------
  560 BIG = 0.0D0 
      IMXER = 1 
      DO 570 I = 1,N 
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1)) 
        IF (BIG .GE. SIZE) GO TO 570 
        BIG = SIZE 
        IMXER = I 
  570   CONTINUE 
      IWORK(16) = IMXER 
! Set Y vector, T, and optional outputs. -------------------------------
  580 DO I = 1,N 
        Y(I) = RWORK(I+LYH-1)
      END DO
      T = TN 
      RWORK(11) = HU 
      RWORK(12) = H 
      RWORK(13) = TN 
      RWORK(15) = TSW 
      IWORK(11) = NST 
      IWORK(12) = NFE 
      IWORK(13) = NJE 
      IWORK(14) = NQU 
      IWORK(15) = NQ 
      IWORK(19) = MUSED 
      IWORK(20) = METH 
      RETURN 
!-----------------------------------------------------------------------
! Block I.                                                              
! The following block handles all error returns due to illegal input    
! (ISTATE = -3), as detected before calling the core integrator.        
! First the error message routine is called.  If the illegal input      
! is a negative ISTATE, the run is aborted (apparent infinite loop).    
!-----------------------------------------------------------------------
! 601  MSG = 'DLSODA-  ISTATE (=I1) illegal.'                           
!      CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)       
  601 common_data%error_message = 'ISTATE has illegal value.' 
      IF (ISTATE .LT. 0) GO TO 800 
      GO TO 700 
! 602  MSG = 'DLSODA-  ITASK (=I1) illegal. '                           
!      CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)        
  602 common_data%error_message = 'ITASK has illegal value.' 
      GO TO 700 
! 603  MSG = 'DLSODA-  ISTATE .gt. 1 but DLSODA not initialized.'       
!      CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)            
  603 common_data%error_message =                                       &
     & 'ISTATE > 1 but LSODA not initialized.'                          
      GO TO 700 
! 604  MSG = 'DLSODA-  NEQ (=I1) .lt. 1     '                           
!      CALL XERRWD (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)       
  604 common_data%error_message =                                       &
     & 'NEQ < 1.'                                                       
      GO TO 700 
! 605  MSG = 'DLSODA-  ISTATE = 3 and NEQ increased (I1 to I2). '       
!      CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)       
  605 common_data%error_message =                                       &
     & 'ISTATE == 3 and NEQ increased.'                                 
      GO TO 700 
! 606  MSG = 'DLSODA-  ITOL (=I1) illegal.  '                           
!      CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)         
  606 common_data%error_message =                                       &
     & 'ITOL has illegal value.'                                        
      GO TO 700 
! 607  MSG = 'DLSODA-  IOPT (=I1) illegal.  '                           
!      CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)         
  607 common_data%error_message =                                       &
     & 'IOPT has illegal value.'                                        
      GO TO 700 
! 608  MSG = 'DLSODA-  JT (=I1) illegal.    '                           
!      CALL XERRWD (MSG, 30, 8, 0, 1, JT, 0, 0, 0.0D0, 0.0D0)           
  608 common_data%error_message =                                       &
     & 'JT has illegal value.'                                          
      GO TO 700 
! 609  MSG = 'DLSODA-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '       
!      CALL XERRWD (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0D0, 0.0D0)      
  609 common_data%error_message =                                       &
     & 'ML has illegal value.'                                          
      GO TO 700 
! 610  MSG = 'DLSODA-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '       
!      CALL XERRWD (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0D0, 0.0D0)     
  610 common_data%error_message =                                       &
     & 'MU has illegal value.'                                          
      GO TO 700 
! 611  MSG = 'DLSODA-  IXPR (=I1) illegal.  '                           
!      CALL XERRWD (MSG, 30, 11, 0, 1, IXPR, 0, 0, 0.0D0, 0.0D0)        
  611 common_data%error_message =                                       &
     & 'IXPR has illegal value.'                                        
      GO TO 700 
! 612  MSG = 'DLSODA-  MXSTEP (=I1) .lt. 0  '                           
!      CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)      
  612 common_data%error_message =                                       &
     & 'MXSTEP < 0.'                                                    
      GO TO 700 
! 613  MSG = 'DLSODA-  MXHNIL (=I1) .lt. 0  '                           
!      CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)      
  613 common_data%error_message =                                       &
     & 'MXHNIL < 0'                                                     
      GO TO 700 
! 614  MSG = 'DLSODA-  TOUT (=R1) behind T (=R2)      '                 
!      CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)                
!      MSG = '      Integration direction is given by H0 (=R1)  '       
!      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)              
  614 common_data%error_message =                                       &
     & 'TOUT and T imply integration direction incompatible with H0'    
      GO TO 700 
! 615  MSG = 'DLSODA-  HMAX (=R1) .lt. 0.0  '                           
!      CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)            
  615 common_data%error_message =                                       &
     & 'HMAX < 0'                                                       
      GO TO 700 
! 616  MSG = 'DLSODA-  HMIN (=R1) .lt. 0.0  '                           
!      CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)            
  616 common_data%error_message =                                       &
     & 'HMIN < 0.0'                                                     
      GO TO 700 
! 617  MSG='DLSODA-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)
!      CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)     
  617 common_data%error_message =                                       &
     & 'RWORK is not long enough.'                                      
      GO TO 700 
! 618  MSG='DLSODA-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)
!      CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)     
  618 common_data%error_message =                                       &
     & 'IWORK is not long enough.'                                      
      GO TO 700 
! 619  MSG = 'DLSODA-  RTOL(I1) is R1 .lt. 0.0        '                 
!      CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)           
  619 common_data%error_message =                                       &
     & 'RTOL < 0.0'                                                     
      GO TO 700 
! 620  MSG = 'DLSODA-  ATOL(I1) is R1 .lt. 0.0        '                 
!      CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)           
  620 common_data%error_message =                                       &
     & 'ATOL < 0.0'                                                     
      GO TO 700 
  621 EWTI = RWORK(LEWT+I-1) 
!      MSG = 'DLSODA-  EWT(I1) is R1 .le. 0.0         '                 
!      CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)            
      common_data%error_message =                                       &
     & 'EWT(i) for some i is < 0.0'                                     
      GO TO 700 
! 622  MSG='DLSODA-  TOUT(=R1) too close to T(=R2) to start integration.
!      CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)                
  622 common_data%error_message =                                       &
     & 'TOUT is too close to T to start integration.'                   
      GO TO 700 
! 623  MSG='DLSODA-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  
!      CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)           
  623 common_data%error_message =                                       &
     & 'TOUT switches the direction of integration.'                    
      GO TO 700 
! 624  MSG='DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   
!      CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)              
  624 common_data%error_message =                                       &
     & 'ITASK = 4 or 5 and TCRIT behind TCUR.'                          
      GO TO 700 
! 625  MSG='DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   
!      CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)            
  625 common_data%error_message =                                       &
     & 'ITASK = 4 or 5 and TCRIT behind TOUT.'                          
      GO TO 700 
! 626  MSG = 'DLSODA-  At start of problem, too much accuracy   '       
!      CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)           
!      MSG='      requested for precision of machine..  See TOLSF (=R1) 
!      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)           
  626 common_data%error_message =                                       &
     & 'At start of problem, too much accuracy is requested.'           
      RWORK(14) = TOLSF 
      GO TO 700 
! 627  MSG = 'DLSODA-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'       
!      CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)        
  627 common_data%error_message =                                       &
     & 'Trouble in DINTDY.'                                             
      GO TO 700 
! 628  MSG = 'DLSODA-  MXORDN (=I1) .lt. 0  '                           
!      CALL XERRWD (MSG, 30, 28, 0, 1, MXORDN, 0, 0, 0.0D0, 0.0D0)      
  628 common_data%error_message =                                       &
     & 'MXORDN < 0'                                                     
      GO TO 700 
! 629  MSG = 'DLSODA-  MXORDS (=I1) .lt. 0  '                           
!      CALL XERRWD (MSG, 30, 29, 0, 1, MXORDS, 0, 0, 0.0D0, 0.0D0)      
  629 common_data%error_message =                                       &
     & 'MXORDS < 0'                                                     
!                                                                       
  700 ISTATE = -3 
      RETURN 
!                                                                       
! 800  MSG = 'DLSODA-  Run aborted.. apparent infinite loop.    '       
!      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)          
  800 BLOCK
        CHARACTER(LEN=200) :: tmp_msg  ! Use a fixed length instead of LEN=:
        tmp_msg = TRIM(common_data%error_message) // ' Run aborted. apparent infinite loop.'
        common_data%error_message = TRIM(tmp_msg)
      END BLOCK                       
      RETURN 
!----------------------- End of Subroutine DLSODA ----------------------
      END                                           
!DECK DLSODAR                                                           
      SUBROUTINE DLSODAR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,  &
     &            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT,        &
     &            G, NG, JROOT, common_data)                            
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_interface, only: DINTDY, DROOTS, DRCHEK, DSTODA 
      use odepack_common 
      type(odepack_common_data), target, intent(inout) :: common_data 
      EXTERNAL F, JAC, G 
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT,      &
     &   NG, JROOT                                                      
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK 
!-----------------------------------------------------------------------
! This is the 12 November 2003 version of                               
! DLSODAR: Livermore Solver for Ordinary Differential Equations, with   
!          Automatic method switching for stiff and nonstiff problems,  
!          and with Root-finding.                                       
!                                                                       
! This version is in double precision.                                  
!                                                                       
! DLSODAR solves the initial value problem for stiff or nonstiff        
! systems of first order ODEs,                                          
!     dy/dt = f(t,y) ,  or, in component form,                          
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).    
! At the same time, it locates the roots of any of a set of functions   
!     g(i) = g(i,t,y(1),...,y(NEQ))  (i = 1,...,ng).                    
!                                                                       
! This a variant version of the DLSODE package.  It differs from it     
! in two ways:                                                          
! (a) It switches automatically between stiff and nonstiff methods.     
! This means that the user does not have to determine whether the       
! problem is stiff or not, and the solver will automatically choose the 
! appropriate method.  It always starts with the nonstiff method.       
! (b) It finds the root of at least one of a set of constraint          
! functions g(i) of the independent and dependent variables.            
! It finds only those roots for which some g(i), as a function          
! of t, changes sign in the interval of integration.                    
! It then returns the solution at the root, if that occurs              
! sooner than the specified stop condition, and otherwise returns       
! the solution according the specified stop condition.                  
!                                                                       
! Authors:       Alan C. Hindmarsh,                                     
!                Center for Applied Scientific Computing, L-561         
!                Lawrence Livermore National Laboratory                 
!                Livermore, CA 94551                                    
! and                                                                   
!                Linda R. Petzold                                       
!                Univ. of California at Santa Barbara                   
!                Dept. of Computer Science                              
!                Santa Barbara, CA 93106                                
!                                                                       
! References:                                                           
! 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE     
!     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),  
!     North-Holland, Amsterdam, 1983, pp. 55-64.                        
! 2.  Linda R. Petzold, Automatic Selection of Methods for Solving      
!     Stiff and Nonstiff Systems of Ordinary Differential Equations,    
!     Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.                 
! 3.  Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined    
!     Output Points for Solutions of ODEs, Sandia Report SAND80-0180,   
!     February 1980.                                                    
!-----------------------------------------------------------------------
! Summary of Usage.                                                     
!                                                                       
! Communication between the user and the DLSODAR package, for normal    
! situations, is summarized here.  This summary describes only a subset 
! of the full set of options available.  See the full description for   
! details, including alternative treatment of the Jacobian matrix,      
! optional inputs and outputs, nonstandard options, and                 
! instructions for special situations.  See also the example            
! problem (with program and output) following this summary.             
!                                                                       
! A. First provide a subroutine of the form:                            
!               SUBROUTINE F (NEQ, T, Y, YDOT)                          
!               DOUBLE PRECISION T, Y(*), YDOT(*)                       
! which supplies the vector function f by loading YDOT(i) with f(i).    
!                                                                       
! B. Provide a subroutine of the form:                                  
!               SUBROUTINE G (NEQ, T, Y, NG, GOUT)                      
!               DOUBLE PRECISION T, Y(*), GOUT(NG)                      
! which supplies the vector function g by loading GOUT(i) with          
! g(i), the i-th constraint function whose root is sought.              
!                                                                       
! C. Write a main program which calls Subroutine DLSODAR once for       
! each point at which answers are desired.  This should also provide    
! for possible use of logical unit 6 for output of error messages by    
! DLSODAR.  On the first call to DLSODAR, supply arguments as follows:  
! F      = name of subroutine for right-hand side vector f.             
!          This name must be declared External in calling program.      
! NEQ    = number of first order ODEs.                                  
! Y      = array of initial values, of length NEQ.                      
! T      = the initial value of the independent variable.               
! TOUT   = first point where output is desired (.ne. T).                
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.       
! RTOL   = relative tolerance parameter (scalar).                       
! ATOL   = absolute tolerance parameter (scalar or array).              
!          the estimated local error in y(i) will be controlled so as   
!          to be less than                                              
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or        
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.           
!          Thus the local error test passes if, in each component,      
!          either the absolute error is less than ATOL (or ATOL(i)),    
!          or the relative error is less than RTOL.                     
!          Use RTOL = 0.0 for pure absolute error control, and          
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error    
!          control.  Caution: actual (global) errors may exceed these   
!          local tolerances, so choose them conservatively.             
! ITASK  = 1 for normal computation of output values of y at t = TOUT.  
! ISTATE = integer flag (input and output).  Set ISTATE = 1.            
! IOPT   = 0 to indicate no optional inputs used.                       
! RWORK  = real work array of length at least:                          
!             22 + NEQ * MAX(16, NEQ + 9) + 3*NG.                       
!          See also Paragraph F below.                                  
! LRW    = declared length of RWORK (in user's dimension).              
! IWORK  = integer work array of length at least  20 + NEQ.             
! LIW    = declared length of IWORK (in user's dimension).              
! JAC    = name of subroutine for Jacobian matrix.                      
!          Use a dummy name.  See also Paragraph F below.               
! JT     = Jacobian type indicator.  Set JT = 2.                        
!          See also Paragraph F below.                                  
! G      = name of subroutine for constraint functions, whose           
!          roots are desired during the integration.                    
!          This name must be declared External in calling program.      
! NG     = number of constraint functions g(i).  If there are none,     
!          set NG = 0, and pass a dummy name for G.                     
! JROOT  = integer array of length NG for output of root information.   
!          See next paragraph.                                          
! Note that the main program must declare arrays Y, RWORK, IWORK,       
! JROOT, and possibly ATOL.                                             
!                                                                       
! D. The output from the first call (or any call) is:                   
!      Y = array of computed values of y(t) vector.                     
!      T = corresponding value of independent variable.  This is        
!          TOUT if ISTATE = 2, or the root location if ISTATE = 3,      
!          or the farthest point reached if DLSODAR was unsuccessful.   
! ISTATE = 2 or 3  if DLSODAR was successful, negative otherwise.       
!           2 means no root was found, and TOUT was reached as desired. 
!           3 means a root was found prior to reaching TOUT.            
!          -1 means excess work done on this call (perhaps wrong JT).   
!          -2 means excess accuracy requested (tolerances too small).   
!          -3 means illegal input detected (see printed message).       
!          -4 means repeated error test failures (check all inputs).    
!          -5 means repeated convergence failures (perhaps bad Jacobian 
!             supplied or wrong choice of JT or tolerances).            
!          -6 means error weight became zero during problem. (Solution  
!             component i vanished, and ATOL or ATOL(i) = 0.)           
!          -7 means work space insufficient to finish (see messages).   
! JROOT  = array showing roots found if ISTATE = 3 on return.           
!          JROOT(i) = 1 if g(i) has a root at t, or 0 otherwise.        
!                                                                       
! E. To continue the integration after a successful return, proceed     
! as follows:                                                           
!  (a) If ISTATE = 2 on return, reset TOUT and call DLSODAR again.      
!  (b) If ISTATE = 3 on return, reset ISTATE to 2, call DLSODAR again.  
! In either case, no other parameters need be reset.                    
!                                                                       
! F. Note: If and when DLSODAR regards the problem as stiff, and        
! switches methods accordingly, it must make use of the NEQ by NEQ      
! Jacobian matrix, J = df/dy.  For the sake of simplicity, the          
! inputs to DLSODAR recommended in Paragraph C above cause DLSODAR to   
! treat J as a full matrix, and to approximate it internally by         
! difference quotients.  Alternatively, J can be treated as a band      
! matrix (with great potential reduction in the size of the RWORK       
! array).  Also, in either the full or banded case, the user can supply 
! J in closed form, with a routine whose name is passed as the JAC      
! argument.  These alternatives are described in the paragraphs on      
! RWORK, JAC, and JT in the full description of the call sequence below.
!                                                                       
!-----------------------------------------------------------------------
! Example Problem.                                                      
!                                                                       
! The following is a simple example problem, with the coding            
! needed for its solution by DLSODAR.  The problem is from chemical     
! kinetics, and consists of the following three rate equations:         
!     dy1/dt = -.04*y1 + 1.e4*y2*y3                                     
!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2                         
!     dy3/dt = 3.e7*y2**2                                               
! on the interval from t = 0.0 to t = 4.e10, with initial conditions    
! y1 = 1.0, y2 = y3 = 0.  The problem is stiff.                         
! In addition, we want to find the values of t, y1, y2, and y3 at which 
!   (1) y1 reaches the value 1.e-4, and                                 
!   (2) y3 reaches the value 1.e-2.                                     
!                                                                       
! The following coding solves this problem with DLSODAR,                
! printing results at t = .4, 4., ..., 4.e10, and at the computed       
! roots.  It uses ITOL = 2 and ATOL much smaller for y2 than y1 or y3   
! because y2 has much smaller values.                                   
! At the end of the run, statistical quantities of interest are         
! printed (see optional outputs in the full description below).         
!                                                                       
!     EXTERNAL FEX, GEX                                                 
!     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y                    
!     DIMENSION Y(3), ATOL(3), RWORK(76), IWORK(23), JROOT(2)           
!     NEQ = 3                                                           
!     Y(1) = 1.                                                         
!     Y(2) = 0.                                                         
!     Y(3) = 0.                                                         
!     T = 0.                                                            
!     TOUT = .4                                                         
!     ITOL = 2                                                          
!     RTOL = 1.D-4                                                      
!     ATOL(1) = 1.D-6                                                   
!     ATOL(2) = 1.D-10                                                  
!     ATOL(3) = 1.D-6                                                   
!     ITASK = 1                                                         
!     ISTATE = 1                                                        
!     IOPT = 0                                                          
!     LRW = 76                                                          
!     LIW = 23                                                          
!     JT = 2                                                            
!     NG = 2                                                            
!     DO 40 IOUT = 1,12                                                 
! 10    CALL DLSODAR(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,      
!    1     IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT,GEX,NG,JROOT)               
!       WRITE(6,20)T,Y(1),Y(2),Y(3)                                     
! 20    FORMAT(' At t =',D12.4,'   Y =',3D14.6)                         
!       IF (ISTATE .LT. 0) GO TO 80                                     
!       IF (ISTATE .EQ. 2) GO TO 40                                     
!       WRITE(6,30)JROOT(1),JROOT(2)                                    
! 30    FORMAT(5X,' The above line is a root,  JROOT =',2I5)            
!       ISTATE = 2                                                      
!       GO TO 10                                                        
! 40    TOUT = TOUT*10.                                                 
!     WRITE(6,60)IWORK(11),IWORK(12),IWORK(13),IWORK(10),               
!    1   IWORK(19),RWORK(15)                                            
! 60  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4,      
!    1   '  No. g-s =',I4/                                              
!    2   ' Method last used =',I2,'   Last switch was at t =',D12.4)    
!     STOP                                                              
! 80  WRITE(6,90)ISTATE                                                 
! 90  FORMAT(///' Error halt.. ISTATE =',I3)                            
!     STOP                                                              
!     END                                                               
!                                                                       
!     SUBROUTINE FEX (NEQ, T, Y, YDOT)                                  
!     DOUBLE PRECISION T, Y, YDOT                                       
!     DIMENSION Y(3), YDOT(3)                                           
!     YDOT(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)                              
!     YDOT(3) = 3.D7*Y(2)*Y(2)                                          
!     YDOT(2) = -YDOT(1) - YDOT(3)                                      
!     RETURN                                                            
!     END                                                               
!                                                                       
!     SUBROUTINE GEX (NEQ, T, Y, NG, GOUT)                              
!     DOUBLE PRECISION T, Y, GOUT                                       
!     DIMENSION Y(3), GOUT(2)                                           
!     GOUT(1) = Y(1) - 1.D-4                                            
!     GOUT(2) = Y(3) - 1.D-2                                            
!     RETURN                                                            
!     END                                                               
!                                                                       
! The output of this program (on a CDC-7600 in single precision)        
! is as follows:                                                        
!                                                                       
!   At t =  2.6400e-01   y =  9.899653e-01  3.470563e-05  1.000000e-02  
!        The above line is a root,  JROOT =    0    1                   
!   At t =  4.0000e-01   Y =  9.851712e-01  3.386380e-05  1.479493e-02  
!   At t =  4.0000e+00   Y =  9.055333e-01  2.240655e-05  9.444430e-02  
!   At t =  4.0000e+01   Y =  7.158403e-01  9.186334e-06  2.841505e-01  
!   At t =  4.0000e+02   Y =  4.505250e-01  3.222964e-06  5.494717e-01  
!   At t =  4.0000e+03   Y =  1.831975e-01  8.941774e-07  8.168016e-01  
!   At t =  4.0000e+04   Y =  3.898730e-02  1.621940e-07  9.610125e-01  
!   At t =  4.0000e+05   Y =  4.936363e-03  1.984221e-08  9.950636e-01  
!   At t =  4.0000e+06   Y =  5.161831e-04  2.065786e-09  9.994838e-01  
!   At t =  2.0745e+07   Y =  1.000000e-04  4.000395e-10  9.999000e-01  
!        The above line is a root,  JROOT =    1    0                   
!   At t =  4.0000e+07   Y =  5.179817e-05  2.072032e-10  9.999482e-01  
!   At t =  4.0000e+08   Y =  5.283401e-06  2.113371e-11  9.999947e-01  
!   At t =  4.0000e+09   Y =  4.659031e-07  1.863613e-12  9.999995e-01  
!   At t =  4.0000e+10   Y =  1.404280e-08  5.617126e-14  1.000000e+00  
!                                                                       
!   No. steps = 361  No. f-s = 693  No. J-s =  64  No. g-s = 390        
!   Method last used = 2   Last switch was at t =  6.0092e-03           
!                                                                       
!-----------------------------------------------------------------------
! Full Description of User Interface to DLSODAR.                        
!                                                                       
! The user interface to DLSODAR consists of the following parts.        
!                                                                       
! 1.   The call sequence to Subroutine DLSODAR, which is a driver       
!      routine for the solver.  This includes descriptions of both      
!      the call sequence arguments and of user-supplied routines.       
!      Following these descriptions is a description of                 
!      optional inputs available through the call sequence, and then    
!      a description of optional outputs (in the work arrays).          
!                                                                       
! 2.   Descriptions of other routines in the DLSODAR package that may be
!      (optionally) called by the user.  These provide the ability to   
!      alter error message handling, save and restore the internal      
!      Common, and obtain specified derivatives of the solution y(t).   
!                                                                       
! 3.   Descriptions of Common blocks to be declared in overlay          
!      or similar environments, or to be saved when doing an interrupt  
!      of the problem and continued solution later.                     
!                                                                       
! 4.   Description of a subroutine in the DLSODAR package,              
!      which the user may replace with his/her own version, if desired. 
!      this relates to the measurement of errors.                       
!                                                                       
!-----------------------------------------------------------------------
! Part 1.  Call Sequence.                                               
!                                                                       
! The call sequence parameters used for input only are                  
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC,       
!     JT, G, and NG,                                                    
! that used only for output is  JROOT,                                  
! and those used for both input and output are                          
!     Y, T, ISTATE.                                                     
! The work arrays RWORK and IWORK are also used for conditional and     
! optional inputs and optional outputs.  (The term output here refers   
! to the return from Subroutine DLSODAR to the user's calling program.) 
!                                                                       
! The legality of input parameters will be thoroughly checked on the    
! initial call for the problem, but not checked thereafter unless a     
! change in input parameters is flagged by ISTATE = 3 on input.         
!                                                                       
! The descriptions of the call arguments are as follows.                
!                                                                       
! F      = the name of the user-supplied subroutine defining the        
!          ODE system.  The system must be put in the first-order       
!          form dy/dt = f(t,y), where f is a vector-valued function     
!          of the scalar t and the vector y.  Subroutine F is to        
!          compute the function f.  It is to have the form              
!               SUBROUTINE F (NEQ, T, Y, YDOT)                          
!               DOUBLE PRECISION T, Y(*), YDOT(*)                       
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)   
!          is output.  Y and YDOT are arrays of length NEQ.             
!          Subroutine F should not alter Y(1),...,Y(NEQ).               
!          F must be declared External in the calling program.          
!                                                                       
!          Subroutine F may access user-defined quantities in           
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array      
!          (dimensioned in F) and/or Y has length exceeding NEQ(1).     
!          See the descriptions of NEQ and Y below.                     
!                                                                       
!          If quantities computed in the F routine are needed           
!          externally to DLSODAR, an extra call to F should be made     
!          for this purpose, for consistent and accurate results.       
!          If only the derivative dy/dt is needed, use DINTDY instead.  
!                                                                       
! NEQ    = the size of the ODE system (number of first order            
!          ordinary differential equations).  Used only for input.      
!          NEQ may be decreased, but not increased, during the problem. 
!          If NEQ is decreased (with ISTATE = 3 on input), the          
!          remaining components of Y should be left undisturbed, if     
!          these are to be accessed in F and/or JAC.                    
!                                                                       
!          Normally, NEQ is a scalar, and it is generally referred to   
!          as a scalar in this user interface description.  However,    
!          NEQ may be an array, with NEQ(1) set to the system size.     
!          (The DLSODAR package accesses only NEQ(1).)  In either case, 
!          this parameter is passed as the NEQ argument in all calls    
!          to F, JAC, and G.  Hence, if it is an array, locations       
!          NEQ(2),... may be used to store other integer data and pass  
!          it to F, JAC, and G.  Each such subroutine must include      
!          NEQ in a Dimension statement in that case.                   
!                                                                       
! Y      = a real array for the vector of dependent variables, of       
!          length NEQ or more.  Used for both input and output on the   
!          first call (ISTATE = 1), and only for output on other calls. 
!          On the first call, Y must contain the vector of initial      
!          values.  On output, Y contains the computed solution vector, 
!          evaluated at T.  If desired, the Y array may be used         
!          for other purposes between calls to the solver.              
!                                                                       
!          This array is passed as the Y argument in all calls to F,    
!          JAC, and G.  Hence its length may exceed NEQ, and locations  
!          Y(NEQ+1),... may be used to store other real data and        
!          pass it to F, JAC, and G.  (The DLSODAR package accesses only
!          Y(1),...,Y(NEQ).)                                            
!                                                                       
! T      = the independent variable.  On input, T is used only on the   
!          first call, as the initial point of the integration.         
!          On output, after each call, T is the value at which a        
!          computed solution y is evaluated (usually the same as TOUT). 
!          If a root was found, T is the computed location of the       
!          root reached first, on output.                               
!          On an error return, T is the farthest point reached.         
!                                                                       
! TOUT   = the next value of t at which a computed solution is desired. 
!          Used only for input.                                         
!                                                                       
!          When starting the problem (ISTATE = 1), TOUT may be equal    
!          to T for one call, then should .ne. T for the next call.     
!          For the initial T, an input value of TOUT .ne. T is used     
!          in order to determine the direction of the integration       
!          (i.e. the algebraic sign of the step sizes) and the rough    
!          scale of the problem.  Integration in either direction       
!          (forward or backward in t) is permitted.                     
!                                                                       
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after    
!          the first call (i.e. the first call with TOUT .ne. T).       
!          Otherwise, TOUT is required on every call.                   
!                                                                       
!          If ITASK = 1, 3, or 4, the values of TOUT need not be        
!          monotone, but a value of TOUT which backs up is limited      
!          to the current internal T interval, whose endpoints are      
!          TCUR - HU and TCUR (see optional outputs, below, for         
!          TCUR and HU).                                                
!                                                                       
! ITOL   = an indicator for the type of error control.  See             
!          description below under ATOL.  Used only for input.          
!                                                                       
! RTOL   = a relative error tolerance parameter, either a scalar or     
!          an array of length NEQ.  See description below under ATOL.   
!          Input only.                                                  
!                                                                       
! ATOL   = an absolute error tolerance parameter, either a scalar or    
!          an array of length NEQ.  Input only.                         
!                                                                       
!             The input parameters ITOL, RTOL, and ATOL determine       
!          the error control performed by the solver.  The solver will  
!          control the vector E = (E(i)) of estimated local errors      
!          in y, according to an inequality of the form                 
!                      max-norm of ( E(i)/EWT(i) )   .le.   1,          
!          where EWT = (EWT(i)) is a vector of positive error weights.  
!          The values of RTOL and ATOL should all be non-negative.      
!          The following table gives the types (scalar/array) of        
!          RTOL and ATOL, and the corresponding form of EWT(i).         
!                                                                       
!             ITOL    RTOL       ATOL          EWT(i)                   
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL        
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)     
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL     
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)  
!                                                                       
!          When either of these parameters is a scalar, it need not     
!          be dimensioned in the user's calling program.                
!                                                                       
!          If none of the above choices (with ITOL, RTOL, and ATOL      
!          fixed throughout the problem) is suitable, more general      
!          error controls can be obtained by substituting a             
!          user-supplied routine for the setting of EWT.                
!          See Part 4 below.                                            
!                                                                       
!          If global errors are to be estimated by making a repeated    
!          run on the same problem with smaller tolerances, then all    
!          components of RTOL and ATOL (i.e. of EWT) should be scaled   
!          down uniformly.                                              
!                                                                       
! ITASK  = an index specifying the task to be performed.                
!          input only.  ITASK has the following values and meanings.    
!          1  means normal computation of output values of y(t) at      
!             t = TOUT (by overshooting and interpolating).             
!          2  means take one step only and return.                      
!          3  means stop at the first internal mesh point at or         
!             beyond t = TOUT and return.                               
!          4  means normal computation of output values of y(t) at      
!             t = TOUT but without overshooting t = TCRIT.              
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to   
!             or beyond TOUT, but not behind it in the direction of     
!             integration.  This option is useful if the problem        
!             has a singularity at or beyond t = TCRIT.                 
!          5  means take one step, without passing TCRIT, and return.   
!             TCRIT must be input as RWORK(1).                          
!                                                                       
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT        
!          (within roundoff), it will return T = TCRIT (exactly) to     
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT, 
!          in which case answers at t = TOUT are returned first).       
!                                                                       
! ISTATE = an index used for input and output to specify the            
!          the state of the calculation.                                
!                                                                       
!          On input, the values of ISTATE are as follows.               
!          1  means this is the first call for the problem              
!             (initializations will be done).  See note below.          
!          2  means this is not the first call, and the calculation     
!             is to continue normally, with no change in any input      
!             parameters except possibly TOUT and ITASK.                
!             (If ITOL, RTOL, and/or ATOL are changed between calls     
!             with ISTATE = 2, the new values will be used but not      
!             tested for legality.)                                     
!          3  means this is not the first call, and the                 
!             calculation is to continue normally, but with             
!             a change in input parameters other than                   
!             TOUT and ITASK.  Changes are allowed in                   
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,        
!             and any optional inputs except H0, MXORDN, and MXORDS.    
!             (See IWORK description for ML and MU.)                    
!             In addition, immediately following a return with          
!             ISTATE = 3 (root found), NG and G may be changed.         
!             (But changing NG from 0 to .gt. 0 is not allowed.)        
!          Note:  A preliminary call with TOUT = T is not counted       
!          as a first call here, as no initialization or checking of    
!          input is done.  (Such a call is sometimes useful for the     
!          purpose of outputting the initial conditions.)               
!          Thus the first call for which TOUT .ne. T requires           
!          ISTATE = 1 on input.                                         
!                                                                       
!          On output, ISTATE has the following values and meanings.     
!           1  means nothing was done; TOUT = t and ISTATE = 1 on input.
!           2  means the integration was performed successfully, and    
!              no roots were found.                                     
!           3  means the integration was successful, and one or more    
!              roots were found before satisfying the stop condition    
!              specified by ITASK.  See JROOT.                          
!          -1  means an excessive amount of work (more than MXSTEP      
!              steps) was done on this call, before completing the      
!              requested task, but the integration was otherwise        
!              successful as far as T.  (MXSTEP is an optional input    
!              and is normally 500.)  To continue, the user may         
!              simply reset ISTATE to a value .gt. 1 and call again     
!              (the excess work step counter will be reset to 0).       
!              In addition, the user may increase MXSTEP to avoid       
!              this error return (see below on optional inputs).        
!          -2  means too much accuracy was requested for the precision  
!              of the machine being used.  This was detected before     
!              completing the requested task, but the integration       
!              was successful as far as T.  To continue, the tolerance  
!              parameters must be reset, and ISTATE must be set         
!              to 3.  The optional output TOLSF may be used for this    
!              purpose.  (Note: If this condition is detected before    
!              taking any steps, then an illegal input return           
!              (ISTATE = -3) occurs instead.)                           
!          -3  means illegal input was detected, before taking any      
!              integration steps.  See written message for details.     
!              Note:  If the solver detects an infinite loop of calls   
!              to the solver with illegal input, it will cause          
!              the run to stop.                                         
!          -4  means there were repeated error test failures on         
!              one attempted step, before completing the requested      
!              task, but the integration was successful as far as T.    
!              The problem may have a singularity, or the input         
!              may be inappropriate.                                    
!          -5  means there were repeated convergence test failures on   
!              one attempted step, before completing the requested      
!              task, but the integration was successful as far as T.    
!              This may be caused by an inaccurate Jacobian matrix,     
!              if one is being used.                                    
!          -6  means EWT(i) became zero for some i during the           
!              integration.  Pure relative error control (ATOL(i)=0.0)  
!              was requested on a variable which has now vanished.      
!              The integration was successful as far as T.              
!          -7  means the length of RWORK and/or IWORK was too small to  
!              proceed, but the integration was successful as far as T. 
!              This happens when DLSODAR chooses to switch methods      
!              but LRW and/or LIW is too small for the new method.      
!                                                                       
!          Note:  Since the normal output value of ISTATE is 2,         
!          it does not need to be reset for normal continuation.        
!          Also, since a negative input value of ISTATE will be         
!          regarded as illegal, a negative output value requires the    
!          user to change it, and possibly other inputs, before         
!          calling the solver again.                                    
!                                                                       
! IOPT   = an integer flag to specify whether or not any optional       
!          inputs are being used on this call.  Input only.             
!          The optional inputs are listed separately below.             
!          IOPT = 0 means no optional inputs are being used.            
!                   Default values will be used in all cases.           
!          IOPT = 1 means one or more optional inputs are being used.   
!                                                                       
! RWORK  = a real array (double precision) for work space, and (in the  
!          first 20 words) for conditional and optional inputs and      
!          optional outputs.                                            
!          As DLSODAR switches automatically between stiff and nonstiff 
!          methods, the required length of RWORK can change during the  
!          problem.  Thus the RWORK array passed to DLSODAR can either  
!          have a static (fixed) length large enough for both methods,  
!          or have a dynamic (changing) length altered by the calling   
!          program in response to output from DLSODAR.                  
!                                                                       
!                       --- Fixed Length Case ---                       
!          If the RWORK length is to be fixed, it should be at least    
!               max (LRN, LRS),                                         
!          where LRN and LRS are the RWORK lengths required when the    
!          current method is nonstiff or stiff, respectively.           
!                                                                       
!          The separate RWORK length requirements LRN and LRS are       
!          as follows:                                                  
!          If NEQ is constant and the maximum method orders have        
!          their default values, then                                   
!             LRN = 20 + 16*NEQ + 3*NG,                                 
!             LRS = 22 + 9*NEQ + NEQ**2 + 3*NG           (JT = 1 or 2), 
!             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ + 3*NG   (JT = 4 or 5). 
!          Under any other conditions, LRN and LRS are given by:        
!             LRN = 20 + NYH*(MXORDN+1) + 3*NEQ + 3*NG,                 
!             LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT + 3*NG,          
!          where                                                        
!             NYH    = the initial value of NEQ,                        
!             MXORDN = 12, unless a smaller value is given as an        
!                      optional input,                                  
!             MXORDS = 5, unless a smaller value is given as an         
!                      optional input,                                  
!             LMAT   = length of matrix work space:                     
!             LMAT   = NEQ**2 + 2              if JT = 1 or 2,          
!             LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5.          
!                                                                       
!                       --- Dynamic Length Case ---                     
!          If the length of RWORK is to be dynamic, then it should      
!          be at least LRN or LRS, as defined above, depending on the   
!          current method.  Initially, it must be at least LRN (since   
!          DLSODAR starts with the nonstiff method).  On any return     
!          from DLSODAR, the optional output MCUR indicates the current 
!          method.  If MCUR differs from the value it had on the        
!          previous return, or if there has only been one call to       
!          DLSODAR and MCUR is now 2, then DLSODAR has switched         
!          methods during the last call, and the length of RWORK        
!          should be reset (to LRN if MCUR = 1, or to LRS if            
!          MCUR = 2).  (An increase in the RWORK length is required     
!          if DLSODAR returned ISTATE = -7, but not otherwise.)         
!          After resetting the length, call DLSODAR with ISTATE = 3     
!          to signal that change.                                       
!                                                                       
! LRW    = the length of the array RWORK, as declared by the user.      
!          (This will be checked by the solver.)                        
!                                                                       
! IWORK  = an integer array for work space.                             
!          As DLSODAR switches automatically between stiff and nonstiff 
!          methods, the required length of IWORK can change during      
!          problem, between                                             
!             LIS = 20 + NEQ   and   LIN = 20,                          
!          respectively.  Thus the IWORK array passed to DLSODAR can    
!          either have a fixed length of at least 20 + NEQ, or have a   
!          dynamic length of at least LIN or LIS, depending on the      
!          current method.  The comments on dynamic length under        
!          RWORK above apply here.  Initially, this length need         
!          only be at least LIN = 20.                                   
!                                                                       
!          The first few words of IWORK are used for conditional and    
!          optional inputs and optional outputs.                        
!                                                                       
!          The following 2 words in IWORK are conditional inputs:       
!            IWORK(1) = ML     These are the lower and upper            
!            IWORK(2) = MU     half-bandwidths, respectively, of the    
!                       banded Jacobian, excluding the main diagonal.   
!                       The band is defined by the matrix locations     
!                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU    
!                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.        
!                       These are required if JT is 4 or 5, and         
!                       ignored otherwise.  ML and MU may in fact be    
!                       the band parameters for a matrix to which       
!                       df/dy is only approximately equal.              
!                                                                       
! LIW    = the length of the array IWORK, as declared by the user.      
!          (This will be checked by the solver.)                        
!                                                                       
! Note: The base addresses of the work arrays must not be               
! altered between calls to DLSODAR for the same problem.                
! The contents of the work arrays must not be altered                   
! between calls, except possibly for the conditional and                
! optional inputs, and except for the last 3*NEQ words of RWORK.        
! The latter space is used for internal scratch space, and so is        
! available for use by the user outside DLSODAR between calls, if       
! desired (but not for use by F, JAC, or G).                            
!                                                                       
! JAC    = the name of the user-supplied routine to compute the         
!          Jacobian matrix, df/dy, if JT = 1 or 4.  The JAC routine     
!          is optional, but if the problem is expected to be stiff much 
!          of the time, you are encouraged to supply JAC, for the sake  
!          of efficiency.  (Alternatively, set JT = 2 or 5 to have      
!          DLSODAR compute df/dy internally by difference quotients.)   
!          If and when DLSODAR uses df/dy, it treats this NEQ by NEQ    
!          matrix either as full (JT = 1 or 2), or as banded (JT =      
!          4 or 5) with half-bandwidths ML and MU (discussed under      
!          IWORK above).  In either case, if JT = 1 or 4, the JAC       
!          routine must compute df/dy as a function of the scalar t     
!          and the vector y.  It is to have the form                    
!               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)          
!               DOUBLE PRECISION T, Y(*), PD(NROWPD,*)                  
!          where NEQ, T, Y, ML, MU, and NROWPD are input and the array  
!          PD is to be loaded with partial derivatives (elements of     
!          the Jacobian matrix) on output.  PD must be given a first    
!          dimension of NROWPD.  T and Y have the same meaning as in    
!          Subroutine F.                                                
!               In the full matrix case (JT = 1), ML and MU are         
!          ignored, and the Jacobian is to be loaded into PD in         
!          columnwise manner, with df(i)/dy(j) loaded into pd(i,j).     
!               In the band matrix case (JT = 4), the elements          
!          within the band are to be loaded into PD in columnwise       
!          manner, with diagonal lines of df/dy loaded into the rows    
!          of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!          ML and MU are the half-bandwidth parameters (see IWORK).     
!          The locations in PD in the two triangular areas which        
!          correspond to nonexistent matrix elements can be ignored     
!          or loaded arbitrarily, as they are overwritten by DLSODAR.   
!               JAC need not provide df/dy exactly.  A crude            
!          approximation (possibly with a smaller bandwidth) will do.   
!               In either case, PD is preset to zero by the solver,     
!          so that only the nonzero elements need be loaded by JAC.     
!          Each call to JAC is preceded by a call to F with the same    
!          arguments NEQ, T, and Y.  Thus to gain some efficiency,      
!          intermediate quantities shared by both calculations may be   
!          saved in a user Common block by F and not recomputed by JAC, 
!          if desired.  Also, JAC may alter the Y array, if desired.    
!          JAC must be declared External in the calling program.        
!               Subroutine JAC may access user-defined quantities in    
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array      
!          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).   
!          See the descriptions of NEQ and Y above.                     
!                                                                       
! JT     = Jacobian type indicator.  Used only for input.               
!          JT specifies how the Jacobian matrix df/dy will be           
!          treated, if and when DLSODAR requires this matrix.           
!          JT has the following values and meanings:                    
!           1 means a user-supplied full (NEQ by NEQ) Jacobian.         
!           2 means an internally generated (difference quotient) full  
!             Jacobian (using NEQ extra calls to F per df/dy value).    
!           4 means a user-supplied banded Jacobian.                    
!           5 means an internally generated banded Jacobian (using      
!             ML+MU+1 extra calls to F per df/dy evaluation).           
!          If JT = 1 or 4, the user must supply a Subroutine JAC        
!          (the name is arbitrary) as described above under JAC.        
!          If JT = 2 or 5, a dummy argument can be used.                
!                                                                       
! G      = the name of subroutine for constraint functions, whose       
!          roots are desired during the integration.  It is to have     
!          the form                                                     
!               SUBROUTINE G (NEQ, T, Y, NG, GOUT)                      
!               DOUBLE PRECISION T, Y(*), GOUT(NG)                      
!          where NEQ, T, Y, and NG are input, and the array GOUT        
!          is output.  NEQ, T, and Y have the same meaning as in        
!          the F routine, and GOUT is an array of length NG.            
!          For i = 1,...,NG, this routine is to load into GOUT(i)       
!          the value at (T,Y) of the i-th constraint function g(i).     
!          DLSODAR will find roots of the g(i) of odd multiplicity      
!          (i.e. sign changes) as they occur during the integration.    
!          G must be declared External in the calling program.          
!                                                                       
!          Caution:  Because of numerical errors in the functions       
!          g(i) due to roundoff and integration error, DLSODAR may      
!          return false roots, or return the same root at two or more   
!          nearly equal values of t.  If such false roots are           
!          suspected, the user should consider smaller error tolerances 
!          and/or higher precision in the evaluation of the g(i).       
!                                                                       
!          If a root of some g(i) defines the end of the problem,       
!          the input to DLSODAR should nevertheless allow integration   
!          to a point slightly past that root, so that DLSODAR can      
!          locate the root by interpolation.                            
!                                                                       
!          Subroutine G may access user-defined quantities in           
!          NEQ(2),... and Y(NEQ(1)+1),... if NEQ is an array            
!          (dimensioned in G) and/or Y has length exceeding NEQ(1).     
!          See the descriptions of NEQ and Y above.                     
!                                                                       
! NG     = number of constraint functions g(i).  If there are none,     
!          set NG = 0, and pass a dummy name for G.                     
!                                                                       
! JROOT  = integer array of length NG.  Used only for output.           
!          On a return with ISTATE = 3 (one or more roots found),       
!          JROOT(i) = 1 if g(i) has a root at T, or JROOT(i) = 0 if not.
!-----------------------------------------------------------------------
! Optional Inputs.                                                      
!                                                                       
! The following is a list of the optional inputs provided for in the    
! call sequence.  (See also Part 2.)  For each such input variable,     
! this table lists its name as used in this documentation, its          
! location in the call sequence, its meaning, and the default value.    
! The use of any of these inputs requires IOPT = 1, and in that         
! case all of these inputs are examined.  A value of zero for any       
! of these optional inputs will cause the default value to be used.     
! Thus to use a subset of the optional inputs, simply preload           
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and   
! then set those of interest to nonzero values.                         
!                                                                       
! Name    Location      Meaning and Default Value                       
!                                                                       
! H0      RWORK(5)  the step size to be attempted on the first step.    
!                   The default value is determined by the solver.      
!                                                                       
! HMAX    RWORK(6)  the maximum absolute step size allowed.             
!                   The default value is infinite.                      
!                                                                       
! HMIN    RWORK(7)  the minimum absolute step size allowed.             
!                   The default value is 0.  (This lower bound is not   
!                   enforced on the final step before reaching TCRIT    
!                   when ITASK = 4 or 5.)                               
!                                                                       
! IXPR    IWORK(5)  flag to generate extra printing at method switches. 
!                   IXPR = 0 means no extra printing (the default).     
!                   IXPR = 1 means print data on each switch.           
!                   T, H, and NST will be printed on the same logical   
!                   unit as used for error messages.                    
!                                                                       
! MXSTEP  IWORK(6)  maximum number of (internally defined) steps        
!                   allowed during one call to the solver.              
!                   The default value is 500.                           
!                                                                       
! MXHNIL  IWORK(7)  maximum number of messages printed (per problem)    
!                   warning that T + H = T on a step (H = step size).   
!                   This must be positive to result in a non-default    
!                   value.  The default value is 10.                    
!                                                                       
! MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff    
!                   (Adams) method.  The default value is 12.           
!                   If MXORDN exceeds the default value, it will        
!                   be reduced to the default value.                    
!                   MXORDN is held constant during the problem.         
!                                                                       
! MXORDS  IWORK(9)  the maximum order to be allowed for the stiff       
!                   (BDF) method.  The default value is 5.              
!                   If MXORDS exceeds the default value, it will        
!                   be reduced to the default value.                    
!                   MXORDS is held constant during the problem.         
!-----------------------------------------------------------------------
! Optional Outputs.                                                     
!                                                                       
! As optional additional output from DLSODAR, the variables listed      
! below are quantities related to the performance of DLSODAR            
! which are available to the user.  These are communicated by way of    
! the work arrays, but also have internal mnemonic names as shown.      
! Except where stated otherwise, all of these outputs are defined       
! on any successful return from DLSODAR, and on any return with         
! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return           
! (ISTATE = -3), they will be unchanged from their existing values      
! (if any), except possibly for TOLSF, LENRW, and LENIW.                
! On any error return, outputs relevant to the error will be defined,   
! as noted below.                                                       
!                                                                       
! Name    Location      Meaning                                         
!                                                                       
! HU      RWORK(11) the step size in t last used (successfully).        
!                                                                       
! HCUR    RWORK(12) the step size to be attempted on the next step.     
!                                                                       
! TCUR    RWORK(13) the current value of the independent variable       
!                   which the solver has actually reached, i.e. the     
!                   current internal mesh point in t.  On output, TCUR  
!                   will always be at least as far as the argument      
!                   T, but may be farther (if interpolation was done).  
!                                                                       
! TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,         
!                   computed when a request for too much accuracy was   
!                   detected (ISTATE = -3 if detected at the start of   
!                   the problem, ISTATE = -2 otherwise).  If ITOL is    
!                   left unaltered but RTOL and ATOL are uniformly      
!                   scaled up by a factor of TOLSF for the next call,   
!                   then the solver is deemed likely to succeed.        
!                   (The user may also ignore TOLSF and alter the       
!                   tolerance parameters in any other way appropriate.) 
!                                                                       
! TSW     RWORK(15) the value of t at the time of the last method       
!                   switch, if any.                                     
!                                                                       
! NGE     IWORK(10) the number of g evaluations for the problem so far. 
!                                                                       
! NST     IWORK(11) the number of steps taken for the problem so far.   
!                                                                       
! NFE     IWORK(12) the number of f evaluations for the problem so far. 
!                                                                       
! NJE     IWORK(13) the number of Jacobian evaluations (and of matrix   
!                   LU decompositions) for the problem so far.          
!                                                                       
! NQU     IWORK(14) the method order last used (successfully).          
!                                                                       
! NQCUR   IWORK(15) the order to be attempted on the next step.         
!                                                                       
! IMXER   IWORK(16) the index of the component of largest magnitude in  
!                   the weighted local error vector ( E(i)/EWT(i) ),    
!                   on an error return with ISTATE = -4 or -5.          
!                                                                       
! LENRW   IWORK(17) the length of RWORK actually required, assuming     
!                   that the length of RWORK is to be fixed for the     
!                   rest of the problem, and that switching may occur.  
!                   This is defined on normal returns and on an illegal 
!                   input return for insufficient storage.              
!                                                                       
! LENIW   IWORK(18) the length of IWORK actually required, assuming     
!                   that the length of IWORK is to be fixed for the     
!                   rest of the problem, and that switching may occur.  
!                   This is defined on normal returns and on an illegal 
!                   input return for insufficient storage.              
!                                                                       
! MUSED   IWORK(19) the method indicator for the last successful step:  
!                   1 means Adams (nonstiff), 2 means BDF (stiff).      
!                                                                       
! MCUR    IWORK(20) the current method indicator:                       
!                   1 means Adams (nonstiff), 2 means BDF (stiff).      
!                   This is the method to be attempted                  
!                   on the next step.  Thus it differs from MUSED       
!                   only if a method switch has just been made.         
!                                                                       
! The following two arrays are segments of the RWORK array which        
! may also be of interest to the user as optional outputs.              
! For each array, the table below gives its internal name,              
! its base address in RWORK, and its description.                       
!                                                                       
! Name    Base Address      Description                                 
!                                                                       
! YH      21 + 3*NG      the Nordsieck history array, of size NYH by    
!                        (NQCUR + 1), where NYH is the initial value    
!                        of NEQ.  For j = 0,1,...,NQCUR, column j+1     
!                        of YH contains HCUR**j/factorial(j) times      
!                        the j-th derivative of the interpolating       
!                        polynomial currently representing the solution,
!                        evaluated at t = TCUR.                         
!                                                                       
! ACOR     LACOR         array of size NEQ used for the accumulated     
!         (from Common   corrections on each step, scaled on output     
!           as noted)    to represent the estimated local error in y    
!                        on the last step.  This is the vector E in     
!                        the description of the error control.  It is   
!                        defined only on a successful return from       
!                        DLSODAR.  The base address LACOR is obtained by
!                        including in the user's program the            
!                        following 2 lines:                             
!                           COMMON /DLS001/ RLS(218), ILS(37)           
!                           LACOR = ILS(22)                             
!                                                                       
!-----------------------------------------------------------------------
! Part 2.  Other Routines Callable.                                     
!                                                                       
! The following are optional calls which the user may make to           
! gain additional capabilities in conjunction with DLSODAR.             
! (The routines XSETUN and XSETF are designed to conform to the         
! SLATEC error handling package.)                                       
!                                                                       
!     Form of Call                  Function                            
!   CALL XSETUN(LUN)          Set the logical unit number, LUN, for     
!                             output of messages from DLSODAR, if       
!                             the default is not desired.               
!                             The default value of LUN is 6.            
!                                                                       
!   CALL XSETF(MFLAG)         Set a flag to control the printing of     
!                             messages by DLSODAR.                      
!                             MFLAG = 0 means do not print. (Danger:    
!                             This risks losing valuable information.)  
!                             MFLAG = 1 means print (the default).      
!                                                                       
!                             Either of the above calls may be made at  
!                             any time and will take effect immediately.
!                                                                       
!   CALL DSRCAR(RSAV,ISAV,JOB) saves and restores the contents of       
!                             the internal Common blocks used by        
!                             DLSODAR (see Part 3 below).               
!                             RSAV must be a real array of length 245   
!                             or more, and ISAV must be an integer      
!                             array of length 55 or more.               
!                             JOB=1 means save Common into RSAV/ISAV.   
!                             JOB=2 means restore Common from RSAV/ISAV.
!                                DSRCAR is useful if one is             
!                             interrupting a run and restarting         
!                             later, or alternating between two or      
!                             more problems solved with DLSODAR.        
!                                                                       
!   CALL DINTDY(,,,,,)        Provide derivatives of y, of various      
!        (see below)          orders, at a specified point t, if        
!                             desired.  It may be called only after     
!                             a successful return from DLSODAR.         
!                                                                       
! The detailed instructions for using DINTDY are as follows.            
! The form of the call is:                                              
!                                                                       
!   LYH = 21 + 3*NG                                                     
!   CALL DINTDY (T, K, RWORK(LYH), NYH, DKY, IFLAG)                     
!                                                                       
! The input parameters are:                                             
!                                                                       
! T         = value of independent variable where answers are desired   
!             (normally the same as the T last returned by DLSODAR).    
!             For valid results, T must lie between TCUR - HU and TCUR. 
!             (See optional outputs for TCUR and HU.)                   
! K         = integer order of the derivative desired.  K must satisfy  
!             0 .le. K .le. NQCUR, where NQCUR is the current order     
!             (see optional outputs).  The capability corresponding     
!             to K = 0, i.e. computing y(t), is already provided        
!             by DLSODAR directly.  Since NQCUR .ge. 1, the first       
!             derivative dy/dt is always available with DINTDY.         
! LYH       = 21 + 3*NG = base address in RWORK of the history array YH.
! NYH       = column length of YH, equal to the initial value of NEQ.   
!                                                                       
! The output parameters are:                                            
!                                                                       
! DKY       = a real array of length NEQ containing the computed value  
!             of the K-th derivative of y(t).                           
! IFLAG     = integer flag, returned as 0 if K and T were legal,        
!             -1 if K was illegal, and -2 if T was illegal.             
!             On an error return, a message is also written.            
!-----------------------------------------------------------------------
! Part 3.  Common Blocks.                                               
!                                                                       
! If DLSODAR is to be used in an overlay situation, the user            
! must declare, in the primary overlay, the variables in:               
!   (1) the call sequence to DLSODAR, and                               
!   (2) the three internal Common blocks                                
!         /DLS001/  of length  255  (218 double precision words         
!                      followed by 37 integer words),                   
!         /DLSA01/  of length  31    (22 double precision words         
!                      followed by  9 integer words).                   
!         /DLSR01/  of length   7  (3 double precision words            
!                      followed by  4 integer words).                   
!                                                                       
! If DLSODAR is used on a system in which the contents of internal      
! Common blocks are not preserved between calls, the user should        
! declare the above Common blocks in the calling program to insure      
! that their contents are preserved.                                    
!                                                                       
! If the solution of a given problem by DLSODAR is to be interrupted    
! and then later continued, such as when restarting an interrupted run  
! or alternating between two or more problems, the user should save,    
! following the return from the last DLSODAR call prior to the          
! interruption, the contents of the call sequence variables and the     
! internal Common blocks, and later restore these values before the     
! next DLSODAR call for that problem.  To save and restore the Common   
! blocks, use Subroutine DSRCAR (see Part 2 above).                     
!                                                                       
!-----------------------------------------------------------------------
! Part 4.  Optionally Replaceable Solver Routines.                      
!                                                                       
! Below is a description of a routine in the DLSODAR package which      
! relates to the measurement of errors, and can be                      
! replaced by a user-supplied version, if desired.  However, since such 
! a replacement may have a major impact on performance, it should be    
! done only when absolutely necessary, and only with great caution.     
! (Note: The means by which the package version of a routine is         
! superseded by the user's version may be system-dependent.)            
!                                                                       
! (a) DEWSET.                                                           
! The following subroutine is called just before each internal          
! integration step, and sets the array of error weights, EWT, as        
! described under ITOL/RTOL/ATOL above:                                 
!     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)              
! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODAR call sequence,  
! YCUR contains the current dependent variable vector, and              
! EWT is the array of weights set by DEWSET.                            
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW), &
     &   JROOT(NG)                                                      
!                                                                       
! If the user supplies this subroutine, it must return in EWT(i)        
! (i = 1,...,NEQ) a positive quantity suitable for comparing errors     
! in y(i) to.  The EWT array returned by DEWSET is passed to the        
! DMNORM routine, and also used by DLSODAR in the computation           
! of the optional output IMXER, and the increments for difference       
! quotient Jacobians.                                                   
!                                                                       
! In the user-supplied version of DEWSET, it may be desirable to use    
! the current values of derivatives of y.  Derivatives up to order NQ   
! are available from the history array YH, described above under        
! optional outputs.  In DEWSET, YH is identical to the YCUR array,      
! extended to NQ + 1 columns with a column length of NYH and scale      
! factors of H**j/factorial(j).  On the first call for the problem,     
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.            
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST       
! can be obtained by including in DEWSET the statements:                
!     DOUBLE PRECISION RLS                                              
!     COMMON /DLS001/ RLS(218),ILS(37)                                  
!     NQ = ILS(33)                                                      
!     NST = ILS(34)                                                     
!     H = RLS(212)                                                      
! Thus, for example, the current value of dy/dt can be obtained as      
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is               
! unnecessary when NST = 0).                                            
!-----------------------------------------------------------------------
!                                                                       
!***REVISION HISTORY  (YYYYMMDD)                                        
! 19811102  DATE WRITTEN                                                
! 19820126  Fixed bug in tests of work space lengths;                   
!           minor corrections in main prologue and comments.            
! 19820507  Fixed bug in RCHEK in setting HMING.                        
! 19870330  Major update: corrected comments throughout;                
!           removed TRET from Common; rewrote EWSET with 4 loops;       
!           fixed t test in INTDY; added Cray directives in STODA;      
!           in STODA, fixed DELP init. and logic around PJAC call;      
!           combined routines to save/restore Common;                   
!           passed LEVEL = 0 in error message calls (except run abort). 
! 19970225  Fixed lines setting JSTART = -2 in Subroutine LSODAR.       
! 20010425  Major update: convert source lines to upper case;           
!           added *DECK lines; changed from 1 to * in dummy dimensions; 
!           changed names R1MACH/D1MACH to RUMACH/DUMACH;               
!           renamed routines for uniqueness across single/double prec.; 
!           converted intrinsic names to generic form;                  
!           removed ILLIN and NTREP (data loaded) from Common;          
!           removed all 'own' variables from Common;                    
!           changed error messages to quoted strings;                   
!           replaced XERRWV/XERRWD with 1993 revised version;           
!           converted prologues, comments, error messages to mixed case;
!           numerous corrections to prologues and internal comments.    
! 20010507  Converted single precision source to double precision.      
! 20010613  Revised excess accuracy test (to match rest of ODEPACK).    
! 20010808  Fixed bug in DPRJA (matrix in DBNORM call).                 
! 20020502  Corrected declarations in descriptions of user routines.    
! 20031105  Restored 'own' variables to Common blocks, to enable        
!           interrupt/restart feature.                                  
! 20031112  Added SAVE statements for data-loaded constants.            
!                                                                       
!-----------------------------------------------------------------------
! Other routines in the DLSODAR package.                                
!                                                                       
! In addition to Subroutine DLSODAR, the DLSODAR package includes the   
! following subroutines and function routines:                          
!  DRCHEK   does preliminary checking for roots, and serves as an       
!           interface between Subroutine DLSODAR and Subroutine DROOTS. 
!  DROOTS   finds the leftmost root of a set of functions.              
!  DINTDY   computes an interpolated value of the y vector at t = TOUT. 
!  DSTODA   is the core integrator, which does one step of the          
!           integration and the associated error control.               
!  DCFODE   sets all method coefficients and test constants.            
!  DPRJA    computes and preprocesses the Jacobian matrix J = df/dy     
!           and the Newton iteration matrix P = I - h*l0*J.             
!  DSOLSY   manages solution of linear system in chord iteration.       
!  DEWSET   sets the error weight vector EWT before each step.          
!  DMNORM   computes the weighted max-norm of a vector.                 
!  DFNORM   computes the norm of a full matrix consistent with the      
!           weighted max-norm on vectors.                               
!  DBNORM   computes the norm of a band matrix consistent with the      
!           weighted max-norm on vectors.                               
!  DSRCAR   is a user-callable routine to save and restore              
!           the contents of the internal Common blocks.                 
!  DGEFA and DGESL   are routines from LINPACK for solving full         
!           systems of linear algebraic equations.                      
!  DGBFA and DGBSL   are routines from LINPACK for solving banded       
!           linear systems.                                             
!  DCOPY    is one of the basic linear algebra modules (BLAS).          
!  DUMACH   computes the unit roundoff in a machine-independent manner. 
!  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all 
!           error messages and warnings.  XERRWD is machine-dependent.  
! Note:  DMNORM, DFNORM, DBNORM, DUMACH, IXSAV, and IUMACH are          
! function routines.  All the others are subroutines.                   
!                                                                       
!-----------------------------------------------------------------------
      EXTERNAL DPRJA, DSOLSY 
      DOUBLE PRECISION DUMACH, DMNORM 
      INTEGER, pointer ::                                               &
     &   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(:),            &
     &   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                     &
     &   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,               &
     &   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU         
      INTEGER, pointer ::                                               &
     &   INSUFR, INSUFI, IXPR, IOWNS2(:), JTYP, MUSED, MXORDN, MXORDS   
      INTEGER, pointer ::                                               &
     &   LG0, LG1, LGX, IOWNR3(:), IRFND, ITASKC, NGC, NGE              
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LENIW,                      &
     &   LENRW, LENWM, LF0, ML, MORD, MU, MXHNL0, MXSTP0                
      INTEGER LEN1, LEN1C, LEN1N, LEN1S, LEN2, LENIWC, LENRWC 
      INTEGER IRFP, IRT, LENYH, LYHNEW 
      DOUBLE PRECISION, pointer :: ROWNS(:),                            &
     &   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND                  
      DOUBLE PRECISION, pointer :: TSW, ROWNS2(:), PDNORM 
      DOUBLE PRECISION, pointer :: ROWNR3(:), T0, TLAST, TOUTC 
      DOUBLE PRECISION ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI, &
     &   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0             
      DIMENSION MORD(2) 
      LOGICAL IHIT 
      CHARACTER*60 MSG 
      SAVE MORD, MXSTP0, MXHNL0 
!-----------------------------------------------------------------------
! The following three internal Common blocks contain                    
! (a) variables which are local to any subroutine but whose values must 
!     be preserved between calls to the routine ("own" variables), and  
! (b) variables which are communicated between subroutines.             
! The block DLS001 is declared in subroutines DLSODAR, DINTDY, DSTODA,  
! DPRJA, and DSOLSY.                                                    
! The block DLSA01 is declared in subroutines DLSODAR, DSTODA, DPRJA.   
! The block DLSR01 is declared in subroutines DLSODAR, DRCHEK, DROOTS.  
! Groups of variables are replaced by dummy arrays in the Common        
! declarations in routines where those variables are not used.          
!-----------------------------------------------------------------------
!      COMMON /DLS001/ ROWNS(209),                                      
!     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,                
!     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),           
!     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                    
!     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,              
!     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU        
!                                                                       
!      COMMON /DLSA01/ TSW, ROWNS2(20), PDNORM,                         
!     1   INSUFR, INSUFI, IXPR, IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS  
!                                                                       
!      COMMON /DLSR01/ ROWNR3(2), T0, TLAST, TOUTC,                     
!     1   LG0, LG1, LGX, IOWNR3(2), IRFND, ITASKC, NGC, NGE             
!                                                                       
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/ 
!                                                                       
!     Common block pointers                                             
      type(DLS001_type), pointer :: DLS001 
      type(DLSA01_type), pointer :: DLSA01 
      type(DLSR01_type), pointer :: DLSR01 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      DLS001 => common_data%DLS001 
                                                                        
      tmp_ptr = c_loc(DLS001%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWNS,[209]) 
      tmp_ptr = c_loc(DLS001%reals(210)) 
      call c_f_pointer(tmp_ptr,CCMAX) 
      tmp_ptr = c_loc(DLS001%reals(211)) 
      call c_f_pointer(tmp_ptr,EL0) 
      tmp_ptr = c_loc(DLS001%reals(212)) 
      call c_f_pointer(tmp_ptr,H) 
      tmp_ptr = c_loc(DLS001%reals(213)) 
      call c_f_pointer(tmp_ptr,HMIN) 
      tmp_ptr = c_loc(DLS001%reals(214)) 
      call c_f_pointer(tmp_ptr,HMXI) 
      tmp_ptr = c_loc(DLS001%reals(215)) 
      call c_f_pointer(tmp_ptr,HU) 
      tmp_ptr = c_loc(DLS001%reals(216)) 
      call c_f_pointer(tmp_ptr,RC) 
      tmp_ptr = c_loc(DLS001%reals(217)) 
      call c_f_pointer(tmp_ptr,TN) 
      tmp_ptr = c_loc(DLS001%reals(218)) 
      call c_f_pointer(tmp_ptr,UROUND) 
                                                                        
      tmp_ptr = c_loc(DLS001%ints(1)) 
      call c_f_pointer(tmp_ptr,INIT) 
      tmp_ptr = c_loc(DLS001%ints(2)) 
      call c_f_pointer(tmp_ptr,MXSTEP) 
      tmp_ptr = c_loc(DLS001%ints(3)) 
      call c_f_pointer(tmp_ptr,MXHNIL) 
      tmp_ptr = c_loc(DLS001%ints(4)) 
      call c_f_pointer(tmp_ptr,NHNIL) 
      tmp_ptr = c_loc(DLS001%ints(5)) 
      call c_f_pointer(tmp_ptr,NSLAST) 
      tmp_ptr = c_loc(DLS001%ints(6)) 
      call c_f_pointer(tmp_ptr,NYH) 
      tmp_ptr = c_loc(DLS001%ints(7)) 
      call c_f_pointer(tmp_ptr,IOWNS,[6]) 
      tmp_ptr = c_loc(DLS001%ints(13)) 
      call c_f_pointer(tmp_ptr,ICF) 
      tmp_ptr = c_loc(DLS001%ints(14)) 
      call c_f_pointer(tmp_ptr,IERPJ) 
      tmp_ptr = c_loc(DLS001%ints(15)) 
      call c_f_pointer(tmp_ptr,IERSL) 
      tmp_ptr = c_loc(DLS001%ints(16)) 
      call c_f_pointer(tmp_ptr,JCUR) 
      tmp_ptr = c_loc(DLS001%ints(17)) 
      call c_f_pointer(tmp_ptr,JSTART) 
      tmp_ptr = c_loc(DLS001%ints(18)) 
      call c_f_pointer(tmp_ptr,KFLAG) 
      tmp_ptr = c_loc(DLS001%ints(19)) 
      call c_f_pointer(tmp_ptr,L) 
      tmp_ptr = c_loc(DLS001%ints(20)) 
      call c_f_pointer(tmp_ptr,LYH) 
      tmp_ptr = c_loc(DLS001%ints(21)) 
      call c_f_pointer(tmp_ptr,LEWT) 
      tmp_ptr = c_loc(DLS001%ints(22)) 
      call c_f_pointer(tmp_ptr,LACOR) 
      tmp_ptr = c_loc(DLS001%ints(23)) 
      call c_f_pointer(tmp_ptr,LSAVF) 
      tmp_ptr = c_loc(DLS001%ints(24)) 
      call c_f_pointer(tmp_ptr,LWM) 
      tmp_ptr = c_loc(DLS001%ints(25)) 
      call c_f_pointer(tmp_ptr,LIWM) 
      tmp_ptr = c_loc(DLS001%ints(26)) 
      call c_f_pointer(tmp_ptr,METH) 
      tmp_ptr = c_loc(DLS001%ints(27)) 
      call c_f_pointer(tmp_ptr,MITER) 
      tmp_ptr = c_loc(DLS001%ints(28)) 
      call c_f_pointer(tmp_ptr,MAXORD) 
      tmp_ptr = c_loc(DLS001%ints(29)) 
      call c_f_pointer(tmp_ptr,MAXCOR) 
      tmp_ptr = c_loc(DLS001%ints(30)) 
      call c_f_pointer(tmp_ptr,MSBP) 
      tmp_ptr = c_loc(DLS001%ints(31)) 
      call c_f_pointer(tmp_ptr,MXNCF) 
      tmp_ptr = c_loc(DLS001%ints(32)) 
      call c_f_pointer(tmp_ptr,N) 
      tmp_ptr = c_loc(DLS001%ints(33)) 
      call c_f_pointer(tmp_ptr,NQ) 
      tmp_ptr = c_loc(DLS001%ints(34)) 
      call c_f_pointer(tmp_ptr,NST) 
      tmp_ptr = c_loc(DLS001%ints(35)) 
      call c_f_pointer(tmp_ptr,NFE) 
      tmp_ptr = c_loc(DLS001%ints(36)) 
      call c_f_pointer(tmp_ptr,NJE) 
      tmp_ptr = c_loc(DLS001%ints(37)) 
      call c_f_pointer(tmp_ptr,NQU) 
                                                                        
      DLSA01 => common_data%DLSA01 
                                                                        
      tmp_ptr = c_loc(DLSA01%reals(1)) 
      call c_f_pointer(tmp_ptr,TSW) 
      tmp_ptr = c_loc(DLSA01%reals(2)) 
      call c_f_pointer(tmp_ptr,ROWNS2,[20]) 
      tmp_ptr = c_loc(DLSA01%reals(22)) 
      call c_f_pointer(tmp_ptr,PDNORM) 
                                                                        
      tmp_ptr = c_loc(DLSA01%ints(1)) 
      call c_f_pointer(tmp_ptr,INSUFR) 
      tmp_ptr = c_loc(DLSA01%ints(2)) 
      call c_f_pointer(tmp_ptr,INSUFI) 
      tmp_ptr = c_loc(DLSA01%ints(3)) 
      call c_f_pointer(tmp_ptr,IXPR) 
      tmp_ptr = c_loc(DLSA01%ints(4)) 
      call c_f_pointer(tmp_ptr,IOWNS2,[2]) 
      tmp_ptr = c_loc(DLSA01%ints(6)) 
      call c_f_pointer(tmp_ptr,JTYP) 
      tmp_ptr = c_loc(DLSA01%ints(7)) 
      call c_f_pointer(tmp_ptr,MUSED) 
      tmp_ptr = c_loc(DLSA01%ints(8)) 
      call c_f_pointer(tmp_ptr,MXORDN) 
      tmp_ptr = c_loc(DLSA01%ints(9)) 
      call c_f_pointer(tmp_ptr,MXORDS) 
                                                                        
      DLSR01 => common_data%DLSR01 
                                                                        
      tmp_ptr = c_loc(DLSR01%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWNR3,[2]) 
      tmp_ptr = c_loc(DLSR01%reals(3)) 
      call c_f_pointer(tmp_ptr,T0) 
      tmp_ptr = c_loc(DLSR01%reals(4)) 
      call c_f_pointer(tmp_ptr,TLAST) 
      tmp_ptr = c_loc(DLSR01%reals(5)) 
      call c_f_pointer(tmp_ptr,TOUTC) 
                                                                        
      tmp_ptr = c_loc(DLSR01%ints(1)) 
      call c_f_pointer(tmp_ptr,LG0) 
      tmp_ptr = c_loc(DLSR01%ints(2)) 
      call c_f_pointer(tmp_ptr,LG1) 
      tmp_ptr = c_loc(DLSR01%ints(3)) 
      call c_f_pointer(tmp_ptr,LGX) 
      tmp_ptr = c_loc(DLSR01%ints(4)) 
      call c_f_pointer(tmp_ptr,IOWNR3,[2]) 
      tmp_ptr = c_loc(DLSR01%ints(6)) 
      call c_f_pointer(tmp_ptr,IRFND) 
      tmp_ptr = c_loc(DLSR01%ints(7)) 
      call c_f_pointer(tmp_ptr,ITASKC) 
      tmp_ptr = c_loc(DLSR01%ints(8)) 
      call c_f_pointer(tmp_ptr,NGC) 
      tmp_ptr = c_loc(DLSR01%ints(9)) 
      call c_f_pointer(tmp_ptr,NGE) 
!-----------------------------------------------------------------------
! Block A.                                                              
! This code block is executed on every call.                            
! It tests ISTATE and ITASK for legality and branches appropriately.    
! If ISTATE .gt. 1 but the flag INIT shows that initialization has      
! not yet been done, an error return occurs.                            
! If ISTATE = 1 and TOUT = T, return immediately.                       
!-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601 
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602 
      ITASKC = ITASK 
      IF (ISTATE .EQ. 1) GO TO 10 
      IF (INIT .EQ. 0) GO TO 603 
      IF (ISTATE .EQ. 2) GO TO 200 
      GO TO 20 
   10 INIT = 0 
      IF (TOUT .EQ. T) RETURN 
!-----------------------------------------------------------------------
! Block B.                                                              
! The next code block is executed for the initial call (ISTATE = 1),    
! or for a continuation call with parameter changes (ISTATE = 3).       
! It contains checking of all inputs and various initializations.       
!                                                                       
! First check legality of the non-optional inputs NEQ, ITOL, IOPT,      
! JT, ML, MU, and NG.                                                   
!-----------------------------------------------------------------------
   20 IF (NEQ(1) .LE. 0) GO TO 604 
      IF (ISTATE .EQ. 1) GO TO 25 
      IF (NEQ(1) .GT. N) GO TO 605 
   25 N = NEQ(1) 
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606 
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607 
      IF (JT .EQ. 3 .OR. JT .LT. 1 .OR. JT .GT. 5) GO TO 608 
      JTYP = JT 
      IF (JT .LE. 2) GO TO 30 
      ML = IWORK(1) 
      MU = IWORK(2) 
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609 
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610 
   30 CONTINUE 
      IF (NG .LT. 0) GO TO 630 
      IF (ISTATE .EQ. 1) GO TO 35 
      IF (IRFND .EQ. 0 .AND. NG .NE. NGC) GO TO 631 
   35 NGC = NG 
! Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40 
      IXPR = 0 
      MXSTEP = MXSTP0 
      MXHNIL = MXHNL0 
      HMXI = 0.0D0 
      HMIN = 0.0D0 
      IF (ISTATE .NE. 1) GO TO 60 
      H0 = 0.0D0 
      MXORDN = MORD(1) 
      MXORDS = MORD(2) 
      GO TO 60 
   40 IXPR = IWORK(5) 
      IF (IXPR .LT. 0 .OR. IXPR .GT. 1) GO TO 611 
      MXSTEP = IWORK(6) 
      IF (MXSTEP .LT. 0) GO TO 612 
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0 
      MXHNIL = IWORK(7) 
      IF (MXHNIL .LT. 0) GO TO 613 
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0 
      IF (ISTATE .NE. 1) GO TO 50 
      H0 = RWORK(5) 
      MXORDN = IWORK(8) 
      IF (MXORDN .LT. 0) GO TO 628 
      IF (MXORDN .EQ. 0) MXORDN = 100 
      MXORDN = MIN(MXORDN,MORD(1)) 
      MXORDS = IWORK(9) 
      IF (MXORDS .LT. 0) GO TO 629 
      IF (MXORDS .EQ. 0) MXORDS = 100 
      MXORDS = MIN(MXORDS,MORD(2)) 
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 614 
   50 HMAX = RWORK(6) 
      IF (HMAX .LT. 0.0D0) GO TO 615 
      HMXI = 0.0D0 
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX 
      HMIN = RWORK(7) 
      IF (HMIN .LT. 0.0D0) GO TO 616 
!-----------------------------------------------------------------------
! Set work array pointers and check lengths LRW and LIW.                
! If ISTATE = 1, METH is initialized to 1 here to facilitate the        
! checking of work space lengths.                                       
! Pointers to segments of RWORK and IWORK are named by prefixing L to   
! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).  
! Segments of RWORK (in order) are denoted  G0, G1, GX, YH, WM,         
! EWT, SAVF, ACOR.                                                      
! If the lengths provided are insufficient for the current method,      
! an error return occurs.  This is treated as illegal input on the      
! first call, but as a problem interruption with ISTATE = -7 on a       
! continuation call.  If the lengths are sufficient for the current     
! method but not for both methods, a warning message is sent.           
!-----------------------------------------------------------------------
   60 IF (ISTATE .EQ. 1) METH = 1 
      IF (ISTATE .EQ. 1) NYH = N 
      LG0 = 21 
      LG1 = LG0 + NG 
      LGX = LG1 + NG 
      LYHNEW = LGX + NG 
      IF (ISTATE .EQ. 1) LYH = LYHNEW 
      IF (LYHNEW .EQ. LYH) GO TO 62 
! If ISTATE = 3 and NG was changed, shift YH to its new location. ------
      LENYH = L*NYH 
      IF (LRW .LT. LYHNEW-1+LENYH) GO TO 62 
      I1 = 1 
      IF (LYHNEW .GT. LYH) I1 = -1 
      CALL DCOPY (LENYH, RWORK(LYH), I1, RWORK(LYHNEW), I1) 
      LYH = LYHNEW 
   62 CONTINUE 
      LEN1N = LYHNEW - 1 + (MXORDN + 1)*NYH 
      LEN1S = LYHNEW - 1 + (MXORDS + 1)*NYH 
      LWM = LEN1S + 1 
      IF (JT .LE. 2) LENWM = N*N + 2 
      IF (JT .GE. 4) LENWM = (2*ML + MU + 1)*N + 2 
      LEN1S = LEN1S + LENWM 
      LEN1C = LEN1N 
      IF (METH .EQ. 2) LEN1C = LEN1S 
      LEN1 = MAX(LEN1N,LEN1S) 
      LEN2 = 3*N 
      LENRW = LEN1 + LEN2 
      LENRWC = LEN1C + LEN2 
      IWORK(17) = LENRW 
      LIWM = 1 
      LENIW = 20 + N 
      LENIWC = 20 
      IF (METH .EQ. 2) LENIWC = LENIW 
      IWORK(18) = LENIW 
      IF (ISTATE .EQ. 1 .AND. LRW .LT. LENRWC) GO TO 617 
      IF (ISTATE .EQ. 1 .AND. LIW .LT. LENIWC) GO TO 618 
      IF (ISTATE .EQ. 3 .AND. LRW .LT. LENRWC) GO TO 550 
      IF (ISTATE .EQ. 3 .AND. LIW .LT. LENIWC) GO TO 555 
      LEWT = LEN1 + 1 
      INSUFR = 0 
      IF (LRW .GE. LENRW) GO TO 65 
      INSUFR = 2 
      LEWT = LEN1C + 1 
      if (common_data%iprint == 1) then 
      MSG='DLSODAR-  Warning.. RWORK length is sufficient for now, but ' 
      CALL XERRWD (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG='      may not be later.  Integration will proceed anyway.   ' 
      CALL XERRWD (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG = '      Length needed is LENRW = I1, while LRW = I2.' 
      CALL XERRWD (MSG, 50, 103, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0) 
      endif 
   65 LSAVF = LEWT + N 
      LACOR = LSAVF + N 
      INSUFI = 0 
      IF (LIW .GE. LENIW) GO TO 70 
      INSUFI = 2 
      if (common_data%iprint == 1) then 
      MSG='DLSODAR-  Warning.. IWORK length is sufficient for now, but ' 
      CALL XERRWD (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG='      may not be later.  Integration will proceed anyway.   ' 
      CALL XERRWD (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG = '      Length needed is LENIW = I1, while LIW = I2.' 
      CALL XERRWD (MSG, 50, 104, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0) 
      endif 
   70 CONTINUE 
! Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1) 
      ATOLI = ATOL(1) 
      DO 75 I = 1,N 
        IF (ITOL .GE. 3) RTOLI = RTOL(I) 
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I) 
        IF (RTOLI .LT. 0.0D0) GO TO 619 
        IF (ATOLI .LT. 0.0D0) GO TO 620 
   75   CONTINUE 
      IF (ISTATE .EQ. 1) GO TO 100 
! if ISTATE = 3, set flag to signal parameter changes to DSTODA. -------
      JSTART = -1 
      IF (N .EQ. NYH) GO TO 200 
! NEQ was reduced.  zero part of yh to avoid undefined references. -----
      I1 = LYH + L*NYH 
      I2 = LYH + (MAXORD + 1)*NYH - 1 
      IF (I1 .GT. I2) GO TO 200 
      DO I = I1,I2 
        RWORK(I) = 0.0D0 
      END DO
      GO TO 200 
!-----------------------------------------------------------------------
! Block C.                                                              
! The next block is for the initial call only (ISTATE = 1).             
! It contains all remaining initializations, the initial call to F,     
! and the calculation of the initial step size.                         
! The error weights in EWT are inverted after being loaded.             
!-----------------------------------------------------------------------
  100 UROUND = DUMACH() 
      TN = T 
      TSW = T 
      MAXORD = MXORDN 
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110 
      TCRIT = RWORK(1) 
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GO TO 625 
      IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0)           &
     &   H0 = TCRIT - T                                                 
  110 JSTART = 0 
      NHNIL = 0 
      NST = 0 
      NJE = 0 
      NSLAST = 0 
      HU = 0.0D0 
      NQU = 0 
      MUSED = 0 
      MITER = 0 
      CCMAX = 0.3D0 
      MAXCOR = 3 
      MSBP = 20 
      MXNCF = 10 
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH 
      CALL F (NEQ, T, Y, RWORK(LF0), common_data%ierr) 
      if (common_data%ierr < 0) then; istate = -8; return; endif 
      NFE = 1 
! Load the initial value vector in YH. ---------------------------------
      DO I = 1,N 
        RWORK(I+LYH-1) = Y(I)
      END DO
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1 
      H = 1.0D0 
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT)) 
      DO I = 1,N 
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621 
        RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1) 
      END DO
!-----------------------------------------------------------------------
! The coding below computes the step size, H0, to be attempted on the   
! first step, unless the user has supplied a value for this.            
! First check that TOUT - T differs significantly from zero.            
! A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))          
! if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted    
! so as to be between 100*UROUND and 1.0E-3.                            
! Then the computed value H0 is given by:                               
!                                                                       
!   H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2                
!                                                                       
! where   w0     = MAX ( ABS(T), ABS(TOUT) ),                           
!         F      = the initial value of the vector f(t,y), and          
!         norm() = the weighted vector norm used throughout, given by   
!                  the DMNORM function routine, and weighted by the     
!                  tolerances initially loaded into the EWT array.      
! The sign of H0 is inferred from the initial values of TOUT and T.     
! ABS(H0) is made .le. ABS(TOUT-T) in any case.                         
!-----------------------------------------------------------------------
      IF (H0 .NE. 0.0D0) GO TO 180 
      TDIST = ABS(TOUT - T) 
      W0 = MAX(ABS(T),ABS(TOUT)) 
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622 
      TOL = RTOL(1) 
      IF (ITOL .LE. 2) GO TO 140 
      DO I = 1,N 
        TOL = MAX(TOL,RTOL(I)) 
      END DO
  140 IF (TOL .GT. 0.0D0) GO TO 160 
      ATOLI = ATOL(1) 
      DO 150 I = 1,N 
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I) 
        AYI = ABS(Y(I)) 
        IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI) 
  150   CONTINUE 
  160 TOL = MAX(TOL,100.0D0*UROUND) 
      TOL = MIN(TOL,0.001D0) 
      SUM = DMNORM (N, RWORK(LF0), RWORK(LEWT)) 
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2 
      H0 = 1.0D0/SQRT(SUM) 
      H0 = MIN(H0,TDIST) 
      H0 = SIGN(H0,TOUT-T) 
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
  180 RH = ABS(H0)*HMXI 
      IF (RH .GT. 1.0D0) H0 = H0/RH 
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0 
      DO I = 1,N 
        RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
      END DO
!                                                                       
! Check for a zero of g at T. ------------------------------------------
      IRFND = 0 
      TOUTC = TOUT 
      IF (NGC .EQ. 0) GO TO 270 
      CALL DRCHEK (1, G, NEQ, Y, RWORK(LYH), NYH,                       &
     &   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT, common_data)   
      if (common_data%ierr < 0) then; istate = -8; return; endif 
      IF (IRT .EQ. 0) GO TO 270 
      GO TO 632 
!-----------------------------------------------------------------------
! Block D.                                                              
! The next code block is for continuation calls only (ISTATE = 2 or 3)  
! and is to check stop conditions before taking a step.                 
! First, DRCHEK is called to check for a root within the last step      
! taken, other than the last root found there, if any.                  
! If ITASK = 2 or 5, and y(TN) has not yet been returned to the user    
! because of an intervening root, return through Block G.               
!-----------------------------------------------------------------------
  200 NSLAST = NST 
!                                                                       
      IRFP = IRFND 
      IF (NGC .EQ. 0) GO TO 205 
      IF (ITASK .EQ. 1 .OR. ITASK .EQ. 4) TOUTC = TOUT 
      CALL DRCHEK (2, G, NEQ, Y, RWORK(LYH), NYH,                       &
     &   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT, common_data)   
      if (common_data%ierr < 0) then; istate = -8; return; endif 
      IF (IRT .NE. 1) GO TO 205 
      IRFND = 1 
      ISTATE = 3 
      T = T0 
      GO TO 425 
  205 CONTINUE 
      IRFND = 0 
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 2) GO TO 400 
!                                                                       
      GO TO (210, 250, 220, 230, 240), ITASK 
  210 IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250 
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG, common_data) 
      IF (IFLAG .NE. 0) GO TO 627 
      T = TOUT 
      GO TO 420 
  220 TP = TN - HU*(1.0D0 + 100.0D0*UROUND) 
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623 
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250 
      T = TN 
      GO TO 400 
  230 TCRIT = RWORK(1) 
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624 
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625 
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245 
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG, common_data) 
      IF (IFLAG .NE. 0) GO TO 627 
      T = TOUT 
      GO TO 420 
  240 TCRIT = RWORK(1) 
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624 
  245 HMX = ABS(TN) + ABS(H) 
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX 
      IF (IHIT) T = TCRIT 
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 5) GO TO 400 
      IF (IHIT) GO TO 400 
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND) 
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250 
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND) 
      IF (ISTATE .EQ. 2 .AND. JSTART .GE. 0) JSTART = -2 
!-----------------------------------------------------------------------
! Block E.                                                              
! The next block is normally executed for all calls and contains        
! the call to the one-step core integrator DSTODA.                      
!                                                                       
! This is a looping point for the integration steps.                    
!                                                                       
! First check for too many steps being taken, update EWT (if not at     
! start of problem), check for too much accuracy being requested, and   
! check for H below the roundoff level in T.                            
!-----------------------------------------------------------------------
  250 CONTINUE 
      IF (METH .EQ. MUSED) GO TO 255 
      IF (INSUFR .EQ. 1) GO TO 550 
      IF (INSUFI .EQ. 1) GO TO 555 
  255 IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500 
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT)) 
      DO I = 1,N 
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510 
        RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1) 
      END DO
  270 TOLSF = UROUND*DMNORM (N, RWORK(LYH), RWORK(LEWT)) 
      IF (TOLSF .LE. 1.0D0) GO TO 280 
      TOLSF = TOLSF*2.0D0 
      IF (NST .EQ. 0) GO TO 626 
      GO TO 520 
  280 IF ((TN + H) .NE. TN) GO TO 290 
      NHNIL = NHNIL + 1 
      IF (NHNIL .GT. MXHNIL) GO TO 290 
      if (common_data%iprint == 1) then 
      MSG = 'DLSODAR-  Warning..Internal T(=R1) and H(=R2) are ' 
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG='      such that in the machine, T + H = T on the next step  ' 
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG = '     (H = step size). Solver will continue anyway.' 
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H) 
      IF (NHNIL .LT. MXHNIL) GO TO 290 
      MSG = 'DLSODAR-  Above warning has been issued I1 times. ' 
      CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      MSG = '     It will not be issued again for this problem.' 
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0) 
      endif 
  290 CONTINUE 
!-----------------------------------------------------------------------
!   CALL DSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPRJA,DSOLSY)
!-----------------------------------------------------------------------
      CALL DSTODA (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),    &
     &   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM),           &
     &   F, JAC, DPRJA, DSOLSY, common_data)                            
      if (common_data%ierr < 0) then; istate = -8; return; endif 
      KGO = 1 - KFLAG 
      GO TO (300, 530, 540), KGO 
!-----------------------------------------------------------------------
! Block F.                                                              
! The following block handles the case of a successful return from the  
! core integrator (KFLAG = 0).                                          
! If a method switch was just made, record TSW, reset MAXORD,           
! set JSTART to -1 to signal DSTODA to complete the switch,             
! and do extra printing of data if IXPR = 1.                            
! Then call DRCHEK to check for a root within the last step.            
! Then, if no root was found, check for stop conditions.                
!-----------------------------------------------------------------------
  300 INIT = 1 
      IF (METH .EQ. MUSED) GO TO 310 
      TSW = TN 
      MAXORD = MXORDN 
      IF (METH .EQ. 2) MAXORD = MXORDS 
      IF (METH .EQ. 2) RWORK(LWM) = SQRT(UROUND) 
      INSUFR = MIN(INSUFR,1) 
      INSUFI = MIN(INSUFI,1) 
      JSTART = -1 
      IF (IXPR .EQ. 0) GO TO 310 
      IF (METH .EQ. 2) THEN 
      MSG='DLSODAR- A switch to the BDF (stiff) method has occurred    ' 
      CALL XERRWD (MSG, 60, 105, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      ENDIF 
      IF (METH .EQ. 1) THEN 
      MSG='DLSODAR- A switch to the Adams (nonstiff) method occurred   ' 
      CALL XERRWD (MSG, 60, 106, 0, 0, 0, 0, 0, 0.0D0, 0.0D0) 
      ENDIF 
      MSG='     at T = R1,  tentative step size H = R2,  step NST = I1 ' 
      CALL XERRWD (MSG, 60, 107, 0, 1, NST, 0, 2, TN, H) 
  310 CONTINUE 
!                                                                       
      IF (NGC .EQ. 0) GO TO 315 
      CALL DRCHEK (3, G, NEQ, Y, RWORK(LYH), NYH,                       &
     &   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT, common_data)   
      if (common_data%ierr < 0) then; istate = -8; return; endif 
      IF (IRT .NE. 1) GO TO 315 
      IRFND = 1 
      ISTATE = 3 
      T = T0 
      GO TO 425 
  315 CONTINUE 
!                                                                       
      GO TO (320, 400, 330, 340, 350), ITASK 
! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
  320 IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250 
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG, common_data) 
      T = TOUT 
      GO TO 420 
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
  330 IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400 
      GO TO 250 
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary. 
  340 IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345 
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG, common_data) 
      T = TOUT 
      GO TO 420 
  345 HMX = ABS(TN) + ABS(H) 
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX 
      IF (IHIT) GO TO 400 
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND) 
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250 
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND) 
      IF (JSTART .GE. 0) JSTART = -2 
      GO TO 250 
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
  350 HMX = ABS(TN) + ABS(H) 
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX 
!-----------------------------------------------------------------------
! Block G.                                                              
! The following block handles all successful returns from DLSODAR.      
! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.        
! ISTATE is set to 2, and the optional outputs are loaded into the      
! work arrays before returning.                                         
!-----------------------------------------------------------------------
  400 DO I = 1,N 
        Y(I) = RWORK(I+LYH-1) 
      END DO
      T = TN 
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420 
      IF (IHIT) T = TCRIT 
  420 ISTATE = 2 
  425 CONTINUE 
      RWORK(11) = HU 
      RWORK(12) = H 
      RWORK(13) = TN 
      RWORK(15) = TSW 
      IWORK(11) = NST 
      IWORK(12) = NFE 
      IWORK(13) = NJE 
      IWORK(14) = NQU 
      IWORK(15) = NQ 
      IWORK(19) = MUSED 
      IWORK(20) = METH 
      IWORK(10) = NGE 
      TLAST = T 
      RETURN 
!-----------------------------------------------------------------------
! Block H.                                                              
! The following block handles all unsuccessful returns other than       
! those for illegal input.  First the error message routine is called.  
! If there was an error test or convergence test failure, IMXER is set. 
! Then Y is loaded from YH and T is set to TN.                          
! The optional outputs are loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
! 500  MSG = 'DLSODA-  At current T (=R1), MXSTEP (=I1) steps   '       
!      CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      taken on this call before reaching TOUT     '       
!      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)        
  500 common_data%error_message =                                       &
     & 'MXSTEP steps taken on this call before reaching TOUT.'          
      ISTATE = -1 
      GO TO 580 
! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
! 510  EWTI = RWORK(LEWT+I-1)                                           
!      MSG = 'DLSODA-  At T (=R1), EWT(I1) has become R2 .le. 0.'       
!      CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)              
  510 common_data%error_message =                                       &
     & 'EWT(i) became zero for some i during the integration.'          
      ISTATE = -6 
      GO TO 580 
! Too much accuracy requested for machine precision. -------------------
! 520  MSG = 'DLSODA-  At T (=R1), too much accuracy requested  '       
!      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      for precision of machine..  See TOLSF (=R2) '       
!      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)             
  520 common_data%error_message =                                       &
     & 'Too much accuracy requested for machine precision.'             
      RWORK(14) = TOLSF 
      ISTATE = -2 
      GO TO 580 
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
! 530  MSG = 'DLSODA-  At T(=R1) and step size H(=R2), the error'       
!      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      test failed repeatedly or with ABS(H) = HMIN'       
!      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)                 
  530 common_data%error_message =                                       &
     & 'Error test failed repeatedly or with ABS(H) = HMIN.'            
      ISTATE = -4 
      GO TO 560 
! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
! 540  MSG = 'DLSODA-  At T (=R1) and step size H (=R2), the    '       
!      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      corrector convergence failed repeatedly     '       
!      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG = '      or with ABS(H) = HMIN   '                           
!      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)                 
  540 common_data%error_message =                                       &
     & 'Convergence failed repeatedly or with ABS(H) = HMIN.'           
      ISTATE = -5 
      GO TO 560 
! RWORK length too small to proceed. -----------------------------------
! 550  MSG = 'DLSODA-  At current T(=R1), RWORK length too small'       
!      CALL XERRWD (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG='      to proceed.  The integration was otherwise successful.
!      CALL XERRWD (MSG, 60, 206, 0, 0, 0, 0, 1, TN, 0.0D0)             
  550 common_data%error_message =                                       &
     & 'RWORK length too small to proceed.'                             
      ISTATE = -7 
      GO TO 580 
! IWORK length too small to proceed. -----------------------------------
! 555  MSG = 'DLSODA-  At current T(=R1), IWORK length too small'       
!      CALL XERRWD (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)          
!      MSG='      to proceed.  The integration was otherwise successful.
!      CALL XERRWD (MSG, 60, 207, 0, 0, 0, 0, 1, TN, 0.0D0)             
  555 common_data%error_message =                                       &
     & 'IWORK length too small to proceed.'                             
      ISTATE = -7 
      GO TO 580 
! Compute IMXER if relevant. -------------------------------------------
  560 BIG = 0.0D0 
      IMXER = 1 
      DO 570 I = 1,N 
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1)) 
        IF (BIG .GE. SIZE) GO TO 570 
        BIG = SIZE 
        IMXER = I 
  570   CONTINUE 
      IWORK(16) = IMXER 
! Set Y vector, T, and optional outputs. -------------------------------
  580 DO I = 1,N 
        Y(I) = RWORK(I+LYH-1)
      END DO
      T = TN 
      RWORK(11) = HU 
      RWORK(12) = H 
      RWORK(13) = TN 
      RWORK(15) = TSW 
      IWORK(11) = NST 
      IWORK(12) = NFE 
      IWORK(13) = NJE 
      IWORK(14) = NQU 
      IWORK(15) = NQ 
      IWORK(19) = MUSED 
      IWORK(20) = METH 
      IWORK(10) = NGE 
      TLAST = T 
      RETURN 
!-----------------------------------------------------------------------
! Block I.                                                              
! The following block handles all error returns due to illegal input    
! (ISTATE = -3), as detected before calling the core integrator.        
! First the error message routine is called.  If the illegal input      
! is a negative ISTATE, the run is aborted (apparent infinite loop).    
!-----------------------------------------------------------------------
! 601  MSG = 'DLSODA-  ISTATE (=I1) illegal.'                           
!      CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)       
  601 common_data%error_message = 'ISTATE has illegal value.' 
      IF (ISTATE .LT. 0) GO TO 800 
      GO TO 700 
! 602  MSG = 'DLSODA-  ITASK (=I1) illegal. '                           
!      CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)        
  602 common_data%error_message = 'ITASK has illegal value.' 
      GO TO 700 
! 603  MSG = 'DLSODA-  ISTATE .gt. 1 but DLSODA not initialized.'       
!      CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)            
  603 common_data%error_message =                                       &
     & 'ISTATE > 1 but LSODA not initialized.'                          
      GO TO 700 
! 604  MSG = 'DLSODA-  NEQ (=I1) .lt. 1     '                           
!      CALL XERRWD (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)       
  604 common_data%error_message =                                       &
     & 'NEQ < 1.'                                                       
      GO TO 700 
! 605  MSG = 'DLSODA-  ISTATE = 3 and NEQ increased (I1 to I2). '       
!      CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)       
  605 common_data%error_message =                                       &
     & 'ISTATE == 3 and NEQ increased.'                                 
      GO TO 700 
! 606  MSG = 'DLSODA-  ITOL (=I1) illegal.  '                           
!      CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)         
  606 common_data%error_message =                                       &
     & 'ITOL has illegal value.'                                        
      GO TO 700 
! 607  MSG = 'DLSODA-  IOPT (=I1) illegal.  '                           
!      CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)         
  607 common_data%error_message =                                       &
     & 'IOPT has illegal value.'                                        
      GO TO 700 
! 608  MSG = 'DLSODA-  JT (=I1) illegal.    '                           
!      CALL XERRWD (MSG, 30, 8, 0, 1, JT, 0, 0, 0.0D0, 0.0D0)           
  608 common_data%error_message =                                       &
     & 'JT has illegal value.'                                          
      GO TO 700 
! 609  MSG = 'DLSODA-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '       
!      CALL XERRWD (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0D0, 0.0D0)      
  609 common_data%error_message =                                       &
     & 'ML has illegal value.'                                          
      GO TO 700 
! 610  MSG = 'DLSODA-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '       
!      CALL XERRWD (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0D0, 0.0D0)     
  610 common_data%error_message =                                       &
     & 'MU has illegal value.'                                          
      GO TO 700 
! 611  MSG = 'DLSODA-  IXPR (=I1) illegal.  '                           
!      CALL XERRWD (MSG, 30, 11, 0, 1, IXPR, 0, 0, 0.0D0, 0.0D0)        
  611 common_data%error_message =                                       &
     & 'IXPR has illegal value.'                                        
      GO TO 700 
! 612  MSG = 'DLSODA-  MXSTEP (=I1) .lt. 0  '                           
!      CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)      
  612 common_data%error_message =                                       &
     & 'MXSTEP < 0.'                                                    
      GO TO 700 
! 613  MSG = 'DLSODA-  MXHNIL (=I1) .lt. 0  '                           
!      CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)      
  613 common_data%error_message =                                       &
     & 'MXHNIL < 0'                                                     
      GO TO 700 
! 614  MSG = 'DLSODA-  TOUT (=R1) behind T (=R2)      '                 
!      CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)                
!      MSG = '      Integration direction is given by H0 (=R1)  '       
!      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)              
  614 common_data%error_message =                                       &
     & 'TOUT and T imply integration direction incompatible with H0'    
      GO TO 700 
! 615  MSG = 'DLSODA-  HMAX (=R1) .lt. 0.0  '                           
!      CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)            
  615 common_data%error_message =                                       &
     & 'HMAX < 0'                                                       
      GO TO 700 
! 616  MSG = 'DLSODA-  HMIN (=R1) .lt. 0.0  '                           
!      CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)            
  616 common_data%error_message =                                       &
     & 'HMIN < 0.0'                                                     
      GO TO 700 
! 617  MSG='DLSODA-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)
!      CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)     
  617 common_data%error_message =                                       &
     & 'RWORK is not long enough.'                                      
      GO TO 700 
! 618  MSG='DLSODA-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)
!      CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)     
  618 common_data%error_message =                                       &
     & 'IWORK is not long enough.'                                      
      GO TO 700 
! 619  MSG = 'DLSODA-  RTOL(I1) is R1 .lt. 0.0        '                 
!      CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)           
  619 common_data%error_message =                                       &
     & 'RTOL < 0.0'                                                     
      GO TO 700 
! 620  MSG = 'DLSODA-  ATOL(I1) is R1 .lt. 0.0        '                 
!      CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)           
  620 common_data%error_message =                                       &
     & 'ATOL < 0.0'                                                     
      GO TO 700 
  621 EWTI = RWORK(LEWT+I-1) 
!      MSG = 'DLSODA-  EWT(I1) is R1 .le. 0.0         '                 
!      CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)            
      common_data%error_message =                                       &
     & 'EWT(i) for some i is < 0.0'                                     
      GO TO 700 
! 622  MSG='DLSODA-  TOUT(=R1) too close to T(=R2) to start integration.
!      CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)                
  622 common_data%error_message =                                       &
     & 'TOUT is too close to T to start integration.'                   
      GO TO 700 
! 623  MSG='DLSODA-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  
!      CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)           
  623 common_data%error_message =                                       &
     & 'TOUT switches the direction of integration.'                    
      GO TO 700 
! 624  MSG='DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   
!      CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)              
  624 common_data%error_message =                                       &
     & 'ITASK = 4 or 5 and TCRIT behind TCUR.'                          
      GO TO 700 
! 625  MSG='DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   
!      CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)            
  625 common_data%error_message =                                       &
     & 'ITASK = 4 or 5 and TCRIT behind TOUT.'                          
      GO TO 700 
! 626  MSG = 'DLSODA-  At start of problem, too much accuracy   '       
!      CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)           
!      MSG='      requested for precision of machine..  See TOLSF (=R1) 
!      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)           
  626 common_data%error_message =                                       &
     & 'At start of problem, too much accuracy is requested.'           
      RWORK(14) = TOLSF 
      GO TO 700 
! 627  MSG = 'DLSODA-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'       
!      CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)        
  627 common_data%error_message =                                       &
     & 'Trouble in DINTDY.'                                             
      GO TO 700 
! 628  MSG = 'DLSODA-  MXORDN (=I1) .lt. 0  '                           
!      CALL XERRWD (MSG, 30, 28, 0, 1, MXORDN, 0, 0, 0.0D0, 0.0D0)      
  628 common_data%error_message =                                       &
     & 'MXORDN < 0'                                                     
      GO TO 700 
! 629  MSG = 'DLSODA-  MXORDS (=I1) .lt. 0  '                           
!      CALL XERRWD (MSG, 30, 29, 0, 1, MXORDS, 0, 0, 0.0D0, 0.0D0)      
  629 common_data%error_message =                                       &
     & 'MXORDS < 0'                                                     
      GO TO 700 
! 630  MSG = 'DLSODAR-  NG (=I1) .lt. 0     '                           
!      CALL XERRWD (MSG, 30, 30, 0, 1, NG, 0, 0, 0.0D0, 0.0D0)          
  630 common_data%error_message =                                       &
     & 'NG < 0'                                                         
      GO TO 700 
! 631  MSG = 'DLSODAR-  NG changed (from I1 to I2) illegally,   '       
!      CALL XERRWD (MSG, 50, 31, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)           
!      MSG = '      i.e. not immediately after a root was found.'       
!      CALL XERRWD (MSG, 50, 31, 0, 2, NGC, NG, 0, 0.0D0, 0.0D0)        
  631 common_data%error_message =                                       &
     & 'NG changed illegally '//                                        &
     & '(i.e. not immediately after a root was found)'                  
      GO TO 700 
! 632  MSG = 'DLSODAR-  One or more components of g has a root  '       
!      CALL XERRWD (MSG, 50, 32, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)           
!      MSG = '      too near to the initial point.    '                 
!      CALL XERRWD (MSG, 40, 32, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)           
  632 common_data%error_message =                                       &
     & 'One or more components of g has a root'//                       &
     & ' too near to the initial point.'                                
!                                                                       
  700 ISTATE = -3 
      RETURN 
!                                                                       
! 800  MSG = 'DLSODAR-  Run aborted.. apparent infinite loop.   '       
!      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)          
  800 BLOCK
        CHARACTER(LEN=200) :: tmp_msg
        tmp_msg = TRIM(common_data%error_message) // ' Run aborted. apparent infinite loop.'
        common_data%error_message = TRIM(tmp_msg)
      END BLOCK
      RETURN
!----------------------- End of Subroutine DLSODAR ---------------------
      END                                           
                                                                        
!=======================================================================
!                 START OF ODEPACK SUB1                                 
!=======================================================================
!DECK DUMACH                                                            
      DOUBLE PRECISION FUNCTION DUMACH () 
!***BEGIN PROLOGUE  DUMACH                                              
!***PURPOSE  Compute the unit roundoff of the machine.                  
!***CATEGORY  R1                                                        
!***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)                     
!***KEYWORDS  MACHINE CONSTANTS                                         
!***AUTHOR  Hindmarsh, Alan C., (LLNL)                                  
!***DESCRIPTION                                                         
! *Usage:                                                               
!        DOUBLE PRECISION  A, DUMACH                                    
!        A = DUMACH()                                                   
!                                                                       
! *Function Return Values:                                              
!     A : the unit roundoff of the machine.                             
!                                                                       
! *Description:                                                         
!     The unit roundoff is defined as the smallest positive machine     
!     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH 
!     in a machine-independent manner.                                  
!                                                                       
!***REFERENCES  (NONE)                                                  
!***ROUTINES CALLED  DUMSUM                                             
!***REVISION HISTORY  (YYYYMMDD)                                        
!   19930216  DATE WRITTEN                                              
!   19930818  Added SLATEC-format prologue.  (FNF)                      
!   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)      
!***END PROLOGUE  DUMACH                                                
!                                                                       
      DOUBLE PRECISION U, COMP 
!***FIRST EXECUTABLE STATEMENT  DUMACH                                  
      U = 1.0D0 
   10 U = U*0.5D0 
      CALL DUMSUM(1.0D0, U, COMP) 
      IF (COMP .NE. 1.0D0) GO TO 10 
      DUMACH = U*2.0D0 
      RETURN 
!----------------------- End of Function DUMACH ------------------------
      END                                           
      SUBROUTINE DUMSUM(A,B,C) 
!     Routine to force normal storing of A + B, for DUMACH.             
      DOUBLE PRECISION A, B, C 
      C = A + B 
      RETURN 
      END                                           
!DECK DCFODE                                                            
      SUBROUTINE DCFODE (METH, ELCO, TESCO) 
!***BEGIN PROLOGUE  DCFODE                                              
!***SUBSIDIARY                                                          
!***PURPOSE  Set ODE integrator coefficients.                           
!***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)                     
!***AUTHOR  Hindmarsh, Alan C., (LLNL)                                  
!***DESCRIPTION                                                         
!                                                                       
!  DCFODE is called by the integrator routine to set coefficients       
!  needed there.  The coefficients for the current method, as           
!  given by the value of METH, are set for all orders and saved.        
!  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.  
!  (A smaller value of the maximum order is also allowed.)              
!  DCFODE is called once at the beginning of the problem,               
!  and is not called again unless and until METH is changed.            
!                                                                       
!  The ELCO array contains the basic method coefficients.               
!  The coefficients el(i), 1 .le. i .le. nq+1, for the method of        
!  order nq are stored in ELCO(i,nq).  They are given by a genetrating  
!  polynomial, i.e.,                                                    
!      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.                   
!  For the implicit Adams methods, l(x) is given by                     
!      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.  
!  For the BDF methods, l(x) is given by                                
!      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,                               
!  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).              
!                                                                       
!  The TESCO array contains test constants used for the                 
!  local error test and the selection of step size and/or order.        
!  At order nq, TESCO(k,nq) is used for the selection of step           
!  size at order nq - 1 if k = 1, at order nq if k = 2, and at order    
!  nq + 1 if k = 3.                                                     
!                                                                       
!***SEE ALSO  DLSODE                                                    
!***ROUTINES CALLED  (NONE)                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   791129  DATE WRITTEN                                                
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)             
!   890503  Minor cosmetic changes.  (FNF)                              
!   930809  Renamed to allow single/double precision versions. (ACH)    
!***END PROLOGUE  DCFODE                                                
!**End                                                                  
      INTEGER METH 
      INTEGER I, IB, NQ, NQM1, NQP1 
      DOUBLE PRECISION ELCO, TESCO 
      DOUBLE PRECISION AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,               &
     &   RQFAC, RQ1FAC, TSIGN, XPIN                                     
      DIMENSION ELCO(13,12), TESCO(3,12) 
      DIMENSION PC(12) 
!                                                                       
!***FIRST EXECUTABLE STATEMENT  DCFODE                                  
      GO TO (100, 200), METH 
!                                                                       
  100 ELCO(1,1) = 1.0D0 
      ELCO(2,1) = 1.0D0 
      TESCO(1,1) = 0.0D0 
      TESCO(2,1) = 2.0D0 
      TESCO(1,2) = 1.0D0 
      TESCO(3,12) = 0.0D0 
      PC(1) = 1.0D0 
      RQFAC = 1.0D0 
      DO 140 NQ = 2,12 
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial          
!     p(x) = (x+1)*(x+2)*...*(x+nq-1).                                  
! Initially, p(x) = 1.                                                  
!-----------------------------------------------------------------------
        RQ1FAC = RQFAC 
        RQFAC = RQFAC/NQ 
        NQM1 = NQ - 1 
        FNQM1 = NQM1 
        NQP1 = NQ + 1 
! Form coefficients of p(x)*(x+nq-1). ----------------------------------
        PC(NQ) = 0.0D0 
        DO IB = 1,NQM1 
          I = NQP1 - IB 
          PC(I) = PC(I-1) + FNQM1*PC(I) 
        END DO
        PC(1) = FNQM1*PC(1) 
! Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        PINT = PC(1) 
        XPIN = PC(1)/2.0D0 
        TSIGN = 1.0D0 
        DO I = 2,NQ 
          TSIGN = -TSIGN 
          PINT = PINT + TSIGN*PC(I)/I 
          XPIN = XPIN + TSIGN*PC(I)/(I+1) 
        END DO
! Store coefficients in ELCO and TESCO. --------------------------------
        ELCO(1,NQ) = PINT*RQ1FAC 
        ELCO(2,NQ) = 1.0D0 
        DO I = 2,NQ 
          ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
        END DO
        AGAMQ = RQFAC*XPIN 
        RAGQ = 1.0D0/AGAMQ 
        TESCO(2,NQ) = RAGQ 
        IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1 
        TESCO(3,NQM1) = RAGQ 
  140   CONTINUE 
      RETURN 
!                                                                       
  200 PC(1) = 1.0D0 
      RQ1FAC = 1.0D0 
      DO 230 NQ = 1,5 
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial          
!     p(x) = (x+1)*(x+2)*...*(x+nq).                                    
! Initially, p(x) = 1.                                                  
!-----------------------------------------------------------------------
        FNQ = NQ 
        NQP1 = NQ + 1 
! Form coefficients of p(x)*(x+nq). ------------------------------------
        PC(NQP1) = 0.0D0 
        DO IB = 1,NQ 
          I = NQ + 2 - IB 
          PC(I) = PC(I-1) + FNQ*PC(I)
        END DO
        PC(1) = FNQ*PC(1) 
! Store coefficients in ELCO and TESCO. --------------------------------
        DO I = 1,NQP1 
          ELCO(I,NQ) = PC(I)/PC(2)
        END DO
        ELCO(2,NQ) = 1.0D0 
        TESCO(1,NQ) = RQ1FAC 
        TESCO(2,NQ) = NQP1/ELCO(1,NQ) 
        TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ) 
        RQ1FAC = RQ1FAC/FNQ 
  230   CONTINUE 
      RETURN 
!----------------------- END OF SUBROUTINE DCFODE ----------------------
      END                                           
!DECK DINTDY                                                            
      SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG, common_data) 
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_common 
      type(odepack_common_data), target, intent(inout) :: common_data 
!***BEGIN PROLOGUE  DINTDY                                              
!***SUBSIDIARY                                                          
!***PURPOSE  Interpolate solution derivatives.                          
!***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)                     
!***AUTHOR  Hindmarsh, Alan C., (LLNL)                                  
!***DESCRIPTION                                                         
!                                                                       
!  DINTDY computes interpolated values of the K-th derivative of the    
!  dependent variable vector y, and stores it in DKY.  This routine     
!  is called within the package with K = 0 and T = TOUT, but may        
!  also be called by the user for any K up to the current order.        
!  (See detailed instructions in the usage documentation.)              
!                                                                       
!  The computed values in DKY are gotten by interpolation using the     
!  Nordsieck history array YH.  This array corresponds uniquely to a    
!  vector-valued polynomial of degree NQCUR or less, and DKY is set     
!  to the K-th derivative of this polynomial at T.                      
!  The formula for DKY is:                                              
!               q                                                       
!   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)      
!              j=K                                                      
!  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR. 
!  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are         
!  communicated by COMMON.  The above sum is done in reverse order.     
!  IFLAG is returned negative if either K or T is out of bounds.        
!                                                                       
!***SEE ALSO  DLSODE                                                    
!***ROUTINES CALLED  XERRWD                                             
!***COMMON BLOCKS    DLS001                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   791129  DATE WRITTEN                                                
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)             
!   890503  Minor cosmetic changes.  (FNF)                              
!   930809  Renamed to allow single/double precision versions. (ACH)    
!   010418  Reduced size of Common block /DLS001/. (ACH)                
!   031105  Restored 'own' variables to Common block /DLS001/, to       
!           enable interrupt/restart feature. (ACH)                     
!   050427  Corrected roundoff decrement in TP. (ACH)                   
!***END PROLOGUE  DINTDY                                                
!**End                                                                  
      INTEGER K, NYH, IFLAG 
      DOUBLE PRECISION T, YH, DKY 
      DIMENSION YH(NYH,*), DKY(*) 
      INTEGER, pointer :: IOWND(:), IOWNS(:),                           &
     &   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                     &
     &   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,               &
     &   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU         
      DOUBLE PRECISION, pointer :: ROWNS(:),                            &
     &   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND                  
!      COMMON /DLS001/ ROWNS(209),                                      
!     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,                
!     2   IOWND(6), IOWNS(6),                                           
!     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                    
!     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,              
!     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU        
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1 
      DOUBLE PRECISION C, R, S, TP 
      CHARACTER*80 MSG 
!                                                                       
!     Common block pointers                                             
      type(DLS001_type), pointer :: DLS001 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      DLS001 => common_data%DLS001 
                                                                        
      tmp_ptr = c_loc(DLS001%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWNS,[209]) 
      tmp_ptr = c_loc(DLS001%reals(210)) 
      call c_f_pointer(tmp_ptr,CCMAX) 
      tmp_ptr = c_loc(DLS001%reals(211)) 
      call c_f_pointer(tmp_ptr,EL0) 
      tmp_ptr = c_loc(DLS001%reals(212)) 
      call c_f_pointer(tmp_ptr,H) 
      tmp_ptr = c_loc(DLS001%reals(213)) 
      call c_f_pointer(tmp_ptr,HMIN) 
      tmp_ptr = c_loc(DLS001%reals(214)) 
      call c_f_pointer(tmp_ptr,HMXI) 
      tmp_ptr = c_loc(DLS001%reals(215)) 
      call c_f_pointer(tmp_ptr,HU) 
      tmp_ptr = c_loc(DLS001%reals(216)) 
      call c_f_pointer(tmp_ptr,RC) 
      tmp_ptr = c_loc(DLS001%reals(217)) 
      call c_f_pointer(tmp_ptr,TN) 
      tmp_ptr = c_loc(DLS001%reals(218)) 
      call c_f_pointer(tmp_ptr,UROUND) 
                                                                        
      tmp_ptr = c_loc(DLS001%ints(1)) 
      call c_f_pointer(tmp_ptr,IOWND,[6]) 
      tmp_ptr = c_loc(DLS001%ints(7)) 
      call c_f_pointer(tmp_ptr,IOWNS,[6]) 
      tmp_ptr = c_loc(DLS001%ints(13)) 
      call c_f_pointer(tmp_ptr,ICF) 
      tmp_ptr = c_loc(DLS001%ints(14)) 
      call c_f_pointer(tmp_ptr,IERPJ) 
      tmp_ptr = c_loc(DLS001%ints(15)) 
      call c_f_pointer(tmp_ptr,IERSL) 
      tmp_ptr = c_loc(DLS001%ints(16)) 
      call c_f_pointer(tmp_ptr,JCUR) 
      tmp_ptr = c_loc(DLS001%ints(17)) 
      call c_f_pointer(tmp_ptr,JSTART) 
      tmp_ptr = c_loc(DLS001%ints(18)) 
      call c_f_pointer(tmp_ptr,KFLAG) 
      tmp_ptr = c_loc(DLS001%ints(19)) 
      call c_f_pointer(tmp_ptr,L) 
      tmp_ptr = c_loc(DLS001%ints(20)) 
      call c_f_pointer(tmp_ptr,LYH) 
      tmp_ptr = c_loc(DLS001%ints(21)) 
      call c_f_pointer(tmp_ptr,LEWT) 
      tmp_ptr = c_loc(DLS001%ints(22)) 
      call c_f_pointer(tmp_ptr,LACOR) 
      tmp_ptr = c_loc(DLS001%ints(23)) 
      call c_f_pointer(tmp_ptr,LSAVF) 
      tmp_ptr = c_loc(DLS001%ints(24)) 
      call c_f_pointer(tmp_ptr,LWM) 
      tmp_ptr = c_loc(DLS001%ints(25)) 
      call c_f_pointer(tmp_ptr,LIWM) 
      tmp_ptr = c_loc(DLS001%ints(26)) 
      call c_f_pointer(tmp_ptr,METH) 
      tmp_ptr = c_loc(DLS001%ints(27)) 
      call c_f_pointer(tmp_ptr,MITER) 
      tmp_ptr = c_loc(DLS001%ints(28)) 
      call c_f_pointer(tmp_ptr,MAXORD) 
      tmp_ptr = c_loc(DLS001%ints(29)) 
      call c_f_pointer(tmp_ptr,MAXCOR) 
      tmp_ptr = c_loc(DLS001%ints(30)) 
      call c_f_pointer(tmp_ptr,MSBP) 
      tmp_ptr = c_loc(DLS001%ints(31)) 
      call c_f_pointer(tmp_ptr,MXNCF) 
      tmp_ptr = c_loc(DLS001%ints(32)) 
      call c_f_pointer(tmp_ptr,N) 
      tmp_ptr = c_loc(DLS001%ints(33)) 
      call c_f_pointer(tmp_ptr,NQ) 
      tmp_ptr = c_loc(DLS001%ints(34)) 
      call c_f_pointer(tmp_ptr,NST) 
      tmp_ptr = c_loc(DLS001%ints(35)) 
      call c_f_pointer(tmp_ptr,NFE) 
      tmp_ptr = c_loc(DLS001%ints(36)) 
      call c_f_pointer(tmp_ptr,NJE) 
      tmp_ptr = c_loc(DLS001%ints(37)) 
      call c_f_pointer(tmp_ptr,NQU) 
!                                                                       
!***FIRST EXECUTABLE STATEMENT  DINTDY                                  
      IFLAG = 0 
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80 
      TP = TN - HU -  100.0D0*UROUND*SIGN(ABS(TN) + ABS(HU), HU) 
      IF ((T-TP)*(T-TN) .GT. 0.0D0) GO TO 90 
!                                                                       
      S = (T - TN)/H 
      IC = 1 
      IF (K .EQ. 0) GO TO 15 
      JJ1 = L - K 
      DO JJ = JJ1,NQ 
        IC = IC*JJ
      END DO
   15 C = IC 
      DO I = 1,N 
        DKY(I) = C*YH(I,L)
      END DO
      IF (K .EQ. NQ) GO TO 55 
      JB2 = NQ - K 
      DO 50 JB = 1,JB2 
        J = NQ - JB 
        JP1 = J + 1 
        IC = 1 
        IF (K .EQ. 0) GO TO 35 
        JJ1 = JP1 - K 
        DO JJ = JJ1,J 
          IC = IC*JJ
        END DO
   35   C = IC 
        DO I = 1,N 
          DKY(I) = C*YH(I,JP1) + S*DKY(I)
        END DO
   50   CONTINUE 
      IF (K .EQ. 0) RETURN 
   55 R = H**(-K) 
      DO I = 1,N 
        DKY(I) = R*DKY(I)
      END DO
      RETURN 
!                                                                       
   80 MSG = 'DINTDY-  K (=I1) illegal      ' 
      if (common_data%iprint == 1) then 
      CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0D0, 0.0D0) 
      endif 
      IFLAG = -1 
      RETURN 
   90 MSG = 'DINTDY-  T (=R1) illegal      ' 
      if (common_data%iprint == 1) then 
      CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0D0) 
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      ' 
      CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN) 
      endif 
      IFLAG = -2 
      RETURN 
!----------------------- END OF SUBROUTINE DINTDY ----------------------
      END                                           
!DECK DSOLSY                                                            
      SUBROUTINE DSOLSY (WM, IWM, X, TEM, common_data) 
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_common 
      type(odepack_common_data), target, intent(inout) :: common_data 
!***BEGIN PROLOGUE  DSOLSY                                              
!***SUBSIDIARY                                                          
!***PURPOSE  ODEPACK linear system solver.                              
!***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)                     
!***AUTHOR  Hindmarsh, Alan C., (LLNL)                                  
!***DESCRIPTION                                                         
!                                                                       
!  This routine manages the solution of the linear system arising from  
!  a chord iteration.  It is called if MITER .ne. 0.                    
!  If MITER is 1 or 2, it calls DGESL to accomplish this.               
!  If MITER = 3 it updates the coefficient h*EL0 in the diagonal        
!  matrix, and then computes the solution.                              
!  If MITER is 4 or 5, it calls DGBSL.                                  
!  Communication with DSOLSY uses the following variables:              
!  WM    = real work space containing the inverse diagonal matrix if    
!          MITER = 3 and the LU decomposition of the matrix otherwise.  
!          Storage of matrix elements starts at WM(3).                  
!          WM also contains the following matrix-related data:          
!          WM(1) = SQRT(UROUND) (not used here),                        
!          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3. 
!  IWM   = integer work space containing pivot information, starting at 
!          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band  
!          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.   
!  X     = the right-hand side vector on input, and the solution vector 
!          on output, of length N.                                      
!  TEM   = vector of work space of length N, not used in this version.  
!  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.  
!          IERSL = 1 if a singular matrix arose with MITER = 3.         
!  This routine also uses the COMMON variables EL0, H, MITER, and N.    
!                                                                       
!***SEE ALSO  DLSODE                                                    
!***ROUTINES CALLED  DGBSL, DGESL                                       
!***COMMON BLOCKS    DLS001                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   791129  DATE WRITTEN                                                
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)             
!   890503  Minor cosmetic changes.  (FNF)                              
!   930809  Renamed to allow single/double precision versions. (ACH)    
!   010418  Reduced size of Common block /DLS001/. (ACH)                
!   031105  Restored 'own' variables to Common block /DLS001/, to       
!           enable interrupt/restart feature. (ACH)                     
!***END PROLOGUE  DSOLSY                                                
!**End                                                                  
      INTEGER IWM 
      DOUBLE PRECISION WM, X, TEM 
      DIMENSION WM(*), IWM(*), X(*), TEM(*) 
      INTEGER, pointer :: IOWND(:), IOWNS(:),                           &
     &   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                     &
     &   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,               &
     &   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU         
      DOUBLE PRECISION, pointer :: ROWNS(:),                            &
     &   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND                  
!      COMMON /DLS001/ ROWNS(209),                                      
!     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,                
!     2   IOWND(6), IOWNS(6),                                           
!     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                    
!     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,              
!     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU        
      INTEGER I, MEBAND, ML, MU 
      DOUBLE PRECISION DI, HL0, PHL0, R 
!     Common block pointers                                             
      type(DLS001_type), pointer :: DLS001 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      IF (TEM(1) == TEM(1)) CONTINUE ! To supress compiler warnings
      DLS001 => common_data%DLS001 
                                                                        
      tmp_ptr = c_loc(DLS001%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWNS,[209]) 
      tmp_ptr = c_loc(DLS001%reals(210)) 
      call c_f_pointer(tmp_ptr,CCMAX) 
      tmp_ptr = c_loc(DLS001%reals(211)) 
      call c_f_pointer(tmp_ptr,EL0) 
      tmp_ptr = c_loc(DLS001%reals(212)) 
      call c_f_pointer(tmp_ptr,H) 
      tmp_ptr = c_loc(DLS001%reals(213)) 
      call c_f_pointer(tmp_ptr,HMIN) 
      tmp_ptr = c_loc(DLS001%reals(214)) 
      call c_f_pointer(tmp_ptr,HMXI) 
      tmp_ptr = c_loc(DLS001%reals(215)) 
      call c_f_pointer(tmp_ptr,HU) 
      tmp_ptr = c_loc(DLS001%reals(216)) 
      call c_f_pointer(tmp_ptr,RC) 
      tmp_ptr = c_loc(DLS001%reals(217)) 
      call c_f_pointer(tmp_ptr,TN) 
      tmp_ptr = c_loc(DLS001%reals(218)) 
      call c_f_pointer(tmp_ptr,UROUND) 
                                                                        
      tmp_ptr = c_loc(DLS001%ints(1)) 
      call c_f_pointer(tmp_ptr,IOWND,[6]) 
      tmp_ptr = c_loc(DLS001%ints(7)) 
      call c_f_pointer(tmp_ptr,IOWNS,[6]) 
      tmp_ptr = c_loc(DLS001%ints(13)) 
      call c_f_pointer(tmp_ptr,ICF) 
      tmp_ptr = c_loc(DLS001%ints(14)) 
      call c_f_pointer(tmp_ptr,IERPJ) 
      tmp_ptr = c_loc(DLS001%ints(15)) 
      call c_f_pointer(tmp_ptr,IERSL) 
      tmp_ptr = c_loc(DLS001%ints(16)) 
      call c_f_pointer(tmp_ptr,JCUR) 
      tmp_ptr = c_loc(DLS001%ints(17)) 
      call c_f_pointer(tmp_ptr,JSTART) 
      tmp_ptr = c_loc(DLS001%ints(18)) 
      call c_f_pointer(tmp_ptr,KFLAG) 
      tmp_ptr = c_loc(DLS001%ints(19)) 
      call c_f_pointer(tmp_ptr,L) 
      tmp_ptr = c_loc(DLS001%ints(20)) 
      call c_f_pointer(tmp_ptr,LYH) 
      tmp_ptr = c_loc(DLS001%ints(21)) 
      call c_f_pointer(tmp_ptr,LEWT) 
      tmp_ptr = c_loc(DLS001%ints(22)) 
      call c_f_pointer(tmp_ptr,LACOR) 
      tmp_ptr = c_loc(DLS001%ints(23)) 
      call c_f_pointer(tmp_ptr,LSAVF) 
      tmp_ptr = c_loc(DLS001%ints(24)) 
      call c_f_pointer(tmp_ptr,LWM) 
      tmp_ptr = c_loc(DLS001%ints(25)) 
      call c_f_pointer(tmp_ptr,LIWM) 
      tmp_ptr = c_loc(DLS001%ints(26)) 
      call c_f_pointer(tmp_ptr,METH) 
      tmp_ptr = c_loc(DLS001%ints(27)) 
      call c_f_pointer(tmp_ptr,MITER) 
      tmp_ptr = c_loc(DLS001%ints(28)) 
      call c_f_pointer(tmp_ptr,MAXORD) 
      tmp_ptr = c_loc(DLS001%ints(29)) 
      call c_f_pointer(tmp_ptr,MAXCOR) 
      tmp_ptr = c_loc(DLS001%ints(30)) 
      call c_f_pointer(tmp_ptr,MSBP) 
      tmp_ptr = c_loc(DLS001%ints(31)) 
      call c_f_pointer(tmp_ptr,MXNCF) 
      tmp_ptr = c_loc(DLS001%ints(32)) 
      call c_f_pointer(tmp_ptr,N) 
      tmp_ptr = c_loc(DLS001%ints(33)) 
      call c_f_pointer(tmp_ptr,NQ) 
      tmp_ptr = c_loc(DLS001%ints(34)) 
      call c_f_pointer(tmp_ptr,NST) 
      tmp_ptr = c_loc(DLS001%ints(35)) 
      call c_f_pointer(tmp_ptr,NFE) 
      tmp_ptr = c_loc(DLS001%ints(36)) 
      call c_f_pointer(tmp_ptr,NJE) 
      tmp_ptr = c_loc(DLS001%ints(37)) 
      call c_f_pointer(tmp_ptr,NQU) 
!                                                                       
!***FIRST EXECUTABLE STATEMENT  DSOLSY                                  
      IERSL = 0 
      GO TO (100, 100, 300, 400, 400), MITER 
! 100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)                          
  100 call dgetrs ('N', n, 1, wm(3), n, iwm(21), x, n, ier) 
      RETURN 
!                                                                       
  300 PHL0 = WM(2) 
      HL0 = H*EL0 
      WM(2) = HL0 
      IF (HL0 .EQ. PHL0) GO TO 330 
      R = HL0/PHL0 
      DO I = 1,N 
        DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2)) 
        IF (ABS(DI) .EQ. 0.0D0) GO TO 390 
        WM(I+2) = 1.0D0/DI 
      END DO
  330 DO I = 1,N 
        X(I) = WM(I+2)*X(I)
      END DO
      RETURN 
  390 IERSL = 1 
      RETURN 
!                                                                       
  400 ML = IWM(1) 
      MU = IWM(2) 
      MEBAND = 2*ML + MU + 1 
!      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)             
      call dgbtrs ('N', n, ml, mu, 1, wm(3), meband, iwm(21), x, n, ier) 
      RETURN 
!----------------------- END OF SUBROUTINE DSOLSY ----------------------
      END                                           
!DECK DEWSET                                                            
      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT) 
!***BEGIN PROLOGUE  DEWSET                                              
!***SUBSIDIARY                                                          
!***PURPOSE  Set error weight vector.                                   
!***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)                     
!***AUTHOR  Hindmarsh, Alan C., (LLNL)                                  
!***DESCRIPTION                                                         
!                                                                       
!  This subroutine sets the error weight vector EWT according to        
!      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,           
!  with the subscript on RTOL and/or ATOL possibly replaced by 1 above, 
!  depending on the value of ITOL.                                      
!                                                                       
!***SEE ALSO  DLSODE                                                    
!***ROUTINES CALLED  (NONE)                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   791129  DATE WRITTEN                                                
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)             
!   890503  Minor cosmetic changes.  (FNF)                              
!   930809  Renamed to allow single/double precision versions. (ACH)    
!***END PROLOGUE  DEWSET                                                
!**End                                                                  
      INTEGER N, ITOL 
      INTEGER I 
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT 
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N) 
!                                                                       
!***FIRST EXECUTABLE STATEMENT  DEWSET                                  
      GO TO (10, 20, 30, 40), ITOL 
   10 CONTINUE 
      DO I = 1,N 
        EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      END DO
      RETURN 
   20 CONTINUE 
      DO I = 1,N 
        EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      END DO
      RETURN 
   30 CONTINUE 
      DO I = 1,N 
        EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      END DO
      RETURN 
   40 CONTINUE 
      DO I = 1,N 
        EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      END DO
      RETURN 
!----------------------- END OF SUBROUTINE DEWSET ----------------------
      END                                           
!DECK DSTODA                                                            
      SUBROUTINE DSTODA (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR,         &
     &   WM, IWM, F, JAC, PJAC, SLVS, common_data)                      
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_common 
      type(odepack_common_data), target, intent(inout) :: common_data 
      EXTERNAL F, JAC, PJAC, SLVS 
      INTEGER NEQ, NYH, IWM 
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM 
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),       &
     &   ACOR(*), WM(*), IWM(*)                                         
      INTEGER, pointer ::                                               &
     &   IOWND(:), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,                 &
     &   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                     &
     &   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,               &
     &   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU         
      INTEGER, pointer ::                                               &
     &   IOWND2(:), ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS         
      DOUBLE PRECISION, pointer ::                                      &
     &   CONIT, CRATE, EL(:), ELCO(:,:), HOLD, RMAX, TESCO(:,:),        &
     &   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND                  
      DOUBLE PRECISION, pointer ::                                      &
     &   ROWND2, CM1(:), CM2(:), PDEST, PDLAST, RATIO,                  &
     &   PDNORM                                                         
!      COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12),               
!     1   HOLD, RMAX, TESCO(3,12),                                      
!     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,                
!     3   IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,                
!     4   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                    
!     5   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,              
!     6   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU        
!      COMMON /DLSA01/ ROWND2, CM1(12), CM2(5), PDEST, PDLAST, RATIO,   
!     1   PDNORM,                                                       
!     2   IOWND2(3), ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS        
      INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ 
      INTEGER LM1, LM1P1, LM2, LM2P1, NQM1, NQM2 
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,&
     &   R, RH, RHDN, RHSM, RHUP, TOLD, DMNORM                          
      DOUBLE PRECISION ALPHA, DM1,DM2, EXM1,EXM2,                       &
     &   PDH, PNORM, RATE, RH1, RH1IT, RH2, RM, SM1(12)                 
      SAVE SM1 
      DATA SM1/0.5D0, 0.575D0, 0.55D0, 0.45D0, 0.35D0, 0.25D0,          &
     &   0.20D0, 0.15D0, 0.10D0, 0.075D0, 0.050D0, 0.025D0/             
!-----------------------------------------------------------------------
! DSTODA performs one step of the integration of an initial value       
! problem for a system of ordinary differential equations.              
! Note: DSTODA is independent of the value of the iteration method      
! indicator MITER, when this is .ne. 0, and hence is independent        
! of the type of chord method used, or the Jacobian structure.          
! Communication with DSTODA is done with the following variables:       
!                                                                       
! Y      = an array of length .ge. N used as the Y argument in          
!          all calls to F and JAC.                                      
! NEQ    = integer array containing problem size in NEQ(1), and         
!          passed as the NEQ argument in all calls to F and JAC.        
! YH     = an NYH by LMAX array containing the dependent variables      
!          and their approximate scaled derivatives, where              
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate       
!          j-th derivative of y(i), scaled by H**j/factorial(j)         
!          (j = 0,1,...,NQ).  On entry for the first step, the first    
!          two columns of YH must be set from the initial values.       
! NYH    = a constant integer .ge. N, the first dimension of YH.        
! YH1    = a one-dimensional array occupying the same space as YH.      
! EWT    = an array of length N containing multiplicative weights       
!          for local error measurements.  Local errors in y(i) are      
!          compared to 1.0/EWT(i) in various error tests.               
! SAVF   = an array of working storage, of length N.                    
! ACOR   = a work array of length N, used for the accumulated           
!          corrections.  On a successful return, ACOR(i) contains       
!          the estimated one-step local error in y(i).                  
! WM,IWM = real and integer work arrays associated with matrix          
!          operations in chord iteration (MITER .ne. 0).                
! PJAC   = name of routine to evaluate and preprocess Jacobian matrix   
!          and P = I - H*EL0*Jac, if a chord method is being used.      
!          It also returns an estimate of norm(Jac) in PDNORM.          
! SLVS   = name of routine to solve linear system in chord iteration.   
! CCMAX  = maximum relative change in H*EL0 before PJAC is called.      
! H      = the step size to be attempted on the next step.              
!          H is altered by the error control algorithm during the       
!          problem.  H can be either positive or negative, but its      
!          sign must remain constant throughout the problem.            
! HMIN   = the minimum absolute value of the step size H to be used.    
! HMXI   = inverse of the maximum absolute value of H to be used.       
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.   
!          HMIN and HMXI may be changed at any time, but will not       
!          take effect until the next change of H is considered.        
! TN     = the independent variable. TN is updated on each step taken.  
! JSTART = an integer used for input only, with the following           
!          values and meanings:                                         
!               0  perform the first step.                              
!           .gt.0  take a new step continuing from the last.            
!              -1  take the next step with a new value of H,            
!                    N, METH, MITER, and/or matrix parameters.          
!              -2  take the next step with a new value of H,            
!                    but with other inputs unchanged.                   
!          On return, JSTART is set to 1 to facilitate continuation.    
! KFLAG  = a completion code with the following meanings:               
!               0  the step was succesful.                              
!              -1  the requested error could not be achieved.           
!              -2  corrector convergence could not be achieved.         
!              -3  fatal error in PJAC or SLVS.                         
!          A return with KFLAG = -1 or -2 means either                  
!          ABS(H) = HMIN or 10 consecutive failures occurred.           
!          On a return with KFLAG negative, the values of TN and        
!          the YH array are as of the beginning of the last             
!          step, and H is the last step size attempted.                 
! MAXORD = the maximum order of integration method to be allowed.       
! MAXCOR = the maximum number of corrector iterations allowed.          
! MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).   
! MXNCF  = maximum number of convergence failures allowed.              
! METH   = current method.                                              
!          METH = 1 means Adams method (nonstiff)                       
!          METH = 2 means BDF method (stiff)                            
!          METH may be reset by DSTODA.                                 
! MITER  = corrector iteration method.                                  
!          MITER = 0 means functional iteration.                        
!          MITER = JT .gt. 0 means a chord iteration corresponding      
!          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is     
!          communicated here as JTYP, but is not used in DSTODA         
!          except to load MITER following a method switch.)             
!          MITER may be reset by DSTODA.                                
! N      = the number of first-order differential equations.            
!-----------------------------------------------------------------------
!     Common block pointers                                             
      type(DLS001_type), pointer :: DLS001 
      type(DLSA01_type), pointer :: DLSA01 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      DLS001 => common_data%DLS001 
                                                                        
      tmp_ptr = c_loc(DLS001%reals(1)) 
      call c_f_pointer(tmp_ptr,CONIT) 
      tmp_ptr = c_loc(DLS001%reals(2)) 
      call c_f_pointer(tmp_ptr,CRATE) 
      tmp_ptr = c_loc(DLS001%reals(3)) 
      call c_f_pointer(tmp_ptr,EL,[13]) 
      tmp_ptr = c_loc(DLS001%reals(16)) 
      call c_f_pointer(tmp_ptr,ELCO,[13,12]) 
      tmp_ptr = c_loc(DLS001%reals(172)) 
      call c_f_pointer(tmp_ptr,HOLD) 
      tmp_ptr = c_loc(DLS001%reals(173)) 
      call c_f_pointer(tmp_ptr,RMAX) 
      tmp_ptr = c_loc(DLS001%reals(174)) 
      call c_f_pointer(tmp_ptr,TESCO,[3,12]) 
      tmp_ptr = c_loc(DLS001%reals(210)) 
      call c_f_pointer(tmp_ptr,CCMAX) 
      tmp_ptr = c_loc(DLS001%reals(211)) 
      call c_f_pointer(tmp_ptr,EL0) 
      tmp_ptr = c_loc(DLS001%reals(212)) 
      call c_f_pointer(tmp_ptr,H) 
      tmp_ptr = c_loc(DLS001%reals(213)) 
      call c_f_pointer(tmp_ptr,HMIN) 
      tmp_ptr = c_loc(DLS001%reals(214)) 
      call c_f_pointer(tmp_ptr,HMXI) 
      tmp_ptr = c_loc(DLS001%reals(215)) 
      call c_f_pointer(tmp_ptr,HU) 
      tmp_ptr = c_loc(DLS001%reals(216)) 
      call c_f_pointer(tmp_ptr,RC) 
      tmp_ptr = c_loc(DLS001%reals(217)) 
      call c_f_pointer(tmp_ptr,TN) 
      tmp_ptr = c_loc(DLS001%reals(218)) 
      call c_f_pointer(tmp_ptr,UROUND) 
                                                                        
      tmp_ptr = c_loc(DLS001%ints(1)) 
      call c_f_pointer(tmp_ptr,IOWND,[6]) 
      tmp_ptr = c_loc(DLS001%ints(7)) 
      call c_f_pointer(tmp_ptr,IALTH) 
      tmp_ptr = c_loc(DLS001%ints(8)) 
      call c_f_pointer(tmp_ptr,IPUP) 
      tmp_ptr = c_loc(DLS001%ints(9)) 
      call c_f_pointer(tmp_ptr,LMAX) 
      tmp_ptr = c_loc(DLS001%ints(10)) 
      call c_f_pointer(tmp_ptr,MEO) 
      tmp_ptr = c_loc(DLS001%ints(11)) 
      call c_f_pointer(tmp_ptr,NQNYH) 
      tmp_ptr = c_loc(DLS001%ints(12)) 
      call c_f_pointer(tmp_ptr,NSLP) 
      tmp_ptr = c_loc(DLS001%ints(13)) 
      call c_f_pointer(tmp_ptr,ICF) 
      tmp_ptr = c_loc(DLS001%ints(14)) 
      call c_f_pointer(tmp_ptr,IERPJ) 
      tmp_ptr = c_loc(DLS001%ints(15)) 
      call c_f_pointer(tmp_ptr,IERSL) 
      tmp_ptr = c_loc(DLS001%ints(16)) 
      call c_f_pointer(tmp_ptr,JCUR) 
      tmp_ptr = c_loc(DLS001%ints(17)) 
      call c_f_pointer(tmp_ptr,JSTART) 
      tmp_ptr = c_loc(DLS001%ints(18)) 
      call c_f_pointer(tmp_ptr,KFLAG) 
      tmp_ptr = c_loc(DLS001%ints(19)) 
      call c_f_pointer(tmp_ptr,L) 
      tmp_ptr = c_loc(DLS001%ints(20)) 
      call c_f_pointer(tmp_ptr,LYH) 
      tmp_ptr = c_loc(DLS001%ints(21)) 
      call c_f_pointer(tmp_ptr,LEWT) 
      tmp_ptr = c_loc(DLS001%ints(22)) 
      call c_f_pointer(tmp_ptr,LACOR) 
      tmp_ptr = c_loc(DLS001%ints(23)) 
      call c_f_pointer(tmp_ptr,LSAVF) 
      tmp_ptr = c_loc(DLS001%ints(24)) 
      call c_f_pointer(tmp_ptr,LWM) 
      tmp_ptr = c_loc(DLS001%ints(25)) 
      call c_f_pointer(tmp_ptr,LIWM) 
      tmp_ptr = c_loc(DLS001%ints(26)) 
      call c_f_pointer(tmp_ptr,METH) 
      tmp_ptr = c_loc(DLS001%ints(27)) 
      call c_f_pointer(tmp_ptr,MITER) 
      tmp_ptr = c_loc(DLS001%ints(28)) 
      call c_f_pointer(tmp_ptr,MAXORD) 
      tmp_ptr = c_loc(DLS001%ints(29)) 
      call c_f_pointer(tmp_ptr,MAXCOR) 
      tmp_ptr = c_loc(DLS001%ints(30)) 
      call c_f_pointer(tmp_ptr,MSBP) 
      tmp_ptr = c_loc(DLS001%ints(31)) 
      call c_f_pointer(tmp_ptr,MXNCF) 
      tmp_ptr = c_loc(DLS001%ints(32)) 
      call c_f_pointer(tmp_ptr,N) 
      tmp_ptr = c_loc(DLS001%ints(33)) 
      call c_f_pointer(tmp_ptr,NQ) 
      tmp_ptr = c_loc(DLS001%ints(34)) 
      call c_f_pointer(tmp_ptr,NST) 
      tmp_ptr = c_loc(DLS001%ints(35)) 
      call c_f_pointer(tmp_ptr,NFE) 
      tmp_ptr = c_loc(DLS001%ints(36)) 
      call c_f_pointer(tmp_ptr,NJE) 
      tmp_ptr = c_loc(DLS001%ints(37)) 
      call c_f_pointer(tmp_ptr,NQU) 
                                                                        
      DLSA01 => common_data%DLSA01 
                                                                        
      tmp_ptr = c_loc(DLSA01%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWND2) 
      tmp_ptr = c_loc(DLSA01%reals(2)) 
      call c_f_pointer(tmp_ptr,CM1,[12]) 
                                                                        
      tmp_ptr = c_loc(DLSA01%reals(14)) 
      call c_f_pointer(tmp_ptr,CM2,[5]) 
      tmp_ptr = c_loc(DLSA01%reals(19)) 
      call c_f_pointer(tmp_ptr,PDEST) 
      tmp_ptr = c_loc(DLSA01%reals(20)) 
      call c_f_pointer(tmp_ptr,PDLAST) 
      tmp_ptr = c_loc(DLSA01%reals(21)) 
      call c_f_pointer(tmp_ptr,RATIO) 
      tmp_ptr = c_loc(DLSA01%reals(22)) 
      call c_f_pointer(tmp_ptr,PDNORM) 
                                                                        
      tmp_ptr = c_loc(DLSA01%ints(1)) 
      call c_f_pointer(tmp_ptr,IOWND2,[3]) 
      tmp_ptr = c_loc(DLSA01%ints(4)) 
      call c_f_pointer(tmp_ptr,ICOUNT) 
      tmp_ptr = c_loc(DLSA01%ints(5)) 
      call c_f_pointer(tmp_ptr,IRFLAG) 
      tmp_ptr = c_loc(DLSA01%ints(6)) 
      call c_f_pointer(tmp_ptr,JTYP) 
      tmp_ptr = c_loc(DLSA01%ints(7)) 
      call c_f_pointer(tmp_ptr,MUSED) 
      tmp_ptr = c_loc(DLSA01%ints(8)) 
      call c_f_pointer(tmp_ptr,MXORDN) 
      tmp_ptr = c_loc(DLSA01%ints(9)) 
      call c_f_pointer(tmp_ptr,MXORDS) 
                                                                        
!                                                                       
      KFLAG = 0 
      TOLD = TN 
      NCF = 0 
      IERPJ = 0 
      IERSL = 0 
      JCUR = 0 
      ICF = 0 
      DELP = 0.0D0 
      IF (JSTART .GT. 0) GO TO 200 
      IF (JSTART .EQ. -1) GO TO 100 
      IF (JSTART .EQ. -2) GO TO 160 
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are     
! initialized.  RMAX is the maximum ratio by which H can be increased   
! in a single step.  It is initially 1.E4 to compensate for the small   
! initial H, but then is normally equal to 10.  If a failure            
! occurs (in corrector convergence or error test), RMAX is set at 2     
! for the next increase.                                                
! DCFODE is called to get the needed coefficients for both methods.     
!-----------------------------------------------------------------------
      LMAX = MAXORD + 1 
      NQ = 1 
      L = 2 
      IALTH = 2 
      RMAX = 10000.0D0 
      RC = 0.0D0 
      EL0 = 1.0D0 
      CRATE = 0.7D0 
      HOLD = H 
      NSLP = 0 
      IPUP = MITER 
      IRET = 3 
! Initialize switching parameters.  METH = 1 is assumed initially. -----
      ICOUNT = 20 
      IRFLAG = 0 
      PDEST = 0.0D0 
      PDLAST = 0.0D0 
      RATIO = 5.0D0 
      CALL DCFODE (2, ELCO, TESCO) 
      DO I = 1,5 
        CM2(I) = TESCO(2,I)*ELCO(I+1,I)
      END DO
      CALL DCFODE (1, ELCO, TESCO) 
      DO I = 1,12 
        CM1(I) = TESCO(2,I)*ELCO(I+1,I)
      END DO
      GO TO 150 
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.    
! IPUP is set to MITER to force a matrix update.                        
! If an order increase is about to be considered (IALTH = 1),           
! IALTH is reset to 2 to postpone consideration one more step.          
! If the caller has changed METH, DCFODE is called to reset             
! the coefficients of the method.                                       
! If H is to be changed, YH must be rescaled.                           
! If H or METH is being changed, IALTH is reset to L = NQ + 1           
! to prevent further changes in H for that many steps.                  
!-----------------------------------------------------------------------
  100 IPUP = MITER 
      LMAX = MAXORD + 1 
      IF (IALTH .EQ. 1) IALTH = 2 
      IF (METH .EQ. MUSED) GO TO 160 
      CALL DCFODE (METH, ELCO, TESCO) 
      IALTH = L 
      IRET = 1 
!-----------------------------------------------------------------------
! The el vector and related constants are reset                         
! whenever the order NQ is changed, or at the start of the problem.     
!-----------------------------------------------------------------------
  150 DO I = 1,L 
        EL(I) = ELCO(I,NQ)
      END DO
      NQNYH = NQ*NYH 
      RC = RC*EL(1)/EL0 
      EL0 = EL(1) 
      CONIT = 0.5D0/(NQ+2) 
      GO TO (160, 170, 200), IRET 
!-----------------------------------------------------------------------
! If H is being changed, the H ratio RH is checked against              
! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to     
! L = NQ + 1 to prevent a change of H for that many steps, unless       
! forced by a convergence or error test failure.                        
!-----------------------------------------------------------------------
  160 IF (H .EQ. HOLD) GO TO 200 
      RH = H/HOLD 
      H = HOLD 
      IREDO = 3 
      GO TO 175 
  170 RH = MAX(RH,HMIN/ABS(H)) 
  175 RH = MIN(RH,RMAX) 
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH) 
!-----------------------------------------------------------------------
! If METH = 1, also restrict the new step size by the stability region. 
! If this reduces H, set IRFLAG to 1 so that if there are roundoff      
! problems later, we can assume that is the cause of the trouble.       
!-----------------------------------------------------------------------
      IF (METH .EQ. 2) GO TO 178 
      IRFLAG = 0 
      PDH = MAX(ABS(H)*PDLAST,0.000001D0) 
      IF (RH*PDH*1.00001D0 .LT. SM1(NQ)) GO TO 178 
      RH = SM1(NQ)/PDH 
      IRFLAG = 1 
  178 CONTINUE 
      R = 1.0D0 
      DO J = 2,L 
        R = R*RH 
        DO I = 1,N 
          YH(I,J) = YH(I,J)*R
        END DO
      END DO
      H = H*RH 
      RC = RC*RH 
      IALTH = L 
      IF (IREDO .EQ. 0) GO TO 690 
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively             
! multiplying the YH array by the Pascal triangle matrix.               
! RC is the ratio of new to old values of the coefficient  H*EL(1).     
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER       
! to force PJAC to be called, if a Jacobian is involved.                
! In any case, PJAC is called at least every MSBP steps.                
!-----------------------------------------------------------------------
  200 IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER 
      IF (NST .GE. NSLP+MSBP) IPUP = MITER 
      TN = TN + H 
      I1 = NQNYH + 1 
      DO 215 JB = 1,NQ 
        I1 = I1 - NYH 
!DIR$ IVDEP                                                             
        DO I = I1,NQNYH 
         YH1(I) = YH1(I) + YH1(I+NYH)
        END DO
  215   CONTINUE 
      PNORM = DMNORM (N, YH1, EWT) 
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is   
! made on the RMS-norm of each correction, weighted by the error        
! weight vector EWT.  The sum of the corrections is accumulated in the  
! vector ACOR(i).  The YH array is not altered in the corrector loop.   
!-----------------------------------------------------------------------
  220 M = 0 
      RATE = 0.0D0 
      DEL = 0.0D0 
      DO I = 1,N 
        Y(I) = YH(I,1)
      END DO
      CALL F (NEQ, TN, Y, SAVF, common_data%ierr) 
      if (common_data%ierr < 0) return 
      NFE = NFE + 1 
      IF (IPUP .LE. 0) GO TO 250 
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - H*EL(1)*J is reevaluated and         
! preprocessed before starting the corrector iteration.  IPUP is set    
! to 0 as an indicator that this has been done.                         
!-----------------------------------------------------------------------
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC,     &
     &           common_data)                                           
      if (common_data%ierr < 0) return 
      IPUP = 0 
      RC = 1.0D0 
      NSLP = NST 
      CRATE = 0.7D0 
      IF (IERPJ .NE. 0) GO TO 430 
  250 DO I = 1,N 
        ACOR(I) = 0.0D0
      END DO
  270 IF (MITER .NE. 0) GO TO 350 
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from           
! the result of the last function evaluation.                           
!-----------------------------------------------------------------------
      DO I = 1,N 
        SAVF(I) = H*SAVF(I) - YH(I,2) 
        Y(I) = SAVF(I) - ACOR(I)
      END DO
      DEL = DMNORM (N, Y, EWT) 
      DO I = 1,N 
        Y(I) = YH(I,1) + EL(1)*SAVF(I) 
        ACOR(I) = SAVF(I)
      END DO
      GO TO 400 
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,         
! and solve the linear system with that as right-hand side and          
! P as coefficient matrix.                                              
!-----------------------------------------------------------------------
  350 DO I = 1,N 
        Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      END DO
      CALL SLVS (WM, IWM, Y, SAVF, common_data) 
      IF (IERSL .LT. 0) GO TO 430 
      IF (IERSL .GT. 0) GO TO 410 
      DEL = DMNORM (N, Y, EWT) 
      DO I = 1,N 
        ACOR(I) = ACOR(I) + Y(I) 
        Y(I) = YH(I,1) + EL(1)*ACOR(I)
      END DO
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence    
! rate constant is stored in CRATE, and this is used in the test.       
!                                                                       
! We first check for a change of iterates that is the size of           
! roundoff error.  If this occurs, the iteration has converged, and a   
! new rate estimate is not formed.                                      
! In all other cases, force at least two iterations to estimate a       
! local Lipschitz constant estimate for Adams methods.                  
! On convergence, form PDEST = local maximum Lipschitz constant         
! estimate.  PDLAST is the most recent nonzero estimate.                
!-----------------------------------------------------------------------
  400 CONTINUE 
      IF (DEL .LE. 100.0D0*PNORM*UROUND) GO TO 450 
      IF (M .EQ. 0 .AND. METH .EQ. 1) GO TO 405 
      IF (M .EQ. 0) GO TO 402 
      RM = 1024.0D0 
      IF (DEL .LE. 1024.0D0*DELP) RM = DEL/DELP 
      RATE = MAX(RATE,RM) 
      CRATE = MAX(0.2D0*CRATE,RM) 
  402 DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT) 
      IF (DCON .GT. 1.0D0) GO TO 405 
      PDEST = MAX(PDEST,RATE/ABS(H*EL(1))) 
      IF (PDEST .NE. 0.0D0) PDLAST = PDEST 
      GO TO 450 
  405 CONTINUE 
      M = M + 1 
      IF (M .EQ. MAXCOR) GO TO 410 
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410 
      DELP = DEL 
      CALL F (NEQ, TN, Y, SAVF, common_data%ierr) 
      if (common_data%ierr < 0) return 
      NFE = NFE + 1 
      GO TO 270 
!-----------------------------------------------------------------------
! The corrector iteration failed to converge.                           
! If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for   
! the next try.  Otherwise the YH array is retracted to its values      
! before prediction, and H is reduced, if possible.  If H cannot be     
! reduced or MXNCF failures have occurred, exit with KFLAG = -2.        
!-----------------------------------------------------------------------
  410 IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430 
      ICF = 1 
      IPUP = MITER 
      GO TO 220 
  430 ICF = 2 
      NCF = NCF + 1 
      RMAX = 2.0D0 
      TN = TOLD 
      I1 = NQNYH + 1 
      DO 445 JB = 1,NQ 
        I1 = I1 - NYH 
!DIR$ IVDEP                                                             
        DO I = I1,NQNYH 
          YH1(I) = YH1(I) - YH1(I+NYH)
        END DO
  445   CONTINUE 
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680 
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670 
      IF (NCF .EQ. MXNCF) GO TO 670 
      RH = 0.25D0 
      IPUP = MITER 
      IREDO = 1 
      GO TO 170 
!-----------------------------------------------------------------------
! The corrector has converged.  JCUR is set to 0                        
! to signal that the Jacobian involved may need updating later.         
! The local error test is made and control passes to statement 500      
! if it fails.                                                          
!-----------------------------------------------------------------------
  450 JCUR = 0 
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ) 
      IF (M .GT. 0) DSM = DMNORM (N, ACOR, EWT)/TESCO(2,NQ) 
      IF (DSM .GT. 1.0D0) GO TO 500 
!-----------------------------------------------------------------------
! After a successful step, update the YH array.                         
! Decrease ICOUNT by 1, and if it is -1, consider switching methods.    
! If a method switch is made, reset various parameters,                 
! rescale the YH array, and exit.  If there is no switch,               
! consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.     
! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for         
! use in a possible order increase on the next step.                    
! If a change in H is considered, an increase or decrease in order      
! by one is considered also.  A change in H is made only if it is by a  
! factor of at least 1.1.  If not, IALTH is set to 3 to prevent         
! testing for that many steps.                                          
!-----------------------------------------------------------------------
      KFLAG = 0 
      IREDO = 0 
      NST = NST + 1 
      HU = H 
      NQU = NQ 
      MUSED = METH 
      DO J = 1,L 
        DO I = 1,N 
          YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
        END DO
      END DO
      ICOUNT = ICOUNT - 1 
      IF (ICOUNT .GE. 0) GO TO 488 
      IF (METH .EQ. 2) GO TO 480 
!-----------------------------------------------------------------------
! We are currently using an Adams method.  Consider switching to BDF.   
! If the current order is greater than 5, assume the problem is         
! not stiff, and skip this section.                                     
! If the Lipschitz constant and error estimate are not polluted         
! by roundoff, go to 470 and perform the usual test.                    
! Otherwise, switch to the BDF methods if the last step was             
! restricted to insure stability (irflag = 1), and stay with Adams      
! method if not.  When switching to BDF with polluted error estimates,  
! in the absence of other information, double the step size.            
!                                                                       
! When the estimates are OK, we make the usual test by computing        
! the step size we could have (ideally) used on this step,              
! with the current (Adams) method, and also that for the BDF.           
! If NQ .gt. MXORDS, we consider changing to order MXORDS on switching. 
! Compare the two step sizes to decide whether to switch.               
! The step size advantage must be at least RATIO = 5 to switch.         
!-----------------------------------------------------------------------
      IF (NQ .GT. 5) GO TO 488 
      IF (DSM .GT. 100.0D0*PNORM*UROUND .AND. PDEST .NE. 0.0D0)         &
     &   GO TO 470                                                      
      IF (IRFLAG .EQ. 0) GO TO 488 
      RH2 = 2.0D0 
      NQM2 = MIN(NQ,MXORDS) 
      GO TO 478 
  470 CONTINUE 
      EXSM = 1.0D0/L 
      RH1 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0) 
      RH1IT = 2.0D0*RH1 
      PDH = PDLAST*ABS(H) 
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQ)/PDH 
      RH1 = MIN(RH1,RH1IT) 
      IF (NQ .LE. MXORDS) GO TO 474 
         NQM2 = MXORDS 
         LM2 = MXORDS + 1 
         EXM2 = 1.0D0/LM2 
         LM2P1 = LM2 + 1 
         DM2 = DMNORM (N, YH(1,LM2P1), EWT)/CM2(MXORDS) 
         RH2 = 1.0D0/(1.2D0*DM2**EXM2 + 0.0000012D0) 
         GO TO 476 
  474 DM2 = DSM*(CM1(NQ)/CM2(NQ)) 
      RH2 = 1.0D0/(1.2D0*DM2**EXSM + 0.0000012D0) 
      NQM2 = NQ 
  476 CONTINUE 
      IF (RH2 .LT. RATIO*RH1) GO TO 488 
! THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
  478 RH = RH2 
      ICOUNT = 20 
      METH = 2 
      MITER = JTYP 
      PDLAST = 0.0D0 
      NQ = NQM2 
      L = NQ + 1 
      GO TO 170 
!-----------------------------------------------------------------------
! We are currently using a BDF method.  Consider switching to Adams.    
! Compute the step size we could have (ideally) used on this step,      
! with the current (BDF) method, and also that for the Adams.           
! If NQ .gt. MXORDN, we consider changing to order MXORDN on switching. 
! Compare the two step sizes to decide whether to switch.               
! The step size advantage must be at least 5/RATIO = 1 to switch.       
! If the step size for Adams would be so small as to cause              
! roundoff pollution, we stay with BDF.                                 
!-----------------------------------------------------------------------
  480 CONTINUE 
      EXSM = 1.0D0/L 
      IF (MXORDN .GE. NQ) GO TO 484 
         NQM1 = MXORDN 
         LM1 = MXORDN + 1 
         EXM1 = 1.0D0/LM1 
         LM1P1 = LM1 + 1 
         DM1 = DMNORM (N, YH(1,LM1P1), EWT)/CM1(MXORDN) 
         RH1 = 1.0D0/(1.2D0*DM1**EXM1 + 0.0000012D0) 
         GO TO 486 
  484 DM1 = DSM*(CM2(NQ)/CM1(NQ)) 
      RH1 = 1.0D0/(1.2D0*DM1**EXSM + 0.0000012D0) 
      NQM1 = NQ 
      EXM1 = EXSM 
  486 RH1IT = 2.0D0*RH1 
      PDH = PDNORM*ABS(H) 
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQM1)/PDH 
      RH1 = MIN(RH1,RH1IT) 
      RH2 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0) 
      IF (RH1*RATIO .LT. 5.0D0*RH2) GO TO 488 
      ALPHA = MAX(0.001D0,RH1) 
      DM1 = (ALPHA**EXM1)*DM1 
      IF (DM1 .LE. 1000.0D0*UROUND*PNORM) GO TO 488 
! The switch test passed.  Reset relevant quantities for Adams. --------
      RH = RH1 
      ICOUNT = 20 
      METH = 1 
      MITER = 0 
      PDLAST = 0.0D0 
      NQ = NQM1 
      L = NQ + 1 
      GO TO 170 
!                                                                       
! No method switch is being made.  Do the usual step/order selection. --
  488 CONTINUE 
      IALTH = IALTH - 1 
      IF (IALTH .EQ. 0) GO TO 520 
      IF (IALTH .GT. 1) GO TO 700 
      IF (L .EQ. LMAX) GO TO 700 
      DO I = 1,N 
        YH(I,LMAX) = ACOR(I)
      END DO
      GO TO 700 
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.       
! Restore TN and the YH array to their previous values, and prepare     
! to try the step again.  Compute the optimum step size for this or     
! one lower order.  After 2 or more failures, H is forced to decrease   
! by a factor of 0.2 or less.                                           
!-----------------------------------------------------------------------
  500 KFLAG = KFLAG - 1 
      TN = TOLD 
      I1 = NQNYH + 1 
      DO 515 JB = 1,NQ 
        I1 = I1 - NYH 
!DIR$ IVDEP                                                             
        DO I = I1,NQNYH 
          YH1(I) = YH1(I) - YH1(I+NYH)
        END DO
  515   CONTINUE 
      RMAX = 2.0D0 
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660 
      IF (KFLAG .LE. -3) GO TO 640 
      IREDO = 2 
      RHUP = 0.0D0 
      GO TO 540 
!-----------------------------------------------------------------------
! Regardless of the success or failure of the step, factors             
! RHDN, RHSM, and RHUP are computed, by which H could be multiplied     
! at order NQ - 1, order NQ, or order NQ + 1, respectively.             
! In the case of failure, RHUP = 0.0 to avoid an order increase.        
! The largest of these is determined and the new order chosen           
! accordingly.  If the order is to be increased, we compute one         
! additional scaled derivative.                                         
!-----------------------------------------------------------------------
  520 RHUP = 0.0D0 
      IF (L .EQ. LMAX) GO TO 540 
      DO I = 1,N 
        SAVF(I) = ACOR(I) - YH(I,LMAX)
      END DO
      DUP = DMNORM (N, SAVF, EWT)/TESCO(3,NQ) 
      EXUP = 1.0D0/(L+1) 
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0) 
  540 EXSM = 1.0D0/L 
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0) 
      RHDN = 0.0D0 
      IF (NQ .EQ. 1) GO TO 550 
      DDN = DMNORM (N, YH(1,L), EWT)/TESCO(1,NQ) 
      EXDN = 1.0D0/NQ 
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0) 
! If METH = 1, limit RH according to the stability region also. --------
  550 IF (METH .EQ. 2) GO TO 560 
      PDH = MAX(ABS(H)*PDLAST,0.000001D0) 
      IF (L .LT. LMAX) RHUP = MIN(RHUP,SM1(L)/PDH) 
      RHSM = MIN(RHSM,SM1(NQ)/PDH) 
      IF (NQ .GT. 1) RHDN = MIN(RHDN,SM1(NQ-1)/PDH) 
      PDEST = 0.0D0 
  560 IF (RHSM .GE. RHUP) GO TO 570 
      IF (RHUP .GT. RHDN) GO TO 590 
      GO TO 580 
  570 IF (RHSM .LT. RHDN) GO TO 580 
      NEWQ = NQ 
      RH = RHSM 
      GO TO 620 
  580 NEWQ = NQ - 1 
      RH = RHDN 
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0 
      GO TO 620 
  590 NEWQ = L 
      RH = RHUP 
      IF (RH .LT. 1.1D0) GO TO 610 
      R = EL(L)/L 
      DO I = 1,N 
        YH(I,NEWQ+1) = ACOR(I)*R
      END DO
      GO TO 630 
  610 IALTH = 3 
      GO TO 700 
! If METH = 1 and H is restricted by stability, bypass 10 percent test. 
  620 IF (METH .EQ. 2) GO TO 622 
      IF (RH*PDH*1.00001D0 .GE. SM1(NEWQ)) GO TO 625 
  622 IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GO TO 610 
  625 IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0) 
!-----------------------------------------------------------------------
! If there is a change of order, reset NQ, L, and the coefficients.     
! In any case H is reset according to RH and the YH array is rescaled.  
! Then exit from 690 if the step was OK, or redo the step otherwise.    
!-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170 
  630 NQ = NEWQ 
      L = NQ + 1 
      IRET = 2 
      GO TO 150 
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more failures have occured.      
! If 10 failures have occurred, exit with KFLAG = -1.                   
! It is assumed that the derivatives that have accumulated in the       
! YH array have errors of the wrong order.  Hence the first             
! derivative is recomputed, and the order is set to 1.  Then            
! H is reduced by a factor of 10, and the step is retried,              
! until it succeeds or H reaches HMIN.                                  
!-----------------------------------------------------------------------
  640 IF (KFLAG .EQ. -10) GO TO 660 
      RH = 0.1D0 
      RH = MAX(HMIN/ABS(H),RH) 
      H = H*RH 
      DO I = 1,N 
        Y(I) = YH(I,1)
      END DO
      CALL F (NEQ, TN, Y, SAVF, common_data%ierr) 
      if (common_data%ierr < 0) return 
      NFE = NFE + 1 
      DO I = 1,N 
        YH(I,2) = H*SAVF(I)
      END DO
      IPUP = MITER 
      IALTH = 5 
      IF (NQ .EQ. 1) GO TO 200 
      NQ = 1 
      L = 2 
      IRET = 3 
      GO TO 150 
!-----------------------------------------------------------------------
! All returns are made through this section.  H is saved in HOLD        
! to allow the caller to change H on the next step.                     
!-----------------------------------------------------------------------
  660 KFLAG = -1 
      GO TO 720 
  670 KFLAG = -2 
      GO TO 720 
  680 KFLAG = -3 
      GO TO 720 
  690 RMAX = 10.0D0 
  700 R = 1.0D0/TESCO(2,NQU) 
      DO I = 1,N 
        ACOR(I) = ACOR(I)*R
      END DO
  720 HOLD = H 
      JSTART = 1 
      RETURN 
!----------------------- End of Subroutine DSTODA ----------------------
      END                                           
!DECK DPRJA                                                             
      SUBROUTINE DPRJA (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM,      &
     &   F, JAC, common_data)                                           
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_common 
      type(odepack_common_data), target, intent(inout) :: common_data 
      EXTERNAL F, JAC 
      INTEGER NEQ, NYH, IWM 
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM 
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),      &
     &   WM(*), IWM(*)                                                  
      INTEGER, pointer :: IOWND(:), IOWNS(:),                           &
     &   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                     &
     &   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,               &
     &   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU         
      INTEGER, pointer ::                                               &
     &   IOWND2(:), IOWNS2(:), JTYP, MUSED, MXORDN, MXORDS              
      DOUBLE PRECISION, pointer :: ROWNS(:),                            &
     &   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND                  
      DOUBLE PRECISION, pointer :: ROWND2, ROWNS2(:), PDNORM 
!      COMMON /DLS001/ ROWNS(209),                                      
!     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,                
!     2   IOWND(6), IOWNS(6),                                           
!     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                    
!     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,              
!     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU        
!      COMMON /DLSA01/ ROWND2, ROWNS2(20), PDNORM,                      
!     1   IOWND2(3), IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS             
      INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP,                      &
     &   MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1                     
      DOUBLE PRECISION CON, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ,         &
     &   DMNORM, DFNORM, DBNORM                                         
!-----------------------------------------------------------------------
! DPRJA is called by DSTODA to compute and process the matrix           
! P = I - H*EL(1)*J , where J is an approximation to the Jacobian.      
! Here J is computed by the user-supplied routine JAC if                
! MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.           
! J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the      
! matrix norm consistent with the weighted max-norm on vectors given    
! by DMNORM) is computed, and J is overwritten by P.  P is then         
! subjected to LU decomposition in preparation for later solution       
! of linear systems with P as coefficient matrix.  This is done         
! by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.           
!                                                                       
! In addition to variables described previously, communication          
! with DPRJA uses the following:                                        
! Y     = array containing predicted values on entry.                   
! FTEM  = work array of length N (ACOR in DSTODA).                      
! SAVF  = array containing f evaluated at predicted y.                  
! WM    = real work space for matrices.  On output it contains the      
!         LU decomposition of P.                                        
!         Storage of matrix elements starts at WM(3).                   
!         WM also contains the following matrix-related data:           
!         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.  
! IWM   = integer work space containing pivot information, starting at  
!         IWM(21).   IWM also contains the band parameters              
!         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.               
! EL0   = EL(1) (input).                                                
! PDNORM= norm of Jacobian matrix. (Output).                            
! IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if              
!         P matrix found to be singular.                                
! JCUR  = output flag = 1 to indicate that the Jacobian matrix          
!         (or approximation) is now current.                            
! This routine also uses the Common variables EL0, H, TN, UROUND,       
! MITER, N, NFE, and NJE.                                               
!-----------------------------------------------------------------------
!     Common block pointers                                             
      type(DLS001_type), pointer :: DLS001 
      type(DLSA01_type), pointer :: DLSA01 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      DLS001 => common_data%DLS001 
                                                                        
      tmp_ptr = c_loc(DLS001%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWNS,[209]) 
      tmp_ptr = c_loc(DLS001%reals(210)) 
      call c_f_pointer(tmp_ptr,CCMAX) 
      tmp_ptr = c_loc(DLS001%reals(211)) 
      call c_f_pointer(tmp_ptr,EL0) 
      tmp_ptr = c_loc(DLS001%reals(212)) 
      call c_f_pointer(tmp_ptr,H) 
      tmp_ptr = c_loc(DLS001%reals(213)) 
      call c_f_pointer(tmp_ptr,HMIN) 
      tmp_ptr = c_loc(DLS001%reals(214)) 
      call c_f_pointer(tmp_ptr,HMXI) 
      tmp_ptr = c_loc(DLS001%reals(215)) 
      call c_f_pointer(tmp_ptr,HU) 
      tmp_ptr = c_loc(DLS001%reals(216)) 
      call c_f_pointer(tmp_ptr,RC) 
      tmp_ptr = c_loc(DLS001%reals(217)) 
      call c_f_pointer(tmp_ptr,TN) 
      tmp_ptr = c_loc(DLS001%reals(218)) 
      call c_f_pointer(tmp_ptr,UROUND) 
                                                                        
      tmp_ptr = c_loc(DLS001%ints(1)) 
      call c_f_pointer(tmp_ptr,IOWND,[6]) 
      tmp_ptr = c_loc(DLS001%ints(7)) 
      call c_f_pointer(tmp_ptr,IOWNS,[6]) 
      tmp_ptr = c_loc(DLS001%ints(13)) 
      call c_f_pointer(tmp_ptr,ICF) 
      tmp_ptr = c_loc(DLS001%ints(14)) 
      call c_f_pointer(tmp_ptr,IERPJ) 
      tmp_ptr = c_loc(DLS001%ints(15)) 
      call c_f_pointer(tmp_ptr,IERSL) 
      tmp_ptr = c_loc(DLS001%ints(16)) 
      call c_f_pointer(tmp_ptr,JCUR) 
      tmp_ptr = c_loc(DLS001%ints(17)) 
      call c_f_pointer(tmp_ptr,JSTART) 
      tmp_ptr = c_loc(DLS001%ints(18)) 
      call c_f_pointer(tmp_ptr,KFLAG) 
      tmp_ptr = c_loc(DLS001%ints(19)) 
      call c_f_pointer(tmp_ptr,L) 
      tmp_ptr = c_loc(DLS001%ints(20)) 
      call c_f_pointer(tmp_ptr,LYH) 
      tmp_ptr = c_loc(DLS001%ints(21)) 
      call c_f_pointer(tmp_ptr,LEWT) 
      tmp_ptr = c_loc(DLS001%ints(22)) 
      call c_f_pointer(tmp_ptr,LACOR) 
      tmp_ptr = c_loc(DLS001%ints(23)) 
      call c_f_pointer(tmp_ptr,LSAVF) 
      tmp_ptr = c_loc(DLS001%ints(24)) 
      call c_f_pointer(tmp_ptr,LWM) 
      tmp_ptr = c_loc(DLS001%ints(25)) 
      call c_f_pointer(tmp_ptr,LIWM) 
      tmp_ptr = c_loc(DLS001%ints(26)) 
      call c_f_pointer(tmp_ptr,METH) 
      tmp_ptr = c_loc(DLS001%ints(27)) 
      call c_f_pointer(tmp_ptr,MITER) 
      tmp_ptr = c_loc(DLS001%ints(28)) 
      call c_f_pointer(tmp_ptr,MAXORD) 
      tmp_ptr = c_loc(DLS001%ints(29)) 
      call c_f_pointer(tmp_ptr,MAXCOR) 
      tmp_ptr = c_loc(DLS001%ints(30)) 
      call c_f_pointer(tmp_ptr,MSBP) 
      tmp_ptr = c_loc(DLS001%ints(31)) 
      call c_f_pointer(tmp_ptr,MXNCF) 
      tmp_ptr = c_loc(DLS001%ints(32)) 
      call c_f_pointer(tmp_ptr,N) 
      tmp_ptr = c_loc(DLS001%ints(33)) 
      call c_f_pointer(tmp_ptr,NQ) 
      tmp_ptr = c_loc(DLS001%ints(34)) 
      call c_f_pointer(tmp_ptr,NST) 
      tmp_ptr = c_loc(DLS001%ints(35)) 
      call c_f_pointer(tmp_ptr,NFE) 
      tmp_ptr = c_loc(DLS001%ints(36)) 
      call c_f_pointer(tmp_ptr,NJE) 
      tmp_ptr = c_loc(DLS001%ints(37)) 
      call c_f_pointer(tmp_ptr,NQU) 
                                                                        
      DLSA01 => common_data%DLSA01 
                                                                        
      tmp_ptr = c_loc(DLSA01%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWND2) 
      tmp_ptr = c_loc(DLSA01%reals(2)) 
      call c_f_pointer(tmp_ptr,ROWNS2,[20]) 
      tmp_ptr = c_loc(DLSA01%reals(22)) 
      call c_f_pointer(tmp_ptr,PDNORM) 
                                                                        
      tmp_ptr = c_loc(DLSA01%ints(1)) 
      call c_f_pointer(tmp_ptr,IOWND2,[3]) 
      tmp_ptr = c_loc(DLSA01%ints(4)) 
      call c_f_pointer(tmp_ptr,IOWNS2,[2]) 
      tmp_ptr = c_loc(DLSA01%ints(6)) 
      call c_f_pointer(tmp_ptr,JTYP) 
      tmp_ptr = c_loc(DLSA01%ints(7)) 
      call c_f_pointer(tmp_ptr,MUSED) 
      tmp_ptr = c_loc(DLSA01%ints(8)) 
      call c_f_pointer(tmp_ptr,MXORDN) 
      tmp_ptr = c_loc(DLSA01%ints(9)) 
      call c_f_pointer(tmp_ptr,MXORDS) 
!                                                                       
      NJE = NJE + 1 
      IERPJ = 0 
      JCUR = 1 
      HL0 = H*EL0 
      GO TO (100, 200, 300, 400, 500), MITER 
! If MITER = 1, call JAC and multiply by scalar. -----------------------
  100 LENP = N*N 
      DO I = 1,LENP 
        WM(I+2) = 0.0D0
      END DO
      CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N, common_data%ierr) 
      if (common_data%ierr < 0) return 
      CON = -HL0 
      DO I = 1,LENP 
        WM(I+2) = WM(I+2)*CON
      END DO
      GO TO 240 
! If MITER = 2, make N calls to F to approximate J. --------------------
  200 FAC = DMNORM (N, SAVF, EWT) 
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC 
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0 
      SRUR = WM(1) 
      J1 = 2 
      DO 230 J = 1,N 
        YJ = Y(J) 
        R = MAX(SRUR*ABS(YJ),R0/EWT(J)) 
        Y(J) = Y(J) + R 
        FAC = -HL0/R 
        CALL F (NEQ, TN, Y, FTEM, common_data%ierr) 
        if (common_data%ierr < 0) return 
        DO I = 1,N 
          WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        END DO
        Y(J) = YJ 
        J1 = J1 + N 
  230   CONTINUE 
      NFE = NFE + N 
  240 CONTINUE 
! Compute norm of Jacobian. --------------------------------------------
      PDNORM = DFNORM (N, WM(3), EWT)/ABS(HL0) 
! Add identity matrix. -------------------------------------------------
      J = 3 
      NP1 = N + 1 
      DO I = 1,N 
        WM(J) = WM(J) + 1.0D0 
        J = J + NP1
      END DO
! Do LU decomposition on P. --------------------------------------------
!      CALL DGEFA (WM(3), N, N, IWM(21), IER)                           
      call dgetrf (n, n, wm(3), n, iwm(21), ier) 
      IF (IER .NE. 0) IERPJ = 1 
      RETURN 
! Dummy block only, since MITER is never 3 in this routine. ------------
  300 RETURN 
! If MITER = 4, call JAC and multiply by scalar. -----------------------
  400 ML = IWM(1) 
      MU = IWM(2) 
      ML3 = ML + 3 
      MBAND = ML + MU + 1 
      MEBAND = MBAND + ML 
      LENP = MEBAND*N 
      DO I = 1,LENP 
        WM(I+2) = 0.0D0
      END DO
      CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND, common_data%ierr) 
      if (common_data%ierr < 0) return 
      CON = -HL0 
      DO I = 1,LENP 
        WM(I+2) = WM(I+2)*CON
      END DO
      GO TO 570 
! If MITER = 5, make MBAND calls to F to approximate J. ----------------
  500 ML = IWM(1) 
      MU = IWM(2) 
      MBAND = ML + MU + 1 
      MBA = MIN(MBAND,N) 
      MEBAND = MBAND + ML 
      MEB1 = MEBAND - 1 
      SRUR = WM(1) 
      FAC = DMNORM (N, SAVF, EWT) 
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC 
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0 
      DO 560 J = 1,MBA 
        DO I = J,N,MBAND 
          YI = Y(I) 
          R = MAX(SRUR*ABS(YI),R0/EWT(I)) 
          Y(I) = Y(I) + R
        END DO
        CALL F (NEQ, TN, Y, FTEM, common_data%ierr) 
        if (common_data%ierr < 0) return 
        DO 550 JJ = J,N,MBAND 
          Y(JJ) = YH(JJ,1) 
          YJJ = Y(JJ) 
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ)) 
          FAC = -HL0/R 
          I1 = MAX(JJ-MU,1) 
          I2 = MIN(JJ+ML,N) 
          II = JJ*MEB1 - ML + 2 
          DO I = I1,I2 
            WM(II+I) = (FTEM(I) - SAVF(I))*FAC
          END DO
  550     CONTINUE 
  560   CONTINUE 
      NFE = NFE + MBA 
  570 CONTINUE 
! Compute norm of Jacobian. --------------------------------------------
      PDNORM = DBNORM (N, WM(ML+3), MEBAND, ML, MU, EWT)/ABS(HL0) 
! Add identity matrix. -------------------------------------------------
      II = MBAND + 2 
      DO I = 1,N 
        WM(II) = WM(II) + 1.0D0 
        II = II + MEBAND
      END DO
! Do LU decomposition of P. --------------------------------------------
!      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)              
      call dgbtrf (n, n, ml, mu, wm(3), meband, iwm(21), ier) 
      IF (IER .NE. 0) IERPJ = 1 
      RETURN 
!----------------------- End of Subroutine DPRJA -----------------------
      END                                           
!DECK DMNORM                                                            
      DOUBLE PRECISION FUNCTION DMNORM (N, V, W) 
!-----------------------------------------------------------------------
! This function routine computes the weighted max-norm                  
! of the vector of length N contained in the array V, with weights      
! contained in the array w of length N:                                 
!   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)                              
!-----------------------------------------------------------------------
      INTEGER N,   I 
      DOUBLE PRECISION V, W,   VM 
      DIMENSION V(N), W(N) 
      VM = 0.0D0 
      DO I = 1,N 
        VM = MAX(VM,ABS(V(I))*W(I))
      END DO
      DMNORM = VM 
      RETURN 
!----------------------- End of Function DMNORM ------------------------
      END                                           
!DECK DFNORM                                                            
      DOUBLE PRECISION FUNCTION DFNORM (N, A, W) 
!-----------------------------------------------------------------------
! This function computes the norm of a full N by N matrix,              
! stored in the array A, that is consistent with the weighted max-norm  
! on vectors, with weights stored in the array W:                       
!   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )  
!-----------------------------------------------------------------------
      INTEGER N,   I, J 
      DOUBLE PRECISION A,   W, AN, SUM 
      DIMENSION A(N,N), W(N) 
      AN = 0.0D0 
      DO 20 I = 1,N 
        SUM = 0.0D0 
        DO J = 1,N 
          SUM = SUM + ABS(A(I,J))/W(J)
        END DO
        AN = MAX(AN,SUM*W(I)) 
   20   CONTINUE 
      DFNORM = AN 
      RETURN 
!----------------------- End of Function DFNORM ------------------------
      END                                           
!DECK DBNORM                                                            
      DOUBLE PRECISION FUNCTION DBNORM (N, A, NRA, ML, MU, W) 
!-----------------------------------------------------------------------
! This function computes the norm of a banded N by N matrix,            
! stored in the array A, that is consistent with the weighted max-norm  
! on vectors, with weights stored in the array W.                       
! ML and MU are the lower and upper half-bandwidths of the matrix.      
! NRA is the first dimension of the A array, NRA .ge. ML+MU+1.          
! In terms of the matrix elements a(i,j), the norm is given by:         
!   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )  
!-----------------------------------------------------------------------
      INTEGER N, NRA, ML, MU 
      INTEGER I, I1, JLO, JHI, J 
      DOUBLE PRECISION A, W 
      DOUBLE PRECISION AN, SUM 
      DIMENSION A(NRA,N), W(N) 
      AN = 0.0D0 
      DO 20 I = 1,N 
        SUM = 0.0D0 
        I1 = I + MU + 1 
        JLO = MAX(I-ML,1) 
        JHI = MIN(I+MU,N) 
        DO J = JLO,JHI 
          SUM = SUM + ABS(A(I1-J,J))/W(J)
        END DO
        AN = MAX(AN,SUM*W(I)) 
   20   CONTINUE 
      DBNORM = AN 
      RETURN 
!----------------------- End of Function DBNORM ------------------------
      END                                           
!DECK DSRCMA                                                            
      SUBROUTINE DSRCMA (RSAV, ISAV, JOB, common_data) 
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_common 
      type(odepack_common_data), target, intent(inout) :: common_data 
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of     
! the Common blocks DLS001, DLSA01, which are used                      
! internally by one or more ODEPACK solvers.                            
!                                                                       
! RSAV = real array of length 240 or more.                              
! ISAV = integer array of length 46 or more.                            
! JOB  = flag indicating to save or restore the Common blocks:          
!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)       
!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)     
!        A call with JOB = 2 presumes a prior call with JOB = 1.        
!-----------------------------------------------------------------------
      INTEGER ISAV, JOB 
      INTEGER, pointer :: ILS(:), ILSA(:) 
      INTEGER I, LENRLS, LENILS, LENRLA, LENILA 
      DOUBLE PRECISION RSAV 
      DOUBLE PRECISION, pointer :: RLS(:), RLSA(:) 
      DIMENSION RSAV(*), ISAV(*) 
      SAVE LENRLS, LENILS, LENRLA, LENILA 
!      COMMON /DLS001/ RLS(218), ILS(37)                                
!      COMMON /DLSA01/ RLSA(22), ILSA(9)                                
      DATA LENRLS/218/, LENILS/37/, LENRLA/22/, LENILA/9/ 
!     Common block pointers                                             
      type(DLS001_type), pointer :: DLS001 
      type(DLSA01_type), pointer :: DLSA01 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      DLS001 => common_data%DLS001 
                                                                        
      tmp_ptr = c_loc(DLS001%reals(1)) 
      call c_f_pointer(tmp_ptr,RLS,[218]) 
                                                                        
      tmp_ptr = c_loc(DLS001%ints(1)) 
      call c_f_pointer(tmp_ptr,ILS,[37]) 
                                                                        
      DLSA01 => common_data%DLSA01 
                                                                        
      tmp_ptr = c_loc(DLSA01%reals(1)) 
      call c_f_pointer(tmp_ptr,RLSA,[22]) 
                                                                        
      tmp_ptr = c_loc(DLSA01%ints(1)) 
      call c_f_pointer(tmp_ptr,ILSA,[9]) 
                                                                        
!                                                                       
      IF (JOB .EQ. 2) GO TO 100 
      DO I = 1,LENRLS 
        RSAV(I) = RLS(I)
      END DO
      DO I = 1,LENRLA 
        RSAV(LENRLS+I) = RLSA(I)
      END DO
!                                                                       
      DO I = 1,LENILS 
        ISAV(I) = ILS(I)
      END DO
      DO I = 1,LENILA 
        ISAV(LENILS+I) = ILSA(I)
      END DO
!                                                                       
      RETURN 
!                                                                       
  100 CONTINUE 
      DO I = 1,LENRLS 
        RLS(I) = RSAV(I) 
      END DO
      DO I = 1,LENRLA 
        RLSA(I) = RSAV(LENRLS+I)
      END DO
!                                                                       
      DO I = 1,LENILS 
        ILS(I) = ISAV(I)
      END DO
      DO I = 1,LENILA 
        ILSA(I) = ISAV(LENILS+I)
      END DO
!                                                                       
      RETURN 
!----------------------- End of Subroutine DSRCMA ----------------------
      END                                           
!DECK DRCHEK                                                            
      SUBROUTINE DRCHEK (JOB, G, NEQ, Y, YH,NYH, G0, G1, GX, JROOT, IRT,&
     &                   common_data)                                   
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_interface, only: DINTDY, DROOTS 
      use odepack_common 
      type(odepack_common_data), target, intent(inout) :: common_data 
      EXTERNAL G 
      INTEGER JOB, NEQ, NYH, JROOT, IRT 
      DOUBLE PRECISION Y, YH, G0, G1, GX 
      DIMENSION NEQ(*), Y(*), YH(NYH,*), G0(*), G1(*), GX(*), JROOT(*) 
      INTEGER, pointer :: IOWND(:), IOWNS(:),                           &
     &   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                     &
     &   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,               &
     &   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU         
      INTEGER, pointer :: IOWND3(:), IOWNR3(:), IRFND, ITASKC, NGC,     &
     &   NGE                                                            
      DOUBLE PRECISION, pointer :: ROWNS(:),                            &
     &   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND                  
      DOUBLE PRECISION, pointer :: ROWNR3(:), T0, TLAST, TOUTC 
!      COMMON /DLS001/ ROWNS(209),                                      
!     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,                
!     2   IOWND(6), IOWNS(6),                                           
!     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,                    
!     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,              
!     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU        
!      COMMON /DLSR01/ ROWNR3(2), T0, TLAST, TOUTC,                     
!     1   IOWND3(3), IOWNR3(2), IRFND, ITASKC, NGC, NGE                 
      INTEGER I, IFLAG, JFLAG 
      DOUBLE PRECISION HMING, T1, TEMP1, TEMP2, X 
      LOGICAL ZROOT 
!-----------------------------------------------------------------------
! This routine checks for the presence of a root in the vicinity of     
! the current T, in a manner depending on the input flag JOB.  It calls 
! Subroutine DROOTS to locate the root as precisely as possible.        
!                                                                       
! In addition to variables described previously, DRCHEK                 
! uses the following for communication:                                 
! JOB    = integer flag indicating type of call:                        
!          JOB = 1 means the problem is being initialized, and DRCHEK   
!                  is to look for a root at or very near the initial T. 
!          JOB = 2 means a continuation call to the solver was just     
!                  made, and DRCHEK is to check for a root in the       
!                  relevant part of the step last taken.                
!          JOB = 3 means a successful step was just taken, and DRCHEK   
!                  is to look for a root in the interval of the step.   
! G0     = array of length NG, containing the value of g at T = T0.     
!          G0 is input for JOB .ge. 2, and output in all cases.         
! G1,GX  = arrays of length NG for work space.                          
! IRT    = completion flag:                                             
!          IRT = 0  means no root was found.                            
!          IRT = -1 means JOB = 1 and a root was found too near to T.   
!          IRT = 1  means a legitimate root was found (JOB = 2 or 3).   
!                   On return, T0 is the root location, and Y is the    
!                   corresponding solution vector.                      
! T0     = value of T at one endpoint of interval of interest.  Only    
!          roots beyond T0 in the direction of integration are sought.  
!          T0 is input if JOB .ge. 2, and output in all cases.          
!          T0 is updated by DRCHEK, whether a root is found or not.     
! TLAST  = last value of T returned by the solver (input only).         
! TOUTC  = copy of TOUT (input only).                                   
! IRFND  = input flag showing whether the last step taken had a root.   
!          IRFND = 1 if it did, = 0 if not.                             
! ITASKC = copy of ITASK (input only).                                  
! NGC    = copy of NG (input only).                                     
!-----------------------------------------------------------------------
!                                                                       
!     Common block pointers                                             
      type(DLS001_type), pointer :: DLS001 
      type(DLSR01_type), pointer :: DLSR01 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      DLS001 => common_data%DLS001 
                                                                        
      tmp_ptr = c_loc(DLS001%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWNS,[209]) 
      tmp_ptr = c_loc(DLS001%reals(210)) 
      call c_f_pointer(tmp_ptr,CCMAX) 
      tmp_ptr = c_loc(DLS001%reals(211)) 
      call c_f_pointer(tmp_ptr,EL0) 
      tmp_ptr = c_loc(DLS001%reals(212)) 
      call c_f_pointer(tmp_ptr,H) 
      tmp_ptr = c_loc(DLS001%reals(213)) 
      call c_f_pointer(tmp_ptr,HMIN) 
      tmp_ptr = c_loc(DLS001%reals(214)) 
      call c_f_pointer(tmp_ptr,HMXI) 
      tmp_ptr = c_loc(DLS001%reals(215)) 
      call c_f_pointer(tmp_ptr,HU) 
      tmp_ptr = c_loc(DLS001%reals(216)) 
      call c_f_pointer(tmp_ptr,RC) 
      tmp_ptr = c_loc(DLS001%reals(217)) 
      call c_f_pointer(tmp_ptr,TN) 
      tmp_ptr = c_loc(DLS001%reals(218)) 
      call c_f_pointer(tmp_ptr,UROUND) 
                                                                        
      tmp_ptr = c_loc(DLS001%ints(1)) 
      call c_f_pointer(tmp_ptr,IOWND,[6]) 
      tmp_ptr = c_loc(DLS001%ints(7)) 
      call c_f_pointer(tmp_ptr,IOWNS,[6]) 
      tmp_ptr = c_loc(DLS001%ints(13)) 
      call c_f_pointer(tmp_ptr,ICF) 
      tmp_ptr = c_loc(DLS001%ints(14)) 
      call c_f_pointer(tmp_ptr,IERPJ) 
      tmp_ptr = c_loc(DLS001%ints(15)) 
      call c_f_pointer(tmp_ptr,IERSL) 
      tmp_ptr = c_loc(DLS001%ints(16)) 
      call c_f_pointer(tmp_ptr,JCUR) 
      tmp_ptr = c_loc(DLS001%ints(17)) 
      call c_f_pointer(tmp_ptr,JSTART) 
      tmp_ptr = c_loc(DLS001%ints(18)) 
      call c_f_pointer(tmp_ptr,KFLAG) 
      tmp_ptr = c_loc(DLS001%ints(19)) 
      call c_f_pointer(tmp_ptr,L) 
      tmp_ptr = c_loc(DLS001%ints(20)) 
      call c_f_pointer(tmp_ptr,LYH) 
      tmp_ptr = c_loc(DLS001%ints(21)) 
      call c_f_pointer(tmp_ptr,LEWT) 
      tmp_ptr = c_loc(DLS001%ints(22)) 
      call c_f_pointer(tmp_ptr,LACOR) 
      tmp_ptr = c_loc(DLS001%ints(23)) 
      call c_f_pointer(tmp_ptr,LSAVF) 
      tmp_ptr = c_loc(DLS001%ints(24)) 
      call c_f_pointer(tmp_ptr,LWM) 
      tmp_ptr = c_loc(DLS001%ints(25)) 
      call c_f_pointer(tmp_ptr,LIWM) 
      tmp_ptr = c_loc(DLS001%ints(26)) 
      call c_f_pointer(tmp_ptr,METH) 
      tmp_ptr = c_loc(DLS001%ints(27)) 
      call c_f_pointer(tmp_ptr,MITER) 
      tmp_ptr = c_loc(DLS001%ints(28)) 
      call c_f_pointer(tmp_ptr,MAXORD) 
      tmp_ptr = c_loc(DLS001%ints(29)) 
      call c_f_pointer(tmp_ptr,MAXCOR) 
      tmp_ptr = c_loc(DLS001%ints(30)) 
      call c_f_pointer(tmp_ptr,MSBP) 
      tmp_ptr = c_loc(DLS001%ints(31)) 
      call c_f_pointer(tmp_ptr,MXNCF) 
      tmp_ptr = c_loc(DLS001%ints(32)) 
      call c_f_pointer(tmp_ptr,N) 
      tmp_ptr = c_loc(DLS001%ints(33)) 
      call c_f_pointer(tmp_ptr,NQ) 
      tmp_ptr = c_loc(DLS001%ints(34)) 
      call c_f_pointer(tmp_ptr,NST) 
      tmp_ptr = c_loc(DLS001%ints(35)) 
      call c_f_pointer(tmp_ptr,NFE) 
      tmp_ptr = c_loc(DLS001%ints(36)) 
      call c_f_pointer(tmp_ptr,NJE) 
      tmp_ptr = c_loc(DLS001%ints(37)) 
      call c_f_pointer(tmp_ptr,NQU) 
                                                                        
      DLSR01 => common_data%DLSR01 
                                                                        
      tmp_ptr = c_loc(DLSR01%reals(1)) 
      call c_f_pointer(tmp_ptr,ROWNR3,[2]) 
      tmp_ptr = c_loc(DLSR01%reals(3)) 
      call c_f_pointer(tmp_ptr,T0) 
      tmp_ptr = c_loc(DLSR01%reals(4)) 
      call c_f_pointer(tmp_ptr,TLAST) 
      tmp_ptr = c_loc(DLSR01%reals(5)) 
      call c_f_pointer(tmp_ptr,TOUTC) 
                                                                        
      tmp_ptr = c_loc(DLSR01%ints(1)) 
      call c_f_pointer(tmp_ptr,IOWND3,[3]) 
      tmp_ptr = c_loc(DLSR01%ints(4)) 
      call c_f_pointer(tmp_ptr,IOWNR3,[2]) 
      tmp_ptr = c_loc(DLSR01%ints(6)) 
      call c_f_pointer(tmp_ptr,IRFND) 
      tmp_ptr = c_loc(DLSR01%ints(7)) 
      call c_f_pointer(tmp_ptr,ITASKC) 
      tmp_ptr = c_loc(DLSR01%ints(8)) 
      call c_f_pointer(tmp_ptr,NGC) 
      tmp_ptr = c_loc(DLSR01%ints(9)) 
      call c_f_pointer(tmp_ptr,NGE) 
                                                                        
!                                                                       
!                                                                       
      IRT = 0 
      DO I = 1,NGC 
        JROOT(I) = 0
      END DO
      HMING = (ABS(TN) + ABS(H))*UROUND*100.0D0 
!                                                                       
      GO TO (100, 200, 300), JOB 
!                                                                       
! Evaluate g at initial T, and check for zero values. ------------------
  100 CONTINUE 
      T0 = TN 
      CALL G (NEQ, T0, Y, NGC, G0, common_data%ierr) 
      if (common_data%ierr < 0) return 
      NGE = 1 
      ZROOT = .FALSE. 
      DO I = 1,NGC 
        IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE. 
      END DO
      IF (.NOT. ZROOT) GO TO 190 
! g has a zero at T.  Look at g at T + (small increment). --------------
      TEMP2 = MAX(HMING/ABS(H), 0.1D0) 
      TEMP1 = TEMP2*H 
      T0 = T0 + TEMP1 
      DO I = 1,N 
        Y(I) = Y(I) + TEMP2*YH(I,2)
      END DO
      CALL G (NEQ, T0, Y, NGC, G0, common_data%ierr) 
      if (common_data%ierr < 0) return 
      NGE = NGE + 1 
      ZROOT = .FALSE. 
      DO I = 1,NGC 
        IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE. 
      END DO
      IF (.NOT. ZROOT) GO TO 190 
! g has a zero at T and also close to T.  Take error return. -----------
      IRT = -1 
      RETURN 
!                                                                       
  190 CONTINUE 
      RETURN 
!                                                                       
!                                                                       
  200 CONTINUE 
      IF (IRFND .EQ. 0) GO TO 260 
! If a root was found on the previous step, evaluate G0 = g(T0). -------
      CALL DINTDY (T0, 0, YH, NYH, Y, IFLAG, common_data) 
      CALL G (NEQ, T0, Y, NGC, G0, common_data%ierr) 
      if (common_data%ierr < 0) return 
      NGE = NGE + 1 
      ZROOT = .FALSE. 
      DO I = 1,NGC 
        IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
      END DO
      IF (.NOT. ZROOT) GO TO 260 
! g has a zero at T0.  Look at g at T + (small increment). -------------
      TEMP1 = SIGN(HMING,H) 
      T0 = T0 + TEMP1 
      IF ((T0 - TN)*H .LT. 0.0D0) GO TO 230 
      TEMP2 = TEMP1/H 
      DO I = 1,N 
        Y(I) = Y(I) + TEMP2*YH(I,2)
      END DO
      GO TO 240 
  230 CALL DINTDY (T0, 0, YH, NYH, Y, IFLAG, common_data) 
  240 CALL G (NEQ, T0, Y, NGC, G0, common_data%ierr) 
      if (common_data%ierr < 0) return 
      NGE = NGE + 1 
      ZROOT = .FALSE. 
      DO 250 I = 1,NGC 
        IF (ABS(G0(I)) .GT. 0.0D0) GO TO 250 
        JROOT(I) = 1 
        ZROOT = .TRUE. 
  250   CONTINUE 
      IF (.NOT. ZROOT) GO TO 260 
! g has a zero at T0 and also close to T0.  Return root. ---------------
      IRT = 1 
      RETURN 
! G0 has no zero components.  Proceed to check relevant interval. ------
  260 IF (TN .EQ. TLAST) GO TO 390 
!                                                                       
  300 CONTINUE 
! Set T1 to TN or TOUTC, whichever comes first, and get g at T1. -------
      IF (ITASKC.EQ.2 .OR. ITASKC.EQ.3 .OR. ITASKC.EQ.5) GO TO 310 
      IF ((TOUTC - TN)*H .GE. 0.0D0) GO TO 310 
      T1 = TOUTC 
      IF ((T1 - T0)*H .LE. 0.0D0) GO TO 390 
      CALL DINTDY (T1, 0, YH, NYH, Y, IFLAG, common_data) 
      GO TO 330 
  310 T1 = TN 
      DO I = 1,N 
        Y(I) = YH(I,1)
      END DO
  330 CALL G (NEQ, T1, Y, NGC, G1, common_data%ierr) 
      if (common_data%ierr < 0) return 
      NGE = NGE + 1 
! Call DROOTS to search for root in interval from T0 to T1. ------------
      JFLAG = 0 
  350 CONTINUE 
      CALL DROOTS (NGC, HMING, JFLAG, T0, T1, G0, G1, GX, X, JROOT,     &
     &             common_data)                                         
      IF (JFLAG .GT. 1) GO TO 360 
      CALL DINTDY (X, 0, YH, NYH, Y, IFLAG, common_data) 
      CALL G (NEQ, X, Y, NGC, GX, common_data%ierr) 
      if (common_data%ierr < 0) return 
      NGE = NGE + 1 
      GO TO 350 
  360 T0 = X 
      CALL DCOPY (NGC, GX, 1, G0, 1) 
      IF (JFLAG .EQ. 4) GO TO 390 
! Found a root.  Interpolate to X and return. --------------------------
      CALL DINTDY (X, 0, YH, NYH, Y, IFLAG, common_data) 
      IRT = 1 
      RETURN 
!                                                                       
  390 CONTINUE 
      RETURN 
!----------------------- End of Subroutine DRCHEK ----------------------
      END                                           
!DECK DROOTS                                                            
      SUBROUTINE DROOTS (NG, HMIN, JFLAG, X0, X1, G0, G1, GX, X, JROOT, &
     &                   common_data)                                   
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_common 
      type(odepack_common_data), target, intent(inout) :: common_data 
      INTEGER NG, JFLAG, JROOT 
      DOUBLE PRECISION HMIN, X0, X1, G0, G1, GX, X 
      DIMENSION G0(NG), G1(NG), GX(NG), JROOT(NG) 
      INTEGER, pointer :: IOWND3(:), IMAX, LAST, IDUM3(:) 
      DOUBLE PRECISION, pointer :: ALPHA, X2, RDUM3(:) 
!      COMMON /DLSR01/ ALPHA, X2, RDUM3(3),                             
!     1   IOWND3(3), IMAX, LAST, IDUM3(4)                               
!-----------------------------------------------------------------------
! This subroutine finds the leftmost root of a set of arbitrary         
! functions gi(x) (i = 1,...,NG) in an interval (X0,X1).  Only roots    
! of odd multiplicity (i.e. changes of sign of the gi) are found.       
! Here the sign of X1 - X0 is arbitrary, but is constant for a given    
! problem, and -leftmost- means nearest to X0.                          
! The values of the vector-valued function g(x) = (gi, i=1...NG)        
! are communicated through the call sequence of DROOTS.                 
! The method used is the Illinois algorithm.                            
!                                                                       
! Reference:                                                            
! Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined        
! Output Points for Solutions of ODEs, Sandia Report SAND80-0180,       
! February 1980.                                                        
!                                                                       
! Description of parameters.                                            
!                                                                       
! NG     = number of functions gi, or the number of components of       
!          the vector valued function g(x).  Input only.                
!                                                                       
! HMIN   = resolution parameter in X.  Input only.  When a root is      
!          found, it is located only to within an error of HMIN in X.   
!          Typically, HMIN should be set to something on the order of   
!               100 * UROUND * MAX(ABS(X0),ABS(X1)),                    
!          where UROUND is the unit roundoff of the machine.            
!                                                                       
! JFLAG  = integer flag for input and output communication.             
!                                                                       
!          On input, set JFLAG = 0 on the first call for the problem,   
!          and leave it unchanged until the problem is completed.       
!          (The problem is completed when JFLAG .ge. 2 on return.)      
!                                                                       
!          On output, JFLAG has the following values and meanings:      
!          JFLAG = 1 means DROOTS needs a value of g(x).  Set GX = g(X) 
!                    and call DROOTS again.                             
!          JFLAG = 2 means a root has been found.  The root is          
!                    at X, and GX contains g(X).  (Actually, X is the   
!                    rightmost approximation to the root on an interval 
!                    (X0,X1) of size HMIN or less.)                     
!          JFLAG = 3 means X = X1 is a root, with one or more of the gi 
!                    being zero at X1 and no sign changes in (X0,X1).   
!                    GX contains g(X) on output.                        
!          JFLAG = 4 means no roots (of odd multiplicity) were          
!                    found in (X0,X1) (no sign changes).                
!                                                                       
! X0,X1  = endpoints of the interval where roots are sought.            
!          X1 and X0 are input when JFLAG = 0 (first call), and         
!          must be left unchanged between calls until the problem is    
!          completed.  X0 and X1 must be distinct, but X1 - X0 may be   
!          of either sign.  However, the notion of -left- and -right-   
!          will be used to mean nearer to X0 or X1, respectively.       
!          When JFLAG .ge. 2 on return, X0 and X1 are output, and       
!          are the endpoints of the relevant interval.                  
!                                                                       
! G0,G1  = arrays of length NG containing the vectors g(X0) and g(X1),  
!          respectively.  When JFLAG = 0, G0 and G1 are input and       
!          none of the G0(i) should be zero.                            
!          When JFLAG .ge. 2 on return, G0 and G1 are output.           
!                                                                       
! GX     = array of length NG containing g(X).  GX is input             
!          when JFLAG = 1, and output when JFLAG .ge. 2.                
!                                                                       
! X      = independent variable value.  Output only.                    
!          When JFLAG = 1 on output, X is the point at which g(x)       
!          is to be evaluated and loaded into GX.                       
!          When JFLAG = 2 or 3, X is the root.                          
!          When JFLAG = 4, X is the right endpoint of the interval, X1. 
!                                                                       
! JROOT  = integer array of length NG.  Output only.                    
!          When JFLAG = 2 or 3, JROOT indicates which components        
!          of g(x) have a root at X.  JROOT(i) is 1 if the i-th         
!          component has a root, and JROOT(i) = 0 otherwise.            
!-----------------------------------------------------------------------
      INTEGER I, IMXOLD, NXLAST 
      DOUBLE PRECISION T2, TMAX, FRACINT, FRACSUB, ZERO,HALF,TENTH,FIVE 
      LOGICAL ZROOT, SGNCHG, XROOT 
      SAVE ZERO, HALF, TENTH, FIVE 
      DATA ZERO/0.0D0/, HALF/0.5D0/, TENTH/0.1D0/, FIVE/5.0D0/ 
!     Common block pointers                                             
      type(DLSR01_type), pointer :: DLSR01 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      DLSR01 => common_data%DLSR01 
                                                                        
      tmp_ptr = c_loc(DLSR01%reals(1)) 
      call c_f_pointer(tmp_ptr,ALPHA) 
      tmp_ptr = c_loc(DLSR01%reals(2)) 
      call c_f_pointer(tmp_ptr,X2) 
      tmp_ptr = c_loc(DLSR01%reals(3)) 
      call c_f_pointer(tmp_ptr,RDUM3,[3]) 
                                                                        
      tmp_ptr = c_loc(DLSR01%ints(1)) 
      call c_f_pointer(tmp_ptr,IOWND3,[3]) 
      tmp_ptr = c_loc(DLSR01%ints(4)) 
      call c_f_pointer(tmp_ptr,IMAX) 
      tmp_ptr = c_loc(DLSR01%ints(5)) 
      call c_f_pointer(tmp_ptr,LAST) 
      tmp_ptr = c_loc(DLSR01%ints(6)) 
      call c_f_pointer(tmp_ptr,IDUM3,[4]) 
                                                                        
!                                                                       
      IF (JFLAG .EQ. 1) GO TO 200 
! JFLAG .ne. 1.  Check for change in sign of g or zero at X1. ----------
      IMAX = 0 
      TMAX = ZERO 
      ZROOT = .FALSE. 
      DO 120 I = 1,NG 
        IF (ABS(G1(I)) .GT. ZERO) GO TO 110 
        ZROOT = .TRUE. 
        GO TO 120 
! At this point, G0(i) has been checked and cannot be zero. ------------
  110   IF (SIGN(1.0D0,G0(I)) .EQ. SIGN(1.0D0,G1(I))) GO TO 120 
          T2 = ABS(G1(I)/(G1(I)-G0(I))) 
          IF (T2 .LE. TMAX) GO TO 120 
            TMAX = T2 
            IMAX = I 
  120   CONTINUE 
      IF (IMAX .GT. 0) GO TO 130 
      SGNCHG = .FALSE. 
      GO TO 140 
  130 SGNCHG = .TRUE. 
  140 IF (.NOT. SGNCHG) GO TO 400 
! There is a sign change.  Find the first root in the interval. --------
      XROOT = .FALSE. 
      NXLAST = 0 
      LAST = 1 
!                                                                       
! Repeat until the first root in the interval is found.  Loop point. ---
  150 CONTINUE 
      IF (XROOT) GO TO 300 
      IF (NXLAST .EQ. LAST) GO TO 160 
      ALPHA = 1.0D0 
      GO TO 180 
  160 IF (LAST .EQ. 0) GO TO 170 
      ALPHA = 0.5D0*ALPHA 
      GO TO 180 
  170 ALPHA = 2.0D0*ALPHA 
  180 X2 = X1 - (X1 - X0)*G1(IMAX) / (G1(IMAX) - ALPHA*G0(IMAX)) 
! If X2 is too close to X0 or X1, adjust it inward, by a fractional ----
! distance that is between 0.1 and 0.5. --------------------------------
      IF (ABS(X2 - X0) < HALF*HMIN) THEN 
        FRACINT = ABS(X1 - X0)/HMIN 
        FRACSUB = TENTH 
        IF (FRACINT .LE. FIVE) FRACSUB = HALF/FRACINT 
        X2 = X0 + FRACSUB*(X1 - X0) 
      ENDIF 
      IF (ABS(X1 - X2) < HALF*HMIN) THEN 
        FRACINT = ABS(X1 - X0)/HMIN 
        FRACSUB = TENTH 
        IF (FRACINT .LE. FIVE) FRACSUB = HALF/FRACINT 
        X2 = X1 - FRACSUB*(X1 - X0) 
      ENDIF 
      JFLAG = 1 
      X = X2 
! Return to the calling routine to get a value of GX = g(X). -----------
      RETURN 
! Check to see in which interval g changes sign. -----------------------
  200 IMXOLD = IMAX 
      IMAX = 0 
      TMAX = ZERO 
      ZROOT = .FALSE. 
      DO 220 I = 1,NG 
        IF (ABS(GX(I)) .GT. ZERO) GO TO 210 
        ZROOT = .TRUE. 
        GO TO 220 
! Neither G0(i) nor GX(i) can be zero at this point. -------------------
  210   IF (SIGN(1.0D0,G0(I)) .EQ. SIGN(1.0D0,GX(I))) GO TO 220 
          T2 = ABS(GX(I)/(GX(I) - G0(I))) 
          IF (T2 .LE. TMAX) GO TO 220 
            TMAX = T2 
            IMAX = I 
  220   CONTINUE 
      IF (IMAX .GT. 0) GO TO 230 
      SGNCHG = .FALSE. 
      IMAX = IMXOLD 
      GO TO 240 
  230 SGNCHG = .TRUE. 
  240 NXLAST = LAST 
      IF (.NOT. SGNCHG) GO TO 250 
! Sign change between X0 and X2, so replace X1 with X2. ----------------
      X1 = X2 
      CALL DCOPY (NG, GX, 1, G1, 1) 
      LAST = 1 
      XROOT = .FALSE. 
      GO TO 270 
  250 IF (.NOT. ZROOT) GO TO 260 
! Zero value at X2 and no sign change in (X0,X2), so X2 is a root. -----
      X1 = X2 
      CALL DCOPY (NG, GX, 1, G1, 1) 
      XROOT = .TRUE. 
      GO TO 270 
! No sign change between X0 and X2.  Replace X0 with X2. ---------------
  260 CONTINUE 
      CALL DCOPY (NG, GX, 1, G0, 1) 
      X0 = X2 
      LAST = 0 
      XROOT = .FALSE. 
  270 IF (ABS(X1-X0) .LE. HMIN) XROOT = .TRUE. 
      GO TO 150 
!                                                                       
! Return with X1 as the root.  Set JROOT.  Set X = X1 and GX = G1. -----
  300 JFLAG = 2 
      X = X1 
      CALL DCOPY (NG, G1, 1, GX, 1) 
      DO 320 I = 1,NG 
        JROOT(I) = 0 
        IF (ABS(G1(I)) .GT. ZERO) GO TO 310 
          JROOT(I) = 1 
          GO TO 320 
  310   IF (SIGN(1.0D0,G0(I)) .NE. SIGN(1.0D0,G1(I))) JROOT(I) = 1 
  320   CONTINUE 
      RETURN 
!                                                                       
! No sign change in the interval.  Check for zero at right endpoint. ---
  400 IF (.NOT. ZROOT) GO TO 420 
!                                                                       
! Zero value at X1 and no sign change in (X0,X1).  Return JFLAG = 3. ---
      X = X1 
      CALL DCOPY (NG, G1, 1, GX, 1) 
      DO 410 I = 1,NG 
        JROOT(I) = 0 
        IF (ABS(G1(I)) .LE. ZERO) JROOT (I) = 1 
  410 END DO 
      JFLAG = 3 
      RETURN 
!                                                                       
! No sign changes in this interval.  Set X = X1, return JFLAG = 4. -----
  420 CALL DCOPY (NG, G1, 1, GX, 1) 
      X = X1 
      JFLAG = 4 
      RETURN 
!----------------------- End of Subroutine DROOTS ----------------------
      END                                           
!DECK DSRCAR                                                            
      SUBROUTINE DSRCAR (RSAV, ISAV, JOB, common_data) 
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc 
      use odepack_common 
      type(odepack_common_data), target, intent(inout) :: common_data 
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of     
! the Common blocks DLS001, DLSA01, DLSR01, which are used              
! internally by one or more ODEPACK solvers.                            
!                                                                       
! RSAV = real array of length 245 or more.                              
! ISAV = integer array of length 55 or more.                            
! JOB  = flag indicating to save or restore the Common blocks:          
!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)       
!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)     
!        A call with JOB = 2 presumes a prior call with JOB = 1.        
!-----------------------------------------------------------------------
      INTEGER ISAV, JOB 
      INTEGER, pointer :: ILS(:), ILSA(:), ILSR(:) 
      INTEGER I, IOFF, LENRLS, LENILS, LENRLA, LENILA, LENRLR, LENILR 
      DOUBLE PRECISION RSAV 
      DOUBLE PRECISION, pointer :: RLS(:), RLSA(:), RLSR(:) 
      DIMENSION RSAV(*), ISAV(*) 
      SAVE LENRLS, LENILS, LENRLA, LENILA, LENRLR, LENILR 
!      COMMON /DLS001/ RLS(218), ILS(37)                                
!      COMMON /DLSA01/ RLSA(22), ILSA(9)                                
!      COMMON /DLSR01/ RLSR(5), ILSR(9)                                 
      DATA LENRLS/218/, LENILS/37/, LENRLA/22/, LENILA/9/ 
      DATA LENRLR/5/, LENILR/9/ 
!     Common block pointers                                             
      type(DLS001_type), pointer :: DLS001 
      type(DLSA01_type), pointer :: DLSA01 
      type(DLSR01_type), pointer :: DLSR01 
      type(c_ptr) :: tmp_ptr 
!-----------------------------------------------------------------------
! This code associates variables with common data                       
!-----------------------------------------------------------------------
      DLS001 => common_data%DLS001 
                                                                        
      tmp_ptr = c_loc(DLS001%reals(1)) 
      call c_f_pointer(tmp_ptr,RLS,[218]) 
                                                                        
      tmp_ptr = c_loc(DLS001%ints(1)) 
      call c_f_pointer(tmp_ptr,ILS,[37]) 
                                                                        
      DLSA01 => common_data%DLSA01 
                                                                        
      tmp_ptr = c_loc(DLSA01%reals(1)) 
      call c_f_pointer(tmp_ptr,RLSA,[22]) 
                                                                        
      tmp_ptr = c_loc(DLSA01%ints(1)) 
      call c_f_pointer(tmp_ptr,ILSA,[9]) 
                                                                        
      DLSR01 => common_data%DLSR01 
                                                                        
      tmp_ptr = c_loc(DLSR01%reals(1)) 
      call c_f_pointer(tmp_ptr,RLSR,[5]) 
                                                                        
      tmp_ptr = c_loc(DLSR01%ints(1)) 
      call c_f_pointer(tmp_ptr,ILSR,[9]) 
                                                                        
!                                                                       
      IF (JOB .EQ. 2) GO TO 100 
      DO I = 1,LENRLS 
        RSAV(I) = RLS(I)
      END DO
       DO I = 1,LENRLA 
         RSAV(LENRLS+I) = RLSA(I)
       END DO
      IOFF = LENRLS + LENRLA 
      DO I = 1,LENRLR 
        RSAV(IOFF+I) = RLSR(I)
      END DO
!                                                                       
      DO I = 1,LENILS 
        ISAV(I) = ILS(I)
      END DO
      DO I = 1,LENILA 
        ISAV(LENILS+I) = ILSA(I)
      END DO
      IOFF = LENILS + LENILA 
      DO I = 1,LENILR 
        ISAV(IOFF+I) = ILSR(I)
      END DO
!                                                                       
      RETURN 
!                                                                       
  100 CONTINUE 
      DO I = 1,LENRLS 
        RLS(I) = RSAV(I)
      END DO
      DO I = 1,LENRLA 
        RLSA(I) = RSAV(LENRLS+I)
      END DO
      IOFF = LENRLS + LENRLA 
      DO I = 1,LENRLR 
        RLSR(I) = RSAV(IOFF+I)
      END DO
!                                                                       
      DO I = 1,LENILS 
        ILS(I) = ISAV(I)
      END DO
      DO I = 1,LENILA 
        ILSA(I) = ISAV(LENILS+I) 
      END DO
      IOFF = LENILS + LENILA 
      DO I = 1,LENILR 
        ILSR(I) = ISAV(IOFF+I)
      END DO
!                                                                       
      RETURN 
!----------------------- End of Subroutine DSRCAR ----------------------
      END                                           
                                                                        
!=======================================================================
!                 START OF ODEPACK SUB1                                 
!=======================================================================
!DECK XERRWD                                                            
      SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2) 
!***BEGIN PROLOGUE  XERRWD                                              
!***SUBSIDIARY                                                          
!***PURPOSE  Write error message with values.                           
!***CATEGORY  R3C                                                       
!***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)                     
!***AUTHOR  Hindmarsh, Alan C., (LLNL)                                  
!***DESCRIPTION                                                         
!                                                                       
!  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,   
!  as given here, constitute a simplified version of the SLATEC error   
!  handling package.                                                    
!                                                                       
!  All arguments are input arguments.                                   
!                                                                       
!  MSG    = The message (character array).                              
!  NMES   = The length of MSG (number of characters).                   
!  NERR   = The error number (not used).                                
!  LEVEL  = The error level..                                           
!           0 or 1 means recoverable (control returns to caller).       
!           2 means fatal (run is aborted--see note below).             
!  NI     = Number of integers (0, 1, or 2) to be printed with message. 
!  I1,I2  = Integers to be printed, depending on NI.                    
!  NR     = Number of reals (0, 1, or 2) to be printed with message.    
!  R1,R2  = Reals to be printed, depending on NR.                       
!                                                                       
!  Note..  this routine is machine-dependent and specialized for use    
!  in limited context, in the following ways..                          
!  1. The argument MSG is assumed to be of type CHARACTER, and          
!     the message is printed with a format of (1X,A).                   
!  2. The message is assumed to take only one line.                     
!     Multi-line messages are generated by repeated calls.              
!  3. If LEVEL = 2, control passes to the statement   STOP              
!     to abort the run.  This statement may be machine-dependent.       
!  4. R1 and R2 are assumed to be in double precision and are printed   
!     in D21.13 format.                                                 
!                                                                       
!***ROUTINES CALLED  IXSAV                                              
!***REVISION HISTORY  (YYMMDD)                                          
!   920831  DATE WRITTEN                                                
!   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)                      
!   930329  Modified prologue to SLATEC format. (FNF)                   
!   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)       
!   930922  Minor cosmetic change. (FNF)                                
!***END PROLOGUE  XERRWD                                                
!                                                                       
!*Internal Notes:                                                       
!                                                                       
! For a different default logical unit number, IXSAV (or a subsidiary   
! routine that it calls) will need to be modified.                      
! For a different run-abort command, change the statement following     
! statement 100 at the end.                                             
!-----------------------------------------------------------------------
! Subroutines called by XERRWD.. None                                   
! Function routine called by XERRWD.. IXSAV                             
!-----------------------------------------------------------------------
!**End                                                                  
!                                                                       
!  Declare arguments.                                                   
!                                                                       
      DOUBLE PRECISION R1, R2 
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR 
      CHARACTER*(*) MSG 
!                                                                       
!  Declare local variables.                                             
!                                                                       
      INTEGER LUNIT, IXSAV, MESFLG 
!                                                                       
!  Get logical unit number and message print flag. 

!  Suppress compiler warnings      
      IF (NMES == NMES .AND. NERR == NERR) CONTINUE                    
!                                                                       
!***FIRST EXECUTABLE STATEMENT  XERRWD                                  
      LUNIT = IXSAV (1, 0, .FALSE.) 
      MESFLG = IXSAV (2, 0, .FALSE.) 
      IF (MESFLG .EQ. 0) GO TO 100 
!                                                                       
!  Write the message.                                                   
!                                                                       
      WRITE (LUNIT,10)  MSG 
   10 FORMAT(1X,A) 
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1 
   20 FORMAT(6X,'In above message,  I1 =',I10) 
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2 
   30 FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10) 
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1 
   40 FORMAT(6X,'In above message,  R1 =',D21.13) 
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2 
   50 FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13) 
!                                                                       
!  Abort the run if LEVEL = 2.                                          
!                                                                       
  100 IF (LEVEL .NE. 2) RETURN 
      STOP 
!----------------------- End of Subroutine XERRWD ----------------------
      END                                           
!DECK IXSAV                                                             
      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET) 
!***BEGIN PROLOGUE  IXSAV                                               
!***SUBSIDIARY                                                          
!***PURPOSE  Save and recall error message control parameters.          
!***CATEGORY  R3C                                                       
!***TYPE      ALL (IXSAV-A)                                             
!***AUTHOR  Hindmarsh, Alan C., (LLNL)                                  
!***DESCRIPTION                                                         
!                                                                       
!  IXSAV saves and recalls one of two error message parameters:         
!    LUNIT, the logical unit number to which messages are printed, and  
!    MESFLG, the message print flag.                                    
!  This is a modification of the SLATEC library routine J4SAVE.         
!                                                                       
!  Saved local variables..                                              
!   LUNIT  = Logical unit number for messages.  The default is obtained 
!            by a call to IUMACH (may be machine-dependent).            
!   MESFLG = Print control flag..                                       
!            1 means print all messages (the default).                  
!            0 means no printing.                                       
!                                                                       
!  On input..                                                           
!    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).          
!    IVALUE = The value to be set for the parameter, if ISET = .TRUE.   
!    ISET   = Logical flag to indicate whether to read or write.        
!             If ISET = .TRUE., the parameter will be given             
!             the value IVALUE.  If ISET = .FALSE., the parameter       
!             will be unchanged, and IVALUE is a dummy argument.        
!                                                                       
!  On return..                                                          
!    IXSAV = The (old) value of the parameter.                          
!                                                                       
!***SEE ALSO  XERRWD, XERRWV                                            
!***ROUTINES CALLED  IUMACH                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   921118  DATE WRITTEN                                                
!   930329  Modified prologue to SLATEC format. (FNF)                   
!   930915  Added IUMACH call to get default output unit.  (ACH)        
!   930922  Minor cosmetic changes. (FNF)                               
!   010425  Type declaration for IUMACH added. (ACH)                    
!***END PROLOGUE  IXSAV                                                 
!                                                                       
! Subroutines called by IXSAV.. None                                    
! Function routine called by IXSAV.. IUMACH                             
!-----------------------------------------------------------------------
!**End                                                                  
      LOGICAL ISET 
      INTEGER IPAR, IVALUE 
!-----------------------------------------------------------------------
      INTEGER IUMACH, LUNIT, MESFLG 
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the    
! listed (local) variables to be saved between calls to this routine.   
!-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG 
      DATA LUNIT/-1/, MESFLG/1/ 
!                                                                       
!***FIRST EXECUTABLE STATEMENT  IXSAV                                   
      IF (IPAR .EQ. 1) THEN 
        IF (LUNIT .EQ. -1) LUNIT = IUMACH() 
        IXSAV = LUNIT 
        IF (ISET) LUNIT = IVALUE 
        ENDIF 
!                                                                       
      IF (IPAR .EQ. 2) THEN 
        IXSAV = MESFLG 
        IF (ISET) MESFLG = IVALUE 
        ENDIF 
!                                                                       
      RETURN 
!----------------------- End of Function IXSAV -------------------------
      END                                           
!DECK IUMACH                                                            
      INTEGER FUNCTION IUMACH() 
!***BEGIN PROLOGUE  IUMACH                                              
!***PURPOSE  Provide standard output unit number.                       
!***CATEGORY  R1                                                        
!***TYPE      INTEGER (IUMACH-I)                                        
!***KEYWORDS  MACHINE CONSTANTS                                         
!***AUTHOR  Hindmarsh, Alan C., (LLNL)                                  
!***DESCRIPTION                                                         
! *Usage:                                                               
!        INTEGER  LOUT, IUMACH                                          
!        LOUT = IUMACH()                                                
!                                                                       
! *Function Return Values:                                              
!     LOUT : the standard logical unit for Fortran output.              
!                                                                       
!***REFERENCES  (NONE)                                                  
!***ROUTINES CALLED  (NONE)                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   930915  DATE WRITTEN                                                
!   930922  Made user-callable, and other cosmetic changes. (FNF)       
!***END PROLOGUE  IUMACH                                                
!                                                                       
!*Internal Notes:                                                       
!  The built-in value of 6 is standard on a wide range of Fortran       
!  systems.  This may be machine-dependent.                             
!**End                                                                  
!***FIRST EXECUTABLE STATEMENT  IUMACH                                  
      IUMACH = 6 
!                                                                       
      RETURN 
!----------------------- End of Function IUMACH ------------------------
      END                                           
