program radiation_map
!use madx_ptc_module
use pointer_lattice
use gauss_dis

implicit none
type(layout), pointer:: ring
real(dp) prec,x0(6)  
real(dp) k1 
type(c_ray) f_ray
type(c_normal_form) normal_form
 integer  no,np,neq,NV
integer i,j,k, case_section
type(probe_8) xs  
!!!!!!!!!!!!!!!!!!!!!
type(fibre), pointer :: p,p2
type(probe) ray_TC 
type(probe_8) ray_8_TC
type(c_damap)  one_turn_map_AP, identity_AP,two_turn_map_AP,lagrange_map_ap,u_c
type(c_damap) a0,as,a1,a2,a_cs,rot,a_tracked
type(c_damap) U,Q,U_0,U_1,U_2
type(c_linear_map) q_cs,q_as,q_rot,q_orb  
real(dp)  pha(3),damp(3)
real(dp)  spin_tune(2)
type(c_vector_field) fc
real(dp)  closed_orbit(6),goal(9),entrance_orbit(6),nu_spin
type(internal_state),target :: bmad_state
type(c_taylor) ct 
type(c_taylor) phase(3),phase_spin
type(c_universal_taylor) tpos(9)  
integer ignore_closed_orbit
integer npmax
integer, allocatable :: jindex(:)
real nave

c_verbose=.false.
prec=1.d-10 ! for printing
use_info=.true.
use_quaternion=.true.
 c_lda_used=2000
 lda_used=2000
 n_cai=-i_
preffile= "C:\document\etienne_programs_new\programs_for_learning\suntao\3nuy\pref.txt"

 
 
call ptc_ini_no_append

 
 call append_empty_layout(m_u)
 
ring=>m_u%start

 call build_lattice_als(ring,mis=.true.,exact=.true., thin=.false. ,onecell=.true.)


call make_node_layout(ring)
 
!call in_bmad_units()

bmad_state=only_4d0+delta0+spin0
bmad_state=nocavity0+spin0

p=>ring%start
call move_to(ring,p,"QF1")   !  Locating QF1 in TC
call make_it_knob(p%magp%bn(2),1)  !  PP : making a Polymorph into a knob
call move_to(ring,p,"QF2")   ! Locating QF2 in TC
call make_it_knob(p%magp%an(2),2)  !  PP : making a Polymorph into a knob

p=>ring%start
closed_orbit=0
call find_orbit_x(closed_orbit,bmad_state, 1.0e-7_dp, fibre1=p)
write(6,*) " closed orbit "
write(6,"(6(1x,g11.4))") closed_orbit
 x0=closed_orbit


write(6,*) " Give order no >=2 for beta "
read(5,*) no


!!!!!!  calling TC !!!!!!
np=2     ! number of free parameters to solve these equations
call init(bmad_state,no,np)   ! Berz's TPSA package is initialized
nV=c_%nv   !     
call alloc(ct,phase_spin)  
call alloc(identity_AP,one_turn_map_AP,a0,as,a1,a2,a_cs,rot,a_tracked,u_c);
call alloc(normal_form);
call alloc(ray_8_TC);call alloc(phase) 
p=>ring%start
   
!  Usual computation of a Taylor maps within TC
ray_TC=x0   !  For TC
identity_AP=1
ray_8_TC=ray_TC + identity_AP ! Connect AP with TC: identity added
 
p=>ring%start
call propagate(ray_8_TC,+bmad_state,fibre1=p) !  TC of PTC proper
 
one_turn_map_AP=ray_8_TC



! compute 2 tunes and the phase slip
call c_normal(one_turn_map_AP,normal_form,phase=phase,dospin=bmad_state%spin) 
 
write(6,*) "Transverse tunes "
write(6,"(3(1x,g17.10,1x))") normal_form%tune(1:c_%nd)
 
call clean(normal_form%H_l,normal_form%H_l,prec=prec)
call clean(normal_form%H_nl,normal_form%H_nl,prec=prec)

 
!
write(6,*) 
write(6,*) "z_1-component of the Linear Vector Field in Phasors Variables  "
call print(normal_form%H_l%v(1));write(6,*);

write(6,*) "z_1-component of the Non-Linear Vector Field in Phasors Variables (2nd order)"
call print(normal_form%H_nl%v(1).cut.4)


if(no>1) then

!# a_t = a_0 o a_1 o a_2 o a_s  for dir=1
call c_full_factorise(normal_form%atot,as,a0,a1,a2,dir=1) 

! computing the "beta function" two different ways
! x^2+p^2 is invariant of normal form n
ct=(1.d0.cmono.'2')+(1.d0.cmono.'02')
! (x^2+p^2) o a^-1 is invariant around the fixed point
ct=ct*a1**(-1)
allocate(jindex(c_%nd2harm))
call clean(ct,ct,prec=prec)
jindex=0
jindex(2)=2
ct=ct.par.jindex
Write(6,*) "Beta via Coefficient of p^2 in invariant"
call print(ct)
!  function x^2 is constructed
ct=(1.d0.cmono.'2')
!  function x^2 is averaged using a_1
call C_AVERAGE(ct,a1,ct) 
jindex=0
jindex(1)=1
jindex(2)=1
!  coefficient of J_1 extracted
ct=ct.par.jindex
call clean(ct,ct,prec=prec)
write(6,*) " Beta via average "
call print(ct)

endif

!!!! Canonize  !!!!
  phase(1)=0.d0
  phase(2)=0.d0
  phase(3)=0.d0
phase_spin=0.d0

!  Only necessary for phase advance 
call c_full_canonise(normal_form%atot,a_cs)
 
p2=>ring%start
call move_to(ring,p2,"SF",reset=.true.)   !  Locating SF in TC
 
!  Usual computation of a Taylor maps within TC
ray_TC=x0   !  For TC
ray_8_TC=ray_TC + a_cs ! Connect AP with TC: identity added
 
p=>ring%start
call propagate(ray_8_TC,+bmad_state,fibre1=p,fibre2=p2) !  TC of PTC proper
a_tracked=ray_8_TC


call c_full_canonise(a_tracked,a_cs,as=as,a0=a0,a1=a1,a2=a2,rotation=rot  &
,phase=phase,nu_spin=phase_spin)
call clean(as,as,prec=prec)

 
call clean(a0,a0,prec=prec)
call clean(a1,a1,prec=prec)
call clean(a2,a2,prec=prec)
call clean(a_cs,a_cs,prec=prec)
call clean(phase,phase,prec=prec)
call clean(phase_spin,phase_spin,prec=prec)

write(6,*) " Spin transformation "
call print(as)
write(6,*) " Fixed point transformation "
call print(a0)
write(6,*) " Linear transformation "
call print(a1)
write(6,*) " Nonlinear transformation "
call print(a2)
write(6,*) " Canonized transformation "
call print(a_cs)
write(6,*) " Phase Advance map"
call print(rot)
write(6,*) " Transverse phase advance "
call print(phase(1:3))
pause 654
write(6,*) " Spin phase advance "
call print(phase_spin)
 
call kill(ct,phase_spin)  
call kill(identity_AP,one_turn_map_AP,a0,as,a1,a2,a_cs,rot,a_tracked,u_c);
call kill(normal_form);
call kill(ray_8_TC);call kill(phase); 
call kill_para(ring)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  calling TC !!!!!!
np=0
no=1
call init(bmad_state,no,np)   ! Berz's TPSA package is initialized
nV=c_%nv   !     
call alloc(ct,phase_spin)  
call alloc(identity_AP,one_turn_map_AP,a0,as,a1,a2,a_cs,rot,a_tracked,u_c);
call alloc(normal_form);
call alloc(ray_8_TC);call alloc(phase); 
p=>ring%start
   
!  Usual computation of a Taylor maps within TC
ray_TC=x0   !  For TC
identity_AP=1
ray_8_TC=ray_TC + identity_AP ! Connect AP with TC: identity added
 
p=>ring%start
call propagate(ray_8_TC,bmad_state,fibre1=p) !  TC of PTC proper
 
one_turn_map_AP=ray_8_TC

 
! compute 2 tunes and the phase slip
call c_normal(one_turn_map_AP,normal_form,phase=phase,dospin=bmad_state%spin) 
 
write(6,*) "Transverse tunes "
write(6,"(3(1x,g17.10,1x))") normal_form%tune(1:c_%nd)
 
 



!!!! Canonize  !!!!
pha=0
damp=0
spin_tune=0

!  Only necessary for phase advance 
!call c_full_canonise(normal_form%atot,a_cs)
call c_fast_canonise(normal_form%atot,a_cs,dospin=bmad_state%spin)  !,phase=pha,damping=damp,q_cs,q_as,q_orb,q_rot,spin_tune ,dospin)
 
 p2=>ring%start
 call move_to(ring,p2,"SF",reset=.true.)   !  Locating SF in TC
 

 !  Usual computation of a Taylor maps within TC
ray_TC=x0   !  For TC
ray_8_TC=ray_TC + a_cs ! Connect AP with TC: identity added
 
p=>ring%start
 
call propagate(ray_8_TC,bmad_state,fibre1=p,fibre2=p2) !  TC of PTC proper

 a_tracked=ray_8_TC
 
 
call c_fast_canonise(a_tracked,a_cs,phase=pha,damping=damp,q_cs=q_cs, &
q_as=q_as,q_orb=q_orb,q_rot=q_rot,spin_tune=spin_tune ,dospin=bmad_state%spin)
 
 
write(6,*) " Spin transformation "
call print(q_as)
write(6,*) " Linear transformation "
call print(q_orb)
 
write(6,*) " Canonized transformation "
call print(q_cs)
write(6,*) " Phase Advance map"
call print(q_rot)
write(6,*) " Transverse phase advance "
write(6,format3) pha
write(6,*) " Spin phase advance "
write(6,format2) spin_tune
 


q_rot=q_as*q_orb

call print(q_rot)

call print(a_cs%q)



call kill(ct,phase_spin)  
call kill(identity_AP,one_turn_map_AP,a0,as,a1,a2,a_cs,rot,a_tracked,u_c);
call kill(normal_form);
call kill(ray_8_TC);call kill(phase); 
 


!!!!!!!!!!!!!!!!!
stop 1000
normal_form%H=0.01d0*normal_form%H
a1=exp(normal_form%H)
a1=c_phasor()*a1*ci_phasor()
call clean(a1,a1,prec=prec)
call alloc(fc)
Fc=-c_logf_spin(a1)
 
a2=exp(fc,a1)
 
call clean(a2,a2,prec=prec)
call print(a2) 
 
!!!!   making map near identity





1000 call ptc_end(graphics_maybe=1,flat_file=.false.)
contains



subroutine  build_lattice_als(ALS,mis,error,exact,sl,thin,onecell)
use madx_ptc_module
use pointer_lattice
implicit none

type(layout), target :: ALS
real(dp),optional :: error(6)
logical, optional :: exact,sl,thin,onecell
real(dp) :: alpha,lbend, cut, ksd, ksf 
type(fibre)  L1,L2,L3,L4,L5,L6,L7,L8,L9,L10 
type(fibre)  L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,CAVM
type(fibre)  L21,L22,L23,L24,L25,L26,L27,L27A,L27B,L27C,L27D,DS
type(fibre)  QF1,QF2,QD1,QD2,QFA1,QFA2,sf,sd,cav,bend,vc5,bend1 
type(layout) :: sfline,sdline,sup1,supb
logical(lp) :: mis,thi=.false.,oneperiod
real(dp) sig(6)
!-----------------------------------
if(present(thin)) thi=thin


call make_states(.true.)
exact_model = .false.;oneperiod = .false.
if(present(exact)) exact_model=exact
if(present(onecell)) oneperiod=onecell
call update_states
madlength = .false.


call set_mad(energy = 1.5d0, method = 2, step = 8)

!madkind2 = matrix_kick_matrix
madkind2 = drift_kick_drift


L1  = drift("L1 ",  2.832695d0);L2  = drift("L2 ",  0.45698d0);
L3  = drift("L3 ",  0.08902d0);L4  = drift("L4 ",  0.2155d0);
L5  = drift("L5 ",  0.219d0);L6  = drift("L6 ",  0.107078d0);
L7  = drift("L7 ",  0.105716d0);L8  = drift("L8 ",  0.135904d0);
L9  = drift("L9 ",  0.2156993d0);L10 = drift("L10",  0.089084d0);
L11= drift("L11",  0.235416d0);L12= drift("L12",  0.1245d0);
L13= drift("L13",  0.511844d0);L14= drift("L14",  0.1788541d0);
L15= drift("L15",  0.1788483d0);L16= drift("L16",  0.511849d0);
L17= drift("L17",  0.1245d0);L18= drift("L18",  0.235405d0);
L19= drift("L19",  0.089095d0);L20= drift("L20",  0.2157007d0);
L21= drift("L21",  0.177716d0);L22= drift("L22",  0.170981d0);
L23= drift("L23",  0.218997d0);L24 = drift ("L24",  0.215503d0);
L25 = drift ("L25",  0.0890187d0);L26 = drift ("L26",  0.45698d0);
L27 = drift ("L27",  2.832696d0);L27a  = drift (" L27a",  0.8596d0);
L27b  = drift (" L27b",  0.1524d0);L27c  = drift (" L27c",  0.04445d0);
L27d  = drift (" L27d",  1.776246d0);ds  = drift (" DS  ", 0.1015d0);

QF1 = QUADRUPOLE(" QF1 ",0.344D0, K1= 2.2474D0+6.447435260914397D-03)
QF2 = QUADRUPOLE(" QF2 ",0.344D0, K1= 2.2474D0)
QD1 = QUADRUPOLE(" QD1 ",0.187D0, K1= -2.3368D0-2.593018157427161D-02); 
QD2 = QUADRUPOLE(" QD2 ",0.187D0, K1= -2.3368D0);  
QFA1= QUADRUPOLE(" QFA1",0.448D0, K1= 2.8856D0);  
QFA2= QUADRUPOLE(" QFA2",0.448D0, K1= 2.8856D0);  

!!! 1/2 mad-x value
ksf=-41.3355516397069748d0;
ksd=56.2564709584745489d0;

sf=sextupole ("sf",2.d0*0.1015d0, K2= ksf);
sd= sextupole("sd", 2.d0*0.1015d0, K2= ksd);

 VC5=marker("vc5");
ALPHA=0.17453292519943295769236907684886d0;
 
LBEND=0.86621d0;
 
BEND = RBEND("BEND", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
BEND1 = RBEND("BEND1", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
 
CAVM=MARK("CAVM");
N_CAV4_F=3
CAV=RFCAVITY("CAV",L=0.0000d0,VOLT=-1000.0d0,REV_FREQ=500.0d6)
!CAV=RFCAVITY("CAV",L=0.0000d0,VOLT=-1.0d0,REV_FREQ=500.0d6)
write(6,*) cav%mag%c4%nf
write(6,*) cav%mag%c4%f
cav%mag%c4%f(1)=1
cav%magp%c4%f(1)=cav%mag%c4%f(1)
!cav%mag%c4%f(3)=-0.3d0*cav%mag%c4%f(1)
!cav%magp%c4%f(3)=cav%mag%c4%f(3)

if(thi) then
 sf=sextupole ("sf",0.d0, K2= ksf*0.203d0);
 sd= sextupole("sd", 0.d0, K2= ksd*0.203d0);
  sfline=(ds+sf+ds);
  sdline=(ds+sd+ds);
else
 sfline=1*sf;
 sdline=1*sd;
endif

SUP1=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27+cavm;

SUPb=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND1+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27+cav;

if(oneperiod) then
 ALS = supb;  !11*sup1+supb;
else
 ALS = 11*sup1+supb;
endif
if(present(sl)) then
L1  = drift("L1 ",  2.832695d0);
 if( sl ) then
  Qf1 = QUADRUPOLE(" QF1 ",L=0.d0, K1= 0.01d0 ); L1  = drift("L1 ",L=0.1d0);
  ALS=L1+QF1;
 endif 
endif

ALS = .ring.ALS

call survey(ALS)

call gaussian_seed(1000)
if(mis) then
 sig=1.d-5; cut=4.d0; 
 if(present(error)) sig=error
 call MESS_UP_ALIGNMENT(ALS,SIG,cut);
endif
end subroutine build_lattice_als



end program radiation_map

