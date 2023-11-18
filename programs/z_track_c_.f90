program radiation_map
!use madx_ptc_module
use pointer_lattice
use gauss_dis

implicit none
type(layout), pointer:: ring
real(dp) prec,x0(6),closed_orbit(6) 
real(dp) k1 ,nu_mod,circ
type(c_ray) f_ray
type(c_normal_form) normal_form
type(probe),target :: xs0
type(probe_8) xs
 integer  no,np,nd,neq,NV,ndpt
integer i,j,k, case_section,ac_modulation
type(real_8) x1
type(complex_8) x2
type(c_damap) ids
!!!!!!!!!!!!!!!!!!!!!
type(fibre), pointer :: p,p_QF1,p_QF2
 
type(internal_state),target :: bmad_state
real(dp), allocatable  :: an(:),bn(:)
real(dp) DC_ac,A_ac,theta_ac, D_ac


c_verbose=.false.
prec=1.d-10 ! for printing
use_info=.true.
use_quaternion=.true.
 c_lda_used=1000
 lda_used=1000
 n_cai=-i_
 

call ptc_ini_no_append

 
 call append_empty_layout(m_u)
 
ring=>m_u%start

 call build_lattice_als(ring,mis=.true.,exact=.true., thin=.false. ,onecell=.true.)


call make_node_layout(ring)

p=>ring%start
call move_to(ring,p,"QF1")   !  Locating QF1 in TC
!  PP : making a Polymorph into a knob
call make_it_knob(p%magp%bn(2),1,s=1.0_dp)  
p_QF1=>p
call move_to(ring,p,"QF2")   ! Locating QF2 in TC
call make_it_knob(p%magp%bn(2),2,s=1.0_dp)  
p_QF2=>p


no=2
nv=7
knob=.true. 
call c_init_all(no,nv)

write(6,*) "Taylors created with two knobs : FPP no phase space maps "

write(6,'(a16,1x,i2)')  " c_%no = ",c_%no
write(6,'(a16,1x,i2)')  " c_%nv = ",c_%nv
write(6,'(a16,1x,i2)') " c_%nd = ",c_%nd
write(6,'(a16,1x,i2)')  " c_%nd2 = ",c_%nd2
write(6,'(a16,1x,i2)')  " c_%nd2harm = ",c_%nd2harm
write(6,'(a16,1x,i2)')  " c_%ndpt = ",c_%ndpt
write(6,'(a16,1x,i2)')  " c_%ndpt_bmad = ",c_%ndpt_bmad
write(6,'(a16,1x,i2)')  " c_%np = ",c_%np
write(6,'(a16,1x,i2)')  " c_%npara_fpp= ",c_%npara_fpp
write(6,'(a16,1x,i2)')  " c_%nd2t= ",c_%nd2t
write(6,'(a16,1x,i2)')  " c_%ndc2t= ",c_%ndc2t
write(6,'(a16,1x,i2)')  " c_%pos_of_delta= ",c_%pos_of_delta


write(6,'(a16,1x,i2)')  " c_%rf= ",c_%rf
write(6,'(a16,1x,L)')  " c_%setknob= ",c_%setknob
write(6,'(a16,1x,L)')  " c_%knob= ",c_%knob


 
no=2
nd=3
np=2
ndpt=0
knob=.true. 
call c_init_all(no,nd,np,ndpt)
  
 

write(6,*) "Taylors created with two knobs : FPP with phase space maps "
write(6,'(a16,1x,i2)')  " c_%no = ",c_%no
write(6,'(a16,1x,i2)')  " c_%nv = ",c_%nv
write(6,'(a16,1x,i2)')  " c_%nd = ",c_%nd
write(6,'(a16,1x,i2)')  " c_%nd2 = ",c_%nd2
write(6,'(a16,1x,i2)')  " c_%nd2harm = ",c_%nd2harm
write(6,'(a16,1x,i2)')  " c_%ndpt = ",c_%ndpt
write(6,'(a16,1x,i2)')  " c_%ndpt_bmad = ",c_%ndpt_bmad
write(6,'(a16,1x,i2)')  " c_%np = ",c_%np
write(6,'(a16,1x,i2)')  " c_%npara_fpp= ",c_%npara_fpp
write(6,'(a16,1x,i2)')  " c_%nd2t= ",c_%nd2t
write(6,'(a16,1x,i2)')  " c_%ndc2t= ",c_%ndc2t
write(6,'(a16,1x,i2)')  " c_%pos_of_delta= ",c_%pos_of_delta


write(6,'(a16,1x,i2)')  " c_%rf= ",c_%rf
write(6,'(a16,1x,L)')  " c_%setknob= ",c_%setknob
write(6,'(a16,1x,L)')  " c_%knob= ",c_%knob


call in_bmad_units()
bmad_state=default+time0      !+nocavity0

call init_all(bmad_state,no,np)
 
write(6,*) "Taylors created with two knobs : PTC with phase space maps "
write(6,*) "Cavities turned on "
write(6,'(a16,1x,i2)')  " c_%no = ",c_%no
write(6,'(a16,1x,i2)')  " c_%nv = ",c_%nv
write(6,'(a16,1x,i2)')  " c_%nd = ",c_%nd
write(6,'(a16,1x,i2)')  " c_%nd2 = ",c_%nd2
write(6,'(a16,1x,i2)')  " c_%nd2harm = ",c_%nd2harm
write(6,'(a16,1x,i2)')  " c_%ndpt = ",c_%ndpt
write(6,'(a16,1x,i2)')  " c_%ndpt_bmad = ",c_%ndpt_bmad
write(6,'(a16,1x,i2)')  " c_%np = ",c_%np
write(6,'(a16,1x,i2)')  " c_%npara_fpp= ",c_%npara_fpp
write(6,'(a16,1x,i2)')  " c_%nd2t= ",c_%nd2t
write(6,'(a16,1x,i2)')  " c_%ndc2t= ",c_%ndc2t
write(6,'(a16,1x,i2)')  " c_%pos_of_delta= ",c_%pos_of_delta


write(6,'(a16,1x,i2)')  " c_%rf= ",c_%rf
write(6,'(a16,1x,L)')  " c_%setknob= ",c_%setknob
write(6,'(a16,1x,L)')  " c_%knob= ",c_%knob


call in_bmad_units()
bmad_state=default+time0+nocavity0

call init_all(bmad_state,no,np)
 
write(6,*) "Taylors created with two knobs : PTC with phase space maps "
write(6,*) "Cavities turned off "
write(6,'(a16,1x,i2)')  " c_%no = ",c_%no
write(6,'(a16,1x,i2)')  " c_%nv = ",c_%nv
write(6,'(a16,1x,i2)')  " c_%nd = ",c_%nd
write(6,'(a16,1x,i2)')  " c_%nd2 = ",c_%nd2
write(6,'(a16,1x,i2)')  " c_%nd2harm = ",c_%nd2harm
write(6,'(a16,1x,i2)')  " c_%ndpt = ",c_%ndpt
write(6,'(a16,1x,i2)')  " c_%ndpt_bmad = ",c_%ndpt_bmad
write(6,'(a16,1x,i2)')  " c_%np = ",c_%np
write(6,'(a16,1x,i2)')  " c_%npara_fpp= ",c_%npara_fpp
write(6,'(a16,1x,i2)')  " c_%nd2t= ",c_%nd2t
write(6,'(a16,1x,i2)')  " c_%ndc2t= ",c_%ndc2t
write(6,'(a16,1x,i2)')  " c_%pos_of_delta= ",c_%pos_of_delta


write(6,'(a16,1x,i2)')  " c_%rf= ",c_%rf
write(6,'(a16,1x,L)')  " c_%setknob= ",c_%setknob
write(6,'(a16,1x,L)')  " c_%knob= ",c_%knob

  




call in_bmad_units()
!call in_ptc_units()
bmad_state=default+time0+nocavity0
bmad_state=bmad_state+modulation0+spin0

call init_all(bmad_state,no,np)
 
 
write(6,*) "Taylors created with 2 knobs : PTC with phase space maps "
write(6,*) " without a cavity and modulation on a magnet "
write(6,'(a16,1x,i2)')  " c_%no = ",c_%no
write(6,'(a16,1x,i2)')  " c_%nv = ",c_%nv
write(6,'(a16,1x,i2)')  " c_%nd = ",c_%nd
write(6,'(a16,1x,i2)')  " c_%nd2 = ",c_%nd2
write(6,'(a16,1x,i2)')  " c_%nd2harm = ",c_%nd2harm
write(6,'(a16,1x,i2)')  " c_%ndpt = ",c_%ndpt
write(6,'(a16,1x,i2)')  " c_%ndpt_bmad = ",c_%ndpt_bmad
write(6,'(a16,1x,i2)')  " c_%np = ",c_%np
write(6,'(a16,1x,i2)')  " c_%npara_fpp= ",c_%npara_fpp
write(6,'(a16,1x,i2)')  " c_%nd2t= ",c_%nd2t
write(6,'(a16,1x,i2)')  " c_%ndc2t= ",c_%ndc2t
write(6,'(a16,1x,i2)')  " c_%pos_of_delta= ",c_%pos_of_delta


write(6,'(a16,1x,i2)')  " c_%rf= ",c_%rf
write(6,'(a16,1x,L)')  " c_%setknob= ",c_%setknob
write(6,'(a16,1x,L)')  " c_%knob= ",c_%knob

!stop
allocate(an(4),bn(4))
an=0
bn=0
DC_ac=1.d0
A_ac=0
theta_ac=0 
 D_ac=1.d0

!!! the modulation with be completely "DA", ie, not really there
an(1)=0.0_dp
call alloc_modulation(p_QF1,an,bn,DC_ac,A_ac,theta_ac, D_ac)




call init_all(bmad_state,no,np)

write(6,*) "Taylors created by a knob : PTC with phase space maps "
write(6,*) "without a cavity but with modulation of a magnet "
write(6,'(a16,1x,i2)')  " c_%no = ",c_%no
write(6,'(a16,1x,i2)')  " c_%nv = ",c_%nv
write(6,'(a16,1x,i2)')  " c_%nd = ",c_%nd
write(6,'(a16,1x,i2)')  " c_%nd2 = ",c_%nd2
write(6,'(a16,1x,i2)')  " c_%nd2harm = ",c_%nd2harm
write(6,'(a16,1x,i2)')  " c_%ndpt = ",c_%ndpt
write(6,'(a16,1x,i2)')  " c_%ndpt_bmad = ",c_%ndpt_bmad
write(6,'(a16,1x,i2)')  " c_%np = ",c_%np
write(6,'(a16,1x,i2)')  " c_%npara_fpp= ",c_%npara_fpp
write(6,'(a16,1x,i2)')  " c_%nd2t= ",c_%nd2t
write(6,'(a16,1x,i2)')  " c_%ndc2t= ",c_%ndc2t
write(6,'(a16,1x,i2)')  " c_%pos_of_delta= ",c_%pos_of_delta


write(6,'(a16,1x,i2)')  " c_%rf= ",c_%rf
write(6,'(a16,1x,L)')  " c_%setknob= ",c_%setknob
write(6,'(a16,1x,L)')  " c_%knob= ",c_%knob

!stop 888

p=>ring%start
closed_orbit=0
call find_orbit_x(closed_orbit,bmad_state, 1.0e-7_dp, fibre1=p)
write(6,*) " closed orbit "
write(6,"(6(1x,g11.4))") closed_orbit

xs0=closed_orbit
 call get_length(ring,circ)
 nu_mod=1.404074285438127d0
!!!! set a modulation clock !!!!!!
xs0%ac(1)%om=twopi*nu_mod/circ ! (B1) ! differs from the first edition   
xs0%ac(1)%x=0.d0 ;       ! (B2) ! differs from the first edition   
write(6,*) " Modulation tune in radians =",circ*xs0%ac(1)%om
xs0%use_q=.true.
xs0%nac=1

!!! for my Windows interface. Not useful in general
my_eprobe=>xs0
!!!

call alloc(xs)
call alloc(ids)
ids=1
xs=xs0+ids

write(6,*) " Probe_8 : before, a PTC structure "
call print(xs)
call propagate(xs,bmad_state,fibre1=p)
write(6,*) " Probe_8 : after, a PTC structure "

call print(xs)
write(6,*) " c_damap : after, a PTC structure "

ids=xs
call print(ids)


pause 5694
 
 
no=2
nd=4
np=2
ndpt=6
ac_modulation=1
knob=.true.
call c_init_all(no,nd,np,ndpt,ac_modulation) 
 

write(6,*) " ***************************************** "
write(6,'(a16,1x,i2)')  " c_%no = ",c_%no
write(6,'(a16,1x,i2)')  " c_%nv = ",c_%nv
write(6,'(a16,1x,i2)')  " c_%nd = ",c_%nd
write(6,'(a16,1x,i2)')  " c_%nd2 = ",c_%nd2
write(6,'(a16,1x,i2)')  " c_%nd2harm = ",c_%nd2harm
write(6,'(a16,1x,i2)')  " c_%ndpt = ",c_%ndpt
write(6,'(a16,1x,i2)')  " c_%ndpt_bmad = ",c_%ndpt_bmad
write(6,'(a16,1x,i2)')  " c_%np = ",c_%np
write(6,'(a16,1x,i2)')  " c_%npara_fpp= ",c_%npara_fpp
write(6,'(a16,1x,i2)')  " c_%nd2t= ",c_%nd2t
write(6,'(a16,1x,i2)')  " c_%ndc2t= ",c_%ndc2t
write(6,'(a16,1x,i2)')  " c_%pos_of_delta= ",c_%pos_of_delta


write(6,'(a16,1x,i2)')  " c_%rf= ",c_%rf
write(6,'(a16,1x,L)')  " c_%setknob= ",c_%setknob
write(6,'(a16,1x,L)')  " c_%knob= ",c_%knob
stop
call ptc_end(graphics_maybe=0,flat_file=.false.)
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

