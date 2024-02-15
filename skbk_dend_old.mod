TITLE BK Ca 2+ -activated K + channel
: Calcium activated K channel.
COMMENT
 Starting from the formulation in De Schutter and Bower, 1994, we
reduced the Ca 2+ dependent activation time to half to account for the larger slow repolarisation at
depolarised states.
Current Model Reference: Karima Ait Ouares , Luiza Filipis , Alexandra Tzilivaki , Panayiota Poirazi , Marco Canepari (2018) Two distinct sets of Ca 2+ and K + channels 
are activated at different membrane potential by the climbing fibre synaptic potential in Purkinje neuron dendrites. 

PubMed link: 

Contact: Filipis Luiza (luiza.filipis@univ-grenoble-alpes.fr)

ENDCOMMENT

UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX bk2
	:USEION cal READ ical
	:USEION cah READ icah
	:USEION car READ icar
	USEION k READ ek WRITE ik
	USEION can READ cani
	USEION car2 READ car2i
	USEION cal2 READ cal2i
	USEION cat READ cati
	USEION capq READ capqi
	RANGE gkbar,gk,zinf,ik,minf,f,minf2, mexp,zexp, a,zinf2,a2
	:GLOBAL bkcoef,bkexp,cahco, carco,tin, tinsh
	GLOBAL bkcoef,bkexp,canco, car2co, calco, catco, capqco, tin, tinsh,alp0,test,cutoffz,kzexp,z0,alp2_0
}


PARAMETER {
	celsius=37	(degC)
	v		(mV)
	gkbar=.08	(mho/cm2)	: Maximum Permeability
	 :cai = .04e-3	(mM)
	 :ical (mA/cm2)
	 :icah (mA/cm2)
	 :icar (mA/cm2)
	cani (mM)
	car2i (mM)
	cati (mM)
	cal2i (mM)
	capqi (mM)
	dt		(ms)
	tin=0.02
	bkcoef=0.023
	bkexp=7.5
	:cahco=2
	:carco=1
	tinsh=400
	canco=1
	car2co=1
	calco=1
	catco=1
	capqco=1
	alp0=4
	a0=1
	f1=1000 	:was 1; factor to decrease Ca-dependent activation
	f2=0.8	:was 1; factor to influence Ca-dependence (both inact and act)
	kb=14.9
	b0=0.11
	vhm=10
	test=1
	cutoffz=10000e-5
	kzexp=4
	z0=0
	alp2_0=4
	vhmexp=20
}


ASSIGNED {
	ik		(mA/cm2)
	minf
	minf2
	mexp
	zinf
	zinf2
	zexp
	gk
	ek 	(mV)
}

STATE {	m z }		: fraction of open channels

BREAKPOINT {
	SOLVE state
:	gk = gkbar*1000*m*z*z
	ik = gkbar*1000*m*z*z*(v - ek)
:	ik = gkbar*1000*m*z*(v - ek)
}
:UNITSOFF
:LOCAL fac

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE state() {	: exact when v held constant; integrates over dt step
	:rate(v, ical, icah, icar)
	rate(v, cani,car2i,cal2i,cati,capqi)
	m = m + mexp*(minf - m)
	z = z + zexp*(zinf - z)
	VERBATIM
	return 0;
	ENDVERBATIM
}

INITIAL {
	:rate(v, ical, icah, icar)
	rate(v, cani,car2i,cal2i,cati,capqi)
	m = minf
	z = zinf
}

FUNCTION alp(v (mV), cani (mM),car2i (mM),cal2i (mM),cati (mM),capqi (mM)) (1/ms) { :callable from hoc
	alp = -alp0/((-cani*canco-car2i*car2co-cal2i*calco-cati*catco-capqi*capqco))
	
}

FUNCTION alp2(v (mV), cani (mM),car2i (mM),cal2i (mM),cati (mM),capqi (mM)) (1/ms) {
alp2 = -alp2_0/((-cani*canco-car2i*car2co-cal2i*calco-cati*catco-capqi*capqco))
}

FUNCTION bet(v (mV)) (1/ms) { :callable from hoc
	bet = b0/exp((v-55)/kb)
}

PROCEDURE rate(v (mV), cani (mM),car2i (mM),cal2i (mM),cati (mM),capqi (mM)) { :callable from hoc
	LOCAL a,a2,b,tinca
	a = alp(v, cani,car2i,cal2i,cati,capqi)
	a2 = alp2(v, cani,car2i,cal2i,cati,capqi)
	:zinf = 1/(1+(a*f2))
	zinf = z0+1/(1+pow(a,test))
	:zinf2=z0+ pow(a2,-1)
	:if (zinf>cutoffz){zinf=cutoffz}
	zexp = (1 - exp(-0.01/kzexp))
	b = bet(v)
	:minf = 8.5/(7.5+b)
	:minf = 1/(1+b)
	minf=1/(1+exp(-(v+vhm)/kb))
	:minf = 1/(b+(a/f1)) : Laila Blomer - to make activation Ca-dependent
	mexp =a0*(1/(1+exp(-(v+vhmexp)/bkexp))) :a0*(1 - exp(-0.01*(bkexp+b))) 
	:printf("v=%15.10g a=%15.10g b=%15.10g mexp=%15.10g zexp=%15.10g\n",v,a,b,mexp,zexp)
	:printf("cani=%15.10g cati=%15.10g cal2i=%15.10g car2i=%15.10g capqi=%15.10g\n",cani,cati,cal2i,car2i,capqi)
}
:UNITSON
