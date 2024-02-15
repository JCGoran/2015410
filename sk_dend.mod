TITLE Calcium activated Potassium channel (SK)  

: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 


:SK channel ( function of Calcium only) Adapted by Luiza Filipis
: undone changes from Luiza by Laila Blomer to normalise channel activation kinetics dependent on Ca

NEURON {
	SUFFIX sk2
	USEION k READ ek WRITE ik
        : USEION cah READ icah
        : USEION car READ icar
		USEION can READ cani
		USEION car2 READ car2i
		USEION cal2 READ cal2i
		USEION cat READ cati
		USEION capq READ capqi
        RANGE  gbar,gkahp,ik, ninf,taun,g
        GLOBAL Cq10,mt, a0, b0
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(pS) = (picosiemens)
	(um) = (micron)
}

PARAMETER {
	gbar = 0.01	(S/cm2)
       : n = 4
        cai = 50.e-6	(mM)
        a0 = 3e4	(1/ms-mM)	:1.3e4	
		b0 = 0.06	(1/ms)		:b0 = 0.06	(1/ms)		keeps inactivation slow	
	    celsius = 37(degC)
	Cq10 = 1		:2
	:mt = 0.2
	mt = 0.8	:1
	:icah (mA/cm2)
	:icar (mA/cm2)
	cani (mM)
	car2i (mM)
	cati (mM)
	cal2i (mM)
	capqi (mM)
	:cahco=0.02
	:carco=1
	canco=1
	car2co=1
	calco=1
	catco=1
	capqco=1
	
}

STATE {	n }

ASSIGNED {
	ik	(mA/cm2)
	g	(S/cm2)
	ninf
	taun	(ms)
	a	(1/ms)
	v	(mV)
	ek	(mV)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gbar*n
	ik = g*(v-ek)
}

INITIAL {
	rate(cani,car2i,cal2i,cati,capqi)
	n=ninf
}

DERIVATIVE state {
	rate(cani,car2i,cal2i,cati,capqi)
	n' = (ninf - n)/taun
}

PROCEDURE rate(cani,car2i,cal2i,cati,capqi) {
	LOCAL q10
	q10 = Cq10^((celsius - 22 (degC))/10 (degC) )
	:a = a0*(-ican*canco-icar2*car2co-ical2*calco-icat*catco-icapq*capqco)/10
	:a = a0*(-cani*canco-car2i*car2co-cal2i*calco-cati*catco-capqi*capqco)^4
	a = a0*(cani*canco+car2i*car2co+cal2i*calco+cati*catco+capqi*capqco)^4
	:taun = b0
	taun = mt * q10* (0.2*(q10/(a + b0)))
	:ninf = a
	ninf = a/(a+b0)
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	       trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	

