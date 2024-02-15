: From 'Dentate gyrus network model' (Santhakumar et al 2005) https://senselab.med.yale.edu/ModelDB/ShowModel?model=51781&file=/dentategyrusnet2005/ccanl.mod#tabs-1

COMMENT
  calcium accumulation into a volume of area*depth next to the
  membrane with a decay (time constant tau) to resting level
  given by the global calcium variable cai0_ca_ion
ENDCOMMENT

NEURON {
  SUFFIX cdp5r_dend
USEION ca WRITE cai
:USEION nca READ inca WRITE ncai VALENCE 2
:USEION lca READ ilca WRITE lcai VALENCE 2
:USEION tca READ itca WRITE tcai VALENCE 2
USEION cal2 READ ical2 WRITE cal2i VALENCE 2
USEION can READ ican WRITE cani VALENCE 2
USEION car2 READ icar2 WRITE car2i VALENCE 2
USEION cat READ icat WRITE cati VALENCE 2
USEION capq READ icapq WRITE capqi VALENCE 2
:RANGE caiinf, caitot, ncai, lcai,tcai, ecatot, elca, enca, etca
RANGE caiinf, caitot, ecatot
RANGE dlcai,dncai,drcai,dtcai,dpqcai
RANGE ecal, ecan, ecar2, ecat, ecapq
GLOBAL catau
}

UNITS {
        (mV) = (millivolt)
  (molar) = (1/liter)
  (mM) = (milli/liter)
  (mA) = (milliamp)
  FARADAY = 96520 (coul)
  R = 8.3134  (joule/degC)
}

INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

PARAMETER {
        celsius = 6.3 (degC)
  depth = 200 (nm)  : assume volume = area*depth
  catau = 0.1 (ms)
  caiinf = 50.e-6 (mM)  : takes precedence over cai0_ca_ion
      : Do not forget to initialize in hoc if different
      : from this default.
  caotot = 2 (mM)
  ica (mA/cm2)
  :inca (mA/cm2)
  :ilca (mA/cm2)
  :itca (mA/cm2)
  ical2 (mA/cm2)
  ican (mA/cm2)
  icar2 (mA/cm2)
  icat (mA/cm2)
  icapq (mA/cm2)
  caitot= 50.e-6 (mM)
  :cai= 50.e-6 (mM) :eg
  n_ions=5
}

ASSIGNED {
  :enca (mV)
  :elca (mV)
  :etca (mV)
  ecal (mV)
  ecan (mV)
  ecar2 (mV)
  ecat (mV) 
  ecapq (mV)
  ecatot (mV)
}

STATE {
  :ncai (mM)
  :lcai (mM)
  :tcai (mM)
  cal2i (mM)
  cani (mM)
  car2i (mM)
  cati (mM)
  capqi (mM)
  cai (mM)
}

INITIAL {
  cai=50.e-6
  VERBATIM
  cal2i = _ion_cal2i;
  cani = _ion_cani;
  car2i = _ion_car2i;
  cati = _ion_cati;
  capqi = _ion_capqi;
  ENDVERBATIM
  
  :ncai=caiinf/n_ions
  :lcai=caiinf/n_ions
  :tcai=caiinf/n_ions
  cal2i = caiinf/n_ions
  cani = caiinf/n_ions
  car2i = caiinf/n_ions
  cati = caiinf/n_ions
  capqi = caiinf/n_ions
  
  :caitot = caiinf
  cai = caiinf  :eg
  ecatot = ktf() * log(caotot/caiinf) 
  
  :enca = ecatot
  :elca = ecatot
  :etca = ecatot
  ecal = ecatot
  ecan = ecatot
  ecar2 = ecatot
  ecat = ecatot
  ecapq = ecatot
}


BREAKPOINT {
  SOLVE integrate METHOD derivimplicit
  :caitot = ncai+lcai+tcai  
  :caitot = cali+cani+car2i+cati+capqi 
  cai= cal2i+cani+car2i+cati+capqi
  :ecatot = ktf() * log(caotot/caitot)
  ecatot = ktf() * log(caotot/cai)  
  :enca = ecatot
  :elca = ecatot
  :etca = ecatot
  ecal = ecatot
  ecan = ecatot
  ecar2 = ecatot
  ecat = ecatot
  ecapq = ecatot
}

DERIVATIVE integrate {

cal2i' = -(ical2)/depth/FARADAY * (1e7) + (caiinf/n_ions - cal2i)/catau
cani' = -(ican)/depth/FARADAY * (1e7) + (caiinf/n_ions - cani)/catau
car2i' = -(icar2)/depth/FARADAY * (1e7) + (caiinf/n_ions - car2i)/catau
cati' = -(icat)/depth/FARADAY * (1e7) + (caiinf/n_ions - cati)/catau
capqi' = -(icapq)/depth/FARADAY * (1e7) + (caiinf/n_ions - capqi)/catau
cai'=-(ical2+ican+icar2+icat+icapq)/depth/FARADAY * (1e-7) + (caiinf - cai)/catau
}

FUNCTION ktf() (mV) {
  ktf = (1000)*R*(celsius +273.15)/(2*FARADAY)
} 