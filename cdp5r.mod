: From 'Dentate gyrus network model' (Santhakumar et al 2005) https://senselab.med.yale.edu/ModelDB/ShowModel?model=51781&file=/dentategyrusnet2005/ccanl.mod#tabs-1

COMMENT
  calcium accumulation into a volume of area*depth next to the
  membrane with a decay (time constant tau) to resting level
  given by the global calcium variable cai0_ca_ion
ENDCOMMENT

NEURON {
  SUFFIX cdp5r
USEION ca WRITE cai
USEION nca READ inca WRITE  ncai VALENCE 2
USEION lca READ ilca WRITE  lcai VALENCE 2 :eg
USEION tca READ itca WRITE tcai VALENCE 2
RANGE caiinf, caitot, ncai, lcai,tcai, ecatot, elca, enca, etca,cai
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
  inca (mA/cm2)
  ilca (mA/cm2)
  itca (mA/cm2)
  caitot= 50.e-6 (mM)
  :cai= 50.e-6 :eg
  n_ions=3
}

ASSIGNED {
  enca (mV)
  elca (mV)
  etca (mV)
  ecatot (mV)
}

STATE {
  ncai (mM)
  lcai (mM)
  tcai (mM)
  cai (mM)
}

INITIAL {
  cai= 50.e-6
  VERBATIM
  ncai = _ion_ncai;
  tcai = _ion_tcai; 
  lcai = _ion_lcai; 
  ENDVERBATIM
  ncai=caiinf/n_ions
  lcai=caiinf/n_ions
  tcai=caiinf/n_ions
  :caitot = caiinf 
  cai = caiinf  
  ecatot = ktf() * log(caotot/caiinf) 
  enca = ecatot
  elca = ecatot
  etca = ecatot
}


BREAKPOINT {
  SOLVE integrate METHOD derivimplicit
  :caitot = ncai+lcai+tcai  
  cai = ncai+tcai+lcai 
  :ecatot = ktf() * log(caotot/caitot) 
  ecatot = ktf() * log(caotot/caiinf)    
  enca = ecatot
  :elca = ecatot
  etca = ecatot
}

DERIVATIVE integrate {
ncai' = -(inca)/depth/FARADAY * (1e7) + (caiinf/n_ions - ncai)/catau
lcai' = -(ilca)/depth/FARADAY * (1e7) + (caiinf/n_ions - lcai)/catau
tcai' = -(itca)/depth/FARADAY * (1e7) + (caiinf/n_ions - tcai)/catau
cai'= -(ilca+inca+itca)/depth/FARADAY * (1e7) + (caiinf -cai)/catau
}

FUNCTION ktf() (mV) {
  ktf = (1000)*R*(celsius +273.15)/(2*FARADAY)
} 