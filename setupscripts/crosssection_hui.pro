;************************************************************************
PRO crosssection_hui, nu, species, cs, eV=eV

; Gives an atom-photon cross-section of given species at given frequency
; nu      => Frequency in s-1
; species => 'H', 'HeI' or 'HeII'
; eV      => If set then nu is interpreted as being electon volts
; cs     <=  Returned cross-section in cm^2
;------------------------------------------------------------------------
  h = 6.62606876d-27                          ; Planck const. in cgs
  erg = 6.2415d11                             ; eV's per erg
  nu_hz = nu
  if keyword_set(eV) then nu_hz = nu/h/erg
  E = h * nu_hz * erg                         ; photon energies in ev
  case species of
     'H': begin
        numin = 3.288e15
        E0=4.298d-1 & cs0=5.475d-14 & P=2.963 
        ya=32.88 & yw=0 & y0=0 & y1=0
     end
     'HeI': begin 
        numin = 5.946e15
        E0=1.361d1 & cs0=9.492d-16 & P=3.188 
        ya=1.469 & yw=2.039 & y0=0.4434 & y1=2.136
     end
     'HeII': begin      
        numin = 1.316d16
        E0=1.720 & cs0=1.369d-14 & P=2.963 
        ya=32.88 & yw=0 & y0=0 & y1=0
     end
  endcase
  x = E/E0 - y0
  y = sqrt(x^2+y1^2)
  cs = cs0 * ((x-1.)^2 + yw^2) * y^(0.5*P-5.5)/(1.+sqrt(y/ya))^P
  underlimit_ind = where(nu_hz lt numin, count)
  ;print,underlimit_ind
  if(count gt 0) then begin 
     cs(underlimit_ind)=0.
     ;cs((where(nu_hz ge numin))(0))=1.d-200
  endif
END

