; HACK OF calc_BB_cse TO PRINT NAMELIST PARAMETERS

;************************************************************************
PRO runbbody
; Hack of bbody.calc_bb_cse to display in namelist format

; Star O6.5III in Sternberg+ 2003
star1e50 = 0
star6e49 = 0
star3e49 = 0
star1e49III = 0
star1e49 = 1
star1e48 = 0
star1e47 = 0
; FREQUENCY BANDS
nus = [2.418e12,2.418e14,3.288e15,5.946e15,13.16e15,100e15] ; in Hz - NOTE: Last value is just ~ inf
nus = double(nus)
; Sternberg 2003 from Vacca 1996
; AGES FROM SCHALLER ET AL 1992
if star1e50 then begin
   ; Fiddled from Schaller+ 1992, 120 Msolar, lifetime ~ 3 Myr
   ;T = 10.0^4.75                     ; K
   ;R  = 22.0                       ; Rsolar
   T = 56000
   R = 22.2
   nus = [3.288,5.946,13.16,1e2]*1e15 ; in Hz - NOTE: Last value is just ~ inf
   nus = double(nus)
endif
if star3e49 then begin
   ; Fiddled from Schaller+ 1992, 120 Msolar, lifetime ~ 3 Myr
   T = 10.0^4.65                     ; K
   R  = 20.8                       ; Rsolar
   ;T = 54000
   ;R = 22.2
   nus = [3.288,5.946,13.16,1e2]*1e15 ; in Hz - NOTE: Last value is just ~ inf
   nus = double(nus)
endif
if star6e49 then begin
   ; Fiddled from Schaller+ 1992, 120 Msolar, lifetime ~ 3 Myr
   T = 10.0^4.7                     ; K
   R  = 22.15                       ; Rsolar
   ;T = 54000
   ;R = 22.2
   nus = [3.288,5.946,13.16,1e2]*1e15 ; in Hz - NOTE: Last value is just ~ inf
   nus = double(nus)
endif
if star1e49III then begin
   T = 41250.0                     ; K
   R  = 14.8                       ; Rsolar
   nus = [3.288,5.946,13.16,1e2]*1e15 ; in Hz - NOTE: Last value is just ~ inf
   nus = double(nus)
endif
; Type 05V (ish) in Sternberg+ 2003 (O5V is 46120, 11.4)
if star1e49 then begin
   ; ~20Msolar, lifetime 9Myr
   ; OLD VALUES FROM GEEN+ 2015 - TOO HIGH T
   ;T = 46000.0                     ; K
   ;R  = 11.2                      ; Rsolar
   T = 39730
   R = 9.6
   nus = [3.288,5.946,13.16,1e2]*1e15 ; in Hz - NOTE: Last value is just ~ inf
   nus = double(nus)
endif
; Type B0V (ish) in Sternberg+ 2003 (B0V in Sternberg is T=33340)
if star1e48 then begin
   ; 15-ish, lifetime 12.75Myr
   T = 33700.0                     ; K
   R  = 8.3                      ; Rsolar
endif
; Type B1V
; omega scorpii is 26550K, 4.51Rsolar
; http://adsabs.harvard.edu/abs/2012ApJ...746..154P also Wikipedia
if star1e47 then begin
   ; omega scorpii is 26550K, 4.51Rsolar
   ; http://adsabs.harvard.edu/abs/2012ApJ...746..154P also Wikipedia
   ; 11Msolar (from the above?), 18Myr lifetime
   T =  28500                ; K
   R  =  4.7                 ; Rsolar
   nus = [3.288,5.946,13.16,1e2]*1e15 ; in Hz - NOTE: Last value is just ~ inf
   nus = double(nus)
endif


; Units
Rsolar = 6.955e10 ; in cm

; Converted stuff
surface = 4.0 * !pi * (R*Rsolar)^2.0

; Compute average Hui et al cross section, for a blackbody at given
; temperature, in the given energy interval, given in frequency (s-1) 
; The ionization frequencies are:
; HI:   3.288d15    s-1
; HeI:  5.946d15    s-1
; HeII: 1.316d16    s-1
; Add photoionization rate calculation with:
;               \int F \sigma d\nu
; Add photo-heating rate calculation with:
;               \int F \sigma (h\nu - ionizationEgy) d\nu
;------------------------------------------------------------------------
  ;init_device,xsize=600,ysize=800
  ;!p.multi=[0,1,3]  &  !p.charsize=2.5
  if not keyword_set(n) then n = 1000L
  ;if not keyword_set(species) then species = 'H'
  hp = 6.62606876d-27                         ; Planck const. in cgs
  erg2ev = 6.24150974d11                      ; Erg to eV conversion
  ;print,nu(0),nu(n-1),bb(0),bb(n-1)

  ; Values for each packet
  ngrp = 3 ; number of groups
  nspe = 3 ; number of species
  nphs = dblarr(ngrp)
  egys = dblarr(ngrp)
  csns = dblarr(ngrp,nspe)
  cses = dblarr(ngrp,nspe)
  spcs = ["H","HeI","HeII"]

  ; Run for each group, species
  for ig=0,ngrp-1 do begin 

     nu0 = nus[ig]
     nu1 = nus[ig+1]
     

     bbody, T, [nu0, nu1], n, nu, bb, /log, /frequency, /photon_count

     ; Calculate average photon energy:
     BB_int = int_tabulated( nu, bb, /DOUBLE )
     avg_egy = int_tabulated( nu, bb*nu  ) / BB_int * hp * erg2ev
     egys[ig] = avg_egy


     ; Calculate photon count per cm-2:
     ph_count = int_tabulated( nu, bb, /double,/sort ) ;;;;weird....
     nphs[ig] = ph_count*surface

     for is=0,nspe-1 do begin
        species = spcs[is]
        
        ; Calculate average cross section:
        crosssection_hui, nu, species, cs
        BB_cs_int = int_tabulated( nu, bb*cs, /DOUBLE )
        avg_csn = BB_cs_int/BB_int
        csns[ig,is] = avg_csn
  
        ; Calculate energy weighted cross section:
        BB_cs_egy_int = int_tabulated( nu, bb*cs*nu, /DOUBLE )
        avg_cse = BB_cs_egy_int / int_tabulated( nu, bb*nu,/double  )
        cses[ig,is] = avg_cse

     endfor
     
  endfor

  ; Print output
  compileme = commalist(nphs)
  print, "!TOTAL NSOURCE:", total(nphs),"photons/s"
  print, "rt_n_source=", commalist(nphs)
  print, "/"
  print, ""
  print, "&RT_GROUPS"
  print, "group_csn(1,:)=",commalist(csns[0,*])
  print, "group_cse(1,:)=",commalist(cses[0,*])
  print, "group_csn(2,:)=",commalist(csns[1,*])
  print, "group_cse(2,:)=",commalist(cses[1,*])
  print, "group_csn(3,:)=",commalist(csns[2,*])
  print, "group_cse(3,:)=",commalist(cses[2,*])
  print, "group_egy     =",commalist(egys)
  print, "spec2group    = 1, 2, 3"
  print, "/"


END
