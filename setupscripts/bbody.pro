;************************************************************************
PRO bbody, T, x_range, n, x, bb, frequency=frequency,                  $
            photon_count=photon_count, log=log
  
; Generates a blackbody spectrum as a function of frequency/wavelength 
; in cgs.
; T        => Temperature in Kelvin
; x_range  => 2-el. vector containing lower and upper limits of
;             x-axis, in cgs units
; n        => number of bins in returned spectrum
; x       <=  Returned x-axis (frequency or wavelength), of length n
; bb      <=  Returned blackbody spectrum, of length n
; frequency: If set, bb is function of frequency, otherwise wavelength
; photon_count: If set return the spectrum in units of a number of
;               photons, otherwise in units of ergs.
; log:       If set, return the spectrum on a logarithmic range.
; The units of the bb-spectrum is then:
; freq=1, ph=0:     erg s-1 cm-2 Hz-1
; freq=1, ph=1:      #  s-1 cm-2 Hz-1
; freq=0, ph=0:     erg s-1 cm-2 cm-1
; freq=0, ph=1:      #  s-1 cm-2 cm-1
;------------------------------------------------------------------------
  if (keyword_set(frequency)) then fr=1 else fr=0
  if (keyword_set(photon_count)) then ph=1 else ph=0
  c = 2.99792458d10                           ; light speed in cm s-2
  h = 6.62606876d-27                          ; Planck const. in cgs
  k = 1.3806504d-16                           ; Boltzmann const. in cgs
  sigma = 5.670d-5
  x=dblarr(n)
  bb=dblarr(n)
  if (keyword_set(log)) then begin
     x = 10.d^(alog10(x_range(0))+dindgen(n)/(n-1)*(alog10(x_range(1)/x_range(0))))
  endif else begin
     x = x_range(0)+dindgen(n)/(n-1)*(x_range(1)-x_range(0))
  endelse
  case fr of 
     0: begin
        case ph of
           0: begin
              bb = $            ; function of wavelength
                 2.d*h*c^2/x^5 * 1.d0/(exp(h/k*c/x*1.d/T)-1.d0)   ;    $
                 ;/  (sigma*T^4/!pi)
           end
           1: begin
              bb = $            ; function of wavelength
                 2.d*c/x^4 * 1.d0/(exp(h/k*c/x*1.d/T)-1.d0)
           end
        endcase
     end
     1: begin
        case ph of
           0: begin
              bb = $            ; function of frequency
                 2.d*h*x^3/c^2 * 1.d0/(exp(h/k*x/T)-1.d0)
           end
           1: begin
              bb = $            ; function of frequency
                 2.d/c^2 * x^2/(exp(h/k*x/T)-1.d0)
;                 2.d*x^2/c^2 * 1.d0/(exp(h/k*x/T)-1.d0)
           end
        endcase
     end     
  endcase
END

;************************************************************************
PRO plot_bbody, T, ps=ps, photon_count=photon_count

; Plot a blackbody spectrum for a given temperature
;------------------------------------------------------------------------
  if N_ELEMENTS(ps) eq 0 then window,0,/retain,xs=600,ys=600
  loadct,0
  !p.multi=[0,2,2]  &  !p.charsize=1.5
  n = 100L                                  ; number of bins
  c = 2.99792458d10                           ; light speed in cm s-2 
  la_range=1.d-8*[9.d1, 1.1d6]                ; wavelengt range (cm.)
  ;la_range=1.d-8*[1000, 15000]               ; wavelengt range (cm.)
  nu_range=c/[la_range(1), la_range(0)]

  title = 'Blackbody spectrum at '+strtrim(string(T))+' Kelvin'
  xtitle='!4k !3 [A]'
  ytitle='B!L!4k!3!N [erg s!U-1!N cm!U-2!N A!U-1!N]'
  bbody, T, la_range, n, la, bb, /log
  la = la * 1.d8 ; conversion to angstrom
  bb = bb / 1.d8
  plot,la,bb,title=title,xtitle=xtitle,ytitle=ytitle
  plot,la,bb,/xlog,/ylog,xtitle=xtitle,ytitle=ytitle;, $
  ;     yrange=[1.d-14, 1.d-2], ystyle=1

  xtitle='!4m!3 [Hz]'
  ytitle='B!L!4m!3!N [erg s!U-1!N cm!U-2!N Hz!U-1!N]'
  bbody, T, nu_range, n, nu, bb, /frequency 
  plot,nu,bb,xtitle=xtitle,ytitle=ytitle
  plot,nu,bb,/xlog,/ylog,xtitle=xtitle,ytitle=ytitle

END

;************************************************************************
PRO calc_BB_cse, T, nu0, nu1, n=n, species=species

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
  if not keyword_set(species) then species = 'H'
  hp = 6.62606876d-27                         ; Planck const. in cgs
  erg2ev = 6.24150974d11                      ; Erg to eV conversion
  bbody, T, [nu0, nu1], n, nu, bb, /log, /frequency, /photon_count
  print,nu(0),nu(n-1),bb(0),bb(n-1)
  ; Calculate average photon energy:
  BB_int = int_tabulated( nu, bb, /DOUBLE )
  avg_egy = int_tabulated( nu, bb*nu  ) / BB_int * hp * erg2ev
  
  ; Calculate average cross section:
  crosssection_hui, nu, species, cs
  BB_cs_int = int_tabulated( nu, bb*cs, /DOUBLE )
  avg_csn = BB_cs_int/BB_int
  
  ; Calculate energy weighted cross section:
  BB_cs_egy_int = int_tabulated( nu, bb*cs*nu, /DOUBLE )
  avg_cse = BB_cs_egy_int / int_tabulated( nu, bb*nu,/double  )

  ; Calculate photon count per cm-2:
  ph_count = int_tabulated( nu, bb, /double,/sort );;;;weird....

  print,'avg_egy=', avg_egy
  print,'avg_csn=', avg_csn
  print,'avg_cse=', avg_cse
  print,'ph_count=',ph_count
  print,'        =',ph_count/3.3985993d25 ;T=1d5 BB
  print,'        =',ph_count/1.2438543d27 ;T=3d5 BB

  ;xrg=[0,3.d16] & yrg=[0,7.d9]
  ;plot,nu,bb,title='BB spectrum',xrange=xrg,yrange=yrg;,/ylog
  ;yrg=[0,7.d-18]
  ;plot,nu,cs,title='Cross sections',xrange=xrg,yrange=yrg;,/ylog
  ;yrg=[0,1.d-7]
  ;plot,nu,bb*cs,title='BB * CS',xrange=xrg,yrange=yrg;,/ylog


END

