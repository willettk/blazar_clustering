Pro karlen_tbtplot, ps = ps, snr = snr

;col_define, 39

dir = '~/Astronomy/proposals/redspirals/'
if keyword_set(ps) then ps_start, filename=dir+'karlen_tbtplot.eps',/color,/landscape, xs=10, ys=5, /encap
openr, rfile, dir+'reds.txt', /get_lun

;Setting SNR
if n_elements(snr) eq 0 then snr = 2
x = snr

line=''
i=0

n=500

;arrays for reading in data
ID=strarr(n)
NII=dblarr(n)
OIII=dblarr(n)
HB=dblarr(n)
HA=dblarr(n)
NEIII=dblarr(n)
OII=dblarr(n)
NEIIIerr=dblarr(n)
OIIerr=dblarr(n)
HAerr=dblarr(n)
HBerr=dblarr(n)
NIIerr=dblarr(n)
OIIIerr=dblarr(n)
gz=dblarr(n)				; (g - z) color

;reading in data
while (not EOF(rfile)) do begin
   readf,rfile,line
   arr=strsplit(line,/extract)
   ID[i]=arr[0]
   HB[i]=arr[7]
   OIII[i]=arr[3]
   HA[i]=arr[9]
   NII[i]=arr[1]
   NEIII[i]=arr[5]
   OII[i]=arr[11]
   HBerr[i]=arr[8]
   OIIIerr[i]=arr[4]
   HAerr[i]=arr[10]
   NIIerr[i]=arr[2]
   NEIIIerr[i]=arr[6]
   OIIerr[i]=arr[12]
   gz[i]=arr[13]
   i++
endwhile

;resizing arrays
temp=where(ID)
HB=HB(temp)
OIII=OIII(temp)
HA=HA(temp)
NII=NII(temp)
OII=OII(temp)
NEIII=NEIII(temp)
HBerr=HBerr(temp)
OIIIerr=OIIIerr(temp)
HAerr=HAerr(temp)
NIIerr=NIIerr(temp)
NEIIIerr=NEIIIerr(temp)
gz=gz(temp)

id = id[temp]

lowzids =[ '587733603730522222', '587738568705507359', '587741488761274453', '587738410327801941', '587736543096734014', '587741602564014184', '587734949666881571', '587725551199780967', '588017978911817806', '587736783606382788', '587738066727534662', '587737827288875214', '587739630096023592', '587733397573009582', '587722982822379684', '587730772799914495', '587734894361706638']

match, lowzids, id, ind1, ind2, count=mcount

;setting NeIII fluxes smaller than 1sigma to upper limit: 1sigma
j=0
while (j lt n_elements(NEIII)) do begin
   if(NeIII[j] lt NeIIIerr[j]) then begin
      NeIII[j]=NeIIIerr[j]
   endif
   j++
endwhile

;calculating line ratios
OHB=alog10(OIII/HB)
NHA=alog10(NII/HA)
NEO=alog10(NEIII/OII)

;calculating nsigma for lines
sigOII=OII/OIIerr
sigNE=NEIII/NEIIIerr
sigHB=HB/HBerr
sigNII=NII/NIIerr


;Define x-sigma sample for BPT & TBT

sig=where((sigOII ge x) and (sigHB ge x) and (sigNII ge x) and (sigNE ge x))
sigOII=where(sigOII ge x)
sigNE=where(sigNE ge x)

ID = ID (sig)
NEO=NEO(sig)
gz=gz(sig)
NHA=NHA(sig)
OHB=OHB(sig)

;for error size calculations

NEIII=NEIII(sig)
NEIIIerr=NEIIIerr(sig)
OII=OII(sig)
OIIerr=OIIerr(sig)

;separating by bpt type, while tracking for null arrays

count=0
counta=0
counts=0

agn=where(OHB gt (0.61/(NHA - 0.47) + 1.19),counta)
sf=where(OHB lt (0.61/(NHA - 0.05) + 1.3) and NHA lt -0.2,counts)
comp=where((OHB le (0.61/(NHA - 0.47) + 1.19)) and (OHB ge (0.61/(NHA - 0.05) + 1.3)),count)

if(counts gt 0) then begin
NEOs = NEO(sf)
gzs = gz(sf)
NHAs = NHA (sf)
OHBs = OHB (sf)
endif

if(counta gt 0) then begin
NEOa = NEO(agn)
gza = gz(agn)
NHAa = NHA (agn)
OHBa = OHB (agn)
endif

if(count gt 0) then begin
   NEOc = NEO(comp)
   gzc = gz(comp)
   NHAc = NHA (comp)
   OHBc = OHB (comp)
endif

;for plotting Kewley and Kauffman lines
i=0
xn=dblarr(1800)
yn=dblarr(1800)
xs=dblarr(1800)
ys=dblarr(1800)
xo=dblarr(1800)
yo=dblarr(1800)
while (i lt 1800) do begin
   ;Kauffman line
   xn[i]=i/1000.-2
   yn[i]=(0.61/(xn[i] - 0.05)) + 1.30
   ;Kewley line
   xs[i] = i/800. -2
   ys[i] = (0.61/(xs[i] - 0.47)) + 1.19
;   ys[i]=(.72/(xs[i]-.32))+1.3 ;Kewley [SII]/Ha line
   ;TBT line
   xo[i] = i/500. -2
   yo[i] = -1.2*xo[i] - 0.4
   i++
endwhile

!p.multi = [0,2,1]

if keyword_set(ps) then cs = 1 else cs = 1
th=4

;Plotting BPT
cgplot, nhaa, ohba, psym=1, xtitle='log ([NII]/H'+greek('alpha')+')',ytitle='log([OIII]/H'+greek('beta')+')', xrange=[-1,0.5],yrange=[-1,1.5], /nodata, charsize = cs, thick=th, ythick=th, xthick=th, color='black', title='BPT'

cgplot, /over, nhaa, ohba, psym=1, color='black' ;AGN
if(counts gt 0) then cgplot, /over, NHAs, OHBs, psym=4,color='purple' ;SF
if(count gt 0) then cgplot, /over, NHAc, OHBc, psym=17, color='green' ;COMP
cgplot, /over, xn, yn, thick = 5, color='red' ;Kauffman
cgplot, /over, xs, ys, thick = 5, color='blue'  ;Kewley
cgtext, -0.9, -0.2, /data, 'Star-forming',charsize=cs, charthick=3
cgtext, -0.1, -0.75, /data, 'Comp.',charsize=cs, charthick=3
cgtext, -0.8, 1.0, /data, 'Seyfert/AGN',charsize=cs, charthick=3

;Plotting TBT
cgplot, NEOa, gza, psym=1, xtitle='log ([NeIII]/[OII])', ytitle='!E0!N(g - z)', xrange=[-2.0, 0.5], yrange=[-0.0, 2.0], /nodata, charsize = cs, thick=th, ythick=th, xthick=th, color='black', title='TBT'

cgplot, /over, NEOa, gza, color='black', psym=1 ;BPT AGN
cgplot, /over, xo, yo, thick=5, color='goldenrod'         ;TBT line
if(counts gt 0) then cgplot, /over, NEOs, gzs, psym=4, color='purple' ;BPT SF
if(count gt 0) then cgplot, /over, NEOc, gzc, psym=17, color='green'  ;BPT COMP
cgtext, -1.8, 0.5, /data, 'Star-forming',charsize=cs, charthick=3
cgtext, -0.5, 1.8, /data, 'Seyfert/AGN',charsize=cs, charthick=3

free_lun, rfile

if keyword_set(ps) then ps_end

;device, /close_file
;set_plot, my_device
end
