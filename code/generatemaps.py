from pcraster import * 
## === FILE IMPORTS ===

setglobaloption("lddin")
setglobaloption("lddfill")

def read_maps(src=''):
    dem = readmap(src+"dem.map")
    n = readmap(src+"n.map")
    lai = readmap(src+"lai.map")
    rr = readmap(src+"rr.map")
    ch = readmap(src+"ch.map")
    per = readmap(src+"per.map")
    thetas1 = readmap(src+"thetas1.map")
    thetai1 = readmap(src+"thetai1.map")
    ksat1 = readmap(src+"ksat1.map")
    landunit = readmap(src+"landunits.map")
    mask = readmap(src+"mask.map")
    psi1 = readmap(src+"psi1.map")
    chanmask = readmap(src+"chanmask.map")

    return dem, n, lai, rr, ch, per, thetas1, thetai1, ksat1, landunit, mask, psi1, chanmask


### === PART 1 Catchment Characteristics ====
def create_catchment(dem, n, lai, rr, ch, per, thetas1, thetai1, ksat1, landunit, mask, psi1, chanmask):
    # lddmap = lddcreate(dem*mask, 1e20, 1e20, 1e20, 1e20)
    dem = dem * mask
    dem = lddcreatedem(dem,999999,999999,999999,999999)
    lddmap = lddcreate(dem*mask, 1e20, 1e20, 1e20, 1e20)
    outlet = pit(lddmap)
    grad = max(sin(atan(slope(dem*mask))), 0.001)
    idmap = nominal(mask)

    print("Reporting maps")

    report(dem*mask, "dem.map")
    report(lddmap, "ldd.map")
    report(landunit*mask, "landunits.map")
    report(outlet, "outlet.map")
    report(grad*mask, "grad.map")
    report(idmap, "id.map")
    report(n*mask, "n.map")
    report(lai*mask, "lai.map")
    report(rr*mask, "rr.map")
    report(ch*mask, "ch.map")
    report(per*mask, "per.map")

    print("Created maps")
# # #

# # ### === PART 2 Channel Data ===
# # # channel properties

def set_chan_properties(mask, chanmask, dem):
    chancoh = 10
    chanman = 0.2
    chanside = 0
    chanwidth = 2
    roadwidth = 6
    chanKsat = 20
    # #
    # # ## map generation
    # chanmask = chanmask/chanmask
    lddchan = lddcreate(chanmask*dem, 1e20,1e20,1e20,1e20)
    changrad = max(0.001, sin(atan(slope(chanmask*dem))))
    chancoh = chanmask * scalar(chancoh)
    chanman = chanmask * scalar(chanman)
    chanside = chanmask * scalar(chanside)
    chanwidt = chanmask * scalar(chanwidth)
    hmxinit = mask * 0
    chanksat = chanmask *chanKsat
    # # # #
    # # # ## map report
    report(lddchan, "lddchan.map")
    report(changrad, "changrad.map")
    report(chancoh, "chancoh.map")
    report(chanman, "chanman.map")
    report(chanside, "chanside.map")
    report(chanwidt, "chanwidt.map")
    report(hmxinit, "hmxinit.map")
    report(chanksat, "chanksat.map")

# # # #
# # # #### === PART 3 Soil Surface ====
# # # ## Already done: r, n
# # # ## pending: roadwidth
def soil_surface(mask, ksat1, thetas1, thetai1, psi1):
    zeromap = mask*0
    # # # #
    # # # ## map report
    report(zeromap, "crust.map")
    report(zeromap, "stone.map")
    report(zeromap, "comp.map")
    report(zeromap, "hard.map")
    report(zeromap, "bufferid.map")
    report(zeromap, "buffervol.map")
    # # #
    # # #
    # # #### === PART 4 Green & Ampt ===
    soildep1 = scalar(200)
    soildep2 = scalar(1000)
    ksat2 = 0.9*ksat1*mask
    # # #

    report(ksat1*mask, "ksat1.map")
    report(thetas1*mask, "thetas1.map")
    report(thetai1*mask, "thetai1.map")
    report(psi1*mask, "psi1.map")
    report(soildep1*mask, "soildep1.map")

    report(ksat2*mask, "ksat2.map")
    report(thetas1*mask, "thetas2.map")
    report(thetai1*mask, "thetai2.map")
    report(psi1*mask, "psi2.map")
    report(soildep2*mask, "soildep2.map")