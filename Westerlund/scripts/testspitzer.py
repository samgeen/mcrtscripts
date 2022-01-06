
def SoundSpeed(cloud,star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = ionisedtemperatures.FindTemperature(star.Teff(t),cloud.metal)
    # Find ci
    # 2.0 factor is because of free electrons
    ci = np.sqrt(2.0 * gamma * kB * Ti * X / mH)
    return ci