import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#import math

def dms2deg(dms): #zamiana minut i sekund na stopnie
    d = dms[0]
    m = dms[1]
    s = dms[2]
    
    deg = d+m/60+s/3600
    return deg

def deg2dms(dd): #zamiana stopni na minuty i sekundy
    deg = int(np.trunc(dd))
    mnt = int(np.trunc((dd-deg) * 60))
    sec = ((dd-deg) * 60 - mnt) * 60
    dms = [deg, abs(mnt), abs(sec)]
    # print(str(deg)+chr(176)+"%0.2d" % abs(mnt)+'\''+"%08.5f" % abs(secq)+'\"')
    return dms

def hms2rad(dms): #zamiana godzin minut sekund na radiany
    d = dms[0]
    m = dms[1]
    s = dms[2]
    
    deg = d+m/60+s/3600
    rad = np.deg2rad(deg*15)
    return rad

def dms2rad(dms): #zamiana stopni minut sekund na radiany
    d = dms[0]
    m = dms[1]
    s = dms[2]
    
    deg = d+m/60+s/3600
    rad = np.deg2rad(deg)
    return rad

def hms2sec(hms): #zamiana godzin minut sekund na sekundy
    sec = hms[0]*3600 + hms[1] * 60 + hms[2]
    return sec

def sec2hms(s):
    hd = s/3600
    h = int(np.trunc(hd))
    m = int(np.trunc((hd-h) * 60))
    s = ((hd-h) * 60 - m) * 60
    hms = [h,abs(m),abs(s)]
    return hms

def rad2hms(rad): #zamiana radianow n a na godziny, minuty i sekundy
    dd = np.rad2deg(rad)
    dd = dd/15
    deg = int(np.trunc(dd))
    mnt = int(np.trunc((dd-deg) * 60))
    sec = ((dd-deg) * 60 - mnt) * 60
    dms = [deg, abs(mnt), abs(sec)]
    return dms

def rad2dms(rad): #zamiana radianow na stopnie, minuty i sekundy
    dd = np.rad2deg(rad)
    dd = dd
    deg = int(np.trunc(dd))
    mnt = int(np.trunc((dd-deg) * 60))
    sec = ((dd-deg) * 60 - mnt) * 60
    dms = [deg, abs(mnt), abs(sec)]
    return dms

def dms2hms(dms): #zamiana stopni minut i sekund na godziny, minuty i sekundy
    sall = dms[0] * (4*60) + dms[1] * 4 + dms[2]/15    
    h = int(sall//3600)
    m = int((sall%3600)//60)
    s = sall%60
    return [h,m,s] 

def julday(y,m,d,h):
    '''
    Simplified Julian Date generator, valid only between
    1 March 1900 to 28 February 2100
    '''
    if m <= 2:
        y = y - 1
        m = m + 12
    # A = np.trunc(y/100)
    # B = 2-A+np.trunc(A/4)
    # C = np.trunc(365.25*y)
    # D = np.trunc(30.6001 * (m+1))
    # jd = B + C + D + d + 1720994.5
    jd = np.floor(365.25*(y+4716))+np.floor(30.6001*(m+1))+d+h/24-1537.5
    return jd


def GMST(jd):
    '''
    calculation of Greenwich Mean Sidereal Time - GMST in hours
    ----------
    jd : TYPE
        julian date
    '''
    T = (jd - 2451545) / 36525
    Tu = jd - 2451545
    g = 280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*T**2-T**3/38710000
    g = (g%360) / 15
    return g




if __name__ == '__main__':
    
    location1_lat: float = 52.0
    location1_lat_rad = np.deg2rad(location1_lat)
    location1_lon: float = 21.0
    location1_lon_rad = np.deg2rad(location1_lon)
    location2_lat: float = 0.0
    location2_lat_rad = np.deg2rad(location2_lat)
    location2_lon: float = 21.0
    location2_lon_rad = np.deg2rad(location2_lon)
 
    # declination =[77, 26, 13.24] 
    # right_ascension = [14, 8, 47.507]

    # declination =[27, 58, 2.980] 
    # right_ascension = [7, 46, 45.037]

    # declination =[-11, 16, 59.73]
    # right_ascension = [13, 26, 26.067]

    #slonce 
    declination =[23, 8, 11.85] 
    right_ascension = [6, 37, 43.973]

    #ksiezyc
    # declination =[-23, 54, 48.69] 
    # right_ascension = [16, 9, 45.978]

    n_hours = 24
    hours = np.arange(1, n_hours + 1)
    jd = julday(2023, 7, 1, hours)
    jd -= 2/24
    GMST0 = GMST(jd)
    #LST = GMST0 + location1_lon_rad / 15
    LST = GMST0 * 15 + location1_lon
    t = np.deg2rad(LST) - hms2rad(right_ascension)
    h = np.arcsin(np.sin(location1_lat_rad) * np.sin(dms2rad(declination)) + np.cos(location1_lat_rad) * np.cos(dms2rad(declination)) * np.cos(t))
    Az = np.arctan2(((-1) * np.cos(dms2rad(declination)) * np.sin(t)), (np.cos(location1_lat_rad) * np.sin(dms2rad(declination)) - np.sin(location1_lat_rad) * np.cos(dms2rad(declination)) * np.cos(t)))

    #SKYPLOT dla jednej gwiazdy
    fig1 = plt.figure(figsize = (8,8))
    ax = fig1.add_subplot(polar = True)
    ax.set_theta_zero_location('N') # ustawienie kierunku północy na górze wykresu
    ax.set_theta_direction(-1)
    ax.set_yticks(range(0, 90+10, 10))                   # Define the yticks
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    ax.set_rlim(0,90) 
    # narysowanie punktu na wykresie 
    ax.scatter(Az, 90-np.rad2deg(h))

   #SKYPLOT dla wielu gwiazd
    # fig = plt.figure(figsize=(8, 8))
    # ax = fig.add_subplot(polar=True)
    # ax.set_theta_zero_location('N')
    # ax.set_theta_direction(-1)
    # ax.set_yticks(range(0, 90 + 10, 10))
    # yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    # ax.set_yticklabels(yLabel)
    # ax.set_rlim(0, 90)
    # ax.scatter(Az, 90 - np.rad2deg(h), label='Star 1', color = 'orange')
    # ax.scatter(Az_sun, 90 - np.rad2deg(h_sun), label='Star 2', color = 'red')
    # ax.legend(loc='upper left')
    # plt.title("Skyplot for Two Stars")


    #Wykres 3D
    fig2 = plt.figure(figsize = (10,10))
    ax2 = fig2.add_subplot(projection = '3d')
    # promień Ziemi
    r = 1
    # siatka wspołrzędnych
    u, v = np.mgrid[0:(2 * np.pi+0.1):0.1, 0:np.pi:0.1]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)
    z[z<0] = 0		# bez tego, narysowalibyśmy całą kulę, a chcemy tylko półkulę
    ax2.plot_surface(x,y,z, alpha = 0.1)
    gx = r * np.sin(Az) * np.cos(h)
    gy = r * np.cos(Az) * np.cos(h)
    gz = r * np.sin(h)
    ax2.plot3D(gx,gy,gz, 'o')


    # Stworzenie wykresu
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    ax3.plot(hours, np.rad2deg(h), color = 'orange')
    
    # Dodanie linii reprezentującej zakres szerokości geograficznej
    ax3.axhline(y=0, color='blue', linestyle='--', label='Horyzont')
    
    ax3.set_xlabel('Godzina')
    ax3.set_ylabel('Wysokość (stopnie)')
    ax3.set_title('Zależność wysokości od czasu')
    ax3.grid(True)
    ax3.legend()
    ax3.set_xticks(hours)
    
    plt.show()

    
