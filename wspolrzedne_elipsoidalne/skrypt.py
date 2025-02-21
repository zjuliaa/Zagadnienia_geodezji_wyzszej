import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import os
from geopy.distance import geodesic

def read_flightradar(file):
    '''
    Parameters
    ----------
    file : .csv file - format as downloaded from fligthradar24
        DESCRIPTION.
    Returns
    -------
    all_data : numpy array
        columns are:
            0 - Timestamp - ?
            1 - year
            2 - month
            3 - day
            4 - hour
            5 - minute
            6 - second
            7 - Latitude [degrees]
            8 - Longitude [degrees]
            9 - Altitude [feet]
            10 - Speed [?]
            11 - Direction [?]
    '''
    with open(file, 'r') as f:
        i = 0
        size= []
        Timestamp = []; date = []; UTC = []; Latitude = []; Longitude = []; 
        Altitude = []; Speed = []; Direction = []; datetime_date = []
        for linia in f:
            if linia[0:1]!='T':
                splited_line = linia.split(',')
                size.append(len(splited_line))
                i+=1
                Timestamp.append(int(splited_line[0]))
                full_date = splited_line[1].split('T')
                date.append(list(map(int,full_date[0].split('-'))))
                UTC.append(list(map(int, full_date[1].split('Z')[0].split(':'))))
                Callsign = splited_line[2]
                Latitude.append(float(splited_line[3].split('"')[1]))
                Longitude.append(float(splited_line[4].split('"')[0]))
                Altitude.append(float(splited_line[5]))
                Speed.append(float(splited_line[6]))
                Direction.append(float(splited_line[7]))
                
    all_data = np.column_stack((np.array(Timestamp), np.array(date), np.array(UTC),
                                np.array(Latitude), np.array(Longitude), np.array(Altitude),
                                np.array(Speed), np.array(Direction)))
    return all_data

def blh2xyz(phi, lam, h):
    
    phi_rad = np.deg2rad(phi)
    lam_rad = np.deg2rad(lam)
    a = 6378137
    e2 = 0.00669438002290
    N = a / (np.sqrt(1 - e2 * np.sin(phi_rad) * np.sin(phi_rad)))
    X= (N + h) * np.cos(phi_rad) * np.cos(lam_rad)
    Y = (N + h) * np.cos(phi_rad) * np.sin(lam_rad)
    Z = (N * (1 - e2) + h) * np.sin(phi_rad)
    return X, Y, Z

def macierz_obrotu(phi, lam):
    R = np.array([[-np.sin(phi)*np.cos(lam), -np.sin(lam), np.cos(phi) * np.cos(lam)], 
                        [-np.sin(phi) * np.sin(lam), np.cos(lam), np.cos(phi) * np.sin(lam)],
                        [np.cos(phi), 0, np.sin(phi)]])
    return R

def calculate_distances(wspolrzedne_lot, lotnisko):
    distances = [geodesic(lotnisko, (lat, lon)).kilometers for lat, lon, alt in wspolrzedne_lot]
    return distances

if __name__ == '__main__':

    azymuty = []
    plik = 'C:/Sem3/GEODEZJA/projekt2/zajecia2/dane/lot25.csv'
    dane = read_flightradar(plik)

    for lot in dane:
        wspolrzedne = dane[ :, [7, 8, 9]]
        lot = np.where(wspolrzedne[:, -1]>0)[0]
        wspolrzedne[:, -1 ] = wspolrzedne[:, -1]*0.3048 + 135.4
        wspolrzedne_lot = wspolrzedne[lot, :]
        wspolrzedne_lotniska = wspolrzedne[lot[0]-1, :]
        Speed = dane[:, 10]
        Timestamp = dane[:, 0] 
        
        xyz_lotniska = blh2xyz(np.deg2rad(wspolrzedne_lotniska[0]), np.deg2rad(wspolrzedne_lotniska[1]), np.deg2rad(wspolrzedne_lotniska[2]) )
        R = macierz_obrotu(np.deg2rad(wspolrzedne_lotniska[0]), np.deg2rad(wspolrzedne_lotniska[1]))

        for flh in wspolrzedne_lot:
            xyz_samolotu = blh2xyz(np.deg2rad(flh[0]), np.deg2rad(flh[1]), np.deg2rad(flh[2]))
            wektor_samolot_lostnisko = np.array(xyz_samolotu) - np.array(xyz_lotniska)
            neu = macierz_obrotu(np.deg2rad(flh[0]), np.deg2rad(flh[1]))
            neu = R.T@wektor_samolot_lostnisko
            az = np.arctan2(neu[1], neu[0])
            azymuty.append(az)

    request = cimgt.GoogleTiles()
    fig = plt.figure(figsize=(10, 5))
    extent = [-12, 22, 50, 57]
    ax = plt.axes(projection=request.crs)
    ax.set_extent(extent)
    ax.add_image(request, 5)
    start_point = (wspolrzedne_lotniska[0], wspolrzedne_lotniska[1])
    end_point = (dane[-1, 7], dane[-1, 8])
    num_points = 100
    lats = np.linspace(start_point[0], end_point[0], num_points)
    longs = np.linspace(start_point[1], end_point[1], num_points)
    line_points = list(zip(lats, longs))
    ax.plot(longs, lats, transform=ccrs.PlateCarree(), color='orange', linewidth=2)

    for flh, az in zip(wspolrzedne_lot, azymuty):
        if 0 <= az < np.pi: 
            color = 'blue'
        else:
            color = 'red'
        ax.plot(flh[1], flh[0], transform=ccrs.PlateCarree(), color = color, marker='o', markersize=5)

    czas_od_startu = (Timestamp - Timestamp[0]) / 60

    fig2 = plt.figure(figsize=(10, 5))
    ax2 = fig2.add_subplot(1, 1, 1)
    ax2.plot(czas_od_startu, Speed)
    ax2.set_xlabel('Czas [minuty]')
    ax2.set_ylabel('Prędkość [węzły]')
    ax2.set_title('Wykres zmian prędkości lotu samolotu')

    fig3 = plt.figure(figsize=(10, 5))
    ax3 = fig3.add_subplot(1, 1, 1)
    ax3.plot(czas_od_startu, wspolrzedne[:, -1 ])
    ax3.set_xlabel('Czas [minuty]')
    ax3.set_ylabel('Wysokość [metry]')
    ax3.set_title('Wykres zmian wysokości lotu samolotu w zależności od czasu')

    distances = calculate_distances(wspolrzedne_lot, start_point)
    czas_trwania_lotu = (Timestamp[-1] - Timestamp[0]) / 60 
    fig4 = plt.figure(figsize=(10, 5))
    ax4 = fig4.add_subplot(1, 1, 1)
    ax4.plot(distances)
    ax4.set_xlabel('Czas [minuty]')
    ax4.set_ylabel('Dystans [kilometry]')
    ax4.set_title('Wykres zależności odległości od lotnsika startu od czasu')
    ax4.grid(True)

    plt.show()
